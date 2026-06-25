from pylab import *
import os
import time
import accretion as setup
import galce as galce
from mock import mock
from mock_migration import mock_migration
from surface_density import get_num

'''
This Galactic chemical evolution model is an updated version of the framework
developed by Lian et al. (2020). The original model is extended to a multi-zone
formulation and includes stellar radial migration following Frankel et al. (2018).
The model adopts seven annuli centered at R = 2, 4, ..., 14 kpc, corresponding
approximately to the disk region from 1 to 15 kpc, and includes two main gas
accretion episodes. The evolutionary history is divided into four phases:
initial, secular, second, and postburst.
'''

def main():
    start_time = time.time()

    #*** Radial and temporal grids ***
    n_radial_bins = 7
    birth_radii_kpc = np.array([2, 4, 6, 8, 10, 12, 14])
    migration_sigma_value = 3
    migration_sigma_matrix = np.full((n_radial_bins, n_radial_bins), migration_sigma_value) #Radial migration width parameter; rows correspond to birth radii and columns to present-day radial bins.

    time_step_gyr = 0.02
    universe_age_gyr = 13.7
    gas_accretion_history = setup.gasAcc(time_step_gyr, universe_age_gyr)

    #*** Gas accretion and star formation efficiency histories ***
    #The seven entries in each array correspond to the annuli at R = 2, 4, ..., 14 kpc.

    #Gas accretion rates in the different evolutionary phases.
    accretion_initial = np.array([0.13, 0.11, 0.07, 0.05, 0.03, 0.006, 0.004])
    accretion_second = np.array([0, 0.03, 0.07, 0.12, 0.15, 0.2, 0.25])
    accretion_postburst = np.array([0, 0, 0, 0, 0, 0, 0])

    #Star formation efficiency coefficients in the different phases; these act as multiplicative factors in the Kennicutt-Schmidt law.
    sfe_initial = np.array([1.4, 1.4, 1.5, 1.6, 1.8, 2, 2.4]) 
    sfe_secular = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
    sfe_second = np.array([0.8, 0.8, 0.8, 0.8, 0.8, 0.6, 0.4])
    sfe_postburst = np.array([0.15, 0.14, 0.13, 0.12, 0.08, 0.08, 0.05])

    #Phase transition times, in Gyr.
    initial_phase_end_gyr = 2.5
    burst_start_gyr = 7.4
    burst_duration_gyr = 0.5

    #Exponential decay timescales.
    tau_sfe = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]) #E-folding timescale for the decline of the SFE coefficient during the secular phase.
    tau_accretion_initial = np.array([100, 100, 100, 100, 100, 100, 100]) #E-folding timescale for the initial-phase gas accretion rate.
    tau_accretion_secular = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]) #E-folding timescale for the secular-phase gas accretion rate.
    tau_accretion_second = np.array([100, 100, 100, 100, 100, 100, 100]) #E-folding timescale for the second-phase gas accretion rate.

    #Generate the gas accretion and SFE histories.
    gas_accretion_history.mw_multi_burst(
        accretion_initial,
        tau_accretion_initial,
        tau_accretion_secular,
        accretion_second,
        tau_accretion_second,
        accretion_postburst,
        sfe_initial,
        tau_sfe,
        sfe_secular,
        sfe_second,
        sfe_postburst,
        initial_phase_end_gyr,
        burst_start_gyr,
        burst_duration_gyr,
    )

  
    #*** Physical parameters for the chemical evolution model ***
    initial_gas_mass = 0 #Initial gas mass.
    imf_type = "kr" #Initial mass function type; "kr" denotes a Kroupa IMF in galCE.
    outflow_mass_loading = 0 #Gas outflow parameter.
    ks_index = 1.4 #Power-law index of the Kennicutt-Schmidt relation.
    ccsn_yield_set = "l18" #Core-collapse supernova yield set; "l18" corresponds to the Limongi & Chieffi (2018) table.
    ia_dtd_power = -1 #Power-law slope of the SN Ia delay-time distribution.
    ia_dtd_min_delay_gyr = 0.15 #Minimum SN Ia delay time, in Gyr.
    burst_infall_abundance = -0.65 #Metallicity [M/H] of the gas accreted during the second phase.

  
    #*** Mock catalog and radial migration parameters ***
    #Different age uncertainties are assigned to young and old populations to better reproduce the observations.
    age_split_gyr = 8 #Age boundary between the young and old samples, in Gyr.
    young_age_error_dex = 0.15 #Age uncertainty width for young stars.
    old_age_error_dex = 0.05 #Age uncertainty width for old stars.

    sample_numbers = get_num("../mw-sbp.txt") #Number of mock stars sampled in each radial bin, based on Lian et al. (2025).

    output_dir = "../model/" #Output directory for the model products.
    os.makedirs(output_dir, exist_ok=True)

    model_filenames = [
        "model_R_2.fits",
        "model_R_4.fits",
        "model_R_6.fits",
        "model_R_8.fits",
        "model_R_10.fits",
        "model_R_12.fits",
        "model_R_14.fits",
    ]

    #*** Run the chemical evolution model ***
    evolution_model = galce.galCE(
        output_dir,
        model_filenames,
        gas_accretion_history,
        initial_gas_mass,
        imf_type,
        outflow_mass_loading,
        time_step_gyr,
        ks_index,
        ccsn_yield_set,
        ia_dtd_power,
        ia_dtd_min_delay_gyr,
        n_radial_bins,
        birth_radii_kpc,
        migration_sigma_matrix,
        burst_infall_abundance,
        burst_start_gyr,
        burst_duration_gyr,
    )
    evolution_model.run()

    #*** Generate mock stellar catalogs from the model tracks ***
    age_error_split = np.full(n_radial_bins, age_split_gyr)
    for radial_index in range(n_radial_bins):
        mock(
            output_dir,
            model_filenames[radial_index],
            sample_numbers[radial_index],
            young_age_error_dex,
            old_age_error_dex,
            age_error_split[radial_index],
        )

    #*** Apply radial migration to the mock stellar catalogs ***
    radial_file_prefixes = ["R_2", "R_4", "R_6", "R_8", "R_10", "R_12", "R_14"]
    mock_migration(output_dir, radial_file_prefixes, migration_sigma_matrix, birth_radii_kpc)

    run_time = time.time() - start_time
    print("Running time:", run_time, "s")
    print("----------------------------->ok \\(^~^)/ ok<-----------------------------")


if __name__ == "__main__":
    main()
