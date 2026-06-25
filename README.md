GCE: Galactic Chemical Evolution Model
=========================================================

This project contains a multi-zone Galactic Chemical Evolution (GCE) model for
the Milky Way disk. The model tracks gas mass, star formation, elemental
abundances, and the chemical enrichment of the interstellar medium by different
nucleosynthesis channels across multiple radial annuli.

The current version is updated from the model framework of Lian et al. (2020).
Compared with the original model, this project extends the calculation to a
multi-zone structure and includes stellar radial migration based on
Frankel et al. (2018).

-----------------------------------------------------------------------------------------
Project Structure:

code/
  run_MW.py
    Main execution script. Most model parameters are set here. Running this
    file generates the chemical evolution tracks and mock stellar catalogs.

  galce.py
    Main class for the chemical evolution model. It reads yield tables, evolves
    gas mass and elemental masses, calculates event rates, and outputs FITS
    files for each radial annulus.

  accretion.py
    Builds the gas accretion rate history and star formation efficiency history.

  imf.py
    Implements the segmented Kroupa initial mass function, used to convert the
    total mass of newly formed stars into the number of stars on different mass
    grids.

  mock.py
    Generates mock stellar samples from the chemical evolution model tracks.

  mock_migration.py
    Applies radial migration to the mock stellar samples.

  migration_list.py
    Calculates the radial migration probability matrix from birth radius to
    present-day radius.

  surface_density.py
    Estimates the number of mock stars to sample in each radial annulus from
    surface density information.

yields/
  Stores nucleosynthetic yield tables for AGB stars, core-collapse supernovae,
  Type Ia supernovae, neutron star mergers, and magneto-rotational supernovae.

model/
  Default output directory for model FITS files and mock stellar catalog text
  files.

Input data files in the project root:
  elements_mass.csv
  lifetime-cube-230226.fits
  lifetime-cube-230313-lowm.fits
  mw-sbp.txt
  Rdisc_and_Reff_birth.csv
  solar_G07.txt

-----------------------------------------------------------------------------------------
To run the GCE model, enter the following commands in the terminal:
   cd code
   python run_MW.py
