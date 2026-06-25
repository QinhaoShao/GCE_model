import numpy as np
import scipy.integrate as spi
import scipy.special as sps
import pandas as pd

def migration(sigma, Rb):
    #Compute the radial migration probability matrix from birth radius Rb to the present-day radial bins.
    #sigma[b, i] controls the migration width from the b-th birth radius to the i-th present-day radial interval.
    pi = np.pi
    age = np.append(np.append(np.arange(10) * 0.02 + 0.01, np.arange(8) * 0.1 + 0.25), np.arange(13) * 1 + 1.5) #Age-bin centers used for the migration probabilities, in Gyr.
    N_bins = len(Rb)
    age_num = len(age)
    P = np.zeros((N_bins, age_num, N_bins)) #P[birth radial bin, age bin, migrated radial bin].
    intervals = [(1, 3), (3, 5), (5, 7), (7, 9), (9, 11), (11, 13), (13, 15)]
    Rd_data = pd.read_csv('../Rdisc_and_Reff_birth.csv')
    Rd_t = np.array(Rd_data['t_lb'])
    Rd_R = np.array(Rd_data['R_disc_birth_50'])
    
    def P_R_R0pm_func(r, b, t, delta):
        #Return the probability density at present-day radius r for a given birth radius, age, and migration strength.
        idv = np.argmin(abs(Rd_t - age[t]))
        Rd = Rd_R[idv]
        
        R0 = Rb[b] - (delta**2 *age[t])/ (16 * Rd) #Drift term of the migration kernel center after accounting for the disk scale length.
        z = R0 / (delta * np.sqrt(2 * age[t] / 8))
        c = delta * np.sqrt(pi / 2 * age[t] / 8) * (sps.erf(z) + 1)
        return c**(-1) * np.exp(-(r - R0)**2 / (2 * delta**2 * age[t] / 8))
    
    for b in range(N_bins):
        for t in range(age_num):
            for i, (n, m) in enumerate(intervals):
                #Integrate over each present-day radial interval to obtain the discrete migration probability.
                delta = sigma[b,i]
                integral_value, _ = spi.quad(P_R_R0pm_func, n, m, args=(b, t, delta))
                P[b, t, i] = integral_value
                
    return P
