import numpy as np

def get_num(data_file):
    #Estimate the number of mock stars to sample in each radial annulus, following Lian et al. (2025).

    data = np.loadtxt(data_file, skiprows=1)
    total_log = data[:, -1]
    sigma_lin = 10 ** total_log

    num_new = []
    for i in range(7):
        idx1 = 2 * i + 1
        idx2 = 2 * i + 2
        avg = (sigma_lin[idx1] + sigma_lin[idx2]) / 2.0 * 1e3
        num_new.append(int(round(avg)))
    
    return num_new
