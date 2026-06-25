from pylab import *
from astropy.io import ascii
from migration_list import *
import random

def migration_model_data(file_path, R, file, sigma, Rb):
    #Read the mock stellar catalog for a given birth radius, group the stars by age, and apply radial migration probabilities.
    model_data = ascii.read(file_path)
    age = model_data['age']
    p_m = migration(sigma, Rb) #Radial migration probability matrix.
    N_bins = len(Rb)

    age_tid = np.append(np.append(np.arange(10)*0.02, np.arange(8)*0.1+0.2), np.arange(14)*1+1)
    age_id = np.digitize(age, age_tid) - 1

    data_grouped = {}

    seven_grouped = {f"R_{i*2+2}": [] for i in range(N_bins)}

    for i, data in enumerate(model_data):
        current_age = age_id[i]
        if current_age not in data_grouped:
            data_grouped[current_age] = []
        data_grouped[current_age].append(data)

    for age_id, data_list in data_grouped.items():
        num = len(data_list)

        distrib = p_m[R, age_id, :] * num
        distrib_n = np.zeros(N_bins, dtype=int)
        for i in range(N_bins):
            distrib_n[i] = int(floor(distrib[i]))
            if random.random() < distrib[i] - distrib_n[i]:
                distrib_n[i] += 1
            samples = random.sample(data_list, distrib_n[i])
            seven_grouped[f"R_{i*2+2}"].extend(samples)
    return seven_grouped

def mock_migration(file_path, filename, sigma, Rb):
    #Combine the mock catalogs from all birth radii and write the migrated present-day radial catalogs R_*.txt.
    all_data = {}
    for i, file in enumerate(filename):
        path = file_path + 'model_' + file + '.txt'
        migration_data = migration_model_data(path, i, file, sigma, Rb)
        for key, value in migration_data.items():
            if key not in all_data:
                all_data[key] = []
            all_data[key].extend(value)

    for key, value in all_data.items():
        file_name = file_path + key + '.txt'
        with open(file_name, 'w') as file:
            file.write("age\tmass\tfe/h\tmg/fe\tsi/fe\tage_0\n")
            for sample in value:
                file.write(f"{sample['age']}\t{sample['mass']}\t{sample['fe/h']}\t{sample['mg/fe']}\t{sample['si/fe']}\t{sample['age_0']}\n")
