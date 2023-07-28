import json
import numpy as np
from ripser import ripser


file = open('mini_example_data.csv', 'rb')
data = np.loadtxt(file, delimiter=",", dtype=np.float64)
print(data.shape, data.dtype)

ph = ripser(data, distance_matrix=True)['dgms']
ph_list = [homology.tolist() for homology in ph]

with open('ph_mini_example_diagram.txt', 'w') as f:
    json.dump(ph_list, f)
