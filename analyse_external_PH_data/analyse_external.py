import json
from persim import plot_diagrams
import numpy as np


with open('ph_vanilla_cubic_diagram.txt') as json_f:
    d = json.load(json_f)
    print(np.array(d[0]).shape)
    plot_diagrams(np.array(d[0]), show=True)
    plot_diagrams(np.array(d[1]), show=True)
