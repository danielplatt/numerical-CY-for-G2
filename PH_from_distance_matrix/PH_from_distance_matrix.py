import json
import numpy as np
from gph import ripser_parallel
import time
from tqdm import tqdm
import psutil

from gtda.plotting import plot_diagram
from gtda.homology._utils import _postprocess_diagrams

from ripser import ripser


def load_as_np():
    file = open('vanilla_modified_distances-generic-deformation.csv', 'rb')
    data = np.loadtxt(file, delimiter=",", dtype=np.float64)
    print(data.shape, data.dtype)
    return data

def explore_data(data):
    print(f'Number of points: {len(data)}')

    closest_point_distances = []
    average_point_distances = []
    check_how_many = 50
    for k in tqdm(range(check_how_many)):
        point_k_distances = []
        for l in range(len(data)):
            if k != l:
                point_k_distances += [np.linalg.norm(data[k] - data[l])]
        closest_point_distances += [np.min(point_k_distances)]
        average_point_distances += [np.average(point_k_distances)]

    print(f'Average distance to closest point: {np.average(closest_point_distances)}')
    print(f'Min distance to closest point: {np.min(closest_point_distances)}')
    print(f'Max distance to closest point: {np.max(closest_point_distances)}')
    print(f'Average distance between two points: {np.average(average_point_distances)}')

def gph_run(data, thresh):
    dgm = ripser_parallel(data, metric="precomputed", thresh=thresh, maxdim=2, n_threads=1)

    dgm_gtda = _postprocess_diagrams([dgm["dgms"]], "ripser", (0, 1), np.inf, True)[0]
    plot_diagram(dgm_gtda, homology_dimensions=(0, 1)).show()


def ripser_run(data, thresh):
    ph = ripser(data, thresh=thresh, distance_matrix=True, maxdim=1)['dgms']
    ph_list = [homology.tolist() for homology in ph]

    with open('ph_vanilla_cubic_diagram.txt', 'w') as f:
        json.dump(ph_list, f)


if __name__ == '__main__':
    print(f'Memory usage in MB: {psutil.Process().memory_info().rss/1024 ** 2} (normally around 116.671875)')
    start_time = time.time()
    data = load_as_np()
    load_time = time.time() - start_time
    print(f'Memory usage in MB: {psutil.Process().memory_info().rss / 1024 ** 2} (normally around 10207.140625 for float64)')
    print(f'Load time: {load_time}')
    thresh = 2*0.11757912072752902-0.0000000001
    ripser_run(data, thresh)
    print(f'Memory usage in MB: {psutil.Process().memory_info().rss / 1024 ** 2}')
    print(f'Threshold: {thresh}. Computation time: {time.time() - start_time - load_time}')
    print(f'Total time: {time.time() - start_time}')
