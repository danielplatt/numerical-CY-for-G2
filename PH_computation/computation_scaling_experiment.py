import numpy as np
from gph import ripser_parallel
import time
from itertools import combinations
from tqdm import tqdm


def load_as_np():
    file = open('pointcloud-80k-pts.csv', 'rb')
    data = np.loadtxt(file, delimiter=",")
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

def main(data, thresh):
    start_time = time.time()
    dgm = ripser_parallel(data, thresh=thresh, maxdim=2, n_threads=-1)

    print(dgm)
    print(f'Threshold: {thresh}. Computation time: {time.time() - start_time}')


if __name__ == '__main__':
    data = load_as_np()
    # explore_data(data)
    main(data, 0.2)

