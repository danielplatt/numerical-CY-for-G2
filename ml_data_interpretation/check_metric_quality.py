import numpy as np
from matplotlib import pyplot as plt


def load_data():
    points = np.load('data/CYmetric/points_complex.npy')
    metric = np.load('data/CYmetric/cymetric.npy')
    losses = np.load('data/CYmetric/unweighted_loss.npy')
    return points, metric, losses


def get_approx_distances_to_near_singularity(points):
    distances = []
    for p in points:
        p = p/np.linalg.norm(p)
        quadric = abs(p[0]**2+p[1]**2+p[2]**2+p[3]**2+p[4]**2)
        cubic = abs(p[0]*(p[1]**2+p[2]**2+p[3]**2-p[4]**2)-(p[1]**2+p[2]**2+p[3]**2-1/2*p[4]**2)-1/2*p[0]**3)
        distances += [max(quadric, cubic)]
    return distances


def plot_distance_versus_loss(points, losses, N=100000):
    distances = get_approx_distances_to_near_singularity(points[:N])
    plt.scatter(distances, losses[:N], s=1, alpha=0.5)
    plt.show()


def plot_distance_to_sing_versus_neckness(points, metric, N=100000):
    distances = get_approx_distances_to_near_singularity(points[:N])

    biggest_evs = []
    smallest_evs = []
    for point_metric in metric[:N]:
        eigenvalues = np.linalg.eigh(point_metric)[0]
        biggest_evs += [eigenvalues[2]]
        smallest_evs += [eigenvalues[0]]

    plt.scatter(distances, biggest_evs[:N], s=1, alpha=0.5)
    plt.show()


def main():
    points, metric, losses = load_data()
    plot_distance_versus_loss(points, losses, N=100000)
    plot_distance_to_sing_versus_neckness(points, metric, N=100000)


if __name__ == '__main__':
    main()
