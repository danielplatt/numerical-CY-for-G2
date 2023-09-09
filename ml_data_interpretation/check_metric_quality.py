import numpy as np
from matplotlib import pyplot as plt
from scipy.linalg import sqrtm


def get_affine_patch_index(points):
    one_entries = np.where(np.abs(points-1) < 0.00000000001)[0]
    assert len(one_entries) == 1

    return one_entries[0]


def get_affine_point_coordinates(point):
    one_entries = np.where(np.abs(point-1) < 0.00000000001)[0]
    assert len(one_entries) == 1
    return np.delete(point, one_entries[0])


def get_four_coord_restriction(five_coord_restriction, affine_patch_index):
    return np.delete(five_coord_restriction, affine_patch_index, axis=1)


def get_real_part_of_metric(matrix):
    return (matrix+np.conj(matrix))/2


def matrix_root(mat):
    return sqrtm(mat)

    eigenpairs = np.linalg.eig(mat)
    eigenvals = eigenpairs[0]
    eigentransformation = np.linalg.inv(eigenpairs[1])
    diagroot = np.diag([np.sqrt(val) for val in eigenvals])
    this_root = np.linalg.inv(eigentransformation).dot(diagroot.dot(eigentransformation))
    return this_root


def load_data():
    points = np.load('data/CYmetric/points_complex.npy')
    metric = np.load('data/CYmetric/cymetric.npy')
    losses = np.load('data/CYmetric/unweighted_loss.npy')
    unrestricted_metric = np.load('data/metrics_and_restrictions/cymetric_cpn.npy')
    restrictions = np.load('data/metrics_and_restrictions/restriction.npy')
    return points, metric, losses, unrestricted_metric, restrictions


def get_approx_distances_to_near_singularity(points):
    singpoints = load_sing_points()
    distances = []
    for n, point in enumerate(points):
        if n % 10 == 0:
            print(f'Computing distance for point number {n}')
        distances += [get_distance_to_sing_locus(point, singpoints)]
    return distances


def get_approx_distances_to_near_singularity_old(points):
    distances = []
    for p in points:
        p = p/np.linalg.norm(p)
        np.linalg
        quadric = abs(p[0]**2+p[1]**2+p[2]**2+p[3]**2+p[4]**2)
        cubic = abs(p[0]*(p[1]**2+p[2]**2+p[3]**2-p[4]**2)-(p[1]**2+p[2]**2+p[3]**2-1/2*p[4]**2)-1/2*p[0]**3)
        distances += [max(quadric, cubic)]
    return distances


def plot_distance_versus_loss(points, losses, N=100000):
    distances = get_approx_distances_to_near_singularity_old(points[:N])
    plt.scatter(distances, losses[:N], s=1, alpha=0.5)
    plt.savefig('output_images/distance_versus_loss.png')
    plt.show()


def load_sing_points():
    return np.genfromtxt('extra_data/singularQuinticPoints.csv', dtype=complex, delimiter=',')


def fubini_study_distance(p, q):
    return np.arccos(np.sqrt((np.abs(np.inner(p, q)))**2/(np.linalg.norm(p)**2*np.linalg.norm(q)**2)))


def get_distance_to_sing_locus(point, singpoints):
    min_distance = 1000
    for singpoint in singpoints:
        min_distance = min(min_distance, fubini_study_distance(singpoint, point))
    return min_distance


def get_real_fubini_study_metric(point):
    '''For points in affine coordinates (i.e. have thrown away the =1 coordinate. For CP^4 this should have length 4.)'''

    hermitian_metric = np.zeros((4,4), dtype=np.complex64)
    for i in [0,1,2,3]:
        for j in [0,1,2,3]:
            if i==j:
                hermitian_metric[i,j] = (1+np.linalg.norm(point, ord=2)**2-np.abs(point[i]))/(1+np.linalg.norm(point, ord=2)**2)**2
            else:
                hermitian_metric[i, j] = (- point[j] * np.conj(point[i]))/(1+np.linalg.norm(point, ord=2)**2)**2
    # print(f'hermitian_metric={hermitian_metric}')
    return hermitian_metric
    #(hermitian_metric+np.conj(hermitian_metric))/2


def get_real_restricted_fubini_study(point, four_to_three_restriction):
    fubini_study_four_by_four = get_real_fubini_study_metric(get_affine_point_coordinates(point))
    restricted_fubini = np.matmul(four_to_three_restriction, np.matmul(fubini_study_four_by_four, np.conj(np.transpose(four_to_three_restriction))))
    return (restricted_fubini+np.conj(restricted_fubini))/2


def get_stretch_factors(point, metric, restriction):
    """Compare the Fubini-Study metric with the learned approximate
    Calabi-Yau metric in a point.

    :param point: Point with 5 coordinates, one of which must be equal to 1
    :param metric: 3x3 metric
    :param restriction: 5x3 restriction map
    :return:
    """
    fubi_metric = get_real_restricted_fubini_study(point, get_four_coord_restriction(restriction, get_affine_patch_index(point)))
    fubi_root = matrix_root(fubi_metric)

    real_metric = get_real_part_of_metric(metric)

    stretch_factors = np.linalg.eigvals(np.matmul(np.linalg.inv(fubi_root), real_metric, np.linalg.inv(fubi_root)))
    # if stretch_factors[-1] < 0:
        # print(f'point: {point}')
        # print(f'fubi_metric: {fubi_metric}')
        # print(f'fubi_root: {fubi_root}')
        # print(f'fubi_root**2: {np.matmul(fubi_root, fubi_root)}')
        # print(f'real_metric: {real_metric}')
        # print(f'stretch_factors: {stretch_factors}')
        # input()
    return stretch_factors


def plot_distance_to_sing_versus_neckness(points, metrics, restrictions, N=100000, distance_computation='crude'):
    if distance_computation=='crude':
        distances = get_approx_distances_to_near_singularity_old(points[:N])
    else:
        distances = get_approx_distances_to_near_singularity(points[:N])

    biggest_evs = []
    smallest_evs = []
    for point, metric, restriction in zip(points[:N], metrics[:N], restrictions[:N]):
        eigenvalues = get_stretch_factors(point, metric, restriction)
        biggest_evs += [eigenvalues[2]]
        smallest_evs += [eigenvalues[0]]

    plt.scatter(distances, biggest_evs[:N], s=1, alpha=0.5)
    #plt.ylim([0, 1000])
    plt.savefig('output_images/distance_versus_neckness.png')
    plt.show()
    #plt.scatter(distances, smallest_evs[:N], s=1, alpha=0.5)
    #plt.show()


def main():
    points, metric, losses, unrestricted_metric, restrictions = load_data()

    plot_distance_versus_loss(points, losses, N=100000)
    plot_distance_to_sing_versus_neckness(points, metric, restrictions, N=100000, distance_computation='crude')




if __name__ == '__main__':
    main()
