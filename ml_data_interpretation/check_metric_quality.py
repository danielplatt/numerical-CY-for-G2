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


def load_data(foldername='data/CYmetric'):
    points = np.load(f'{foldername}/points_complex.npy')
    metric = np.load(f'{foldername}/cymetric.npy')
    losses = np.load(f'{foldername}/unweighted_loss.npy')
    unrestricted_metric = None # np.load('data/metrics_and_restrictions/cymetric_cpn.npy')
    restrictions = np.load(f'{foldername}/restriction.npy')
    return points, metric, losses, unrestricted_metric, restrictions


def get_approx_distances_to_near_singularity(points):
    singpoints = load_sing_points()
    distances = []
    for n, point in enumerate(points):
        if n % 10 == 0:
            print(f'Computing distance for point number {n}')
        distances += [get_distance_to_sing_locus(point, singpoints)]
    return distances


def get_approx_distances_to_near_singularity_old(points, which_variety='quintic'):
    distances = []
    for p in points:
        p = p/np.linalg.norm(p)
        if which_variety=='quintic':
            quadric = abs(p[0]**2+p[1]**2+p[2]**2+p[3]**2+p[4]**2)
            cubic = abs(p[0]*(p[1]**2+p[2]**2+p[3]**2-p[4]**2)-(p[1]**2+p[2]**2+p[3]**2-1/2*p[4]**2)-1/2*p[0]**3)
            distances += [max(quadric, cubic)]
        elif which_variety=='quadric_quartic':
            dist_to_zero_dim_sing = max(abs(p[4]**2+p[5]**2), abs(p[0]), abs(p[1]), abs(p[2]), abs(p[3]))
            dist_to_one_dim_sing = max(abs(p[1]**4+p[2]**4+p[3]**4), abs(p[0]), abs(p[4]), abs(p[5]))
            dist_to_any_sing = min(dist_to_one_dim_sing, dist_to_zero_dim_sing)
            distances += [dist_to_any_sing]
    return distances


def plot_distance_versus_loss(points, losses, N=100000, which_variety='quintic', output_filename_suffix=None):
    distances = get_approx_distances_to_near_singularity_old(points[:N], which_variety=which_variety)
    plt.scatter(distances, losses[:N], s=1, alpha=0.5)
    if output_filename_suffix is not None:
        plt.savefig(f'output_images/distance_versus_loss_{output_filename_suffix}.png')
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

    dimension = len(point)

    hermitian_metric = np.zeros((dimension,dimension), dtype=np.complex64)
    for i in range(dimension):
        for j in range(dimension):
            if i==j:
                hermitian_metric[i,j] = (1+np.linalg.norm(point, ord=2)**2-np.abs(point[i]))/(1+np.linalg.norm(point, ord=2)**2)**2
            else:
                hermitian_metric[i, j] = (- point[j] * np.conj(point[i]))/(1+np.linalg.norm(point, ord=2)**2)**2
    # print(f'hermitian_metric={hermitian_metric}')
    return hermitian_metric
    #(hermitian_metric+np.conj(hermitian_metric))/2


def get_real_restricted_fubini_study(point, restrict_by_poly_equations):
    affine_fubini = get_real_fubini_study_metric(get_affine_point_coordinates(point))
    restricted_fubini = np.matmul(restrict_by_poly_equations, np.matmul(affine_fubini, np.conj(np.transpose(restrict_by_poly_equations))))
    return (restricted_fubini+np.conj(restricted_fubini))/2


def get_stretch_factors(point, metric, restriction):
    """Compare the Fubini-Study metric with the learned approximate
    Calabi-Yau metric in a point.

    :param point: Point with any number of coordinates, one of which must be equal to 1
    :param metric: metric (can be 3x3 i.e. already restricted, or unrestricted, i.e. of dimension of ambient vector space, i.e. 6x6 on CP^5)
    :param restriction: 5x3 restriction map
    :return:
    """
    if not metric.shape==(3,3):
        # print(f'metric before restrict: {metric}')
        # print(np.linalg.eigvals((metric+np.conj(metric))/2))
        # need to restrict metric
        metric=np.matmul(restriction, np.matmul(metric, np.conj(np.transpose(restriction))))
        # print(np.linalg.eigvals((metric + np.conj(metric)) / 2))
    fubi_metric = get_real_restricted_fubini_study(point, get_four_coord_restriction(restriction, get_affine_patch_index(point)))
    fubi_root = matrix_root(fubi_metric)

    real_metric = get_real_part_of_metric(metric)

    stretch_factors = np.linalg.eigvals(np.matmul(np.linalg.inv(fubi_root), real_metric, np.linalg.inv(fubi_root)))
    # if stretch_factors[-1] < 0:
    #     print(f'point: {point}')
    #     print(f'fubi_metric: {fubi_metric}')
    #     print(f'fubi_root: {fubi_root}')
    #     print(f'fubi_root**2: {np.matmul(fubi_root, fubi_root)}')
    #     print(f'real_metric: {real_metric}. ev: {np.linalg.eigvals(real_metric)}')
    #     print(f'Hermitian metric: {metric}')
    #     print(f'stretch_factors: {stretch_factors}')
    #     input()
    return stretch_factors


def plot_distance_to_sing_versus_neckness(points, metrics, restrictions, which_variety='quadric_quartic', N=100000, distance_computation='crude', output_filename_suffix=None):
    if distance_computation=='crude':
        distances = get_approx_distances_to_near_singularity_old(points[:N], which_variety=which_variety)
    else:
        distances = get_approx_distances_to_near_singularity(points[:N])

    biggest_evs = []
    smallest_evs = []
    for point, metric, restriction in zip(points[:N], metrics[:N], restrictions[:N]):
        eigenvalues = get_stretch_factors(point, metric, restriction)
        biggest_evs += [eigenvalues[2]]
        smallest_evs += [eigenvalues[0]]

    plt.scatter(distances, biggest_evs[:N], s=1, alpha=0.5)
    if output_filename_suffix is not None:
        plt.savefig(f'output_images/distance_versus_neckness_{output_filename_suffix}.png')
    plt.show()
    #plt.scatter(distances, smallest_evs[:N], s=1, alpha=0.5)
    #plt.show()


def main():
    points, metric, losses, _, restrictions = load_data(foldername='CICY/cicy1/CYmetric') #'CICY/cicy2/CYmetric'   ////   'data/CYmetric'

    # plot_distance_versus_loss(points, losses, N=100000, which_variety='quadric_quartic', output_filename_suffix='cicy2')

    plot_distance_to_sing_versus_neckness(points, metric, restrictions, which_variety='quadric_quartic', N=10000, distance_computation='crude', output_filename_suffix='cicy1')

    points, metric, losses, _, restrictions = load_data(foldername='CICY/cicy2/CYmetric')
    plot_distance_to_sing_versus_neckness(points, metric, restrictions, which_variety='quadric_quartic', N=10000,
                                          distance_computation='crude', output_filename_suffix='cicy2')



if __name__ == '__main__':
    main()
