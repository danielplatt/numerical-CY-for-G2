import numpy as np
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn_extra.cluster import KMedoids


def load_data():
    points_real = np.load('data/One_form/points_real.npy')
    omega_norm = np.load('data/One_form/norm_square.npy')
    return points_real, omega_norm


def show_histograms(points_real, omega_norm):
    max_norm = max(omega_norm)
    min_norm = min(omega_norm)
    avg_norm = np.average(omega_norm)

    plotmax = max_norm
    for i in range(5):
        omega_selected = [element for element in omega_norm if element <= plotmax]
        fig, ax = plt.subplots(figsize=(10, 7))
        ax.hist(omega_selected, bins=100)
        plt.show()
        plotmax /= 10


def fubini_study_distance(p, q):
    return np.abs(np.inner(p, q))**2/(np.abs(np.inner(p, p))**2 *np.abs(np.inner(q, q))**2)


def fubini_study_distance_tuple(tup):
    return fubini_study_distance(tup[0], tup[1])


def average_nearest_point(points_real, dist, howmany):
    nearest_distances = []
    for i, pointi in enumerate(points_real[:howmany]):
        print(f'Finding nearest point for point {i}')
        nearest_dist = np.inf
        for k, pointk in enumerate(points_real[:howmany]):
            if i==k:
                continue
            nearest_dist = min(nearest_dist, dist(pointi, pointk))
        nearest_distances += [nearest_dist]
    return nearest_distances


def normalize_points_by_first_coordinate(points):
    return np.transpose(np.transpose(points)/np.transpose(points[:,0]))


def get_normalized_small_points(points_real, omega_norm):
    small_points = []
    small_norms = []
    for point, norm in zip(points_real, omega_norm):
        if norm <= 2.0e-7:
            small_points += [point]
            small_norms += [norm]

    small_points = np.array(small_points)
    small_points = normalize_points_by_first_coordinate(small_points)

    return small_points, small_norms


def find_agglomerate_cluster(points_real, omega_norm):
    # this is meaningless because I made a particular, not motivated choice for distance_threshold
    small_points, _ = get_normalized_small_points(points_real, omega_norm)

    cutoff = 1.9899897075452778e-07 # np.average(average_nearest_point(small_points, fubini_study_distance, len(small_points)))

    distance_matrix = np.zeros((len(small_points), len(small_points)))
    for i, pointi in enumerate(small_points):
        for k, pointk in enumerate(small_points):
            distance_matrix[i,k] = fubini_study_distance(pointi, pointk)

    ACluster = AgglomerativeClustering(metric='precomputed', n_clusters=None, linkage='complete', distance_threshold=445000*cutoff)
    AClusterResult = ACluster.fit(distance_matrix)
    print(AClusterResult.n_clusters_)


def find_medoids_cluster(points_real, omega_norm):
    small_points, _ = get_normalized_small_points(points_real, omega_norm)

    k_values = []
    scores = []

    cluster_centers = None

    for k in range(10):
        if k == 0:
            continue

        print(f'Now clustering for k={k}')
        kmedoids = KMedoids(n_clusters=k, random_state=0, metric=fubini_study_distance).fit(small_points)

        k_values += [k]
        scores += [kmedoids.inertia_]

        if k == 4:
            cluster_centers = kmedoids.cluster_centers_

    plt.plot(k_values, scores)
    plt.savefig('output_images/medoids_scores.png')
    plt.show()

    print('Cluster centers:')
    print(cluster_centers)


if __name__ == '__main__':
    points_real, omega_norm = load_data()

    #avg_nearest_distances = average_nearest_point(points_real, fubini_study_distance, 1000)
    #print(avg_nearest_distances)
    #print(np.average(avg_nearest_distances), np.std(avg_nearest_distances))

    #find_agglomerate_cluster(points_real, omega_norm)

    find_medoids_cluster(points_real, omega_norm)
