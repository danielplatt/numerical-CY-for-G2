import numpy as np
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn_extra.cluster import KMedoids
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.decomposition import PCA
import sys

np.set_printoptions(threshold=sys.maxsize)


def load_data(points_path='data/One_form/points_real.npy', omega_norm_path='data/One_form/norm_square.npy'):
    points_real = np.load(points_path)
    omega_norm = np.load(omega_norm_path)
    print(f'{points_path} points_real.shape: {points_real.shape}')
    print(f'{omega_norm_path} omega_norm.shape: {omega_norm.shape}')
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
    return np.arccos(np.sqrt(np.abs(np.inner(p, q))**2/(np.abs(np.inner(p, p)) *np.abs(np.inner(q, q)))))


def euclidean_distance(p, q):
    return np.linalg.norm(p-q)


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


def get_normalized_small_points(points_real, omega_norm, return_all_points='False', num_points=500):

    indexes = sorted(range(len(omega_norm)), key=omega_norm.__getitem__)
    omega_norm_sorted = np.array(list(map(omega_norm.__getitem__, indexes)))
    points_real_sorted = np.array(list(map(points_real.__getitem__, indexes)))

    # ind = np.where(omega_norm<1e-6)
    # return points_real[ind], None

    if return_all_points:
        return normalize_points_by_first_coordinate(points_real_sorted[:num_points]), normalize_points_by_first_coordinate(points_real_sorted[num_points:])
    else:
        return normalize_points_by_first_coordinate(points_real_sorted[:num_points]), omega_norm_sorted[:num_points]

    # small_points = []
    # small_norms = []
    # omega_L1_norm = sum(omega_norm)/len(omega_norm)
    # for point, norm in zip(points_real, omega_norm):
    #     if norm <= 0.001*omega_L1_norm:
    #         small_points += [point]
    #         small_norms += [norm]
    #
    # small_points = np.array(small_points)
    # small_points = normalize_points_by_first_coordinate(small_points)
    #
    # return small_points, small_norms


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


def find_medoids_cluster(points_real, omega_norm, label='', silhouette=False, num_points='500'):
    small_points, _ = get_normalized_small_points(points_real, omega_norm, num_points=num_points)
    np.set_printoptions(threshold=np.inf)
    # print(small_points)

    print(f'small_points.shape: {small_points.shape}')

    k_values = []
    scores = []
    silhouette_scores = []

    cluster_centers = None

    for k in range(100):
        if k <= 0:
            continue

        print(f'Now clustering for k={k}')
        kmedoids = KMedoids(n_clusters=k, random_state=0, metric=fubini_study_distance).fit(small_points)
        if silhouette:
            cluster_labels = KMedoids(n_clusters=k, random_state=0, metric=fubini_study_distance).fit_predict(small_points)

        k_values += [k]
        scores += [kmedoids.inertia_]
        if k>=2 and silhouette:
            silhouette_scores += [silhouette_score(small_points, cluster_labels)]

        if k == 4:
            cluster_centers = kmedoids.cluster_centers_

    print('list:')
    print(list(zip(k_values, scores)))
    plt.plot(k_values, scores)
    plt.title(f'Medoids scores for {label} (num_points={len(small_points)})')
    plt.savefig(f'output_images/medoids_scores_{label}_{len(small_points)}.png')
    plt.show()
    if silhouette:
        plt.title(f'Average silhouette scores for {label} (num_points={len(small_points)})')
        plt.plot(k_values[1:], silhouette_scores)
        plt.savefig(f'output_images/averagesilhouette_scores_{label}_{len(small_points)}.png')
        plt.show()

    # print('k=4 Cluster centers:')
    # print(cluster_centers)

def plot_points(points_real, omega_norm, label='', clustering='full_points', num_points=100, distance='fubini'):
    small_points, large_points = get_normalized_small_points(points_real, omega_norm, return_all_points='True', num_points=num_points)
    print(f'np.mean(small_points, axis=0): {np.mean(small_points, axis=0)}')
    print(f'np.std(small_points, axis=0): {np.std(small_points, axis=0)}')

    # small_points = np.multiply(small_points, np.array([1,1, 1, 1, 1/5, 1/5]))

    if clustering == 'full_points':
        print('Start clustering...')
        if distance=='fubini':
            kmedoids = KMedoids(n_clusters=2, random_state=0, metric=fubini_study_distance, method='pam').fit(small_points)
        else:
            # use Euclidean distance
            kmedoids = KMedoids(n_clusters=2, random_state=0, metric=euclidean_distance, method='pam').fit(small_points)

        point_labels = kmedoids.fit_predict(small_points)
        print(small_points)
        print('Clustering ended.')

    # 2D plot
    print('Start PCA...')
    pca = PCA(n_components=2)
    pca.fit(small_points)
    print(f'mean: {pca.mean_}')
    print(f'components: {pca.components_}')
    projected_points = pca.transform(small_points)
    print('PCA ended.')

    if clustering == 'projected_points':
        if distance=='fubini':
            kmedoids = KMedoids(n_clusters=2, random_state=0, metric=fubini_study_distance, method='pam').fit(projected_points)
        else:
            # use Euclidean distance
            kmedoids = KMedoids(n_clusters=2, random_state=0, metric=euclidean_distance, method='pam').fit(projected_points)
        point_labels = kmedoids.fit_predict(small_points)
        print(small_points)
        print('Clustering ended.')

    plt.title(f'Points with small norm for {label} (number points={num_points})\n distance={distance}')
    if clustering:
        plt.scatter(projected_points[point_labels == 0][:,0], projected_points[point_labels == 0][:,1], color='red')
        plt.scatter(projected_points[point_labels == 1][:, 0], projected_points[point_labels == 1][:, 1], color='blue')
        plt.scatter(projected_points[point_labels == 2][:,0], projected_points[point_labels == 2][:,1], color='green')
        plt.scatter(projected_points[point_labels == 3][:, 0], projected_points[point_labels == 3][:, 1], color='black')
    else:
        print(projected_points.shape)
        plt.scatter(projected_points[:][:,0], projected_points[:][:,1], color='red')
    plt.show()

    #3D plot
    # PCA
    pca = PCA(n_components=3)
    pca.fit(small_points)
    print(f'mean: {pca.mean_}')
    print(f'components: {pca.components_}')
    projected_points = pca.transform(small_points)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    if clustering == 'projected_points':
        if distance == 'fubini':
            kmedoids = KMedoids(n_clusters=2, random_state=0, metric=fubini_study_distance, method='pam').fit(projected_points)
        else:
            # use Euclidean distance
            kmedoids = KMedoids(n_clusters=2, random_state=0, metric=euclidean_distance, method='pam').fit(projected_points)
        point_labels = kmedoids.fit_predict(small_points)
        print(small_points)
        print('Clustering ended.')


    if clustering:
        ax.scatter(projected_points[point_labels == 0][:, 0], projected_points[point_labels == 0][:, 1], projected_points[point_labels == 0][:, 2], marker='o')
        ax.scatter(projected_points[point_labels == 1][:, 0], projected_points[point_labels == 1][:, 1], projected_points[point_labels == 1][:, 2], marker='^')
    else:
        ax.scatter(projected_points[:][:, 0], projected_points[:][:, 1], projected_points[:][:, 2], marker='o')
    plt.title(f'Points with small norm for {label} (number points={num_points})\n distance={distance}')
    plt.show()

    # hardcoded projection vectors
    if label == 'cicy1' or label == 'cicy2':
        z2 = np.array([0,0,1,0,0,0])
        z4 = np.array([0,0,0,0,1,0])
        z5 = np.array([0, 0, 0, 0, 0, 1])
        print(small_points.shape, z2.shape)
        component1 = np.inner(small_points, z2)
        component2 = np.inner(small_points, z4)
        component3 = np.inner(small_points, z5)
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(component1, component2, component3, marker='o')
        plt.title(f'Points with small norm for {label} (number points={num_points}) [proj onto z2,z4,z5]\n distance={distance}')
        plt.show()


def points_statistical_data(points_real, omega_norm, label=''):
    points_normalised = normalize_points_by_first_coordinate(points_real)
    print(points_normalised)
    print(points_normalised.shape)
    print(np.mean(points_normalised, axis=0))
    print(np.std(points_normalised, axis=0))

def eyeballing_cicy1(points_real, omega_norm):
    zIcoords = points_real[:,2]
    indices = sorted(range(len(points_real)), key=zIcoords.__getitem__)
    for i in indices:
        print(f'zI={zIcoords[i]} -> omega_norm={omega_norm[i]}')
    plt.scatter(zIcoords, omega_norm)
    plt.show()



if __name__ == '__main__':
    load_paths = [
        # ['data/One_form/points_real.npy', 'data/One_form/norm_square.npy'],
        ['data/quintic2/one_form/points_real.npy', 'data/quintic2/one_form/norm_square.npy'],
        # ['data/cicy1/one_form/points_real.npy', 'data/cicy1/one_form/norm_square.npy'],
        # ['data/cicy2/one_form/points_real.npy', 'data/cicy2/one_form/norm_square.npy'], # no small points here...
    ]
    labels = [
        # 'oldquintic',
        'newquintic',
        # 'cicy1',
        # 'cicy2'
    ]
    for label, load_path_pair in zip(labels, load_paths):
        points_real, omega_norm = load_data(points_path=load_path_pair[0], omega_norm_path=load_path_pair[1])

        random_indices = np.random.choice(points_real.shape[0], size=2000, replace=False)
        # eyeballing_cicy1(points_real[:], omega_norm[:])
        # points_statistical_data(points_real, omega_norm, label=label)
        # plot_points(points_real, omega_norm, label=label, clustering='full_points', num_points=100, distance='fubini')
        find_medoids_cluster(points_real, omega_norm, label=label, num_points=20000)
