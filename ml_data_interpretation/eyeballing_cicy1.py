import numpy as np
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn_extra.cluster import KMedoids
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.decomposition import PCA
import sys


def load_data(points_path='data/One_form/points_real.npy', omega_path='data/One_form/omega.npy'):
    points_real = np.load(points_path)
    omega_norm = np.load(omega_path)
    print(f'{points_path} points_real.shape: {points_real.shape}')
    print(f'{omega_path} omega.shape: {omega_norm.shape}')
    return points_real, omega_norm


def main():
    points_path = 'data/cicy1/one_form/points_real.npy'
    omega_path = 'data/cicy1/one_form/omega.npy'
    points_real, omega = load_data(points_path=points_path, omega_path=omega_path)
    print(omega)


if __name__ == '__main__':
    main()
