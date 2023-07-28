import numpy as np


def main():
    np.random.seed(10)
    random_matrix = np.random.rand(20,20)
    random_distance_matrix = random_matrix + np.transpose(random_matrix)
    for i in range(20):
        random_distance_matrix[i,i] = 0
    print(random_distance_matrix)
    np.savetxt('mini_example_data.csv', random_distance_matrix, delimiter=',')


if __name__ == '__main__':
    main()
