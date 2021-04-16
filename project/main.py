import argparse
from sklearn.datasets import make_blobs
import random
import struct
import numpy as np
import matplotlib.pyplot as plt

from prj_lib import prj_main

MAX_N = 500
MAX_K = 128

def generate_data(k, n, d):
    return make_blobs(n_samples=n, n_features=d, centers=k)

def random_n():
    return random.randint(int(MAX_N / 2), MAX_N)

def random_k(n):
    max_k = min(n, MAX_K)
    return random.randint(int(max_k / 2), max_k)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", help="number of clusters", type=int, required=True)
    parser.add_argument("-n", help="number of data points", type=int, required=True)
    parser.add_argument("--Random", help="random data generation", action="store_true")

    args = parser.parse_args()

    initial_k = k = args.k
    n = args.n
    r = args.Random
    d = random.choice([2,3])
    m = 300

    print("Starting")
    
    if r:
        n = random_n()
        initial_k = random_k(n)
        k = 0
        
    data, clusters = generate_data(initial_k, n, d)
    print(f"Generated data with paramters n={n}, k={initial_k}, d={d}")

    print(f"Clustering data with Normalized Spectral Clustering and K-means++")
    k, binary_kmpp_labels, kmpp_jaccard, binary_nsc_labels, nsc_jaccard = prj_main(k, n, d, m, data.flatten().tolist(), clusters.flatten().tolist())

    kmpp_labels = np.array(zip(*struct.iter_unpack("Q", binary_kmpp_labels)).__next__())
    nsc_labels = np.array(zip(*struct.iter_unpack("Q", binary_nsc_labels)).__next__())

    if d == 3:
        fig = plt.figure()
        add_3d_subplot(fig, 1, "Normalized Spectral Clustering", data, nsc_labels, nsc_jaccard)
        add_3d_subplot(fig, 2, "K-means++", data, kmpp_labels, kmpp_jaccard)

    elif d == 2:
        fig = plt.figure()
        add_2d_subplot(fig, 1, "Normalized Spectral Clustering", data, nsc_labels, nsc_jaccard)
        add_2d_subplot(fig, 2, "K-means++", data, kmpp_labels, kmpp_jaccard)

    desc = f"""
    Clustered data was generated with parameters n={n}, k={initial_k}, d={d}
    Both algorithms created K={k} clusters
    """
    plt.figtext(0.5, 0.01, desc, wrap=True, horizontalalignment='center', fontsize=12)
    
    fig.set_size_inches(7, 8, forward=True)   
    plt.savefig("clusters.pdf")
    

def add_2d_subplot(fig, index, title, data, labels, jaccard):
    ax = fig.add_subplot(1, 2, index)
    cmap = plt.cm.get_cmap("hsv")
    for cluster in set(labels):
        ax.scatter(data[:, 0][labels == cluster],
                   data[:, 1][labels == cluster],
                   color=cmap(cluster * 30))
    ax.set_title(title + "\nJaccard Measure=%.3f" % jaccard) 
    ax.set_xlabel('X')
    ax.set_aspect('equal')


def add_3d_subplot(fig, index, title, data, labels, jaccard):
    ax = fig.add_subplot(1, 2, index, projection='3d')
    cmap = plt.cm.get_cmap("hsv")
    for cluster in set(labels):
        ax.scatter(data[:, 0][labels == cluster],
                   data[:, 1][labels == cluster],
                   data[:, 2][labels == cluster],
                   color=cmap(cluster * 30))
    ax.set_title(title + "\nJaccard Measure=%.3f" % jaccard) 
    ax.set_xlabel('X')
    ax.set_ylabel('Y')


if __name__ == "__main__":
    main()
