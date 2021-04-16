# main.py
# Clustering comparison main entry point

import argparse
from sklearn.datasets import make_blobs
import random
import struct
import numpy as np
import matplotlib.pyplot as plt
import time

from clustering import nsc_and_kmpp

KMPP_MAX_ITER = 300

# Estimates of Maximum Capacity within 5 minutes
MAX_N = 550
MAX_K = 128


def generate_data(k, n, d):
    """
    Returns n random samples of d featues centered around k random centroids
    """
    return make_blobs(n_samples=n, n_features=d, centers=k)


def random_n():
    """
    Returns a reasonable number of samples for clustering
    """
    return random.randint(int(MAX_N / 2), MAX_N)


def random_k(n):
    """
    Returns a reasonable number of clusters for n samples
    """
    max_k = min(n, MAX_K)
    return random.randint(int(max_k / 2), max_k)


def output_clusters_pdf(data, n, d, nsc_labels, nsc_jaccard, kmpp_labels, kmpp_jaccard, desc):
    """
    Writes a PDF file with visualization and description of the given data to ./clusters.pdf
    """
    fig = plt.figure()

    if d == 3:
        add_3d_subplot(fig, 1, "Normalized Spectral Clustering", data, nsc_labels, nsc_jaccard)
        add_3d_subplot(fig, 2, "K-means++", data, kmpp_labels, kmpp_jaccard)

    elif d == 2:
        add_2d_subplot(fig, 1, "Normalized Spectral Clustering", data, nsc_labels, nsc_jaccard)
        add_2d_subplot(fig, 2, "K-means++", data, kmpp_labels, kmpp_jaccard)

    plt.figtext(0.5, 0.01, desc, wrap=True, horizontalalignment='center', fontsize=12)

    fig.set_size_inches(7, 8, forward=True)
    plt.savefig("clusters.pdf")


def add_2d_subplot(fig, index, title, data, labels, jaccard):
    """
    Adds a 2d subplot based on data and labels to fig
    """
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
    """
    Adds a 3d subplot based on data and labels to fig
    """
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


def valid_parameters(k, n, r):
    return 0 < k < n

def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        description="""Clustering Algorithms Comparison
This program generates random data based on the given parameters and clusters the data
with both K-means++ algorithm and Normalized Spectral Clustering

For completion within 5 minutes, the Maximum Recommended Capacity of this program is
n = 512, k = 128.
"""
    )
    parser.add_argument("-k", help="number of clusters", type=int, required=True)
    parser.add_argument("-n", help="number of data points", type=int, required=True)
    parser.add_argument("--Random", help="randomize parameters", action="store_true")

    args = parser.parse_args()

    initial_k = k = args.k
    n = args.n
    r = args.Random
    d = random.choice([2,3])
    m = KMPP_MAX_ITER

    if r:
        n = random_n()
        initial_k = random_k(n)
        k = 0

    if not valid_parameters(initial_k, n, r):
        print("invalid parameters")
        return

    print(f"Generating data with paramters n={n}, k={initial_k}, d={d}")
    data, clusters = generate_data(initial_k, n, d)

    print(f"Writing data to ./data.txt")
    with open("./data.txt", "w") as f:
        for i in range(n):
            f.write(", ".join(list(map(lambda d : "%.8f" % d, data[i, :])) + [str(clusters[i])]))
            f.write("\n")

    print(f"Clustering data with Normalized Spectral Clustering and K-means++")
    (k,
     binary_kmpp_labels, kmpp_jaccard,
     binary_nsc_labels, nsc_jaccard) = nsc_and_kmpp(k, n, d, m,
                                                    data.flatten().tolist(),
                                                    clusters.flatten().tolist())
    kmpp_labels = np.array(zip(*struct.iter_unpack("Q", binary_kmpp_labels)).__next__())
    nsc_labels = np.array(zip(*struct.iter_unpack("Q", binary_nsc_labels)).__next__())
    print("Done in %.3f seconds" % (time.time() - start_time))

    print(f"Writing calculated clusters to ./clusters.txt")
    with open("./clusters.txt", "w") as f:
        f.write(f"{k}\n")

        cluster_members = {cluster : [] for cluster in range(k)}
        for i in range(n):
            cluster_members[nsc_labels[i]].append(i)
        for cluster in range(k):
            f.write(", ".join(map(str, cluster_members[cluster])))
            f.write("\n")

        cluster_members = {cluster : [] for cluster in range(k)}
        for i in range(n):
            cluster_members[kmpp_labels[i]].append(i)
        for cluster in range(k):
            f.write(", ".join(map(str, cluster_members[cluster])))
            f.write("\n")

    desc = f"""
    Clustered data was generated with parameters n={n}, k={initial_k}, d={d}
    Both algorithms created K={k} clusters
    """
    output_clusters_pdf(data, n, d, nsc_labels, nsc_jaccard, kmpp_labels, kmpp_jaccard, desc);


if __name__ == "__main__":
    main()
