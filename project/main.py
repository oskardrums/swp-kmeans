# main.py
# Clustering comparison main entry point

import argparse
from sklearn.datasets import make_blobs
import random
import struct
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from clustering import nsc, kmpp

KMPP_MAX_ITER = 300

# Estimates of Maximum Capacity within 5 minutes
MAX_N = 550
# We can handle a much larger K value within 5 minutes,
# but we don't want too many clusters if we don't have a lot of samples.
MAX_K = 8


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

    else:
        return

    plt.figtext(0.5, 0.01, desc, wrap=True, horizontalalignment='center', fontsize=12)

    fig.set_size_inches(7, 8, forward=True)
    plt.savefig("clusters.pdf")


def write_labels(n, k, labels, f):
    """
    Writes clusters labels to file f
    """
    cluster_members = {cluster : [] for cluster in range(k)}
    
    for i in range(n):
        cluster_members[labels[i]].append(i)
        
    for cluster in range(k):
        f.write(b", ".join(map(lambda x: str(x).encode(), cluster_members[cluster])))
        f.write(b"\n")


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


def valid_parameters(k, n):
    """
    Validates input parameters
    """
    return 0 < k < n

def normalized_spectral_clustering(samples, labels, k=0):
    """
    Takes an array of n d-dimensional observations and
    returns 3-tuple where the 1st elemenet the number of 
    calculated clusters, the 2nd elemenet is an array of n
    labels such that the i'th label is the number of the cluster 
    matched to the i'th observation, and the 3nd elemenet is the 
    Jaccard Measure of the given labels and the calculated labels.

    Hence, labels is an array of labels considered the source of truth for
    the given samples, and is used to assess the accuracy of the calculated
    labels.

    The parameter k is the number of clusters to create.
    If k is 0, the Eigengap Heuristic is used to determine a viable candidate.
    """
    n, d = samples.shape
    k, binary_labels, jaccard_measure =  nsc(k, n, d, 300,
                                             samples.flatten().tolist(),
                                             labels.flatten().tolist())
    calculated_labels = np.array(zip(*struct.iter_unpack("Q", binary_labels)).__next__())
    
    return k, calculated_labels, jaccard_measure


def k_means_pp(samples, labels, k):
    """
    Takes an array of n d-dimensional observations and
    returns 2-tuple where the 1st elemenet is an array of n
    labels such that the i'th label is the number of the cluster 
    matched to the i'th observation, and the 2nd elemenet is the 
    Jaccard Measure of the given labels and the calculated labels.

    Hence, labels is an array of labels considered the source of truth for
    the given samples, and is used to assess the accuracy of the calculated
    labels.

    The parameter k is the number of clusters to create.
    """
    n, d = samples.shape
    binary_labels, jaccard_measure =  kmpp(k, n, d, 300,
                                           samples.flatten().tolist(),
                                           labels.flatten().tolist())
    calculated_labels = np.array(zip(*struct.iter_unpack("Q", binary_labels)).__next__())
    
    return calculated_labels, jaccard_measure
    

def main():
    parser = argparse.ArgumentParser(
        description="""
        Clustering Algorithms Comparison
        This program generates random data based on the given parameters and clusters the data
        with both K-means++ algorithm and Normalized Spectral Clustering

        For completion within 5 minutes, the Maximum Recommended Capacity of this program is
        n = 512, k = 128.
        """
    )
    parser.add_argument("-k", help="number of clusters", type=int)
    parser.add_argument("-n", help="number of data points", type=int)
    parser.add_argument("-d", help="dimension of data points", type=int)
    parser.add_argument("-m", help="maximum number of iterations to perform for kmeans", type=int)
    parser.add_argument("-r", help="randomize parameters", action="store_true")

    args = parser.parse_args()

    initial_k = k = args.k or 0
    n = args.n or 0
    r = args.r or False
    d = args.d or random.choice([2,3])
    m = args.m or KMPP_MAX_ITER

    if r:
        n = random_n()
        initial_k = random_k(n)
        k = 0

    if not valid_parameters(initial_k, n):
        print("invalid parameters")
        parser.parse_args(["-h"])
        return

    print(f"Generating data with paramters n={n}, k={initial_k}, d={d}")
    data, clusters = generate_data(initial_k, n, d)

    print(f"Writing data to ./data.txt")
    with open("./data.txt", "wb") as f:
        for i in range(n):
            f.write(b", ".join(list(map(lambda d : b"%.8f" % d, data[i, :])) + [str(clusters[i]).encode()]))
            f.write(b"\n")

    print(f"Clustering data with Normalized Spectral Clustering")
    k, nsc_labels, nsc_jaccard = normalized_spectral_clustering(data, clusters, k)
    
    print(f"Clustering data with K-means++")
    kmpp_labels, kmpp_jaccard = k_means_pp(data, clusters, k)
    
    print(f"Writing calculated clusters to ./clusters.txt")
    with open("./clusters.txt", "wb") as f:
        f.write(b"%u\n" % k)
        write_labels(n, k, nsc_labels, f)
        write_labels(n, k, kmpp_labels, f)

    print(f"Drawing visualization at ./clusters.pdf")
    desc = f"""
    Clustered data was generated with parameters n={n}, k={initial_k}, d={d}
    Both algorithms created K={k} clusters
    """
    output_clusters_pdf(data, n, d, nsc_labels, nsc_jaccard, kmpp_labels, kmpp_jaccard, desc);


if __name__ == "__main__":
    main()
