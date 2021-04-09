import argparse
from sklearn.datasets import make_blobs
import random

from main_lib import doit

def generate_data(k, n, d):
    return make_blobs(n_samples=n, n_features=d, centers=k)

def random_k():
    # TODO - eigengap heuristic
    return random.randint(5,10)

def random_n():
    return random.randint(50,100)

def normalized_spectral_clustering(k, x):
    w = weighted_adjacency_matrix(x)
    l = normalized_graph_laplacian(w)
    k = k or eigengap_heuristic(l)
    u = trim_columns(k, l)
    t = renormalize(u)
    return k, kmeans(k, t)

def kmeans(k, x)
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", help="number of clusters", type=int, required=True)
    parser.add_argument("-n", help="number of data points", type=int, required=True)
    parser.add_argument("--Random", help="random data generation", action="store_true")

    args = parser.parse_args()

    k = args.k
    n = args.n
    r = args.Random
    d = random.choice((2,3))
    m = 300

    if r:
        k = random_k()
        n = random_n()

    data, clusters = generate_data(k, n, d)
    doit(k, n, d, m, data, clusters)
    k, nsc_clusters = normalized_spectral_clustering(k, data)
    
    kmeans_clusters = kmeans(k, data)    

    print(data, nsc_clusters, kmeans_clusters)

if __name__ == "__main__":
    main()
