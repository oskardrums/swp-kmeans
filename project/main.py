import argparse
from sklearn.datasets import make_blobs
import random
import struct

from prj_lib import prj_main

def generate_data(k, n, d):
    return make_blobs(n_samples=n, n_features=d, centers=k)

def random_n():
    n = 1<<7
    return random.randint(int(n / 2), n)

def random_k(n):
    k = min(n, 1<<4)
    return random.randint(int(k / 2), k)

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
        k = 0
        n = random_n()
        data, clusters = generate_data(random_k(n), n, d)
    else:
        data, clusters = generate_data(k, n, d)
    
    binary_kmpp_labels, binary_nsc_labels = prj_main(k, n, d, m, data.flatten().tolist(), clusters.flatten().tolist())

    kmpp_labels = list(zip(*struct.iter_unpack("Q", binary_kmpp_labels)).__next__())
    nsc_labels = list(zip(*struct.iter_unpack("Q", binary_nsc_labels)).__next__())

    print(kmpp_labels[:16])
    print(nsc_labels[:16])
    print(list(clusters)[:16])
    
if __name__ == "__main__":
    main()
