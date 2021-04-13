import argparse
from sklearn.datasets import make_blobs
import random
import struct

from prj_lib import prj_main

def generate_data(k, n, d):
    return make_blobs(n_samples=n, n_features=d, centers=k)
    
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
    print(data)
    
    binary_labels = prj_main(k, n, d, m, data.tobytes(), clusters.tobytes())

    labels = zip(*struct.iter_unpack("Q", binary_labels)).__next__()

    print(labels)
    
if __name__ == "__main__":
    main()
