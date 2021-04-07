import argparse
from sklearn.datasets import make_blobs
import random

def generate_data(k, n, d):
    print(k, n, d)
    x, _ = make_blobs(n_samples=n, n_features=d, centers=k)
    return x

def random_k():
    return random.randint(5,10)

def random_n():
    return random.randint(50,100)

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

    if r:
        k = random_k()
        n = random_n()

    data = generate_data(k, n, d)

    print(data)

if __name__ == "__main__":
    main()
