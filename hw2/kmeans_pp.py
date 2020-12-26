#!/bin/env python3

import argparse

import pandas as pd
import numpy as np

from mykmeanspp import fit


def k_means_pp(df, k, n, d, m):
    np.random.seed(0)

    centroid_indices = [0 for _ in range(k)]
    centroids = np.full((k, d), 0)

    centroid_indices[0] = np.random.choice(n, 1)[0]
    centroids[0] = df.iloc[centroid_indices[0], :]

    for j in range(1, k):
        min_distances = [np.min([np.sum(np.square(df.iloc[i, :] - centroids[k])) for k in range(j)]) for i in range(n)]
        sum_distances = np.sum(min_distances)
        centroid_indices[j] = np.random.choice(n, 1, p=min_distances / sum_distances)[0]
        centroids[j] = df.iloc[centroid_indices[j], :]

    print(",".join(map(str,centroid_indices)))
    fit(k, n, d, m, df.to_numpy().flatten().tolist(), list(map(int, centroid_indices)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("K", type=int)
    parser.add_argument("N", type=int)
    parser.add_argument("d", type=int)
    parser.add_argument("MAX_ITER", type=int)
    parser.add_argument("filename")
    args = parser.parse_args()

    k = args.K
    n = args.N
    d = args.d
    m = args.MAX_ITER
    df = pd.read_csv(args.filename, header=None)

    k_means_pp(df, k, n, d, m)

    return


if __name__ == "__main__":
    main()

