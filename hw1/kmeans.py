#!/bin/env python3
import argparse

class Observation(object):
    def __init__(self, dimension):
        self.data = None
        self.cluster = 0

    def set_data(self, data):
        self.data = data

class Cluster(object):
    def __init__(self, dimension):
        self.centroid = None
        self.cardinality = 0

    def set_data(self, data):
        self.data = data


class KmContext(object):
    def __init__(self, k, n, d, m):
        self.num_observed = 0
        self.num_clusters = k
        self.num_observations = n
        self.dimension = d
        self.max_iterations = m
        self.clusters = [Cluster(self.dimension)] * self.num_clusters
        self.observations = [Observation(self.dimension)] * self.num_observations

    def observe(self, data):
        if self.num_observed < self.num_clusters:
            self.clusters[self.num_observed].set_data(data)
        self.observations[self.num_observed].set_data(data)

    def converge(self):
        for _ in range(self.max_iterations):
            for i in range(self.num_observations):
                
def scanner():
    while True:
        try:
            yield map(float, input().split(","))
        except EOFError:
            break


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("K", type=int)
    parser.add_argument("N", type=int)
    parser.add_argument("d", type=int)
    parser.add_argument("MAX_ITER", type=int)
    args = parser.parse_args()

    k = args.K
    n = args.N
    d = args.d
    m = args.MAX_ITER

    ctx = KmContext(k, n, d, m)

    for observation in scanner():
        ctx.observe(observation)

    ctx.converge()

    ctx.dump()

    return


if __name__ == "__main__":
    main()

