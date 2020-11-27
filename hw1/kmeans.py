#!/bin/env python3

import argparse

class Observation(object):
    def __init__(self, dimension):
        self.vector = [0.] * dimension
        self.cluster = 0

    def set_vector(self, vector):
        self.vector = vector

class Cluster(object):
    def __init__(self, dimension):
        self.dimension = dimension
        self.centroid = [0.] * self.dimension
        self.cardinality = 0

    def set_center(self, vector):
        self.centroid = vector

    def __eq__(self, other):
        return self.centroid == other.centroid

    def distance_squared(self, ob):
        vector = ob.vector
        return sum(((self.centroid[i] - vector[i])**2 for i in range(self.dimension)))

    def add(self, ob):
        vector = ob.vector
        self.centroid = [
                (self.centroid[i] * self.cardinality + vector[i]) / (self.cardinality + 1)
                for i in range(self.dimension)
                ]
        self.cardinality += 1


class KmContext(object):
    def __init__(self, k, n, d, m):
        self.num_observed = 0
        self.num_clusters = k
        self.num_observations = n
        self.dimension = d
        self.max_iterations = m
        self.clusters = [Cluster(self.dimension) for _ in range(self.num_clusters)]
        self.observations = [Observation(self.dimension) for _ in range(self.num_observations)]

    def observe(self, vector):
        if self.num_observed < self.num_clusters:
            self.clusters[self.num_observed].set_center(vector)
        self.observations[self.num_observed].set_vector(vector)
        self.num_observed += 1

    def choose_cluster(self, ob):
        res = 0
        dis = self.clusters[0].distance_squared(ob)
        for i, cluster in enumerate(self.clusters[1:]):
            d = cluster.distance_squared(ob)
            if d < dis:
                dis = d
                res = i + 1
        return res

    def converge(self):
        for _ in range(self.max_iterations):
            next_clusters = [Cluster(self.dimension) for _ in range(self.num_clusters)]
            for i in range(self.num_observations):
                curr_ob = self.observations[i]
                curr_cluster_index = self.choose_cluster(curr_ob)
                next_clusters[curr_cluster_index].add(curr_ob)
            if all(map(lambda i: self.clusters[i] == next_clusters[i], range(self.num_clusters))):
                return 1
            else:
                self.clusters = next_clusters
        return 0

    def dump(self):
        for cluster in self.clusters:
            print(",".join(map(lambda f: str(round(f, 2)), cluster.centroid)))


def scanner():
    while True:
        try:
            yield list(map(float, input().split(",")))
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

