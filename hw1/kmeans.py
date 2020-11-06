#!/bin/env python3
import argparse

class KmContext(object):
    def __init__(self, k, n, d, m):
        self.k = k
        self.n = n
        self.d = d
        self.m = m
        self.centroids = ...

    def observe(o):


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

