#!/bin/env python
import argparse
from random import random, choice, randint


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n", type=int)
    parser.add_argument("d", type=int)
    args = parser.parse_args()

    n = args.n
    d = args.d

    for i in range(n):
#        for j in range(n):
#            print("%.2f" % random(), end=",")
#        print()
       # print("hello")
        print(",".join(("%.2f" % (random() * choice([1, -1]) * randint(1, 2**16)) for _ in range(d))))
#        print(",".join(("%.2f" % random() for _ in range(d))))
    
if __name__ == "__main__":
    main()

