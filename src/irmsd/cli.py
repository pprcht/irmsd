import argparse
import numpy as np
import irmsd

def main():
    p = argparse.ArgumentParser(prog="irmsd")
    p.add_argument("--a", type=float, required=True)
    p.add_argument("--x", nargs="+", type=float, required=True)
    p.add_argument("--y", nargs="+", type=float, required=True)
    args = p.parse_args()

    x = np.array(args.x, dtype=float)
    y = np.array(args.y, dtype=float)

    # call the public API
    out = irmsd.saxpy(args.a, x, y)
    print(" ".join(map(str, out)))

