"""Generate 10k train + 500 test samples -> data/{train,test,stats}.npz."""
import os
import numpy as np
from sbw import (transmission_scenery, sample_params, build_scene,
                 signal_from_scene, NL, N)

OUT = "data"
N_TRAIN, N_TEST = 10000, 500


def make(n, T_flat, rng):
    X = np.empty((n, NL, N), np.float32)
    Y = np.empty((n, 6), np.float32)
    for i in range(n):
        p = sample_params(rng)
        X[i] = signal_from_scene(build_scene(p), T_flat)
        Y[i] = p
        if (i + 1) % 1000 == 0:
            print(f"  {i + 1}/{n}")
    return X, Y


def main():
    os.makedirs(OUT, exist_ok=True)
    T_flat = transmission_scenery()
    rng = np.random.default_rng(0)
    print("train:")
    Xtr, Ytr = make(N_TRAIN, T_flat, rng)
    print("test:")
    Xte, Yte = make(N_TEST, T_flat, rng)
    np.savez(f"{OUT}/train.npz", X=Xtr, Y=Ytr)
    np.savez(f"{OUT}/test.npz", X=Xte, Y=Yte)
    # global stats preserve absolute amplitude (needed to recover amplitudes)
    np.savez(f"{OUT}/stats.npz", mean=Xtr.mean(), std=Xtr.std())
    print(f"done. mean={Xtr.mean():.3g} std={Xtr.std():.3g}")


if __name__ == "__main__":
    main()
