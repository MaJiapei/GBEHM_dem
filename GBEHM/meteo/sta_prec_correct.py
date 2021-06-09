import numpy as np

tl = 2
ts = -2


def cspg(prec, w, tmax, tmin, td):
    if np.isnan(np.asarray([prec, w, tmax, tmin, td], dtype=float)).any():
        return np.asarray([prec, w, tmax, tmin, td], dtype=float)[0]

    if 0.01 < prec < 0.1:
        return 0.1

    if prec < 0.01:
        return 0

    if tmin >= tl:
        CR = np.exp(-0.041 * w) * 100
        Dpw = 0.29

    elif ts >= tmax:
        CR = np.exp(-0.056 * w) * 100
        Dpw = 0.3

    else:
        CRsnow = np.exp(-0.056 * w) * 100
        CRrain = np.exp(-0.041 * w) * 100
        CR = CRsnow - (CRsnow - CRrain) * (td - tl) / (tl - ts)
        Dpw = 0.29

    return 100 * (1 / CR) * (prec + Dpw)


def tretyakov(prec, w, tmax, tmin):
    if np.isnan(np.asarray([prec, w, tmax, tmin], dtype=float)).any():
        return np.asarray([prec, w, tmax, tmin], dtype=float)[0]

    if 0.01 < prec < 0.1:
        return 0.1
    if prec < 0.01:
        return 0

    ws = w * (np.log(0.4 / 0.01) / np.log(10 / 0.01))
    if tmin >= tl:
        CR = 100 - 4.77 * ws ** 0.56
        Dpw = 0.2

    elif tmax <= ts:
        CR = 103.1 - 8.67 * ws + 0.30 * tmax
        Dpw = 0.15

    else:
        CR = 96.99 - 4.46 * ws + 0.88 * tmax + 0.22 * tmin
        Dpw = 0.15

    return 100 * (1 / CR) * (prec + +Dpw)


def india(prec, w, tmax, tmin):
    if np.isnan(np.asarray([prec, w, tmax, tmin], dtype=float)).any():
        return np.asarray([prec, w, tmax, tmin], dtype=float)[0]

    if 0.01 < prec < 0.1:
        return 0.1
    if prec < 0.01:
        return 0

    ws = w * (np.log(0.3 / 0.01) / np.log(10 / 0.01))
    if tmin >= tl:
        CR = 100 - 4.77 * ws ** 0.56
        Dpw = 0.2

    elif tmax <= ts:
        CR = 103.1 - 8.67 * ws + 0.30 * tmax
        Dpw = 0.15

    else:
        CR = 96.99 - 4.46 * ws + 0.88 * tmax + 0.22 * tmin
        Dpw = 0.15

    return 100 * (1 / CR) * (prec + +Dpw)


def mk2(prec, tmax, tmin):
    if np.isnan(np.asarray([prec, tmax, tmin], dtype=float)).any():
        return np.asarray([prec, tmax, tmin], dtype=float)[0]

    if 0.01 < prec < 0.1:
        return 0.1
    if prec < 0.01:
        return 0

    if tmin >= tl:
        CR = 100 / 1.05
        Dpw = 0.2

    elif tmax <= ts:
        CR = 100 / 1.13
        Dpw = 0.15

    else:
        CR = 100 / 1.05
        Dpw = 0.15

    return 100 * (1 / CR) * (prec + +Dpw)


def uswb8(prec, w, tmax, tmin):
    if np.isnan(np.asarray([prec, w, tmax, tmin], dtype=float)).any():
        return np.asarray([prec, w, tmax, tmin], dtype=float)[0]

    if 0.01 < prec < 0.1:
        return 0.1
    if prec < 0.01:
        return 0

    ws = w * (np.log(0.7 / 0.01) / np.log(10 / 0.01))
    if tmin >= tl:
        CR = np.exp(4.605 - 0.062 * ws ** 0.58)
        Dpw = 0.03

    elif tmax <= ts:
        CR = np.exp(4.606 - 0.157 * ws ** 1.28)
        Dpw = 0.15

    else:
        CR = 100.77 - 8.35 * ws
        Dpw = 0.15

    return 100 * (1 / CR) * (prec + +Dpw)


def nepal(prec, w, tmax, tmin):
    if np.isnan(np.asarray([prec, w, tmax, tmin], dtype=float)).any():
        return np.asarray([prec, w, tmax, tmin], dtype=float)[0]

    if 0.01 < prec < 0.1:
        return 0.1
    if prec < 0.01:
        return 0

    ws = w * (np.log(1 / 0.01) / np.log(10 / 0.01))
    if tmin >= tl:
        CR = np.exp(4.605 - 0.062 * ws ** 0.58)
        Dpw = 0.2

    elif tmax <= ts:
        CR = np.exp(4.606 - 0.157 * ws ** 1.28)
        Dpw = 0.15

    else:
        CR = 100.77 - 8.35 * ws
        Dpw = 0.15

    return 100 * (1 / CR) * (prec + +Dpw)


if __name__ == "__main__":
    import pandas as pd

    df = pd.read_pickle(r"D:\data\MeteoData\XJ\Precip_corrected_update.pik")
