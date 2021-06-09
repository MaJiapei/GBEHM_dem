import numpy as np


def value_check(value, varname):
    # waiting for implementation.
    valueRange = {'wsp10m': [0, 100],
                  'U10m': [0, 100],
                  'U10': [0, 100],
                  'V10m': [0, 100],
                  'temp2m': [200, 340],
                  'T2': [200, 340],
                  'pres': [1000, 120000],
                  'rh2m': [0.01, 1],
                  'q2m': [0.001, 0.1],
                  'swr': [0, 1500],
                  'lwr': [10, 1000],
                  'prec': [0, 500]}
    if valueRange[varname][0] <= value <= valueRange[varname][1]:
        return True
    elif np.isnan(value):
        return -1
    else:
        return False
