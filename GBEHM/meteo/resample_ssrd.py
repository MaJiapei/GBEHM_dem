from netCDF4 import Dataset, num2date, date2num
import pandas as pd
import numpy as np
from datetime import timedelta
from multiprocessing import Pool
from GBEHM.Tool import Sun


def intpssrd(file, oname):
    onf = Dataset(file)
    ossrd = onf['ssrd'][:]
    otime = onf['time']
    lon = onf['longitude'][:]
    lat = onf['latitude'][:]
    oTime = num2date(otime[:], units=otime.units, calendar=otime.calendar)
    index = pd.date_range(oTime[0] - timedelta(hours=6), oTime[-1], freq='H')
    rows, cols = ossrd.shape[1], ossrd.shape[2]
    df = pd.DataFrame(index=index, columns=np.arange(rows * cols))
    for i in range(rows * cols):
        pass

    nf = Dataset(oname, "w")


if __name__ == "__main__":
    intpssrd(r'C:\Users\jiape\Desktop\interim_ssrd.nc', 'test.nc')
