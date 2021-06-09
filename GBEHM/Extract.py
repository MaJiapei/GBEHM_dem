import numpy as np
import pandas as pd

from GBEHM.ncTool import CoordVar, PosIndex
from GBEHM.Tool import GetGLT, Grid


def _Pos1D(Lat, Lon, lat, lon, Var):
    extent = np.array([np.max(lat), np.min(lat), np.min(lon), np.max(lon)])
    VarDim = Var.ndim
    index_lat = np.where((Lat > extent[1]) & (Lat < extent[0]))
    index_lon = np.where((Lon > extent[2]) & (Lon < extent[3]))
    lat_range = [index_lat[0][0] - 5, index_lat[0][-1] + 5]
    lon_range = [index_lon[0][0] - 5, index_lon[0][-1] + 5]

    lon_nc, lat_nc = np.meshgrid(Lon[lon_range[0]:lon_range[1]], \
                                 Lat[lat_range[0]:lat_range[1]])
    points = np.array([lat_nc.flatten(), lon_nc.flatten()]).T

    if VarDim == 2:
        var = Var[lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]]
    if VarDim == 3:
        var = Var[:, lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]]

    return points, var


def byRaster(Fnc, Flon, Flat, var, method="nearest"):
    Lat, Lon, Time, Var = CoordVar(Fnc, varName=var)
    lon, lat = GetGLT(Flon, Flat)

    CoordDim = Lat.ndim

    if CoordDim == 1:
        points, var = _Pos1D(Lat, Lon, lat, lon, Var)
        Nvar = Grid(var, (points[:, 0], points[:, 1]), (lat, lon), method=method)
    if CoordDim == 2:
        pass

    return Nvar


def byPoints(point, fnc, varName):
    '''
    Using  points to extract the corresponding value in the netCDF file.

    :point: list
        contains latitude and longitude of the point[[lon,lat]].
    :fnc:str
        name of the input netCDF file.
    :varName: str
        name of the variable.
    :return:a panda dataframe contains the value in the position of the points.
    '''
    lat, lon, time, Var = CoordVar(fnc, varName=varName)

    VarDim, CoordDim = Var.ndim, lat.ndim

    npoints = len(point)

    if VarDim == 2:
        df = pd.DataFrame(columns=np.arange(npoints))
        for i, p in enumerate(point):
            index = PosIndex(p, lat=lat, lon=lon)
            Value = Var[index[0], index[1]]
            df.iloc[:, i] = Value
    if VarDim == 3:
        if isinstance(time, np.ndarray):
            df = pd.DataFrame(index=time, columns=np.arange(npoints))
            for i, p in enumerate(point):
                index = PosIndex(p, lat=lat, lon=lon)
                Value = Var[:, index[0], index[1]]
                df.iloc[:, i] = Value
        else:
            df = pd.DataFrame(columns=np.arange(npoints))
            for i, p in enumerate(point):
                index = PosIndex(p, lat=lat, lon=lon)
                Value = Var[:, index[0], index[1]]
                df.iloc[:, i] = Value
    else:
        raise ValueError('Unsupport Dimensions more than 3.')

    return df


if __name__ == "__main__":
    sitePath = r"D:\work\QilianPrecipCorrection\TibetPlateauSites.csv"
    dataPath = r"D:\work\QilianPrecipCorrection\out.nc"

    sdf = pd.read_csv(sitePath)
    point = np.array(sdf[['Longitude', 'Latitude']]).tolist()
    df = byPoints(point, dataPath, "prcp")
    print('a')
