import re
import warnings

import numpy as np
import numpy.ma as ma

from netCDF4 import Dataset, num2date, date2num
from GBEHM.Tool import MaskArrayFillNaN

import xarray as xr
import rasterio


def CoordVar(fnc, varName=None, deNaN=False, fillValue=None):
    '''

    :param fnc:
    :param varName: default none
    :return:
    '''
    nf = Dataset(fnc)
    vars = nf.variables.keys()
    tag = 0
    time = 0
    mark = False
    for var in vars:
        if re.search(r"lat", var, re.I):
            lat = nf[var][:]
            tag += 1
        if re.search(r'lon', var, re.I):
            lon = nf[var][:]
            tag += 1

        if re.search(r'time', var, re.I):
            try:
                time = nf[var]
                time = num2date(time[:], units=time.units)
            except:
                time = -2
                warnings.warn('Error in reading variable time!')

        if varName == var:
            Var = nf[var][:]
            # test if the specified variable is masked.
            # if masked,the variable will be interpolated to remove the NaN value.
            if deNaN:
                if isinstance(Var, np.ma.core.MaskedArray):
                    Var = MaskArrayFillNaN(Var)
                else:
                    if fillValue != None:
                        Var = ma.array(Var, mask=Var == fillValue)
                        Var = MaskArrayFillNaN(Var)

            mark = True

    nf.close()

    if tag != 2:
        raise IndexError("Can't Find Coordinate Infomation!")

    if isinstance(time, int) and time != -2:
        warnings.warn("Time dimension not found!")

    if lat.ndim != lon.ndim:
        raise IndexError(str(lat.shape) + 'unequal to' + str(lon.shape))

    if varName:
        if not mark:
            raise ValueError("No Variable:%s" % varName)
        else:
            return lat, lon, time, Var
    else:
        return lat, lon, time


def PosIndex(point, fnc=None, lat=None, lon=None):
    '''

    :param point: list contains latitude and longitude of the point.
    :param fnc: name of the input netCDF file.
    :param lat:
    :param lon:
    :return: Index (Y,X) of the point in the coordinate.
    '''
    if fnc:
        lat, lon, time = CoordVar(fnc)

    CoordDim = lat.ndim
    if CoordDim == 1:
        Lon, Lat = np.meshgrid(lon, lat)
        dis = (Lat - point[1]) ** 2 + (Lon - point[0]) ** 2
        pos = np.where(dis == np.min(dis))
        Ypos, Xpos = pos[0][0], pos[1][0]
        return Ypos, Xpos

    if CoordDim == 2:
        dis = (lat - point[1]) ** 2 + (lon - point[0]) ** 2
        pos = np.where(dis == np.min(dis))
        Ypos, Xpos = pos[0][0], pos[1][0]
        return Ypos, Xpos


def WriteArray2NC(arr, outName, lon=None, lat=None, time=None,
                  format='NETCDF3_CLASSIC', arrUnits=None, arrName=None, fillvalue=-9999.):
    '''

    :param arr:
    :param outName:
    :param lon:
    :param lat:
    :param time:
    :param format:
    :param arrUnits:
    :param arrName:
    :param fillvalue:
    :return:
    '''

    # data dimensions check
    if arr.ndim == 2:
        rows, cols, nbands = arr.shape, 1
    elif arr.ndim == 3:
        rows, cols, nbands = arr.shape
    else:
        raise NotImplementedError('Data having dimensions more than 3 is not supported now.')

    if time:
        if time.__len__() != nbands:
            raise ValueError('Dimensions of time is not equal to the third dimensions of data.')

    # coordinate check
    if lon and lat:
        if lon.ndim == 1:
            if lon.__len__() != cols or lat.__len__() != rows:
                raise ValueError('Dimensions of coord is not equal to that of data.')
        elif lon.ndim == 2:
            if lon.shape != [rows, cols] or lat.shape != [rows, cols]:
                raise ValueError('Dimensions of coord is not equal to that of data.')
    #
    #
    ncf = Dataset(outName, "w", format=format)
    ncf.createDimension('lon', cols)
    ncf.createDimension('lat', rows)
    if time:
        ncf.createDimension('time', None)
        TIME = ncf.createVariable('time', 'f8', ('time'))
        TIME[:] = [date2num(i, units='hours since 2000-01-01 00:00:00') for i in time]
        TIME.long_name = 'time'
        TIME.uints = 'hours since 2000-01-01 00:00:00'
        TIME.axis = 'T'

        var = ncf.createVariable('var', 'f4', ('time', 'lat', 'lon'), fill_value=fillvalue)
        var[:] = arr
        var.long_name = arrName
        var.units = arrUnits

    else:
        if arr.ndim == 3:
            var = ncf.createVariable('var', 'f4', ('time', 'lat', 'lon'), fill_value=fillvalue)
            var[:] = arr
            var.long_name = arrName
            var.units = arrUnits
        elif arr.ndim == 2:
            var = ncf.createVariable('var', 'f4', ('lat', 'lon'), fill_value=fillvalue)
            var[:] = arr
            var.long_name = arrName
            var.units = arrUnits

    if lon and lat:
        if lon.ndim == 1 and lat.ndim == 1:
            LON = ncf.createVariable('lon', 'f4', ('lon'))
            LAT = ncf.createVariable('lat', 'f4', ('lat'))
            LON[:] = Lon
            LAT[:] = Lat

        elif lon.ndim == 2 and lat.ndim == 2:
            LON = ncf.createVariable('lon', 'f4', ('lat', 'lon'))
            LAT = ncf.createVariable('lat', 'f4', ('lat', 'lon'))
            LON[:] = Lon
            LAT[:] = Lat

        LON.long_name = 'longitude'
        LON.units = 'degree_east'
        LON.axis = 'X'

        LAT.long_name = 'latitude'
        LAT.units = 'degrees_north'
        LAT.axis = 'Y'

    ncf.close()


if __name__ == "__main__":
    Fname = r"C:\Users\jpma\Downloads\wrfout_heihe_2012-06-01.nc"
    Oname = "../test.tiff"
