import os
import platform
import re

import math
import numpy as np
import numpy.ma as ma
import pandas as pd
from scipy.interpolate import griddata
from datetime import datetime

import ogr, osr

__sysname = platform.system()
__version = "v1.0.0"


class UnsupportedDimensionsError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def _PatternConvert(pattern):
    pattern = pattern.replace(".", "\.")
    pattern = pattern.replace("*", ".*")
    pattern = pattern + "$"
    return pattern


# --------------------------------------

def GetGLT(flon, flat, flipup=True):
    Lon = np.loadtxt(flon, skiprows=6)
    Lat = np.loadtxt(flat, skiprows=6)
    if flipup:
        Lon = np.flipud(Lon)
        Lat = np.flipud(Lat)
    return Lon, Lat


# --------------------------------------
def Grid(value, Opoints, Npoints, method='nearest'):
    '''

    :param value: numpy.ndarray,
    :param Opoints:tuple that contains lat and long array of the value
    :param Npoints:tuple that contains lat and long array of the target value.
    :param method:default for nearest
    :return:
    '''
    points = np.array([Opoints[0].flatten(), Opoints[1].flatten()]).T
    nr, nc = Npoints[0].shape[0], Npoints[0].shape[1]

    if value.ndim == 2:
        Nvar = griddata(points, \
                        value.flatten(), \
                        Npoints, method=method)
    elif value.ndim == 3:
        Nvar = np.zeros((value.shape[0], nr, nc))
        for i, valuei in enumerate(value):
            Nvar[i] = griddata(points, \
                               valuei.flatten(), \
                               Npoints, method=method)
    else:
        raise UnsupportedDimensionsError('Dimension %d is unsupported.' % value.ndim)

    return Nvar


def AsciiCoord(AscFiles):
    '''

    :param AscFiles:
    :param r:
    :param c:
    :return:
    '''

    df = pd.read_table(AscFiles, nrows=6, sep='\s+', names=['name', 'value'])

    GeoTransform = (df.iloc[2, 1], df.iloc[4, 1], 0.0,
                    df.iloc[3, 1], 0.0, df.iloc[4, 1])

    r = int(df.iloc[1, 1])
    c = int(df.iloc[0, 1])

    coord = np.zeros((2, r, c))
    for i in range(r):
        for j in range(c):
            coord[1, i, j] = GeoTransform[0] + (j + 0.5) * GeoTransform[1]  # x
            coord[0, i, j] = GeoTransform[3] + (i + 0.5) * GeoTransform[5]  # y

    coord[0] = np.flipud(coord[0])

    return coord


def AsciiMask(MaskFile, attri=False):
    try:
        mask = np.loadtxt(MaskFile, skiprows=6, dtype=np.int)
        Info = pd.read_table(MaskFile, sep='\s+', nrows=6, names=['name', 'values'])
        mask = ma.array(mask, mask=(mask == Info.iloc[5, 1]))

        if attri:
            return mask, Info
        else:
            return mask
    except:
        raise ValueError('Unrecognized format!')


# ------------------------------------------
def IntpNaN2d(X, axis=0, method='linear'):
    if axis >= 2 or X.ndim != 2:
        raise ValueError('dimension or aixs error.')
    df = pd.DataFrame(X)
    X = df.interpolate(axis=axis, method=method)
    return X.get_values()


def IntpNaN3d(X, axis=0, method='linear'):
    if axis >= 3 or X.ndim != 3:
        raise ValueError('dimension or aixs error.')

    Y = X[:]
    if axis == 1 or axis == 2:
        for i, Xi in enumerate(X[:]):
            Xi = IntpNaN2d(Xi, axis=axis, method=method)
            Y[i] = Xi
    else:
        X = X.transpose((2, 1, 0))
        Y = X[:]
        for i, Xi in enumerate(X[:]):
            Xi = IntpNaN2d(Xi, axis=axis, method=method)
            Y[i] = Xi
        Y = Y.transpose((2, 1, 0))
    return Y


def MaskArrayFillNaN(Var):
    fv = Var.fill_value
    Var = Var.data
    Var[Var == fv] = np.nan

    if Var.ndim == 2:
        Var = IntpNaN2d(Var)
    if Var.ndim == 3:
        Var = IntpNaN3d(Var)
    return Var


def coord_conv(Point, in_proj4, out_proj4):
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromProj4(out_proj4)

    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromProj4(in_proj4)

    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(Point[0], Point[1])

    # create coordinate transformation

    coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    point.Transform(coordTransform)

    return point.GetX(), point.GetY()


def check_within(Point, proj4, shp):
    '''
        To check whether a point is within a shape file.
    :param Point:list or like
        [Lon,Lat]
    :param shp:str.
        shape file name,must have only one feature.
    :return:True or False
        if within return true,else false.
    '''

    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromProj4(proj4)

    driver = ogr.GetDriverByName('ESRI Shapefile')
    polyshp = driver.Open(shp)
    polylyr = polyshp.GetLayer(0)
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(Point[0], Point[1])

    point.AssignSpatialReference(inSpatialRef)

    poly_feature = polylyr.GetNextFeature()
    if point.Within(poly_feature.geometry()):
        return True
    else:
        return False


# ------------------------------------------

class Sun:

    def __init__(self, coords):
        self.coords = coords

    def getSunriseTime(self, time=False):
        return self.calcSunTime(isRiseTime=True, time=time)

    def getSunsetTime(self, time=False):
        return self.calcSunTime(isRiseTime=False, time=time)

    def getCurrentUTC(self):
        now = datetime.now()
        return [now.day, now.month, now.year]

    def calcSunTime(self, isRiseTime, zenith=90.8, time=False):

        # isRiseTime == False, returns sunsetTime
        if not time:
            time = self.getCurrentUTC()
        else:
            time = datetime.strptime(time, "%Y-%m-%d").day, \
                   datetime.strptime(time, "%Y-%m-%d").month, \
                   datetime.strptime(time, "%Y-%m-%d").year

        day, month, year = time[0], time[1], time[2]

        longitude = self.coords['longitude']
        latitude = self.coords['latitude']

        TO_RAD = math.pi / 180

        # 1. first calculate the day of the year
        N1 = math.floor(275 * month / 9)
        N2 = math.floor((month + 9) / 12)
        N3 = (1 + math.floor((year - 4 * math.floor(year / 4) + 2) / 3))
        N = N1 - (N2 * N3) + day - 30

        # 2. convert the longitude to hour value and calculate an approximate time
        lngHour = longitude / 15

        if isRiseTime:
            t = N + ((6 - lngHour) / 24)
        else:  # sunset
            t = N + ((18 - lngHour) / 24)

        # 3. calculate the Sun's mean anomaly
        M = (0.9856 * t) - 3.289

        # 4. calculate the Sun's true longitude
        L = M + (1.916 * math.sin(TO_RAD * M)) + (0.020 * math.sin(TO_RAD * 2 * M)) + 282.634
        L = self.forceRange(L, 360)  # NOTE: L adjusted into the range [0,360)

        # 5a. calculate the Sun's right ascension

        RA = (1 / TO_RAD) * math.atan(0.91764 * math.tan(TO_RAD * L))
        RA = self.forceRange(RA, 360)  # NOTE: RA adjusted into the range [0,360)

        # 5b. right ascension value needs to be in the same quadrant as L
        Lquadrant = (math.floor(L / 90)) * 90
        RAquadrant = (math.floor(RA / 90)) * 90
        RA = RA + (Lquadrant - RAquadrant)

        # 5c. right ascension value needs to be converted into hours
        RA = RA / 15

        # 6. calculate the Sun's declination
        sinDec = 0.39782 * math.sin(TO_RAD * L)
        cosDec = math.cos(math.asin(sinDec))

        # 7a. calculate the Sun's local hour angle
        cosH = (math.cos(TO_RAD * zenith) - (sinDec * math.sin(TO_RAD * latitude))) / (
                    cosDec * math.cos(TO_RAD * latitude))

        if cosH > 1:
            return {'status': False, 'msg': 'the sun never rises on this location (on the specified date)'}

        if cosH < -1:
            return {'status': False, 'msg': 'the sun never sets on this location (on the specified date)'}

        # 7b. finish calculating H and convert into hours

        if isRiseTime:
            H = 360 - (1 / TO_RAD) * math.acos(cosH)
        else:  # setting
            H = (1 / TO_RAD) * math.acos(cosH)

        H = H / 15

        # 8. calculate local mean time of rising/setting
        T = H + RA - (0.06571 * t) - 6.622

        # 9. adjust back to UTC
        UT = T - lngHour
        UT = self.forceRange(UT, 24)  # UTC time in decimal format (e.g. 23.23)

        # 10. Return
        hr = self.forceRange(int(UT), 24)
        min = round((UT - int(UT)) * 60, 0)

        return {
            'status': True,
            'decimal': UT,
            'hr': hr,
            'min': min
        }

    def forceRange(self, v, max):
        # force v to be >= 0 and < max
        if v < 0:
            return v + max
        elif v >= max:
            return v - max

        return v


if __name__ == "__main__":
    print(check_within([90, 30], r"D:\code\indp_val\shp\test_extent.shp"))
