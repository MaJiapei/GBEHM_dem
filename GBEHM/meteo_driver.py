import numpy as np
from datetime import datetime
from netCDF4 import Dataset, num2date, date2num
import pandas as pd
import warnings
from osgeo import gdal, osr, ogr


class Update():
    def __init__(self, ncfile, outfile, txtfile, spref_file):
        '''

        :param ncfile:str
            input netCDF filename
        :param outfile:str
            output filename
        :param txtfile: str
            measurement filename with a format of csv,having a form like below.
            varname,time,x,y,value
            T2,2017-01-01 16:00:00,104.23,30.25,290
            T2,2017-10-21 18:00:00,85.597,43.427,292
            T2,2017-10-21 19:00:00,85.597,43.427,292
            T2,2017-10-21 20:00:00,85.597,43.427,292
            T2,2017-10-21 21:00:00,85.597,43.427,292
            U10,2017-10-22 18:00:00,85.597,43.427,5.2
            hgt,2017-10-02 18:00:00,104.23,30.25,7.2
        :param spref_file:str
            reference raster,coming from basin extraction.
            It must have same extent and spatial projection as input netCDF file

        '''
        self.ncfile = ncfile
        self.ofile = outfile
        self.txtfile = txtfile
        self.extentf = spref_file
        self._nc_check()

    def _nc_check(self):
        f = Dataset(self.ncfile)
        self.f = f

        self.varname = list(f.variables.keys())

        self.times = num2date(f["time"][:], units=f['time'].units)

        ras = gdal.Open(self.extentf)
        Xsize, Ysize = ras.RasterXSize, ras.RasterYSize
        X0, xres, _, Y0, _, yres = ras.GetGeoTransform()
        self.res = xres
        self.x = np.linspace(X0 + xres / 2., X0 + (Xsize - 1) * xres, num=Xsize)
        self.y = np.linspace(Y0 + yres / 2., Y0 + (Ysize - 1) * yres, num=Ysize)

        self.spref = ras.GetProjection()

    def _var_check(self, varname):

        if varname in self.varname:
            return True
        else:
            return False

    def _time_check(self, itime):

        if self.times[0] <= itime <= self.times[-1]:
            return abs((itime - self.times)).argmin()
        else:
            return False

    def _proj_trans(self, Point):
        inputEPSG = 4326

        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(Point[0], Point[1])

        # create coordinate transformation
        inSpatialRef = osr.SpatialReference()
        inSpatialRef.ImportFromEPSG(inputEPSG)

        outSpatialRef = osr.SpatialReference()
        outSpatialRef.ImportFromWkt(self.spref)
        coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
        point.Transform(coordTransform)

        return [point.GetX(), point.GetY()]

    def _coord_check(self, point):
        # checking whether the point in the measurement is within the netCDF extent.
        X, Y = self._proj_trans(point)

        # if the minimum distance is witin a half pixel resolution,
        # the measurement will be thought out of the netCDF extent and will be ignored.
        mdis = np.sqrt(((X - self.x) ** 2).min() + ((Y - self.y) ** 2).min())
        if mdis > 0.5 * 2 ** 0.5 * self.res:
            return False
        else:
            index = [((X - self.x) ** 2).argmin(), ((Y - self.y) ** 2).argmin()]
            return index

    def _value_check(self, varname, value):
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
                      'prec': [0, 10]}
        if valueRange[varname][0] <= value <= valueRange[varname][1]:
            return True
        else:
            return False

    def _txt_check(self):
        df = pd.read_csv(self.txtfile, parse_dates=['time'])

        df = df[df['varname'].apply(self._var_check)]
        df['time_index'] = df['time'].apply(self._time_check)

        f = lambda x: self._coord_check([x['x'], x['y']])
        df['xy_index'] = df.apply(f, axis=1)

        df = df.replace(False, np.nan).dropna().drop(['x', 'y', 'time'], axis=1)

        f = lambda x: self._value_check(x['varname'], x['value'])
        df['value_qc'] = df.apply(f, axis=1)
        df = df[df['value_qc']].drop('value_qc', axis=1)

        df = df.astype({'time_index': 'int'})

        print('%d valid records found in the measurement file.' % df.__len__())
        self.df = df.set_index('varname')

    def direct_replace(self):
        self._txt_check()
        nf = Dataset(self.ofile, "w")
        nf.setncatts(self.f.__dict__)

        # copy dimensions of the input netCDF file.
        for dim in self.f.dimensions.values():
            nf.createDimension(dim.name, dim.size)

        # copy variables of the input netCDF file.
        for var in self.f.variables.values():
            nf.createVariable(var.name, var.dtype.str, var.dimensions)
            # copy attributes
            nf[var.name].setncatts(var.__dict__)

            if var.name not in self.df.index.unique():
                nf[var.name][:] = self.f[var.name][:]
            else:
                print("%s is updating..." % var.name)
                tarr = self.f[var.name][:]
                tdf = self.df.loc[[var.name]]

                for irec in range(tdf.__len__()):
                    # direct replace.
                    tarr[tdf.ix[irec, 'time_index'], tdf.ix[irec,
                                                            'xy_index'][1], tdf.ix[irec, 'xy_index'][0]] = tdf.iloc[
                        irec, 0]
                nf[var.name][:] = tarr[:]
        nf.close()
        self.f.close()


if __name__ == "__main__":
    dup = Update(r"D:\code\update_driver\driver_data\Manas_drivers_XJ_201710.nc",
                 r"D:\code\update_driver\driver_new\temp.nc",
                 r"D:\code\update_driver\measurement.dat",
                 r"D:\code\pythowork\workspace1\station\station0\deminput\bedslope.asc")
    dup.direct_replace()
