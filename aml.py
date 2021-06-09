import arcpy
from arcpy import env
from arcpy.sa import *
import os
import warnings
from numpy import ma
from tqdm import tqdm

from GBEHM.DataStdr import dataCheck


class deminput:
    def __init__(self, dem100, net100, dir100, dem1k, dir1k, basin1k, res=100, workspace='.'):
        '''

        :param dem100:
        :param net100:
        :param dir100:
        :param dem1k:
        :param dir1k:
        :param basin1k:
        :param workspace:
        '''
        self.dem100 = dem100
        self.net100 = net100
        self.dir100 = dir100
        self.dem1k = dem1k
        self.dir1k = dir1k
        self.basin1k = basin1k

        self._org_res = res
        self.res = 10 * self._org_res

        env.workspace = workspace
        env.overwriteOutput = True

        self._Slope()
        self._SlopeLength()
        self._Elevation()
        self._BedSlope()

    def _Slope(self):
        print('        Slope Calculation Beginning...')
        outSlope = Slope(self.dem100, "PERCENT_RISE")
        outBlockStat = BlockStatistics(outSlope, "Rectangle 10 10 cell", "MEAN")
        arcpy.Resample_management(outBlockStat, "Tslope", "{}".format(self.res), "NEAREST")
        outExtractByMask = ExtractByMask("Tslope", self.basin1k)
        RasterdataCheck(outExtractByMask, [0, 400])
        arcpy.RasterToASCII_conversion(outExtractByMask, "slope.asc")

        # arcpy.Delete_management('Tslope')

    def _SlopeLength(self):
        print("        Slope Length Calculation Beginning...")
        Raster1 = Raster(self.dir100)
        Raster2 = Raster(self.net100)
        rl1 = Con(Raster1 == 1, Raster2)
        rl2 = Con(Raster1 == 2, 1.4142 * Raster2)
        rl4 = Con(Raster1 == 4, Raster2)
        rl8 = Con(Raster1 == 8, 1.4142 * Raster2)
        rl16 = Con(Raster1 == 16, Raster2)
        rl32 = Con(Raster1 == 32, 1.4142 * Raster2)
        rl64 = Con(Raster1 == 64, Raster2)
        rl128 = Con(Raster1 == 128, 1.4142 * Raster2)
        arcpy.MosaicToNewRaster_management([rl1, rl2, rl4, rl8, rl16, rl32, rl64, rl128], env.workspace,
                                           "rla", "", "32_BIT_FLOAT", "", 1, "LAST", "FIRST")

        rla11 = BlockStatistics("rla", NbrRectangle(10, 10, 'CELL'), "SUM")
        rla12 = BlockStatistics("rla", NbrRectangle(100, 100, 'CELL'), "SUM") / (10 * 10)
        arcpy.MosaicToNewRaster_management([rla12, rla11], env.workspace,
                                           "rla1", "", "32_BIT_FLOAT", "", 1, "LAST", "FIRST")
        arcpy.Resample_management("rla1", "rla2", "{}".format(self.res), "NEAREST")
        slope_length1 = ExtractByMask("rla2", self.basin1k)
        slope_length2 = self.res * self.res / (2 * slope_length1 * self._org_res)
        # L = A(x)/2W(X)
        slope_length = Con(slope_length2 > self.res, self.res, slope_length2)

        RasterdataCheck(slope_length, [0, self.res])
        arcpy.RasterToASCII_conversion(slope_length, "slope_length.asc")
        arcpy.Delete_management("rla")
        arcpy.Delete_management("rla1")
        arcpy.Delete_management("rla2")

    def _Elevation(self):
        print("        Elevation Calculation Beginning...")
        elev = ExtractByMask(self.dem1k, self.basin1k)

        RasterdataCheck(elev, [-200, 9000])
        arcpy.RasterToASCII_conversion(elev, "elevation.asc")

    def _BedSlope(self):
        print("        Bed Slope Caculation Beginning...")
        rb1 = Slope(self.dem100, "PERCENT_RISE")
        rb2 = Con(Raster(self.net100) == 1, rb1)
        rb31 = BlockStatistics(rb2, "Rectangle 10 10 cell", "MEAN")
        rb32 = BlockStatistics(rb2, "Rectangle 100 100 cell", "MEAN")
        arcpy.MosaicToNewRaster_management([rb32, rb31], env.workspace,
                                           "rb3", "", "32_BIT_FLOAT", "", 1, "LAST", "FIRST")
        arcpy.Resample_management("rb3", "rb4", "{}".format(self.res), "NEAREST")
        rb4 = Raster("rb4")
        rb5 = Con(rb4 <= 0, 0.001, rb4)
        bedslope = ExtractByMask(rb5, self.basin1k)

        RasterdataCheck(bedslope, [0, 400])
        arcpy.RasterToASCII_conversion(bedslope, "bedslope.asc")
        arcpy.Delete_management('rb3')
        arcpy.Delete_management('rb4')


def RasterdataCheck(raster, range):
    '''

    :param raster:
    :param name:
    :return:
    '''

    arr = arcpy.RasterToNumPyArray(raster)
    dataCheck(arr, range, fillValue=raster.noDataValue)


class morphx:
    def __init__(self, pbasin, dir1k, bedslope, fpram, workspace='.'):
        self.home = workspace
        self.pbasin = pbasin
        self.bedslope = bedslope
        self.dir1k = dir1k
        self.fpram = fpram
        env.workspace = self.home
        env.overwriteOutput = True

        self._Split()

    def _Split(self):
        print("Splitting Parameters into Subbasin ...")
        f = open(self.fpram, "r")
        lines = f.readlines()
        f.close()
        for line in tqdm(lines):
            id, name, numsub = line.split()[0], line.split()[1], int(line.split()[2])
            path = self.home + "/subs/sub_" + name
            if not os.path.exists(path):
                os.makedirs(path)

            # --------NEXT STEP----------------------

            '''
            distance.asc
            watershed.asc
            bedslope.asc
            
            '''

            wsheds = Con(Raster(self.pbasin) == numsub, Raster(self.pbasin))
            arcpy.RasterToASCII_conversion(wsheds, path + "/watershed.asc")

            bedslope = Con(Raster(self.pbasin) == numsub, Raster(self.bedslope))
            arcpy.RasterToASCII_conversion(bedslope, path + "/bedslope.asc")

            dir = Con(Raster(self.pbasin) == numsub, Raster(self.dir1k))
            flowlen = FlowLength(dir)
            arcpy.RasterToASCII_conversion(flowlen, path + "/distance.asc")


if __name__ == "__main__":
    pass
