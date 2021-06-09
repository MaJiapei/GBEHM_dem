'''
Author: jpma
date:2017-8-10
email:jiapeima@lzb.ac.cn
description:
    This program  focuses on extracting the subbasins by using DEM and the point of hydrological station.
'''
import datetime
import shutil
from functools import wraps

import arcpy
import numpy as np
from netCDF4 import Dataset
from tqdm import trange

from aml import *
from Fhorton import priver, morph, parameters

env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")


def ShowName(func):
    def wapper(*args, **kwargs):
        print("...Runing %s ...  " % func.__name__)
        return func(*args, **kwargs)

    return wapper


def CheckEnv(func):
    @wraps(func)
    def wapper(*args, **kwargs):
        if arcpy.env.workspace == None:
            arcpy.env.workspace = os.getcwd()
        return func(*args, **kwargs)

    return wapper


def SaveData(func):
    @wraps(func)
    def wapper(*args, **kwargs):
        if 'savePath' in kwargs:
            value = func(*args, **kwargs)
            for valuei in value:
                valuei.save(kwargs['savePath'] + "\\" + str(valuei).split('\\')[-1])
            return value
        else:
            return func(*args, **kwargs)

    return wapper


@ShowName
@CheckEnv
@SaveData
def GetRivers(DEMPath, AccThreshold, pF=""):
    env.parallelProcessingFactor = pF
    FillRaster = Fill(DEMPath)
    FlowdirRaster = FlowDirection(FillRaster, "NORMAL")
    FlowAccRaster = FlowAccumulation(FlowdirRaster, data_type="INTEGER")
    FlowAccThRaster = Con(FlowAccRaster > AccThreshold, 1)
    env.parallelProcessingFactor = ""
    return FillRaster, FlowdirRaster, FlowAccRaster, FlowAccThRaster


@ShowName
@CheckEnv
@SaveData
def ReconstructDem(DEM100Path):
    raster100m = arcpy.Raster(DEM100Path)
    res = raster100m.meanCellWidth
    ext = raster100m.extent

    Dxsize = raster100m.width % 10
    Dysize = raster100m.height % 10

    cellSizeW, cellSizeH = raster100m.meanCellWidth, raster100m.meanCellHeight

    extent = [ext.XMin, ext.YMin + Dysize * cellSizeH, ext.XMax - Dxsize * cellSizeW, ext.YMax]
    extent = " ".join(str(i) for i in extent)
    arcpy.Clip_management(raster100m, extent, 'dem1h.tif')
    dem1h = arcpy.Raster(env.workspace + '\\dem1h.tif')
    dem1k = Aggregate(dem1h, 10, "MEAN")

    return dem1h, dem1k, res


@ShowName
def _GetBasin(basinDir, shpFile, flowFill, flowAcc, flowDir, flowAccTh, flowFill100, flowDir100, flowAccTh100,
              SearchRadius):
    shpPath = basinDir + "shpPath/"

    outSnapPour = SnapPourPoint(shpFile, flowAcc, SearchRadius, "FID")
    AllWatershed = Watershed(flowDir, outSnapPour)
    arcpy.RasterToPolygon_conversion(AllWatershed,
                                     shpPath + "watershed.shp", "NO_SIMPLIFY", "VALUE")

    arcpy.Buffer_analysis(shpPath + "watershed.shp", \
                          shpPath + "watershed_buffer.shp", \
                          "4000", "FULL", "ROUND", "NONE")
    WatershedRaster = ExtractByMask(AllWatershed, shpPath + "watershed_buffer.shp")

    extentR = Con(WatershedRaster > -1, 1)

    xsize = extentR.width
    ysize = extentR.height

    # Getting the river net and classify them using STRHLER method
    accth = ExtractByMask(flowAccTh, extentR)
    outStreamOrder = StreamOrder(accth, flowDir, "STRAHLER")
    StreamToFeature(accth, flowDir, shpPath + "river.shp",
                    "NO_SIMPLIFY")

    # Getting the subbasin inside the basin and creatting the polygon
    outStreamLink = StreamLink(accth, flowDir)
    subbaisns = Watershed(flowDir, outStreamLink)
    arcpy.RasterToPolygon_conversion(subbaisns,
                                     shpPath + "subbasin.shp",
                                     "SIMPLIFY", "VALUE")

    # --------------Cut those raster to extent of the basin--------------
    AsciiPath = basinDir + "ascii/"
    os.mkdir(AsciiPath)
    FullAsciiPath = basinDir + "full_ascii/"
    os.mkdir(FullAsciiPath)
    path100 = basinDir + "dem100/"
    os.mkdir(path100)

    flowAccCut = ExtractByMask(flowAcc, extentR)
    arcpy.RasterToASCII_conversion(flowAccCut, AsciiPath + "acc1k.asc")

    FlowAccThCut = ExtractByMask(flowAccTh, extentR)
    arcpy.RasterToASCII_conversion(FlowAccThCut, AsciiPath + "net1k.asc")

    subbaisnsCut = ExtractByMask(subbaisns, extentR)
    arcpy.RasterToASCII_conversion(subbaisnsCut, AsciiPath + "wshed1k.asc")

    # dir1k.asc need boundry in Priver.f90.

    flowDirCut = ExtractByMask(flowDir, shpPath + "watershed_buffer.shp")
    arcpy.RasterToASCII_conversion(flowDirCut, AsciiPath + "dir1k.asc")

    FillRasterNoCut = ExtractByRectangle(flowFill, extentR, 'INSIDE')
    flowDirNoCut = ExtractByRectangle(flowDir, extentR, 'INSIDE')
    FillRasterNoCut.save(FullAsciiPath + "dem1k")
    flowDirNoCut.save(FullAsciiPath + "dir1k")
    extentR.save(FullAsciiPath + "basin1k")

    FillRaster100 = ExtractByRectangle(flowFill100, extentR, 'INSIDE')
    FlowDir100 = ExtractByRectangle(flowDir100, extentR, 'INSIDE')
    FlowNet100 = ExtractByRectangle(flowAccTh100, extentR, 'INSIDE')
    FillRaster100.save(path100 + "dem100")
    FlowDir100.save(path100 + "dir100")
    FlowNet100.save(path100 + "net100")

    arcpy.RasterToNetCDF_md(AsciiPath + "wshed1k.asc", "coord.nc", "basin", "m", "x", "y")

    return AsciiPath, FullAsciiPath, path100, extentR, xsize, ysize


@ShowName
def PostProcess(basinDir, AsciiPath, path100, fullAsciiPath, extentR, subDat=None, landusePath=None, **kwargs):
    # -----------------------------------------------------
    # Produce those ASCII file needed by GBHM,which
    # includes accumulation, river net, watershed,direction.
    # acc1k.asc
    # net1k.asc
    # wshed1k.asc
    # dir1k.asc   need boundry
    # ------------------------------------------------------

    # -----------------Priver.f90------------------------
    print("        Running Priver.f90...")
    _GetHortonPara(AsciiPath)
    # ---------------------------------------------------

    # ----------------deminput.aml----------------------
    print("        Running deminput.aml...")
    dem100 = path100 + "dem100"
    net100 = path100 + "net100"
    dir100 = path100 + "dir100"
    dem1k = fullAsciiPath + "dem1k"
    dir1k = fullAsciiPath + "dir1k"
    basin1k = fullAsciiPath + "basin1k"
    wk = basinDir + "deminput/"
    os.mkdir(wk)
    deminput(dem100, net100, dir100, dem1k, dir1k, basin1k, workspace=wk, **kwargs)
    # -----------------------------------------------------

    # --------------morphx.aml-----------------------------
    print("Running morphx.aml...")
    pbasin = AsciiPath + "para/pbasin"
    dir1k = fullAsciiPath + "dir1k"
    bedslope = basinDir + "deminput/bedslope.asc"
    fpram = AsciiPath + "para/subcatchment.txt"
    wk = basinDir + "morphx/"
    os.mkdir(wk)
    morphx(pbasin, dir1k, bedslope, fpram, workspace=wk)
    shutil.copy(AsciiPath + "para/subcatchment.txt", wk)
    shutil.copy(AsciiPath + "para/subbasin.dat", wk)
    # -----------------------------------------------------

    # --------------morph.f&parameter.f------------------------
    if subDat:
        print('        Running morph.f and parameter.f... ')
        shutil.copy(subDat, wk)
        morph_para(wk)

    # ----------------------GLT Build---------------------------
    print('        Building Geographic Lookup Table...')
    os.mkdir(basinDir + "GLTMask")
    GLTBuild(basinDir + "ascii/wshed1k.asc", basinDir + "GLTMask/", MaskFile=True)

    # -----------------------land use----------------------------
    if landusePath:
        print('...Extracting landuse...')
        arcpy.env.snapRaster = extentR
        landuse = ExtractByMask(landusePath, extentR)
        arcpy.env.snapRaster = None
        arcpy.RasterToASCII_conversion(landuse, AsciiPath + "landuse.asc")


@ShowName
def SaveNetCDF(xsize, ysize, basinDir, AsciiPath, workspace):
    # ----------------------saving as netCDF format----------------------
    # bedslope.asc, slope.asc, slope_length.asc
    # landuse_lucc.asc, soil_depth,watershed_lon.asc,watershed_lat.asc
    # watershed_mask.asc

    nf = Dataset(workspace + "/GisPara.nc", "w")
    nf.createDimension("x", xsize)
    nf.createDimension("y", ysize)

    sp = nf.createVariable("slope", "f4", ("y", "x"), fill_value=-9999)
    splength = nf.createVariable("slope_length", "f4", ("y", "x"), fill_value=-9999)

    elev = nf.createVariable("elevation", "f4", ("y", "x"), fill_value=-9999)

    lon = nf.createVariable("longitude", 'f4', ('y', 'x'))
    lat = nf.createVariable("latitude", 'f4', ('y', 'x'))
    mk = nf.createVariable("mask", 'f4', ('y', 'x'))

    sp[:] = np.flipud(arcpy.RasterToNumPyArray(basinDir + "deminput/" + "slope.asc"))
    splength[:] = np.flipud(arcpy.RasterToNumPyArray(basinDir + "deminput/" + "slope_length.asc"))
    elev[:] = np.flipud(arcpy.RasterToNumPyArray(basinDir + "deminput/" + "elevation.asc"))

    if os.path.exists(AsciiPath + "landuse.asc"):
        ld = nf.createVariable("landuse", "f4", ("y", "x"), fill_value=-9999)
        ld[:] = np.flipud(arcpy.RasterToNumPyArray(AsciiPath + "landuse.asc")).astype('int')
        ld.units = ''
        ld.longname = 'landuse using classification standards of USGS'

    lon[:] = np.flipud(arcpy.RasterToNumPyArray(basinDir + "GLTMask/" + "wshed1k.asc_lon.asc"))
    lat[:] = np.flipud(arcpy.RasterToNumPyArray(basinDir + "GLTMask/" + "wshed1k.asc_lat.asc"))
    mk[:] = np.flipud(arcpy.RasterToNumPyArray(basinDir + "GLTMask/" + "wshed1k.asc_mask.asc"))

    sp.units = 'percent'
    sp.longname = 'slope in gradient'
    sp.scale_factor = 0.01
    sp.add_offset = 0.

    splength.units = 'm'
    splength.longname = "slope length"

    elev.units = "m"
    elev.longname = "elevation"

    lon.units = "degrees_east"
    lon.longname = "longitude"

    lat.units = "degrees_north"
    lat.longname = "latitude"

    mk.units = ''
    mk.longname = 'mask'

    nf.creattime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    nf.author = 'jpma'
    nf.program = "jpma's subbasin extracting program."
    nf.close()


def Main(DEM100Path, StationPath, WorkSpace, SearchRadius, AccThreshold, landusePath=None, subdat=None):
    env.workspace = WorkSpace

    dem1h, dem1k, res = ReconstructDem(DEM100Path)

    FillRaster100, FlowdirRaster100, FlowAccRaster100, FlowAccThRaster100 = \
        GetRivers(dem1h, AccThreshold, pF='8')

    flowFill, flowDir, flowAcc, flowAccTh = \
        GetRivers(dem1k, AccThreshold)

    with arcpy.da.SearchCursor(StationPath, ["FID", "SHAPE@"]) as cursor:
        for row in cursor:
            # disasssemble the mutil points shpfile into single point shpfile
            basinDir = WorkSpace + "/station" + str(row[0]) + "/"
            shpPath = basinDir + "shpPath/"
            os.makedirs(basinDir)
            os.makedirs(shpPath)
            shpFile = shpPath + str(row[0]) + ".shp"
            arcpy.CopyFeatures_management(row[1], shpFile)

            # Using the single point to define the basin and get the polygon
            # of the basin

            AsciiPath, fullAsciiPath, path100, extentR, xsize, ysize = \
                _GetBasin(basinDir, shpFile, flowFill, flowAcc,
                          flowDir, flowAccTh, FillRaster100, FlowdirRaster100, FlowAccThRaster100, SearchRadius)
            PostProcess(basinDir, AsciiPath, path100, fullAsciiPath, extentR,
                        res=res, subDat=subdat, landusePath=landusePath)
            SaveNetCDF(xsize, ysize, basinDir, AsciiPath, WorkSpace)


def _GetHortonPara(AsciiPath):
    obj = priver
    # ---------------------Notice----------------------
    # This step integrates the Priver.f90 as a
    # python module. It should be noted that the two
    # parameters 'gis_dir', 'para_dir' are
    # strings both having 100 characters, thus the rest
    # characters should be replaced by blank.
    # -------------------------------------------------

    obj.gis_dir = AsciiPath.ljust(100)
    os.mkdir(AsciiPath + "para/")
    obj.para_dir = (AsciiPath + "para/").ljust(100)
    obj.getpara()

    os.chdir(AsciiPath + "para\\")
    # If don't change the pwd path, there will be a bug.
    arcpy.ASCIIToRaster_conversion("pbasin.asc",
                                   "./pbasin", "INTEGER")


def morph_para(model_dir):
    if not os.path.exists(model_dir + "subcatchment.dat"):
        raise IOError("subcatchment.dat doesn't exist.")
    dem_home = model_dir + "subs/".ljust(200)
    nsub = len(os.listdir(model_dir + "subs/"))
    fileList = os.listdir(model_dir + "subs/")
    temp = Raster(model_dir + "subs/" + fileList[0] + "/bedslope.asc")
    nc, nr = temp.width, temp.height
    print("nsub:%d   nc:%d   nr%d" % (nsub, nc, nr))

    model_dir = model_dir.ljust(200)

    # ------------Run morph_horton.f--------------
    morph(nsub, dem_home, model_dir)

    # ------------Run parameter_horton.f-------------
    obj = parameters
    obj.nc = nc
    obj.nr = nr
    obj.nsub = nsub
    obj.home = model_dir
    obj.para()


@CheckEnv
def GLTBuild(Image, Outpath, MaskFile=False):
    '''
        Using a raster file to get the GLT file.
    :param Image: str
        path of the file used to creat the GLT file.
    :param Outpath: str
        path to save the GLT file.
    :return:
    '''
    raster = arcpy.Raster(Image)
    spRef = raster.spatialReference
    ext = raster.extent
    ll = arcpy.Point(ext.XMin, ext.YMin)

    ResInfo = [raster.extent.upperLeft.X, raster.meanCellWidth,
               raster.extent.upperLeft.Y, raster.meanCellHeight]

    xsize = raster.width
    ysize = raster.height

    print("xsize: %f" % (xsize))
    print("ysize: %f" % (ysize))

    LonLatGLT = np.zeros((2, ysize, xsize))

    for j in trange(ysize, ncols=75):
        for i in range(xsize):
            xCoord = ResInfo[0] + (i + 0.5) * ResInfo[1]
            yCoord = ResInfo[2] - (j + 0.5) * ResInfo[3]
            LonLat = _ProTrans(spRef, [xCoord, yCoord])
            LonLatGLT[0, j, i] = LonLat[0]
            LonLatGLT[1, j, i] = LonLat[1]

    LonGLT = arcpy.NumPyArrayToRaster(LonLatGLT[0], ll, raster.meanCellWidth, raster.meanCellHeight)
    LatGLT = arcpy.NumPyArrayToRaster(LonLatGLT[1], ll, raster.meanCellWidth, raster.meanCellHeight)

    arcpy.DefineProjection_management(LonGLT, spRef)
    arcpy.DefineProjection_management(LatGLT, spRef)

    print('Saving longitude file...')
    arcpy.RasterToASCII_conversion(LonGLT, Outpath + "/" + Image.split('/')[-1] + "_LON.asc")
    print('Saving latitude file...')
    arcpy.RasterToASCII_conversion(LatGLT, Outpath + "/" + Image.split('/')[-1] + "_LAT.asc")

    if MaskFile:
        print('Saving mask file...')
        MASK = arcpy.RasterToNumPyArray(Image)
        MASK[MASK != raster.noDataValue] = 1
        MASK[MASK == raster.noDataValue] = 0
        MASK = arcpy.NumPyArrayToRaster(MASK, ll, raster.meanCellWidth, raster.meanCellHeight)
        arcpy.DefineProjection_management(MASK, spRef)
        arcpy.RasterToASCII_conversion(MASK, Outpath + "/" + Image.split('/')[-1] + "_MASK.asc")


def _ProTrans(PrjInfo, XY, TargetCoord=4326):
    '''
        coordinate convert
    :param PrjInfo: arcpy geoprocessing spatial reference
        source spatialreference coming from a raster
    :param XY:list or array-like
        point that will be converted to target coordinate.
    :param TargetCoord:int
        coordinate system's factory code
    :return:
    '''
    point = arcpy.Point()
    point.X = XY[0]
    point.Y = XY[1]

    TargetSRS = arcpy.SpatialReference(TargetCoord)
    ptgeo = arcpy.PointGeometry(point, PrjInfo)
    ptgeo = ptgeo.projectAs(TargetSRS)

    return [ptgeo.firstPoint.X, ptgeo.firstPoint.Y]


if __name__ == "__main__":
    DEM100Path = r"G:\GBEHM_training_materials\data\data\dem15.tif"
    StationPath = r"G:\GBEHM_training_materials\data\data\hy_station.shp"
    subdata = r"G:\GBEHM_training_materials\data\data\subcatchment.dat"
    SearchRadius = 1000
    AccThreshold = 100
    if os.path.exists(r"D:\wsp2"):
        shutil.rmtree(r"D:\wsp2")
        os.mkdir(r"D:\wsp2")
    else:
        os.mkdir(r"D:\wsp2")

    workspace = r"D:\wsp2"

    Landuse = r"G:\GBEHM_training_materials\data\data\landuse50.tif"
    # raster100marr = arcpy.RasterToNumPyArray(DEM100Path)
    Main(DEM100Path, StationPath, workspace, SearchRadius, AccThreshold, landusePath=Landuse, subdat=subdata)
