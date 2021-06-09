import numpy as np
import pandas as pd

from netCDF4 import Dataset
from GBEHM.Extract import byRaster
from GBEHM.Tool import GetGLT




def etorh(qe, airtk, pres):
    '''
    This code comes from forcing_class.py by Li Hongyi,it is used
    to convert the specific humidity to relative humidity

    :param self:
    :param qe:
    :param airtk:
    :param pres:
    :return:
    '''
    airtc = airtk - 273.15
    es = 611. * np.exp(17.27 * airtc / (237.3 + airtc))
    es1 = 611. * np.exp(21.87 * airtc / (265.5 + airtc))
    # es[airtc < 0.] = es1[airtc < 0.]cdo

    rh = qe * pres / (0.622 * es)
    qe = np.where(rh > 1., 0.622 * es * 0.95 / pres, qe)
    rh = np.where(rh > 1., 0.95, rh)
    qe = np.where(rh < 0., 0.622 * es * 0.05 / pres, qe)
    rh = np.where(rh < 0., 0.05, rh)

    return rh, qe



def wrfDriver(fname,GLTFile):
    pass

def CommonDriver(fname,GLTFile,varName,Opath):

    '''
    :param fname:
    :param GLTFile:
    :param varName:
    :return:
    '''
    Flon,Flat=GLTFile[0], GLTFile[1]
    var=byRaster(fname,Flon,Flat,varName,method="nearest")


def XJDriver(Ipath,Oname):
    '''
    
    :return: 
    '''
    '''
    Unit 
    10MWS10                  m/s
    2mtemperature            C
    PSFC                     hpa
    2mrh2                    %
    DLWRF_P8_L1_GLL0_avg     W/m^2
    DSWRF_P8_L1_GLL0_avg     W/m^2
    RAINC                    mm/h
    
    '''
    index=pd.date_range('20170910','20180330',freq='h')


def _unitCovert():
    '''
cond
    :return:
    '''
    #-----------------UNIT IN GBHM--------------------
    # wind at 10m:                     U10 m/s V10 m/s
    # Temperature at 2m:               T2 K
    # Pressure at surface:             P0 Pa
    # Relative humidity:               RH %
    # Specific humidity:               Q  Kg/Kg
    # Downward long wave radiance:     w/m^2
    # Downward short wave radiance:    w/m^2
    # Precipitation:                   mm/s
    #-------------------------------------------------


    pass



def ForcingBadValueRemove(var,name):
    '''

    :param var:
    :param name:
    :return:
    '''
    RangeTable={'U10':[-100,100],'V10':[-100,100],'T2':[180,340],\
                'P0':[1000,101325],'RH':[0,1],\
                'LR':[10,1000],'SR':[0,1500],'PREC':[0,0.028]}

    if var.min() < RangeTable[name][0] or var.max > RangeTable[name][1]:
        print('-------------VALUE REPORT---------------------')
        print('      Bad Value in %s Found!!!      '%name)
        print(' Max Value: %f,Min Value: %f' % (var.max, var.min))
        print(' Excuting Process of Bad Value Remove.')

        var[var < RangeTable[name][0]] = RangeTable[name][0]
        var[var > RangeTable[name][1]] = RangeTable[name][1]

    else:
        print('----------VALUE REPORT------------')
        print('      Bad Value in %s not Found.      '%name)

    return var



if __name__=="__main__":
    fname=r"D:\data\XJsimu\drivers\XJmeteo_2mtemperature_HOURS.nc"
    Flon=r"../../data/watershed_lon.asc"
    Flat=r"../../data/watershed_lat.asc"