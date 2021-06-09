import numpy as np
import pandas as pd
import warnings


def byTime(pikData, freq='M', outPath='./', nrc=(), format='npy', prefix='', suffix='', ):
    '''
    Splitting a dataFrame to a specified period.

    :param pikData: dataFrame created by pandas of pickle file.
    :param freq: frequency the out file will use.
    :param outPath: path where the out file will be saved.
    :param nrc: tuple contains rows and columns of the array that dataFrame will be saved.
    :param format: format that will be saved.
    :param prefix: prefix that will be added to the front the out name.
    :param suffix: suffix that will be added to the tail of the out name.
    :return:
    '''
    if isinstance(pikData, str):
        df = pd.read_pickle(pikData)
    else:
        df = pikData

    try:
        Mtag = df.to_period(freq).index.unique()
    except:
        raise ValueError('Error in creating time tag.')

    if format == 'npy':
        if nrc:
            for iMtag in list(map(str, Mtag)):
                print('Dealing with %s' % iMtag)
                oName = outPath + prefix + iMtag + suffix + '.npy'
                tdf = df.loc[iMtag]
                tarray = np.array(tdf).reshape(tdf.shape[0], nrc[0], nrc[1])
                np.save(oName, tarray)
        else:
            for iMtag in list(map(str, Mtag)):
                print('Dealing with %s' % iMtag)
                oName = outPath + prefix + iMtag + suffix + '.npy'
                tarray = np.array(df.loc[iMtag])
                np.save(oName, tarray)
    elif format == 'pik':
        for iMtag in list(map(str, Mtag)):
            print('Dealing with %s' % iMtag)
            oName = outPath + prefix + iMtag + suffix + '.pik'
            df.loc[iMtag].to.pickle(oName)
    else:
        raise ValueError('Unknown format ' + format)


def dataCheck(arr, range, fillValue=None):
    '''
        check if the arr contains invalid value and drop those value if drop=True
    :param arr:numpy.ndarray
        data array
    :param range:list
        two elements list like [0,10],first value is the minimum
        second value is maximum
    :param fillValue:none,int,float
        fillvalue in the data array
    :return:
    '''

    minValue = range[0]
    maxValue = range[1]

    arr = arr.astype(np.float)

    if fillValue != None:
        arr[arr == fillValue] = np.nan
        if (np.nanmax(arr) > maxValue) | (np.nanmin(arr) < minValue):
            warnings.warn("data is beyond the valid range.")

    else:
        if np.isnan(arr).any():
            if (np.nanmax(arr) > maxValue) | (np.nanmin(arr) < minValue):
                warnings.warn("data is beyond the valid range.")
        else:
            if (arr.max() > maxValue) | (arr.min() < minValue):
                warnings.warn("data is beyond the valid range.")


if __name__ == "__main__":
    print(byTime.__doc__)
