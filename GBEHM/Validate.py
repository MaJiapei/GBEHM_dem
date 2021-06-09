import numpy as np
import pandas as pd

'''
This program is used to validate the preciptation between simulation and stations.
'''

__rule = [[0.1, 4.9], [5, 9.9], [10, 14.9], [15., 24.9], [25., 49.9], [50., 99.9], [100., 250.]]
__col = ['0.1-4.9', '5.-9.9', '10.-14.9', '15.-24.9', '25.-49.9', '50.-99.9', '100.-250.']


class QualityTable:

    def __BasciSTA(self, mea, simu):
        t_mea, m_mea, s_mea, a_mea = mea.sum(), mea.max(), mea.min(), mea.mean()
        t_simu, m_simu, s_simu, a_simu = simu.sum(), simu.max(), simu.min(), simu.mean()
        Binfo = pd.DataFrame(columns=[mea.name, simu.name], index=['sum', 'max', 'min', 'mean'])
        Binfo.ix[:, 1] = [t_simu, m_simu, s_simu, a_simu]
        Binfo.ix[:, 0] = [t_mea, m_mea, s_mea, a_mea]
        return Binfo

    def __CapRate(self, mea, simu):

        nday = float(mea.__len__())

        Cinfo0 = pd.DataFrame(columns=['N', 'Y'], index=[mea.name, 'precnt', simu.name, 'precnt'])
        Cinfo1 = pd.DataFrame(columns=['--', 'NA', 'NC', 'NB'], index=['times', 'precnt'])
        Cinfo2 = pd.Series(index=['TS', 'PO', 'FAR'])

        Tmea_Y = mea[mea != 0].__len__()
        Tmea_N = mea[mea == 0].__len__()
        Tsim_Y = simu[simu != 0].__len__()
        Tsim_N = simu[simu == 0].__len__()

        Cinfo0['N'] = Tmea_N, float(Tmea_N) / (Tmea_N + Tmea_Y), Tsim_N, float(Tsim_N) / (Tsim_N + Tsim_Y)
        Cinfo0['Y'] = Tmea_Y, float(Tmea_Y) / (Tmea_N + Tmea_Y), Tsim_Y, float(Tsim_Y) / (Tsim_N + Tsim_Y)

        # ----------------------------------------------------------------------------------------------------

        p1 = pd.DataFrame(mea[mea == 0]). \
            join(pd.DataFrame(simu[simu == 0]), how='inner').__len__()
        p2 = pd.DataFrame(mea[mea != 0]). \
            join(pd.DataFrame(simu[simu != 0]), how='inner').__len__()
        p3 = pd.DataFrame(mea[mea != 0]). \
            join(pd.DataFrame(simu[simu == 0]), how='inner').__len__()
        p4 = pd.DataFrame(mea[mea == 0]). \
            join(pd.DataFrame(simu[simu != 0]), how='inner').__len__()

        Cinfo1.ix['times'] = p1, p2, p3, p4
        Cinfo1.ix['precnt'] = p1 / nday, p2 / nday, p3 / nday, p4 / nday

        Cinfo2['TS'] = 100 * float(p2) / (p2 + p3 + p4)
        Cinfo2['PO'] = 100 * float(p3) / (p2 + p3)
        Cinfo2['FAR'] = 100 * float(p4) / (p2 + p4)

        return Cinfo0, Cinfo1, Cinfo2

    def __FreqVal(self, mea, simu):
        Finfo2 = pd.DataFrame(columns=__col, index=['TS', 'PO', 'FAR'])
        for i, ru in enumerate(__rule):
            try:
                a = mea[(mea > ru[0]) & (mea < ru[1])]
                b = simu.ix[a.index]
                Cinfo0, Cinfo1, Cinfo2 = __CapRate(a, b)
                Finfo2.iloc[:, i] = Cinfo2
            except:
                pass
        return Finfo2

    def __CoorRMSE(self, mea, simu):
        r1 = mea.corr(simu)
        a = pd.DataFrame(mea[mea != 0])
        b = pd.DataFrame(simu[simu != 0])
        c = a.join(b, how='inner')
        r2 = c.ix[:, 0].corr(c.ix[:, 1])

        rmse = ((c.ix[:, 0] - c.ix[:, 1]) ** 2).mean()
        rmse = np.sqrt(rmse)

        CRinfo = pd.Series(index=['r1', 'r2', 'rmse'])
        CRinfo[:] = r1, r2, rmse

        return CRinfo

    def __TimesInDifRange(self, mea, simu):
        Tr = pd.DataFrame(index=[mea.name, simu.name], columns=__col)
        for i, ru in enumerate(__rule):
            try:
                a = mea[(mea > ru[0]) & (mea < ru[1])].__len__()
                b = simu[(simu > ru[0]) & (simu < ru[1])].__len__()
                Tr.ix[:, i] = a, b
            except:
                pass
        return Tr

    def show(self, mea, simu):
        if not (isinstance(mea, pd.core.series.Series) and \
                isinstance(simu, pd.core.series.Series)):
            raise TypeError('Wrong Input type,must be pandas.core.series.Series.')

        print('#-------------Basic Infomation-----------------------#')
        print(self.__BasciSTA(mea, simu))
        print("-----------------------------------------------------#")
        print('\n')

        print('#-------------Correlation and RMSE-------------------#')
        print(self.__CoorRMSE(mea, simu))
        print('#----------------------------------------------------#')
        print('\n')

        print('#-------------Capture Proportion---------------------#')
        Cinfo0, Cinfo1, Cinfo2 = self.__CapRate(mea, simu)
        print('Dry and Wet:')
        print(Cinfo0)
        print('\n')
        print('Catch Prop:')
        print(Cinfo1)
        print('\n')
        print('TS Marks in All Ranges:')
        print(Cinfo2)
        print('#----------------------------------------------------#')
        print('\n')

        print('#------------------------TS in Different Range---------------------------------#')
        print(self.__FreqVal(mea, simu))
        print("-------------------------------------------------------------------------------#")
        print('\n')

        print('#------------------------Times in Different Range-------------------------------#')
        print(self.__TimesInDifRange(mea, simu))
        print("--------------------------------------------------------------------------------#")


def NSE(mea, simu):
    '''

    :param mea:
    :param simu:
    :return:
    '''

    return 1 - ((mea - simu) ** 2).sum() / \
           ((mea - mea.mean()) ** 2).sum()


def PBias(mea, simu):
    '''

    :param mea:
    :param simu:
    :return:
    '''
    a1 = (mea - simu).sum()
    a2 = mea.sum()
    value = 100 * a1 / a2
    return value


def RSR(mea, simu):
    return np.sqrt(((mea - simu) ** 2).mean()) / mea.std(ddof=0)


if __name__ == "__main__":
    pass
