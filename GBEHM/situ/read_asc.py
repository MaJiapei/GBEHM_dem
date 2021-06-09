import pandas as pd
import numpy as np
from GBEHM.meteo import sta_prec_correct
import warnings


class gsod():

    def __init__(self, fname):
        self.fname = fname
        site = pd.read_csv("site_info.csv")
        f = lambda x: str(x['USAF']) + "-" + str(x['WBAN'])
        site['ID'] = site.apply(f, axis=1)

        self.site = site

    def getpik(self):

        colspecs = [(0, 6), (7, 12), (14, 22),
                    (24, 30), (31, 33), (35, 41), (42, 44),
                    (46, 52), (53, 55), (57, 63), (64, 66),
                    (68, 73), (74, 76), (78, 83), (84, 86),
                    (88, 93), (95, 100), (102, 108), (108, 109),
                    (110, 116), (116, 117), (118, 123), (123, 124),
                    (125, 130), (132, 138)]

        head = ['STN', 'WBAN', 'TIME', 'TEMP',
                'Temp_c', 'DEWP', 'Dewp_c', 'SLP',
                'Slp_c', 'STP', 'Stp_c', 'VISIB',
                'Visib_c', 'WSDP', 'Wsdp_c', 'MXSPD',
                'GUST', 'MAX', 'Max_f', 'MIN',
                'Min_f', 'PRCP', 'Prcp_f', 'SNDP',
                'FRSHTT']

        df = pd.read_fwf(self.fname, colspecs, skiprows=1,
                         names=head, parse_dates=[2], index_col=[2])
        return df

    def _correct(self, x):

        if x.isna().any():
            return x['PRCP']

        prec = x['PRCP']
        w = x['WSDP']
        tmax = x['MAX']
        tmin = x['MIN']
        td = x['TEMP']
        ID = x['STN'] + "-" + x['WBAN']
        country = self.site[self.site['ID'] == ID]['CTRY'].values[0]

        if country in ["KZ", "KG", "UZ", "TI", "TX"]:
            return sta_prec_correct.tretyakov(prec, w, tmax, tmin)
        elif country == "PK":
            return sta_prec_correct.mk2(prec, tmax, tmin)
        elif country == "BG":
            return sta_prec_correct.uswb8(prec, w, tmax, tmin)
        elif country == "IN":
            return sta_prec_correct.india(prec, w, tmax, tmin)
        elif country == "NP":
            return sta_prec_correct.nepal(prec, w, tmax, tmin)
        elif country == 'CH':
            return sta_prec_correct.cspg(prec, w, tmax, tmin, td)
        else:
            return -1
            warnings.warn("Unsupported precipitation guages.%s" % country)

    def get_prec(self, loss_corc=True):
        '''

        :param df:
        :param loss_corc:
        :return:
        '''

        df = self.getpik()
        df = df[['STN', 'WBAN', 'TEMP', 'PRCP', 'MAX', 'MIN', 'WSDP']]
        df = df.astype({'STN': str, 'WBAN': str, 'TEMP': float,
                        'PRCP': float, 'MAX': float, 'MIN': float, 'WSDP': float})

        df.replace({'TEMP': 9999.9, 'PRCP': 99.99, 'MAX': 9999.9, 'MIN': 9999.9,
                    'WSDP': 999.9}, np.nan, inplace=True)

        df['TEMP'] = (df['TEMP'] - 32) / 18
        df['MAX'] = (df['MAX'] - 32) / 18
        df['MIN'] = (df['MIN'] - 32) / 18
        df['WSDP'] = df['WSDP'] * 0.54
        df['PRCP'] = df['PRCP'] * 25.4

        if not loss_corc:
            return df[['STN', 'WBAN', 'PRCP']]
        else:
            df['PRCP'] = df.apply(self._correct, axis=1)
            return df[['STN', 'WBAN', 'PRCP']]


if __name__ == "__main__":
    obj = gsod(r"D:\data\MeteoData\world\gosd\studyarea\all_site.op")
    df = obj.get_prec()
    df.to_pickle(r"D:\data\MeteoData\world\gosd\studyarea\all_site.pik")
