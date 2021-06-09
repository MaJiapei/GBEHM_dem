import os, shutil
import numpy as np
import pandas as pd
from datetime import datetime
from GBEHM.Tool import check_within, coord_conv
from GBEHM.meteo import pyh_qual


class Mea():
    def __init__(self, df, site, id, lon, lat, elev=None,
                 varname="prec", proj4="+proj=longlat +datum=WGS84 +no_defs",
                 lack_th=10):
        '''
        :param df:
        :param site:
        :param id:
        :param lon:
        :param lat:
        '''
        self.data_df = df
        self.site_df = site

        self.lon = lon
        self.lat = lat
        self.id = id
        self.proj4 = proj4
        self.elev = elev
        self.varname = varname
        self.lack_th = lack_th

    def update_coord(self, out_proj4):

        site_df = self.site_df
        f = lambda x: [x[self.lon], x[self.lat]]
        self.site_df[[self.lon, self.lat]] = site_df.apply(f, axis=1). \
            apply(coord_conv, args=(self.proj4, out_proj4)).tolist()

        self.proj4 = out_proj4

    def value_check(self, drop_bad=True):
        df = self.data_df.copy()
        vfun = np.vectorize(pyh_qual.value_check)
        mask = vfun(np.array(df), self.varname)
        print("----------BAD VALUES SUMMARY------------------------")
        bad_row, bad_col = np.where(~mask)
        print("Find %d bad values(s)." % len(bad_row))
        print("Time                         Id            Value    ")
        for i in range(len(bad_row)):
            print("%s        %s        %f"
                  % (df.index[bad_row[i]].strftime("%Y-%m-%d %H:%M:%S"),
                     str(df.columns[bad_col[i]]),
                     df.iloc[bad_row[i], bad_col[i]]))

        if drop_bad:
            df[~mask] = np.nan
            self.data_df = df

    def select_within(self, extend_shp, inplace=False, ):
        '''

        :param extend_shp: str
            filename of shape file
        :param id: str
        :param lon: str
        :param lat: str
        :param site_info:
        :param inplace:
        :return:
        '''
        site_df = self.site_df
        lon = self.lon
        lat = self.lat
        id = self.id

        ID = self.data_df.columns
        site_df = site_df[site_df.apply(lambda x: x[id] in ID, axis=1)]
        f = lambda x: check_within([x[lon], x[lat]], self.proj4, extend_shp)
        site_df = site_df[site_df.apply(f, axis=1)]

        if inplace:
            self.data_df = self.data_df[site_df[id]]
        else:
            return self.data_df[site_df[id]]

    def to_daily_csv(self, outpath):

        if os.path.exists(outpath):
            shutil.rmtree(outpath)
            os.mkdir(outpath)
        else:
            os.mkdir(outpath)

        site = self.site_df
        df = self.data_df

        for itime in df.index:
            print('Dealing with %s...' % itime.strftime("%Y-%m-%d"))
            mdf = site.merge(df.loc[itime, :].dropna().to_frame(self.varname), right_index=True, left_on=self.id,
                             how='inner')
            oname_csv = "%s/" % outpath + itime.strftime("%Y-%m-%d") + ".csv"
            mdf = mdf.round(3)
            mdf.to_csv(oname_csv, index=False)

    def to_anomaly_csv(self, outpath):

        if os.path.exists(outpath):
            shutil.rmtree(outpath)
            os.mkdir(outpath)
            os.mkdir(outpath + "/daily_mean")
            os.mkdir(outpath + "/daily_anomaly")
        else:
            os.mkdir(outpath)
            os.mkdir(outpath + "/daily_mean")
            os.mkdir(outpath + "/daily_anomaly")

        site = self.site_df
        df = self.data_df

        davg_df = df.groupby(df.index.dayofyear).mean()
        print('Saving daily mean...')
        for itime in davg_df.index:
            mdf = site.merge(davg_df.loc[itime, :].dropna().to_frame('daily_mean'), right_index=True, left_on=self.id,
                             how='inner')
            em_array = np.array(mdf[[self.id, self.lon, self.lat, self.elev, 'daily_mean']])
            np.savetxt(outpath + "/daily_mean/" + "%03d" % (itime) + ".dat",
                       em_array, fmt="%6d%12.2f%12.2f%9.2f%7.2f")

        print('Saving daily anomaly...')
        for iindex in df.index:
            tdf = df.loc[iindex] - davg_df.loc[iindex.dayofyear]
            mdf = site.merge(tdf.dropna().to_frame('ano'), right_index=True, left_on=self.id,
                             how='inner')
            oname_csv = "%s/daily_anomaly/" % outpath + iindex.strftime("%Y-%m-%d") + ".csv"
            mdf = mdf.round(3)
            mdf.to_csv(oname_csv, index=False)

    def month_sum(self):
        site = self.site_df
        df = self.data_df
        mdf = df.resample("m").sum().to_period()
        qul_ctr = df.isna().resample('m').sum().to_period()
        mdf[qul_ctr > self.lack_th] = -999
        return mdf

    def day_month_ratio(self):
        mdf = self.month_sum().replace(-999, np.nan)
        df = self.data_df.to_period('d')
        md_df = mdf.resample('d').ffill()

        rdf = df / md_df
        rdf[md_df == 0] = 0
        return rdf

    def to_anusplin_dat(self, outpath, dem=True):

        if os.path.exists(outpath):
            shutil.rmtree(outpath)
            os.mkdir(outpath)
            os.mkdir(outpath + "/month_prec")
            os.mkdir(outpath + "/ratio")
        else:
            os.mkdir(outpath)
            os.mkdir(outpath + "/month_prec")
            os.mkdir(outpath + "/ratio")

        mdf = self.month_sum()
        rdf = self.day_month_ratio()
        site = self.site_df

        print("Saving month precipitation...")
        for imonth in mdf.index:
            tmdf = mdf.loc[imonth].replace(-999, np.nan).dropna().to_frame('month_prec')
            emdf = site.merge(tmdf, left_on=self.id, right_index=True, how='inner')

            if dem:
                em_array = np.array(emdf[[self.id, self.lon, self.lat, self.elev, 'month_prec']])
                np.savetxt(outpath + "/month_prec/" + imonth.strftime("%Y-%m") + ".dat",
                           em_array, fmt="%6d%12.2f%12.2f%9.2f%7.2f")
            else:
                em_array = np.array(emdf[[self.id, self.lon, self.lat, 'month_prec']])
                np.savetxt(outpath + "/month_prec/" + imonth.strftime("%Y-%m") + ".dat",
                           em_array, fmt="%6d%12.2f%12.2f%7.2f")

        print("Saving ratio precipitation...")
        for itime in rdf.index:
            trdf = rdf.loc[itime].dropna().to_frame('ratio')
            edf = site.merge(trdf, left_on=self.id, right_index=True, how="inner")
            oname = outpath + "/ratio/%s.csv" % itime.strftime("%Y-%m-%d")
            edf = edf.round(3)
            edf.to_csv(oname, index=False)


if __name__ == "__main__":
    df = pd.read_pickle(r"D:\code\XJ_intp\station_data\combine_sta_prec.pik")
    # site_list=[517050,518260,421310,423790,563070,560430,566410,561510,561520,420710]
    # df.drop(site_list,axis=1,inplace=True)
    site_info = r"D:\code\XJ_intp\station_data\all_sta_dem.csv"
    site_df = pd.read_csv(site_info)
    obj = Mea(df, site_df, id='SITE_NO', lon='Xpos', lat='Ypos', elev='SRTM_1km', varname='prec')
    # out_proj4=r"+proj=lcc +lat_1=30 +lat_2=35 +lat_0=30 +lon_0=87 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    # obj.update_coord(out_proj4)
    obj.to_anomaly_csv(r"D:\code\XJ_intp\ano_csv")
