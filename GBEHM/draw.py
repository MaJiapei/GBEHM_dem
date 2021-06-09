import pandas as pd
from matplotlib import pyplot as pt
import numpy as np
from GBHM.Validate import NSE, PBias


def _draw(SimRunoff, **kwargs):
    fig = pt.figure()
    pt.plot(SimRunoff.index, SimRunoff.values, label='Simulated')
    pt.ylabel(r'Streamflow' + '$(m^3/s)$')
    pt.xlabel(r'Time')
    if 'precp' in kwargs:
        precp = kwargs['precp']
        precp = precp.resample('D').sum().loc[SimRunoff.index]
        ax = pt.gca()
        ax1 = ax.twinx()
        ax1.invert_yaxis()
        ax1.set_ylim(50, 0)
        ax1.bar(precp.index, precp.values, color='red', label='Precipitation')
        ax1.set_ylabel('Daily Precipitation(mm/d)')

    if 'valRunoff' in kwargs:
        valrunoff = kwargs['valRunoff']
        valrunoff = valrunoff.loc[SimRunoff.index]
        ax = pt.gca()
        ax.plot(valrunoff.index, valrunoff.values, label='Measurement')

        nash = NSE(SimRunoff.values, valrunoff.values)
        pbias = PBias(SimRunoff.values, valrunoff.values)
        ax.text(0.8, 0.8, 'NSE=%5.2f' % nash, transform=ax.transAxes)
        ax.text(0.8, 0.9, 'PBIAS=%5.2f' % pbias, transform=ax.transAxes)
    fig.legend(loc=2, ncol=1 + kwargs.__len__())

    if kwargs['save']:
        pt.savefig(kwargs['save'])
    else:
        pt.show()


def draw_runoff(simu_runoff, val_runoff=None, prec=None, save=None):
    df_simu = pd.read_table(simu_runoff, sep='\s+', names=['year', 'month', 'day', 'simu'], \
                            parse_dates={'time': [0, 1, 2]}, \
                            keep_date_col=False, dtype={'simu': np.float})
    df_simu.set_index('time', inplace=True, drop=True)
    if val_runoff:
        df_mea = pd.read_csv(val_runoff, parse_dates={'time': [0, 1, 2]}, keep_date_col=False)
        df_mea.set_index('time', inplace=True)
        df_mea.columns = ['site']
        _draw(df_simu, valRunoff=df_mea, save=save)
    else:
        _draw(df_simu, save=save)


if __name__ == "__main__":
    pass
