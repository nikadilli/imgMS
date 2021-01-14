import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
from decimal import Decimal
import datetime
import logging

from imgMS.side_functions import *


class DataReader():
    """
    Reads data into pandas Dataframe from a file.

    Parameters
    ----------
    filetype : str
        Type of the file to read. If not specified, csv is used. 
        Possible options are csv, xlsx and asc.
    instrument : str
        Type of the instrument used for measurement. If not specified, raw data is expected. 
        Possible options are Agilent, and Element.

    Returns
    -------
    data : dataframe
        `data` as a dataframe, which can be passed to MSData.

    """

    def __init__(self, filename, filetype=None, instrument=None):

        self.filename = filename
        self.filetype = filetype
        self.instrument = instrument
        self.data = self.read(filename, filetype, instrument)

    def __call__(self):
        return self.data

    def read(self, filename, filetype, instrument):
        if instrument == 'Element':
            skipfooter = 4
            header = 1
            drop = 9
        elif instrument == 'Agilent':
            skipfooter = 4
            header = 3
            drop = 3
        else:
            skipfooter = 0
            header = 0
            drop = 0

        if filetype == 'xlsx':
            imported = pd.ExcelFile(filename)
            data = imported.parse(
                0, index_col=0, skipfooter=skipfooter, header=header)
            data = data.drop(data.index[:drop], axis=0)

        elif filetype == 'csv':
            data = pd.read_csv(filename, sep=',', index_col=0, skipfooter=skipfooter,
                               header=header, engine='python')

        elif filetype == 'asc':
            data = pd.read_csv(filename, sep='\t', index_col=0, skipfooter=skipfooter,
                               header=header, engine='python')
            data = data.drop(data.index[:drop], axis=0)
            data.dropna(axis=1, how='all', inplace=True)
            data = data.apply(pd.to_numeric, errors='coerce')

        else:
            warnings.warn('File type not supported.')

        return data


class Background():
    def __init__(self, isotope, laser_on, laser_off, width=0.8, offset=0.15):
        assert isotope.data is not None
        self.bcg_all = [remove_outliers(
            isotope.data[laser_off[i][0]:laser_off[i][1]], offset, width) for i in range(len(laser_off))]
        self.bcg_means = [x.mean() for x in self.bcg_all]

        bcg_mskd = np.copy(isotope.data)
        for (s, e) in laser_on:
            bcg_mskd[s-2:e+2] = np.nan
        not_nan = np.logical_not(np.isnan(bcg_mskd))
        indices = np.arange(len(bcg_mskd))
        self.bcg_interp = np.interp(
            indices, indices[not_nan], bcg_mskd[not_nan])

    def __call__(self):
        plt.plot(self.bcg_means)

    def __repr__(self):
        return f'{self.__class__.__name__}: {self.bcg_means} '


class Iolite():
    """
    Class holding Iolite data.

    Parameters
    ----------
    path : str
        Path to Iolite .csv file.

    """

    def __init__(self, path):
        self.data = self.read_iolite(path)
        self.peak_names = self.names_from_iolite()

    def read_iolite(self, path):
        iolite = pd.read_csv(path, sep=",", engine='python')
        return iolite

    def names_from_iolite(self):
        names = list(self.data[' Comment'].dropna())
        return names


class Lase():
    pass


class Param():
    """
    Class holding meta data and parameters of analysis.

    Parameters
    ----------
    path : str
        Path to param .xlsx file.

    """

    def __init__(self, path, logger=None):
        self.path = path
        self.peak_names = self.read_names()
        self.is_coef = self.read_is_coef()
        self.ts_coef = self.read_ts_coef()
        self.logger = logger

    def read_names(self):
        xl = pd.ExcelFile(self.path)
        if 'names' in xl.sheet_names:
            names = list(xl.parse('names', header=None)[0])
            return names
        else:
            if self.logger is not None:
                self.logger.error(
                    'Can not read names of peak from Param file.')

    def read_is_coef(self):
        xl = pd.ExcelFile(self.path)
        if 'internal standard' in xl.sheet_names:
            internal_std = xl.parse(
                'internal standard', index_col=0, header=0)
            return internal_std
        else:
            if self.logger is not None:
                self.logger.error(
                    'Can not read coefficients of internal standards from Param file.')

    def read_ts_coef(self):
        xl = pd.ExcelFile(self.path)
        if 'total sum' in xl.sheet_names:
            sum_koeficients = xl.parse(
                'total sum', index_col=0, header=None).to_dict()[1]
            return sum_koeficients
        else:
            if self.logger is not None:
                self.logger.error(
                    'Can not read coefficients for total sum from Param file.')


class Selector():
    def __init__(self, ms_data, s=60, sdmul=10, iolite=None, logger=None):
        self.possible_methods = ('treshold', 'iolite')
        self.ms_data = ms_data
        self.logger = logger
        self.filter_line = self.ms_data.data.sum(1)
        self.method = 'treshold'
        self.iolite = iolite
        self.start = self.ms_data.data.index[get_index(self.ms_data.data, s)]
        if self.logger is not None:
            self.logger.info(f'Starting point found at {self.start}s.')
        self.sdmul = sdmul
        self.skip = {'bcg_start': 0,
                     'bcg_end': 0,
                     'sample_start': 0,
                     'sample_end': 0}    # time in seconds to skip from each bcg and sample

    def __call__(self):
        if self.method == 'treshold':
            return self.create_selector_treshold()
        elif self.method == 'iolite':
            return self.create_selector_iolite()

    def create_selector_iolite(self):
        """select starts and ends of ablation using iolite file """
        if self.logger is not None:
            self.logger.info('Selecting peak bounds by iolite.')

        lst = [x for x in self.iolite.data.loc[:6,
                                               ' Comment'] if isinstance(x, str)]

        if len(lst) == 2:
            if self.logger is not None:
                self.logger.info('Selecting spots.')
            difflst = get_diff_lst(self.iolite.data)
        elif len(lst) == 1:
            if self.logger is not None:
                self.logger.info('Selecting lines.')
            difflst = get_diff_lst_line(self.iolite.data)
        else:
            if self.logger is not None:
                self.logger.error('Iolite selection failed.')

        timeindex = []
        for i in range(0, len(difflst)+1):
            timeindex.append(sum(difflst[:i])+self.start)
        index = [get_index(self.ms_data.data, x) for x in timeindex]

        starts = [index[i] for i in range(len(index)) if i % 2 == 0]
        ends = [index[i]+1 for i in range(len(index)) if i % 2 != 0]

        return starts, ends

    def create_selector_treshold(self):
        """
        select starts and ends of ablation based on selected element or sum of all using treshold
        calculated from background
        """
        if self.logger is not None:
            self.logger.info('Selecting peak bounds by setting treshold.')
        bcg_nr = self.ms_data.time_to_number(self.start)
        bcg = self.filter_line[0:bcg_nr].mean()
        std = self.filter_line[0:bcg_nr].std()
        ind = [True if value > bcg+self.sdmul *
               std else False for value in self.filter_line]
        ind2 = ind[1:]
        ind2.append(False)
        index = [i for i in range(0, len(ind)) if ind[i] != ind2[i]]

        starts = [index[i] for i in range(len(index)) if i % 2 == 0]
        ends = [index[i] for i in range(len(index)) if i % 2 != 0]

        return starts, ends

    def set_skip(self, bcg_s=None, bcg_e=None, sig_s=None, sig_e=None):
        """
        set time skipped on start and end of background and ablation in seconds
        """
        if bcg_s is not None:
            self.skip['bcg_start'] = bcg_s
        if bcg_e is not None:
            self.skip['bcg_end'] = bcg_e
        if sig_s is not None:
            self.skip['sample_start'] = sig_s
        if sig_e is not None:
            self.skip['sample_end'] = sig_e

    def create_on_off(self, starts, ends):
        """
        from starts and ends of ablation create laser_on and laser_off with skipped values
        """
        laser_off = []
        laser_on = []

        laser_off.append((0+self.ms_data.time_to_number(
            self.skip['bcg_start']), starts[0]-self.ms_data.time_to_number(self.skip['bcg_end'])))

        for i in range(len(starts)-1):
            laser_off.append((ends[i]+self.ms_data.time_to_number(self.skip['bcg_start']),
                              starts[i+1]-self.ms_data.time_to_number(self.skip['bcg_end'])))
            laser_on.append((starts[i]+self.ms_data.time_to_number(self.skip['sample_start']),
                             ends[i]-self.ms_data.time_to_number(self.skip['sample_end'])))

        laser_off.append((ends[-1]+self.ms_data.time_to_number(self.skip['bcg_start']), len(
            self.ms_data.time)-2-self.ms_data.time_to_number(self.skip['bcg_end'])))
        laser_on.append((starts[-1]+self.ms_data.time_to_number(self.skip['sample_start']),
                         ends[-1]-self.ms_data.time_to_number(self.skip['sample_end'])))
        return laser_on, laser_off
