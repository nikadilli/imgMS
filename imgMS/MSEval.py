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
    filename: str
        Path to file to read.
    filetype : str
        Type of the file to read. If not specified, csv is used. 
        Possible options are csv, xlsx and asc.
    instrument : str
        Type of the instrument used for measurement. If not specified, raw data is expected. 
        Possible options are Agilent, and Element.
    """

    def __init__(self, filename, filetype=None, instrument=None):

        self.filename = filename
        self.filetype = filetype
        self.instrument = instrument
        self.data = self.read(filename, filetype, instrument)

    def __call__(self):
        return self.data

    def read(self, filename, filetype, instrument):
        """
        Reads data into pandas Dataframe from a file.

        Parameters
        ----------
        filename: str
            Path to file to read.
        filetype : str
            Type of the file to read. If not specified, csv is used. 
            Possible options are csv, xlsx and asc.
        instrument : str
            Type of the instrument used for measurement. If not specified, raw data is expected. 
            Possible options are Agilent, MC Nu-Sapphire, and Element.

        Returns
        -------
        data : dataframe
            `data` as a dataframe, which can be passed to MSData.
        """
        cols_to_drop = []

        if instrument == 'Element':
            skipfooter = 4
            header = 1
            drop = 9
            index_col = 0
        elif instrument == 'Agilent':
            skipfooter = 4
            header = 3
            drop = 3
            index_col = 0
        elif instrument =='MC Nu-Sapphire':
            skipfooter = 0
            header = 14
            drop = 0
            index_col = 1
            cols_to_drop = ['Cycle', 'Section', 'Type', 'Trigger Status']
        else:
            skipfooter = 0
            header = 0
            drop = 0
            index_col = 0

        if filetype == 'xlsx':
            imported = pd.ExcelFile(filename)
            data = imported.parse(
                0, index_col=index_col, skipfooter=skipfooter, header=header)
            data = data.drop(data.index[:drop], axis=0)

        elif filetype == 'csv':
            data = pd.read_csv(filename, sep=',', index_col=index_col, skipfooter=skipfooter,
                               header=header, engine='python')

        elif filetype == 'asc':
        	try:
            	data = pd.read_csv(filename, sep='\t', index_col=index_col, skipfooter=skipfooter, 
            	                   header=header, engine='python')
            except:
            	data = pd.read_csv(filename, sep=',', index_col=index_col, skipfooter=skipfooter,
            	                   header=header, engine='python')
                               
            data = data.drop(data.index[:drop], axis=0)
            data.dropna(axis=1, how='all', inplace=True)
            data = data.apply(pd.to_numeric, errors='coerce')

        else:
            warnings.warn('File type not supported.')

        if cols_to_drop:
            data = data.drop(cols_to_drop, axis=1)
            

        return data


class Background():
    """
    Class for calculating and visualising background.

    Parameters
    ----------
    isotope : MSData.Isotope
        Isotope from which background is calculated
    laser_on : list
        List of sets, where each set has 2 values (start of the laser and
        end of the laser.) indicating when the laser was fireing.
    laser_off : list
        List of sets, where each set has 2 values (end of the laser and
        start of the laser.) indicating when the laser wasn't fireing.    
    offset : float (Optional)
        Parameter for outlier removal. 
        1-offset = upper treshold for percentile filtering. 
        Accepts values between 0 and 1. Default is 0.15.
    width : float (Optional)
        Parameter for outlier removal.  Width od returned values, 
        where 1-offset-width = lower treshold for percentile filtering.
        Accepts values between 0 and 1. Default is 0.8.
    """

    def __init__(self, isotope, laser_on, laser_off, width=0.8, offset=0.15):
        assert isotope.data is not None

        # background by means
        self.bcg_all = []
        for i in range(len(laser_off)):
            data_points = isotope.data[laser_off[i][0]:laser_off[i][1]]
            if len(data_points) == 0:
                data_points = [0]
            if len(data_points) > 0:
                data_points = remove_outliers(data_points, offset, width)
            self.bcg_all.append(data_points)

        self.bcg_means = [x.mean() for x in self.bcg_all]

        # background by interpolation
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
        """
        Read Iolite file.

        Parameters
        ----------
        path : str
            Path to Iolite .csv file.
        """
        iolite = pd.read_csv(path, sep=",", engine='python')
        return iolite

    def names_from_iolite(self):
        """
        Return peak names from iolite. 
        """

        names = list(self.data[' Comment'].dropna())
        return names

    def on_and_off_times(self):
        times = self.data['Timestamp'][self.data.ne(self.data.shift()).apply(
            lambda x: x.index[x].tolist())[' Laser State']]
        return times


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
        """
        Read sheet with peak names from Param file.

        Returns
        ----------
        names : list
            List of peak names.
        """

        xl = pd.ExcelFile(self.path)
        if 'names' in xl.sheet_names:
            names = list(xl.parse('names', header=None)[0])
            return names
        else:
            if self.logger is not None:
                self.logger.error(
                    'Can not read names of peak from Param file.')

    def read_is_coef(self):
        """
        Read sheet with internal standard coefficients from Param file.

        Returns
        ----------
        internal_std : DataFrame
            DataFrame with internal standard coefficients.
        """

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
        """
        Read sheet with total sum correction coefficients from Param file.

        Returns
        ----------
        internal_std : dict
            DataFrame with total sum correction coefficients .
        """

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
    """
    Class for peak and background identification.

    Parameters
    ----------
    ms_data : MSData
        Mass spectrometry data class, for which the selector will be used.
    s : float (Optional)
        Start of first peak in seconds for synchronisation with Iolite data if using 
        method Iolite or for calculation  of backround standard deviation if using 
        method treshold. Default is 60s.
    sdmul: int (Optional)
        Multiplier for backround standard deviation in method treshold. Default is 10. 
    iolite: MSEval.Iolite (Optional)
        Iolite class, only necessary if method is Iolite. Default is None.
    logger: logger class (optional)
        If logger is pssed all methods of Selector will log in the activity.
    """

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
        """
        Select starts and ends of ablation using iolite file.
        """

        if self.logger is not None:
            self.logger.info('Selecting peak bounds by iolite.')

        times = self.iolite.on_and_off_times()

        starts = []
        ends = []
        for i, t in enumerate(times[1:]):
            timeindex = get_difference(get_timestamp(
                times[0]), get_timestamp(t)) + self.start
            if i % 2 == 0:
                starts.append(get_index(self.ms_data.data, timeindex))
            if i % 2 != 0:
                ends.append(get_index(self.ms_data.data, timeindex)+1)

        return starts, ends

    def create_selector_treshold(self):
        """
        Select starts and ends of ablation based on selected element or sum of all using treshold
        calculated from background.
        """

        if self.logger is not None:
            self.logger.info('Selecting peak bounds by setting treshold.')
        bcg_nr = self.ms_data.time_to_number(self.start)
        bcg = self.filter_line.iloc[0:bcg_nr].mean()
        std = self.filter_line.iloc[0:bcg_nr].std()
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
        Set time skipped on start and end of background and ablation in seconds.
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
        From starts and ends of ablation create laser_on and laser_off with skipped values.
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
