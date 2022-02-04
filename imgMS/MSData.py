import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
from scipy import interpolate
import warnings
from decimal import Decimal
import datetime
import logging
import itertools

from imgMS.side_functions import *
from imgMS.MSEval import *


class MSData():
    """
    LA-ICP-MS data structure. Primary object for LA-ICP-MS data reduction with imgMS.

    Parameters
    ----------
    datareader: MSEval.DataReader
        Object for reading and import of LA-ICP-MS data.
    logger: logger class (optional)
        If logger is pssed all methods of MSData will log in the activity.
    """

    def __init__(self, datareader=None, logger=None):
        self.logger = logger
        if datareader is not None:
            self.datareader = datareader
            self.data = self.datareader()
            self.data.index = list(map(float, self.data.index))
            self.time = np.array(self.data.index)
            self.matrix = self.data.values
            self.isotope_names = np.array(self.data.columns)
            if self.logger is not None:
                self.logger.info(f'Reading data {self.datareader.filename}.')
        else:
            self.datareader = None
            self.data = None
            self.time = None
            self.matrix = None
            self.isotope_names = []
            if self.logger is not None:
                self.logger.info(f'Creating empty MSData.')

        self.isotopes = {}
        self.param = None
        self.selector = None
        self.starts, self.ends = None, None
        self.laser_off, self.laser_on = None, None
        self.names = None
        self.srms = None
        self.means = None
        self.quantified = None
        self.lod = None
        self.corrected_IS = None
        self.corrected_TS = None

        self.create_isotopes()

    def __call__(self, isotopes=None, *args, **kwargs):
        """
        Plots data in time dependant graph.

        Parameters
        ----------
        isotopes: list (optional)
            List of isotopes to be shown in the plot. |If not specified plots all isotopes in measurement.
        *args, **kwargs:
            All other plotting arguments to be passed to matplotlib.pyplot.plot.
        """

        plot_data(self.data, isotopes=isotopes, *args, **kwargs)

    def __repr__(self):
        res = f'{self.__class__.__name__}\n Analysis time: {self.time.max()}s\n File path: {self.datareader.filename}\n Measured isotopes: {self.isotope_names}'
        return res

    def time_to_number(self, time, integration_time=None):
        """
        Takes time in seconds returns number of measured values.
        The result depends on integration time of MS method.

        Parameters
        ----------
        time: float
            Time in seconds to be converted into number of values in data.
        integration_time: float (optional)
            Integration time of LA-ICP-MS measurement. If not specified will be calculated from data.
        """

        if not integration_time:
            integration_time = (self.time[2]-self.time[1])
        val = int(abs(time)//integration_time)
        if time < 0:
            val = -val
        return val

    def read_param(self, path):
        """
        Import excel file with additional parameters (names of peaks, internal standard values,
        total sum correction coefficients) for data reduction. A sample PARAM file can be found in data folder.

        Parameters
        ----------
        path: str
            Path to excel param file.
        """

        self.param = Param(path, self.logger)
        if self.param.peak_names:
            self.names = self.param.peak_names

    def read_srms(self, path='./SRM.xlsx'):
        """
        Import excel file with standard reference values. Default file is part of the imgMS package
        and contains values for NIST610, NIST612 and NIST 614.

        Parameters
        ----------
        path: str
            Path to excel SRM file.
        """

        self.srms = pd.ExcelFile(path).parse(index_col=0)

    def set_names(self, names):
        """
        Sets list of names of peaks. The name of SRM must be equal to the name in SRM file.

        Parameters
        ----------
        names: list
            List of names.
        """

        self.names = names

    def create_isotopes(self):
        for key in self.isotope_names:
            self.isotopes[key] = Isotope(key, self, logger=self.logger)

    def select(self, method='treshold', selector=None, s=60, sdmul=10, iolite=None):
        """
        Selects starts and ends of peaks using imgMS.Selector.

        Parameters
        ----------
        method: str (optional)
            Name of the method to be used for identifying peaks. Possible options are 'treshold' and 'iolite'.
            Default is 'treshold'.
        selector: MSEval.Selector (optional)
            Class for identifying peaks. If Selector is passed none of the other parameters are necessary. If not,
            Selector is created by selected settings.
        s: float (optional)
            Start of the first peak in seconds from the start of analysis. Default is 60. Necessary if Selector is not
            passed and used for synchronisation of data with iolite.
        sdmul: float (optional)
            Coeficient by which a standard deviation of background is multiplied to calculate treshold. Only used if
            method = treshold.
        iolite: MSEval.Iolite (optional)
            Iolite class holding data from .Iolite file. Necessary if method = Iolite and Selector not passed.
        """

        if selector is None:
            self.selector = Selector(
                self, s=s, sdmul=sdmul, iolite=iolite, logger=self.logger)
        else:
            self.selector = selector

        self.selector.method = method
        self.starts, self.ends = self.selector()
        self.laser_on, self.laser_off = self.selector.create_on_off(
            self.starts, self.ends)

    def graph(self, ax=None, logax=False, el=None, *args, **kwargs):
        """
        Create matplotlib graph of intensity in time for ablation and highlights peaks and background signal
        if the peaks are already identifyied.

        Parameters
        ----------
        ax: matplotlib axes (optional)
            Axes to plot in, if not specified, create new ax.
        logax: bool (optional)
            If True use logarythmic x axes. Default False.
        el: str (optional)
            Element to plot. If not specified plot all measured elements.
        *args, **kwargs:
            All other plotting arguments to be passed to matplotlib.pyplot.plot.
        """

        if ax == None:
            fig, ax = plt.subplots()

        ax.cla()
        ax.clear()
        # if element is defined, plot only one element, otherwise all
        if el:
            self.data.plot(ax=ax, y=el, kind='line', legend=False)
        else:
            self.data.plot(ax=ax, kind='line', legend=True)
            ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.075),
                      fancybox=True, shadow=True, ncol=10)

        if logax:
            ax.set_yscale('log')

        if self.starts and self.ends:
            # create lines for start and end of each ablation
            for i in range(0, len(self.starts)):
                ax.axvline(x=self.time[self.starts[i]],
                           color='blue', linewidth=2)
            for i in range(0, len(self.ends)):
                ax.axvline(x=self.time[self.ends[i]],
                           color='blue', linewidth=2)

        if self.laser_off:
            # higlights bacground
            for off in self.laser_off:
                try:
                    ax.axvspan(
                        self.time[off[0]], self.time[off[1]], alpha=0.2, color='red')
                except:
                    warnings.warn('something is wrong')

        if self.laser_on:
            # higlihts ablation
            for on in self.laser_on:
                ax.axvspan(self.time[on[0]], self.time[on[1]],
                           alpha=0.2, color='green')

        plt.show()

    def average_isotopes(self, despiked=True, bcgcor_method='all', method='intensity'):
        """
        Calculate average value for all peaks for all isotopes.

        Parameters
        ----------
        bcgcor_method : str
            Method of background calculation. Possible options are [None, 'all', 'beginning', 'end']. Default is 'all'.
        despiked : bool
            If True use despiked data, else use original. Default is True.
        method : str
            Method to get average peaks. Possible options are ['intensity', 'integral']. Default is 'intensity'.

        Returns
        -------
        data : DataFrame
            DF where index are peak names and columns are isotopes
        """

        if self.logger is not None:
            self.logger.info(
                f'Averaging of signal using despiking: {despiked}, bcgcor: {bcgcor_method}, method: {method}.')

        self.means = pd.DataFrame()
        for el, isotope in self.isotopes.items():
            isotope.average_intensity(
                despiked=despiked, bcgcor_method=bcgcor_method, method=method)
            self.means[el] = isotope.means

        if self.names:
            self.means.index = self.names
        else:
            if self.logger is not None:
                self.logger.error('Missing peak names.')
        return self.means

    def quantify_isotopes(self, srm_name='NIST610', drift_correction=True):
        """
        Calculate quantified value for all peaks for all isotopes.

        Parameters
        ----------
        srm_name : str
            Standard reference material used for quantification. The name must be at least one of the
            peaks and listed in SRM file.
        Returns
        -------
        data : DataFrame
            Quantified data in DF where index are peak names and columns are isotopes
        """

        self.quantified = pd.DataFrame()
        for el, isotope in self.isotopes.items():

            if element_strip(el) not in self.srms.columns:
                if self.logger is not None:
                    self.logger.error(f'Missing srm {el}.')

            isotope.quantify(srm_name=srm_name,
                             drift_correction=drift_correction)
            self.quantified[el] = isotope.quantified

        if self.logger is not None:
            self.logger.info(
                f'Quantification of signal using SRM: {srm_name}.')

        self.quantified.index = [
            name for name in self.names if name != srm_name]
        return self.quantified

    def IS_correction(self):
        """
        Calculates correction for each element given in internal standard correction
        from PARAM file.

        Returns
        -------
        corrected data : dict
            dict of internal standards used for correction as keys and DataFrames
            where index are peak names and columns are isotopes with values in ppm.
        """
        self.corrected_IS = {}

        if self.param is not None:
            if self.param.is_coef.empty:
                if self.logger is not None:
                    self.logger.error(
                        'Missing coefficients of internal standards.')
        else:
            if self.logger is not None:
                self.logger.error(
                    'Missing param file for IS correction.')

        self.correction_elements = list(self.param.is_coef.columns)

        for el in self.correction_elements:
            if el in list(map(elem_resolution, self.isotope_names)):
                self.corrected_IS[el] = (correction(
                    self.quantified, el, self.param.is_coef))
        return self.corrected_IS

    def TS_correction(self, suma=1000000, skip_isotopes=[], return_oxides=False):
        """
        Calculates total sum correction [1] using coefficients given in PARAM file.
        If coefficients in PARAM file are not given, uses default values. Default values
        assume all elements are in most common oxide form.

        [1] Liu, Y., Hu, Z., Gao, S., Günther, D., Xu, J., Gao, C. and Chen, H., 2008.
        In situ analysis of major and trace elements of anhydrous minerals by LA-ICP-MS
        without applying an internal standard. Chemical Geology, 257(1-2), pp.34-43.

        Parameters
        ----------
        suma : float
            Total sum of measured elements in ppm used for correction.
            Normally is equal to 100% (default).
        skip_isotopes : list
            List of isotopes to be skipped, in total sum correction one element
            can't be measured on multiple isotopes.
        return_oxides : bool
            If True return data in oxide form as was used for total sum correction.
            Default is False.

        Returns
        -------
        corrected data : DataFrame
            DataFrame where index are peak names and columns are isotopes with values in ppm.
        """
        if self.param is not None:
            if not self.param.ts_coef:
                if self.logger is not None:
                    self.logger.error(
                        'Missing coeficients for total sum correction.')
            else:
                ts_coef = self.param.ts_coef
        else:
            if self.logger is not None:
                self.logger.error(
                    'Using default TS coef.')
            ts_coef = pd.ExcelFile(
                './default_sum_koef.xlsx').parse(0, index_col=0, header=None).to_dict()[1]

        if self.logger is not None:
            self.logger.info(
                f'Using these coeficients for total sum: {ts_coef}.')

        self.corrected_TS = self.quantified.copy()
        self.corrected_TS.drop(skip_isotopes, axis=1, inplace=True)

        elem_lst = [element_strip(el) for el in self.corrected_TS.columns]
        if len(elem_lst) != len(set(elem_lst)):
            if self.logger is not None:
                self.logger.warning(
                    'Multiple isotopes of one element measured. Recomended to use only one isotope for total sum correction.')

        # calculate oxide form of specified elements
        for key in ts_coef:
            elem = element_formater(key, self.corrected_TS.columns)
            if not elem:
                continue
            self.corrected_TS[elem] = self.corrected_TS[elem] / \
                ts_coef[key] * 100

        # calculate coefficient for correction and correct data
        koef = suma/self.corrected_TS.sum(1)
        self.corrected_TS = self.corrected_TS.multiply(np.array(koef), axis=0)

        # Return oxide form back to isotopes if not specified otherwise
        if not return_oxides:
            for key in ts_coef:
                elem = element_formater(key, self.corrected_TS.columns)
                if not elem:
                    continue
                self.corrected_TS[elem] = self.corrected_TS[elem] * \
                    ts_coef[key] / 100
        else:
            ox_names = {el: el + ' oxide' for el in ts_coef.keys()}
            self.corrected_TS.rename(columns=ox_names, inplace=True)

        return self.corrected_TS

    def detection_limit(self, method='intensity', scale='all', ablation_time=60):
        """
        Calculate detection limit for all isotopes.

        Parameters
        ----------
        method : str
            Possible methods ['integral','intensity']. LoD must be calculates same way as average of isotopes.
        scale : str
            Possible methods ['beginning', 'all']. LoD must be calculates same way as average of isotopes.
        ablation_time : int
            Time in seconds of one ablation. Necessary only for method 'integral'.
        """
        if self.logger is not None:
            self.logger.info(
                f'Calculating detection limit by method: {method}, scale: {scale}')
        self.lod = pd.Series()
        for el, isotope in self.isotopes.items():
            isotope.detection_limit(method, scale, ablation_time)
            self.lod[el] = isotope.lod
        self.lod.name = 'lod'
        return self.lod

    def report(self):
        """
        Clean all data in MSData. Replace values lower than limit of detection.
        If the value is above LoD, round to specific decimal place.
        """
        if self.corrected_TS is not None:
            ts_lod = pd.Series({el: self.lod[el] if el in self.lod.index else self.lod[el.replace(
                ' oxide', '')] for el in self.corrected_TS.columns}, name='LoD')
            self.corrected_TS = self.corrected_TS.append(ts_lod)
            for column in self.corrected_TS:
                elem = element_formater(column, list(self.lod.index))
                self.corrected_TS[column] = [
                    report(value, self.lod, elem) for value in self.corrected_TS[column]]

        if self.corrected_IS is not None:
            self.corrected_IS = {key: df.append(self.lod)
                                 for key, df in self.corrected_IS.items()}
            for el, df in self.corrected_IS.items():
                for column in df:
                    df[column] = [report(value, self.lod, column)
                                  for value in df[column]]

        if self.quantified is not None:
            self.quantified = self.quantified.append(self.lod)
            for column in self.quantified:
                self.quantified[column] = [
                    report(value, self.lod, column) for value in self.quantified[column]]

    def export(self, path):
        """
        Export all data to excel file.

        Parameters
        ----------
        path : str
            Path to excel file where to save data.
        """

        if self.logger is not None:
            self.logger.info(
                f'Saving data to: {path}.')

        writer = pd.ExcelWriter(path, engine='xlsxwriter')
        data = {'raw': self.data,
                'average': self.means,
                'quantified': self.quantified}

        if self.corrected_IS is not None:
            data.update(self.corrected_IS)

        if self.corrected_TS is not None:
            data['corrected_TS'] = self.corrected_TS

        for key, df in data.items():
            df.to_excel(writer, sheet_name=key)

        writer.save()

    def formated_export(self, path, table='quantified'):
        """
        Export all data to excel file in a formated style, containing means, sd, max and min for repeated samples.

        Parameters
        ----------
        path : str
            Path to excel file where to save data.
        """

        if self.logger is not None:
            self.logger.info(
                f'Saving data to: {path}.')

        if table == 'quantified':
            frame = self.quantified
            formatted_export(frame, path)
        elif table == 'TScorrected':
            frame = self.corrected_TS
            formatted_export(frame, path)
        elif 'IScorrected' in table:
            for elem in self.corrected_IS.keys():
                formatted_export(self.corrected_IS[elem], "{0}_{2}.{1}".format(
                    *path.rsplit('.', 1) + [elem]))
        else:
            raise ValueError(
                f'{table} not in possible otions for table. Use quantified, TScorrected or IScorrected')

    def create_maps(self, despiked=False, bcgcor_method='all', dx=1, dy=1):
        """
        Create elemental distribution matrix for all isotopes.

        Parameters
        ----------
        despiked : bool
            If True use despiked data, else use original. Default is False.
        bcgcor_method : str
            Method of background calculation. Possible options are [None,
            'all', 'beginning', 'end']. Default is 'all'.
        dx : float
            X-axis distance between two values, usually in μm. Can be
            calculated by (scan speed [μm/s]/ integration time [s])
        dy : float
            Y-axis distance between two values, usually in μm. Distance
            between two lines, usually equal to ablation spot size.)
        """

        for el, isotope in self.isotopes.items():
            isotope.elemental_distribution(despiked, bcgcor_method, dx, dy)

    def quantify_maps(self, slopes, intercepts):
        """
        Calculate quantification for elemental distribution matrix for all isotopes.

        Parameters
        ----------
        intercepts : dict
            Dict of intercepts for each isotope.
        slopes : dict
            Dict of slopes for each isotope.
        """

        for el, isotope in self.isotopes.items():
            isotope.elmap.quantify_map(slopes[el], intercepts[el])

    def export_matrices(self, path, quantified=False):
        """
        Export all elemental distribution data to excel file. Each isotope
        will be a matrix on one sheet.

        Parameters
        ----------
        path : str
            Path to excel file where to save data.
        """

        writer = pd.ExcelWriter(path, engine='xlsxwriter')
        for el, isotope in self.isotopes.items():
            if quantified is False:
                df = pd.DataFrame(isotope.elmap.matrix,
                                  index=isotope.elmap.y, columns=isotope.elmap.x)
            else:
                df = pd.DataFrame(isotope.elmap.qmap,
                                  index=isotope.elmap.y, columns=isotope.elmap.x)
            df.to_excel(writer, sheet_name=el)
        writer.save()

    def import_matrices(self, path):
        """
        Import all elemental distribution data from excel file. Each isotope
        will be a matrix on one sheet.

        Parameters
        ----------
        path : str
            Path to the excel file.
        """
        for el, isotope in self.isotopes.items():
            isotope.elmap.read_matrix(path, el)


class Isotope():
    """
    Extract specific Isotope data from MSData object.

    Parameters
    ----------
    isotope_name : str
        Isotope to extract from MSData (i.e. Li7)
    ms_data : MSData
        MS data to extract from.

    Returns
    -------
    data : array
        `data` interpreted as an array.
    """

    def __init__(self, isotope_name, ms_data=None, logger=None):
        self.isotope_name = elem_resolution(isotope_name)
        if ms_data is not None:
            self.ms_data = ms_data
            self.data = ms_data.data[isotope_name].values
            self.time = ms_data.data.index
        self.isotope_number = int(
            ''.join([c for c in self.isotope_name if c.isnumeric()]))
        self.despiked = None
        self.bcg = None
        self.bcg_corrected = None
        self.bcgcor_method = None
        self.means = None
        self.ratio = None
        self.quantified = None
        self.elmap = None
        self.logger = logger

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        if isinstance(idx, (int, slice)):
            return self.data[idx]
        return [self.data[i] for i in idx]

    def __call__(self):
        plt.plot(self.data)

    def __repr__(self):
        res = f'{self.__class__.__name__} {self.isotope_name} ({len(self)} values)\n{self.data[:10]}'
        if len(self) > 10:
            res = res[:-1] + '...]'
        return res

    def __iter__(self): return iter(self.data)
    def __setitem__(self, i, o): self.data[i] = o
    def __delitem__(self, i): del(self.data[i])

    def element(self):
        """
        Returns name of the element.
        """

        return ''.join([c for c in self.isotope_name if c.isalpha()])

    def despike(self, win=3, treshold=10):
        """
        Removes high outliers of the data using rolling mean and standard deviation.
        Replace values that are greater than n standard deviations above the mean
        with the mean adjacent values.

        Parameters
        ----------
        win: int
            Window for rolling mean.
        tresholds: int
            Multiplier of std above which values will be replaced.
        """

        sig = pd.Series(self.data)
        skip = int((win-1)/2)
        outliers = np.zeros(len(sig), dtype=bool)
        outliers[skip:-skip] = True
        rm = np.array(sig.rolling(win).mean().dropna())
        rs = np.array(sig.rolling(win).std().dropna())**0.5
        outliers[skip:-skip] = sig[skip:-skip] > rm+rs*treshold
        sig[outliers] = rm[outliers[skip:-skip]]
        self.despiked = sig

    def bcg_correction(self, method='all', despiked=True):
        """
        Calculate bacground correction for selected Isotope.

        Parameters
        ----------
        method : str
            Method of background calculation. Possible options are ['all', 'beginning', 'end']. Default is 'all'.
        despiked : bool
            State of data to be corrected. If True use despiked data, else use original. Default is True.
        """

        if not self.ms_data.laser_on:
            if self.logger is not None:
                self.logger.error(
                    'Bacground correction requires selection of peaks first.')
            return

        self.bcg = Background(
            self, self.ms_data.laser_on, self.ms_data.laser_off)

        if not despiked:
            data = np.copy(self.data)
        elif despiked:
            if self.despiked is None:
                self.despike()
            data = np.copy(self.despiked)

        if method == 'beginning':
            data_out = data-self.bcg.bcg_means[0]
        elif method == 'end':
            data_out = data-self.bcg.bcg_means[-1]
        elif method == 'all':
            data_out = data-self.bcg.bcg_interp
        else:
            if self.logger is not None:
                self.logger.warning(
                    'Select correct method  for background correction. Possible options: [\'beginning\', \'end\', \'all\']')
            return

        data_out[data_out < 0] = 0
        self.bcg_corrected = data_out
        self.bcgcor_method = method

    def average_intensity(self, despiked=True, bcgcor_method='all', method='intensity'):
        """
        Calculate average intensity for peaks in selected Isotope.

        Parameters
        ----------
        bcgcor_method : str
            Method of background calculation. Possible options are [None, 'all', 'beginning', 'end']. Default is 'all'.
        despiked : bool
            State of data to be corrected. If True use despiked data, else use original. Default is True.
        method : str
            Method to get average peaks. Possible options are ['intensity', 'integral']. Default is 'intensity'.
        """
        laser_on = self.ms_data.laser_on

        if bcgcor_method is not None:
            if self.bcg_corrected is None or self.bcgcor_method != bcgcor_method:
                self.bcg_correction(method=bcgcor_method, despiked=despiked)
                data = np.copy(self.bcg_corrected)
            else:
                data = np.copy(self.bcg_corrected)
        else:
            if not despiked:
                data = np.copy(self.data)
            elif despiked:
                if self.despiked is None:
                    self.despike()
                data = np.copy(self.despiked)

        means = []
        for on in laser_on:
            if method == 'intensity':
                means.append(np.mean(data[on[0]:on[1]]))
            elif method == 'integral':
                sample_y = data[on[0]:on[1]]
                sample_x = self.time[on[0]:on[1]]
                means.append(np.trapz(sample_y, sample_x))
            else:
                if self.logger is not None:
                    self.logger.error(
                        'Select correct method for peak averaging. Possible options: [\'intensity\', \'integral\']')
                return
        self.means = np.array(means)

    def quantify(self, srm_name='NIST610', drift_correction=True):
        """
        Calculate quantification of average peaks in Isotope.

        Parameters
        ----------
        srm_name : str
            Name of the standard reference material to be used for quantification. Must be in SRM file.
        """
        if self.ms_data.srms is not None:
            if srm_name in self.ms_data.srms.index:
                srm = self.ms_data.srms.loc[srm_name, :]
            else:
                if self.logger is not None:
                    self.logger.error(
                        f'Standard reference material {srm_name} not in file.')
                return
        else:
            if self.logger is not None:
                self.logger.error('Missing standard reference materials file.')
            return

        if self.means is None:
            if self.logger is not None:
                self.logger.error('Missing average values for peaks.')
            return

        if self.ms_data.names is None:
            if self.logger is not None:
                self.logger.error('Missing names of peaks.')
            return

        if element_strip(self.isotope_name) not in self.ms_data.srms.columns:
            if self.logger is not None:
                self.logger.error(f'Unknown isotope: {self.isotope_name}.')
            return

        if drift_correction:
            if self.ms_data.names.count(srm_name) < 2:
                stdsig = np.mean(
                    [[m for m, s in zip(self.means, self.ms_data.names) if s == srm_name]])
            else:
                # stdsig_mskd = np.array(
                #     [m if s == srm_name else np.nan for m, s in zip(self.means, self.ms_data.names)])
                # not_nan = np.logical_not(np.isnan(stdsig_mskd))
                # indices = np.arange(len(stdsig_mskd))
                # stdsig = np.interp(indices, indices[not_nan], stdsig_mskd[not_nan])

                x = np.array([i for i, s in enumerate(
                    self.ms_data.names) if s == srm_name])
                y = np.array(
                    [m for m, s in zip(self.means, self.ms_data.names) if s == srm_name])
                f = interpolate.interp1d(
                    x, y, kind='slinear', fill_value="extrapolate")

                xnew = np.arange(len(self.means))
                stdsig = f(xnew)

                xnew = np.arange(len(self.means))
                stdsig = f(xnew)
        else:
            stdsig = np.mean(
                [[m for m, s in zip(self.means, self.ms_data.names) if s == srm_name]])

        # TODO: interpolate standard signal by time, not every spot has to be the same length?

        self.ratio = srm[element_strip(self.isotope_name)]/stdsig
        quantified = self.means * self.ratio

        self.quantified = np.array(
            [q for q, s in zip(quantified, self.ms_data.names) if s != srm_name])

    def detection_limit(self, method='intensity', scale='all', ablation_time=60):
        """
        Calculate limit of detection for isotope.

        Parameters
        ----------
        method : str
            Possible methods ['integral','intensity']. LoD must be calculates same way as average of isotopes.
        scale : str
            Possible methods ['beginning', 'all']. LoD must be calculates same way as average of isotopes.
        ablation_time: int
            Time in seconds of one ablation. Necessary only for method 'integral'.
        """
        if scale == 'all':
            bcg = np.array(
                list(itertools.chain.from_iterable(self.bcg.bcg_all)))
        elif scale == 'beginning':
            bcg = np.array(self.bcg.bcg_all[0])
        elif scale == 'end':
            bcg = np.array(self.bcg.bcg_all[-1])
        else:
            if self.logger is not None:
                self.logger.error(f'Unknown scale: {scale}.')

        if self.ratio is not None:
            if method == 'integral':
                self.lod = (bcg.std()*ablation_time*self.ratio)
            elif method == 'intensity':
                self.lod = (bcg.std()*3*self.ratio)
        else:
            self.lod = 0

    def elemental_distribution(self, despiked=False, bcgcor_method='all', dx=1, dy=1):
        """
        Creates ElementalMap from Isotope data.

        Parameters
        ----------
        despiked : bool (Optional)
            State of data to be corrected. If True use despiked data, false use original. Default is False.
        bcgcor_method : str (Optional)
            Method of background calculation. Possible options are [None, 'all', 'beginning', 'end']. Default is 'all'.
        dx : float
            X-axis distance between two values, usually in μm. Can be
            calculated by (scan speed [μm/s]/ integration time [s])
        dy : float
            Y-axis distance between two values, usually in μm. Distance
            between two lines, usually equal to ablation spot size.)
        """

        laser_on = self.ms_data.laser_on

        if bcgcor_method is not None:
            if self.bcgcor_method != bcgcor_method:
                self.bcg_correction(method=bcgcor_method, despiked=despiked)
            data = np.copy(self.bcg_corrected)
        else:
            if not despiked:
                data = np.copy(self.data)
            elif despiked:
                if self.despiked is None:
                    self.despike()
                data = np.copy(self.despiked)

        l = []
        for on in laser_on:
            arr = data[on[0]:on[1]]
            arr[arr < 0] = 0
            l.append(arr)

        lens = [len(r) for r in l]
        maxlen = max(lens)
        mat = np.empty([len(l), maxlen])
        mat[:] = np.nan
        mask = np.arange(maxlen) < np.array(lens)[:, None]
        mat[mask] = np.concatenate(l)

        self.elmap = ElementalMap(mat, dx, dy)


class Peak():
    """
    Separate one peak of isotope data.

    Parameters
    ----------
    isotope: MSData.Isotope
        Isotope to use for peak extraction.
    bounds: set
        Set of 2 values, start and end of peak.
    """

    def __init__(self, isotope, bounds):
        if isotope.data is None:
            raise ValueError(
                'Missing data for isotope: {}'.format(isotope.isotope_name))
        self.data = isotope.data[bounds[0]:bounds[1]]

    def __len__(self):
        return len(self.data)

    def __call__(self, *args, **kwargs):
        plt.plot(self.data, *args, **kwargs)

    def __repr__(self):
        pass


class ElementalMap():
    def __init__(self, matrix=None, dx=None, dy=None):
        """
        Class of Elemental map holding data from MSData representing elemental distribution of one isotope.

        Parameters
        ----------
        matrix : numpy matrix
            Matrix of z values of image.
        dx: float
            Step between z values in x dirextion.
        dy: float
            Step between z values in y dirextion.
        """
        self.matrix = matrix
        self.dx = dx
        self.dy = dy
        self.qmap = None
        if self.matrix is not None and self.dx is not None and self.dy is not None:
            self.create_xy()

    def __call__(self, fig=None, ax=None, units='', axis=True, clb=True, quantified=False, title='', *args, **kwargs):
        """
        Show image of elemental map.
        Parameters
        ----------
        ax : matplotlib axes
            Axes where plot will be drawn, if None new axes will be created.
        units : str
            Label for colorbar
        axis : boolean (Optional)
            If true then show axis with labels in the image. Default is True.
        clb : boolean (Optional)
            Wheather to show colorbal in image. Default is True.
        quantified : bool (optional)
            Whether to plot intensities or quantified map. Default is False.
        title: str (Optional)
            Figure title to be used. Default is no title.
        """
        if fig is None or ax is None:
            fig, ax = plt.subplots()
        ax.cla()
        if quantified is False:
            im = ax.imshow(self.matrix, extent=[
                0, self.x[-1], self.y[-1], 0], *args, **kwargs)
        else:
            im = ax.imshow(self.qmap, extent=[
                0, self.x[-1], self.y[-1], 0], *args, **kwargs)

        if not axis:
            ax.axis('off')

        if clb:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            clb = fig.colorbar(im, cax=cax)
            clb.ax.set_title(units)

        fig.suptitle(title)

    def __repr__(self):
        return f'{self.__class__.__name__}: {self.matrix.shape}'

    def create_xy(self):
        """
        Create values of x and y axis from dx and dy.
        """

        self.x = [i*self.dx for i in range(self.matrix.shape[1])]
        self.y = [i*self.dy for i in range(self.matrix.shape[0])]

    def read_matrix(self, path, el):
        """
        Import created map from excel, where each isotope is on one list.

        Parameters
        ----------
        path : str
            Path to excel file.
        el: str
            Name of isotope to import.
        """

        file = pd.ExcelFile(path)
        df = file.parse(el, index_col=0)
        self.matrix = np.matrix(df)
        self.dx = df.columns[1]-df.columns[0]
        self.dy = df.index[1]-df.index[0]
        self.create_xy()

    def write_matrix(self, path, el):
        """
        Save created map to excel, where each isotope is on one list.

        Parameters
        ----------
        path : str
            Path to excel file.
        el: str
            Name of isotope to export. Will be used as sheet name.
        """

        writer = pd.ExcelWriter(path, engine='xlsxwriter')
        df = pd.DataFrame(self.matrix, index=self.y, columns=self.x)
        df.to_excel(writer, sheet_name=el)
        writer.save()

    def rotate(self):
        """
        Rotate map by 90 degrees.
        """

        self.matrix = np.rot90(self.matrix)
        self.x, self.y = self.y, self.x

    def flip(self, axis=None):
        """
        Flip matrix of image.

        Parameters
        ----------
        axis : str
            Axis along which to flip over. The default, axis=None,
        will flip over all of the axes, axis=0 flip matrix vertically,
        axis=1 flip matrixhorizontally.
        """

        self.matrix = np.flip(self.matrix, axis)

    def quantify_map(self, slope, intercept=0):
        """
        Quantify elemental map.

        Parameters
        ----------
        slope: float
            Slope of linear regression for quantification.
        intercept: float
            Intercept of linear regression for quantification.
        """

        self.qmap = (self.matrix-intercept)/slope
