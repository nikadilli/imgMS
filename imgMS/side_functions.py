import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from scipy import stats
import warnings
from decimal import Decimal
import datetime
import logging


def get_logger(pathname):
    """
    Create log file similar to Ilaps-GUI log file

    Parameters
    ----------
    pathname : str
    Path to a log file(.txt). If file doesn't exists it is created.
    """
    logging.basicConfig(level=logging.INFO, format='%(levelname)s (%(name)s): %(asctime)s >>> %(message)s',
                        datefmt='%d-%b-%y %H:%M:%S', filename=pathname)
    logger = logging.getLogger("main")
    return logger


def get_timestamp(strTime):
    """format string time from iolite to timestamp"""
    return datetime.datetime.strptime(strTime, '%Y-%m-%d %H:%M:%S.%f')


def get_difference(start, now):
    """return time in seconds between 2 timestamps"""
    diff = now - start
    return diff.total_seconds()


def get_index(data, time):
    """return closest index of MS time given time in seconds"""
    for i in range(len(data.index)-1):
        if (data.index[i] <= time) and (data.index[i+1] > time):
            return i


def get_diff_lst(iolite):
    """return list of times in seconds from start to every start and end of laser ablation for spots"""
    lst = []
    for i in range(1, len(iolite['Timestamp'])-1):
        if (i-4) % 5 == 0:
            lst.append(get_difference(get_timestamp(
                iolite.loc[i-2, 'Timestamp']), get_timestamp(iolite.loc[i, 'Timestamp'])))
            lst.append(get_difference(get_timestamp(
                iolite.loc[i, 'Timestamp']), get_timestamp(iolite.loc[i+1, 'Timestamp'])))
    lst.append(get_difference(get_timestamp(
        iolite.loc[i, 'Timestamp']), get_timestamp(iolite.loc[i+1, 'Timestamp'])))
    return lst


def get_diff_lst_line(iolite):
    """return list of times in seconds from start to every start and end of laser ablation for lines"""
    lst = []
    for i in range(1, len(iolite['Timestamp'])-1):
        if (i-6) % 7 == 0:
            lst.append(get_difference(get_timestamp(
                iolite.loc[i-2, 'Timestamp']), get_timestamp(iolite.loc[i, 'Timestamp'])))
            lst.append(get_difference(get_timestamp(
                iolite.loc[i, 'Timestamp']), get_timestamp(iolite.loc[i+1, 'Timestamp'])))
    lst.append(get_difference(get_timestamp(
        iolite.loc[i-2, 'Timestamp']), get_timestamp(iolite.loc[i, 'Timestamp'])))
    return lst


def element_strip(elem):
    """formatting of element name to match colnames of reference material"""
    elem = elem.replace('(LR)', '').replace('(MR)', '').replace('(HR)', '')
    elem = ''.join([c for c in elem if c.isalpha()])
    return elem


def elem_resolution(elem):
    """formatting of element name to remove resolution for ELEMENT2 instrument"""
    elem = elem.replace('(LR)', '').replace('(MR)', '').replace('(HR)', '')
    return elem


def remove_outliers(data, offset, width):
    """
    Function to filter data by percentile value.

    Parameters
    ----------
    data: nparray
        array of data
    offset: float
        1-offset = upper treshold for percentile filtering. 
        Accepts values between 0 and 1.
    width: float
        Width od returned values, where 1-offset-width = lower treshold for percentile filtering.
        Accepts values between 0 and 1.

    Returns
    -------
    data_out : nparray
        Filtered data.
    """

    if offset < 0 or 1 < offset:
        raise ValueError(f'offset={offset}: Offset must be between 0 and 1')

    if width < 0 or 1 < width:
        raise ValueError(f'Width={width}: Width must be between 0 and 1')

    upper, lower = np.percentile(data, [round(
        100*(1-offset), 4), round(100*(1-offset-width), 4)])  # round to clear float error

    data_out = np.array([x for x in data if x >= lower and x <= upper])

    return data_out


def plot_data(data, isotopes=None, ax=None, *args, **kwargs):
    """
    Create a plot of MS time dependant data.

    Parameters
    ----------
    data : DataFrame
        Dataframe of ms data where column names are names of measured isotopes and index is time.
    isotopes : list
        List of isotopes to plot. (Optional)
    ax : matplotlib axes
        Matplotlib axes to show the plot. (Optional)

    """
    if ax is None:
        fig, ax = plt.subplots()
    if isotopes is not None:
        ax.plot(data[isotopes], *args, **kwargs)
    else:
        ax.plot(data, *args, **kwargs)


def correction(data, elem, internal_std):
    """
    Calculates internal standard correction.

    Parameters
    ----------
    data: DataFrame
        DataFrame with quantified values where columns are measured isotopes
    el: str
        Element used as internal standard
    internal_std: DataFrame
        DataFrame of values of internal standard  where columns are elements for 
        correction and each row represents one measurement e.g. Spot.
    """
    ratio = data[element_formater(elem, data.columns)].div(
        list(internal_std[elem]))
    return data.apply(lambda x: x/ratio)


def element_formater(elem, lst_of_elems):
    """matches the given element format to the one used in list"""
    if elem in lst_of_elems:
        return elem
    elif elem not in lst_of_elems:
        elem = elem.replace(' oxide', '')
        if elem in lst_of_elems:
            return elem
        elif elem not in lst_of_elems:
            elem = elem.replace('(LR)', '').replace(
                '(MR)', '').replace('(HR)', '')
            if elem in lst_of_elems:
                return elem
            elif elem not in lst_of_elems:
                elem = elem + '(LR)'
                if elem in lst_of_elems:
                    return elem
                elif elem not in lst_of_elems:
                    elem = elem.replace('(LR)', '')
                    elem = elem + '(MR)'
                    if elem in lst_of_elems:
                        return elem
                    elif elem not in lst_of_elems:
                        elem = elem.replace('(MR)', '')
                        elem = elem + '(HR)'
                        if elem in lst_of_elems:
                            return elem
                        elif elem not in lst_of_elems:
                            elem = elem.replace('(HR)', '')
                            elem = ''.join([c for c in elem if c.isalpha()])
                            if elem in lst_of_elems:
                                return elem
                            else:
                                return


def report(x, LoD, elem):
    """
    Replace values lower than limit of detection.
    If the value is above LoD, round to specific decimal place.
    """
    if x < LoD[elem]:
        return '< LoD'
    else:
        if x >= 100:
            return float(Decimal(x).quantize(Decimal('1')))
        else:
            if x >= 10:
                return float(Decimal(x).quantize(Decimal('0.1')))
            else:
                if x >= 1:
                    return float(Decimal(x).quantize(Decimal('0.01')))
                else:
                    if x >= 0:
                        return float(Decimal(x).quantize(Decimal('0.001')))


ox_name = {'Na23': 'Na2O (%)',
           'Mg24': 'MgO (%)',
           'Al27': 'Al2O3 (%)',
           'Si28': 'SiO2 (%)',
           'P31': 'P2O5 (%)',
           'K39': 'K2O (%)',
           'Ca44': 'CaO (%)',
           'Ti47': 'TiO2',
           'Mn55': 'MnO (%)',
           'Fe56': 'Fe2O3 (%)',
           'Co59': 'CoO',
           'Cu63': 'CuO (%)',
           'Sn118': 'SnO2',
           'Sb121': 'Sb2O3',
           'Pb208': 'PbO'}

ox_names = {el + ' oxide': ox for el, ox in ox_name.items()}

cols_to_perc = [x for x in ox_names.values() if '%' in x]


def z_score_method(df, variable_name, threshold=3):
    # Takes two parameters: dataframe & variable of interest as string
    columns = df.columns
    z = np.abs(sp.stats.zscore(df))

    outlier = []
    index = 0
    for item in range(len(columns)):
        if columns[item] == variable_name:
            index = item
    for i, v in enumerate(z[:, index]):
        if v > threshold:
            outlier.append(i)
        else:
            continue
    return outlier


def multivariate_outliers(df, cols=cols_to_perc, threshold=3):
    outliers = []
    for column in cols:
        outliers.extend(z_score_method(df, column, threshold))

    out_dict = {i: outliers.count(i) for i in range(len(df))}
    multi_out = [key for key, val in out_dict.items() if val > 7]

    return multi_out


def formatted_export(frame, path):
    writer = pd.ExcelWriter(path, engine='xlsxwriter')

    # Exporting analysis with outliers marked by red colour

    outliers = {}

    frame.replace('< LoD', np.nan, inplace=True)
    grouped = frame.groupby(frame.index)

    outliers_list = []
    for sample, group in grouped:
        if len(group) < 3:
            continue
        outliers[sample] = multivariate_outliers(group, threshold=1)

    frame.replace(np.nan, '< LoD', inplace=True)
    frame.to_excel(writer, sheet_name='Outliers')

    # formatting of excel
    workbook = writer.book
    worksheet = writer.sheets['Outliers']

    outlier_fmt = workbook.add_format({'color': 'red'})

    counter = 0
    tmp = ''
    for i, (sample, row) in enumerate(frame.iterrows()):
        if sample not in outliers.keys():
            continue

        if sample == tmp:
            counter += 1
        else:
            counter = 0
            tmp = sample
            idxs = outliers[sample]

        if counter in idxs:
            worksheet.set_row(row=i+1, cell_format=outlier_fmt)

    # exporting stats without outliers

    frame.replace('< LoD', np.nan, inplace=True)
    grouped = frame.groupby(frame.index)
    sumary = pd.DataFrame(columns=pd.concat(
        [pd.Series(['sample', 'stat']), pd.Series(frame.columns)]))
    if 'LoD' in frame.index:
        sumary = sumary.append(pd.concat(
            [pd.Series({'sample': 'LoD', 'stat': '-'}), frame.loc['LoD']]), ignore_index=True)

    # calculating stats
    for name, group in grouped:
        if name == 'LoD':
            continue

        # remove outliers
        # if name in outliers.keys():
        #     group.reset_index(inplace=True)
        #     group.drop(outliers[name], inplace=True)
        #     group.index = group['index']
        #     group.drop(['index'], axis=1, inplace=True)

        priemer = pd.concat(
            [pd.Series({'sample': name, 'stat': 'mean'}), group.mean()])
        sd = pd.concat(
            [pd.Series({'sample': name, 'stat': '2*std'}), group.std()*2])
        maxim = pd.concat(
            [pd.Series({'sample': name, 'stat': 'max'}), group.max()])
        minim = pd.concat(
            [pd.Series({'sample': name, 'stat': 'min'}), group.min()])
        med = pd.concat(
            [pd.Series({'sample': name, 'stat': 'median'}), group.median()])

        sumary = pd.concat(
            [sumary, pd.concat([priemer, sd, maxim, minim, med], axis=1).transpose()])

    sumary.replace(np.nan, '< LoD', inplace=True)

    sumary.to_excel(writer, sheet_name='Stats', index=False)

    # formatting of excel
    workbook = writer.book
    worksheet = writer.sheets['Stats']

    lod_fmt = workbook.add_format({'bold': True})
    bold_fmt = [workbook.add_format({'bg_color': '#E4FAEE', 'top': 1, 'bold': True}),
                workbook.add_format({'bg_color': '#FAF3E4', 'top': 1, 'bold': True})]
    normal_fmt = [workbook.add_format({'bg_color': '#E4FAEE', 'bold': False}),
                  workbook.add_format({'bg_color': '#FAF3E4', 'bold': False})]

    worksheet.set_row(row=1, cell_format=lod_fmt)

    tmp = False

    for i, (label, s) in enumerate(zip(sumary.index, sumary['sample'])):

        if s == 'LoD':
            continue

        if label == 0:
            worksheet.set_row(row=i+1, cell_format=bold_fmt[int(tmp)])
        else:
            worksheet.set_row(row=i+1, cell_format=normal_fmt[int(tmp)])

        if label == 4:
            tmp = not tmp

    writer.save()
