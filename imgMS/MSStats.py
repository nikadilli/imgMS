from imgMS import MSData
from imgMS import MSEval
from imgMS.side_functions import *

import numpy as np
import xlsxwriter
import io
import os

from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
import matplotlib.pyplot as plt

import itertools

stats_dict = {'mean':None, 'std':None, 'sum':None, 'min':None, 'max':None, 'med':None,}

class interactive_average():
    """
    Creates interactive graph for the selection of the area of interest and
    calculates basic stats as mean, std, median, min and max. The values can
    be printed or saved as an excel file.  

    Parameters
    ----------
    msdata: MSData object from imgMS
        Data of elemental imaging where stats will be calculated from a selected area of the elemental distribution image. 
    """
    def __init__(self, msdata, elem=None):
        
        if elem is None: 
            self.elem = msdata.isotope_names[0]
        else: 
            if elem in msdata.isotope_names:
                self.elem = elem
            else:
                print(f'{elem} not in isotopes')
        self.msdata = msdata
        self.elmap = msdata.isotopes[self.elem].elmap
        self.get_attr()
        self.stats_dict = { e:{'mean':None, 'std':None, 'sum':None, 'min':None, 'max':None, 'med':None,} for e in self.msdata.isotope_names }
        

    def get_attr(self):
        """
        get all the data required for interactive mean calculation for given map
        """
        self.matrix = self.elmap.matrix
        self.shape = self.matrix.shape
        self.dx = self.elmap.dx
        self.dy = self.elmap.dy
        
        self.pixx = np.arange(self.shape[1])
        self.pixy = np.arange(self.shape[0])
        self.xv, self.yv = np.meshgrid(self.pixx,self.pixy)
        self.pix = np.vstack( (self.xv.flatten(), self.yv.flatten()) ).T

        self.array = np.array(self.matrix)
        self.array[np.isnan(self.array)] = 0
        
    def switch_elem(self, elem):
        """
        Change the element for interactive calculation
        
        Parameters
        ----------
        elem: str
        Name of the element, which will be used for selecting the boundaries of area where stats are calculated.
        """
        self.elem = elem
        self.elmap = self.msdata.isotopes[self.elem].elmap
        self.get_attr()
        
    def updateArray(self, array, indices):
        lin = np.arange(array.size)
        newArray = np.array(array.flatten())
        newArray[lin[~indices]] = 0.25
        newArray[lin[indices]] = 1
        return newArray.reshape(array.shape)
    
    def stats(self, indices):
        """
        Calculates the statistical values (mean, std, sum, minimum, 
        maximum and median) from selected area of image.
        
        Parameters
        ----------
        indices: list
        List of bool to select which pixels are in the area.
        """
        for e in self.msdata.isotopes:
            array = data.isotopes[e].elmap.matrix
            lin = np.arange(array.size)
            newArray = np.array(array.flatten())
            mean = newArray[lin[indices]].mean()
            std = newArray[lin[indices]].std()
            suma = newArray[lin[indices]].sum()
            mina = newArray[lin[indices]].min()
            maxa = newArray[lin[indices]].max()
            meda = np.median(newArray[lin[indices]])
            
            self.stats_dict[e]['mean'] = newArray[lin[indices]].mean()
            self.stats_dict[e]['std'] = newArray[lin[indices]].std()
            self.stats_dict[e]['sum'] = newArray[lin[indices]].sum()
            self.stats_dict[e]['min'] = newArray[lin[indices]].min()
            self.stats_dict[e]['max'] = newArray[lin[indices]].max()
            self.stats_dict[e]['med'] = np.median(newArray[lin[indices]])
            
            print(e)
            print( f'mean: {mean}')
            print( f'std: {std}')
            print( f'sum: {suma}')
            print( f'min: {mina}')
            print( f'max: {maxa}')
            print( f'med: {meda}')

    def onSelect(self, verts):
        """
        Selection of area in the image.
        """
        p = Path(verts)
        ind = p.contains_points(self.pix, radius=1)
        alpha = self.updateArray(self.array, ind)
        self.im.set_alpha(alpha)
        self.fig.canvas.draw_idle()
        self.stats(ind)

    def __call__(self, ax=None, vmax=None):
        self.ax = ax
        if self.ax == None:
            self.fig, self.ax = plt.subplots()

        self.im = self.ax.imshow(self.array, vmax=vmax)
        self.ax.set_xlim([0, self.shape[1]])
        self.ax.set_ylim([0, self.shape[0]])

        lsso = LassoSelector(ax=self.ax, onselect=self.onSelect)

        plt.show()

    def export_stats(self, filename, sheetname='Sheet1'):
        """
        Export calculated stat for every isotope in the selected area 
        into excel.
        
        Parameters
        ----------
        filename: str
        Excel file name to save data.
        """
        workbook = xlsxwriter.Workbook(filename)
            
        wks1=workbook.add_worksheet(sheetname)
        
        wks1.write(0,0,'elem')
        wks1.write(1,0,'mean')
        wks1.write(2,0,'std')
        wks1.write(3,0,'sum')
        wks1.write(4,0,'min')
        wks1.write(5,0,'max')
        wks1.write(6,0,'med')
        
        for i, e in enumerate(self.msdata.isotopes):
            
            wks1.write(0,i+1,e)
            wks1.write(1,i+1,self.stats_dict[e]['mean'])
            wks1.write(2,i+1,self.stats_dict[e]['std'])
            wks1.write(3,i+1,self.stats_dict[e]['sum'])
            wks1.write(4,i+1,self.stats_dict[e]['min'])
            wks1.write(5,i+1,self.stats_dict[e]['max'])
            wks1.write(6,i+1,self.stats_dict[e]['med'])

        imgdata=io.BytesIO()
        self.fig.savefig(imgdata, format='png')
        wks1.insert_image(8,3, '', {'image_data': imgdata})

        workbook.close()
        






















        
