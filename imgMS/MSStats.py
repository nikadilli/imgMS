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

class interactive_average():
    """
    Creates interactive graph for the selection of the area of interest and
    calculates basic stats as mean, std, median, min and max. The values can
    be printed or saved as an excel file.  

    Parameters
    ----------
    elmap: elmap class from imgMS
        Data of elemental image of one isotope stored as numpy array. 
    """
    def __init__(self, elmap):
        #TODO add possibility to pass multiple images / eg all isotopes of a
        # map at once 
        self.elmap = elmap

        self.get_attr()
        

    def get_attr(self):
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
        
    def updateArray(self, array, indices):
        lin = np.arange(array.size)
        newArray = np.array(array.flatten())
        newArray[lin[~indices]] = 0.25
        newArray[lin[indices]] = 1
        return newArray.reshape(array.shape)
    
    def stats(self, array, indices):
        lin = np.arange(array.size)
        newArray = np.array(array.flatten())
        self.mean = newArray[lin[indices]].mean()
        self.std = newArray[lin[indices]].std()
        self.sum = newArray[lin[indices]].sum()
        self.min = newArray[lin[indices]].min()
        self.max = newArray[lin[indices]].max()
        self.med = np.median(newArray[lin[indices]])
        print( f'mean: {self.mean}')
        print( f'std: {self.std}')
        print( f'sum: {self.sum}')
        print( f'min: {self.min}')
        print( f'max: {self.max}')
        print( f'med: {self.med}')

    def onSelect(self, verts):
        p = Path(verts)
        ind = p.contains_points(self.pix, radius=1)
        alpha = self.updateArray(self.array, ind)
        self.im.set_alpha(alpha)
        self.fig.canvas.draw_idle()
        self.stats(self.array, ind)

    def __call__(self, vmax=None):
        self.fig, self.ax = plt.subplots()

        self.im = self.ax.imshow(self.array, vmax=vmax)
        self.ax.set_xlim([0, self.shape[1]])
        self.ax.set_ylim([0, self.shape[0]])

        lsso = LassoSelector(ax=self.ax, onselect=self.onSelect)

        plt.show()

    def export_stats(self, filename, sheetname='Sheet1'):
        
        workbook = xlsxwriter.Workbook(filename)
            
        wks1=workbook.add_worksheet(sheetname)

        wks1.write(0,0,'mean')
        wks1.write(1,0,'std')
        wks1.write(2,0,'sum')
        wks1.write(3,0,'min')
        wks1.write(4,0,'max')
        wks1.write(5,0,'med')

        wks1.write(0,1,self.mean)
        wks1.write(1,1,self.std)
        wks1.write(2,1,self.sum)
        wks1.write(3,1,self.min)
        wks1.write(4,1,self.max)
        wks1.write(5,1,self.med)

        imgdata=io.BytesIO()
        self.fig.savefig(imgdata, format='png')
        wks1.insert_image(0,3, '', {'image_data': imgdata})

        workbook.close()




















        
