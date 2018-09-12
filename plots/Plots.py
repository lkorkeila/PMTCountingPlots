import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

class DataPlotter(object):
    def __init__(self, dataclass=None):
        self.dataclass = dataclass

    def StackPlot_FirstDataFiles(self,ratecorrect=True):
        '''Plots a stack plot of the first two datafiles in
        the dataclass' open signal and background files on the
        same plot'''
        sigx = np.array(self.dataclass.opendatafiles[0].x)
        sigy = np.array(self.dataclass.opendatafiles[0].y)
        bkgx = np.array(self.dataclass.openbkgfiles[0].x)
        bkgy = np.array(self.dataclass.openbkgfiles[0].y)
        if ratecorrect==True:
            sigrate = self.dataclass.opendatafiles[0].livetime
            bkgrate = self.dataclass.openbkgfiles[0].livetime
            sigy = sigy/sigrate
            bkgy = bkgy/bkgrate
        bg_subtract_sigy = sigy - bkgy
        plt.plot(sigx, sigy, color='b', alpha=0.7, linewidth=4,label='PMT bulbdown')
        plt.plot(bkgx, bkgy, color='g', alpha=0.7, linewidth=4,label='Background')
        plt.xlabel("Energy (keV)",fontsize=26)
        plt.ylabel("Counts/second",fontsize=26)
        plt.title("Comparison of Watchboy PMT and background radioactivity",fontsize=32)
        plt.legend(loc=1, fontsize=22)
        plt.show()

    def BkgSubtract_FirstDataFiles(self,ratecorrect=True,use_uncal=False):
        '''Plots a background-subtracted plot of the first two datafiles in
        the dataclass' open signal and background files on the
        same plot'''
        sigx, sigy = [], []
        bkgx, bkgy = [], []
        if use_uncal is False: 
            bkgx = np.array(self.dataclass.openbkgfiles[0].x)
            sigx = np.array(self.dataclass.opendatafiles[0].x)
        else:
            bkgx = np.array(self.dataclass.openbkgfiles[0].x_uncal)
            sigx = np.array(self.dataclass.opendatafiles[0].x_uncal)
        sigy = np.array(self.dataclass.opendatafiles[0].y)
        bkgy = np.array(self.dataclass.openbkgfiles[0].y)
        sigrate = None
        sigrate = None
        if ratecorrect==True:
            siglivetime = self.dataclass.opendatafiles[0].livetime
            bkglivetime = self.dataclass.openbkgfiles[0].livetime
            sigyrate = sigy/siglivetime
            bkgyrate = bkgy/bkglivetime
        bkg_subtract_sigy = sigyrate - bkgyrate
        bkg_subtract_sigy_unc = np.sqrt((sigy/(siglivetime**2)) + (bkgy/(bkglivetime**2)))
        plt.errorbar(x=sigx,y=bkg_subtract_sigy,yerr=bkg_subtract_sigy_unc, 
            elinewidth=4, markersize=7, color='black',marker='o',linestyle='none')
        if use_uncal is False:
            plt.xlabel("Energy (keV)",fontsize=26)
        else:
            plt.xlabel("ADC counts", fontsize=26)
        plt.ylabel("Counts/second",fontsize=26)
        plt.title("Background-subtracted Watchboy PMT Bulbdown count rate",fontsize=32)
        plt.legend(loc=1, fontsize=22)
        plt.grid(True)
        plt.show()

    
