import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.optimize as scp

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
    
    def BkgSubtract_RoughShift(self,ratecorrect=True,use_uncal=False):
        '''Plots a background-subtracted plot of the first two datafiles in
        the dataclass' open signal and background files on the
        same plot.  Shifts the data set ith the lower mean by the difference
        in the means of fitted gaussians'''
        gauss = lambda x, m, s,C: C*(((s**2)*2*np.pi)**(-1/2)*np.exp(-1/2*(x-m)**2/s**2))

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
        #Fit a gaussian to signal and background data
        p0 = [1460, 7.,0.0001] #Initial guess on mean and sigma
        #Define range to fit gaussian over
        [xleft, xright] = [1450., 1470.]
        #Get indices of sigx and bkgx that are between these values
        xmin_s=np.where(sigx<xleft)[0][-1] #array index of point just below xleft
        xmax_s=np.where(sigx>xright)[0][0] #arr index of point just above xright
        xmin_b=np.where(bkgx<xleft)[0][-1] #array index of point just below xleft
        xmax_b=np.where(bkgx>xright)[0][0] #arr index of point just above xright
        #Cut ybins and xbins down to our peak ROI
        sigx = sigx[xmin_s:xmax_s]
        bkgx = bkgx[xmin_b:xmax_b]
        sigy = sigy[xmin_s:xmax_s]
        bkgy = bkgy[xmin_b:xmax_b]
        print("BKGX, STRIPPED: " + str(bkgx))
        sigyrate = sigyrate[xmin_s:xmax_s]
        bkgyrate = bkgyrate[xmin_s:xmax_s]
        sigyrate_unc = np.sqrt(sigy/(siglivetime**2))
        bkgyrate_unc = np.sqrt(bkgy/(bkglivetime**2))
        popt_s, pcov_s = scp.curve_fit(gauss, sigx,sigyrate, p0=p0,sigma=sigyrate_unc)
        popt_b, pcov_b = scp.curve_fit(gauss, bkgx,bkgyrate, p0=p0,sigma=bkgyrate_unc)
        print(popt_s)
        print(popt_b)
        sigyrate_bestfit = gauss(sigx, popt_s[0], popt_s[1],popt_s[2])
        plt.plot(sigx, sigyrate_bestfit,color='b')
        plt.plot(sigx, sigyrate, color='r')
        plt.show()
        #First, shift the x values of the signal upward
        sigx_shift = (popt_b[0] -popt_s[0])
        #Now fit again with x-shifted in the signal data
        #Now, find the number of indices in the x-axis corresponding to this shift
        firstx = sigx[0]
        ind_shift = np.where((sigx - firstx)>sigx_shift)[0][0]
        #Now, shift the signal distribution up by this many indices
        #add_shift = list(np.zeros(ind_shift))
        add_shift = list(np.zeros(1))
        sigyrate = list(add_shift) + list(sigyrate)
        sigy = list(add_shift) + list(sigy)
        del sigyrate[len(sigyrate)-len(add_shift)-1:len(sigyrate)-1]
        del sigy[len(sigy)-len(add_shift)-1:len(sigy)-1]
        sigyrate = np.array(sigyrate)
        sigy = np.array(sigy)
        sigyrate_unc = np.sqrt(sigy/(siglivetime**2))
        bkg_subtract_sigy = sigyrate - bkgyrate
        bkg_subtract_sigy_unc = np.sqrt(sigyrate_unc**2 + \
                bkgyrate_unc**2) 
        print(len(sigx))
        #bkg_subtract_sigy = sigyrate - bkgyrate
        #bkg_subtract_sigy_unc = np.sqrt((sigy/(siglivetime**2)) + (bkgy/(bkglivetime**2)))
        plt.errorbar(x=bkgx,y=bkg_subtract_sigy,yerr=bkg_subtract_sigy_unc, 
            elinewidth=4, markersize=7, color='black',marker='o',linestyle='none')
        #plt.plot(sigx, bkg_subtract_sigbestfit, markersize=7, color='black',
        #       marker='o', linestyle='none')
        if use_uncal is False:
            plt.xlabel("Energy (keV)",fontsize=26)
        else:
            plt.xlabel("ADC counts", fontsize=26)
        plt.ylabel("Counts/second",fontsize=26)
        #plt.title("Background-subtracted Watchboy PMT Bulbdown count rate",fontsize=32)
        plt.title("Difference of best fit 40K peaks \n "+\
                "(Signal data shifted by 1 bin towards background mean)",fontsize=32)
        #plt.legend(loc=1, fontsize=22)
        plt.grid(True)
        plt.show()
  
    def BkgSubtract_FDF_Peakshift(self,ratecorrect=True,use_uncal=False):
        '''Plots a background-subtracted plot of the first two datafiles in
        the dataclass' open signal and background files on the
        same plot.  Shifts the data set ith the lower mean by the difference
        in the means of fitted gaussians'''
        gauss = lambda x, m, s,C,b: C*(((s**2)*2*np.pi)**(-1/2)*np.exp((-1./2)*(x-m)**2/s**2)) + b
        #gauss = lambda x, m, s,C: C*(((s**2)*2*np.pi)**(-1/2)*np.exp(-1/2*(x-m)**2/s**2)) 
        gauss_unc = lambda x,m,s,C,m_unc,s_unc,C_unc: np.exp((-1./2)*(x-m)**2/s**2)*np.sqrt(\
                (C_unc*((s**2)*2*np.pi)**(-1/2))**2 + \
                (C*m_unc*(x-m)*((s**6)*2*np.pi)**(-1./2.))**2 +  \
                (C*s_unc*((((x-m)**2)*((s**8)*2*np.pi)**(-1./2.)) - \
                    (((s**4)*2*np.pi)**(-1/2))))**2)

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
        #Fit a gaussian to signal and background data
        p0 = [1460, 8.,0.015,0.0005] #Initial guess on mean and sigma
        #p0 = [609., 3., 0.005, 0.0015]
        #Define range to fit gaussian over
        [xleft, xright] = [1454., 1468.]
        #[xleft, xright] = [606., 614.]
        #Get indices of sigx and bkgx that are between these values
        xmin_s=np.where(sigx<xleft)[0][-1] #array index of point just below xleft
        xmax_s=np.where(sigx>xright)[0][0] #arr index of point just above xright
        xmin_b=np.where(bkgx<xleft)[0][-1] #array index of point just below xleft
        xmax_b=np.where(bkgx>xright)[0][0] #arr index of point just above xright
        #Cut ybins and xbins down to our peak ROI
        sigx = sigx[xmin_s:xmax_s]
        bkgx = bkgx[xmin_b:xmax_b]
        sigy = sigy[xmin_s:xmax_s]
        bkgy = bkgy[xmin_b:xmax_b]
        sigyrate = sigyrate[xmin_s:xmax_s]
        bkgyrate = bkgyrate[xmin_s:xmax_s]
        sigyrate_unc = np.sqrt(sigy/(siglivetime**2))
        bkgyrate_unc = np.sqrt(bkgy/(bkglivetime**2))
        popt_s, pcov_s = scp.curve_fit(gauss, sigx,sigyrate, p0=p0,sigma=sigyrate_unc)
        popt_b, pcov_b = scp.curve_fit(gauss, bkgx,bkgyrate, p0=p0,sigma=bkgyrate_unc)
        print(popt_s)
        print(popt_b)
        #First, shift the x values of the signal to match the background 
        sigx_shift = (popt_b[0] -popt_s[0])
        print("SIGX_SHIFT: " + str(sigx_shift))
        sigx_shifted = sigx - sigx_shift #This is super naughty
        #Now fit again with x-shifted in the signal data
        #popt_s, pcov_s = scp.curve_fit(gauss, sigx_shifted,sigyrate, p0=p0,sigma=sigyrate_unc)
        print("ERRORS: " + str(np.sqrt(np.diag(pcov_s))))
        muerr_s = np.sqrt(np.diag(pcov_s))[0]
        sigerr_s = np.sqrt(np.diag(pcov_s))[1]
        counterr_s = np.sqrt(np.diag(pcov_s))[2]
        #flaterr_s = np.sqrt(np.diag(pcov_s))[3]
        muerr_b = np.sqrt(np.diag(pcov_b))[0]
        sigerr_b = np.sqrt(np.diag(pcov_b))[1]
        counterr_b = np.sqrt(np.diag(pcov_b))[2]
        #flaterr_b = np.sqrt(np.diag(pcov_b))[3]
        sigyrate_bestfit = gauss(sigx_shifted, popt_s[0], popt_s[1],popt_s[2],popt_s[3])
        plt.plot(sigx, sigyrate_bestfit,color='b')
        plt.plot(sigx, sigyrate, color='r')
        plt.show()
        bkgyrate_bestfit = gauss(bkgx, popt_b[0], popt_b[1],popt_b[2],popt_b[3])
        sigyrate_bestfit_unc = gauss_unc(sigx_shifted, popt_s[0], popt_s[1], popt_s[2],
                muerr_s, sigerr_s,counterr_s)
        plt.errorbar(x=sigx_shifted, y=sigyrate_bestfit, yerr=sigyrate_bestfit_unc,
                linestyle='none', marker='o')
        plt.show()
        bkgyrate_bestfit_unc = gauss_unc(bkgx, popt_b[0], popt_b[1], popt_b[2],
                muerr_b, sigerr_b,counterr_b)
        print("BEST FIT SIGNAL MU AND SIG: " + str(popt_s))
        print("BEST FIT BKG MU AND SIG: " + str(popt_b))
        #Shift mean of signal curve upward (it's the lower one)
        bkg_subtract_sigbestfit = sigyrate_bestfit - bkgyrate_bestfit 
        bkg_subtract_sigbestfit_unc = np.sqrt(sigyrate_bestfit_unc**2 + \
                bkgyrate_bestfit_unc**2) 
        print(len(sigx))
        print(len(bkg_subtract_sigbestfit))
        #bkg_subtract_sigy = sigyrate - bkgyrate
        #bkg_subtract_sigy_unc = np.sqrt((sigy/(siglivetime**2)) + (bkgy/(bkglivetime**2)))
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.errorbar(x=sigx,y=bkg_subtract_sigbestfit,yerr=bkg_subtract_sigbestfit_unc, 
            elinewidth=4, markersize=7, color='black',marker='o',linestyle='none')
        #plt.plot(sigx, bkg_subtract_sigbestfit, markersize=7, color='black',
        #       marker='o', linestyle='none')
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(16)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(16)
        if use_uncal is False:
            plt.xlabel("Energy (keV)",fontsize=26)
        else:
            plt.xlabel("ADC counts", fontsize=26)
        plt.ylabel("Counts/second",fontsize=26)
        #plt.title("Background-subtracted Watchboy PMT Bulbdown count rate",fontsize=32)
        plt.title("Difference of best fit signal & background 40K peak (Signal-shifted)\n"+\
                "Watchboy PMT Bulbdown count rate",fontsize=32)
        #plt.legend(loc=1, fontsize=22)
        plt.grid(True)
        plt.show()
  
