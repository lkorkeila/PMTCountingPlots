import os,sys
import glob
import plots.Plots as p
import lib.datastruct as ds

#Define your experiment data and background counting data file locations
DEBUG=True

MAINDIR = os.path.dirname(__file__)
DATADIR = os.path.abspath(os.path.join(MAINDIR, "data", "signal"))
BKGDATADIR = os.path.abspath(os.path.join(MAINDIR, "data","background"))

if __name__=='__main__':
    print("MAKING PMT COUNTING STACK PLOTS")

    #Get list of all signal and background data
    listofdatafiles = glob.glob(DATADIR + '/*.npz')
    listofbkgfiles = glob.glob(BKGDATADIR + '/*.npz')

    if DEBUG is True:
        print("PRINTING LIST OF DATA FILES, THEN BACKGROUND FILES")
        print(listofdatafiles)
        print(listofbkgfiles)

    # First, initialize your data structure with all desired data file
    # and background file locations
    alldatafiles = ds.datacollection(listofdatafiles, listofbkgfiles)

    # Load the current list of data files into data to be analyzed
    alldatafiles.load_datafiles()
    #Loads in bkg files for use. Uncomment to calculate bkgs w/ peak sidebands
    alldatafiles.load_bkgdatafiles()

    ourplotter = p.DataPlotter(dataclass=alldatafiles)
    ourplotter.StackPlot_FirstDataFiles()
