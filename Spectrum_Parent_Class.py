# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 12:08:45 2019

@author: Molly Lockart / Joe Butler
"""

##############################################################################
##############################################################################
#Imports

import os

try:
    import Tkinter as tk
except ImportError:
    import tkinter as tk

#from tkinter import messagebox
from tkinter import filedialog
#from tkinter import simpledialog
#from tkinter import Label, Frame, OptionMenu, StringVar
#from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.colors import LinearSegmentedColormap
#from scipy.optimize import minimize
#from scipy.signal import savgol_filter

# Our processing classes
from MimsENDOR import MimsENDOR
from HYSCORE import HYSCORE

##############################################################################
##############################################################################


class Spectrum:
    """Parent class for EPR spectrum read in from file"""

    def __init__(self, full_file_name):
        """initializes a spectrum from file based on metadata file or a dialog;
        determines child class, reads and does initial processing"""

        self.filename = full_file_name
        # variable full_file_name comes from FILENAMELIST[j] in executed code
        (self.path, self.file) = os.path.split(self.filename)
        self.metadata_file = self.filename.replace(".DTA", ".DSC")

        def elexsys_dict(self):
            """builds a searchable dictionary out of Elexsys .DSC file.
            only gets experimental parameters, ignores pulse programs."""

            bruker_file = open(self.metadata_file, 'r')
            edict = {}
            edict['FullFileName'] = self.filename  # file for metadata

            lookup = 'PlsSPELPrg'
            # the last line that needs to be in dictionary; the pulse sequence starts

            for line in bruker_file:
                line = line.strip()

                if lookup in line:
                    break

                else:
                    if ("#" not in line and ".DVC" not in line and "begin" not in line
                            and "end" not in line and not line.startswith("\\n\\")
                            and not line.startswith("\\n") and not line.startswith("*")
                            and not line.startswith(";") and line != ''):

                        if '=' in line:
                            splitline = line.partition('=')
                            edict[splitline[0]] = splitline[2]
                        else:
                            splitline = line.split(maxsplit=1)
                            if len(splitline) == 2:
                                edict[splitline[0]] = splitline[1]

            return edict

        def build_metadata(self):
            """determine type of spectrum and build dictionary of
            spectral metadata"""

            (spectrumname, ext) = os.path.splitext(os.path.split(self.file)[1])
            if ext == '.DTA':
                # maybe Bruker Elexsys file
                if os.path.isfile(self.filename.replace(".DTA", ".DSC")):
                    # assign type and build metadata dictionary
                    file_type = 'Elexsys'
                    metadata = elexsys_dict(self)
                    metadata['SpecName'] = spectrumname
                    metadata['FileType'] = file_type
                else:
                    file_type = "Unknown"
                    metadata = {}
                    metadata['FileType'] = file_type
                    # query for metadata
            else:
                file_type = "Unknown"
                metadata = {}
                metadata['FileType'] = file_type
                # query for metadata

            exp = metadata.get('PlsSPELEXPSlct', "Unknown")

            # now select child class of Spectrum
            if 'Mims ENDOR ESE' or 'Mims ENDOR' in exp:
                # an ENDOR spectrum
                spectype = "ENDOR"
            elif 'HYSCORE' or 'HYSCORE Split' in exp:
                # HYSCORE spectrum
                spectype = "HYSCORE"
            elif exp == 'Log T1 Recovery':
                # T1e recovery, on log scale, integrated echo
                spectype = "log T1e"
            else:
                # unknown
                spectype = "Unknown"
                print("Pulse sequence not recognised")

            return spectype, metadata

##############################################################################
# initialization varaibles

        # Creates a searchable dictionary with the experimental parameters
        (self.spec_type, self.metadata) = build_metadata(self)

        # finds the correct child class and uses it to process the data
        if self.spec_type == "ENDOR":
            self.spect = MimsENDOR(self.metadata, self.filename)  # contains Amps and Coords
            # self.spect = MimsENDOR(amp, x_axis, x_unit)  # contains Amps and Coords

            # this is the part that doesn't work so far. self.amp is the raw y
            # data that needs to be processed

        elif self.spec_type == "HYSCORE":
            self.spect = HYSCORE(self.metadata, self.ffilename)
            #newest attempted addition, untested as of 1/31/19 am

        elif self.spec_type == "Log T1e":
            pass

        else:
            print("Unrecognised type of spectrum: ", self.spec_type)



##############################################################################
# Class Methods

    def lineplot(self):
        """makes a line plot from x and y data and saves it as a .tif figure"""

        plt.figure()
        #plt.plot(self.spect.calc_endorfreq(), self.spect.smoothed, 'k', linewidth=1.5)
        # plt.plot(self.spect.return_axis('x'), self.spect.phased_complex.imag, 'r', linewidth=1.5)
        plt.plot(self.spect.return_axis('x'), self.spect.phased, 'b', linewidth=1.5)

        # plt.xlim([self.spect.return_axis('x')[0], self.spect.return_axis('x')[-1]])
        plt.xlabel('ENDOR Frequency ('+self.spect.return_units('x')+')', fontsize=14)
        plt.ylabel('Absolute ENDOR Effect', fontsize=14)
        plt.tick_params(labelsize=12)
        save = [self.filename.replace(".DTA", ".tiff")]
        plt.savefig("".join(save), dpi=100, bbox_inches='tight', format="tiff")

# SpecData.return_xaxis
#############################################################################
#############################################################################
# Executed Code

def choose_files():
    """choose files from file dialog box;
    either the Bruker Elexsys .DTA;
        or any other type"""

    root = tk.Tk()
    root.withdraw()
    filenamelist = filedialog.askopenfilenames(
        parent=root,
        defaultextension='.DTA',
        title="Select Data files",
        filetypes=[('Bruker', '.DTA'), ('all', '.*')])
    return filenamelist


# Test what is here so far
FILENAMELIST = choose_files()
THE_SPECTRUM = list()

for j in range(0, len(FILENAMELIST)):
    THE_SPECTRUM.append(Spectrum(FILENAMELIST[j]))
    THE_SPECTRUM[j].lineplot()
