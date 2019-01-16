# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 12:08:45 2019

@author: Student
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
import MimsENDOR

##############################################################################
##############################################################################
#Class definition

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

        def read_bruker_file(self):
            """Reads in Bruker .DTA file containing only the y values (1D)
            or the z values (if 2D experiment) and x,y axes are created"""

            # Read in binary data as Elexsys floats
            rawbin = np.fromfile(self.filename, dtype='>f8')
                #big endian, f8 is 64-bit floating point number

            if self.complex:
                # Complex data as pairs of real, imag numbers to be combined
                size = int(rawbin.size / 2)
                if size != int(self.metadata['XPTS']):
                    print('Size of data is inconsistent')
                amp_data = np.zeros((size, 1), dtype=complex)
                amp_data = rawbin[0::2]+1j* rawbin[1::2]
            else:
                amp_data = rawbin

            if not self.onedim:
                # Not 1D data; assume 2D square as for HYSCORE
                dimx = int(self.metadata['XPTS'])
                dimy = int(self.metadata['YPTS'])
                amp_data = np.reshape(amp_data, (dimx, dimy))

            return amp_data

        def buildaxis(self, axis):
            """builds the x axis. Gets the number of points, the x minimum, and
            the x maximum from the dictionary created at the beginning of the
            script. Uses those values to create the x axis."""

            axdim = int(self.metadata[axis+'PTS'])
            axmin = float(self.metadata[axis+'MIN'])
            axrange = float(self.metadata[axis+'WID'])
            # xstep = xrange / (xdim_pad - 1)
            ax_n = (np.linspace(
                axmin, axmin + axrange, num=axdim, endpoint=False))
            ax_unit = self.metadata[axis+'UNI']

            return ax_n, ax_unit
##############################################################################
# initialization varaibles

        # Creates a searchable dictionary with the experimental parameters
        (self.spec_type, self.metadata) = build_metadata(self)

        # Decides what kind of data it is
        if self.metadata['IKKF'] == 'CPLX':
            self.complex = True
        else:
            self.complex = False

        if self.metadata['YTYP'] == 'NODATA':
            self.onedim = True
        else:
            self.onedim = False

        # Reads the raw amplitude data (either 1D or 2D)
        self.amp = read_bruker_file(self)

        # Builds the x axis if one dimensional (ENDOR)
        (self.xaxis, self.xunit) = buildaxis(self, 'X')

        # Builds the y axis too if two dimsntional (HYSCORE)
        if not self.onedim:
            (self.yaxis, self.yunit) = buildaxis(self, 'Y')

        # finds the correct child class and uses it to process the data
        if self.spec_type == "ENDOR":
            self.spect = MimsENDOR  # contains Amps and Coords

        elif self.spec_type == "HYSCORE":
            pass

        elif self.spec_type == "Log T1e":
            pass

        else:
            print("Unrecognised type of spectrum: ", self.spec_type)
##############################################################################
# Class Methods

    def lineplot(self):
        """makes a line plot from x and y data and saves it as a .tif figure"""

        plt.figure()
        plt.plot(self.xaxis, self.amp.real, 'k', linewidth=1.5)
        plt.xlim([self.xaxis[0], self.xaxis[-1]])
        plt.xlabel('RF Frequency ('+self.xunit+')', fontsize=14)
        plt.ylabel('Absolute ENDOR Effect', fontsize=14)
        plt.tick_params(labelsize=12)
        save = [self.filename.replace(".DTA", ".tiff")]
        plt.savefig("".join(save), dpi=100, bbox_inches='tight', format="tiff")


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
