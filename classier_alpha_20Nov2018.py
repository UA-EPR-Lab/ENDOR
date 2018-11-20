# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 13:16:19 2018

@author: Molly Lockart
"""
# for data import
import os, sys

try:
    import Tkinter as tk
except ImportError:
    import tkinter as tk

from tkinter import messagebox
from tkinter import filedialog
from tkinter import simpledialog

# for data processing
import numpy as np
import matplotlib.pyplot as plt
# from scipy import interpolate
from scipy.optimize import minimize
from scipy.signal import savgol_filter

# Parent class; will add subclasses for processing Mims ENDOR, HYSCORE, etc.

def choose_files():

    """choose files from file dialog box; can choose either the .DSC or the .DTA"""

    root = tk.Tk()
    root.withdraw()
#    messagebox.showinfo("For ENDOR Subtractions",
#                        "First choose file, then choose file to subtract from first file.")
    filenamelist = filedialog.askopenfilenames(parent=root, defaultextension='.DTA',
                                               title="Select Data files",
                                               filetypes=[('Bruker', '.DTA'),
                                                          ('all', '.*')])


    return filenamelist


class ProcessSpectra:

    """this is the parent class. It will choose a file, crease a file dictionar
    y,and phase the raw data"""

    def __init__(self, file_dictionary=0, rawdata=0, phased=0, phase_angle=0):

        """initializes the class Process_Spectra. This will load a file, create a dictionary, and
        read the .DTA file every time this class is invoked."""

        self.filename = filename

        def create_dict(self):

            """this creates a searchable dictionary out of the description file.
            it only takes experimental parameters; it ignores the pulse programs
            at the end of the description file."""


            if self.filename.endswith(".DSC"):
                self.dict_filename = self.filename

            else:

                try:
                    self.dict_filename = filename.replace(".DTA", ".DSC")

                except ValueError:
                    messagebox.showinfo('there is no valid .DSC file')

            dictionary = {}

            bruker_file = open(self.dict_filename, 'r')

            for line in bruker_file:
                line = line.strip()
                lookup = 'PlsSPELPrg'  # the last line that needs to be in the
                # dictionary the pulse sequence starts

                if lookup in line:
                    break

                else:
                    if ("#" not in line and ".DVC" not in line
                            and "begin" not in line and "end" not in line
                            and not line.startswith("\\n\\")
                            and not line.startswith("\\n")
                            and not line.startswith("*")
                            and not line.startswith(";") and line != ''):

                        line = line.split()

                        if "=" in line:
                            dictionary[line[0]] = line[2]

                        else:
                            dictionary[line[0]] = line[1:]

            return dictionary

        def read_dta_file(self):

            """this reads in the data file. The Bruker .DTA file only contains the
            y values (if one-dimensional experiment) or the z values (if it's a two
            dimensional experiment). The x axis is created later"""

            rawdata = np.fromfile(self.filename, dtype='>f8')

            return rawdata

        def phase(self):

            """this phases the data by minimizing the imaginary component of the
            data."""

            realdata = np.ones([int((len(self.rawdata)) / 2)])
            imagdata = np.ones([int((len(self.rawdata)) / 2)])

            for i in range(0, int((len(self.rawdata)) / 2)):
                realdata[i] = self.rawdata[i * 2]
                imagdata[i] = self.rawdata[i * 2 + 1]
                complexdata = realdata + 1j * imagdata

            def phase(phi):

                """phase function"""

                return np.sum(abs((np.sin(phi) * complexdata.real -
                                   np.cos(phi) * complexdata.imag)))

            start_pos = np.array(0)
            res = minimize(phase, start_pos, method='Powell')
            phi = res.x
            complexdataraw = np.exp(-phi * 1j) * complexdata
            phaseddata = complexdataraw.real

            return phaseddata, phi

        self.file_dictionary = create_dict(self)
        self.rawdata = read_dta_file(self)
        self.phased = phase(self)[0]
        self.phase_angle = phase(self)[1]


    def get_from_dict(self, key):

        """this gets something from the dictionary. the key input is
        what you want to pull from the dictionary. Ex: get_from_dict('XPTS')
        returns the number of points in the measurement."""

        value = self.file_dictionary[key]
        if (key != 'XPTS') and (key != 'XMIN') and (key != 'XWID'):
            value = " ".join(value) #takes string out of tuple so you can int() it
            value = value.strip("'")
        return value


    def lineplot(self):

        """makes a line plot from x and y data and saves it as a .tif figure"""

        plt.figure()
        plt.plot(spectrum.xdata, spectrum.data, 'k', linewidth=1.5)
        xmin = float("".join(ProcessSpectra.get_from_dict(self, 'XMIN')))
        xwid = float("".join(ProcessSpectra.get_from_dict(self, 'XWID')))
        plt.xlim([xmin, xmin + xwid])
        plt.xlabel('RF Frequency (MHz)', fontsize=14)
        plt.ylabel('Absolute ENDOR Effect', fontsize=14)
        plt.xticks(np.arange(xmin, xmin + xwid + 1, step=1))
        plt.tick_params(labelsize=12)
        save = [spectrum.filename, ".tiff"]
        plt.savefig("".join(save), dpi=100, bbox_inches='tight', format="tiff")




class MimsENDOR(ProcessSpectra):

    """This is a subclass that processes Mims ENDOR spectra."""

    def __init__(self):
        super().__init__(self)

        self.lndata = np.log1p(self.phased)

        def baseline_correct(self):

            """baseline corrects the data by fitting it to a zero-order polynomial.
            It chooses the maximum 10% of y values to fit, which won't contain
            any signal. This eliminates the need to choose regions to fit the
            poolynomial to."""

            basex = np.arange(0, len(self.lndata))
            lndata1 = np.transpose(self.lndata)
            indexarray = np.transpose(np.argsort(lndata1))
            ind = len(indexarray)
            baselineind = indexarray[int(.9 * ind):ind - 1]
            polyx1 = np.arange(0, len(baselineind))
            polyy1 = np.arange(0, len(baselineind))
            polyx1 = basex[baselineind]
            polyy1 = lndata1[baselineind]
            polynomial = np.polyfit(polyx1, polyy1, 1)
            baseline = np.polyval(polynomial, basex)
            baseline_corrected = self.lndata - baseline

            return baseline_corrected

        def exp(self):

            """takes the exp of the data"""

            expdata = np.exp(self.baseline_corrected)

            return expdata

        def flipdata(self):

            """flip data to get absolute ENDOR effect"""

            flipdata = np.subtract(self.expdata, 1)
            flipdata = flipdata * -1

            return flipdata

        def smooth(self):

            """smoothes using a Savitsky-Golay filter. the default is to fit a 4th
            order polynomial over an odd number of points. this can be changed
            depending on how much you want to smooth. Increase the number of points
            to smooth."""

            smoothed = savgol_filter(self.processed, 25, 6)

            return smoothed

        def buildxy(self):

            """builds the x axis. Gets the number of points, the x minimum, and
            the x maximum from the dictionary created at the beginning of the
            script. Uses those values to create the x axis."""

            xdim = float("".join(ProcessSpectra.get_from_dict(self, 'XPTS')))
            # xdim_string = str(xdim)
            # pad = 0
            # xdim_pad = np.pad(xdim, (pad, pad), 'constant')
            # return xdim_pad
            xmin = float("".join(ProcessSpectra.get_from_dict(self, 'XMIN')))
            xrange = float("".join(ProcessSpectra.get_from_dict(self, 'XWID')))
            # xstep = xrange / (xdim_pad - 1)
            freqx_n = (np.linspace(xmin, xmin + xrange, num=xdim, endpoint=False))

            return freqx_n

        self.baseline_corrected = baseline_correct(self)
        self.expdata = exp(self)
        self.processed = flipdata(self)
        self.data = smooth(self)
        self.xdata = buildxy(self)

    def calc_endorfreq(self):

        """calculates the ENDOR frequency for a proton. Uses the dictionary
        created in the beginning of the script to find the magnetic field to
        use in the calulation. Returns a frequency x-axis that is centered
        at the nucelar ENDOR frequency."""

        try:
            b0vl = float("".join(ProcessSpectra.get_from_dict(self, 'B0VL')))
        except KeyError:
            root = tk.Tk()
            root.withdraw()
            try:
                b0vl = float(
                    simpledialog.askstring(
                        "Bruker Error",
                        "Enter Magnetic Field in Gauss",
                        parent=root)
                ) / 10000  # enter field in Gauss and convert to Tesla
            except NameError:
                endorfreqx = spectrum.xdata
                print("Magnetic field should be a number")
                return endorfreqx
        endorfreq = (300 * b0vl) / 7.046
        endorfreqx = (endorfreq - spectrum.xdata) * (-1)

        return endorfreqx

    def plot(self):

        """inherits the lineplot method from ProcessSpectra"""

        ProcessSpectra.lineplot(self)

FILENAMELIST = choose_files()



for j in range(0, len(FILENAMELIST)):
    filename = FILENAMELIST[j]
    spectrum = MimsENDOR()
    figure = MimsENDOR.plot(spectrum)
    spectrumnam = os.path.split(filename)[1]
    plt.title(spectrumname)
    np.savetxt(filename.replace(".DTA", ".txt"),
               np.c_[spectrum.xdata, spectrum.data], delimiter=",")
