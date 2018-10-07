# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 10:32:01 2018

@author: Others
"""

"""processes ENDOR files from the raw data file. Takes Bruker .DTA files and
processes them."""

# %matplotlib tk
# this is necessary to get the tkinter dialog box to close on a Mac.

# Import all necessary packages
import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.signal import savgol_filter

# Choose Files
root = tk.Tk()
messagebox.showinfo("For ENDOR Subtractions","First choose file, then choose file to subtract from first file.")
        #root = tk.Tk()
filenames = filedialog.askopenfilenames(parent=root)
filenamelist=[0,len(filenames)]
for i in range(0,len(filenames)):
    filenamelist[i]="".join(filenames[i])
    if filenamelist[i].endswith('.DTA'):
        filenamelist[i] = filenamelist[i]
    if not filenamelist[i].endswith('.DTA'):
        if filenamelist[i].endswith('.DSC'):
            filenamelist[i] = filenamelist[i].replace('.DSC', '.DTA')
        else:
            messagebox.showinfo("Error", "must be a .DTA file")

# Create data arrays for output
processed=[0,len(filenames)]
endorfreqx=[0,len(filenames)]
freqx=[0,len(filenames)]

# Parent function
def endor_process():
    """parent function; this is all you need to call to process the spectra"""

    def create_dict(filename):
        """this creates a searchable dictionary out of the description file.
        it only takes experimental parameters; it ignores the pulse programs
        at the end of the description file."""
        if filename.endswith(".DSC"):
            filename = filename
        else:
            try:
                filename = filename.replace(".DTA", ".DSC")  # just in case the
    # .DTA file is chosen instead
            except:
                messagebox.showinfo('this is not a valid .DSC file')

        dictionary = {}

        bruker_file = open(filename, 'r')
        for line in bruker_file:
            line = line.strip()
            lookup = 'PlsSPELPrg'  # the last line that needs to be in the
    # dictionary the pulse sequence starts

            if lookup in line:
                break
            else:
                if ("#" not in line and ".DVC" not in line and "begin" not in line
                        and "end" not in line and not line.startswith("\\n\\")
                        and not line.startswith("\\n") and not line.startswith("*")
                        and not line.startswith(";") and line != ''):
                    line = line.split()
                    if "=" in line:
                        dictionary[line[0]] = line[2]
                    else:
                        dictionary[line[0]] = line[1:]
        return dictionary

    def get_from_dict(key):
        """this gets something from the dictionary. the key input is
        what you want to pull from the dictionary. Ex: get_from_dict('XPTS')
        returns the number of points in the measurement."""
        value = file_dictionary[key]
        if (key != 'XPTS') and (key != 'XMIN') and (key != 'XWID'):
            value = " ".join(value)
            value = value.strip("'")
        return value


    def read_dta_file(filename):
        """this reads in the data file. The Bruker .DTA file only contains the
        y values (if one-dimensional experiment) or the z values (if it's a two
        dimensional experiment). The x axis is created later"""
        data = np.fromfile(filename, dtype='>f8')
        return data


    def phase(data):
        """this phases the data by minimizing the imaginary component of the
        data."""
        realdata = np.ones([int((len(data))/2)])
        imagdata = np.ones([int((len(data))/2)])
        for i in range(0, int((len(data))/2)):
            realdata[i] = data[i*2]
            imagdata[i] = data[i*2+1]
        complexdata = realdata+1j*imagdata

        def phase(phi):
            """phase function"""
            return np.sum(abs((np.sin(phi) * complexdata.real -
                               np.cos(phi) * complexdata.imag)))
        start_pos = np.array(0)
        res = minimize(phase, start_pos, method='Powell')
        phi = res.x
        complexdataraw = np.exp(-phi*1j)*complexdata
        phaseddata = complexdataraw.real
        # imagdataraw=complexdataraw.imag
        return phaseddata


    def naturallog(phaseddata):
        """takes ln of the data"""
        lndata = np.log1p(phaseddata)
        return lndata


    def baseline_correct(lndata):
        """baseline corrects the data by fitting it to a zero-order polynomial.
        It chooses the maximum 10% of y values to fit, which won't contain
        any signal. This eliminates the need to choose regions to fit the
        poolynomial to."""
        basex = np.arange(0, len(lndata))
        lndata1 = np.transpose(lndata)
        indexarray = np.transpose(np.argsort(lndata1))
        ind = len(indexarray)
        baselineind = indexarray[int(.9*ind):ind-1]
        polyx1 = np.arange(0, len(baselineind))
        polyy1 = np.arange(0, len(baselineind))
        polyx1 = basex[baselineind]
        polyy1 = lndata1[baselineind]
        # plt.plot(x1,y1)
        polynomial = np.polyfit(polyx1, polyy1, 1)
        baseline = np.polyval(polynomial, basex)
        # plt.plot(x,lndata,x,baseline)
        baseline_corrected = lndata-baseline
        # plt.plot(x,baseline_corrected)
        return baseline_corrected


    def exp(baseline_corrected):
        """takes the exp of the data"""
        expdata = np.exp(baseline_corrected)
        return expdata


    def flipdata(expdata):
        """flip data to get absolute ENDOR effect"""
        flipdata = np.subtract(expdata, 1)
        flipdata = flipdata*-1
        return flipdata


    def smooth(processed):
        """smoothes using a Savitsky-Golay filter. the default is to fit a 4th
        order polynomial over 6 points. this can be changed depending on how
        much you want to smooth. Increase the number of points to smooth more
        """
        smoothed = savgol_filter(processed, 7, 4)
        # For future this could be a window that you type the order and the
        # number of points into, and then it will plot it to show you the
        #smooth before moving on
        return smoothed


    def buildxy():
        """builds the x axis. Gets the number of points, the x minimum, and
        the x maximum from the dictionary created at the beginning of the
        script. Uses those values to create the x axis."""
        pad = 0
        xdim = list(map(float, get_from_dict('XPTS')))
        xdim = float(xdim[0])
        xdim_string = str(xdim)
        """should serve to get the largest dimension and pad the smaller data
        set with zeroes before subtracting them. Should be able to graph this
        new, subtrated data set (just one line, only works for 2 inputs)."""
        for i in range(0,len(filenames)):
            length_0 = len(xdim_string[0])
            length_i = len(xdim_string[i])
            if length_0 / length_i == 1:
                pad = 0
            else:
                if length_0 - length_i != 0:
                    pad = abs(length_0 - length_i)
                else:
                    0
                #run as for I in the range of the # of files, and then re-pull xmin?
        xdim_pad = np.pad(xdim, (pad, pad), 'constant')
        #return xdim_pad
        xmin = list(map(float, get_from_dict('XMIN')))
        xmin = float(xmin[0])
        xrange = list(map(float, get_from_dict('XWID')))
        xrange = float(xrange[0])
        xstep = xrange/(xdim_pad-1)
        freqx = (np.arange(xmin, xmin+xdim_pad*xstep, xstep))
        return freqx


    def calc_endorfreq():
        """calculates the ENDOR frequency for a proton. Uses the dictionary
        created in the beginning of the script to find the magnetic field to
        use in the calulation. Returns a frequency x-axis that is centered
        at the nucelar ENDOR frequency."""
        b0vl = float(str(get_from_dict('B0VL')))
        endorfreq = (300*b0vl)/7.046
        endorfreqx = (endorfreq-freqx) * (-1)
        return endorfreqx


    # def pad_axis():
    #     """should serve to get the largest dimension and pad the smaller data
    #     set with zeroes before subtracting them. Should be able to graph this
    #     new, subtrated data set (just one line, only works for 2 inputs)."""
    #     xdim = list(map(float, get_from_dict('XPTS')))
    #     xdim = float(xdim[0])
    #     xdim_pad = (if (xdim(0) / xdim(i) != 1:
    #                 pad = if (xdim(0) - xdim(i) > 0
    #                           abs(xdim(0) - xdim(i))
    #                           else = 0)
    #                 np.pad(xdim, (pad, pad) 'constant')))
    #     return xdim_pad


    
    file_dictionary = create_dict(filename)
    data = read_dta_file(filename)
    phaseddata = phase(data)
    lndata = naturallog(phaseddata)
    baseline_corrected = baseline_correct(lndata)
    expdata = exp(baseline_corrected)
    flipdata = flipdata(expdata)
    smoothed = smooth(flipdata)
    freqx = buildxy()
    endorfreqx = calc_endorfreq()
    processed = smoothed

    expx = np.arange(0, len(processed))

    plt.figure(1)
    plt.plot(endorfreqx, processed, linewidth=2)
    plt.figure(2)
    plt.plot(freqx, processed, linewidth=2)
    plt.show()

    return processed, endorfreqx, freqx

for i in range(0,len(filenames)):
    filename = filenamelist[i]
    processed[i], endorfreqx[i], freqx[i] = endor_process()


root.destroy()


#this works well; it doesn't subtract things yet though.
# this needs to be tested for actual use; can we run from github, or have to copy to python?
