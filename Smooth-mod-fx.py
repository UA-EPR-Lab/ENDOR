# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 10:56:57 2018

@author: Joe Butler
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
from scipy import interpolate
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
filtered = [0, len(filenames)]
endorfreqx = [0, len(filenames)]
freqx = [0, len(filenames)]
x_spline = [0, len(filenames)]
y_spline = [0, len(filenames)]

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

    # def filter_data(flipdata):
    #     """should serve to filter outliers of +/- 3 std dev using a moving average of a set number of points around 'j' """
    #     std_dev = np.std(flipdata)
    #     out_low = - 3 * std_dev
    #     out_hi = 3 * std_dev
    #     win_pct = 0.01
    #     flip_pct = int(win_pct * len(flipdata))
    #     # percent of total data set observed by the window for moving average,
    #     # to be set as a GUI input variable later, = 1/2 the actual desired window size (1/2 on each side of center)
    #     j = -1
    #     while j < (len(flipdata) - 1):
    #         j = j+1
    #         if flipdata[j] < out_low or flipdata[j] > out_hi:
    #             for i in range ((-1 * flip_pct), flip_pct, 1):
    #                 flipdata[j] = flipdata[j+i] / (win_pct *len(flipdata))
    #         else: flipdata[j] = flipdata[j]
    #         #currently giving argument error for int(), 11/26/18
    #     filtered = flipdata
    #     return filtered

    def filter_data(flipdata):
        """this attempt tries to rectify the standard deviation to be within the
        window, and not over the entire graph using the same +/- 3 std dev
        currently need to set up a fix--index error on line defining std_dev"""
        
        win_pct = 0.01
        window = int(win_pct * len(flipdata))
        deviation = 2
        avg = 3
        
        j = -1
        while  j < (len(flipdata) - (avg +1)):
            j = j + 1
            
            try:
                for i in range ((-1 * window), window):
                    if j < (len(flipdata) - window - 1) and j > (len(flipdata) + window):
                        std_dev = np.std(range(flipdata[j - window], flipdata[j + window]))
                    else: std_dev = np.std(flipdata)

                out_low = -deviation * std_dev
                out_hi = deviation * std_dev
            except TypeError as e1:
                print(e1)
                
            try:
                if flipdata[j] < out_low or flipdata[j] > out_hi:
                    for k in range (-avg, avg + 1):
                        flipdata[j] = np.average(flipdata[j + k])
                else: flipdata[j] = flipdata[j]
            except TypeError as e2:
                print(e2)
                    
        filtered = flipdata
        return filtered
    

    def smooth(filtered):
        """smoothes using a Savitsky-Golay filter. the default is to fit a 4th
        order polynomial over an odd number of points. this can be changed depending on how
        much you want to smooth. Increase the number of points to smooth more
        """
        points = 45
        poly = 6
        
        smoothed = savgol_filter(filtered, points, poly)
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
        pad = 0
        xdim_pad = np.pad(xdim, (pad, pad), 'constant')
        #return xdim_pad
        xmin = list(map(float, get_from_dict('XMIN')))
        xmin = float(xmin[0])
        xrange = list(map(float, get_from_dict('XWID')))
        xrange = float(xrange[0])
        xstep = xrange/(xdim_pad-2)
        freqx_n = (np.arange(xmin, xmin+xdim_pad*xstep, xstep))
        return freqx_n


    def calc_endorfreq():
        """calculates the ENDOR frequency for a proton. Uses the dictionary
        created in the beginning of the script to find the magnetic field to
        use in the calulation. Returns a frequency x-axis that is centered
        at the nucelar ENDOR frequency."""
        b0vl = float(str(get_from_dict('B0VL')))
        endorfreq = (300*b0vl)/7.046
        endorfreqx = (endorfreq-freqx) * (-1)
        #endor_max = np.where(endorfreqx == np.max(endorfreqx))
        return endorfreqx
    

#    def pad_axis():
#         """should serve to get the largest dimension and pad the smaller data
#         set with zeroes before subtracting them. Should be able to graph this
#         new, subtrated data set (just one line, only works for 2 inputs)."""
#         xdim = list(map(float, get_from_dict('XPTS')))
#         xdim = float(xdim[0])
#         xstep = xrange/(xdim_pad-1)
#         endor_max0 = np.where(endorfreqx[0] == np.max(endorfreqx[0]))
#         endor_max1 = np.where(endorfreqx[1] == np.max(endorfreqx[1]))
#         #return xdim_pad
    
    
    file_dictionary = create_dict(filename)
    data = read_dta_file(filename)
    phaseddata = phase(data)
    lndata = naturallog(phaseddata)
    baseline_corrected = baseline_correct(lndata)
    expdata = exp(baseline_corrected)
    flipdata = flipdata(expdata)
    filtered = filter_data(flipdata)
    smoothed = smooth(filtered)
    freqx = buildxy()
    endorfreqx = calc_endorfreq()
    #processed = smoothed
    
    def spline_interpolation():
        """using a cubic spline interpolation to create a curve between each 
        set of 2 points. if no smoothing is desired, s = 0. usually
        s = m - sqrt(2m) where m = # datapoints (xpts). Order of the spline
        = k (cubic k = 3)."""
    
        x_pre_spline = freqx
        y_pre_spline = smoothed
        xdim = float(get_from_dict('XPTS')[0])
        # xdim = float(xdim[0])
        xmin = float(get_from_dict('XMIN')[0])
        # xmin = float(xmin)
        xrange = float(get_from_dict('XWID')[0])
        # xrange = float(xrange)
        xstep = xrange / xdim
        tck = interpolate.splrep(x_pre_spline, y_pre_spline)
        #s = ((xdim) - np.sqrt(2 * (xdim)))
        x_spline = np.arange(xmin, xmin +xrange, xstep)
        y_spline = interpolate.splev(x_spline, tck, der=0)
                
        return x_spline, y_spline
    
    
   # spline = spline_interpolation()
    #x_spline = spline[0]
    #y_spline = spline[1]
    
#    pad_def = pad_axis()

    expx = np.arange(0, len(smoothed))

    plt.figure(1)
    plt.plot(endorfreqx, smoothed, linewidth=1)
    plt.title('smoothed')
    plt.figure(2)
    plt.plot(endorfreqx, flipdata, linewidth=1)
    plt.title('flipdata')
    plt.figure(3)
    plt.plot(endorfreqx, filtered, linewidth=1)
    plt.title('filtered')
    
    plt.show() #will want to put ths entire plotting section after the 
    # different file sizes function

    return smoothed, endorfreqx, freqx, filtered
#, x_spline, y_spline

for i in range(0,len(filenames)):
    filename = filenamelist[i]
    smoothed[i], endorfreqx[i], freqx[i], filtered[i] = endor_process()
    # x_spline[i], y_spline[i]

# def different_file_sizes ():
#     for i in range(0,len(filenames)):
#         endor_max = np.where(processed == np.max(processed[i]))
#         return endor_max #and endor_max[i]
# endor_max = different_file_sizes()
# # do this out here, but use acubic spline to create ideal data then generate based on user input where you want to graph


#processed_subtracted = [(processed[0]) - (processed[1])]
#plt.figure(4)
#plt.plot(freqx[0], processed_subtracted[0], linewidth=1)
#plt.title('Subtraction from RF')

# def optimize_smooth():
#     """this should, ideally, compare the integral of the selection of points
#     [i,i+1] to the total integral of the range [0, end] and generate a ratio.
#     this ratio should then be able to be used to determine whether the new 
#     smoothing should be applied to the selection. Hopefully this will extra-smooth
#     the edges of the signal where there is little emphasis, and can make those
#     edge ranges the same so that when they are subtracted, there is a straight
#     line at 0 which should make for 'pretty' comparison spectra. We'll see how this goes"""
#     for i in range(0, len(processed)):  
#         xstep = float((max(freqx[0]) - min(freqx[0]))) / float(len(freqx[0]))
#         int_point = integrate.simps(y_spline[0], x_spline[0], dx = xstep, even = 'avg')
# optimize_smoothed = optimize_smooth()
root.destroy()



# needs to be tested with files of different sizes to really see if it works
# I had to comment something so that I could change the description
