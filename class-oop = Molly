# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 13:16:19 2018

@author: Molly Lockart
"""
# for data import
import os

try:
    import Tkinter as tk
except ImportError:
    import tkinter as tk

from tkinter import messagebox
from tkinter import filedialog
from tkinter import simpledialog
from tkinter import Label, Frame, OptionMenu, StringVar
from tkinter import ttk

# for data processing
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from scipy.optimize import minimize
from scipy.signal import savgol_filter


# Parent class; will add subclasses for processing Mims ENDOR, HYSCORE, etc.


def choose_files():

    """choose files from file dialog box; choose either the .DSC or the .DTA"""

    root = tk.Tk()
    root.withdraw()
#    messagebox.showinfo("For ENDOR Subtractions",
#                        "First choose file, then choose file to subtract
                            #from first file.")
    filenamelist = filedialog.askopenfilenames(parent=root, defaultextension='.DTA',
                                               title="Select Data files",
                                               filetypes=[('Bruker', '.DTA'),
                                                          ('all', '.*')])
    return filenamelist


def create_dict():

    """this creates a searchable dictionary out of the description file.
    it only takes experimental parameters; it ignores the pulse programs
    at the end of the description file."""


    if filename.endswith(".DSC"):
        dict_filename = filename

    else:

        try:
            dict_filename = filename.replace(".DTA", ".DSC")

        except ValueError:
            messagebox.showinfo('there is no valid .DSC file')

    dictionary = {}

    bruker_file = open(dict_filename, 'r')

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


def get_from_dict(key):

    """this gets something from the dictionary. the key input is
    what you want to pull from the dictionary. Ex: get_from_dict('XPTS')
    returns the number of points in the measurement."""

    value = file_dictionary[key]
    if (key != 'XPTS') and (key != 'XMIN') and (key != 'XWID'):
        value = " ".join(value) #takes string out of tuple so you can int() it
        value = value.strip("'")
    return value


class ProcessSpectra:

    """this is the parent class. It will choose a file, crease a file dictionar
    y,and phase the raw data"""

    def __init__(self, file_dictionary=0, rawdata=0, phased=0, phase_angle=0):

        """initializes the class Process_Spectra. This will load a file,
        create a dictionary, and read the .DTA file every time this class
        is invoked."""

        self.filename = filename

        def read_dta_file(self):

            """this reads in the data file. The Bruker .DTA file only contains
            the y values (if one-dimensional experiment) or the z values
            (if 2D experiment). The x axis is created later"""

            if exp == 'Mims ENDOR ESE':

                rawdata = np.fromfile(self.filename, dtype='>f8')

            else:
                if exp =='HYSCORE' or exp =='HYSCORE Split':
                    #big endian, f8 is 64-bit floating point number
                    rawbin = np.fromfile(self.filename, dtype='>f8')
                    rawbin = np.reshape(rawbin, (len(rawbin), 1,)) #transopse
                    value = get_from_dict('IKKF') #if complex data

                    if value == 'CPLX':
                        size = rawbin.size/2
                        size = int(size)
                        rawdata = np.zeros((size, 1), dtype=complex)
                    for i in range(0, len(rawdata)):
                        rawdata[i, 0] = np.complex(rawbin[((2*i)-1)+1],
                                                   rawbin[(2*i)+1])

                    if value != 'CPLX':
                        rawdata = rawdata
                    #gets number of points to reshape data into a square matrix
                    dim = int(np.sqrt(len(rawdata)))
                    rawdata = np.reshape(rawdata, (dim, dim))

            return rawdata

        def phase(self):

            """this phases the data by minimizing the imaginary component of the
            data."""

            if exp == "Mims ENDOR ESE":
                realdata = np.ones([int((len(self.rawdata)) / 2)])
                imagdata = np.ones([int((len(self.rawdata)) / 2)])

                for i in range(0, int((len(self.rawdata)) / 2)):
                    realdata[i] = self.rawdata[i * 2]
                    imagdata[i] = self.rawdata[i * 2 + 1]
                    complexdata = realdata + 1j * imagdata

                def f_phase(phi):

                    """phase function"""

                    return np.sum(abs((np.sin(phi) * complexdata.real -
                                       np.cos(phi) * complexdata.imag)))

                start_pos = np.array(0)
                res = minimize(f_phase, start_pos, method='Powell')
                phi = res.x
                complexdataraw = np.exp(-phi * 1j) * complexdata
                phased = complexdataraw.real

            else:

                if exp == 'HYSCORE' or exp =='HYSCORE Split':

                    def f_phase(x_a):

                        """this defines the phasing function"""

                        return np.sum(abs((np.sin(x_a)*
                                           self.rawdata.real-np.cos(x_a)*
                                           self.rawdata.imag)))

                start_pos = np.array(0)
                res = minimize(f_phase, start_pos, method='Powell')
                phi = res.x
                complexdataraw = np.exp(-phi * 1j) * self.rawdata
                phased = complexdataraw.real

            return phased


        self.rawdata = read_dta_file(self)
        self.phased = phase(self)
        self.phase_angle = phase(self)


    def lineplot(self):

        """makes a line plot from x and y data and saves it as a .tif figure"""

        plt.figure()
        plt.plot(spectrum.endorxdata, spectrum.data, 'k', linewidth=1.5)
        #xmin = float("".join(get_from_dict('XMIN')))
        #xwid = float("".join(get_from_dict('XWID')))
        #plt.xlim([xmin, xmin + xwid])
        plt.xlabel('Frequency (MHz)', fontsize=14)
        plt.ylabel('Absolute ENDOR Effect', fontsize=14)
        #plt.xticks(np.arange(xmin, xmin + xwid + 1, step=1))
        plt.tick_params(labelsize=12)
        save = [spectrum.filename, ".tiff"]
        plt.savefig("".join(save), dpi=100, bbox_inches='tight', format="tiff")

    def contour_plot(self):

        """makes a contour plot from x and y data and saves it as tif figure"""

        plt.figure()
        colors = ['w', 'royalblue', 'limegreen', 'gold', 'coral', 'firebrick']
        density_cmap = LinearSegmentedColormap.from_list('my_cmap', colors)

        try:
            vmin = np.ndarray.max(spectrum.zdata[1]) + 0.3
            plt.contourf(spectrum.xdata, spectrum.ydata, spectrum.zdata, 60,
                     cmap=density_cmap, vmin=vmin, vmax=0)

        except ValueError:
            vmin = np.ndarray.min(spectrum.zdata[1])
            plt.contourf(spectrum.xdata, spectrum.ydata, spectrum.zdata, 60,
                     cmap=density_cmap, vmin=vmin, vmax=0)

        #plt.colorbar()
        # more color contours = lower vmax
        # more noise = lower vmin
        #title = get_from_dict('TITL')
        #plt.title(title)
        plt.gca().set_aspect('equal')
        plt.ylabel('Frequency (MHz)', fontsize=12)
        plt.xlabel('Frequency (MHz)', fontsize=12)
        plt.tick_params(labelsize=12)
        save = [spectrum.filename, ".tiff"]
        plt.savefig("".join(save), dpi=100, bbox_inches='tight', format="tiff")



class MimsENDOR(ProcessSpectra):

    """This is a subclass that processes Mims ENDOR spectra."""

    def __init__(self):
        super().__init__(self)

        self.lndata = np.log1p(self.phased)

        def baseline_correct(self):

            """baseline corrects data by fitting to a zero-order polynomial.
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

        def flip(self):

            """flip data to get absolute ENDOR effect"""

            flipdata = np.subtract(self.expdata, 1)
            flipdata = flipdata * -1

            return flipdata

        def filter_data(self):

            """this attempt tries to rectify the standard deviation to be within the
            window, and not over the entire graph using the same +/- 3 std dev
            currently need to set up a fix--index error on line defining std_dev"""

            fixeddata = np.copy(self.flipdata)
            win_pct = 0.01
            window = int(win_pct * len(fixeddata))
            deviation = 5

            j = -1

            while  j < (len(fixeddata) - window - 1):

                j = j + 1

                # detect the glitch
                if j <  window - 1 or j > (len(fixeddata) - window):
                    fixeddata[j] = fixeddata[j]

                else:

                    arrays = np.arange((j - window),j), np.arange((j + 1) , (j + window + 1))
                    dev_range = np.append(arrays[0], arrays[1])
                    std_dev_compare = np.std(fixeddata[dev_range])
                    average = np.average(fixeddata[dev_range])

                    if abs(fixeddata[j] - average) > deviation * std_dev_compare:
                        fixeddata[j] = average
            return fixeddata

        def smooth(self):

            """smoothes using a Savitsky-Golay filter. the default is to fit
            a 4th order polynomial over an odd number of points. this can be
            changed depending on how much you want to smooth. Increase the
            number of points to smooth."""

            application_window = tk.Tk()
            application_window.withdraw()

            order_pts = simpledialog.askstring("Savitsky-Golay Filter Options",
                                               "Enter: polynomial order (0-6), number of smoothing points (e.g.: 4,15)",
                                parent=application_window)
            try:
                order = int(order_pts[0])
                pts = int(order_pts[2:])

            except TypeError:
                order = 4
                pts = 15

            smoothed = savgol_filter(self.deglitched, pts, order)

            return smoothed

        def buildxy(self):

            """builds the x axis. Gets the number of points, the x minimum, and
            the x maximum from the dictionary created at the beginning of the
            script. Uses those values to create the x axis."""

            xdim = int("".join(get_from_dict('XPTS')))
            # xdim_string = str(xdim)
            # pad = 0
            # xdim_pad = np.pad(xdim, (pad, pad), 'constant')
            # return xdim_pad
            xmin = float("".join(get_from_dict('XMIN')))
            xrange = float("".join(get_from_dict('XWID')))
            # xstep = xrange / (xdim_pad - 1)
            freqx_n = (np.linspace(xmin, xmin + xrange, num=xdim,
                                   endpoint=False))

            return freqx_n

        def calc_endorfreq(self):

            """calculates the ENDOR frequency for a proton. Uses the dictionary
            created in the beginning of the script to find the magnetic field to
            use in the calulation. Returns a frequency x-axis that is centered
            at the nucelar ENDOR frequency."""

            try:
                b0vl = float("".join(get_from_dict('B0VL')))
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
            endorfreqx = (endorfreq - self.xdata) * (-1)

            return endorfreqx


        self.baseline_corrected = baseline_correct(self)
        self.expdata = exp(self)
        self.flipdata = flip(self)
        self.deglitched = filter_data(self)
        self.data = smooth(self)
        self.xdata = buildxy(self)
        self.endorxdata = calc_endorfreq(self)



    def plot(self):

        """inherits the lineplot method from ProcessSpectra"""

        ProcessSpectra.lineplot(self)

class Hyscore(ProcessSpectra):

    """This is a subclass that processes Mims ENDOR spectra."""

    def __init__(self):
        super().__init__(self)

        def baseline_correct(self):
            """SAMPLE DOCSTRING this baseline corrects both the x and y axis"""
            x_a = np.arange(0, len(self.phased))
            # vector the same length as the data in one direction

            base = np.zeros((dim_matrix, dim_matrix))
            for i in range(0, len(self.phased)):
                phased_baseline = np.polyfit(x_a, self.phased[i, :], 3)
                # baseline correct in y direction
                baseline = np.polyval(phased_baseline, x_a)
                self.phased[i, :] = self.phased[i, :]-baseline
                base[i, :] = baseline

            for i in range(0, len(self.phased)):
                phased_baseline = np.polyfit(x_a, self.phased[:, i], 3)
                # baseline correct in x direction
                baseline = np.polyval(phased_baseline, x_a)
                self.phased[:, i] = self.phased[:, i]-baseline
                base[:, i] = base[:, i]+baseline
                self.baselined = self.phased
            return self.baselined


        def Butler_mize_window(self):
            """SAMPLE DOCSTRING this is windowing the data from an improved
            diagonal blackman window that can change the width of the window"""
            wind = np.ones([1, dim_matrix*2])
            a_0 = (1-alpha_0)/2
            a_1 = 1/2
            a_2 = a_1-a_0
            for i in range(0, dim_matrix*2):
                wind[0, i] = (a_0-a_1*np.cos(2*np.pi*(i)/(2*dim_matrix-1)) +
                              a_2*np.cos(4*np.pi*(i)/(2*dim_matrix-1)))
                # normal blackman function
            for i in range(1+x_m+2*x_0+dim_matrix, x_m+2*x_0+2*dim_matrix):
                # creates the right half
                wind[0, i] = wind[0, i]
            for i in range(0, x_m):  # define left side
                wind[0, i] = x_m+1
            wind[0, x_m+dim_matrix] = 1
            wind[0, x_m+dim_matrix-1] = 1
            # these two lines define the center;
            # they make positions 127, 128 both 1

            wind[0, x_m] = 0  # makes left side zero
            wind[0, 2*dim_matrix-1] = 0  # makes right side zero
            dwind = np.ones([dim_matrix, dim_matrix])
            # create the array for the next step
            for i in range(0, dim_matrix):
                dwind[:, i] = wind[0, (dim_matrix-i):2*dim_matrix-i]
            wind2 = np.ones([1, dim_matrix])
            for i_2 in range(dim_matrix-round(dim_matrix/4), dim_matrix):
                wind2[0, i_2] = abs(np.sin(np.pi/2*((i_2)/
                                                    (round(dim_matrix/4)))))
                # Taper
            wind2[0, dim_matrix-1] = 0
            wind3 = (wind2*(np.ones([dim_matrix, dim_matrix])))
            windowed = self.baselined*wind3*np.transpose(wind3)*dwind
            return windowed


        def fourier_transform2d(self):
            """ SAMPLE DOCSTRING 2D Fourier transform. Data is zerofilled
            (a matrix of all zeros is created and data is dropped in it) first.
            The transform is shifted so zero frequency is in the middle."""
            zerofill = np.zeros([512, 512])
            zerofill[:len(self.windowed), :len(self.windowed)] = self.windowed
            transform = np.fft.fft2(zerofill)
            transform = np.fft.fftshift(transform)
            transformed = np.absolute(transform)
            tmax = transformed.max()
            zdata = (transformed - 0)/(tmax - 0)
            # if you want a minimum at zero for your normalized plot,
            # then replace zero with the min.
            return zdata


        def buildxy(self):

            """ SAMPLE DOCSTRING this constructs x and y axis based
            off of description file"""

            x_dim = float("".join(get_from_dict('XPTS')))
            xmin = float("".join(get_from_dict('XMIN')))
            xrange = float("".join(get_from_dict('XWID')))

            d_x = xrange/(x_dim-1)
            x_axis = (np.arange(xmin, xmin+x_dim*d_x, d_x))

            # y_dim = float("".join(ProcessSpectra.get_from_dict('YPTS')))
            # ymin = list(map(float, get_from_dict('YMIN')))
            # yrange = float("".join(ProcessSpectra.get_from_dict('YWID')))

            frwidth = 1000/(x_axis[0])
            frinc = frwidth/(len(self.zdata))
            freq = np.arange(-frwidth, frwidth, frinc*2)
            xdata = freq
            ydata = freq
            return xdata, ydata

        dim_matrix = len(self.phased)
        x_m = 0
        x_0 = 0
        alpha_0 = 0.16

        self.baselined = baseline_correct(self)
        self.windowed = Butler_mize_window(self)
        # self.zdata = fourier_transform2d(self)
        self.zdata = np.log(fourier_transform2d(self))
        self.xdata = buildxy(self)[0]
        self.ydata = buildxy(self)[1]

    def contour_plot(self):

        """inherits the lineplot method from ProcessSpectra"""

        ProcessSpectra.contour_plot(self)


FILENAMELIST = choose_files()


for j in range(0, len(FILENAMELIST)):

    filename = FILENAMELIST[j]
    file_dictionary = create_dict()
    # spectrum  = np.arange(0,len(FILENAMELIST))
    exp = get_from_dict('PlsSPELEXPSlct')

    if exp == 'Mims ENDOR ESE' or exp == 'Mims ENDOR':
        spectrum = MimsENDOR()

        figure = MimsENDOR.plot(spectrum)
        spectrumname = os.path.split(filename)[1]
        #plt.title(spectrumname)
        np.savetxt(filename.replace(".DTA", ".txt"),
                   np.c_[spectrum.xdata, spectrum.data], delimiter=",")

    if  exp == 'HYSCORE' or exp =='HYSCORE Split':
        spectrum = Hyscore()
        figure = Hyscore.contour_plot(spectrum)
       # figure = Hyscore.plot()
        spectrum.data = [spectrum.xdata, spectrum.ydata, spectrum.zdata]

        spectrumname = os.path.split(filename)[1]
        #plt.title(spectrumname)
        np.savetxt(filename.replace(".DTA", ".txt"),
                   np.c_[spectrum.xdata, spectrum.ydata, spectrum.zdata],
                   delimiter=",")
# this saves a .txt where the first colum is the x data, the second column is
#y,and the next 512 columns are z
