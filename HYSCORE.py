# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 10:06:15 2019

@author: Joe Butler / Molly Lockart
"""
from Spec_Data import Spec_Data
import numpy as np


class Hyscore(Spec_Data):

    """This is a subclass that processes HYSCORE spectra."""

    def __init__(self, metadata, full_file_name):
        super().__init__( metadata, full_file_name)

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
                self.baseline_corrected = self.phased

            return self.baseline_corrected


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
            windowed = self.baseline_corrected*wind3*np.transpose(wind3)*dwind

            return windowed


        def fourier_transform2d(self):

            """ SAMPLE DOCSTRING 2D Fourier transform. Data is zerofilled
            (a matrix of all zeros is created and data is dropped in it) first.
            The transform is shifted so zero frequency is in the middle."""

            zerofill = np.zeros(1024 * np.array([1,1])) #so it will always be square
            zerofill[:len(self.windowed), :len(self.windowed)] = self.windowed
            transform = np.fft.fft2(zerofill)
            transform = np.fft.fftshift(transform) # shift center to zero
            transformed = np.absolute(transform)
            tmax = transformed.max()
            zdata = (transformed)/(tmax) # normalize to maximum value

            return zdata


        def buildxy(self):

            """ SAMPLE DOCSTRING this constructs x and y axis based
            off of description file"""

            x_dim = float(self.metadata['XPTS'])
            xmin = float(self.metadata['XMIN'])
            xrange = float(self.metadata['XWID'])

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

        self.baseline_corrected = baseline_correct(self)
        self.windowed = Butler_mize_window(self)
        # self.zdata = fourier_transform2d(self)
        self.zdata = np.log(fourier_transform2d(self))
        self.xdata = buildxy(self)[0]
        self.ydata = buildxy(self)[1]

