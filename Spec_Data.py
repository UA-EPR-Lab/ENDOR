# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 11:24:07 2019

@author: Molly Lockart / Joe Butler
"""
import numpy as np
from scipy.optimize import minimize
# import matplotlib.pyplot as plt

##############################################################################
##############################################################################

class Spec_Data:

    def __init__(self, metadata, full_file_name):

        # print(metadata['XPTS'])

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

        def buildaxis(self, axis): # future add a different # of points option
            """builds the x axis. Gets the number of points, the x minimum, and
            the x maximum from the dictionary created at the beginning of the
            script. Uses those values to create the x axis."""

            axdim = int(self.metadata[axis+'PTS'])
            axmin = float(self.metadata[axis+'MIN'])
            axrange = float(self.metadata[axis+'WID'])
            # xstep = xrange / (xdim_pad - 1)
            ax_n = (np.linspace(
                axmin, axmin + axrange, num=axdim, endpoint=False))
            ax_unit = self.metadata[axis+'UNI'][1:-1]

            return ax_n, ax_unit

        def phase(self):

            """this phases the data by minimizing the imaginary component of the
            data."""


            def f_phase(phi):

                """phase function"""

                return np.sum(abs((np.sin(phi) * self.amp.real -
                                   np.cos(phi) * self.amp.imag)))

            start_pos = np.array(0)
            res = minimize(f_phase, start_pos, method='Powell')
            phi = res.x
            phased_complex = np.exp(-phi * 1j) * self.amp
            phased=phased_complex.real

            return phased, phi, phased_complex


##############################################################################
# initialization varaibles

        self.metadata = metadata
        self.filename = full_file_name

        if self.metadata['IKKF'] == 'CPLX':
            self.complex = True
        else:
            self.complex = False

        if self.metadata['YTYP'] == 'NODATA':
            self.onedim = True
        else:
            self.onedim = False

        self.amp = read_bruker_file(self)

        # Decides what kind of data it is
        if self.metadata['IKKF'] == 'CPLX':
            self.complex = True
        else:
            self.complex = False

        if self.metadata['YTYP'] == 'NODATA':
            self.onedim = True
        else:
            self.onedim = False

        (self.x_axis, self.x_unit) = buildaxis(self, 'X')

        # Builds the y axis too if two dimsntional (HYSCORE)
        if self.onedim == False:
            (self.y_axis, self.yunit) = buildaxis(self, 'Y')

        (self.phased, self.phi, self.phased_complex) = phase(self)


##############################################################################
# Methods

    def return_axis(self, axis):
        if axis == 'x':
            return self.x_axis
        else:
            return self.y_axis

    def return_units(self, axis):
        if axis == 'x':
            return self.x_unit
        else:
            return self.y_unit

# methods for plotting here, interpolation, smoothing, phasing


#    def lineplot(self):
#        """makes a line plot from x and y data and saves it as a .tif figure"""
#
#        plt.figure()
#        plt.plot(self.x_axis, self.amp.real, 'k', linewidth=1.5)
#        plt.xlim([self.x_axis[0], self.x_axis[-1]])
#        plt.xlabel('RF Frequency ('+self.x_unit+')', fontsize=14)
#        plt.ylabel('Absolute ENDOR Effect', fontsize=14)
#        plt.tick_params(labelsize=12)
#        save = [self.filename.replace(".DTA", ".tiff")]
#        plt.savefig("".join(save), dpi=100, bbox_inches='tight', format="tiff")
## methods for plotting here, interpolation, smoothing, phasing, baseline correction with polynomial input for multiple corrections

# Also a method for returning an axis to feed plotting later.