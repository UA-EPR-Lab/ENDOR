# -*- coding: utf-8 -*-
# @author: Joe Butler
#child class of SpecData. Initializes SpecData to create the amplitudes and
# axes.

# this is where we will add the Mims ENDOR methods like ln, baseline, exp
# deglitch, and flip to the initialization.

# we can add a method for the calculation of ENDOR freq x.
from Spec_Data import Spec_Data
import numpy as np
import tkinter as tk

##############################################################################
##############################################################################

class MimsENDOR(Spec_Data):


    def __init__(self, metadata, full_file_name):
        super().__init__( metadata, full_file_name)

        #this knows self.phased from SpecData. It is the phased data that is ready
        # to be processed further.
        print(len(self.phased))
        # this is just to check that it knows what self.amp is.
        
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

            order_pts = np.simpledialog.askstring("Savitsky-Golay Filter Options",
                                               "Enter: polynomial order (0-6), number of smoothing points (e.g.: 4,15)",
                                parent=application_window)
            try:
                order = int(order_pts[0])
                pts = int(order_pts[2:])

            except TypeError:
                order = 4
                pts = 15

            smoothed = np.savgol_filter(self.deglitched, pts, order)

            return smoothed
        
    def calc_endorfreq(self):

            """calculates the ENDOR frequency for a proton. Uses the dictionary
            created in the beginning of the script to find the magnetic field to
            use in the calulation. Returns a frequency x-axis that is centered
            at the nucelar ENDOR frequency."""

            try:
                b0vl = float("".join(self.metadata('B0VL')))
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
                    endorfreqx = self.spectrum.xdata
                    print("Magnetic field should be a number")
                    return endorfreqx
            endorfreq = (300 * b0vl) / 7.046
            endorfreqx = (endorfreq - self.xdata) * (-1)

            return endorfreqx
        

