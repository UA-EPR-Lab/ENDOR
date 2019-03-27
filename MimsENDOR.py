# -*- coding: utf-8 -*-
# @author: Molly Lockart / Joe Butler
#child class of SpecData. Initializes SpecData to create the amplitudes and
# axes.

# this is where we will add the Mims ENDOR methods like ln, baseline, exp
# deglitch, and flip to the initialization.

# we can add a method for the calculation of ENDOR freq x.
from Spec_Data import Spec_Data
import numpy as np
import tkinter as tk
from tkinter import simpledialog
from scipy.signal import savgol_filter

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

        def deglitch(self):

            """this attempt tries to rectify the standard deviation to be within the
            window, and not over the entire graph using the same +/- 3 std dev
            currently need to set up a fix--index error on line defining std_dev"""

            deglitched = np.copy(self.flipdata)
            win_pct = 0.01
            window = int(win_pct * len(deglitched))
            deviation = 5

            j = -1

            while  j < (len(deglitched) - window - 1):

                j = j + 1

                # detect the glitch
                if j <  window - 1 or j > (len(deglitched) - window):
                    deglitched[j] = deglitched[j]

                else:

                    arrays = np.arange((j - window),j), np.arange((j + 1) , (j + window + 1))
                    dev_range = np.append(arrays[0], arrays[1])
                    std_dev_compare = np.std(deglitched[dev_range])
                    average = np.average(deglitched[dev_range])

                    if abs(deglitched[j] - average) > deviation * std_dev_compare:
                        deglitched[j] = average
            return deglitched

        def smooth(self):

            """smoothes using a Savitsky-Golay filter. the default is to fit
            a 4th order polynomial over an odd number of points. this can be
            changed depending on how much you want to smooth. Increase the
            number of points to smooth."""

            # errors and exceptions mean to keep user input within the bounds
            # of [0-6, 1-49 (odd)] and if ooutside, to return to the defualt.
            
            application_window = tk.Tk()
            application_window.withdraw()
            order_pts = simpledialog.askstring("Savitsky-Golay Filter Options",
                                               """Enter: Polynomial order (6 = max), and ODD number of smoothing points (49 = max)
            No space accepted in input.
            Default (4,15)""",
                                parent=application_window)
            if order_pts is None:
                smoothed = self.deglitched
            else:
                try:
                    order = int(order_pts[0]) 
                    if order > 6:
                        order = 6
                        application_window = tk.Tk()
                        application_window.withdraw()
                        simpledialog.messagebox.showerror("Polynomial Error",
                                                      """Invalid Sav-Gol Polynomial Order input detected. 
                Max (6) selected""",
                                                      parent = application_window)
                    pts = int(order_pts[1])
                    
                    if (pts/2) - np.round((pts / 2)) > 0:
                        pts = pts
                    else:
                        pts = pts -1
                        application_window = tk.Tk()
                        application_window.withdraw()
                        simpledialog.messagebox.showerror("Point Error",
                                                      """Invalid Sav-Gol Point input detected. 
                Rounded down to next odd number""",
                                                      parent = application_window)
                        
                    if pts > 49:
                        pts = 49
                        application_window = tk.Tk()
                        application_window.withdraw()
                        simpledialog.messagebox.showerror("Point Error",
                                                      """Invalid Sav-Gol Point input detected. 
                Max (49) selected""",
                                                      parent = application_window)
                    
                except ValueError:
                    order_pts = [4,15]
                    application_window = tk.Tk()
                    application_window.withdraw()
                    
                    simpledialog.messagebox.showerror("ValueError",
                                                      """Invalid Sav-Gol parameter input detected.
                No space accepted in input.
                Defualt (4,15) selected""",
                                                      parent = application_window)
                    
                except IndexError:
                    order_pts = [4,15]
                    application_window = tk.Tk()
                    application_window.withdraw()
                    
                    simpledialog.messagebox.showerror("IndexError",
                                                      """Invalid Sav-Gol parameter input detected.
                Defualt (4,15) selected""",
                                                      parent = application_window)            
                    
                order = order_pts[0]
                pts = order_pts[1]
                
                
                smoothed = savgol_filter(self.deglitched, pts, order)

            return smoothed
##############################################################################
##############################################################################
# initialization variables

        self.baseline_corrected = baseline_correct(self)
        self.expdata = exp(self)
        self.flipdata = flip(self)
        self.deglitched = deglitch(self)
        self.smoothed = smooth(self)

##############################################################################
##############################################################################
# class methods

    def calc_endorfreq(self):

            """calculates the ENDOR frequency for a proton. Uses the dictionary
            created in the beginning of the script to find the magnetic field to
            use in the calulation. Returns a frequency x-axis that is centered
            at the nucelar ENDOR frequency."""

            try:
                # b0vl = float("".join(self.metadata['B0VL'])
                b0vl = float(self.metadata['B0VL'])
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
                    endorfreqx = self.x_axis
                    print("Magnetic field should be a number")
                    return endorfreqx
            endorfreq = (300 * b0vl) / 7.046
            endorfreqx = (endorfreq - self.x_axis) * (-1)

            return endorfreqx


# test for uploading changes to GitHub