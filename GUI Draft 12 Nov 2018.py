# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 12:40:28 2018

@author: Molly Lockart
"""
try:
    import Tkinter as tk
except ImportError:
    import tkinter as tk

from tkinter import messagebox
from tkinter import filedialog
from tkinter import simpledialog
import matplotlib.backends.tkagg as tkagg
from matplotlib.backends.backend_agg import FigureCanvasAgg
import numpy as np
import matplotlib.pyplot as plt
#from scipy import interpolate
from scipy.optimize import minimize
from scipy.signal import savgol_filter

LARGE_FONT = ("Verdana", 12)


# Choose Files
#root = tk.Tk()

def load_file():
    
    global data
    global file_dictionary
    global filename

    def choose_files():
        messagebox.showinfo("For ENDOR Subtractions","First choose file, then choose file to subtract from first file.")
        filenames = filedialog.askopenfilenames()
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
        filename = filenamelist[0]
        return filename
    
    
    
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
    
#    def get_from_dict(key):
#        """this gets something from the dictionary. the key input is
#        what you want to pull from the dictionary. Ex: get_from_dict('XPTS')
#        returns the number of points in the measurement."""
#        value = file_dictionary[key]
#        if (key != 'XPTS') and (key != 'XMIN') and (key != 'XWID'):
#            value = " ".join(value)
#            value = value.strip("'")
#        return value
    
    
    def read_dta_file():
        """this reads in the data file. The Bruker .DTA file only contains the
        y values (if one-dimensional experiment) or the z values (if it's a two
        dimensional experiment). The x axis is created later"""
        data = np.fromfile(filename, dtype='>f8')
        return data

    
    filename = choose_files()
    file_dictionary = create_dict(filename)
    data = read_dta_file()
     
    
    return filename, file_dictionary, data

def process():
    
    global processed 
    global freqx
    global endorfreqx
    
    def get_from_dict(key):
        """this gets something from the dictionary. the key input is
        what you want to pull from the dictionary. Ex: get_from_dict('XPTS')
        returns the number of points in the measurement."""
        value = file_dictionary[key]
        if (key != 'XPTS') and (key != 'XMIN') and (key != 'XWID'):
            value = " ".join(value)
            value = value.strip("'")
        return value
    
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
        order polynomial over an odd number of points. this can be changed depending on how
        much you want to smooth. Increase the number of points to smooth more
        """
        smoothed = savgol_filter(processed, 45, 6)
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
        # xdim_string = str(xdim)
        pad = 0
        xdim_pad = np.pad(xdim, (pad, pad), 'constant')
        #return xdim_pad
        xmin = list(map(float, get_from_dict('XMIN')))
        xmin = float(xmin[0])
        xrange = list(map(float, get_from_dict('XWID')))
        xrange = float(xrange[0])
        xstep = xrange/(xdim_pad-1)
        freqx_n = (np.arange(xmin, xmin+xdim_pad*xstep, xstep))
        return freqx_n


    def calc_endorfreq():
        """calculates the ENDOR frequency for a proton. Uses the dictionary
        created in the beginning of the script to find the magnetic field to
        use in the calulation. Returns a frequency x-axis that is centered
        at the nucelar ENDOR frequency."""
        try:
            b0vl = float(str(get_from_dict('B0VL')))
        except:
            root = tk.Tk()
            root.withdraw()
            b0vl = float(simpledialog.askstring("Bruker Error", "Enter Magnetic Field in Gauss",
                                parent=root))/10000 #enter field in Gauss and convert to Tesla
        endorfreq = (300*b0vl)/7.046
        endorfreqx = (endorfreq-freqx) * (-1)
        #endor_max = np.where(endorfreqx == np.max(endorfreqx))
        return endorfreqx 
    
    phaseddata = phase(data)
    lndata = naturallog(phaseddata)
    baseline_corrected = baseline_correct(lndata)
    expdata = exp(baseline_corrected)
    flipdata = flipdata(expdata)
    smoothed = smooth(flipdata)
    freqx = buildxy()
    endorfreqx = calc_endorfreq()
    processed = smoothed
    
    return processed 

class endorprocess(tk.Tk):
    def __init__(self, *args, **kwargs): #defines initialization function for the class
        tk.Tk.__init__(self, *args, **kwargs) #initializes tkinter
        window = tk.Frame(self) #makes your opening window
        window.pack(side = "top", fill = "both", expand = True) #fills the window 
        
        window.grid_columnconfigure(0, weight = 1) #sets a grid so yo can organize your window later
        window.grid_rowconfigure(0, weight = 1)
        
        self.frames = {} #empty dictionary for now--will populate it later
        
        frame = StartPage(window, self)
        
        self.frames[StartPage] = frame #not sure what this bit does
        frame.grid(row = 0, column = 0, sticky = "nsew") #nsew direction can be changed        
        
        self.show_frame(StartPage)
        
    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()
            
class StartPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent, background = 'white')
        label = tk.Label(self, text = "ENDOR Processing!!!", font = LARGE_FONT)
        label.grid(row = 0, column = 1, sticky = "n") #can use grid here somehow
           
        
        button1 = tk.Button(self, text = "load file!!!", bg = 'blue', command = load_file) 
        button1.grid(row = 1, column = 0)
        
        button2 = tk.Button(self, text = "process!!!", command = process) 
        button2.grid(row = 1, column = 2)
     
        try:
            f = plt.plot(endorfreqx, processed, linewidth=2)
            canvas = FigureCanvasTkAgg(f, master=self)
            canvas.show()
        except: 
            pass
       
app  = endorprocess()
app.mainloop() #tkinter code; necessary
    

