"""
NMAT: Nano-micromotor Analysis Tool

v0.7

05/01/2022

Changes from past version:

    -) Functions that compute particle-specific information, like MSD, MSAD, instvel,
        have been moved inside the Particle class.
    -) The script now computes the MSD directly from the x,y data. The time
        displacement (timeD) vector is now computed and not read from file, normalising
        the input time.
    -) Third and fourth degree fittings and full equation fitting added.
    -) R^2 added to the fittings.
    -) Added an option: "use only trackings longer than X s".

    
TO DO:
    
    -) Extract error from the fitting of the average MSD and add it to the results
    -) Add summary file of ALL the MSDs, trajectories, MSADs, etc., in the same
        file, inside the "individual" folder.
    

@author: Rafael Mestre; rmestre@ibecbarcelona.eu; rafmescas1@gmail.com
@contributions: certain functions come from the tracking software of Albert Miguel, v1.6.2

This code was written to be compatible for both Python 2.6+ and 3.3+, as well
as both Windows, Linux and Mac OS.

Additions to make the code compatible to Python 2 and 3
Check the page: https://python-future.org/compatible_idioms.html
for details
"""


from __future__ import print_function #To use print()
from __future__ import division #Division between integers is a float

import subprocess
def install(package):
    try:
        subprocess.check_call(["python", '-m', 'pip', 'install', package])
    except subprocess.CalledProcessError as e:
        print('Error installing package '+package+': ')
        raise

try:
    from future.utils import raise_with_traceback #To raise exceptions with traceback as:
    from future.utils import raise_from #For exception chaining
except:
    install('future')
    from future.utils import raise_with_traceback #To raise exceptions with traceback as:
    from future.utils import raise_from #For exception chaining    

from builtins import int #Log integers are just int. "builtins" module needs to be installed
from builtins import object
from builtins import range
import sys
version = sys.version_info[0]

import warnings

#ignoring warnings that result in nan's
warnings.filterwarnings("ignore", message="divide by zero encountered") 
warnings.filterwarnings("ignore", message="invalid value encountered")
warnings.filterwarnings("ignore", message="Covariance of the parameters")
warnings.filterwarnings("ignore", message="overflow encountered")
warnings.filterwarnings("ignore", message="Degrees of freedom <= 0 for slice")
warnings.filterwarnings("ignore", message="Mean of empty slice")


import tkinter
from tkinter import ttk
from tkinter import filedialog

import numpy as np
import os
import glob
import sys
from matplotlib import pyplot as plt
plt.switch_backend('agg') #doesn't show plots
import csv
from scipy.optimize import curve_fit

try:
    import seaborn as sns
except:
    install('seaborn')
    import seaborn as sns
    
try:
    import tidynamics
except:
    install('tidynamics')
    import tidynamics

try:
    from pathlib import Path
except:
    install('pathlib')
    from pathlib import Path

sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 3})
sns.set_style("ticks")
sns.set_palette(sns.color_palette("hls", 8))

def linear(x,a,b):
    return a + b*x

def quadratic(x,a,b):
    return a*x + b*x**2

def cubic(x,a,b,c):
    return a*x + b*x**2 + c*x**3

def fourthOrder(x,a,b,c,d):
    return a*x + b*x**2 + c*x**3 + d*x**4

def fullEquation(x,v,tau,D):
    return 4*D*x + 2*(v**2)*(tau**2)*(x/tau + np.exp(-x/tau) - 1)

def linearZero(x,a):
    return a*x
    
def control(x,a):
    return a*x

def exponential(x,a,b):
    return np.exp(a)*x**b

def powerLaw(x,a,b):
    return a*x**b

def powerLaw2(x,a,b):
    return x**a
    
def exponential2(x,a,b):
    return a*np.exp(-x/b)

def calculate_rsquared(xdata, ydata, popt, function):
    #From: https://stackoverflow.com/a/37899817
    residuals = ydata - function(np.asarray(xdata), *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared
        

plt.close("all")

def calculate_D(order):
    #Taken from tracking code v 1.6.2
    #Calculates the differentiation matrix D of arbitrary order n
    #For equidistant points x0, x1, ..., xn you can calculate the derivative of a set of points f(x) with the formula [f'(x0), f'(x1),...,f'(xn)] = 1/h * D * [f(x0), f(x1),...,f(xn)] 
    
    #In our case h = 1 (one frame of difference)
    order = order + 1
    nodes = np.arange(1,order+1)
    nodes_mask = np.ma.array(nodes, mask = False)

    #Baricentric weights
    weights = np.ndarray.astype(nodes, float)
    for i in range(order):
        nodes_mask.mask[i] = True
        weights[i] = 1.0/np.prod(nodes[i] - nodes_mask)
        nodes_mask.mask[i] = False

    #Calculation of Dij
    D = np.empty((order,order), dtype = float)
    for i in range(order):
        for j in range(order):
            if i != j:
                D[i,j] = weights[j]/(weights[i] * (nodes[i]-nodes[j]))
            else:
                nodes_mask.mask[i] = True
                D[i,j] = np.sum(1.0/ (nodes[i]-nodes_mask)[np.invert(nodes_mask.mask)])
                nodes_mask.mask[i] = False

    return D

'''For some reason I can't understand, this function needs to be outside...
If I put it inside centered_differences, the D calculated is not the same...
Need to figure out why'''
#Store D values
Ddict = {}
for order in range(2,10,2):
    Ddict[order] = calculate_D(order)
'''Up to here'''


def centered_difference(points, order):
    #Computes the first derivative of a set of points by applying the method of centered differences. Order can be 2, 4, 6 and 8
    #Taken from Tracking code v 1.6.2
    
    Dij = Ddict[order]
    coef = Dij[int(order/2), :][::-1]
    acc = int((len(coef)-1)/2)

    #Apply central difference
    deriv = np.convolve(points,coef,'valid')
    deriv = np.hstack((acc*[0],deriv,acc*[0]))

    #Fill initial and final velocities
    deriv[0] = np.array([points[1]-points[0]])
    deriv[-1] = np.array([points[-1]-points[-2]])
    if order > 2:
        deriv[1] = np.array([points[2]-points[0]])/2.0
        deriv[-2] = np.array([points[-1]-points[-3]])/2.0
        if order > 4:
            coef4CD = Ddict[4][2,:]
            deriv[2] = np.dot(coef4CD, points[:5])
            deriv[-3] = np.dot(coef4CD, points[-5:])
            if order > 6:
                coef6CD = Ddict[6][3,:]
                deriv[3] = np.dot(coef6CD, points[:7])
                deriv[-4] = np.dot(coef6CD, points[-7:])

    return deriv

def extend_angles(angles, period = 360):
    #Function to extend the angles to more than 360 degrees
    #Uses unwrap function from numpy
    #Taken from tracking software v 1.6.2
    #Modified so that it returns ALWAYS angles in rads
    
    if period == 360:
        return np.unwrap(np.deg2rad(angles))
    elif period == 2*np.pi:
        return np.unwrap(angles)    

def angular_auto_correlation(vector):
    l = len(vector)
    intList = np.arange(1,l)
    AAC = np.arange(1,l, dtype = float)
    AAC_dev = np.arange(1,l, dtype =float)
    for interval in intList:
        intVector = [1]+[0]*(interval-1)+[-1]
        cosDif = [np.cos(np.deg2rad(vector[i])) for i in range(len(vector))]
        convolution = np.convolve(cosDif, intVector, 'valid')
        AACList = convolution
        AAC[interval-1] = np.average(AACList)
        AAC_dev[interval-1] = np.std(AACList)
    return intList, AAC,AAC_dev


def velocity_auto_correlation(vx,vy):
    #Adapted from https://stackoverflow.com/questions/48844295/computing-autocorrelation-of-vectors-with-numpy
    
#    vx = np.asarray(vx) #Needs to be an array for the following to work
#    vy = np.asarray(vy) #Needs to be an array for the following to work
    n = len(vx)
    # correlate each component indipendently
    acorrx = np.correlate(vx,vx,'full')[n-1:]
    acorry = np.correlate(vy,vy,'full')[n-1:]
    # sum the correlations for each component
    acorr = np.sum([acorrx,acorry], axis = 0)
    # print(np.asarray(acorr).shape)
    # print(acorry)
    # divide by the number of values actually measured and return
    acorr = acorr/(n - np.arange(n))
    
    return acorr

class Particle:
    def __init__(self):
        
        self.particleLabel = None
        self.time = list()
        self.X = list()
        self.Y = list()
        self.MSD = list()
        self.time = list()
        self.v = list()
        self.vx = list()
        self.vy = list()
        self.w = list()
        self.angle = list()
        self.angleExtended = list()
        self.timeD = list()
        self.MSAD = list()
        self.autoCor = list()
        self.FPS = None
        self.fileName = ''
        
        self.alpha = list()
        self.alphaFitting = None
        self.diffusionAlphaFitting = None
        self.r_squaredAlphaFitting = None
        
        self.speedQuadraticFitting = None
        self.speedSquaredQuadraticFitting = None
        self.diffusionQuadraticFitting = None
        self.r_squaredQuadraticFitting = None
        
        self.diffusionLinearFitting = None
        self.r_squaredLinearFitting = None

        self.speedCubicFitting = None
        self.speedSquaredCubicFitting = None
        self.tauCubicFitting = None
        self.diffusionCubicFitting = None
        self.r_squaredCubicFitting = None

        self.speedFourthOrderFitting = None
        self.speedSquaredFourthOrderFitting = None
        self.tauFourthOrderFitting = None
        self.diffusionFourthOrderFitting = None
        self.r_squaredFourthOrderFitting = None

        self.speedFullFitting = None
        self.tauFullFitting = None
        self.diffusionFullFitting = None
        self.r_squaredFullFitting = None

        self.rotSpeedQuadraticFitting = None
        self.rotSpeedSquaredQuadraticFitting = None
        self.rotDiffusionQuadraticFitting = None
        self.rotDiffusionLinearFitting = None
        self.tauQuadraticFitting = None
        self.tauLinearFitting = None
        self.r_squaredRotLinearFitting = None
        self.r_squaredRotQuadraticFitting = None

        
        self.valid = False

    def calculateLocalAlpha(self):
        
        if self.valid:
            self.alpha = list()
            for i in range(len(self.MSD)-1):
                num = np.log(self.MSD[i+1])-np.log(self.MSD[i])
                den = np.log(self.timeD[i+1])-np.log(self.timeD[i])
                self.alpha.append(num/den)

    def calculateInstVel(self):
        
        if self.valid and len(self.v) == 0:
            #Calculates the instantenous velocity from the X,Y coordinates
            #From the centered differences of second order
            #Adapted from tracking software v 1.6.2
            order = 2
            X = self.X
            Y = self.Y
            
            if len(X) > 2*order:
                vx = centered_difference(X, order)
                vy = centered_difference(Y, order)
            
            vel_mod = np.sqrt(vx**2+vy**2) * self.FPS
            self.v = vel_mod
            self.vx = vx
            self.vy = vy


    def calculateAngles(self):

        if self.valid and len(self.angle) == 0:
            #Angle taken from velocity
            anglesVel = np.around(np.arctan2(self.vy, 
                                self.vx) * (180/np.pi)) 
            self.angle = anglesVel
            angles = list(extend_angles(anglesVel)) #Continuous angles
            self.angleExtended = angles
            #Calculates also rotational velocity
            w = centered_difference(angles, 2) * self.FPS
            self.w = w

    def calculateMSAD(self):
        
        if self.valid and len(self.MSAD) == 0:
            self.MSAD = self.mean_square(self.angleExtended)

    def mean_square(self,vector):
        '''Taken from tracking software v. 1.6.2'''
        #Input: vector with data
        #Output: mean square displacement given by MSD(p) = sum( (f(i+p)-f(i))**2)/total
        length = len(vector)
        intList = np.arange(1,length) #intervals
        MSD = np.arange(1,length, dtype = float) #To save the MSD values
        for interval in intList:
            intVector = [1]+[0]*(interval-1)+[-1] #Ex: for int = 3 you have [1,0,0,-1]
            #With "valid" you only take the overlapping points of the convolution
            convolution = np.convolve(vector,intVector,'valid')
            MSDlist = convolution**2
            MSD[interval-1] = np.average(MSDlist)
        return MSD

    def calculateAutoCorrelation(self):
        
        if self.valid and len(self.autoCor) == 0:
            vector = self.angleExtended
            l = len(vector)
            intList = np.arange(1,l)
            AAC = np.arange(1,l, dtype = float)
            AAC_dev = np.arange(1,l, dtype =float)
            for interval in intList:
                intVector = [1]+[0]*(interval-1)+[-1]
                convolution = np.convolve(vector, intVector, 'valid')
                AACList = np.cos(convolution)
                AAC[interval-1] = np.average(AACList)
                AAC_dev[interval-1] = np.std(AACList)
            self.autoCor = AAC


    def autocorrFFT(self,x):
        '''Calculates the autocorrelation FFT of a list of numbers.
        It's needed by the method MSD_fft'''
        
        N=len(x)
        F = np.fft.fft(x, n=2*N)  #2*N because of zero-padding
        PSD = F * F.conjugate()
        res = np.fft.ifft(PSD)
        res = (res[:N]).real   #now we have the autocorrelation in convention B
        n=np.arange(N, 0, -1) #divide res(m) by (N-m)
        return res/n #Normalized auto-correlation


    def MSD_fft(self,xvector,yvector,dt=1):
        r'''Performs the MSD very efficiently using FFT. The result is time averaged.    
        The discrete MSD is separated in S1 and 2*S2.
        Based on https://www.neutron-sciences.org/articles/sfn/abs/2011/01/sfn201112010/sfn201112010.html
        Code adapted from https://stackoverflow.com/questions/34222272/computing-mean-square-displacement-using-python-and-fft
        Scales with O(NlogN).
        '''
        
        pos = np.array([xvector, yvector]).T
        N=len(pos)
        
        time_list = np.arange(0,N)*dt
        
        
        D=np.square(pos).sum(axis=1) #x(i)**2 + y(i)**2
        D=np.append(D,0) #To make S1[0] equal to D[0]
        Q=2*D.sum()
        S1=np.zeros(N)
        
        for m in range(N):
            Q=Q-D[m-1]-D[N-m]
            S1[m]=Q/(N-m)
            
    #    print 's1', time.time() - s1t
    
        S2=sum([self.autocorrFFT(pos[:, i]) for i in range(pos.shape[1])])
        
    #    print 's2', time.time()-s2t
        
        msd = S1 - 2*S2
        return time_list, msd[0:]


    def calculateMSD(self):
        
        if self.valid and len(self.MSD) == 0:

            _, MSD = self.MSD_fft(xvector=self.X,
                                      yvector=self.Y,dt=1/self.FPS)
            self.MSD = MSD[1:]


class GUI:
    def __init__(self, master):
        
        self.master = master
        
        self.data = False
        
        self.dataOrdering = tkinter.StringVar(master,'Horizontal')
        self.dataType = tkinter.StringVar(master,'csv')
        self.delimeter = tkinter.StringVar(master,'\t')
        
        self.calculateMSD = tkinter.IntVar()
        self.calculateTrajectory = tkinter.IntVar()
        self.instVelocity = tkinter.IntVar()
        self.calculateMSAD = tkinter.IntVar()
        self.autocorrelation = tkinter.IntVar()
        
        self.autoCorFit = tkinter.IntVar()

        self.localExponent = tkinter.IntVar() #Calculate the local exponent
        self.MSDFitting = tkinter.IntVar() #Fit the MSD
        self.MSADFitting = tkinter.IntVar() #Firt the MSAD
        self.fitAlpha = tkinter.IntVar() #Perform logarithmic fitting
        
        self.fittingInterval = tkinter.StringVar(root, value="")
        self.useOnlyInterval = tkinter.StringVar(root, value="")
        
        self.nbPoints = -1
        
        self.calculateAverage = tkinter.IntVar() #Make averages over particles

        self.doAll = tkinter.IntVar()
        
        self.useOnly = tkinter.IntVar()
        
        self.particles = list()
        
        self.dn = Path(os.path.dirname(os.path.realpath(__file__)))
        

        self.folderPath = tkinter.StringVar(root, value=self.dn)
        self.pathEntry = tkinter.Entry(master,width=25,textvariable=self.folderPath)
        self.browseButton = tkinter.Button(master,text="Browse",command=self.browseFolder)

        self.dataOrderHorizontal = tkinter.Radiobutton(master, text='Horizontal',
                                                   value='Horizontal',variable=self.dataOrdering,
                                                   command=self.updateValuesAndEntries)
        self.dataOrderVertical = tkinter.Radiobutton(master, text='Vertical',
                                                   value='Vertical',variable=self.dataOrdering,
                                                   command=self.updateValuesAndEntries)


        self.dataTypeCsv = tkinter.Radiobutton(master, text='.csv',
                                                   value='csv',variable=self.dataType,
                                                   command=self.updateValuesAndEntries)
        self.dataTypeTxt = tkinter.Radiobutton(master, text='.txt',
                                                   value='txt',variable=self.dataType,
                                                   command=self.updateValuesAndEntries)

        self.delimeterTab = tkinter.Radiobutton(master, text=r'\t',
                                                   value='\t',variable=self.delimeter,
                                                   command=self.updateValuesAndEntries)
        self.delimeterComma = tkinter.Radiobutton(master, text=',',
                                                   value=',',variable=self.delimeter,
                                                   command=self.updateValuesAndEntries)

        # self.dataOrderHorizontal.select()

        self.firstText = tkinter.StringVar(master, value ="")
        
        self.text1 = tkinter.Label(master, text = self.firstText.get(),wraplength = 270, height=2)
        self.text2 = tkinter.Label(master,text = "What do you want to analyze?")
        self.text3 = tkinter.Label(master,text = "Data distribution:")
        self.text4 = tkinter.Label(master,text = "File type:")
        self.text5 = tkinter.Label(master,text = "Delimeter:")

        self.MSDButton = tkinter.Checkbutton(master, text="MSD", 
                                             variable=self.calculateMSD, 
                                             onvalue = True, offvalue = False)
        self.trajectoryButton = tkinter.Checkbutton(master, text="Trajectory", 
                                                    variable=self.calculateTrajectory, 
                                                    onvalue = True, offvalue = False)
        self.vButton = tkinter.Checkbutton(master, text="Instantaneous velocity", 
                                                    variable=self.instVelocity, 
                                                    onvalue = True, offvalue = False)
        self.autoCorrButton = tkinter.Checkbutton(master, text="Velocity autocorrelation", 
                                                    variable=self.autocorrelation, 
                                                    onvalue = True, offvalue = False,
                                                    command = self.autoCorChecked)
        self.MSADButton = tkinter.Checkbutton(master, text="MSAD", 
                                                    variable=self.calculateMSAD, 
                                                    onvalue = True, offvalue = False)

        self.autoCorFittingButton = tkinter.Checkbutton(master, text="Fit?", 
                                                    variable=self.autoCorFit, 
                                                    onvalue = True, offvalue = False)
        
        
        self.MSDOptions = tkinter.Label(master,text = "MSD options")
        
        
        self.exponentButton = tkinter.Checkbutton(master, text="Calculate local exponent of MSD?", 
                                                    variable=self.localExponent, 
                                                    onvalue = True, offvalue = False)       
        self.fittingButton = tkinter.Checkbutton(master, text="Fit MSD?", 
                                                    variable=self.MSDFitting, 
                                                    onvalue = True, offvalue = False,
                                                    command = self.fittingChecked)
        self.fittingCombo = ttk.Combobox(master, values = ["Linear", "Quadratic",
                                                           "Cubic", "Fourth order", 
                                                           "Full equation", "All"],width=14)
        self.fittingButtonMSAD = tkinter.Checkbutton(master, text="Fit MSAD?", 
                                                    variable=self.MSADFitting, 
                                                    onvalue = True, offvalue = False,
                                                    command = self.fittingMSADChecked)
        self.fittingComboMSAD = ttk.Combobox(master, values = ["Linear", "Quadratic","Both"],width=14)

        self.alphaButton = tkinter.Checkbutton(master,text="Perform logarithmic fitting to MSD?",
                                               variable=self.fitAlpha, 
                                               onvalue = True, offvalue = False)
        
        
        
        
        self.fitUpTo = tkinter.Label(master, text="Interval of interest:")
        self.fitIntervalEntry = tkinter.Entry(master,width=5,textvariable = self.fittingInterval)
        self.textSeconds = tkinter.StringVar(root, value="s out of ? s")
        self.seconds = tkinter.Label(master, text = self.textSeconds.get(),width=13)
        self.fittingInterval.trace_add("write", self.intervalsModified)
        
        self.useOnlyButton = tkinter.Checkbutton(master, text="Only trajectories longer than ", 
                                                    variable=self.useOnly, 
                                                    onvalue = True, offvalue = False, 
                                                    command = self.useOnlyChecked)
        self.useOnlyEntry = tkinter.Entry(master,width=5,textvariable = self.useOnlyInterval)
        self.textS = tkinter.Label(master, text = "s")
        self.useOnlyEntry.config(state="disabled")
        self.useOnlyInterval.trace_add("write", self.intervalsModified)

        
        self.AverageButton = tkinter.Checkbutton(master, text="Do averages between all particles?", 
                                                    variable=self.calculateAverage, 
                                                    onvalue = True, offvalue = False)

        
        self.analysisButton = tkinter.Button(master,text="Analyze",
                                             command=self.clickAnalysis)
        self.analysisButton.config(state="disabled")
        
        self.doAllButton = tkinter.Checkbutton(master,text="Do all",
                                               variable=self.doAll, 
                                               onvalue = True, offvalue = False,
                                               command = self.doEverything)
        self.doAllButton.config(state="disabled")
    
        '''Grid'''
        
        rowAnalyses = 6
        rowOptions = rowAnalyses + 6
        
        self.pathEntry.grid(row=0,column=0,padx=25,pady=20,columnspan=2,sticky = 'W')
        self.browseButton.grid(row=0,column=2,pady=20,columnspan=1)
        
        self.text1.grid(row=5,column=0,ipadx = 10, ipady = 5, columnspan=3, rowspan=1)

        self.text3.grid(row=1,column=0,ipadx = 10, ipady = 0, columnspan=1, rowspan=1, sticky='W')
        self.text4.grid(row=1,column=1,ipadx = 10, ipady = 0, columnspan=1, rowspan=1, sticky='W')
        self.text5.grid(row=1,column=2,ipadx = 10, ipady = 0, columnspan=1, rowspan=1, sticky='W')
        
        self.dataOrderHorizontal.grid(row=2,ipadx=10,ipady=0,column=0,columnspan=1, sticky='W')
        self.dataOrderVertical.grid(row=3,ipadx=10,ipady=0,column=0,columnspan=1, sticky='W')

        self.dataTypeCsv.grid(row=2,ipadx=10,ipady=0,column=1,columnspan=1, sticky='W')
        self.dataTypeTxt.grid(row=3,ipadx=10,ipady=0,column=1,columnspan=1, sticky='W')
        
        self.delimeterTab.grid(row=2,ipadx=10,ipady=0,column=2,columnspan=1, sticky='W')
        self.delimeterComma.grid(row=3,ipadx=10,ipady=0,column=2,columnspan=1, sticky='W')
        
        self.text2.grid(row=rowAnalyses,column=0,ipady = 5, columnspan=3)
        
        self.MSDButton.grid(row=rowAnalyses+1,column=0,sticky = 'W', ipadx = 10, columnspan=3)
        self.trajectoryButton.grid(row=rowAnalyses+2,column=0,sticky = 'W', ipadx = 10, columnspan=3)
        self.vButton.grid(row=rowAnalyses+3,column=0,sticky = 'W', ipadx = 10, columnspan=3)   
        self.MSADButton.grid(row=rowAnalyses+4,column=0,sticky = 'W', ipadx = 10, columnspan=3)   
        self.autoCorrButton.grid(row=rowAnalyses+5,column=0,sticky = 'W', ipadx = 10, columnspan=3)   
        
        self.MSDOptions.grid(row=rowOptions,column=0,ipady = 15, columnspan=3)
        
        self.exponentButton.grid(row=rowOptions+1,column=0,sticky = 'W',ipadx = 10, columnspan=3)
        self.fittingButton.grid(row=rowOptions+2,column=0,sticky = 'W', ipadx = 10, columnspan=3)
        self.fittingButtonMSAD.grid(row=rowOptions+3,column=0,sticky = 'W', ipadx = 10, columnspan=3)

        self.alphaButton.grid(row=rowOptions+4,column=0,sticky = 'W', ipadx = 10, columnspan=3)
    
        self.fitUpTo.grid(row=rowOptions+5,column=0,sticky = 'W',padx = 10,pady = 8)
        self.fitIntervalEntry.grid(row=rowOptions+5,column=1,sticky = 'W',pady = 8, columnspan=2)
        self.seconds.grid(row=rowOptions+5,column=1,sticky = 'W',pady = 8,padx = 50, columnspan=2)
        
        self.useOnlyButton.grid(row=rowOptions+6,column=0,sticky = 'W',pady = 8, padx = 10, columnspan=2)
        self.useOnlyEntry.grid(row=rowOptions+6,column=2,sticky = 'W',pady = 8, columnspan=1)
        self.textS.grid(row=rowOptions+6,column=2,sticky = 'W',pady = 8,padx = 50, columnspan=1)
        
        self.AverageButton.grid(row=rowOptions+7,column=0,sticky = 'W', ipadx = 10,ipady = 12, columnspan=3)
        self.analysisButton.grid(row=rowOptions+8,column=0, columnspan=3)

        self.doAllButton.grid(row=rowOptions+9,column=0,columnspan=3)

        #Initialize in the current folder
        self.updateValuesFolder()
        
        
#
    def browseFolder(self):
        self.dn = Path(filedialog.askdirectory(initialdir=self.dn))
                
        self.updateValuesFolder()
        
        if self.nbParticles > 0:
            self.fittingInterval.set("%.2f" % (self.longestTime/10))
            self.useOnlyInterval.set("%.2f" % (self.longestTime/10))
        else:
            self.fittingInterval.set("")
            self.useOnlyInterval.set("")


    def updateValuesAndEntries(self):
        
        self.updateValuesFolder()
        
        self.intervalsModified()

    def updateValuesFolder(self):
        
        self.readFiles()
        
        #Get valid particles
        self.validParticles= [p for p in self.particles if p.valid == True]
        self.nbParticles = len(self.validParticles)
        if self.nbParticles > 0:
            self.firstText.set("A total number of %d files were found with a total of %d valid particles" % (self.nbFiles,self.nbParticles))
            #Find particle with longest time
            totalTime = list()
            for part in range(self.nbParticles):
                totalTime.append(self.validParticles[part].timeD[-1])
            self.longestTime = max(totalTime)
            self.textSeconds.set("s out of %.2f s" % (self.longestTime))
            self.data = True
            self.analysisButton.config(state="normal")
            # self.fittingInterval.set("%.2f" % (self.longestTime/10))
            # self.useOnlyInterval.set("%.2f" % (self.longestTime/10))
            self.doAllButton.config(state="normal")
            
        else:
            self.firstText.set("No data was found in that folder")
            self.textSeconds.set("s out of ? s")
            self.data = False
            self.analysisButton.config(state="disabled")
            self.doAllButton.config(state="disabled")
            # self.fittingInterval.set("")
            # self.useOnlyInterval.set("")
        self.seconds.config(text=self.textSeconds.get())
        self.text1.config(text=self.firstText.get())
        self.folderPath.set(self.dn)
        self.pathEntry.config(text=self.folderPath)  
        self.fitIntervalEntry.config(text=self.fittingInterval)
        self.useOnlyEntry.config(text=self.useOnlyInterval)

    def doEverything(self):

        if self.doAll.get():
            self.MSDButton.select()
            self.trajectoryButton.select()
            self.vButton.select()  
            self.MSADButton.select() 
            self.autoCorrButton.select()  
                        
            self.exponentButton.select()
            self.fittingButton.select()
            self.fittingButtonMSAD.select()
            self.alphaButton.select()
                    
            self.AverageButton.select()
            
            self.fittingCombo.set("All")
            self.fittingComboMSAD.set("Both")
            self.fittingChecked()
            self.fittingMSADChecked()
        elif not self.doAll.get():
            self.MSDButton.deselect()
            self.trajectoryButton.deselect()
            self.vButton.deselect()  
            self.MSADButton.deselect() 
            self.autoCorrButton.deselect()  
                        
            self.exponentButton.deselect()
            self.fittingButton.deselect()
            self.fittingButtonMSAD.deselect()
            self.alphaButton.deselect()
                    
            self.AverageButton.deselect()
            
            self.fittingChecked()
            self.fittingMSADChecked()
                
        
    def fittingChecked(self):
        if self.MSDFitting.get():
            lv_x = self.fittingButton.winfo_x()
            lv_y = self.fittingButton.winfo_y()
            self.fittingCombo.place(x=lv_x+120,y=lv_y+2)
        else:
            self.fittingCombo.place_forget()

    def fittingMSADChecked(self):
        if self.MSADFitting.get():
            lv_x = self.fittingButtonMSAD.winfo_x()
            lv_y = self.fittingButtonMSAD.winfo_y()
            self.fittingComboMSAD.place(x=lv_x+120,y=lv_y+2)
        else:
            self.fittingComboMSAD.place_forget()


    def autoCorChecked(self):
        if self.autocorrelation.get():
            lv_x = self.autoCorrButton.winfo_x()
            lv_y = self.autoCorrButton.winfo_y()
            # self.autoCorFittingButton.place(x=lv_x+220,y=lv_y)
        else:
            self.autoCorFittingButton.place_forget()

    def useOnlyChecked(self):
        if self.useOnly.get():
            self.useOnlyEntry.config(state="normal")
            self.updateValuesAndEntries()
        else:
            self.useOnlyEntry.config(state="disabled")
            self.updateValuesAndEntries()

    def intervalsModified(self,var=None, index=None, mode=None):
         
        self.updateValuesFolder()
        
        if self.nbParticles > 0:
            disableAnalysis = False
        else:
            disableAnalysis = True
        
        try: 
            if self.useOnly.get():
                float(self.useOnlyInterval.get())
                self.useOnlyEntry.config(fg="black")
        except:
            self.useOnlyEntry.config(fg="red")
            disableAnalysis = True
            
        try:
            float(self.fittingInterval.get())
            if self.nbParticles > 0:
                if float(self.fittingInterval.get()) > self.longestTime:
                    self.fitIntervalEntry.config(fg="red")
                    disableAnalysis = True
                else:
                    self.fitIntervalEntry.config(fg="black")
            else:
                self.fitIntervalEntry.config(fg="black")
        except:
            self.fitIntervalEntry.config(fg="red")
            disableAnalysis = True
        
        
        if disableAnalysis:
            self.analysisButton.config(state="disabled")
        # else:
            # self.analysisButton.config(state="normal")


    def readFiles(self):
        
        self.particles = list()
        self.files = [i for i in self.dn.glob('*.{}'.format(self.dataType.get()))]
        self.nbFiles = len(self.files)
        
        for f in range(self.nbFiles):
                             
            fileName = self.files[f]
                        
            self.firstParticle = 0
            
            self.readData(fileName)
            
        for p in range(len(self.particles)):
            self.checkParticleValidity(self.particles[p])
            



    def readData(self,fileName):
        '''Reads the data from the tracking file'''


        with open(str(fileName),'r') as f:
            
            fileReader = csv.reader(f, delimiter=self.delimeter.get())
            
            if self.dataOrdering.get() == 'Horizontal':
                #These need to start with True, otherwise they might overwrite
                #information from the previous particle. They should be turned on 
                #when the particle is detected.
                timeRead = True
                XRead = True
                YRead = True

                for row in fileReader:
                                        
                    if len(row) == 0:
                        continue
                    
                    if 'particle' in row[0].lower():
                        #The algorithm reads data about a certain particle
                        #until the next instance of "particle" is found
                        try:
                            p = Particle()
                            p.particleLabel = int(row[1])
                            p.fileName = fileName
                            self.particles.append(p)
                            timeRead = False
                            XRead = False
                            YRead = False
                        except:
                            continue

                    if 'time' in row[0].lower() and not timeRead:
                        try:
                            self.particles[-1].time = [float(i) for i in row[1:] if i != '']
                            deltaT = self.particles[-1].time[1] - self.particles[-1].time[0]
                            self.particles[-1].timeD = [deltaT*i for i in range(1,len(self.particles[-1].time))]
                            self.particles[-1].FPS = 1/self.particles[-1].timeD[0]
                            timeRead = True
                        except:
                            continue
                    if 'x' in row[0].lower() and not XRead:
                        try:
                            self.particles[-1].X = [float(i) for i in row[1:] if i != '']
                            XRead = True
                        except:
                            continue
                    if 'y' in row[0].lower() and not YRead:
                        try:
                            self.particles[-1].Y = [float(i) for i in row[1:] if i != '']
                            YRead = True
                        except:
                            continue
                    
            if self.dataOrdering.get() == 'Vertical':
                
                particleFound = False
                timeFound = False
                xFound = False
                yFound = False
                
                time = list()
                x = list()
                y = list()
                
                timeColumn = None
                xColumn = None
                yColumn = None
                
                for row in fileReader:
                                        
                    if len(row) == 0:
                        continue

                    if 'particle' in row[0].lower() and not particleFound:
                        try:
                            p = Particle()
                            p.particleLabel = int(row[1])
                            p.fileName = fileName
                            self.particles.append(p)
                            particleFound = True
                        except:
                            continue
                    
                    
                        
                    if timeFound and xFound and yFound:
                        if not particleFound:
                            p = Particle()
                            p.particleLabel = 0
                            p.fileName = fileName
                            self.particles.append(p)
                            particleFound = True
                        try:
                            time.append(float(row[timeColumn]))
                            x.append(float(row[xColumn]))
                            y.append(float(row[yColumn]))
                        except:
                            break
                        
                    if not next((s for s in row if 'time' in s.lower()), None) is None and not timeFound:
                        timeColumn = row.index(next((s for s in row if 'time' in s.lower()), None))
                        timeFound = True
                    if not next((s for s in row if 'x' in s.lower()), None) is None and not xFound:
                        xColumn = row.index(next((s for s in row if 'x' in s.lower()), None))
                        xFound = True
                    if not next((s for s in row if 'y' in s.lower()), None) is None and not yFound:
                        yColumn = row.index(next((s for s in row if 'y' in s.lower()), None))
                        yFound = True
                    
                if timeFound and xFound and yFound:
                    
                    self.particles[-1].time = time
                    self.particles[-1].X = x
                    self.particles[-1].Y = y
                    if len(time)>1:
                        deltaT = time[1] - time[0]
                        self.particles[-1].timeD = [deltaT*i for i in range(1,len(self.particles[-1].time))]
                        self.particles[-1].FPS = 1/self.particles[-1].timeD[0]

    
            
    def checkParticleValidity(self,p):
        #A particle is only valid if it has X, Y, timeD and time information
        #Instantaneous velocity is note necessary
        
        if (not p.X) or (not p.Y) or  (not p.timeD) or (not p.time):
            p.valid = False
        else:
            if self.useOnly.get():
                try:
                    if p.timeD[-1] >= float(self.useOnlyInterval.get()):
                        p.valid = True
                    else:
                        p.valid = False
                except:
                    p.valid = False
            else:
                p.valid = True

    def clickAnalysis(self):
        
        print("Starting analysis")
        #Calculate the number of points for the fitting interval
        totalTime = list()
        for part in range(self.nbParticles):
            totalTime.append(self.validParticles[part].timeD[-1])
        maxParticle = np.argmax(totalTime)
        self.nbPoints = np.abs(np.asarray(self.validParticles[maxParticle].timeD) - float(self.fitIntervalEntry.get())).argmin()
        
        if self.nbPoints <= 1:
            self.nbPoints = -1

        # if self.useOnly.get():
        #     self.nbPointsMax = np.abs(np.asarray(self.validParticles[maxParticle].timeD) - float(self.useOnlyInterval.get())).argmin()
        # else:
        #     self.nbPointsMax = self.nbPoints
        

        if self.data:
            self.doDataAnalysis()

    def doDataAnalysis(self):
        
        if not os.path.exists(str(self.dn / Path('Plots'))):
            os.mkdir(str(self.dn/ Path('Plots')))
                
        #If something related to angles is selected
        #The angles are calculated, not taken from the tracking file
        #(because sometimes they don't select that option)
        if self.calculateMSAD.get() or self.autocorrelation.get() or self.MSADFitting.get():
            self.calculateAngles() #The angles are calculated, extended and expressed in rad

        #The MSAD and autocorrelation are calculated before computing the averages
        if self.MSADFitting.get() or self.calculateMSAD.get():
            self.computeMSAD()            

        if self.autocorrelation.get():
            self.calculateAutocorrelation()

        if self.MSDFitting.get() or self.calculateMSD.get():
            self.computeMSD()            
        
        if self.calculateAverage.get():
            self.computeAverages()
        

        if self.calculateMSD.get():
            self.doMSD()
            if self.calculateAverage.get():
                self.doAverageMSD()
        
        if self.calculateTrajectory.get():  
            self.doTrajectory()
            
        if self.instVelocity.get():
            #Instantenous velocity is not a requirement, it should be calculated
            #if it doesn't exist
            # self.calculateInstVel()
            self.doInstVel()
        
            
        if self.calculateMSAD.get():
            self.doMSAD()
            if self.calculateAverage.get():
                self.doAverageMSAD()
            
        if self.autocorrelation.get():
            self.doAutoCor()
            if self.calculateAverage.get():
                self.doAverageAutoCor()

        #Do MSD fittings
            
        if self.MSDFitting.get() and self.fittingCombo.get() == "Quadratic":
            self.doQuadraticFitting()
            if self.calculateAverage.get():
                self.doQuadraticFittingAverage()
            
        if self.MSDFitting.get() and self.fittingCombo.get() == "Linear":
            self.doLinearFitting()
            if self.calculateAverage.get():
                self.doLinearFittingAverage()

        if self.MSDFitting.get() and self.fittingCombo.get() == "Cubic":
            self.doCubicFitting()
            if self.calculateAverage.get():
                self.doCubicFittingAverage()

        if self.MSDFitting.get() and self.fittingCombo.get() == "Fourth order":
            self.doFourthOrderFitting()
            if self.calculateAverage.get():
                self.doFourthOrderFittingAverage()

        if self.MSDFitting.get() and self.fittingCombo.get() == "Full equation":
            self.doFullEquationFitting()
            if self.calculateAverage.get():
                self.doFullEquationFittingAverage()
                
        if self.MSDFitting.get() and self.fittingCombo.get() == "All":
            self.doLinearFitting()
            self.doQuadraticFitting()
            self.doCubicFitting()
            self.doFourthOrderFitting()
            self.doFullEquationFitting()
            if self.calculateAverage.get():
                self.doLinearFittingAverage()
                self.doQuadraticFittingAverage()
                self.doCubicFittingAverage()
                self.doFourthOrderFittingAverage()
                self.doFullEquationFittingAverage()
                
        #Do MSAD fittings
        
        if self.MSADFitting.get() and self.fittingComboMSAD.get() == "Quadratic":
            self.doQuadraticFittingMSAD()
            if self.calculateAverage.get():
                self.doQuadraticFittingMSADAverage()
            
        if self.MSADFitting.get() and self.fittingComboMSAD.get() == "Linear":
            self.doLinearFittingMSAD()
            if self.calculateAverage.get():
                self.doLinearFittingMSADAverage()
                
        if self.MSADFitting.get() and self.fittingComboMSAD.get() == "Both":
            self.doLinearFittingMSAD()
            self.doQuadraticFittingMSAD()
            if self.calculateAverage.get():
                self.doLinearFittingMSADAverage()
                self.doQuadraticFittingMSADAverage()

            
        if self.localExponent.get():
            self.doLocalAlpha()
            if self.calculateAverage.get():
                self.doAverageAlpha()
                
        if self.fitAlpha.get():
            self.doAlphaFitting()
            if self.calculateAverage.get():
                self.doAlphaFittingAverage()        

        print("Analysis finished")


                                                          
    def calculateAngles(self):
        #Calculates the angles from the instantaneous velocity
        #They are calculated in degrees
        #Adapted from tracking software v 1.6.2
        
        #Checks if information about instantaneous velocity is there

        
        for part in range(self.nbParticles):
            #Calculate instantaneous velocity if it hasn't been done
            self.validParticles[part].calculateInstVel()
            self.validParticles[part].calculateAngles()


#     def calculateAutocorrelation(self):
#         #Calculates the velocity autocorrelation function
#         for part in range(self.nbParticles):
#             vx = self.validParticles[part].vx
#             vy = self.validParticles[part].vy
#             vel = [[vx[i],vy[i]] for i in range(len(vx))]
# #            vac = velocity_auto_correlation(vx,vy)
#             vac = tidynamics.acf(vel) * self.validParticles[part].v[0]**2
# #            vac = vac * self.validParticles[part].v**2
# #            _, AAC, AAC_dev = angular_auto_correlation(self.validParticles[part].angleExtended)
# #            print(vac[:10])
# #            print(AAC[:10])
#             self.validParticles[part].autoCor = vac[:-1]

    def calculateAutocorrelation(self):
        
        for part in range(self.nbParticles):
            self.validParticles[part].calculateAutoCorrelation()

    
    def computeAverages(self):
        
        #Take all MSDs and the total length of video
        listMSD = list()
        totalTime = list()
        listMSAD = list()
        listAutoCor = list()
        totalTimeMSAD = list()
        totalTimeAutoCor = list()
        for part in range(self.nbParticles):
            
            totalTime.append(self.validParticles[part].timeD[-1])
            listMSD.append(self.validParticles[part].MSD)
            if len(self.validParticles[part].MSAD) > 0:
                listMSAD.append(self.validParticles[part].MSAD)
                totalTimeMSAD.append(self.validParticles[part].timeD[-1])
            if len(self.validParticles[part].autoCor) > 0:
                listAutoCor.append(self.validParticles[part].autoCor)
                totalTimeAutoCor.append(self.validParticles[part].timeD[-1])
                
        #Average all MSDs
        lengthsMSD = [len(listMSD[i]) for i in range(len(listMSD))]
        longestMSD = np.max(lengthsMSD)
        #Need to find a way to make this below better
        self.nMSD = [len([listMSD[i][j] for i in range(len(listMSD)) if j < len(listMSD[i])]) for j in range(longestMSD)]
        self.averageMSD = [np.mean([listMSD[i][j] for i in range(len(listMSD)) if j < len(listMSD[i])]) for j in range(longestMSD)]                    
        self.stdMSD = [np.std([listMSD[i][j] for i in range(len(listMSD)) if j < len(listMSD[i])],ddof=1) for j in range(longestMSD)]
        self.semMSD = [self.stdMSD[i]/np.sqrt(self.nMSD[i]) for i in range(len(self.nMSD))]
        
        #Average all MSADs
        if len(listMSAD) > 0:
            lengthsMSAD = [len(listMSAD[i]) for i in range(len(listMSAD))]
            longestMSAD = np.max(lengthsMSAD)
            self.averageMSAD = [np.mean([listMSAD[i][j] for i in range(len(listMSAD)) if j < len(listMSAD[i])]) for j in range(longestMSAD)]                    
            longestParticleMSAD = totalTime.index(max(totalTimeMSAD))
            self.averageTimeDMSAD = self.validParticles[longestParticleMSAD].timeD   
            self.nMSAD = [len([listMSAD[i][j] for i in range(len(listMSAD)) if j < len(listMSAD[i])]) for j in range(longestMSAD)]
            self.stdMSAD = [np.std([listMSAD[i][j] for i in range(len(listMSAD)) if j < len(listMSAD[i])],ddof=1) for j in range(longestMSAD)]
            self.semMSAD = [self.stdMSAD[i]/np.sqrt(self.nMSAD[i]) for i in range(len(self.nMSAD))]
        
        #Average all autocorrelations
        if len(listAutoCor) > 0:
            lengthsAutoCor = [len(listAutoCor[i]) for i in range(len(listAutoCor))]
            longestAutoCor = np.max(lengthsAutoCor)
            self.averageAutoCor = [np.mean([listAutoCor[i][j] for i in range(len(listAutoCor)) if j < len(listAutoCor[i])]) for j in range(longestAutoCor)]                    
            longestParticleAutoCor = totalTime.index(max(totalTimeAutoCor))
            self.averageTimeDAutoCor = self.validParticles[longestParticleAutoCor].timeD        
        
        #Find particle with the longest time
        longestParticle = totalTime.index(max(totalTime))
        self.averageTimeD = self.validParticles[longestParticle].timeD
        
    

    def doLinearFitting(self):
        
        #Creates the folder fittings if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD')))
            
        #Creates the folder for alpha fitting if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD/Linear'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD/Linear')))
        
        directorySave = self.dn / Path('Plots/FittingsMSD/Linear') 
        
        #Plot MSD with fitting info
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            popt, pcov = curve_fit(linearZero, timeD[:self.nbPoints], MSD[:self.nbPoints])  
            
            self.validParticles[part].diffusionLinearFitting = popt[0]/4
            
            self.validParticles[part].r_squaredLinearFitting = calculate_rsquared(timeD[:self.nbPoints],MSD[:self.nbPoints],popt,linearZero)
            
            plt.plot(timeD[:self.nbPoints],MSD[:self.nbPoints])
            
            label = "Fitting equation: MSD = $4D_{t}t$\n$D_t$ = %.4f $\mu$m$^2$/s\n$R^2$ = %.4f" % (popt[0]/4,self.validParticles[part].r_squaredLinearFitting)
            plt.plot(timeD[:self.nbPoints],linearZero(np.asarray(timeD[:self.nbPoints]),*popt),'--',
                     label = label)
        
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSD ($\mu$m$^2$)')
            plt.axis('tight')
            plt.legend()
            plt.title('MSD linear fitting for particle ' + str(part))
            
            fig0.savefig(str(directorySave/('MSD_linearFitting_P'+str(part)+'.png')))
            fig0.savefig(str(directorySave/('MSD_linearFitting_P'+str(part)+'.svg')),format='svg',dpi=1200)
            plt.close()

        #Save data
        for p in range(self.nbParticles):
            with open(str(directorySave/('MSD_linearFitting_P'+str(p)+'.csv')), 'w') as textfile:
                textfile.write('File,'+str(self.validParticles[p].fileName)+'\n')
                textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                textfile.write("Fitting equation,MSD = 4D*t\n")
                textfile.write("D (um^2/s),%.6f\n" % (self.validParticles[p].diffusionLinearFitting))
                textfile.write("\nR^2,%.6f\n" % (self.validParticles[p].r_squaredLinearFitting))
        
        #Save summary data
        with open(str(directorySave/'MSD_linearFitting_Summary.csv'), 'w') as textfile:
            Dlist = list()
            textfile.write("Fitting equation,MSD = 4D*t\n")
            textfile.write('D (um^2/s),')
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].diffusionLinearFitting))
                Dlist.append(self.validParticles[p].diffusionLinearFitting)
            textfile.write('\n')
            textfile.write('\n')
            textfile.write('Average diffusion,%.6f\n' % (np.mean(Dlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(Dlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(Dlist,ddof=1)/np.sqrt(len(Dlist))))

    def doLinearFittingAverage(self):
        #The doAverageMSD function should have been run first
        #so that the averageMSD variable exists
        
        #Creates the folder if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD')))
        
        directorySave = self.dn / Path('Plots/FittingsMSD')          

        #Plot MSD  fitting info
        fig0 = plt.figure(0,figsize=(12, 10))
                
        popt, pcov = curve_fit(linearZero, self.averageTimeD[:self.nbPoints], self.averageMSD[:self.nbPoints])  
        
        r_squared = calculate_rsquared(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints],popt,linearZero)

        plt.plot(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints])
        
        label = "Fitting equation: MSD = $4D_{t}t$\n$D_t$ = %.4f $\mu$m$^2$/s\n$R^2$ = %.4f" % (popt[0]/4, r_squared)
        plt.plot(self.averageTimeD[:self.nbPoints],linearZero(np.asarray(self.averageTimeD[:self.nbPoints]),*popt),'--',
                 label = label)
    
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$)')
        plt.axis('tight')
        plt.legend()
        plt.title('MSD linear fitting for all particles')
        
        fig0.savefig(str(directorySave/'MSD_linearFitting_average.png'))
        fig0.savefig(str(directorySave/'MSD_linearFitting_average.svg'),format='svg',dpi=1200)
        plt.close()
    
        #Save data
        with open(str(directorySave/'MSD_linearFitting_average.csv'), 'w') as textfile:
            textfile.write("Fitting equation,MSD = 4D*t\n")
            textfile.write("D (um^2/s),%.6f\n" % (popt[0]/4))
            textfile.write("\nR^2,%.6f\n" % (r_squared))


    def doLinearFittingMSAD(self):
        
        #Creates the folder fittings if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSAD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSAD')))
            
        #Creates the folder for alpha fitting if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSAD/Linear'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSAD/Linear')))
        
        directorySave = self.dn / Path('Plots/FittingsMSAD/Linear')
        
        #Plot MSD with fitting info
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSAD = self.validParticles[part].MSAD
            
            popt, pcov = curve_fit(linearZero, timeD[:self.nbPoints], MSAD[:self.nbPoints])  
            
            self.validParticles[part].rotDiffusionLinearFitting = popt[0]/4
            self.validParticles[part].tauLinearFitting = 1/self.validParticles[part].rotDiffusionLinearFitting
            
            self.validParticles[part].r_squaredRotLinearFitting = calculate_rsquared(timeD[:self.nbPoints],MSAD[:self.nbPoints],popt,linearZero)

            plt.plot(timeD[:self.nbPoints],MSAD[:self.nbPoints])
            
            label = "Fitting equation: MSD = $4D_{r}t$\n$D_t$ = %.4f rad$^2$/s\n$\\tau_r$ = %.4f s\n$R^2$ = %.4f" % (popt[0]/4,4/popt[0],self.validParticles[part].r_squaredRotLinearFitting)
            plt.plot(timeD[:self.nbPoints],linearZero(np.asarray(timeD[:self.nbPoints]),*popt),'--',
                     label = label)
        
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSAD (rad$^2$)')
            plt.axis('tight')
            plt.legend()
            plt.title('MSAD linear fitting for particle ' + str(part))
            
            fig0.savefig(str(directorySave/('MSAD_linearFitting_P'+str(part)+'.png')))
            fig0.savefig(str(directorySave/('MSAD_linearFitting_P'+str(part)+'.svg')),format='svg',dpi=1200)
            plt.close()

        #Save data
        for p in range(self.nbParticles):
            with open(str(directorySave/('MSAD_linearFitting_P'+str(p)+'.csv')), 'w') as textfile:
                textfile.write('File,'+str(self.validParticles[p].fileName)+'\n')
                textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                textfile.write("Fitting equation,MSD = 4Dr*t\n")
                textfile.write("Dr (rad^2/s),%.6f\n" % (self.validParticles[p].rotDiffusionLinearFitting))
                textfile.write("tau (s),%.6f\n" % (self.validParticles[p].tauLinearFitting))
                textfile.write("\nR^2,%.6f\n" % (self.validParticles[p].r_squaredRotLinearFitting))

        #Save summary data
        with open(str(directorySave/('MSAD_linearFitting_Summary.csv')), 'w') as textfile:
            Dlist = list()
            taulist = list()
            textfile.write("Fitting equation,MSAD = 4Dr*t\n")
            textfile.write('Dr (rad^2/s),')
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].rotDiffusionLinearFitting))
                Dlist.append(self.validParticles[p].rotDiffusionLinearFitting)
            textfile.write('\n')
            textfile.write("tau (s),")
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].tauLinearFitting))
                taulist.append(self.validParticles[p].tauLinearFitting)
            textfile.write('\n')
            textfile.write('\n')
            textfile.write('Average diffusion,%.6f\n' % (np.mean(Dlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(Dlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(Dlist,ddof=1)/np.sqrt(len(Dlist))))
            textfile.write('\n')
            textfile.write('Average tau,%.6f\n' % (np.mean(taulist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(taulist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(taulist,ddof=1)/np.sqrt(len(taulist))))

    def doLinearFittingMSADAverage(self):
        #The doAverageMSAD function should have been run first
        #so that the averageMSAD variable exists
        
        #Creates the folder if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSAD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSAD')))
        
        directorySave = self.dn / Path('Plots/FittingsMSAD')        

        #Plot MSD  fitting info
        fig0 = plt.figure(0,figsize=(12, 10))
                
        popt, pcov = curve_fit(linearZero, self.averageTimeD[:self.nbPoints], self.averageMSAD[:self.nbPoints])  
                
        r_squared = calculate_rsquared(self.averageTimeD[:self.nbPoints],self.averageMSAD[:self.nbPoints],popt,linearZero)

        plt.plot(self.averageTimeD[:self.nbPoints],self.averageMSAD[:self.nbPoints])
        
        label = "Fitting equation: MSAD = $4D_{r}t$\n$D_r$ = %.4f rad$^2$/s\n$\\tau_r$ = %.4f s\n$R^2$ = %.4f" % (popt[0]/4,4/popt[0],r_squared)
        plt.plot(self.averageTimeD[:self.nbPoints],linearZero(np.asarray(self.averageTimeD[:self.nbPoints]),*popt),'--',
                 label = label)
    
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSAD (rad$^2$/s)')
        plt.axis('tight')
        plt.legend()
        plt.title('MSAD linear fitting for all particles')
        
        fig0.savefig(str(directorySave/'MSAD_linearFitting_average.png'))
        fig0.savefig(str(directorySave/'MSAD_linearFitting_average.svg'),format='svg',dpi=1200)
        plt.close()
    
        #Save data
        with open(str(directorySave/'MSAD_linearFitting_average.csv'), 'w') as textfile:
            textfile.write("Fitting equation,MSAD = 4Dr*t\n")
            textfile.write("Dr (rad^2/s),%.6f\n" % (popt[0]/4))
            textfile.write("tau (s),%.6f\n" % (4/popt[0]))
            textfile.write("\nR^2,%.6f\n" % (r_squared))




    def doQuadraticFitting(self):
        
        #Creates the folder fittings if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD')))
            
        #Creates the folder for alpha fitting if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD/Quadratic'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD/Quadratic')))
        
        directorySave = self.dn / Path('Plots/FittingsMSD/Quadratic')
        
        #Plot MSD with fitting info
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            popt, pcov = curve_fit(quadratic, timeD[:self.nbPoints], MSD[:self.nbPoints])  
            
            self.validParticles[part].speedQuadraticFitting = np.sqrt(popt[1])
            self.validParticles[part].speedSquaredQuadraticFitting = popt[1]
            self.validParticles[part].diffusionQuadraticFitting = popt[0]/4
            
            self.validParticles[part].r_squaredQuadraticFitting = calculate_rsquared(timeD[:self.nbPoints],MSD[:self.nbPoints],popt,quadratic)

            plt.plot(timeD[:self.nbPoints],MSD[:self.nbPoints])
            
            label = "Fitting equation: MSD = $4D_{t}t + v^2t^2$\n$D_t$ = %.4f $\mu$m$^2$/s\n$v$ = %.4f $\mu$m/s\n$v^2$ = %.4f $\mu$m$^2/s^2$\n$R^2$ = %.4f" % (popt[0]/4,np.sqrt(popt[1]),popt[1],self.validParticles[part].r_squaredQuadraticFitting)
            plt.plot(timeD[:self.nbPoints],quadratic(np.asarray(timeD[:self.nbPoints]),*popt),'--',
                     label = label)
        
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSD ($\mu$m$^2$)')
            plt.axis('tight')
            plt.legend()
            plt.title('MSD quadratic fitting for particle ' + str(part))
            
            fig0.savefig(str(directorySave/('MSD_quadraticFitting_P'+str(part)+'.png')))
            fig0.savefig(str(directorySave/('MSD_quadraticFitting_P'+str(part)+'.svg')),format='svg',dpi=1200)
            plt.close()

        #Save data
        for p in range(self.nbParticles):
            with open(str(directorySave/('MSD_quadraticFitting_P'+str(p)+'.csv')), 'w') as textfile:
                textfile.write('File,'+str(self.validParticles[p].fileName)+'\n')
                textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                textfile.write("Fitting equation,MSD = 4D*t + v^2*t^2\n")
                textfile.write("D (um^2/s),%.6f\n" % (self.validParticles[p].diffusionQuadraticFitting))
                textfile.write("v (um/s),%.6f\n" % (self.validParticles[p].speedQuadraticFitting))
                textfile.write("v^2 (um^2/s^2),%.6f\n" % (self.validParticles[p].speedSquaredQuadraticFitting))
                textfile.write("\nR^2,%.6f\n" % (self.validParticles[p].r_squaredQuadraticFitting))

        #Save summary data
        with open(str(directorySave/('MSD_quadraticFitting_Summary.csv')), 'w') as textfile:
            Dlist = list()
            vlist = list()
            vsquarelist = list()
            textfile.write("Fitting equation,MSD = 4D*t + v^2*t^2\n")
            textfile.write('D (um^2/s),')
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].diffusionQuadraticFitting))
                Dlist.append(self.validParticles[p].diffusionQuadraticFitting)
            textfile.write('\n')
            textfile.write("v (um/s),")
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].speedQuadraticFitting))
                vlist.append(self.validParticles[p].speedQuadraticFitting)
            textfile.write('\n')
            textfile.write("v^2 (um^2/s^2),")
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].speedSquaredQuadraticFitting))
                vsquarelist.append(self.validParticles[p].speedSquaredQuadraticFitting)
            textfile.write('\n')
            textfile.write('\n')
            textfile.write('Average diffusion,%.6f\n' % (np.mean(Dlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(Dlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(Dlist,ddof=1)/np.sqrt(len(Dlist))))
            textfile.write('\n')
            textfile.write('Average speed,%.6f\n' % (np.mean(vlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(vlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(vlist,ddof=1)/np.sqrt(len(vlist))))
            textfile.write('\n')
            textfile.write('Average speed squared,%.6f\n' % (np.mean(vsquarelist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(vsquarelist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(vsquarelist,ddof=1)/np.sqrt(len(vsquarelist))))

    def doQuadraticFittingAverage(self):
        #The doAverageMSD function should have been run first
        #so that the averageMSD variable exists
        
        #Creates the folder if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD')))
        
        directorySave = self.dn / Path('Plots/FittingsMSD')         

        #Plot MSD  fitting info
        fig0 = plt.figure(0,figsize=(12, 10))
                
        popt, pcov = curve_fit(quadratic, self.averageTimeD[:self.nbPoints], self.averageMSD[:self.nbPoints])  
        
        r_squared = calculate_rsquared(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints],popt,quadratic)

        plt.plot(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints])
        
        label = "Fitting equation: MSD = $4D_{t}t + v^2t^2$\n$D_t$ = %.4f $\mu$m$^2$/s\n$v$ = %.4f $\mu$m/s\n$v^2$ = %.4f $\mu$m$^2/s^2$\n$R^2$ = %.4f" % (popt[0]/4,np.sqrt(popt[1]),popt[1], r_squared)
        plt.plot(self.averageTimeD[:self.nbPoints],quadratic(np.asarray(self.averageTimeD[:self.nbPoints]),*popt),'--',
                 label = label)
    
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$)')
        plt.axis('tight')
        plt.legend()
        plt.title('MSD quadratic fitting for all particles')
        
        fig0.savefig(str(directorySave/('MSD_quadraticFitting_average.png')))
        fig0.savefig(str(directorySave/('MSD_quadraticFitting_average.svg')),format='svg',dpi=1200)
        plt.close()
    
        #Save data
        with open(str(directorySave/('MSD_quadraticFitting_average.csv')), 'w') as textfile:
            textfile.write("Fitting equation,MSD = 4D*t + v^2*t^2\n")
            textfile.write("D (um^2/s),%.6f\n" % (popt[0]/4))
            textfile.write("v (um/s),%.6f\n" % (np.sqrt(popt[1])))
            textfile.write("v^2 (um^2/s^2),%.6f\n" % (popt[1]))
            textfile.write("\nR^2,%.6f\n" % (r_squared))


    def doQuadraticFittingMSAD(self):
        '''Does quadratic fitting to the MSAD function'''
        
        #Creates the folder fittings if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSAD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSAD')))
            
        #Creates the folder for alpha fitting if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSAD/Quadratic'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSAD/Quadratic')))
        
        directorySave = self.dn / Path('Plots/FittingsMSAD/Quadratic')
        
        #Plot MASD with fitting info
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSAD = self.validParticles[part].MSAD
            popt, pcov = curve_fit(quadratic, timeD[:self.nbPoints], MSAD[:self.nbPoints])  
            
            self.validParticles[part].rotSpeedQuadraticFitting = np.sqrt(popt[1])
            self.validParticles[part].rotSpeedSquaredQuadraticFitting = popt[1]
            self.validParticles[part].rotDiffusionQuadraticFitting = popt[0]/4
            self.validParticles[part].tauQuadraticFitting = 1/self.validParticles[part].rotDiffusionQuadraticFitting
            
            self.validParticles[part].r_squaredRotQuadraticFitting = calculate_rsquared(timeD[:self.nbPoints],MSAD[:self.nbPoints],popt,quadratic)

            plt.plot(timeD[:self.nbPoints],MSAD[:self.nbPoints])
            
            label = "Fitting equation: MSAD = $4D_{r}t + \omega^2t^2$\n$D_r$ = %.4f rad$^2$/s\n$\omega$ = %.4f rad/s\n$\\tau_r$ = %.4f s\n$R^2$ = %.4f" % (popt[0]/4,np.sqrt(popt[1]),4/popt[0],self.validParticles[part].r_squaredRotQuadraticFitting)
            plt.plot(timeD[:self.nbPoints],quadratic(np.asarray(timeD[:self.nbPoints]),*popt),'--',
                     label = label)
        
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSAD (rad$^2$)')
            plt.axis('tight')
            plt.legend()
            plt.title('MSAD quadratic fitting for particle ' + str(part))
            
            fig0.savefig(str(directorySave/('MSAD_quadraticFitting_P'+str(part)+'.png')))
            fig0.savefig(str(directorySave/('MSAD_quadraticFitting_P'+str(part)+'.svg')),format='svg',dpi=1200)
            plt.close()

        #Save data
        for p in range(self.nbParticles):
            with open(str(directorySave/('MSAD_quadraticFitting_P'+str(p)+'.csv')), 'w') as textfile:
                textfile.write('File,'+str(self.validParticles[p].fileName)+'\n')
                textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                textfile.write("Fitting equation,MSAD = 4Dr*t + w^2*t^2\n")
                textfile.write("Dr (um^2/s),%.6f\n" % (self.validParticles[p].rotDiffusionQuadraticFitting))
                textfile.write("w (um/s),%.6f\n" % (self.validParticles[p].rotSpeedQuadraticFitting))
                textfile.write("tau (s),%.6f\n" % (self.validParticles[p].tauQuadraticFitting))
                textfile.write("\nR^2,%.6f\n" % (self.validParticles[p].r_squaredRotQuadraticFitting))

        #Save summary data
        with open(str(directorySave/'MSAD_quadraticFitting_Summary.csv'), 'w') as textfile:
            Dlist = list()
            wlist = list()
            wsquarelist = list()
            taulist = list()
            textfile.write("Fitting equation,MSAD = 4Dr*t + w^2*t^2\n")
            textfile.write('Dr (rad^2/s),')
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].rotDiffusionQuadraticFitting))
                Dlist.append(self.validParticles[p].rotDiffusionQuadraticFitting)
            textfile.write('\n')
            textfile.write("w^2 (rad^2/s^2),")
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].rotSpeedSquaredQuadraticFitting))
                wsquarelist.append(self.validParticles[p].rotSpeedSquaredQuadraticFitting)
            textfile.write('\n')
            textfile.write("w (rad/s),")
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].rotSpeedQuadraticFitting))
                wlist.append(self.validParticles[p].rotSpeedQuadraticFitting)
            textfile.write('\n')
            textfile.write("tau (s),")
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].tauQuadraticFitting))
                taulist.append(self.validParticles[p].tauQuadraticFitting)
            textfile.write('\n')
            textfile.write('\n')
            textfile.write('Average diffusion,%.6f\n' % (np.mean(Dlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(Dlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(Dlist,ddof=1)/np.sqrt(len(Dlist))))
            textfile.write('\n')
            textfile.write('Average rotational squared speed,%.6f\n' % (np.mean(wsquarelist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(wsquarelist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(wsquarelist,ddof=1)/np.sqrt(len(wsquarelist))))
            textfile.write('\n')
            textfile.write('Average rotational speed,%.6f\n' % (np.mean(wlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(wlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(wlist,ddof=1)/np.sqrt(len(wlist))))
            textfile.write('\n')
            textfile.write('Average tau,%.6f\n' % (np.mean(taulist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(taulist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(taulist,ddof=1)/np.sqrt(len(taulist))))

    def doQuadraticFittingMSADAverage(self):
        #The doAverageMASD function should have been run first
        #so that the averageMASD variable exists
        
        #Creates the folder if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSAD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSAD')))
        
        directorySave = self.dn / Path('Plots/FittingsMSAD')         

        #Plot MASD  fitting info
        fig0 = plt.figure(0,figsize=(12, 10))
                
        popt, pcov = curve_fit(quadratic, self.averageTimeD[:self.nbPoints], self.averageMSAD[:self.nbPoints])  
                
        r_squared = calculate_rsquared(self.averageTimeD[:self.nbPoints],self.averageMSAD[:self.nbPoints],popt,quadratic)

        plt.plot(self.averageTimeD[:self.nbPoints],self.averageMSAD[:self.nbPoints])
        
        label = "Fitting equation: MSAD = $4D_{r}t + \omega^2t^2$\n$D_r$ = %.4f rad$^2$/s\n$\omega$ = %.4f rad/s\n$\\tau_r$ = %.4f s\n$R^2$ = %.4f" % (popt[0]/4,np.sqrt(popt[1]),4/popt[0],r_squared)
        plt.plot(self.averageTimeD[:self.nbPoints],quadratic(np.asarray(self.averageTimeD[:self.nbPoints]),*popt),'--',
                 label = label)
    
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSAD (rad$^2$)')
        plt.axis('tight')
        plt.legend()
        plt.title('MSAD quadratic fitting for all particles')
        
        fig0.savefig(str(directorySave/('MSAD_quadraticFitting_average.png')))
        fig0.savefig(str(directorySave/('MSAD_quadraticFitting_average.svg')),format='svg',dpi=1200)
        plt.close()
    
        #Save data
        with open(str(directorySave/('MSAD_quadraticFitting_average.csv')), 'w') as textfile:
            textfile.write("Fitting equation,MSAD = 4Dr*t + w^2*t^2\n")
            textfile.write("Dr (rad^2/s),%.6f\n" % (popt[0]/4))
            textfile.write("w^2 (rad^2/s^2),%.6f\n" % (popt[1]))            
            textfile.write("w (rad/s),%.6f\n" % (np.sqrt(popt[1])))
            textfile.write("tau (s),%.6f\n" % (4/popt[0]))
            textfile.write("\nR^2,%.6f\n" % (r_squared))





    def doCubicFitting(self):
        
        #Creates the folder fittings if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD')))
            
        #Creates the folder for alpha fitting if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD/Cubic'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD/Cubic')))
        
        directorySave = self.dn / Path('Plots/FittingsMSD/Cubic')
        
        #Plot MSD with fitting info
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            popt, pcov = curve_fit(cubic, timeD[:self.nbPoints], MSD[:self.nbPoints])  
            
            
            self.validParticles[part].speedCubicFitting = np.sqrt(popt[1])
            self.validParticles[part].speedSquaredCubicFitting = popt[1]
            self.validParticles[part].tauCubicFitting = -self.validParticles[part].speedSquaredCubicFitting/(3*popt[2])
            self.validParticles[part].diffusionCubicFitting = popt[0]/4
                        
            self.validParticles[part].r_squaredCubicFitting = calculate_rsquared(timeD[:self.nbPoints],MSD[:self.nbPoints],popt,cubic)

            plt.plot(timeD[:self.nbPoints],MSD[:self.nbPoints])
            
            label = "Fitting equation: MSD = $4D_{t}t + v^2t^2 - \\frac{v^2}{3\\tau_{r}}t^3$"
            label += "\n$D_t$ = %.4f $\mu$m$^2$/s" % (self.validParticles[part].diffusionCubicFitting)
            label += "\n$v$ = %.4f $\mu$m/s" % (self.validParticles[part].speedCubicFitting)
            label += "\n$v^2$ = %.4f $\mu$m$^2/s^2$" % (self.validParticles[part].speedSquaredCubicFitting)
            label += "\n$\\tau_{r}$ = %.4f $s$" % (self.validParticles[part].tauCubicFitting)
            label += "\n$R^2$ = %.4f" % (self.validParticles[part].r_squaredCubicFitting)
            
            plt.plot(timeD[:self.nbPoints],cubic(np.asarray(timeD[:self.nbPoints]),*popt),'--',
                     label = label)
        
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSD ($\mu$m$^2$)')
            plt.axis('tight')
            plt.legend()
            plt.title('MSD cubic fitting for particle ' + str(part))
            
            fig0.savefig(str(directorySave/('MSD_cubicFitting_P'+str(part)+'.png')))
            fig0.savefig(str(directorySave/('MSD_cubicFitting_P'+str(part)+'.svg')),format='svg',dpi=1200)
            plt.close()

        #Save data
        for p in range(self.nbParticles):
            with open(str(directorySave/('MSD_cubicFitting_P'+str(p)+'.csv')), 'w') as textfile:
                textfile.write('File,'+str(self.validParticles[p].fileName)+'\n')
                textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                textfile.write("Fitting equation,MSD = 4D*t + v^2*t^2 + -(v^2/3*tau)*t^3\n")
                textfile.write("D (um^2/s),%.6f\n" % (self.validParticles[p].diffusionCubicFitting))
                textfile.write("v (um/s),%.6f\n" % (self.validParticles[p].speedCubicFitting))
                textfile.write("v^2 (um^2/s^2),%.6f\n" % (self.validParticles[p].speedSquaredCubicFitting))
                textfile.write("tau (s),%.6f\n" % (self.validParticles[p].tauCubicFitting))
                textfile.write("\nR^2,%.6f\n" % (self.validParticles[p].r_squaredCubicFitting))

        #Save summary data
        with open(str(directorySave/('MSD_cubicFitting_Summary.csv')), 'w') as textfile:
            Dlist = list()
            vlist = list()
            vsquarelist = list()
            taulist = list()
            textfile.write("Fitting equation,MSD = 4D*t + v^2*t^2 + -(v^2/3*tau)*t^3\n")
            textfile.write('D (um^2/s),')
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].diffusionCubicFitting))
                Dlist.append(self.validParticles[p].diffusionCubicFitting)
            textfile.write('\n')
            textfile.write("v (um/s),")
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].speedCubicFitting))
                vlist.append(self.validParticles[p].speedCubicFitting)
            textfile.write('\n')
            textfile.write("v^2 (um^2/s^2),")
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].speedSquaredCubicFitting))
                vsquarelist.append(self.validParticles[p].speedSquaredCubicFitting)
            textfile.write('\n')
            textfile.write("tau (s),")
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].tauCubicFitting))
                taulist.append(self.validParticles[p].tauCubicFitting)
            textfile.write('\n')
            textfile.write('\n')
            textfile.write('Average diffusion,%.6f\n' % (np.mean(Dlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(Dlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(Dlist,ddof=1)/np.sqrt(len(Dlist))))
            textfile.write('\n')
            textfile.write('Average speed,%.6f\n' % (np.mean(vlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(vlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(vlist,ddof=1)/np.sqrt(len(vlist))))
            textfile.write('\n')
            textfile.write('Average speed squared,%.6f\n' % (np.mean(vsquarelist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(vsquarelist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(vsquarelist,ddof=1)/np.sqrt(len(vsquarelist))))
            textfile.write('\n')
            textfile.write('Average rotational time,%.6f\n' % (np.mean(taulist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(taulist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(taulist,ddof=1)/np.sqrt(len(taulist))))

    def doCubicFittingAverage(self):
        #The doAverageMSD function should have been run first
        #so that the averageMSD variable exists
        
        #Creates the folder if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD')))
        
        directorySave = self.dn / Path('Plots/FittingsMSD')         

        #Plot MSD  fitting info
        fig0 = plt.figure(0,figsize=(12, 10))
                
        popt, pcov = curve_fit(cubic, self.averageTimeD[:self.nbPoints], self.averageMSD[:self.nbPoints])  
        
        r_squared = calculate_rsquared(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints],popt,cubic)

        plt.plot(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints])
        
        label = "Fitting equation: MSD = $4D_{t}t + v^2t^2 - \\frac{v^2}{3\\tau_{r}}t^3$"
        label += "\n$D_t$ = %.4f $\mu$m$^2$/s" % (popt[0]/4)
        label += "\n$v$ = %.4f $\mu$m/s" % (np.sqrt(popt[1]))
        label += "\n$v^2$ = %.4f $\mu$m$^2/s^2$" % (popt[1])
        label += "\n$\\tau_{r}$ = %.4f $s$" % (-popt[1]/(3*popt[2]))
        label += "\n$R^2$ = %.4f" % (r_squared)

        plt.plot(self.averageTimeD[:self.nbPoints],cubic(np.asarray(self.averageTimeD[:self.nbPoints]),*popt),'--',
                 label = label)
    
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$)')
        plt.axis('tight')
        plt.legend()
        plt.title('MSD cubic fitting for all particles')
        
        fig0.savefig(str(directorySave/('MSD_cubicFitting_average.png')))
        fig0.savefig(str(directorySave/('MSD_cubicFitting_average.svg')),format='svg',dpi=1200)
        plt.close()
    
        #Save data
        with open(str(directorySave/('MSD_cubicFitting_average.csv')), 'w') as textfile:
            textfile.write("Fitting equation,MSD = 4D*t + v^2*t^2 + -(v^2/3*tau)*t^3\n")
            textfile.write("D (um^2/s),%.6f\n" % (popt[0]/4))
            textfile.write("v (um/s),%.6f\n" % (np.sqrt(popt[1])))
            textfile.write("v^2 (um^2/s^2),%.6f\n" % (popt[1]))
            textfile.write("tau (s),%.6f\n" % (-popt[1]/(3*popt[2])))
            textfile.write("\nR^2,%.6f\n" % (r_squared))



    def doFourthOrderFitting(self):
        
        #Creates the folder fittings if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD')))
            
        #Creates the folder for alpha fitting if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD/Fourth order'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD/Fourth order')))
        
        directorySave = self.dn / Path('Plots/FittingsMSD/Fourth order')
        
        #Plot MSD with fitting info
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            popt, pcov = curve_fit(fourthOrder, timeD[:self.nbPoints], MSD[:self.nbPoints])  
            
            
            self.validParticles[part].speedFourthOrderFitting = np.sqrt(popt[1])
            self.validParticles[part].speedSquaredFourthOrderFitting = popt[1]
            self.validParticles[part].tauFourthOrderFitting = -self.validParticles[part].speedSquaredFourthOrderFitting/(3*popt[2])
            self.validParticles[part].diffusionFourthOrderFitting = popt[0]/4
            
            self.validParticles[part].r_squaredFourthOrderFitting = calculate_rsquared(timeD[:self.nbPoints],MSD[:self.nbPoints],popt,fourthOrder)

            plt.plot(timeD[:self.nbPoints],MSD[:self.nbPoints])
            
            label = "Fitting equation: MSD = $4D_{t}t + v^2t^2 - \\frac{v^2}{3\\tau_{r}}t^3 + \\frac{v^2}{12\\tau_{r}^2}t^4$"
            label += "\n$D_t$ = %.4f $\mu$m$^2$/s" % (self.validParticles[part].diffusionFourthOrderFitting)
            label += "\n$v$ = %.4f $\mu$m/s" % (self.validParticles[part].speedFourthOrderFitting)
            label += "\n$v^2$ = %.4f $\mu$m$^2/s^2$" % (self.validParticles[part].speedSquaredFourthOrderFitting)
            label += "\n$\\tau_{r}$ = %.4f $s$" % (self.validParticles[part].tauFourthOrderFitting)
            label += "\n$R^2$ = %.4f" % (self.validParticles[part].r_squaredFourthOrderFitting)

            
            plt.plot(timeD[:self.nbPoints],fourthOrder(np.asarray(timeD[:self.nbPoints]),*popt),'--',
                     label = label)
        
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSD ($\mu$m$^2$)')
            plt.axis('tight')
            plt.legend()
            plt.title('MSD fourth order fitting for particle ' + str(part))
            
            fig0.savefig(str(directorySave/('MSD_fourthOrderFitting_P'+str(part)+'.png')))
            fig0.savefig(str(directorySave/('MSD_fourthOrderFitting_P'+str(part)+'.svg')),format='svg',dpi=1200)
            plt.close()

        #Save data
        for p in range(self.nbParticles):
            with open(str(directorySave/('MSD_fourthOrderFitting_P'+str(p)+'.csv')), 'w') as textfile:
                textfile.write('File,'+str(self.validParticles[p].fileName)+'\n')
                textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                textfile.write("Fitting equation,MSD = 4D*t + v^2*t^2 + -(v^2/3*tau)*t^3 + (v^2/12*tau^2)*t^4\n")
                textfile.write("D (um^2/s),%.6f\n" % (self.validParticles[p].diffusionFourthOrderFitting))
                textfile.write("v (um/s),%.6f\n" % (self.validParticles[p].speedFourthOrderFitting))
                textfile.write("v^2 (um^2/s^2),%.6f\n" % (self.validParticles[p].speedSquaredFourthOrderFitting))
                textfile.write("tau (s),%.6f\n" % (self.validParticles[p].tauFourthOrderFitting))
                textfile.write("\nR^2,%.6f\n" % (self.validParticles[p].r_squaredFourthOrderFitting))

        #Save summary data
        with open(str(directorySave/('MSD_fourthOrderFitting_Summary.csv')), 'w') as textfile:
            Dlist = list()
            vlist = list()
            vsquarelist = list()
            taulist = list()
            textfile.write("Fitting equation,MSD = 4D*t + v^2*t^2 + -(v^2/3*tau)*t^3 + (v^2/12*tau^2)*t^4\n")
            textfile.write('D (um^2/s),')
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].diffusionFourthOrderFitting))
                Dlist.append(self.validParticles[p].diffusionFourthOrderFitting)
            textfile.write('\n')
            textfile.write("v (um/s),")
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].speedFourthOrderFitting))
                vlist.append(self.validParticles[p].speedFourthOrderFitting)
            textfile.write('\n')
            textfile.write("v^2 (um^2/s^2),")
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].speedSquaredFourthOrderFitting))
                vsquarelist.append(self.validParticles[p].speedSquaredFourthOrderFitting)
            textfile.write('\n')
            textfile.write("tau (s),")
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].tauFourthOrderFitting))
                taulist.append(self.validParticles[p].tauFourthOrderFitting)
            textfile.write('\n')
            textfile.write('\n')
            textfile.write('Average diffusion,%.6f\n' % (np.mean(Dlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(Dlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(Dlist,ddof=1)/np.sqrt(len(Dlist))))
            textfile.write('\n')
            textfile.write('Average speed,%.6f\n' % (np.mean(vlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(vlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(vlist,ddof=1)/np.sqrt(len(vlist))))
            textfile.write('\n')
            textfile.write('Average speed squared,%.6f\n' % (np.mean(vsquarelist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(vsquarelist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(vsquarelist,ddof=1)/np.sqrt(len(vsquarelist))))
            textfile.write('\n')
            textfile.write('Average rotational time,%.6f\n' % (np.mean(taulist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(taulist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(taulist,ddof=1)/np.sqrt(len(taulist))))

    def doFourthOrderFittingAverage(self):
        #The doAverageMSD function should have been run first
        #so that the averageMSD variable exists
        
        #Creates the folder if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD')))
        
        directorySave = self.dn / Path('Plots/FittingsMSD')         

        #Plot MSD  fitting info
        fig0 = plt.figure(0,figsize=(12, 10))
                
        popt, pcov = curve_fit(fourthOrder, self.averageTimeD[:self.nbPoints], self.averageMSD[:self.nbPoints])  
        
        r_squared = calculate_rsquared(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints],popt,fourthOrder)

        plt.plot(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints])
        
        label = "Fitting equation: MSD = $4D_{t}t + v^2t^2 - \\frac{v^2}{3\\tau_{r}}t^3 + \\frac{v^2}{12\\tau_{r}^2}t^4$"
        label += "\n$D_t$ = %.4f $\mu$m$^2$/s" % (popt[0]/4)
        label += "\n$v$ = %.4f $\mu$m/s" % (np.sqrt(popt[1]))
        label += "\n$v^2$ = %.4f $\mu$m$^2/s^2$" % (popt[1])
        label += "\n$\\tau_{r}$ = %.4f $s$" % (-popt[1]/(3*popt[2]))
        label += "\n$R^2$ = %.4f" % (r_squared)

        plt.plot(self.averageTimeD[:self.nbPoints],fourthOrder(np.asarray(self.averageTimeD[:self.nbPoints]),*popt),'--',
                 label = label)
    
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$)')
        plt.axis('tight')
        plt.legend()
        plt.title('MSD fourth order fitting for all particles')
        
        fig0.savefig(str(directorySave/('MSD_fourthOrderFitting_average.png')))
        fig0.savefig(str(directorySave/('MSD_fourthOrderFitting_average.svg')),format='svg',dpi=1200)
        plt.close()
    
        #Save data
        with open(str(directorySave/('MSD_fourthOrderFitting_average.csv')), 'w') as textfile:
            textfile.write("Fitting equation,MSD = 4D*t + v^2*t^2 + -(v^2/3*tau)*t^3 + (v^2/12*tau^2)*t^4\n")
            textfile.write("D (um^2/s),%.6f\n" % (popt[0]/4))
            textfile.write("v (um/s),%.6f\n" % (np.sqrt(popt[1])))
            textfile.write("v^2 (um^2/s^2),%.6f\n" % (popt[1]))
            textfile.write("tau (s),%.6f\n" % (-popt[1]/(3*popt[2])))
            textfile.write("\nR^2,%.6f\n" % (r_squared))




    def doFullEquationFitting(self):
        
        #Creates the folder fittings if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD')))
            
        #Creates the folder for alpha fitting if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD/Full equation'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD/Full equation')))
        
        directorySave = self.dn / Path('Plots/FittingsMSD/Full equation')
        
        #Plot MSD with fitting info
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            try:
                popt, pcov = curve_fit(fullEquation, timeD[:self.nbPoints], MSD[:self.nbPoints])  
            except:
                print('Full equation fitting failed.')
                popt = [np.nan]*4
            
            self.validParticles[part].speedFullEquationFitting = popt[0]
            self.validParticles[part].tauFullEquationFitting = popt[1]
            self.validParticles[part].diffusionFullEquationFitting = popt[2]
            
            self.validParticles[part].r_squaredFullEquationFitting = calculate_rsquared(timeD[:self.nbPoints],MSD[:self.nbPoints],popt,fullEquation)

            plt.plot(timeD[:self.nbPoints],MSD[:self.nbPoints])
            
            label = "Fitting equation: MSD = $4D_{t}t + 2v^2\\tau_{r}^2 (\\frac{t}{\\tau_{r}} + \\exp(-t/\\tau_{r}) - 1)$"
            label += "\n$D_t$ = %.4f $\mu$m$^2$/s" % (self.validParticles[part].diffusionFullEquationFitting)
            label += "\n$v$ = %.4f $\mu$m/s" % (self.validParticles[part].speedFullEquationFitting)
            label += "\n$\\tau_{r}$ = %.4f $s$" % (self.validParticles[part].tauFullEquationFitting)
            label += "\n$R^2$ = %.4f" % (self.validParticles[part].r_squaredFullEquationFitting)

          
            plt.plot(timeD[:self.nbPoints],fullEquation(np.asarray(timeD[:self.nbPoints]),*popt),'--',
                     label = label)
        
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSD ($\mu$m$^2$)')
            plt.axis('tight')
            plt.legend()
            plt.title('MSD full equation fitting for particle ' + str(part))
            
            fig0.savefig(str(directorySave/('MSD_fullEquationFitting_P'+str(part)+'.png')))
            fig0.savefig(str(directorySave/('MSD_fullEquationFitting_P'+str(part)+'.svg')),format='svg',dpi=1200)
            plt.close()

        #Save data
        for p in range(self.nbParticles):
            with open(str(directorySave/('MSD_fullEquationFitting_P'+str(p)+'.csv')), 'w') as textfile:
                textfile.write('File,'+str(self.validParticles[p].fileName)+'\n')
                textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                textfile.write("Fitting equation,MSD = 4D*t + 2(v^2)*(tau^2)*(t/tau + exp(-t/tau) - 1)\n")
                textfile.write("D (um^2/s),%.6f\n" % (self.validParticles[p].diffusionFullEquationFitting))
                textfile.write("v (um/s),%.6f\n" % (self.validParticles[p].speedFullEquationFitting))
                textfile.write("tau (s),%.6f\n" % (self.validParticles[p].tauFullEquationFitting))
                textfile.write("\nR^2,%.6f\n" % (self.validParticles[p].r_squaredFullEquationFitting))

        #Save summary data
        with open(str(directorySave/('MSD_fullEquationFitting_Summary.csv')), 'w') as textfile:
            Dlist = list()
            vlist = list()
            taulist = list()
            textfile.write("Fitting equation,MSD = 4D*t + 2(v^2)*(tau^2)*(t/tau + exp(-t/tau) - 1)\n")
            textfile.write('D (um^2/s),')
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].diffusionFullEquationFitting))
                Dlist.append(self.validParticles[p].diffusionFullEquationFitting)
            textfile.write('\n')
            textfile.write("v (um/s),")
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].speedFullEquationFitting))
                vlist.append(self.validParticles[p].speedFullEquationFitting)
            textfile.write('\n')
            textfile.write("tau (s),")
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].tauFullEquationFitting))
                taulist.append(self.validParticles[p].tauFullEquationFitting)
            textfile.write('\n')
            textfile.write('\n')
            textfile.write('Average diffusion,%.6f\n' % (np.mean(Dlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(Dlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(Dlist,ddof=1)/np.sqrt(len(Dlist))))
            textfile.write('\n')
            textfile.write('Average speed,%.6f\n' % (np.mean(vlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(vlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(vlist,ddof=1)/np.sqrt(len(vlist))))
            textfile.write('\n')
            textfile.write('Average rotational time,%.6f\n' % (np.mean(taulist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(taulist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(taulist,ddof=1)/np.sqrt(len(taulist))))

    def doFullEquationFittingAverage(self):
        #The doAverageMSD function should have been run first
        #so that the averageMSD variable exists
        
        #Creates the folder if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD')))
        
        directorySave = self.dn / Path('Plots/FittingsMSD')         

        #Plot MSD  fitting info
        fig0 = plt.figure(0,figsize=(12, 10))
                
        try:
            popt, pcov = curve_fit(fullEquation, self.averageTimeD[:self.nbPoints], self.averageMSD[:self.nbPoints])  
        except:
            print('Full equation fitting failed.')
            popt = [np.nan]*4
        
        r_squared = calculate_rsquared(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints],popt,fullEquation)

        plt.plot(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints])
        
        label = "Fitting equation: MSD = $4D_{t}t + 2v^2\\tau_{r}^2 (\\frac{t}{\\tau_{r}} + \\exp(-t/\\tau_{r}) - 1)$"
        label += "\n$D_t$ = %.4f $\mu$m$^2$/s" % (popt[2])
        label += "\n$v$ = %.4f $\mu$m/s" % (popt[0])
        label += "\n$\\tau_{r}$ = %.4f $s$" % (popt[1])
        label += "\n$R^2$ = %.4f" % (r_squared)

        plt.plot(self.averageTimeD[:self.nbPoints],fullEquation(np.asarray(self.averageTimeD[:self.nbPoints]),*popt),'--',
                 label = label)
    
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$)')
        plt.axis('tight')
        plt.legend()
        plt.title('MSD full equation fitting for all particles')
        
        fig0.savefig(str(directorySave/('MSD_fullEquationFitting_average.png')))
        fig0.savefig(str(directorySave/('MSD_fullEquationFitting_average.svg')),format='svg',dpi=1200)
        plt.close()
    
        #Save data
        with open(str(directorySave/('MSD_fullEquationFitting_average.csv')), 'w') as textfile:
            textfile.write("Fitting equation,MSD = 4D*t + 2(v^2)*(tau^2)*(t/tau + exp(-t/tau) - 1)\n")
            textfile.write("D (um^2/s),%.6f\n" % (popt[2]))
            textfile.write("v (um/s),%.6f\n" % (popt[0]))
            textfile.write("tau (s),%.6f\n" % (popt[1]))
            textfile.write("\nR^2,%.6f\n" % (r_squared))







    def doAlphaFitting(self):
        
        
        #Creates the folder fittings if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD')))
            
        #Creates the folder for alpha fitting if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD/Logarithmic'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD/Logarithmic')))
        
        directorySave = self.dn / Path('Plots/FittingsMSD/Logarithmic')
        
        #Plot MSD in log-log with fitting info
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            popt, pcov = curve_fit(linear, np.log(timeD[:self.nbPoints]), np.log(MSD[:self.nbPoints]))  
            
            self.validParticles[part].alphaFitting = popt[1]
            self.validParticles[part].diffusionAlphaFitting = np.exp(popt[0])
            
            self.validParticles[part].r_squaredAlphaFitting = calculate_rsquared(np.log(timeD[:self.nbPoints]),np.log(MSD[:self.nbPoints]),popt,linear)

            plt.plot(timeD[:self.nbPoints],MSD[:self.nbPoints])
            
            label = "Fitting equation: MSD = $Dt^{\\alpha}$\n$\\alpha$ = %.2f\n$D$ = %.6f\n$R^2$ = %.4f" % (popt[1],np.exp(popt[0]),self.validParticles[part].r_squaredAlphaFitting)
            plt.plot(timeD[:self.nbPoints],powerLaw(np.asarray(timeD[:self.nbPoints]),np.exp(popt[0]),popt[1]),'--',
                     label = label)
        
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSD ($\mu$m$^2$)')
            plt.yscale('log')
            plt.xscale('log')
            plt.axis('tight')
            plt.legend()
            plt.title('MSD exponent for particle ' + str(part))
            
            fig0.savefig(str(directorySave/('MSD_logFitting_P'+str(part)+'.png')))
            fig0.savefig(str(directorySave/('MSD_logFitting_P'+str(part)+'.svg')),format='svg',dpi=1200)
            plt.close()

        #Save data
        for p in range(self.nbParticles):
            with open(str(directorySave/('MSD_logFitting_P'+str(p)+'.csv')), 'w') as textfile:
                textfile.write('File,'+str(self.validParticles[p].fileName)+'\n')
                textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                textfile.write("Fitting equation,MSD = D*t^alpha\n")
                textfile.write("D (um^2/s^alpha),%.6f\n" % (self.validParticles[p].diffusionAlphaFitting))
                textfile.write("Alpha,%.6f\n" % (self.validParticles[p].alphaFitting))
                textfile.write("\nR^2,%.6f\n" % (self.validParticles[p].r_squaredAlphaFitting))

        #Save summary data
        with open(str(directorySave/'MSD_logFitting_Summary.csv'), 'w') as textfile:
            Dlist = list()
            alphalist = list()
            textfile.write("Fitting equation,MSD = D*t^alpha\n")
            textfile.write('D (um^2/s^alpha),')
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].diffusionAlphaFitting))
                Dlist.append(self.validParticles[p].diffusionAlphaFitting)
            textfile.write('\n')
            textfile.write("Alpha,")
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].alphaFitting))
                alphalist.append(self.validParticles[p].alphaFitting)
            textfile.write('\n')
            textfile.write('\n')
            textfile.write('Average diffusion,%.6f\n' % (np.mean(Dlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(Dlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(Dlist,ddof=1)/np.sqrt(len(Dlist))))
            textfile.write('\n')
            textfile.write('Average alpha,%.6f\n' % (np.mean(alphalist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(alphalist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(alphalist,ddof=1)/np.sqrt(len(alphalist))))


    def doAlphaFittingAverage(self):
        #The doAverageMSD function should have been run first
        #so that the averageMSD variable exists
        
        #Creates the folder if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/FittingsMSD'))):
            os.mkdir(str(self.dn / Path('Plots/FittingsMSD')))
        
        directorySave = self.dn / Path('Plots/FittingsMSD')         

        #Plot MSD in log-log with fitting info
        fig0 = plt.figure(0,figsize=(12, 10))
                
        popt, pcov = curve_fit(linear, np.log(self.averageTimeD[:self.nbPoints]), np.log(self.averageMSD[:self.nbPoints]))  
        
        r_squared = calculate_rsquared(np.log(self.averageTimeD[:self.nbPoints]),np.log(self.averageMSD[:self.nbPoints]),popt,linear)

        plt.plot(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints])
        
        label = "Fitting equation: MSD = $Dt^{\\alpha}$\n$\\alpha$ = %.2f\n$D$ = %.6f\n$R^2$ = %.4f" % (popt[1],np.exp(popt[0]),r_squared)
        plt.plot(self.averageTimeD[:self.nbPoints],powerLaw(np.asarray(self.averageTimeD[:self.nbPoints]),np.exp(popt[0]),popt[1]),'--',
                 label = label)
    
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$)')
        plt.yscale('log')
        plt.xscale('log')
        plt.axis('tight')
        plt.legend()
        plt.title('MSD exponent for all particles')
        
        fig0.savefig(str(directorySave/'MSD_logFitting_average.png'))
        fig0.savefig(str(directorySave/'MSD_logFitting_average.svg'),format='svg',dpi=1200)
        plt.close()

        #Save data
        with open(str(directorySave/'MSD_logFitting_average.csv'), 'w') as textfile:
            textfile.write("Fitting equation,MSD = D*t^alpha\n")
            textfile.write("D (um^2/s^alpha),%.6f\n" % (popt[1]))
            textfile.write("Alpha,%.6f\n" % (np.exp(popt[0])))        
            textfile.write("\nR^2,%.6f\n" % (r_squared))



                
                
    def doLocalAlpha(self):
        
        if not os.path.exists(str(self.dn / Path('Plots/Exponent'))):
            os.mkdir(str(self.dn / Path('Plots/Exponent')))
        
        directorySave = self.dn / Path('Plots/Exponent')


        #Plot all alphas (long version)
        fig0 = plt.figure(0,figsize=(12, 10))
        for part in range(self.nbParticles):
        
            self.validParticles[part].calculateLocalAlpha()
            
            timeD = self.validParticles[part].timeD[:-1] #By default, alpha has one less value
            alpha = self.validParticles[part].alpha
            
            plt.plot(timeD,alpha)
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD exponent')
        plt.axis('tight')
        plt.title('MSD exponent for all particles (long)')
        
        fig0.savefig(str(directorySave/'alpha_long_all.png'))
        fig0.savefig(str(directorySave/'alpha_long_all.svg'),format='svg',dpi=1200)
        plt.close()


        #Plot all alphas (short version, 1/10 of data points)
        fig0 = plt.figure(0,figsize=(12, 10))
        for part in range(self.nbParticles):
        
            timeD = self.validParticles[part].timeD[:-1] #By default, alpha has one less value
            alpha = self.validParticles[part].alpha
            
            plt.plot(timeD[:self.nbPoints],alpha[:self.nbPoints])
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD exponent')
        plt.axis('tight')
        plt.title('MSD exponent for all particles (short)')
        
        fig0.savefig(str(directorySave/'alpha_short_all.png'))
        fig0.savefig(str(directorySave/'alpha_short_all.svg'),format='svg',dpi=1200)
        plt.close()
        
        #Creates the folder for individual particles if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/Exponent/Individual'))):
            os.mkdir(str(self.dn / Path('Plots/Exponent/Individual')))
        
        directorySave = self.dn / Path('Plots/Exponent/Individual')
        
        #Plots all alphas one by one (short version)
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD[:-1] #By default, alpha has one less value
            alpha = self.validParticles[part].alpha
            
            plt.plot(timeD[:self.nbPoints],alpha[:self.nbPoints])
        
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSD exponent')
            plt.axis('tight')
            plt.title('MSD exponent for particle ' + str(part))
            
            fig0.savefig(str(directorySave/('alpha_short_P'+str(part)+'.png')))
            fig0.savefig(str(directorySave/('alpha_short_P'+str(part)+'.svg')),format='svg',dpi=1200)
            plt.close()
                

        #Plots all alphas one by one (long version)
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD[:-1] #By default, alpha has one less value
            alpha = self.validParticles[part].alpha
            
            plt.plot(timeD,alpha)
                
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSD exponent')
            plt.axis('tight')
            plt.title('MSD exponent for particle ' + str(part))
            
            fig0.savefig(str(directorySave/('alpha_long_P'+str(part)+'.png')))
            fig0.savefig(str(directorySave/('alpha_long_P'+str(part)+'.svg')),format='svg',dpi=1200)
            plt.close()
        
        #Save data
        for p in range(self.nbParticles):
            alpha = self.validParticles[p].alpha
            timeD = self.validParticles[p].timeD[:-1] #By default, alpha has one less value
            if len(alpha) > 0:
                with open(str(directorySave/('alpha_P'+str(p)+'.csv')), 'w') as textfile:
                    textfile.write('File,'+str(self.validParticles[p].fileName)+'\n')
                    textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                    textfile.write('TimeD (s),')
                    for t in timeD:
                        textfile.write("%.2f," % (t))
                    textfile.write('\n')
                    textfile.write('MSD exponent,')
                    for m in alpha:
                        textfile.write(("%.6f," % (m)))

    def doAverageAlpha(self):
        #The doAverageMSD function should have been run first
        #so that the averageMSD variable exists
        
        #Creates the folder alpha if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/Exponent'))):
            os.mkdir(str(self.dn / Path('Plots/Exponent')))
        
        directorySave = self.dn / Path('Plots/Exponent')
                
        
        self.averageAlpha = list()
        for i in range(len(self.averageMSD)-1):
            num = np.log(self.averageMSD[i+1])-np.log(self.averageMSD[i])
            den = np.log(self.averageTimeD[i+1])-np.log(self.averageTimeD[i])
            self.averageAlpha.append(num/den)

        #Plot average alpha, long version
        fig0 = plt.figure(0,figsize=(12, 10))
        
        plt.plot(self.averageTimeD[:-1],self.averageAlpha)
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD exponent')
        plt.axis('tight')
        plt.title('Average MSD exponent (long)')
        
        fig0.savefig(str(directorySave/'alpha_average_long.png'))
        fig0.savefig(str(directorySave/'alpha_average_long.svg'),format='svg',dpi=1200)
        plt.close()
        
        #Plot average alpha, short version
        fig0 = plt.figure(0,figsize=(12, 10))
                
        plt.plot(self.averageTimeD[:self.nbPoints],self.averageAlpha[:self.nbPoints])
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD exponent')
        plt.axis('tight')
        plt.title('Average MSD exponent (short)')
        
        fig0.savefig(str(directorySave/'alpha_average_short.png'))
        fig0.savefig(str(directorySave/'alpha_average_short.svg'),format='svg',dpi=1200)
        plt.close()  
        
        with open(str(directorySave/'alpha_average.csv'), 'w') as textfile:
            textfile.write('TimeD (s),')
            for t in self.averageTimeD:
                textfile.write("%.2f," % (t))
            textfile.write('\n')
            textfile.write('MSD exponent,')
            for m in self.averageAlpha:
                textfile.write(("%.6f," % (m)))     
                


    def doInstVel(self):
        
        if not os.path.exists(str(self.dn / Path('Plots/InstantaneousVelocity'))):
            os.mkdir(str(self.dn / Path('Plots/InstantaneousVelocity')))
        
        directorySave = self.dn / Path('Plots/InstantaneousVelocity')
        
        #Plot all velocities
        fig0 = plt.figure(0,figsize=(12, 10))
        
        for part in range(self.nbParticles):
            
            self.validParticles[part].calculateInstVel()
            timeD = self.validParticles[part].timeD
            v = self.validParticles[part].v[:-1] #The speed has one value less
            
            if len(v) > 0:
                plt.plot(timeD,v)
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('Instantaneous velocity ($\mu$m/s)')
        plt.axis('tight')
        plt.title('Instantaneous velocity for all particles')
        
        fig0.savefig(str(directorySave/'instVel_all.png'))
        fig0.savefig(str(directorySave/'instVel_all.svg'),format='svg',dpi=1200)
        plt.close()

        #Creates the folder for individual particles if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/InstantaneousVelocity/Individual'))):
            os.mkdir(str(self.dn / Path('Plots/InstantaneousVelocity/Individual')))
        
        directorySave = self.dn / Path('Plots/InstantaneousVelocity/Individual')
        
        #Plots all InstantaneousVelocity one by one
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            v = self.validParticles[part].v[:-1] #The speed has one value less
            
            if len(v) > 0: 
                averageV = np.mean(v)
                stdV = np.std(v)
                
                label = "v$_{average}$ = %.2f $\pm$ %.2f (SD)" % (averageV,stdV)
                
                plt.plot(timeD,v,label = label)
            
                plt.xlabel('$\Delta$t (s)')
                plt.ylabel('Instantaneous velocity ($\mu$m/s)')
                plt.axis('tight')
                plt.title('Instantaneous velocity for particle ' + str(part))
                plt.legend()
                
                fig0.savefig(str(directorySave/('instVel_P'+str(part)+'.png')))
                fig0.savefig(str(directorySave/('instVel_P'+str(part)+'.svg')),format='svg',dpi=1200)
                plt.close()
        
        #Save data
        for p in range(self.nbParticles):
            v = self.validParticles[p].v[:-1]
            timeD = self.validParticles[p].timeD
            if len(v) > 0:
                averageV = np.mean(v)
                stdV = np.std(v)
                with open(str(directorySave/('instVel_P'+str(p)+'.csv')), 'w') as textfile:
                    textfile.write('File,'+str(self.validParticles[p].fileName)+'\n')
                    textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                    textfile.write('Average velocity,%.6f\n' % (averageV))
                    textfile.write('STD velocity,%.6f\n' % (stdV))
                    textfile.write('Time (s),')
                    for t in timeD:
                        textfile.write("%.2f," % (t))
                    textfile.write('\n')
                    textfile.write('Instantaneous velocity (um/s),')
                    for m in v:
                        textfile.write(("%.6f," % (m)))

    def doAverageAutoCor(self):
        
        #Creates the folder MSD if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/Autocorrelation'))):
            os.mkdir(str(self.dn / Path('Plots/Autocorrelation')))
        
        directorySave = self.dn / Path('Plots/Autocorrelation') 
                
        #Plot average MSD, long version
        fig0 = plt.figure(0,figsize=(12, 10))
        
        plt.plot(self.averageTimeDAutoCor,self.averageAutoCor)
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('Angular autocorrelation')
        plt.axis('tight')
        plt.title('Average autocorrelation function (long)')
        
        fig0.savefig(str(directorySave/'AutoCor_average_long.png'))
        fig0.savefig(str(directorySave/'AutoCor_average_long.svg'),format='svg',dpi=1200)
        plt.close()
        
        #Plot average MSD, short version
        fig0 = plt.figure(0,figsize=(12, 10))
                
        plt.plot(self.averageTimeDAutoCor[:self.nbPoints],self.averageAutoCor[:self.nbPoints])
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('Angular autocorrelation')
        plt.axis('tight')
        plt.title('Average autocorrelation function (short)')
        
        fig0.savefig(str(directorySave/'AutoCor_average_short.png'))
        fig0.savefig(str(directorySave/'AutoCor_average_short.svg'),format='svg',dpi=1200)
        plt.close()  
        
        with open(str(directorySave/'AutoCor_average.csv'), 'w') as textfile:
            textfile.write('TimeD (s),')
            for t in self.averageTimeDAutoCor:
                textfile.write("%.2f," % (t))
            textfile.write('\n')
            textfile.write('Angular autocorrelation,')
            for m in self.averageAutoCor:
                textfile.write(("%.6f," % (m)))        


    def doAutoCor(self):
                
        if not os.path.exists(str(self.dn / Path('Plots/Autocorrelation'))):
            os.mkdir(str(self.dn / Path('Plots/Autocorrelation')))
        
        directorySave = self.dn / Path('Plots/Autocorrelation')
        
        #Plot all autocorrelations (long version)
        fig0 = plt.figure(0,figsize=(12, 10))
        for part in range(self.nbParticles):
        
            timeD = self.validParticles[part].timeD
            autoCor = self.validParticles[part].autoCor
            
            if len(autoCor) > 0:
                plt.plot(timeD,autoCor)
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('Angular autocorrelation')
        plt.axis('tight')
        plt.title('Autocorrelation function for all particles (long)')
        
        fig0.savefig(str(directorySave/'autocor_long_all.png'))
        fig0.savefig(str(directorySave/'autocor_long_all.svg'),format='svg',dpi=1200)
        plt.close()


        #Plot all autocorrelations (short version, 1/10 of data points)
        fig0 = plt.figure(0,figsize=(12, 10))
        for part in range(self.nbParticles):
        
            timeD = self.validParticles[part].timeD
            autoCor = self.validParticles[part].autoCor
            
            if len(autoCor) > 0:
                plt.plot(timeD[:self.nbPoints],autoCor[:self.nbPoints])
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('Angular autocorrelation')
        plt.axis('tight')
        plt.title('Autocorrelation function for all particles (short)')
        
        fig0.savefig(str(directorySave/'autocor_short_all.png'))
        fig0.savefig(str(directorySave/'autocor_short_all.svg'),format='svg',dpi=1200)
        plt.close()
        
        
        #Creates the folder for individual particles if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/Autocorrelation/Individual'))):
            os.mkdir(str(self.dn / Path('Plots/Autocorrelation/Individual')))
        
        directorySave = self.dn / Path('Plots/Autocorrelation/Individual')
        
        #Plots all Autocorrelation one by one (short version)
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            autoCor = self.validParticles[part].autoCor
            
            if len(autoCor) > 0:
                plt.plot(timeD[:self.nbPoints],autoCor[:self.nbPoints])
            
                plt.xlabel('$\Delta$t (s)')
                plt.ylabel('Angular autocorrelation')
                plt.axis('tight')
                plt.title('Autocorrelation for particle ' + str(part))
                
                fig0.savefig(str(directorySave/('autocor_short_P'+str(part)+'.png')))
                fig0.savefig(str(directorySave/('autocor_short_P'+str(part)+'.svg')),format='svg',dpi=1200)
                plt.close()
                

        #Plots all Autocorrelation one by one (long version)
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            autoCor = self.validParticles[part].autoCor
            
            if len(autoCor) > 0:
                plt.plot(timeD,autoCor)
                
                plt.xlabel('$\Delta$t (s)')
                plt.ylabel('Angular autocorrelation')
                plt.axis('tight')
                plt.title('Autocorrelation for particle ' + str(part))
                
                fig0.savefig(str(directorySave/('autocor_long_P'+str(part)+'.png')))
                fig0.savefig(str(directorySave/('autocor_long_P'+str(part)+'.svg')),format='svg',dpi=1200)
                plt.close()
        
        #Save data
        for p in range(self.nbParticles):
            autoCor = self.validParticles[p].autoCor
            timeD = self.validParticles[p].timeD
            if len(autoCor) > 0:
                with open(str(directorySave/('autocor_P'+str(p)+'.csv')), 'w') as textfile:
                    textfile.write('File,'+str(self.validParticles[p].fileName)+'\n')
                    textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                    textfile.write('TimeD (s),')
                    for t in timeD:
                        textfile.write("%.2f," % (t))
                    textfile.write('\n')
                    textfile.write('Angular autocorrelation,')
                    for m in autoCor:
                        textfile.write(("%.6f," % (m)))

    def doAverageMSAD(self):
        '''Plots the average MSAD with error bars in long and short format'''

        #Creates the folder MSAD if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/MSAD'))):
            os.mkdir(str(self.dn / Path('Plots/MSAD')))
        
        directorySave = self.dn / Path('Plots/MSAD')
                
        #Plot average MSAD, long version
        fig0 = plt.figure(0,figsize=(12, 10))
        
        plt.errorbar(self.averageTimeDMSAD,self.averageMSAD,yerr=self.semMSAD)
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSAD (rad$^2$)')
        plt.axis('tight')
        plt.title('Average MSAD (long)')
        
        fig0.savefig(str(directorySave/'MSAD_average_long.png'))
        fig0.savefig(str(directorySave/'MSAD_average_long.svg'),format='svg',dpi=1200)
        plt.close()
        
        #Plot average MSAD, short version
        fig0 = plt.figure(0,figsize=(12, 10))
        
        plt.errorbar(self.averageTimeDMSAD[:self.nbPoints],self.averageMSAD[:self.nbPoints],yerr=self.semMSAD[:self.nbPoints])
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSAD (rad$^2$)')
        plt.axis('tight')
        plt.title('Average MSAD (short)')
        
        fig0.savefig(str(directorySave/'MSAD_average_short.png'))
        fig0.savefig(str(directorySave/'MSAD_average_short.svg'),format='svg',dpi=1200)
        plt.close()  
        
        with open(str(directorySave/'MSAD_average.csv'), 'w') as textfile:
            textfile.write('TimeD (s),')
            for t in self.averageTimeDMSAD:
                textfile.write("%.2f," % (t))
            textfile.write('\n')
            textfile.write('MSAD (rad^2),')
            for m in self.averageMSAD:
                textfile.write(("%.6f," % (m)))    
            textfile.write('\n')
            textfile.write('STD (rad^2),')
            for m in self.stdMSAD:
                textfile.write(("%.6f," % (m)))  
            textfile.write('\n')
            textfile.write('SEM (rad^2),')
            for m in self.semMSAD:
                textfile.write(("%.6f," % (m)))  
            textfile.write('\n')
            textfile.write('n,')
            for m in self.nMSAD:
                textfile.write(("%s," % (str(m))))     







    def computeMSAD(self):
        '''The MSAD from Albert's tracking code is wrong.
        Even if the data is in the Excel file, this function calculates the 
        MSAD again. For that, it needs angle (continuous!!) and times data'''
        
        for part in range(self.nbParticles):
            self.validParticles[part].calculateMSAD()
            # self.validParticles[part].MSAD = self.mean_square(self.validParticles[part].angleExtended)
        


    def doMSAD(self):
        
        if not os.path.exists(str(self.dn / Path('Plots/MSAD'))):
            os.mkdir(str(self.dn / Path('Plots/MSAD')))
        
        directorySave = self.dn / Path('Plots/MSAD')
                      
        #Plot all MSAD (long version)
        fig0 = plt.figure(0,figsize=(12, 10))
        for part in range(self.nbParticles):
        
            timeD = self.validParticles[part].timeD
            MSAD = self.validParticles[part].MSAD
            
            if len(MSAD) > 0:
                plt.plot(timeD,MSAD)
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSAD (rad$^2$)')
        plt.axis('tight')
        plt.title('MSAD for all particles (long)')
        
        fig0.savefig(str(directorySave/'MSAD_long_all.png'))
        fig0.savefig(str(directorySave/'MSAD_long_all.svg'),format='svg',dpi=1200)
        plt.close()


        #Plot all MSAD (short version, 1/10 of data points)
        fig0 = plt.figure(0,figsize=(12, 10))
        for part in range(self.nbParticles):
        
            timeD = self.validParticles[part].timeD
            MSAD = self.validParticles[part].MSAD
            
            if len(MSAD) > 0:
                plt.plot(timeD[:self.nbPoints],MSAD[:self.nbPoints])
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSAD (rad$^2$)')
        plt.axis('tight')
        plt.title('MSAD for all particles (short)')
        
        fig0.savefig(str(directorySave/'MSAD_short_all.png'))
        fig0.savefig(str(directorySave/'MSAD_short_all.svg'),format='svg',dpi=1200)
        plt.close()
        
        
        #Creates the folder for individual particles if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/MSAD/Individual'))):
            os.mkdir(str(self.dn / Path('Plots/MSAD/Individual')))
        
        directorySave = self.dn / Path('Plots/MSAD/Individual')
        
        #Plots all MSAD one by one (short version)
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSAD = self.validParticles[part].MSAD
            
            if len(MSAD) > 0:
                plt.plot(timeD[:self.nbPoints],MSAD[:self.nbPoints])
            
                plt.xlabel('$\Delta$t (s)')
                plt.ylabel('MSAD (rad$^2$)')
                plt.axis('tight')
                plt.title('MSAD for particle ' + str(part))
                
                fig0.savefig(str(directorySave/('MSAD_short_P'+str(part)+'.png')))
                fig0.savefig(str(directorySave/('MSAD_short_P'+str(part)+'.svg')),format='svg',dpi=1200)
                plt.close()
                

        #Plots all MSAD one by one (long version)
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSAD = self.validParticles[part].MSAD
            
            if len(MSAD) > 0:
                plt.plot(timeD,MSAD)
                
                plt.xlabel('$\Delta$t (s)')
                plt.ylabel('MSAD (rad$^2$)')
                plt.axis('tight')
                plt.title('Autocorrelation for particle ' + str(part))
                
                fig0.savefig(str(directorySave/('MSAD_long_P'+str(part)+'.png')))
                fig0.savefig(str(directorySave/('MSAD_long_P'+str(part)+'.svg')),format='svg',dpi=1200)
                plt.close()
        
        #Save data
        for p in range(self.nbParticles):
            MSAD = self.validParticles[p].MSAD
            timeD = self.validParticles[p].timeD
            if len(MSAD) > 0:
                with open(str(directorySave/('MSAD_P'+str(p)+'.csv')), 'w') as textfile:
                    textfile.write('File,'+str(self.validParticles[p].fileName)+'\n')
                    textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                    textfile.write('TimeD (s),')
                    for t in timeD:
                        textfile.write("%.2f," % (t))
                    textfile.write('\n')
                    textfile.write('MSAD (rad^2),')
                    for m in MSAD:
                        textfile.write(("%.6f," % (m)))
            
            
    def doTrajectory(self):
        
        if not os.path.exists(str(self.dn / Path('Plots/Trajectory'))):
            os.mkdir(str(self.dn / Path('Plots/Trajectory')))
        
        directorySave = self.dn / Path('Plots/Trajectory')
        
        Xmax,Ymax,Xmin,Ymin = [list() for i in range(4)]
        
        #Normalize all trajectories
        for part in range(self.nbParticles):
                    
            X = self.validParticles[part].X
            Y = self.validParticles[part].Y
            
            #Normalize to (0,0)
            self.validParticles[part].X = [x - X[0] for x in X]
            self.validParticles[part].Y = [y - Y[0] for y in Y]
            
            Xmax.append(max(self.validParticles[part].X))
            Xmin.append(min(self.validParticles[part].X))
            Ymax.append(max(self.validParticles[part].Y))
            Ymin.append(min(self.validParticles[part].Y))
            
        Xmax = max(Xmax)
        Xmin = max(Xmin)
        Ymax = max(Ymax)
        Ymin = max(Ymin)
        
        totalMax = max([np.abs(Xmax),np.abs(Xmin),np.abs(Ymax),np.abs(Ymin)])
        
        lims = (-totalMax,totalMax)
        
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
            
            X = self.validParticles[part].X
            Y = self.validParticles[part].Y
            
            plt.plot(X,Y)
            
#            '''The next part is done to equalize the limits of x and y'''
#            x1,x2 = plt.xlim()
#            diffx = x2 - x1
#            y1,y2 = plt.ylim()
#            diffy = y2- y1
#            diff = diffx-diffy
#            if diff > 0:
#                plt.ylim(y1-(diff/2),y2+(diff/2))
#            elif diff < 0:
#                plt.xlim(x1-(-diff/2),x2+(-diff/2))
#            '''Done'''

            plt.xlim(lims)
            plt.ylim(lims)
            
            plt.xlabel('X position ($\mu$m)')
            plt.ylabel('Y position ($\mu$m)')
#            plt.axis('tight')
            plt.title('Trajectory for particle ' + str(part))
            
            fig0.savefig(str(directorySave/('trajectory_P'+str(part)+'.png')))
            fig0.savefig(str(directorySave/('trajectory_P'+str(part)+'.svg')),format='svg',dpi=1200)
            plt.close()
            
        for p in range(self.nbParticles):
            with open(str(directorySave/('trajectory_P'+str(p)+'.csv')), 'w') as textfile:
                textfile.write('File,'+str(self.validParticles[p].fileName)+'\n')
                textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                textfile.write('TimeD (s),')
                timeD = self.validParticles[p].timeD
                for t in timeD:
                    textfile.write("%.2f," % (t))
                textfile.write('\n')
                textfile.write('X (um),')
                X = self.validParticles[p].X
                for x in X:
                    textfile.write(("%.6f," % (x)))            
                textfile.write('\n')
                textfile.write('Y (um),')
                Y = self.validParticles[p].Y
                for y in Y:
                    textfile.write(("%.6f," % (y)))       

    def doAverageMSD(self):
        '''Plots the average MSD with error bars in long and short format'''
        
        #Creates the folder MSD if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/MSD'))):
            os.mkdir(str(self.dn / Path('Plots/MSD')))
        
        directorySave = self.dn / Path('Plots/MSD') 
                
        #Plot average MSD, long version
        fig0 = plt.figure(0,figsize=(12, 10))
        
        plt.errorbar(self.averageTimeD,self.averageMSD,yerr=self.semMSD)
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$)')
        plt.axis('tight')
        plt.title('Average MSD (long)')
        
        fig0.savefig(str(directorySave/'MSD_average_long.png'))
        fig0.savefig(str(directorySave/'MSD_average_long.svg'),format='svg',dpi=1200)
        plt.close()
        
        #Plot average MSD, short version
        fig0 = plt.figure(0,figsize=(12, 10))
                
        plt.errorbar(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints],yerr=self.semMSD[:self.nbPoints])
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$)')
        plt.axis('tight')
        plt.title('Average MSD (short)')
        
        fig0.savefig(str(directorySave/'MSD_average_short.png'))
        fig0.savefig(str(directorySave/'MSD_average_short.svg'),format='svg',dpi=1200)
        plt.close()  
        
        with open(str(directorySave/'MSD_average.csv'), 'w') as textfile:
            textfile.write('TimeD (s),')
            for t in self.averageTimeD:
                textfile.write("%.2f," % (t))
            textfile.write('\n')
            textfile.write('MSD (um^2/s),')
            for m in self.averageMSD:
                textfile.write(("%.6f," % (m))) 
            textfile.write('\n')
            textfile.write('STD (um^2/s),')
            for m in self.stdMSD:
                textfile.write(("%.6f," % (m)))  
            textfile.write('\n')
            textfile.write('SEM (um^2/s),')
            for m in self.semMSD:
                textfile.write(("%.6f," % (m)))  
            textfile.write('\n')
            textfile.write('n,')
            for m in self.nMSD:
                textfile.write(("%s," % (str(m))))     
            
    def doMSD(self):
        
        #Creates the folder MSD if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/MSD'))):
            os.mkdir(str(self.dn / Path('Plots/MSD')))
        
        directorySave = self.dn / Path('Plots/MSD')
        
        #Plot all MSDs (short version, 1/10 of data points)
        fig0 = plt.figure(0,figsize=(12, 10))
        for part in range(self.nbParticles):
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            plt.plot(timeD[:self.nbPoints],MSD[:self.nbPoints])
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$)')
        plt.axis('tight')
        plt.title('Short MSD for all particles')
        
        fig0.savefig(str(directorySave/'MSD_short_all.png'))
        fig0.savefig(str(directorySave/'MSD_short_all.svg'),format='svg',dpi=1200)
        plt.close()

        #Plots all MSDs (all points)
        fig0 = plt.figure(0,figsize=(12, 10))
        for part in range(self.nbParticles):
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            plt.plot(timeD,MSD)
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$)')
        plt.axis('tight')
        plt.title('Long MSD for all particles')
        
        fig0.savefig(str(directorySave/'MSD_long_all.png'))
        fig0.savefig(str(directorySave/'MSD_long_all.svg'),format='svg',dpi=1200)
        plt.close()

        #Creates the folder for individual particles if it's the first time
        if not os.path.exists(str(self.dn / Path('Plots/MSD/Individual'))):
            os.mkdir(str(self.dn / Path('Plots/MSD/Individual')))
        
        directorySave = self.dn / Path('Plots/MSD/Individual')
        
        #Plots all MSDs one by one (short version)
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            plt.plot(timeD[:self.nbPoints],MSD[:self.nbPoints])
            
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSD ($\mu$m$^2$)')
            plt.axis('tight')
            plt.title('Short MSD for particle ' + str(part))
            
            fig0.savefig(str(directorySave/('MSD_short_P'+str(part)+'.png')))
            fig0.savefig(str(directorySave/('MSD_short_P'+str(part)+'.svg')),format='svg',dpi=1200)
            plt.close()

        #Plots all MSDs one by one (long version)
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            plt.plot(timeD,MSD)
            
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSD ($\mu$m$^2$)')
            plt.axis('tight')
            plt.title('Long MSD for particle ' + str(part))
            
            fig0.savefig(str(directorySave/('MSD_long_P'+str(part)+'.png')))
            fig0.savefig(str(directorySave/('MSD_long_P'+str(part)+'.svg')),format='svg',dpi=1200)
            plt.close()
        
        for p in range(self.nbParticles):
            with open(str(directorySave/('MSD_P'+str(p)+'.csv')), 'w') as textfile:
                textfile.write('File,'+str(self.validParticles[p].fileName)+'\n')
                textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                textfile.write('TimeD (s),')
                timeD = self.validParticles[p].timeD
                for t in timeD:
                    textfile.write("%.2f," % (t))
                textfile.write('\n')
                textfile.write('MSD (um^2/s),')
                MSD = self.validParticles[p].MSD
                for m in MSD:
                    textfile.write(("%.6f," % (m)))
        

    def computeMSD(self):

        for part in range(self.nbParticles):
            self.validParticles[part].calculateMSD()

if __name__ == '__main__':
    root = tkinter.Tk()
    gui = GUI(root)
    root.title("Motion Visualization Tool v0.6")
    w = 350
    h = 800
    x = 200
    y = 150
    root.geometry('%dx%d+%d+%d' % (w, h, x, y))
#    root.geometry("300x280+200+200")
    root.mainloop()



