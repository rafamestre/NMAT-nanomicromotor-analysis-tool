"""
Motion Visualiation Tool

v0.2

15/04/2019

Changes from past version:
    
    -) Fixed issues reading .csv files:
        -) Happening to Shidong: If the files were opened and then saved, 
        Excel adds extra empty spaces in the empty cells. The reading was creating
        mistakes.
        -) If there is a row without any columns, it was giving error. A check-point
        to check for columns has been added.
    -) Added automatic installation of packages "future" and "seaborn"
    -) The std and sem of the MSD are computed. SEM is added as error in the
    average MSD plot. STD, SEM and n are saved in the .csv file.
    -) Corrected error of averaging MSAD when no MSAD is in the files.
    -) Trajectories for all particles have the same axes

@author: Rafael Mestre

This code was written to be compatible for both Python 2.6+ and 3.3+


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
import seaborn as sns

sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 3})
sns.set_style("ticks")
sns.set_palette(sns.color_palette("hls", 8))

def linear(x,a,b):
    return a + b*x

def quadratic(x,a,b):
    return a*x + b*x**2

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

plt.close("all")

class Particle:
    def __init__(self):
        
        self.particleLabel = 0
        self.time = list()
        self.X = list()
        self.Y = list()
        self.MSD = list()
        self.time = list()
        self.v = list()
        self.timeD = list()
        self.MSAD = list()
        self.autoCor = list()
        self.autoCorstd = list()
        self.fileName = ''
        
        self.alpha = list()
        self.alphaFitting = 0
        self.diffusionAlphaFitting = 0
        self.speedQuadraticFitting = 0
        self.speedSquaredQuadraticFitting = 0
        self.diffusionQuadraticFitting = 0
        self.diffusionLinearFitting = 0

        
        self.valid = False

    def calculateLocalAlpha(self):
        
        if self.valid:
            self.alpha = list()
            for i in range(len(self.MSD)-1):
                num = np.log(self.MSD[i+1])-np.log(self.MSD[i])
                den = np.log(self.timeD[i+1])-np.log(self.timeD[i])
                self.alpha.append(num/den)


class GUI:
    def __init__(self, master):
        
        self.master = master
        
        self.data = False
        
        self.calculateMSD = tkinter.IntVar()
        self.calculateTrajectory = tkinter.IntVar()
        self.instVelocity = tkinter.IntVar()
        self.calculateMSAD = tkinter.IntVar()
        self.autocorrelation = tkinter.IntVar()
        
        self.autoCorFit = tkinter.IntVar()

        self.localExponent = tkinter.IntVar() #Calculate the local exponent
        self.MSDFitting = tkinter.IntVar() #Fit the MSD to a polynomial
        self.fitAlpha = tkinter.IntVar() #Perform logarithmic fitting
        
        self.fittingInterval = tkinter.StringVar(root, value="")
        self.nbPoints = -1
        
        self.calculateAverage = tkinter.IntVar() #Make averages over particles

        self.doAll = tkinter.IntVar()
        
        self.particles = list()
        
        self.dn = os.path.dirname(os.path.realpath(__file__))+'\\'
        

        self.folderPath = tkinter.StringVar(root, value=self.dn)
        self.pathEntry = tkinter.Entry(master,width=28,textvariable=self.folderPath)
        self.browseButton = tkinter.Button(master,text="Browse",command=self.browseFolder)

        self.firstText = tkinter.StringVar(master, value ="")
        
        self.text1 = tkinter.Label(master, text = self.firstText.get())

        self.text2 = tkinter.Label(master,text = "What do you want to analyze?")

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
        self.fittingCombo = ttk.Combobox(master, values = ["Quadratic", "Linear","Both"],width=14)
        self.alphaButton = tkinter.Checkbutton(master,text="Perform logarithmic fitting?",
                                               variable=self.fitAlpha, 
                                               onvalue = True, offvalue = False)
        
        self.AverageButton = tkinter.Checkbutton(master, text="Do averages between all particles?", 
                                                    variable=self.calculateAverage, 
                                                    onvalue = True, offvalue = False)
        
        
        
        self.fitUpTo = tkinter.Label(master, text="Interval of interest:")
        self.fitIntervalEntry = tkinter.Entry(master,width=10,textvariable = self.fittingInterval)
        self.textSeconds = tkinter.StringVar(root, value="s out of ? s")
        self.seconds = tkinter.Label(master, text = self.textSeconds.get())
        
        self.analysisButton = tkinter.Button(master,text="Analyze",
                                             command=self.clickAnalysis)
        self.analysisButton.config(state="disabled")
        
        self.doAllButton = tkinter.Checkbutton(master,text="Do all",
                                               variable=self.doAll, 
                                               onvalue = True, offvalue = False,
                                               command = self.doEverything)
        self.doAllButton.config(state="disabled")
    
        self.pathEntry.grid(row=0,column=0,padx=25,pady=20,columnspan=2,sticky = 'W')
        self.browseButton.grid(row=0,column=2,pady=20,columnspan=1,sticky = 'W')
        
        self.text1.grid(row=1,column=0,ipadx = 25, ipady = 25, columnspan=3)
        self.text2.grid(row=2,column=0,ipady = 5, columnspan=3)
        
        self.MSDButton.grid(row=3,column=0,sticky = 'W', ipadx = 10, columnspan=3)
        self.trajectoryButton.grid(row=4,column=0,sticky = 'W', ipadx = 10, columnspan=3)
        self.vButton.grid(row=5,column=0,sticky = 'W', ipadx = 10, columnspan=3)   
        self.MSADButton.grid(row=6,column=0,sticky = 'W', ipadx = 10, columnspan=3)   
        self.autoCorrButton.grid(row=7,column=0,sticky = 'W', ipadx = 10, columnspan=3)   
        
        self.MSDOptions.grid(row=8,column=0,ipady = 15, columnspan=3)
        
        self.exponentButton.grid(row=9,column=0,sticky = 'W',ipadx = 10, columnspan=3)
        self.fittingButton.grid(row=10,column=0,sticky = 'W', ipadx = 10, columnspan=3)
        self.alphaButton.grid(row=11,column=0,sticky = 'W', ipadx = 10, columnspan=3)
    
        self.fitUpTo.grid(row=12,column=0,sticky = 'W',padx = 10,pady = 8)
        self.fitIntervalEntry.grid(row=12,column=1,sticky = 'W',pady = 8, columnspan=2)
        self.seconds.grid(row=12,column=1,sticky = 'W',pady = 8,padx = 50, columnspan=2)
        
        self.AverageButton.grid(row=13,column=0,sticky = 'W', ipadx = 10,ipady = 12, columnspan=3)
        self.analysisButton.grid(row=14,column=0, columnspan=3)

        self.doAllButton.grid(row=15,column=0,columnspan=3)

        #Initialize in the current folder
        self.updateValuesFolder()
        
        
#
    def browseFolder(self):
        self.dn = filedialog.askdirectory(initialdir=self.dn)
        self.dn += "\\"
        
        self.updateValuesFolder()

    def updateValuesFolder(self):
        
        self.readFiles()
        
        #Get valid particles
        self.validParticles= [p for p in self.particles if p.valid == True]
        self.nbParticles = len(self.validParticles)
        if self.nbParticles > 0:
            self.firstText.set("A total number of %d files were found\nwith a total of %d valid particles" % (self.nbFiles,self.nbParticles))
            #Find particle with longest time
            totalTime = list()
            for part in range(self.nbParticles):
                totalTime.append(self.validParticles[part].timeD[-1])
            self.longestTime = max(totalTime)
            self.textSeconds.set("s out of %.2f s" % (self.longestTime))
            self.data = True
            self.analysisButton.config(state="normal")
            self.fittingInterval.set("%.2f" % (self.longestTime/10))
            self.doAllButton.config(state="normal")
        else:
            self.firstText.set("No data was found in that folder")
            self.textSeconds.set("s out of ? s")
            self.data = False
            self.analysisButton.config(state="disabled")
            self.doAllButton.config(state="disabled")
            self.fittingInterval.set("")
        self.seconds.config(text=self.textSeconds.get())
        self.text1.config(text=self.firstText.get())
        self.folderPath.set(self.dn[:-1])
        self.pathEntry.config(text=self.folderPath)  
        self.fitIntervalEntry.config(text=self.fittingInterval)

    def doEverything(self):

        if self.doAll.get():
            self.MSDButton.select()
            self.trajectoryButton.select()
            self.vButton.select()  
            self.MSADButton.select() 
            self.autoCorrButton.select()  
                        
            self.exponentButton.select()
            self.fittingButton.select()
            self.alphaButton.select()
                    
            self.AverageButton.select()
            
            self.fittingCombo.set("Both")
            self.fittingChecked()
        elif not self.doAll.get():
            self.MSDButton.deselect()
            self.trajectoryButton.deselect()
            self.vButton.deselect()  
            self.MSADButton.deselect() 
            self.autoCorrButton.deselect()  
                        
            self.exponentButton.deselect()
            self.fittingButton.deselect()
            self.alphaButton.deselect()
                    
            self.AverageButton.deselect()
            
            self.fittingChecked()
                
        
    def fittingChecked(self):
        if self.MSDFitting.get():
            lv_x = self.fittingButton.winfo_x()
            lv_y = self.fittingButton.winfo_y()
            self.fittingCombo.place(x=lv_x+120,y=lv_y+2)
        else:
            self.fittingCombo.place_forget()


    def autoCorChecked(self):
        if self.autocorrelation.get():
            lv_x = self.autoCorrButton.winfo_x()
            lv_y = self.autoCorrButton.winfo_y()
#            self.autoCorFittingButton.place(x=lv_x+220,y=lv_y)
        else:
            self.autoCorFittingButton.place_forget()

    def readFiles(self):
        
        self.particles = list()
        self.files = [i for i in glob.glob(self.dn+'*.{}'.format('csv'))]
        self.nbFiles = len(self.files)
        
        for f in range(self.nbFiles):
                             
            fileName = self.files[f]
            
#            directorySave = directory+'\\Plots\\'
            
            self.firstParticle = 0
            
            self.readData(fileName)
            
        for p in range(len(self.particles)):
            self.checkParticleValidity(self.particles[p])

    def readData(self,f):
        
        summaryParticle = True
        
        with open(f,'r') as csvfile:
            fileReader = csv.reader(csvfile, delimiter=',')
            for row in fileReader:
                if len(row) == 0:
                    continue
                if row[0] == 'Particle':
                    if not summaryParticle:
                        p = Particle()
                        p.particleLabel = row[1]
                        p.fileName = f
                        self.particles.append(p)
                    if summaryParticle:
                        summaryParticle = False
                        '''The first tile the word "Particle" appears is in the
                        summary table. Therefore, when Particle appears the first
                        time, it's ignored and summaryParticle is turned to 1'''
                if row[0] == 'Time (seconds)':
                    self.particles[-1].time = [float(i) for i in row[1:] if i != '']
                if row[0] == 'X (microm)':
                    self.particles[-1].X = [float(i) for i in row[1:] if i != '']       
                if row[0] == 'Y (microm)':
                    self.particles[-1].Y = [float(i) for i in row[1:] if i != '']            
                if row[0] == 'Instantaneous velocity (microm/seconds)':
                    self.particles[-1].v = [float(i) for i in row[1:] if i != '']
#                    p.velocityLongList.append([float(i) for i in row[1:]])
                if row[0] == 'Time displacement (seconds)':
                    self.particles[-1].timeD = [float(i) for i in row[1:] if i != '']  
                if row[0] == 'MSD (microm**2)':
                    self.particles[-1].MSD = [float(i) for i in row[1:] if i != '']
                if row[0] == "Mean square angular displacement from vel vector (degrees**2)":
                    self.particles[-1].MSAD = [float(i) for i in row[1:] if i != '']
                if row[0] == "Angular auto-correlation from vel vector ":
                    self.particles[-1].autoCor = [float(i) for i in row[1:] if i != '']
                if row[0] == "Angular auto-correlation from vel vector standard deviation":
                    self.particles[-1].autoCorstd = [float(i) for i in row[1:] if i != '']
    
    def checkParticleValidity(self,p):
        
        if (not p.X) or (not p.Y) or (not p.MSD) or (not p.timeD) or (not p.time):
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

        if self.data:
            self.doDataAnalysis()

    def doDataAnalysis(self):
        
        if not os.path.exists(self.dn+'\\Plots\\'):
            os.mkdir(self.dn+'\\Plots\\')
                
        
        if self.calculateAverage.get():
            self.computeAverages()
        
        if self.calculateMSD.get():
            self.doMSD()
            if self.calculateAverage.get():
                self.doAverageMSD()
        
        if self.calculateTrajectory.get():  
            self.doTrajectory()
            
        if self.instVelocity.get():
            self.doInstVel()
            
        if self.calculateMSAD.get():
            self.doMSAD()
            if self.calculateAverage.get():
                self.doAverageMSAD()
            
        if self.autocorrelation.get():
            self.doAutoCor()
            if self.calculateAverage.get():
                self.doAverageAutoCor()
            
        if self.MSDFitting.get() and self.fittingCombo.get() == "Quadratic":
            self.doQuadraticFitting()
            if self.calculateAverage.get():
                self.doQuadraticFittingAverage()
            
        if self.MSDFitting.get() and self.fittingCombo.get() == "Linear":
            self.doLinearFitting()
            if self.calculateAverage.get():
                self.doLinearFittingAverage()
                
        if self.MSDFitting.get() and self.fittingCombo.get() == "Both":
            self.doLinearFitting()
            self.doQuadraticFitting()
            if self.calculateAverage.get():
                self.doLinearFittingAverage()
                self.doQuadraticFittingAverage()
            
        if self.localExponent.get():
            self.doLocalAlpha()
            if self.calculateAverage.get():
                self.doAverageAlpha()
                
        if self.fitAlpha.get():
            self.doAlphaFitting()
            if self.calculateAverage.get():
                self.doAlphaFittingAverage()        

        print("Analysis finished")

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
        
        
    def doLinearFittingAverage(self):
        #The doAverageMSD function should have been run first
        #so that the averageMSD variable exists
        
        #Creates the folder if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\Fittings\\'):
            os.mkdir(self.dn+'\\Plots\\Fittings\\')
        
        directorySave = self.dn+'\\Plots\\Fittings\\'          

        #Plot MSD  fitting info
        fig0 = plt.figure(0,figsize=(12, 10))
                
        popt, pcov = curve_fit(linearZero, self.averageTimeD[:self.nbPoints], self.averageMSD[:self.nbPoints])  
        
        plt.plot(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints])
        
        label = "Fitting equation: MSD = $4D_{t}t$\n$D_t$ = %.4f $\mu$m$^2$/s" % (popt[0]/4)
        plt.plot(self.averageTimeD[:self.nbPoints],linearZero(np.asarray(self.averageTimeD[:self.nbPoints]),*popt),'--',
                 label = label)
    
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$/s)')
        plt.axis('tight')
        plt.legend()
        plt.title('MSD linear fitting for all particles')
        
        fig0.savefig(directorySave+'MSD_linearFitting_average.png')
        fig0.savefig(directorySave+'MSD_linearFitting_average.svg',format='svg',dpi=1200)
        plt.close()
    
        #Save data
        with open(directorySave+'MSD_linearFitting_average.csv', 'w') as textfile:
            textfile.write("Fitting equation,MSD = 4D*t\n")
            textfile.write("D (um^2/s),%.6f\n" % (popt[0]/4))


    def doLinearFitting(self):
        
        
        #Creates the folder fittings if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\Fittings'):
            os.mkdir(self.dn+'\\Plots\\Fittings')
            
        #Creates the folder for alpha fitting if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\Fittings\\Linear\\'):
            os.mkdir(self.dn+'\\Plots\\Fittings\\Linear\\')
        
        directorySave = self.dn+'\\Plots\\Fittings\\Linear\\\\'  
        
        #Plot MSD with fitting info
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            popt, pcov = curve_fit(linearZero, timeD[:self.nbPoints], MSD[:self.nbPoints])  
            
            self.validParticles[part].diffusionLinearFitting = popt[0]/4
            
            plt.plot(timeD[:self.nbPoints],MSD[:self.nbPoints])
            
            label = "Fitting equation: MSD = $4D_{t}t$\n$D_t$ = %.4f $\mu$m$^2$/s" % (popt[0]/4)
            plt.plot(timeD[:self.nbPoints],linearZero(np.asarray(timeD[:self.nbPoints]),*popt),'--',
                     label = label)
        
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSD ($\mu$m$^2$/s)')
            plt.axis('tight')
            plt.legend()
            plt.title('MSD linear fitting for particle ' + str(part))
            
            fig0.savefig(directorySave+'MSD_linearFitting_P'+str(part)+'.png')
            fig0.savefig(directorySave+'MSD_linearFitting_P'+str(part)+'.svg',format='svg',dpi=1200)
            plt.close()

        #Save data
        for p in range(self.nbParticles):
            with open(directorySave+'MSD_linearFitting_P'+str(p)+'.csv', 'w') as textfile:
                textfile.write('File,'+self.validParticles[p].fileName+'\n')
                textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                textfile.write("Fitting equation,MSD = 4D*t\n")
                textfile.write("D (um^2/s),%.6f\n" % (self.validParticles[p].diffusionLinearFitting))
        
        #Save summary data
        with open(directorySave+'MSD_linearFitting_Summary.csv', 'w') as textfile:
            Dlist = list()
            textfile.write("Fitting equation,MSD = 4D*t\n")
            textfile.write('D (um^2/s),')
            for p in range(self.nbParticles):
                textfile.write("%.6f," % (self.validParticles[p].diffusionLinearFitting))
                Dlist.append(self.validParticles[p].diffusionLinearFitting)
            textfile.write('\n')
            textfile.write('Average diffusion,%.6f\n' % (np.mean(Dlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(Dlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(Dlist,ddof=1)/np.sqrt(len(Dlist))))



    def doQuadraticFittingAverage(self):
        #The doAverageMSD function should have been run first
        #so that the averageMSD variable exists
        
        #Creates the folder if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\Fittings\\'):
            os.mkdir(self.dn+'\\Plots\\Fittings\\')
        
        directorySave = self.dn+'\\Plots\\Fittings\\'          

        #Plot MSD  fitting info
        fig0 = plt.figure(0,figsize=(12, 10))
                
        popt, pcov = curve_fit(quadratic, self.averageTimeD[:self.nbPoints], self.averageMSD[:self.nbPoints])  
        
        plt.plot(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints])
        
        label = "Fitting equation: MSD = $4D_{t}t + v^2t^2$\n$D_t$ = %.4f $\mu$m$^2$/s\n$v$ = %.4f $\mu$m/s\n$v^2$ = %.4f $\mu$m$^2/s^2$" % (popt[0]/4,np.sqrt(popt[1]),popt[1])
        plt.plot(self.averageTimeD[:self.nbPoints],quadratic(np.asarray(self.averageTimeD[:self.nbPoints]),*popt),'--',
                 label = label)
    
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$/s)')
        plt.axis('tight')
        plt.legend()
        plt.title('MSD quadratic fitting for all particles')
        
        fig0.savefig(directorySave+'MSD_quadraticFitting_average.png')
        fig0.savefig(directorySave+'MSD_quadraticFitting_average.svg',format='svg',dpi=1200)
        plt.close()
    
        #Save data
        with open(directorySave+'MSD_quadraticFitting_average.csv', 'w') as textfile:
            textfile.write("Fitting equation,MSD = 4D*t + v^2*t^2\n")
            textfile.write("D (um^2/s),%.6f\n" % (popt[0]/4))
            textfile.write("v (um/s),%.6f\n" % (np.sqrt(popt[1])))
            textfile.write("v^2 (um^2/s^2),%.6f\n" % (popt[1]))
        

    def doQuadraticFitting(self):
        
        #Creates the folder fittings if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\Fittings'):
            os.mkdir(self.dn+'\\Plots\\Fittings')
            
        #Creates the folder for alpha fitting if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\Fittings\\Quadratic\\'):
            os.mkdir(self.dn+'\\Plots\\Fittings\\Quadratic\\')
        
        directorySave = self.dn+'\\Plots\\Fittings\\Quadratic\\\\'  
        
        #Plot MSD with fitting info
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            popt, pcov = curve_fit(quadratic, timeD[:self.nbPoints], MSD[:self.nbPoints])  
            
            self.validParticles[part].speedQuadraticFitting = np.sqrt(popt[1])
            self.validParticles[part].speedSquaredQuadraticFitting = popt[1]
            self.validParticles[part].diffusionQuadraticFitting = popt[0]/4
            
            plt.plot(timeD[:self.nbPoints],MSD[:self.nbPoints])
            
            label = "Fitting equation: MSD = $4D_{t}t + v^2t^2$\n$D_t$ = %.4f $\mu$m$^2$/s\n$v$ = %.4f $\mu$m/s\n$v^2$ = %.4f $\mu$m$^2/s^2$" % (popt[0]/4,np.sqrt(popt[1]),popt[1])
            plt.plot(timeD[:self.nbPoints],quadratic(np.asarray(timeD[:self.nbPoints]),*popt),'--',
                     label = label)
        
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSD ($\mu$m$^2$/s)')
            plt.axis('tight')
            plt.legend()
            plt.title('MSD quadratic fitting for particle ' + str(part))
            
            fig0.savefig(directorySave+'MSD_quadraticFitting_P'+str(part)+'.png')
            fig0.savefig(directorySave+'MSD_quadraticFitting_P'+str(part)+'.svg',format='svg',dpi=1200)
            plt.close()

        #Save data
        for p in range(self.nbParticles):
            with open(directorySave+'MSD_quadraticFitting_P'+str(p)+'.csv', 'w') as textfile:
                textfile.write('File,'+self.validParticles[p].fileName+'\n')
                textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                textfile.write("Fitting equation,MSD = 4D*t + v^2*t^2\n")
                textfile.write("D (um^2/s),%.6f\n" % (self.validParticles[p].diffusionQuadraticFitting))
                textfile.write("v (um/s),%.6f\n" % (self.validParticles[p].speedQuadraticFitting))
                textfile.write("v^2 (um^2/s^2),%.6f\n" % (self.validParticles[p].speedSquaredQuadraticFitting))

        #Save summary data
        with open(directorySave+'MSD_quadraticFitting_Summary.csv', 'w') as textfile:
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
            textfile.write('Average diffusion,%.6f\n' % (np.mean(Dlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(Dlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(Dlist,ddof=1)/np.sqrt(len(Dlist))))
            textfile.write('Average speed,%.6f\n' % (np.mean(vlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(vlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(vlist,ddof=1)/np.sqrt(len(vlist))))
            textfile.write('Average speed squared,%.6f\n' % (np.mean(vsquarelist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(vsquarelist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(vsquarelist,ddof=1)/np.sqrt(len(vsquarelist))))

    def doAlphaFittingAverage(self):
        #The doAverageMSD function should have been run first
        #so that the averageMSD variable exists
        
        #Creates the folder if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\Fittings\\'):
            os.mkdir(self.dn+'\\Plots\\Fittings\\')
        
        directorySave = self.dn+'\\Plots\\Fittings\\'          

        #Plot MSD in log-log with fitting info
        fig0 = plt.figure(0,figsize=(12, 10))
                
        popt, pcov = curve_fit(linear, np.log(self.averageTimeD[:self.nbPoints]), np.log(self.averageMSD[:self.nbPoints]))  
        
        plt.plot(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints])
        
        label = "Fitting equation: MSD = $Dt^{\\alpha}$\n$\\alpha$ = %.2f\n$D$ = %.6f" % (popt[1],np.exp(popt[0]))
        plt.plot(self.averageTimeD[:self.nbPoints],powerLaw(np.asarray(self.averageTimeD[:self.nbPoints]),np.exp(popt[0]),popt[1]),'--',
                 label = label)
    
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$/s)')
        plt.yscale('log')
        plt.xscale('log')
        plt.axis('tight')
        plt.legend()
        plt.title('MSD exponent for all particles')
        
        fig0.savefig(directorySave+'MSD_logFitting_average.png')
        fig0.savefig(directorySave+'MSD_logFitting_average.svg',format='svg',dpi=1200)
        plt.close()

        #Save data
        with open(directorySave+'MSD_logFitting_average.csv', 'w') as textfile:
            textfile.write("Fitting equation,MSD = D*t^alpha\n")
            textfile.write("D (um^2/s^alpha),%.6f\n" % (popt[1]))
            textfile.write("Alpha,%.6f\n" % (np.exp(popt[0])))        

    def doAlphaFitting(self):
        
        
        #Creates the folder fittings if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\Fittings'):
            os.mkdir(self.dn+'\\Plots\\Fittings')
            
        #Creates the folder for alpha fitting if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\Fittings\\Logarithmic\\'):
            os.mkdir(self.dn+'\\Plots\\Fittings\\Logarithmic\\')
        
        directorySave = self.dn+'\\Plots\\Fittings\\Logarithmic\\\\'  
        
        #Plot MSD in log-log with fitting info
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            popt, pcov = curve_fit(linear, np.log(timeD[:self.nbPoints]), np.log(MSD[:self.nbPoints]))  
            
            self.validParticles[part].alphaFitting = popt[1]
            self.validParticles[part].diffusionAlphaFitting = np.exp(popt[0])
            
            plt.plot(timeD[:self.nbPoints],MSD[:self.nbPoints])
            
            label = "Fitting equation: MSD = $Dt^{\\alpha}$\n$\\alpha$ = %.2f\n$D$ = %.6f" % (popt[1],np.exp(popt[0]))
            plt.plot(timeD[:self.nbPoints],powerLaw(np.asarray(timeD[:self.nbPoints]),np.exp(popt[0]),popt[1]),'--',
                     label = label)
        
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSD ($\mu$m$^2$/s)')
            plt.yscale('log')
            plt.xscale('log')
            plt.axis('tight')
            plt.legend()
            plt.title('MSD exponent for particle ' + str(part))
            
            fig0.savefig(directorySave+'MSD_logFitting_P'+str(part)+'.png')
            fig0.savefig(directorySave+'MSD_logFitting_P'+str(part)+'.svg',format='svg',dpi=1200)
            plt.close()

        #Save data
        for p in range(self.nbParticles):
            with open(directorySave+'MSD_logFitting_P'+str(p)+'.csv', 'w') as textfile:
                textfile.write('File,'+self.validParticles[p].fileName+'\n')
                textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                textfile.write("Fitting equation,MSD = D*t^alpha\n")
                textfile.write("D (um^2/s^alpha),%.6f\n" % (self.validParticles[p].diffusionAlphaFitting))
                textfile.write("Alpha,%.6f\n" % (self.validParticles[p].alphaFitting))

        #Save summary data
        with open(directorySave+'MSD_logFitting_Summary.csv', 'w') as textfile:
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
            textfile.write('Average diffusion,%.6f\n' % (np.mean(Dlist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(Dlist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(Dlist,ddof=1)/np.sqrt(len(Dlist))))
            textfile.write('Average alpha,%.6f\n' % (np.mean(alphalist)))
            textfile.write('Standard deviation,%.6f\n' % (np.std(alphalist,ddof=1)))
            textfile.write('Standard error of the mean,%.6f\n' % (np.std(alphalist,ddof=1)/np.sqrt(len(alphalist))))



    def doAverageAlpha(self):
        #The doAverageMSD function should have been run first
        #so that the averageMSD variable exists
        
        #Creates the folder alpha if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\Exponent'):
            os.mkdir(self.dn+'\\Plots\\Exponent')
        
        directorySave = self.dn+'\\Plots\\Exponent\\'  
                
        
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
        
        fig0.savefig(directorySave+'alpha_average_long.png')
        fig0.savefig(directorySave+'alpha_average_long.svg',format='svg',dpi=1200)
        plt.close()
        
        #Plot average alpha, short version
        fig0 = plt.figure(0,figsize=(12, 10))
                
        plt.plot(self.averageTimeD[:self.nbPoints],self.averageAlpha[:self.nbPoints])
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD exponent')
        plt.axis('tight')
        plt.title('Average MSD exponent (short)')
        
        fig0.savefig(directorySave+'alpha_average_short.png')
        fig0.savefig(directorySave+'alpha_average_short.svg',format='svg',dpi=1200)
        plt.close()  
        
        with open(directorySave+'alpha_average.csv', 'w') as textfile:
            textfile.write('TimeD (s),')
            for t in self.averageTimeD:
                textfile.write("%.2f," % (t))
            textfile.write('\n')
            textfile.write('MSD exponent,')
            for m in self.averageAlpha:
                textfile.write(("%.6f," % (m)))     
                
                
                
    def doLocalAlpha(self):
        
        if not os.path.exists(self.dn+'\\Plots\\Exponent'):
            os.mkdir(self.dn+'\\Plots\\Exponent')
        
        directorySave = self.dn+'\\Plots\\Exponent\\'


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
        
        fig0.savefig(directorySave+'alpha_long_all.png')
        fig0.savefig(directorySave+'alpha_long_all.svg',format='svg',dpi=1200)
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
        
        fig0.savefig(directorySave+'alpha_short_all.png')
        fig0.savefig(directorySave+'alpha_short_all.svg',format='svg',dpi=1200)
        plt.close()
        
        #Creates the folder for individual particles if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\Exponent\\Individual\\'):
            os.mkdir(self.dn+'\\Plots\\Exponent\\Individual\\')
        
        directorySave = self.dn+'\\Plots\\Exponent\\Individual\\'
        
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
            
            fig0.savefig(directorySave+'alpha_short_P'+str(part)+'.png')
            fig0.savefig(directorySave+'alpha_short_P'+str(part)+'.svg',format='svg',dpi=1200)
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
            
            fig0.savefig(directorySave+'alpha_long_P'+str(part)+'.png')
            fig0.savefig(directorySave+'alpha_long_P'+str(part)+'.svg',format='svg',dpi=1200)
            plt.close()
        
        #Save data
        for p in range(self.nbParticles):
            alpha = self.validParticles[p].alpha
            timeD = self.validParticles[p].timeD[:-1] #By default, alpha has one less value
            if len(alpha) > 0:
                with open(directorySave+'alpha_P'+str(p)+'.csv', 'w') as textfile:
                    textfile.write('File,'+self.validParticles[p].fileName+'\n')
                    textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                    textfile.write('TimeD (s),')
                    for t in timeD:
                        textfile.write("%.2f," % (t))
                    textfile.write('\n')
                    textfile.write('MSD exponent,')
                    for m in alpha:
                        textfile.write(("%.6f," % (m)))



    def doInstVel(self):
        
        if not os.path.exists(self.dn+'\\Plots\\InstantaneousVelocity'):
            os.mkdir(self.dn+'\\Plots\\InstantaneousVelocity')
        
        directorySave = self.dn+'\\Plots\\InstantaneousVelocity\\'
        
        #Plot all velocities
        fig0 = plt.figure(0,figsize=(12, 10))
        for part in range(self.nbParticles):
        
            timeD = self.validParticles[part].timeD
            v = self.validParticles[part].v[:-1] #The speed has one value less
            
            if len(v) > 0:
                plt.plot(timeD,v)
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('Instantaneous velocity ($\mu$m/s)')
        plt.axis('tight')
        plt.title('Instantaneous velocity for all particles')
        
        fig0.savefig(directorySave+'instVel_all.png')
        fig0.savefig(directorySave+'instVel_all.svg',format='svg',dpi=1200)
        plt.close()

        #Creates the folder for individual particles if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\InstantaneousVelocity\\Individual\\'):
            os.mkdir(self.dn+'\\Plots\\InstantaneousVelocity\\Individual\\')
        
        directorySave = self.dn+'\\Plots\\InstantaneousVelocity\\Individual\\'
        
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
                plt.ylabel('Velocity autocorrelation ($\mu$m$^2$/s$^2$)')
                plt.axis('tight')
                plt.title('Instantaneous velocity for particle ' + str(part))
                plt.legend()
                
                fig0.savefig(directorySave+'instVel_P'+str(part)+'.png')
                fig0.savefig(directorySave+'instVel_P'+str(part)+'.svg',format='svg',dpi=1200)
                plt.close()
        
        #Save data
        for p in range(self.nbParticles):
            v = self.validParticles[p].v[:-1]
            timeD = self.validParticles[p].timeD
            if len(v) > 0:
                averageV = np.mean(v)
                stdV = np.std(v)
                with open(directorySave+'instVel_P'+str(p)+'.csv', 'w') as textfile:
                    textfile.write('File,'+self.validParticles[p].fileName+'\n')
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
        if not os.path.exists(self.dn+'\\Plots\\Autocorrelation'):
            os.mkdir(self.dn+'\\Plots\\Autocorrelation')
        
        directorySave = self.dn+'\\Plots\\Autocorrelation\\'  
                
        #Plot average MSD, long version
        fig0 = plt.figure(0,figsize=(12, 10))
        
        plt.plot(self.averageTimeDAutoCor,self.averageAutoCor)
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('Velocity autocorrelation ($\mu$m$^2$/s$^2$)')
        plt.axis('tight')
        plt.title('Average autocorrelation function (long)')
        
        fig0.savefig(directorySave+'AutoCor_average_long.png')
        fig0.savefig(directorySave+'AutoCor_average_long.svg',format='svg',dpi=1200)
        plt.close()
        
        #Plot average MSD, short version
        fig0 = plt.figure(0,figsize=(12, 10))
                
        plt.plot(self.averageTimeDAutoCor[:self.nbPoints],self.averageAutoCor[:self.nbPoints])
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('Velocity autocorrelation ($\mu$m$^2$/s$^2$)')
        plt.axis('tight')
        plt.title('Average autocorrelation function (short)')
        
        fig0.savefig(directorySave+'AutoCor_average_short.png')
        fig0.savefig(directorySave+'AutoCor_average_short.svg',format='svg',dpi=1200)
        plt.close()  
        
        with open(directorySave+'AutoCor_average.csv', 'w') as textfile:
            textfile.write('TimeD (s),')
            for t in self.averageTimeDAutoCor:
                textfile.write("%.2f," % (t))
            textfile.write('\n')
            textfile.write('MSD (um^2/s),')
            for m in self.averageAutoCor:
                textfile.write(("%.6f," % (m)))        


    def doAutoCor(self):
                
        if not os.path.exists(self.dn+'\\Plots\\Autocorrelation'):
            os.mkdir(self.dn+'\\Plots\\Autocorrelation')
        
        directorySave = self.dn+'\\Plots\\Autocorrelation\\'
        
        #Plot all autocorrelations (long version)
        fig0 = plt.figure(0,figsize=(12, 10))
        for part in range(self.nbParticles):
        
            timeD = self.validParticles[part].timeD
            autoCor = self.validParticles[part].autoCor
            
            if len(autoCor) > 0:
                plt.plot(timeD,autoCor)
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('Velocity autocorrelation ($\mu$m$^2$/s$^2$)')
        plt.axis('tight')
        plt.title('Autocorrelation function for all particles (long)')
        
        fig0.savefig(directorySave+'autocor_long_all.png')
        fig0.savefig(directorySave+'autocor_long_all.svg',format='svg',dpi=1200)
        plt.close()


        #Plot all autocorrelations (short version, 1/10 of data points)
        fig0 = plt.figure(0,figsize=(12, 10))
        for part in range(self.nbParticles):
        
            timeD = self.validParticles[part].timeD
            autoCor = self.validParticles[part].autoCor
            
            if len(autoCor) > 0:
                plt.plot(timeD[:self.nbPoints],autoCor[:self.nbPoints])
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('Velocity autocorrelation ($\mu$m$^2$/s$^2$)')
        plt.axis('tight')
        plt.title('Autocorrelation function for all particles (short)')
        
        fig0.savefig(directorySave+'autocor_short_all.png')
        fig0.savefig(directorySave+'autocor_short_all.svg',format='svg',dpi=1200)
        plt.close()
        
        
        #Creates the folder for individual particles if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\Autocorrelation\\Individual\\'):
            os.mkdir(self.dn+'\\Plots\\Autocorrelation\\Individual\\')
        
        directorySave = self.dn+'\\Plots\\Autocorrelation\\Individual\\'
        
        #Plots all Autocorrelation one by one (short version)
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            autoCor = self.validParticles[part].autoCor
            
            if len(autoCor) > 0:
                plt.plot(timeD[:self.nbPoints],autoCor[:self.nbPoints])
            
                plt.xlabel('$\Delta$t (s)')
                plt.ylabel('Velocity autocorrelation ($\mu$m$^2$/s$^2$)')
                plt.axis('tight')
                plt.title('Autocorrelation for particle ' + str(part))
                
                fig0.savefig(directorySave+'autocor_short_P'+str(part)+'.png')
                fig0.savefig(directorySave+'autocor_short_P'+str(part)+'.svg',format='svg',dpi=1200)
                plt.close()
                

        #Plots all Autocorrelation one by one (long version)
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            autoCor = self.validParticles[part].autoCor
            
            if len(autoCor) > 0:
                plt.plot(timeD,autoCor)
                
                plt.xlabel('$\Delta$t (s)')
                plt.ylabel('Velocity autocorrelation ($\mu$m$^2$/s$^2$)')
                plt.axis('tight')
                plt.title('Autocorrelation for particle ' + str(part))
                
                fig0.savefig(directorySave+'autocor_long_P'+str(part)+'.png')
                fig0.savefig(directorySave+'autocor_long_P'+str(part)+'.svg',format='svg',dpi=1200)
                plt.close()
        
        #Save data
        for p in range(self.nbParticles):
            autoCor = self.validParticles[p].autoCor
            timeD = self.validParticles[p].timeD
            if len(autoCor) > 0:
                with open(directorySave+'autocor_P'+str(p)+'.csv', 'w') as textfile:
                    textfile.write('File,'+self.validParticles[p].fileName+'\n')
                    textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                    textfile.write('TimeD (s),')
                    for t in timeD:
                        textfile.write("%.2f," % (t))
                    textfile.write('\n')
                    textfile.write('Autocorrelation (um^2/s^2),')
                    for m in autoCor:
                        textfile.write(("%.6f," % (m)))

    def doAverageMSAD(self):
        
        #Creates the folder MSD if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\MSAD'):
            os.mkdir(self.dn+'\\Plots\\MSAD')
        
        directorySave = self.dn+'\\Plots\\MSAD\\'  
                
        #Plot average MSD, long version
        fig0 = plt.figure(0,figsize=(12, 10))
        
        plt.plot(self.averageTimeDMSAD,self.averageMSAD)
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSAD (degree$^2$)')
        plt.axis('tight')
        plt.title('Average MSAD (long)')
        
        fig0.savefig(directorySave+'MSAD_average_long.png')
        fig0.savefig(directorySave+'MSAD_average_long.svg',format='svg',dpi=1200)
        plt.close()
        
        #Plot average MSD, short version
        fig0 = plt.figure(0,figsize=(12, 10))
                
        plt.plot(self.averageTimeDMSAD[:self.nbPoints],self.averageMSAD[:self.nbPoints])
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSAD (degree$^2$)')
        plt.axis('tight')
        plt.title('Average MSAD (short)')
        
        fig0.savefig(directorySave+'MSAD_average_short.png')
        fig0.savefig(directorySave+'MSAD_average_short.svg',format='svg',dpi=1200)
        plt.close()  
        
        with open(directorySave+'MSAD_average.csv', 'w') as textfile:
            textfile.write('TimeD (s),')
            for t in self.averageTimeDMSAD:
                textfile.write("%.2f," % (t))
            textfile.write('\n')
            textfile.write('MSD (um^2/s),')
            for m in self.averageMSAD:
                textfile.write(("%.6f," % (m)))        

            
    def doMSAD(self):
                
        if not os.path.exists(self.dn+'\\Plots\\MSAD'):
            os.mkdir(self.dn+'\\Plots\\MSAD')
        
        directorySave = self.dn+'\\Plots\\MSAD\\'
        
        #Plot all MSAD (long version)
        fig0 = plt.figure(0,figsize=(12, 10))
        for part in range(self.nbParticles):
        
            timeD = self.validParticles[part].timeD
            MSAD = self.validParticles[part].MSAD
            
            if len(MSAD) > 0:
                plt.plot(timeD,MSAD)
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSAD (degree$^2$)')
        plt.axis('tight')
        plt.title('MSAD for all particles (long)')
        
        fig0.savefig(directorySave+'MSAD_long_all.png')
        fig0.savefig(directorySave+'MSAD_long_all.svg',format='svg',dpi=1200)
        plt.close()


        #Plot all MSAD (short version, 1/10 of data points)
        fig0 = plt.figure(0,figsize=(12, 10))
        for part in range(self.nbParticles):
        
            timeD = self.validParticles[part].timeD
            MSAD = self.validParticles[part].MSAD
            
            if len(MSAD) > 0:
                plt.plot(timeD[:self.nbPoints],MSAD[:self.nbPoints])
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSAD (degree$^2$)')
        plt.axis('tight')
        plt.title('MSAD for all particles (short)')
        
        fig0.savefig(directorySave+'MSAD_short_all.png')
        fig0.savefig(directorySave+'MSAD_short_all.svg',format='svg',dpi=1200)
        plt.close()
        
        
        #Creates the folder for individual particles if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\MSAD\\Individual\\'):
            os.mkdir(self.dn+'\\Plots\\MSAD\\Individual\\')
        
        directorySave = self.dn+'\\Plots\\MSAD\\Individual\\'
        
        #Plots all MSAD one by one (short version)
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSAD = self.validParticles[part].MSAD
            
            if len(MSAD) > 0:
                plt.plot(timeD[:self.nbPoints],MSAD[:self.nbPoints])
            
                plt.xlabel('$\Delta$t (s)')
                plt.ylabel('MSAD (degree$^2$)')
                plt.axis('tight')
                plt.title('MSAD for particle ' + str(part))
                
                fig0.savefig(directorySave+'MSAD_short_P'+str(part)+'.png')
                fig0.savefig(directorySave+'MSAD_short_P'+str(part)+'.svg',format='svg',dpi=1200)
                plt.close()
                

        #Plots all MSAD one by one (long version)
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSAD = self.validParticles[part].MSAD
            
            if len(MSAD) > 0:
                plt.plot(timeD,MSAD)
                
                plt.xlabel('$\Delta$t (s)')
                plt.ylabel('MSAD (degree$^2$)')
                plt.axis('tight')
                plt.title('Autocorrelation for particle ' + str(part))
                
                fig0.savefig(directorySave+'MSAD_long_P'+str(part)+'.png')
                fig0.savefig(directorySave+'MSAD_long_P'+str(part)+'.svg',format='svg',dpi=1200)
                plt.close()
        
        #Save data
        for p in range(self.nbParticles):
            MSAD = self.validParticles[p].MSAD
            timeD = self.validParticles[p].timeD
            if len(MSAD) > 0:
                with open(directorySave+'MSAD_P'+str(p)+'.csv', 'w') as textfile:
                    textfile.write('File,'+self.validParticles[p].fileName+'\n')
                    textfile.write('Particle label,'+str(self.validParticles[p].particleLabel)+'\n')
                    textfile.write('TimeD (s),')
                    for t in timeD:
                        textfile.write("%.2f," % (t))
                    textfile.write('\n')
                    textfile.write('MSAD (degree^2),')
                    for m in MSAD:
                        textfile.write(("%.6f," % (m)))
            
            
    def doTrajectory(self):
        
        if not os.path.exists(self.dn+'\\Plots\\Trajectory'):
            os.mkdir(self.dn+'\\Plots\\Trajectory')
        
        directorySave = self.dn+'\\Plots\\Trajectory\\'
        
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
            
            fig0.savefig(directorySave+'trajectory_P'+str(part)+'.png')
            fig0.savefig(directorySave+'trajectory_P'+str(part)+'.svg',format='svg',dpi=1200)
            plt.close()
            
        for p in range(self.nbParticles):
            with open(directorySave+'trajectory_P'+str(p)+'.csv', 'w') as textfile:
                textfile.write('File,'+self.validParticles[p].fileName+'\n')
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
        
        #Creates the folder MSD if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\MSD'):
            os.mkdir(self.dn+'\\Plots\\MSD')
        
        directorySave = self.dn+'\\Plots\\MSD\\'  
                
        #Plot average MSD, long version
        fig0 = plt.figure(0,figsize=(12, 10))
        
        plt.errorbar(self.averageTimeD,self.averageMSD,yerr=self.semMSD)
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$/s)')
        plt.axis('tight')
        plt.title('Average MSD (long)')
        
        fig0.savefig(directorySave+'MSD_average_long.png')
        fig0.savefig(directorySave+'MSD_average_long.svg',format='svg',dpi=1200)
        plt.close()
        
        #Plot average MSD, short version
        fig0 = plt.figure(0,figsize=(12, 10))
                
        plt.errorbar(self.averageTimeD[:self.nbPoints],self.averageMSD[:self.nbPoints],yerr=self.semMSD[:self.nbPoints])
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$/s)')
        plt.axis('tight')
        plt.title('Average MSD (short)')
        
        fig0.savefig(directorySave+'MSD_average_short.png')
        fig0.savefig(directorySave+'MSD_average_short.svg',format='svg',dpi=1200)
        plt.close()  
        
        with open(directorySave+'MSD_average.csv', 'w') as textfile:
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
        if not os.path.exists(self.dn+'\\Plots\\MSD'):
            os.mkdir(self.dn+'\\Plots\\MSD')
        
        directorySave = self.dn+'\\Plots\\MSD\\'
        
        #Plot all MSDs (short version, 1/10 of data points)
        fig0 = plt.figure(0,figsize=(12, 10))
        for part in range(self.nbParticles):
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            plt.plot(timeD[:self.nbPoints],MSD[:self.nbPoints])
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$/s)')
        plt.axis('tight')
        plt.title('Short MSD for all particles')
        
        fig0.savefig(directorySave+'MSD_short_all.png')
        fig0.savefig(directorySave+'MSD_short_all.svg',format='svg',dpi=1200)
        plt.close()

        #Plots all MSDs (all points)
        fig0 = plt.figure(0,figsize=(12, 10))
        for part in range(self.nbParticles):
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            plt.plot(timeD,MSD)
            
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('MSD ($\mu$m$^2$/s)')
        plt.axis('tight')
        plt.title('Long MSD for all particles')
        
        fig0.savefig(directorySave+'MSD_long_all.png')
        fig0.savefig(directorySave+'MSD_long_all.svg',format='svg',dpi=1200)
        plt.close()

        #Creates the folder for individual particles if it's the first time
        if not os.path.exists(self.dn+'\\Plots\\MSD\\Individual\\'):
            os.mkdir(self.dn+'\\Plots\\MSD\\Individual\\')
        
        directorySave = self.dn+'\\Plots\\MSD\\Individual\\'
        
        #Plots all MSDs one by one (short version)
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            plt.plot(timeD[:self.nbPoints],MSD[:self.nbPoints])
            
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSD ($\mu$m$^2$/s)')
            plt.axis('tight')
            plt.title('Short MSD for particle ' + str(part))
            
            fig0.savefig(directorySave+'MSD_short_P'+str(part)+'.png')
            fig0.savefig(directorySave+'MSD_short_P'+str(part)+'.svg',format='svg',dpi=1200)
            plt.close()

        #Plots all MSDs one by one (long version)
        for part in range(self.nbParticles):
            
            fig0 = plt.figure(0,figsize=(12, 10))
        
            timeD = self.validParticles[part].timeD
            MSD = self.validParticles[part].MSD
            
            plt.plot(timeD,MSD)
            
            plt.xlabel('$\Delta$t (s)')
            plt.ylabel('MSD ($\mu$m$^2$/s)')
            plt.axis('tight')
            plt.title('Long MSD for particle ' + str(part))
            
            fig0.savefig(directorySave+'MSD_long_P'+str(part)+'.png')
            fig0.savefig(directorySave+'MSD_long_P'+str(part)+'.svg',format='svg',dpi=1200)
            plt.close()
        
        for p in range(self.nbParticles):
            with open(directorySave+'MSD_P'+str(p)+'.csv', 'w') as textfile:
                textfile.write('File,'+self.validParticles[p].fileName+'\n')
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
        

if __name__ == '__main__':
    root = tkinter.Tk()
    gui = GUI(root)
    root.title("Motion Visualization Tool v0.3")
    w = 350
    h = 700
    x = 200
    y = 200
    root.geometry('%dx%d+%d+%d' % (w, h, x, y))
#    root.geometry("300x280+200+200")
    root.mainloop()



