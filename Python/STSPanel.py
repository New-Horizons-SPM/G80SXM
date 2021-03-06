# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 18:09:36 2022

@author: jced0001
"""

from Panel import Panel
import tkinter as tk
import numpy as np
import os
import nanonispy as nap
import math
from scipy.signal import savgol_filter as savgol
class STSPanel(Panel):
    datFile = []; stsPos = []; stsOffset = True                                 # list of .dat filenames. stspos: location of xy pos 
    logScale = False
    datFileCustom = []; customSTSPos = []
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width, height, dpi, mainPanel):
        super().__init__(master, width, height, dpi, mainPanel=mainPanel)
        self.buttons()
    ###########################################################################
    # Panel
    ###########################################################################
    def buttons(self):
        self.btn = {
            "Multi":    tk.Button(self.master, text="Add Multi",  command=self._browseMulti),
            "Single":   tk.Button(self.master, text="Add Single", command=self._browseSingle),
            "Custom":   tk.Button(self.master, text="Add Custom", command=self._browseCustom),
            "Offset":   tk.Button(self.master, text="Offset",     command=self._offset),
            "Scale":    tk.Button(self.master, text="Linear",     command=self._scale),
            "Undo":     tk.Button(self.master, text="Undo Last",  command=self._undo),
            "Reset":    tk.Button(self.master, text="Reset",      command=self._reset),
            "Inset":    tk.Button(self.master, text="Inset",      command=super().addInset),
            "Imprint":  tk.Button(self.master, text="Imprint",    command=super()._imprint),
            "Close":    tk.Button(self.master, text="Close",      command=self.destroy)
            }
    ###########################################################################
    # Update and Plotting
    ###########################################################################
    def update(self):
        if(not self.mainPanel.init): return
        
        self.ax.cla()                                                           # Clear the axis
        self._plotSTS()                                                         # Loops through .dat files, takes the derivative of IV curves and plot dI/dV
        self.ax.set_position([0.13, 0.1, 0.83, 0.83])                           # Leave room for axis labels and title
        
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()                                                      # Redraw the canvas with the updated figure
    
    def _plotSTS(self):
        dat_xchannel = 'Bias calc (V)'
        dat_ychannel = 'Current (A)'
        sg_pts  = 7; sg_poly = 1; num_offset = 3                                # These will eventually be user input #todo
        max_val = 0; offset = 0; cnt = 0                                        # Loop variables
        
        datFiles = self.datFile.copy()
        if(self.datFileCustom): datFiles.append(self.datFileCustom)             # Weird way to loop through combined list but only way I could figure out how to do it
        for df in datFiles:                                                     # Loop through each .dat file, get the IV curve, take the derivative and plot a filtered version of dI/dV
            dat = nap.read.Spec(df)                                             # dat Spec object
            
            V = dat.signals[dat_xchannel]
            I = dat.signals[dat_ychannel]
            dV = V[1] - V[0]
            
            didv = savgol(I,sg_pts,sg_poly,deriv=1,delta=dV)
            if(self.logScale): didv = np.log(didv); didv = didv - np.min(didv)
           
            self.ax.plot(V,didv + cnt*offset,linewidth=1.3)
           
            max_val = max(max_val,didv.max())
            if cnt == 0:                                                        # Only do this on the first iteration
               offset = num_offset*0.25*max_val*self.stsOffset                  # offset for the next curve
            cnt += 1
            
        self.ax.set_xlabel("Bias (V)")
        self.ax.set_ylabel(["dI/dV (arb)","log(dI/dV) (arb)"][self.logScale]);
        self.ax.set_title("Point Spectroscopy")
    ###########################################################################
    # Browsing STS Files
    ###########################################################################
    def _browseMulti(self):                                                     # Select folder containing STS spectra to display. Returns False if no files found, otherwise true
        stsPath      = super()._browseFolder()
        if(stsPath):
            file_list    = os.listdir(stsPath)
            self.stsPos  = []
            self.datFile = [stsPath + "/" + f for f in file_list if f.endswith(".dat")] # Get .dat filenames in selected directory
            if(self.datFile):
                for df in range(len(self.datFile)):
                    dat = nap.read.Spec(self.datFile[df])                       # dat Spec object
                    
                    # Location of STS on image
                    x = (np.array(dat.header['X (m)']).astype(float)*self.mainPanel.pixelCalibration[0] - self.mainPanel.im_offset[0])
                    y = (np.array(dat.header['Y (m)']).astype(float)*self.mainPanel.pixelCalibration[1] - self.mainPanel.im_offset[1])
                    
                    ox,oy = (0,0)                                               # Origin before rotation is at 0,0 (bottom left)
                    rotAbout = self.mainPanel.lxy/2                             # Nanonis rotates frame about centre
                    angle = -math.pi*self.mainPanel.scanAngle/180.0             # Angle of the frame in nanonis (convert to rad)
                    oX,oY = super().rotate(rotAbout,(ox,oy),angle)              # Rotate the origin about the centre of the frame by angle
                    X = x - oX                                                  # x w.r.t new origin is X
                    Y = y - oY                                                  # y w.r.t new origin is Y
                    Xb = X*math.cos(angle) + Y*math.cos(math.pi/2 - angle)      # Xb is x w.r.t new basis
                    Yb = Y*math.cos(angle) + X*math.cos(math.pi/2 + angle)      # Yb is y w.r.t new basis
                    
                    self.stsPos.append([Xb,Yb])
                    self.mainPanel.update(upd=[0,3])
                return True
            else:
                print("No .dat files in folder")
                return False
    
    def _browseSingle(self):
        stsPath      = super()._browseFile()
        if(not stsPath.endswith(".dat")):                                       # Needs to be dat file
            print("Expecting .dat file")
            return                                                              # Return if no SXM file chosen
        self.datFile.append(stsPath)
        
        dat = nap.read.Spec(stsPath)                                            # dat Spec object
        
        # Location of STS on image
        x = (np.array(dat.header['X (m)']).astype(float)*self.mainPanel.pixelCalibration[0] - self.mainPanel.im_offset[0])
        y = (np.array(dat.header['Y (m)']).astype(float)*self.mainPanel.pixelCalibration[1] - self.mainPanel.im_offset[1])
        
        ox,oy = (0,0)                                                           # Origin before rotation is at 0,0 (bottom left)
        rotAbout = self.mainPanel.lxy/2                                         # Nanonis rotates frame about centre
        angle = -math.pi*self.mainPanel.scanAngle/180.0                         # Angle of the frame in nanonis (convert to rad)
        oX,oY = super().rotate(rotAbout,(ox,oy),angle)                          # Rotate the origin about the centre of the frame by angle
        X = x - oX                                                              # x w.r.t new origin is X
        Y = y - oY                                                              # y w.r.t new origin is Y
        Xb = X*math.cos(angle) + Y*math.cos(math.pi/2 - angle)                  # Xb is x w.r.t new basis
        Yb = Y*math.cos(angle) + X*math.cos(math.pi/2 + angle)                  # Yb is y w.r.t new basis
        
        self.stsPos.append([Xb,Yb])
        self.mainPanel.update(upd=[0,3])
    ###########################################################################
    # Custom STS Files
    ###########################################################################
    def _browseCustom(self):
        stsPath      = self._browseFile()
        if(not stsPath.endswith(".dat")):                                       # Needs to be dat file
            print("Expecting .dat file")
            return                                                              # Return if no SXM file chosen
        self.datFileCustom = stsPath
        self.mainPanel.customSTSBind()                                          # Follow logic for cursor bind in LineProfilePanel
        self.update()
        
    def setMarkSTS(self,stsPos,setMarker=False):
        self.customSTSPos = stsPos
        if(setMarker):
            self.datFile.append(self.datFileCustom)
            self.stsPos.append(stsPos)
            self.datFileCustom = []
            self.customSTSPos  = []
    
    def cancelMarkSTS(self):
        self.datFileCustom = []
        self.customSTSPos  = []
        self.update()
    ###########################################################################
    # Misc Button Functions
    ###########################################################################
    def _undo(self):
        self.stsPos  = self.stsPos[0:-1]
        self.datFile = self.datFile[0:-1]
        self.mainPanel.update(upd=[0,3])
        
    def _reset(self):
        self.stsPos = []
        self.datFile = []
        self.mainPanel.update(upd=[0,3])
    
    def _scale(self):
        self.logScale = not self.logScale
        text = ["Linear","Log"]
        self.btn['Scale'].configure(text=text[self.logScale])
        self.update()
    
    def _offset(self):
        self.stsOffset = not self.stsOffset
        self.update()
    ###########################################################################
    # Save (WIP)
    ###########################################################################
    def buildSaveString(self):
        saveString = "#STSPanel\n"                                              # Line 1: Header
        
        saveString += str(self.stsOffset) + "\n"                                # Line 2: offset curves
        
        saveString += str(self.logScale)  + "\n"                                # Line 3: log scale = 1. linear scale = 0
        
        stsPos   = []
        for i in range(len(self.datFile)):
            stsPos.append(','.join(str(pos) for pos in self.stsPos[i]))
        saveString += ('|'.join(i for i in self.datFile) + "\n")                # Line 4: List of sts .dat files
        saveString += ('|'.join(i for i in stsPos) + "\n")                      # Line 5: Their position on the sxm
    ###########################################################################
    # Load (WIP)
    ###########################################################################