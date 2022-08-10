# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 17:52:07 2022

@author: jced0001
"""

from Panel import Panel
import tkinter as tk
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
import numpy as np

class AIMLPanel(Panel):
    cb = []
    labels = []
    loadedIm = []
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width, height, dpi, mainPanel):
        super().__init__(master, width, height, dpi, mainPanel=mainPanel)
        self.exportPath = "C:/Users/jced0001/Development/Data/Highlights/hBN/Dep8/Dep8b/LabelledData/"
        self.buttons()
    
    ###########################################################################
    # Initialisation
    ###########################################################################
    def init(self):
        self.nextSXM(0)
        
    ###########################################################################
    # Panel
    ###########################################################################
    def buttons(self):
        self.btn = {
            "Close"         : tk.Button(self.master, text="Close",      command=self.destroy),
            "Previous"      : tk.Button(self.master, text="Prev",       command=lambda: self.nextSXM(-1)),
            "Next"          : tk.Button(self.master, text="Next",       command=lambda: self.nextSXM(1)),
            "Save"          : tk.Button(self.master, text="Save",       command=self.save),
            "+1"            : tk.Button(self.master, text="+1",         command=lambda: self.addLabel("+1")),
            "white_flower"  : tk.Button(self.master, text="interesting",command=lambda: self.addLabel("white_flower")),
            "duel"          : tk.Button(self.master, text="close dbl",  command=lambda: self.addLabel("duel")),
            "bow_and_arrow" : tk.Button(self.master, text="far dbl",    command=lambda: self.addLabel("bow_and_arrow")),
            "crayon"        : tk.Button(self.master, text="dull tip",   command=lambda: self.addLabel("crayon")),
            "striped_pole"  : tk.Button(self.master, text="unstable",   command=lambda: self.addLabel("striped_pole")),
            "chart"         : tk.Button(self.master, text="bad area",   command=lambda: self.addLabel("chart")),
            "poop"          : tk.Button(self.master, text="fire",       command=lambda: self.addLabel("poop")),
            "AugFlip"       : tk.Button(self.master, text="AugFlip",    command=lambda: self.augmentType("AugFlip")),
            "AugRot"        : tk.Button(self.master, text="AugRot",     command=lambda: self.augmentType("AugRot")),
            "AugPlane"      : tk.Button(self.master, text="AugPlane",   command=lambda: self.augmentType("AugPlane")),
            "AugNoise"      : tk.Button(self.master, text="AugRot",     command=lambda: self.augmentType("AugNoise")),
            "ExportFolder"  : tk.Button(self.master, text="Save Path",  command=self.savePath)
            }
    
    ###########################################################################
    # Update and Plotting
    ###########################################################################
    def update(self):
        if(not self.mainPanel.init): return
        
        self.ax.cla()                                                           # Clear the axis
        if(self.cb): self.cb.remove()
        self.fig.axes.append(self.ax)
        
        image = self.loadedIm
        if(not len(image)): image = self.mainPanel.tiltedim/1e-9
        im = self.ax.imshow(image)
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        self.cb = self.fig.colorbar(im, cax=cax)
        
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()                                                      # Redraw the canvas with the updated figure
    
    def nextSXM(self,d):
        path,filename = os.path.split(self.mainPanel.filename)                  # Split the path and filename
        files = os.listdir(path)                                                # Get all files in current directory
        sxmFiles = [path + "/" + f for f in files if f.endswith(".sxm")]        # Get all the .sxm files
        sxmFiles.sort(key=os.path.getctime)                                     # Sort files by asc datetime
        
        idx = sxmFiles.index(path + "/" + filename) + d                         # Index of the next (or prev if d=-1) sxmFile to look at
        if(idx < 0 or idx > len(sxmFiles) - 1): return                          # Do nothing if we've hit the end of the line
        
        nextFile = sxmFiles[idx]                                                # Next sxm file to load
        
        self.labels = []
        self.loadedIm = []
        pickledFiles = os.listdir(self.exportPath)
        pickledFile = os.path.split(nextFile)[1] + ".pkl"
        if(pickledFile in pickledFiles):
            with open(self.exportPath + pickledFile, 'rb') as handle:
                labelledData = pickle.load(handle)
                self.labels  = labelledData['labels']
                self.loadedIm = labelledData['Z (m)']
        
        for btn in self.btn:
            self.btn[btn].configure(bg='SystemButtonFace')
        for btn in self.labels:
            self.btn[btn].configure(bg='red')
        
        if(not d): return
        
        self.mainPanel.loadSXM(nextFile)                                        # Load it in the main panel
    
    ###########################################################################
    # Data Augmentation
    ###########################################################################
    def augmentType(self,augType):
        if(augType in self.augTypes):
            del self.augTypes[self.augTypes.index(augType)]
            self.btn[augType].configure(bg="SystemButtonInterface")
            return
        self.augTypes.append(augType)
        self.btn[augType].configure(bg="Red")
    
    def augment(self,im):
        pass
    
    def augFlip(self,im):
        flipVert = np.flipud(im)
        flipHorz = np.fliplr(im)
        flipBoth = np.flipud(flipHorz)
        return [flipVert, flipHorz, flipBoth]
    
    def augRot(self,im):
        rot90  = np.rot90(im)
        rot180 = np.rot90(rot90)
        rot270 = np.rot90(rot180)
        return [rot90, rot180, rot270]
    
    def augPlane(self,im):
        pass
    def augNoise(self,im):
        pass
        
    ###########################################################################
    # Labelling and Saving
    ###########################################################################
    def addLabel(self,label):
        if(label in self.labels):
            del self.labels[self.labels.index(label)]
            self.btn[label].configure(bg='SystemButtonFace')
            return
        
        self.labels.append(label)
        self.btn[label].configure(bg='red')
    
    def savePath(self):
        path = super()._browseFolder()
        if(path):
            self.exportPath = path
            if(not path.endswith('/')): self.exportPath += '/'
    
    def save(self):
        if(not len(self.labels)):
            print("Cannot save unlabled data")
            return
        
        _,filename = os.path.split(self.mainPanel.filename)                     # Split the path and filename
        imdata  = self.mainPanel.tiltedim.copy()
        imdata -= np.min(imdata)
        imdata /= np.max(imdata)
        labelledData = {"filename" : filename,
                        "labels"   : self.labels,
                        "size"     : self.mainPanel.lxy,
                        "Z (m)"    : imdata
                        }
        pickleFilename = self.exportPath + filename + ".pkl"
        with open(pickleFilename, 'wb') as handle:
            pickle.dump(labelledData,handle)
        
        self.nextSXM(1)