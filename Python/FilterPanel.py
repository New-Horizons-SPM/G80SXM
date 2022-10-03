# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 18:09:36 2022

@author: jced0001
"""

from Panel import Panel
import customtkinter as ctk
import numpy as np
from scipy.signal import convolve2d
from scipy import ndimage

class FilterPanel(Panel):
    activeFilters = {"rollV" : 0,
                     "rollH" : 0}
    filterOrder = []; activeFilterOrder = []                                    # Implement this at some point
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width, height, dpi, mainPanel):
        super().__init__(master, width, height, dpi, mainPanel=mainPanel)
        self.buttons()
        self.buildFilters()
    ###########################################################################
    # Panel
    ###########################################################################
    def buttons(self):
        self.btn = {
            "RollV+": ctk.CTkButton(self.master, text="Roll Vert +", command=lambda: self.updateFilter("rollV", 1)),
            "RollV-": ctk.CTkButton(self.master, text="Roll Vert -", command=lambda: self.updateFilter("rollV",-1)),
            "RollH+": ctk.CTkButton(self.master, text="Roll Horz +", command=lambda: self.updateFilter("rollH", 1)),
            "RollH-": ctk.CTkButton(self.master, text="Roll Horz -", command=lambda: self.updateFilter("rollH",-1)),
            "HPF+":   ctk.CTkButton(self.master, text="High Pass +", command=lambda: self.updateFilter("HP",1)),
            "HPF-":   ctk.CTkButton(self.master, text="High Pass -", command=lambda: self.updateFilter("HP",-1)),
            "SetFlt": ctk.CTkButton(self.master, text="Set Filter",  command=self.setFilter),
            "Close":  ctk.CTkButton(self.master, text="Close",       command=self.destroy)
            }
    ###########################################################################
    # Update and Plotting (WIP)
    ###########################################################################
    def update(self):
        if(not self.mainPanel.init): return
        
        self.ax.cla()                                                           # Clear the axis
        self._previewFilters()                                                  # Filter the image in this panel to show a preview before setting
        self._addPlotCaption()                                                  # Add a caption to the plot in the upper left to show what filters are previewing
        self.ax.set_position([0, 0, 1, 1])                                      # Make it take up the whole canvas
        self.ax.axis('off')                                                     # Hide the axis
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()                                                      # Redraw the canvas with the updated figure
    
    def _previewFilters(self):                                                  # Take a snapshot of the current sxm image and preview the filters being applied to it before setting them
        im = np.copy(self.mainPanel.unfilteredIm)                               # The snapshot of the current sxm image
        
        filteredIm = self.applyFilters(im)
        
        extent = self.mainPanel.extent
        cmap = self.cmaps[self.mainPanel.cmap][1]
        vmin = self.mainPanel.vmin; vmax = self.mainPanel.vmax
        self.ax.imshow(filteredIm,extent=extent,cmap=cmap(),vmin=vmin,vmax=vmax)
    
    def _addPlotCaption(self):
        plotCaption  = 'Filter List:'                                           # List of filters in the order they are applied (text box will display these)
        offset = 0                                                              # Used to offset lines in the informative textbox at the top left of the image
        for f in self.filterOrder:
            offset += 0.032                                                     # Appropriate offset between each line of new text in the text box
            plotCaption += "\n" + f + ": " + str(self.filters[f][0])
            
        extent = self.mainPanel.extent
        props = dict(boxstyle='round',facecolor='white',alpha=0.5)
        self.ax.text(0.025*extent[1],(0.95-offset)*extent[3],plotCaption,bbox=props)
    ###########################################################################
    # Filter Control
    ###########################################################################
    def buildFilters(self):
        self.filters = {"rollV" : [0, lambda im: self.rollVert(im)],
                        "rollH" : [0, lambda im: self.rollHorz(im)],
                        "HP"    : [0, lambda im: self.highPass(im)]}
        self.activeFilters = self.filters.copy()
        
    def updateFilter(self,filtType,inc):
        self.filters[filtType][0] += inc
        
        if(self.filters[filtType][0] == 0):
            del self.filterOrder[self.filterOrder.index(filtType)]              # Remove this filter from the order of applied filters if it has no contribution
            
        if(self.filters[filtType][0] < 0):
            self.filters[filtType][0] = 0                                       # Bottom out at zero
        
        if(self.filters[filtType][0] > 0 and not (filtType in self.filterOrder)):
            self.filterOrder.append(filtType)                                   # Append if this filter is in use
            
        self.update()
    
    def setFilter(self):
        self.activeFilters     = self.filters.copy()
        self.activeFilterOrder = self.filterOrder.copy()
        self.mainPanel.update()
    ###########################################################################
    # Filters
    ###########################################################################
    def applyFilters(self,im,active=False):
        filters    = self.filters.copy()
        filterList = self.filterOrder.copy()
        
        if(active):
            filters    = self.activeFilters.copy()
            filterList = self.activeFilterOrder.copy()
        
        filteredIm = im
        for i in filterList:                                                    # Loop through the filters we're previewing and apply them in the correct order
            for p in range(filters[i][0]):
                filteredIm = filters[i][1](filteredIm)
            
        # Some filters can change the number of pixels... easier to just pad out the image to what it was
        # pad = np.zeros(self.mainPanel.unfilteredIm.shape)
        # pad[:filteredIm.shape[0],:filteredIm.shape[1]] = filteredIm
        
        return filteredIm
    
    def rollVert(self,im,p=1):
        return convolve2d(im, np.ones((2*p,1)) / 2 / p, mode='same')
       
    def rollHorz(self,im,p=1):
        return convolve2d(im, np.ones((1,2*p)) / 2 / p, mode='same')
    
    def highPass(self,im):                                                      # Gaussian hp filter taken from https://stackoverflow.com/questions/6094957/high-pass-filter-for-image-processing-in-python-by-using-scipy-numpy
        lowpass = ndimage.gaussian_filter(im, 1)
        gauss_highpass = im - lowpass
        return gauss_highpass
    ###########################################################################
    # Save (WIP)
    ###########################################################################
    def buildSaveString(self):
        saveString = "#FilterPanel\n"                                           # Line 1: Header
        
        # Save self.activeFilterOrder
        # Save self.activeFilters[f][0] for f in self.activeFilterOrder
        
        return saveString