# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 15:16:02 2022

@author: jced0001
"""

import numpy as np
import nanonispy as nap
import ntpath
from ase.visualize.plot import plot_atoms
import copy
from Panel import Panel
import tkinter as tk
import customtkinter as ctk
import math

class GridPanel(Panel):
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width, height, dpi, mainPanel):
        super().__init__(master, width, height, dpi, mainPanel=mainPanel)
        self.init()
        self.buttons()
        
    ###########################################################################
    # Initialisation
    ###########################################################################
    def init(self):
        self.gridData = []                                                      # 3ds grid data goes here
        self.bbox = []                                                          # Box parameters for plotting outline of the grid on main panel overlay
        self.ychannel = 'LI Demod 1 X (A)'                                      # Default channel is the demod channel for now
        self.currentExtractPos = []                                             # x,y coorinates of the point spectra we with to extract while moving the mouse
        self.extractPos = []                                                    # Array of x,y coordinates within the grid to extract point spectra from
        self.currentAveragedPointsPos = []                                      # Array of x,y coordinates within the grid to extract and average spectra from
        self.averagedPointsPos = []                                             # List of arrays containing x,y coordinates within the grid to extract and average spectra from
        self.showAllAveragedPoints = True
        self.bias = []                                                          # Contains the index of the bias/sweep signal slice
        self.atoms = []                                                         # Initialise atoms variable (for plotting molecules/atoms on overlay)
        self.plotCaption = True                                                 # Show the plot caption by default
        self.scaleBar = True                                                    # Show a scale bar by default
        self.bound = False                                                      # Flag to prevent binding two functions at once
        
    ###########################################################################
    # Panel
    ###########################################################################
    def buttons(self):
        self.btn = {
            "load3ds":  ctk.CTkButton(self.master, text="Load 3ds",   command=self.load3ds),            # Load a 3ds file
            "Channel":  ctk.CTkButton(self.master, text="Current (A)",command=self.cycleChannel),
            "cmap":     ctk.CTkButton(self.master, text="viridis",    command=super()._cmap),           # Button to cycle through colour maps
            "Inset":    ctk.CTkButton(self.master, text="Inset",      command=super().addInset),        # Inset this into the main panel
            "Imprint":  ctk.CTkButton(self.master, text="Imprint",    command=super()._imprint),        # Imprint the spectra on STS panel
            "PullMol":  ctk.CTkButton(self.master, text="Copy Atoms", command=self.pullMolecules),      # Pull in molecules plotted on the main panel
            "ClearMol": ctk.CTkButton(self.master, text="Clear Atoms",command=self.clearMolecules),     # Clear all molecules from grid overlay
            "Caption":  ctk.CTkButton(self.master, text="Caption",    command=self._toggleCaption),     # Toggle the plot caption
            "ScaleBar": ctk.CTkButton(self.master, text="Scale Bar",  command=self._toggleScalebar),    # Toggle the plot scalebar
            "Close":    ctk.CTkButton(self.master, text="Close",      command=self.destroy)             # Close the panel
            }
            
    def buttonHelp(self):
        helpStr = "Load in a .3ds file containing grid data.\nFor cloud measurements, use the STS Panel"
        self.btn['load3ds'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Change the colour map"
        self.btn['cmap'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Copy atoms drawn on the main figure. Only those within the grid region will be displayed"
        self.btn['PullMol'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Clear atoms drawn on this figure"
        self.btn['ClearMol'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Show/Hide the plot caption"
        self.btn['Caption'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Add the above plot as an inset on the main figure. Double click a location in the main figure to repoisition the inset and use the scroll wheel to change its size"
        self.btn['Inset'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Imprint the overlay drawn by this panel on the main figure so it persits after closing this panel"
        self.btn['Imprint'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Close this panel"
        self.btn['Close'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Show/Hide the scalebar"
        self.btn['ScaleBar'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Adjust the slider to scan bias"
        self.slider.bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
    
    def special(self):                                                          # Special canvas UI
        self.slider = ctk.CTkSlider(self.master, orient=tk.HORIZONTAL, from_=0, to=1, width=420, command=self.changeBias) # Slider to select which bias/sweep signal slice to look show
        self.slider.grid(row=10,column=self.pos,columnspan=4,rowspan=2)          # Make it take up the entire length of the panel
        if(self.gridData):                                                      # Figure out how many ticks the slider needs. 1 tick = 1 dV
            to = len(self.gridData.signals['sweep_signal']) - 1
            self.slider.configure(to=to)
    
    def removeSpecial(self):
        self.slider.grid_forget()                                               # Called when panel is closed
        
    def changeBias(self,event):
        self.bias = int(event)                                                  # Change the bias on a slider event
        self.mainPanel.update(upd=[3,5])                                        # Update this panel and the STS panel (to show the vertical dashed line at selected bias)
        
    ###########################################################################
    # Update and Plotting
    ###########################################################################
    def update(self):
        if(not self.mainPanel.init): return
        if(not self.gridData):
            self.updateHelpLabel("Start by loading a 3ds file")
            return
        
        if(self.helpText == ""):
            if(not self.mainPanel.stsPanel.active):
                self.updateHelpLabel("Open the STS Panel to further interact with the grid")
            
        self.ax.cla()                                                           # Clear the axis
        self.plot()
        self.updateOverlay()
        self.ax.axis('off')                                                     # Hide the axis
        self.ax.set_position([0, 0, 1, 1])                                      # Make it take up the whole canvas
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()                                                      # Redraw the canvas with the updated figure
    
    def plot(self):
        Y = np.flipud(self.data[1])
        cmap = self.cmaps[self.cmap][1]
        im = Y[:,:,self.bias]
        vmin = np.min(im[np.nonzero(im)])
        vmax = np.max(im[np.nonzero(im)])
        self.ax.imshow(im, extent=self.extent,cmap=cmap(),vmin=vmin,vmax=vmax)
    
    def updateOverlay(self):
        left, right = self.ax.get_xlim();                                       # Remember plot extent
        bottom, top = self.ax.get_ylim();                                       # Remember plot extent
        self.plotMolecules()
        self.markExtractions()                                                  # Mark locations where point spectra have been extracted
        self.markAveragedExtractions()                                          # Mark locations where the point spectra are being averaged
        if(self.scaleBar): super().addPlotScalebar()                            # Add a scale bar to the plot
        if(self.plotCaption):                                                   # Caption the image with Vbias and Iset
            plotCaption  = r'V$_{bias}$ = '  + "{:.3f}".format(self.getBias()) + ' V' # Show bias in a box in the top left of the image
            offset  = self.gridData.header['dim_px'][1]
            offset /= self.gridData.header['dim_px'][0]
            offset  = 1 - offset
            offset /= 10                                                        # Add a smal offset to the plot caption to keep it within the image, when the image is wider than it is tall
            if(offset < 0): offset = 0                                          # Don't want the plot caption to start creeping over the top of the image
            pos = [0.025*self.extent[1],(0.94-offset)*self.extent[3]]           # Put the caption in the top left
            super().addPlotCaption(plotCaption, pos)                            # Add the caption to the plot
        self.ax.set_xlim((left,right)); self.ax.set_ylim((bottom,top))          # Put back extent
        
    def plotMolecules(self):
        for m in self.atoms:
            pop = []
            for idx,a in enumerate(m.positions):
                if(a[0] < 0 or a[1] < 0):                                       # Weird thing happens when plotting a molecule and some/all atoms are out of the image extent
                    pop.append(idx)                                             # so find any atoms outside image area
            
            for p in np.flip(pop):
                m.pop(p)                                                        # and get rid of them (in reverse order so indexes don't change)
            if(len(m)): plot_atoms(m, self.ax, radii=0.5e-10)                   # Plot this molecule if there are any atoms left
                
    def markExtractions(self):
        for marker in self.extractPos:
            self.ax.plot(marker[0],marker[1],'.',markersize=12)                 # Plot a marker at each point we're extracting a spectrum
        
        if(len(self.currentExtractPos)):
            self.ax.plot(self.currentExtractPos[0],self.currentExtractPos[1],'x',markersize=12) # Plot a cross at the cursor position when selecting a point
    
    def markAveragedExtractions(self):
        idx = 0
        for idx,points in enumerate(self.averagedPointsPos):
            c = self.mainPanel.mplibColours[idx+1]                              #+1 because I don't wanna start from black
            for marker in points:
                self.ax.plot(marker[0],marker[1],'x',markersize=12,c=c)
                if(not self.showAllAveragedPoints): break                       # Quit the loop after the first point if we're not meant to be marking all locations
        
        c = self.mainPanel.mplibColours[idx+2]                                  # Go to the next colour for points currently being placed
        for marker in self.currentAveragedPointsPos:
            self.ax.plot(marker[0],marker[1],'x',markersize=12,c=c)
            
    def load3ds(self):
        path,filename = ntpath.split(self._browseFile())                        # Open file dialog to browse for a 3ds file
        if(not filename.endswith('.3ds')):
            print("Expecting 3ds file")
            return                                                              # 3ds files only
            
        header_override = {'Delay before measuring (s)':'0',                    # Not used, just stop nap from complaining
                           'Start time':'0',
                           'End time':'0',
                           'Comment':'none'}
        
        self.gridData = nap.read.Grid(path + '/' + filename, header_override=header_override) # Read the file. see gridData.header.keys() and gridData.signals.keys() for structure
        
        lxy = np.array(self.gridData.header['size_xy'])                         # Length by width in real space (m)
        self.extent = np.array([0,lxy[0],0,lxy[1]])                             # Extent of the data for plotting
        
        # Trig shit
        scanAngle = math.pi*self.mainPanel.scanAngle/180                        # Angle of the scan frame in radians
                                                                                # Take the bottom left corner of the scan frame as the 'origin'
        OA = self.mainPanel.lxy/2                                               # Vector from origin to centre of scan frame
        
        AB = np.array(self.gridData.header['pos_xy']) - np.array(self.mainPanel.sxm.header['scan_offset']) # Vector from centre of scan frame to centre of grid
        AB = np.array(self.rotate([0,0],AB,scanAngle))                          # Rotate this vector about 0,0 to get it in the same basis as scan frame
        
        C = -lxy/2                                                              # Vector from centre of grid to bottom left corner of grid
        
        V = OA + AB + C                                                         # Vector from bottom left corner of scan frame to bottom left corner of grid
        ox,oy = V
        origin = OA + AB                                                        # Vector from bottom left of scan frame to centre of grid
        self.bbox = {"bbox"   : [V,lxy[0],lxy[1]],
                     "origin" : origin,
                     "angle"  : -self.gridData.header['angle']}                 # Outline of the grid w.r.t the sxm TOPO
        
        self.bias = 0                                                           # Reset the bias/sweep signal index upon load
        
        to = len(self.gridData.signals['sweep_signal']) - 1                     # The index of the last sweep signal point
        self.slider.configure(to=to)                                            # Update the slider selector on the gui to cover index range of loaded 3ds file
        
        self.updateHelpLabel("")
        self.data = []
        self.data.append(self.gridData.signals['sweep_signal'])
        self.data.append(copy.deepcopy(self.gridData.signals[self.ychannel]))
        self.mainPanel.update(upd=[0,3,5])                                      # Update the main panel, sts panel, and this panel
        
    ###########################################################################
    # Extract point STS from grid
    ###########################################################################
    def extractBind(self):
        if(self.bound): return
        self.bound = True
        self.lcCursorBind     = self.canvas.get_tk_widget().bind('<Button-1>', self.addExtract)     # Bind the left click for adding a point
        self.rcCursorBind     = self.canvas.get_tk_widget().bind('<Button-3>', self.cancelExtract)  # Bind the right click for cancelling during selection
        self.motionCursorBind = self.canvas.get_tk_widget().bind('<Motion>', self.motionExtract)    # Bind the motion of the mouse for placing a point
        
        self.updateHelpLabel("Select a pixel within the grid to plot its spectrum on the STS panel\nLeft click to select, right click to cancel")
        
    def extractUnbind(self):                                                    # Unbind all controls
        self.canvas.get_tk_widget().unbind('<Button-1>', self.lcCursorBind)
        self.canvas.get_tk_widget().unbind('<Button-3>', self.rcCursorBind)
        self.canvas.get_tk_widget().unbind('<Motion>', self.motionCursorBind)
        self.bound = False
        
        self.updateHelpLabel("")
        
    def addExtract(self,event):
        self.extractUnbind()                                                    # Unbind, we've just selected a point
        self.extractPos.append(self.currentExtractPos)                          # Add the location of the click to the list of extraction locations
        self.currentExtractPos = []                                             # Clear the buffer used when placing a point at the cursor location
        self.mainPanel.update(upd=[3,5])                                        # Update the sts panel and this panel
    
    def motionExtract(self,event):
        size = self.fig.get_size_inches()*self.fig.dpi                          # size in pixels
        x = event.x
        y = size[1] - event.y
        X = super()._getX(x)
        Y = super()._getY(y)
        
        self.currentExtractPos = np.array([X,Y])                                # Refresh with current cursor location
        self.mainPanel.update(upd=[3,5])                                        # Update this panel to reposition x marker, update sts panel to replot sts from cursor location
    
    def cancelExtract(self,event):
        self.extractUnbind()                                                    # Unbind, we've cancelled placement of this marker
        self.currentExtractPos = []                                             # Clear buffer
        self.mainPanel.update(upd=[3,5])                                        # Update sts panel and this panel
        
    def getPointSpectra(self):
        if(not self.gridData): return [],[]                                     # If we're here, we haven't loaded a 3ds file yet
        indexes = self.getGridIndexes()                                         # Convert x,y coordinates to the x,y indexes of point spectra in the 3D array
        spectra = []                                                            # Array to store the extracted point spectra
        for idx in indexes:
            spectra.append(self.data[1][idx[1],idx[0],:])                       # Append point spectra from desired locations
        
        return self.data[0],copy.deepcopy(spectra)
    
    def getAveragedPointSpectra(self):
        if(not self.gridData): return [],[]                                     # If we're here, we haven't loaded a 3ds file yet
        indexList = self.getAveragedGridIndexes()
        spectra = []
        sweep = self.gridData.signals['sweep_signal']                           # The signal swept during acquisition (typically bias)
        for indexes in indexList:
            averagedSpectrum = 0*sweep
            for idx in indexes:
                averagedSpectrum += self.data[1][idx[1],idx[0],:]
            spectra.append(averagedSpectrum/len(indexes))
        
        return copy.deepcopy(sweep),copy.deepcopy(spectra)
    
    def getGridIndexes(self):                                                   # This function converts real-space x,y coords into x,y array indicies
        lxy = np.array(self.gridData.header['size_xy'])
        indexes = []
        for pos in self.extractPos:                                             # Do this for all the previously selected points
            x_idx = int((pos[0]/lxy[0]) * self.gridData.header['dim_px'][0])
            y_idx = int((pos[1]/lxy[1]) * self.gridData.header['dim_px'][1])
            indexes.append(np.array([x_idx,y_idx]))
        
        if(len(self.currentExtractPos)):                                        # And also the mouse is hovering over, if we're currently choosing a point
            pos = self.currentExtractPos
            x_idx = int((pos[0]/lxy[0]) * self.gridData.header['dim_px'][0])
            y_idx = int((pos[1]/lxy[1]) * self.gridData.header['dim_px'][1])
            indexes.append(np.array([x_idx,y_idx]))
            
        return indexes
    
    def getAveragedGridIndexes(self):                                           # This function converts real-space x,y coords into x,y array indicies
        lxy = np.array(self.gridData.header['size_xy'])
        indexList = []
        for points in self.averagedPointsPos:                                   # Do this for all the previously selected points
            indexes = []
            for pos in points:
                x_idx = int((pos[0]/lxy[0]) * self.gridData.header['dim_px'][0])
                y_idx = int((pos[1]/lxy[1]) * self.gridData.header['dim_px'][1])
                indexes.append(np.array([x_idx,y_idx]))
            indexList.append(copy.deepcopy(indexes))    
        return indexList
        
    ###########################################################################
    # Average several spectra from points within grid
    ###########################################################################
    def averageGridPointsBind(self):
        if(self.bound): return
        self.bound = True
        self.lcCursorBind     = self.canvas.get_tk_widget().bind('<Button-1>', self.addGridPointToAverage)  # Bind the left click for adding a point
        self.rcCursorBind     = self.canvas.get_tk_widget().bind('<Button-3>', self.finishAveraging)        # Bind the right click for cancelling during selection
        self.motionCursorBind = self.canvas.get_tk_widget().bind('<Motion>',   self.motionExtract)          # Bind the motion of the mouse for placing a point (uses the function from extracting single grid points)
        
        self.updateHelpLabel("Select a group of pixels within the grid to plot their averaged spectra on the STS panel\nLeft click to select, right click within the grid to finish")
        
    def averageGridPointsUnbind(self):                                             # Unbind all controls
        self.canvas.get_tk_widget().unbind('<Button-1>', self.lcCursorBind)
        self.canvas.get_tk_widget().unbind('<Button-3>', self.rcCursorBind)
        self.canvas.get_tk_widget().unbind('<Motion>', self.motionCursorBind)
        self.bound = False
        
        self.updateHelpLabel("")
        
    def addGridPointToAverage(self,event):
        self.currentAveragedPointsPos.append(self.currentExtractPos)            # Add the location of the click to the list of extraction locations we're currently building
        self.currentExtractPos = []                                             # Clear the buffer used when placing a point at the cursor location
        self.mainPanel.update(upd=[3,5])                                        # Update the sts panel and this panel
    
    def finishAveraging(self,event):
        self.extractUnbind()                                                    # Unbind, we've cancelled placement of this marker
        if(len(self.currentAveragedPointsPos)):
            self.averagedPointsPos.append(self.currentAveragedPointsPos)        # This collection of points become averaged into a single spectra
        self.currentExtractPos = []                                             # Clear buffer
        self.currentAveragedPointsPos = []                                      # Clear buffer
        self.mainPanel.update(upd=[3,5])                                        # Update sts panel and this panel
        
    ###########################################################################
    # Pull Molecules from main panel
    ###########################################################################
    def pullMolecules(self):                                                    # This function pulls the atoms from the main panel and plots them over the grid data
        self.atoms = copy.deepcopy(self.mainPanel.atoms)                        # Deep copy because we're playing with the atoms' positions
        for m in self.atoms:
            X,Y    = self.bbox['bbox'][0]
            origin = self.bbox['origin']
            scanAngle = self.mainPanel.scanAngle
            angle  = math.pi*(self.bbox['angle'] + scanAngle)/180
            for pos in m.positions:
                pos[0:2] = self.rotate(origin,pos[0:2],-angle)                  # Rotate
                pos -= np.array([X,Y,0])                                        # And shift the atoms according to rotation and offset of grid frame w.r.t scan frame
        self.update()                                                           # Just update this panel
    
    def clearMolecules(self):
        self.atoms = []                                                         # Clear all the atoms
        self.update()                                                           # Just update this panel
        
    ###########################################################################
    # Misc
    ###########################################################################
    def cycleChannel(self):
        if(not self.gridData): return
        
        channels = list(self.gridData.signals.keys())
        
        if(not self.ychannel in channels):
            self.ychannel = channels[0]
            return
        
        idx = channels.index(self.ychannel) + 1
        if(idx == len(channels)): idx = 0
        self.ychannel = channels[idx]
        
        shape = list(self.gridData.signals[self.ychannel].shape)
        
        if(idx > 0):
            if(not len(shape) == 3):
                self.cycleChannel()
            elif(not shape[2] == len(self.data[0])):
                print("shape2:",shape[2])
                self.cycleChannel()
        
        self.smooth()
        self.btn['Channel'].configure(text=self.ychannel)
        self.update()
        
    def smooth(self):
        self.data    = []
        if(not self.gridData): return
        
        sweep = self.gridData.signals['sweep_signal']
        self.data.append(sweep)
        self.data.append(copy.deepcopy(self.gridData.signals[self.ychannel]))
        for idx,i in enumerate(self.data[1]):
            for jdx,c in enumerate(i):
                curve = self.mainPanel.stsPanel.getDIDV(curve=[sweep,c],ychannel=self.ychannel)
                self.data[1][idx][jdx] = curve[1]
        self.update()
        
    def getBias(self):
        if(not self.gridData): return np.nan
        return self.gridData.signals['sweep_signal'][self.bias]                 # Convert bias index to actual bias value
    
    def _toggleCaption(self):
        self.plotCaption = not self.plotCaption                                 # Toggles the display of the plot caption
        self.update()
        
    def _toggleScalebar(self):
        self.scaleBar = not self.scaleBar                                       # Toggles the display of the scale bar
        self.update()
        
    def _undo(self):
        if(len(self.extractPos)):
            self.extractPos.pop()                                               # Pop the last spectra off the back of the list to undo
            self.mainPanel.update(upd=[3,5])
            
    def undoAverage(self):
        if(len(self.averagedPointsPos)):
            self.averagedPointsPos.pop()
            self.mainPanel.update(upd=[3,5])
            
    def _reset(self):
        self.extractPos = []                                                    # Clear all the spectra we've selected
        self.averagedPointsPos = []
        self.mainPanel.update(upd=[3,5])