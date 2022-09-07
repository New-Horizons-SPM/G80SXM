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
            "load3ds":  tk.Button(self.master, text="Load 3ds",   command=self.load3ds),            # Load a 3ds file
            "cmap":     tk.Button(self.master, text="viridis",    command=super()._cmap),           # Button to cycle through colour maps
            "Inset":    tk.Button(self.master, text="Inset",      command=super().addInset),        # Inset this into the main panel
            "Imprint":  tk.Button(self.master, text="Imprint",    command=super()._imprint),        # Imprint the spectra on STS panel
            "Undo":     tk.Button(self.master, text="Undo Last",  command=self._undo),              # Undo the last point spectra on STS panel
            "Reset":    tk.Button(self.master, text="Reset",      command=self._reset),             # Remove all point spectra from STS panel
            "PullMol":  tk.Button(self.master, text="Pull Mol",   command=self.pullMolecules),      # Pull in molecules plotted on the main panel
            "ClearMol": tk.Button(self.master, text="Clear Mol",  command=self.clearMolecules),     # Clear all molecules from grid overlay
            "Caption":  tk.Button(self.master, text="Caption",    command=self._toggleCaption),     # Toggle the plot caption
            "ScaleBar": tk.Button(self.master, text="Scale Bar",  command=self._toggleScalebar),    # Toggle the plot scalebar
            "Close":    tk.Button(self.master, text="Close",      command=self.destroy)             # Close the panel
            }
    
    def special(self):                                                          # Special canvas UI
        self.slider = tk.Scale(self.master, orient=tk.HORIZONTAL, from_=0, to=1, length=420, command=self.changeBias) # Slider to select which bias/sweep signal slice to look show
        self.slider.grid(row=8,column=self.pos,columnspan=4,rowspan=2)          # Make it take up the entire length of the panel
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
        Y = self.gridData.signals[self.ychannel]
        Y = np.flipud(Y)
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
        if(self.scaleBar): super().addPlotScalebar()                            # Add a scale bar to the plot
        if(self.plotCaption):                                                   # Caption the image with Vbias and Iset
            plotCaption  = r'V$_{bias}$ = '  + "{:.3f}".format(self.getBias()) + ' V' # Show bias in a box in the top left of the image
            offset  = self.gridData.header['dim_px'][1]
            offset /= self.gridData.header['dim_px'][0]
            offset  = 1 - offset
            offset /= 10
            if(offset < 0): offset = 1
            pos = [0.025*self.extent[1],(0.94-offset)*self.extent[3]]           # Put the caption in the top left
            super().addPlotCaption(plotCaption, pos)
        self.ax.set_xlim((left,right)); self.ax.set_ylim((bottom,top))          # Put back extent
        
    def plotMolecules(self):
        for m in self.atoms:
            plot_atoms(m, self.ax, radii=0.5e-10)                               # Plot this molecule
                
    def markExtractions(self):
        for marker in self.extractPos:
            self.ax.plot(marker[0],marker[1],'.',markersize=12)                 # Plot a marker at each point we're extracting a spectrum
        
        if(len(self.currentExtractPos)):
            self.ax.plot(self.currentExtractPos[0],self.currentExtractPos[1],'x',markersize=12) # Plot a cross at the cursor position when selecting a point
            
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
        
        ox,oy = np.array(self.gridData.header['pos_xy']) - self.mainPanel.im_offset # Origin of the grid w.r.t the origin of the main panel in real space (m)
        origin = [ox,oy]
        ox = ox - lxy[0]/2                                                      # x origin of the bbox starts from the bottom left corner of the rectangle, not the middle
        oy = oy - lxy[1]/2                                                      # y origin of the bbox starts from the bottom left corner of the rectangle, not the middle
        self.bbox = {"bbox"   : [(ox,oy),lxy[0],lxy[1]],
                     "angle"  : -self.gridData.header['angle'],
                     "origin" : origin}                                         # Outline of the grid w.r.t the sxm TOPO
        
        self.bias = 0                                                           # Reset the bias/sweep signal index upon load
        
        to = len(self.gridData.signals['sweep_signal']) - 1                     # The index of the last sweep signal point
        self.slider.configure(to=to)                                            # Update the slider selector on the gui to cover index range of loaded 3ds file
        
        self.updateHelpLabel("")
        
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
            spectra.append(self.gridData.signals[self.ychannel][idx[1],idx[0],:]) # Append point spectra from desired locations
        
        sweep = self.gridData.signals['sweep_signal']                           # The signal swept during acquisition (typically bias)
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
            
    ###########################################################################
    # Pull Molecules from main panel
    ###########################################################################
    def pullMolecules(self):                                                    # This function pulls the atoms from the main panel and plots them over the grid data
        self.atoms = copy.deepcopy(self.mainPanel.atoms)                        # Deep copy because we're playing with the atoms' positions
        for m in self.atoms:
            X,Y    = self.bbox['bbox'][0]
            origin = self.bbox['origin']
            angle  = math.pi*self.bbox['angle']/180
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
    def getBias(self):
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
    
    def _reset(self):
        self.extractPos = []                                                    # Clear all the spectra we've selected
        self.mainPanel.update(upd=[3,5])