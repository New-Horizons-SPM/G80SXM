# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 15:00:50 2022

@author: jced0001
"""

from Panel import Panel
import tkinter as tk
from   tkinter import filedialog
from LineProfilePanel import LineProfilePanel
from FFTPanel import FFTPanel
from STSPanel import STSPanel
from FilterPanel import FilterPanel
# from AIMLPanel import AIMLPanel
from GridPanel import GridPanel
import numpy as np
import nanonispy as nap
import nanonispyfit as napfit
import math
import os
import io
import pickle
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib
from ase.io import read
from ase.visualize.plot import plot_atoms
matplotlib.use("TkAgg")

class MainPanel(Panel):
    init = False
    
    # Main Figure
    plotChannel = 'Z'
    scaleBar = True; plotCaption = True
    mplibColours = ['black'] + plt.rcParams['axes.prop_cycle'].by_key()['color'] + ['white']
    
    # Inset
    moveInsetActive = False
    insetPos = np.array([0.65, 0.65, 0.3, 0.3]);
    insetColours = ['black','white'] + plt.rcParams['axes.prop_cycle'].by_key()['color']
    insetCmap = 0
    
    # Tilt
    tiltActive= False
    # tiltFactor= 1e-14                                                           # tiltFactor will be user input eventually.
    tiltFactor= 100                                                             # tiltFactor will be user input eventually.
    curTilt   = np.array([[1/2, 0, 0],[0,1,0],[1,1,0]]);                        # Tilt plane is defined by three points. curTilt keeps track of a tilt correction in progress
    tiltPlane = np.array([[1/2, 0, 0],[0,1,0],[1,1,0]]);
    
    # Plane Fit
    planeFitArea = []; curPlaneFitArea = [];
    planeFitCursor = -1
    
    # Draw Atoms
    atoms = []; curMol = []                                                     # curMol is the filename of the molecule being placed
    courseRot = 10; fineRot = 1                                                 # course and fine rotation precision in degrees
    molFiles = []; molPos = []; molRot = []; molRotX = []; molRotY = []         # list of filenames and corresponding positions of the molecules on the canvas.
    curRot = 0; curRotX = 0; curRotY = 0;
    
    # Alternate function buttons
    shiftL = False
    
    # Misc
    piezoFactor=np.array([1,1])
    bound = False
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width=512, height=512, dpi=96):
        super().__init__(master, width=width, height=height, dpi=dpi)
        self.buildSubPanels()
        self.buttons()
        super().create()
        self.canvas.figure = self.fig
        self.canvas.draw()
        
        # Always bound
        self.canvas.get_tk_widget().bind("<MouseWheel>", self._mouseWheel)      # Used for adjusting colour range on sxm image
    
    ###########################################################################
    # Panel
    ###########################################################################
    def buildSubPanels(self):                                                   # Sub panels are build in here. add new ones here first
        commonParams = [self.master, self.width, self.height, self.dpi, self]
        
        self.linePanel = LineProfilePanel(*commonParams)                        # 1D Line Plot: Plots 1D lines through a cursor
        self.fftPanel  = FFTPanel(*commonParams)
        self.stsPanel  = STSPanel(*commonParams)
        self.fltPanel  = FilterPanel(*commonParams)
        # self.aimlPanel = AIMLPanel(*commonParams)
        self.gridPanel = GridPanel(*commonParams)
        
    def buttons(self):
        self.btn = {
            "Plot":     tk.Button(self.master, text="Load SXM", command=self.loadSXM),              # Button to plot a new sxm
            "DrawHAT":  tk.Button(self.master, text="HAT.xyz",  command=lambda: self.placeMolecule("../xyz/HAT.xyz")), # Button to place a HAT molecule
            "DrawDCA":  tk.Button(self.master, text="DCA.xyz",  command=lambda: self.placeMolecule("../xyz/DCA.xyz")), # Button to place a DCA molecule
            "DrawCu":   tk.Button(self.master, text="Cu.xyz",   command=lambda: self.placeMolecule("../xyz/Cu.xyz")),  # Button to place a Cu atom
            "Undo":     tk.Button(self.master, text="Undo",     command=self.undoMolecule),         # Button to place a Cu atom
            "cmap":     tk.Button(self.master, text="viridis",  command=super()._cmap),             # Button to cycle through colour maps
            "Tilt":     tk.Button(self.master, text="Tilt",     command=self.tilt),                 # Button to tilt correct an image using arrow keys
            "PlaneFit": tk.Button(self.master, text="Plane Fit",command=self.planeFit),             # Button to plane fit an area and subtract from image (like tilt but auto)
            "Shift":    tk.Button(self.master, text="Shift",    command=self.toggleShiftL),         # Button to toggle shift mode
            "Caption":  tk.Button(self.master, text="Caption",  command=self._toggleCaption),       # Toggle the plot caption
            "Profiles": tk.Button(self.master, text="Profiles", command=self.linePanel.create),     # Button to activate 1D Profiles panel
            "FFT":      tk.Button(self.master, text="FFT",      command=self.fftPanel.create),      # Button to activate FFT panel
            "STS":      tk.Button(self.master, text="STS",      command=self.stsPanel.create),      # Button to activate STS panel
            "Filter":   tk.Button(self.master, text="Filter",   command=self.fltPanel.create),      # Button to activate Filter panel
            # "AIML":     tk.Button(self.master, text="LabelMode",command=self.aimlPanel.create),     # Button to activate AI machine learning data labeler
            "Grid":     tk.Button(self.master, text="STSGrid",  command=self.gridPanel.create),     # Button to activate AI machine learning data labeler
            "InsetCol": tk.Button(self.master, text="Inset Col",command=self._insetCmap),           # Change the inset font and line colours
            "RemInset": tk.Button(self.master, text="Rem Inset",command=self._removeInset),         # Remove the inset from main panel
            "FlipIm":   tk.Button(self.master, text="Flip scan",command=self._flipScan),            # Flip the current scan
            "RotIm":    tk.Button(self.master, text="Rotate",   command=self._rotateScan),          # Flip the current scan
            "Channel":  tk.Button(self.master, text="Z (m)",    command=self._channel),             # Flip the current scan
            "Save":     tk.Button(self.master, text="Save",     command=self._save),                # Save all active panels to a .g80 file
            "PNG":      tk.Button(self.master, text="Exp PNG",  command=self._exportPNG),           # Export the canvas to png
            "Load":     tk.Button(self.master, text="Load",     command=self._load),                # Load a .g80 file
            "Quit":     tk.Button(self.master, text="Quit",     command=self.quit)                  # Button to quit the program
            }
    ###########################################################################
    # Main Panel Updates
    ###########################################################################
    def update(self,upd=[-1]):
        if(not self.init):
            self.updateHelpLabel("Click 'Load SXM' to get started")
            return                                                              # Get outta here if no sxm is loaded yet
        
        if(self.helpText == ""):
            self.updateHelpLabel("Scroll over image to adjust contrast.\n"
                                 + "Shift inactive: vmax\nShift active: vmin")
            
        if(-1 in upd or 0 in upd):                                              # upd=0 selects mainPanel
            self.ax.cla()                                                       # Clear the axis
            self._updateSXM()                                                   # Show the processed (tilt-corrected, filtered, etc.) sxm image
            self.ax.set_position([0, 0, 1, 1])                                  # Make it take up the whole canvas
            self.ax.axis('off')                                                 # Hide the axis
            
        if((-1 in upd or 1 in upd) and self.linePanel.active):                  # upd=1 selects lineProfile panel
            self.linePanel.update()
            
        if((-1 in upd or 2 in upd) and self.fftPanel.active):                   # upd=2 selects fftPanel
            self.fftPanel.update()
        
        if((-1 in upd or 3 in upd) and self.stsPanel.active):                   # upd=3 selects stsPanel
            self.stsPanel.update()
            
        if(-1 in upd or 4 in upd):                                              # upd=4 selects aimlPanel
            # self.aimlPanel.update()
            pass
            
        if((-1 in upd or 5 in upd) and self.gridPanel.active):                  # upd=5 selects gridPanel
            self.gridPanel.update()
            
        if(-1 in upd or 0 in upd):                                              # upd=0 selects mainPanel. Update the overlay at the end, since other panels may contribute to the overlay
            self._updateOverlay()                                               # Add things to the foreground (e.g. scalebar and plot label)
        
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()                                                      # Redraw the canvas with the updated figure
        
    def _updateSXM(self):
        self.planeFitIm = self._planeFit(self.im)
        
        self.tiltedim = self.planeFitIm + self._tiltPlane()                     # Adjust the raw image according to tilt correction
        
        self.unfilteredIm = self.tiltedim
        
        self.finalim = self.fltPanel.applyFilters(self.unfilteredIm,True)
        
        cmap = self.cmaps[self.cmap][1]
        self.ax.imshow(self.finalim,extent=self.extent,cmap=cmap(),
                       vmin=self.vmin,vmax=self.vmax)                           # Show the processed sxm image
    
    def _updateOverlay(self):
        left, right = self.ax.get_xlim();                                       # Remember plot extent
        bottom, top = self.ax.get_ylim();                                       # Remember plot extent
        self._markSTS()
        self._drawCursor()
        self._drawFitArea()
        self._drawMolecules()
        self._drawGridBox()
        if(self.scaleBar): super().addPlotScalebar()                            # Add a scale bar to the plot
        if(self.plotCaption):                                                   # Caption the image with Vbias and Iset
            plotCaption  = r'V$_{bias}$ = '  + str(self.bias)    + ' V'         # Show bias in a box in the top left of the image
            plotCaption += r'; I$_{set}$ = ' + str(self.current) + ' pA'        # Show setpoint current in top left
            pos = [0.025*self.extent[1],0.94*self.extent[3]]                    # Put the caption in the top left
            super().addPlotCaption(plotCaption, pos)
    
        self.ax.set_xlim((left,right)); self.ax.set_ylim((bottom,top))          # Put back extent
    
    def loadSXM(self, filename=None, load=False):
        if(not filename):
            filename = self._browseFile()
            if(not filename.endswith(".sxm")):                                  # Needs to be an SXM
                print("Expecting .sxm file")
                return                                                          # Return if no SXM file chosen
        
        self.filename = filename                                                # Only get here if a .sxm file is chosen
        self.sxm     = nap.read.Scan(self.filename)                             # Everything about the sxm
        self.bias    = np.array(self.sxm.header['bias'])                        # Image bias (V)
        self.current = np.array(self.sxm.header['z-controller']['Setpoint'][0]. # Image current (1e12: convert to pA)
                                split()[0]).astype(float)*1e12                  # Piezo/tip position within the scan window. Used to align STS and SXM in real space, since sxm image axes are plotted with 0 in the corner.
        self.scanAngle = float(self.sxm.header['scan_angle'])
        
        theta = -math.pi*self.scanAngle/180
        R = np.array([[math.cos(theta)**2,math.sin(theta)**2],
                      [math.sin(theta)**2,math.cos(theta)**2]])
        self.pixelCalibration = np.matmul(R,self.piezoFactor.T)
        
        self.lxy  = self.sxm.header['scan_range']                               # Real length of scan (m)
        self.lxy *= self.pixelCalibration                                       # Piezo calibration
        self.pxy  = self.sxm.header['scan_pixels']                              # Num pixels [x,y]
        self.im_offset  = np.array(self.sxm.header['scan_offset'])
        self.im_offset *= self.pixelCalibration
        self.im_offset -= self.lxy/2                         
        self.size = np.array(self.fig.get_size_inches()*self.fig.dpi)           # Figure size in pixels
        self.dxy = self.lxy/self.pxy                                            # Real size of each pixel on the figure (this is not the resolution of the actual data)
        
        channelList = self._channel(getChannels=True)                           # Get the list of recorded channels for this SXM
        if(not self.plotChannel in channelList):                                # If the current selected channel is not recorded in this file...
            self.plotChannel = channelList[0]                                   # Default the channel to the first one in the list if...
            for c in channelList:
                if('Z' in c):                                                   # ...we don't find one containing Z
                    self.plotChannel = c                                        # If we do find a channel containing Z, default to this one
                    break
        
        self.btn['Channel'].configure(text=self.plotChannel)                    # Set the button text to channel name
        
        rawim = np.array(self.sxm.signals[self.plotChannel]['forward'])         # Raw sxm image. Take the forward scan by convention
        rawim = np.nan_to_num(rawim,nan=np.nanmin(rawim))                       # nan to zero
        self.im = rawim - rawim.min()                                           # Set min value to zero
        
        way  = self.sxm.header['scan_dir']                                      # Scan direction (up/down)
        if way == 'up':                                                         # Invert image according to scan direction
            self.im = np.flipud(self.im)
        
        self.finalim = self.im
        
        if(not load):
            self.vmin = np.min(self.finalim)
            self.vmax = np.max(self.finalim)
            
            self.tiltPlane = np.array([[1/2, 0, 0],[0,1,0],[1,1,0]])            # Init tilt correct
        
        self.planeFitActive = False                                             # Init plane fit (need to move this into if(not load) statement after putting this in the g80 file)
        
        self.extent = (0, self.lxy[0], 0, self.lxy[1])                          # Real-space image boundaries
        
        self.init = True
        
        ## Re-initialise panels that need it
        self.linePanel.init()
        # self.aimlPanel.init()
        self.gridPanel.init()
        
        self.updateHelpLabel("")
        
        self.update(upd=[-1])
    ###########################################################################
    # Inset
    ###########################################################################
    def addInset(self,fig):
        self.moveInsetBind = self.canvas.get_tk_widget().bind('<Double-Button-1>', self._moveInsetBind) # Bind the double click so we can move around the inset
        
        buf = io.BytesIO()
        pickle.dump(fig, buf)                                                   # pickle dump the figure
        buf.seek(0)                                                             # go to byte 0
        insetFig = pickle.load(buf)                                             # reload the figure into a temporary figure
        insetAx = insetFig.axes[0]                                              # Take the axes
        insetAx.remove()                                                        # Remove from the temporary figure
        insetAx.figure = self.fig                                               # Point the axes to the main figure
        if(len(self.fig.axes) > 1):                                             # Remove the previous inset from the sxm figure if there was one
            self.fig.axes[1].remove()
        
        insetAx.set_title("")                                                   # Get rid of the title
        
        self.fig.axes.append(insetAx)                                           # Append the axis to the sxm figure
        self.fig.add_axes(insetAx)
        
        plt.close(insetFig)
            
        self.fig.axes[1].set_position(self.insetPos)                            # Adjust the size and position of the inset
        
        self.setInsetCmap()
        
    def _moveInsetBind(self,event):
        if(self.bound): return
        self.bound = True
        self.lcInsetBind     = self.canvas.get_tk_widget().bind('<Button-1>', self._moveInsetUnbind)
        self.motionInsetBind = self.canvas.get_tk_widget().bind('<Motion>', self._moveInsetPos)
        self.moveInsetActive = True
        
    def _moveInsetUnbind(self,event=[]):
        self.canvas.get_tk_widget().unbind('<Button-1>', self.lcInsetBind)
        self.canvas.get_tk_widget().unbind('<Motion>', self.motionInsetBind)
        self.moveInsetActive = False
        self.bound = False
      
    def _moveInsetPos(self, event):
        size = self.fig.get_size_inches()*self.fig.dpi                          # size in pixels
        x = event.x
        y = size[1] - event.y
        X = super()._getX(x)/self.lxy[0]                                        # Position from 0 to 1 instead of real coordinate
        Y = super()._getY(y)/self.lxy[1] - 0.05                                 # Position from 0 to 1 instead of real coordinate
        
        self.insetPos[0:2] = [X, Y]                                             # Change the position to x,y
        self.fig.axes[1].set_position(self.insetPos)                            # Don't redraw the figure, just reposition it
        self.update(upd=[])                                                     # Update the canvas
    
    def _insetCmap(self):
        self.insetCmap += 1                                                     # Cycle through the colours
        if(self.insetCmap == len(self.insetColours)): self.insetCmap = 0
        self.setInsetCmap()
        self.update(upd=[])
    
    def setInsetCmap(self):
        # Also see Temporary Styling
        # https://matplotlib.org/1.5.3/users/style_sheets.html
        if(len(self.fig.axes) > 1):
            self.fig.axes[1].spines['bottom'].set_color(self.insetColours[self.insetCmap])
            self.fig.axes[1].spines['top'].set_color(self.insetColours[self.insetCmap])
            self.fig.axes[1].spines['right'].set_color(self.insetColours[self.insetCmap])
            self.fig.axes[1].spines['left'].set_color(self.insetColours[self.insetCmap])
            self.fig.axes[1].tick_params(axis='x', colors=self.insetColours[self.insetCmap])
            self.fig.axes[1].tick_params(axis='y', colors=self.insetColours[self.insetCmap])
            self.fig.axes[1].yaxis.label.set_color(self.insetColours[self.insetCmap])
            self.fig.axes[1].xaxis.label.set_color(self.insetColours[self.insetCmap])
        
    def _removeInset(self):
        if(len(self.fig.axes) > 1):                                             # If there is a secondary set of axes on the sxm figure, remove it
            self.fig.axes[1].remove()
            self.update(upd=[])
            self.canvas.get_tk_widget().unbind('<Double-Button-1>', self.moveInsetBind) # Unbind the double click
    ###########################################################################
    # Colourmap Saturation
    ###########################################################################
    def _mouseWheel(self,event):
        if(self.moveInsetActive):
            self.insetPos[2:4] += np.array([event.delta,event.delta])/12000     # Change the inset size, leaving position as is
            self.fig.axes[1].set_position(self.insetPos)                        # Don't redraw the figure, just resize it
            self.update(upd=[])                                                 # Update the canvas
            return
        
        if(self.shiftL):     self._vmin(event.delta)
        if(not self.shiftL): self._vmax(event.delta)
        self.update(upd=[0])
        
    def _vmin(self,delta):
        std = np.std(self.finalim)                                              # This block adjusts the colour map saturation boundaries based on...
        med = np.median(self.finalim)
        
        sxmStdMin = -(self.vmin - med)/std - delta/2400
        
        vmin = med - sxmStdMin*std                                              # sxmStdMin and max are set by mouse scroll wheel on sxm canvas.
        
        if(vmin < np.min(self.finalim)): vmin = np.min(self.finalim)            # Cap vmin
        
        if(vmin > self.vmax): vmin = self.vmax - abs(delta/2400)
        
        self.vmin = vmin
        
    def _vmax(self,delta):
        std = np.std(self.finalim)                                              # This block adjusts the colour map saturation boundaries based on...
        med = np.median(self.finalim)
        
        sxmStdMax =  (self.vmax - med)/std - delta/2400
        
        vmax = med + sxmStdMax*std
        
        if(vmax > 2*np.max(self.finalim)): vmax = 2*np.max(self.finalim)        # Cap vmax
        
        if(vmax < self.vmin): vmax = self.vmin + abs(delta/2400)
        
        self.vmax = vmax
    ###########################################################################
    # Drawing Cursor for 1D Profiles
    ###########################################################################
    def _drawCursor(self):
        if(not self.linePanel.active and not self.linePanel.imprint): return
        
        if(self.linePanel.plotMode == 1):
            linewidth = 2; linestyle = 'dashed'
            for idx,cPos in enumerate(self.linePanel.cPos):
                segInfo  = self.linePanel.segInfo[idx]
                
                x = cPos[:,0]
                y = cPos[:,1]
                c = self.mplibColours[segInfo[3]]                               # Get the default matplotlib colour for this line so it matches the colour on profile panel
                self.ax.plot(x*self.extent[1],y*self.extent[3],c=c,             # Horizontal 1D line through current cursor position
                                    linewidth=linewidth,linestyle=linestyle)
            
                showInfo = segInfo[2]
                if(not showInfo): continue
                
                length = segInfo[0]
                angle  = segInfo[1]
                self.ax.annotate("{:.2f} nm, {:.1f}$^\circ$".format(length,angle),
                                    xy=(cPos[0]*self.lxy),fontsize=13,color=c)
            
            # props = {'ha': 'center', 'va': 'center'}
            # self.ax.text(*(cPos[0]*self.lxy),
            #              "{:.2f} nm, {:.1f}$^\circ$".format(segInfo[0],segInfo[1]),
            #              props,rotation=segInfo[1],fontsize=10,color='black')
            
        if(self.linePanel.plotMode == 0):
            idx = self.linePanel.activeCursor[1]
            linewidth = [0.5,0.75][self.linePanel.activeCursor[0] > -1]             # Cursor is thin when set and thick when placing
            linestyle = ['dashed','solid'][self.linePanel.activeCursor[0] > -1]     # Cursor is dotted when set and solid when placing
            self.ax.axhline(y=self.linePanel.cPos[idx][0][1]*self.extent[3],c='b',    # Horizontal 1D line through current cursor position
                               linewidth=linewidth,linestyle=linestyle)     
            self.ax.axvline(x=self.linePanel.cPos[idx][0][0]*self.extent[1],c='r',    # Vertical 1D line through current cursor position
                               linewidth=linewidth,linestyle=linestyle)
        
    def cursorBind(self):
        if(self.bound): return
        self.bound = True
        self.lcCursorBind     = self.canvas.get_tk_widget().bind('<Button-1>', self._setCursor)
        self.motionCursorBind = self.canvas.get_tk_widget().bind('<Motion>', self._motionCursor)
        
        self.updateHelpLabel("Left click the image to position the cursor")
        
    def _cursorUnbind(self):
        self.canvas.get_tk_widget().unbind('<Button-1>', self.lcCursorBind)
        self.canvas.get_tk_widget().unbind('<Motion>', self.motionCursorBind)
        self.bound = False
        
        self.updateHelpLabel("")
        
    def _setCursor(self, event):
        size = self.fig.get_size_inches()*self.fig.dpi                          # size in pixels
        x = event.x
        y = size[1] - event.y
        X = super()._getX(x)
        Y = super()._getY(y)
        
        cPos = np.array([X,Y])/self.lxy
        self.linePanel.setCursor(cPos, True)
        
        self._cursorUnbind()
        self.update(upd=[0,1])
    
    def _motionCursor(self, event):
        size = self.fig.get_size_inches()*self.fig.dpi                          # size in pixels
        x = event.x
        y = size[1] - event.y
        X = super()._getX(x)
        Y = super()._getY(y)
        
        cPos = np.array([X,Y])/self.lxy
        cPos[cPos > 1] = 1
        cPos[cPos < 0] = 0
        self.linePanel.setCursor(cPos)
        
        self.update(upd=[0,1])
    ###########################################################################
    # Drawing STS Markers
    ###########################################################################
    def _markSTS(self):                                                         # Colour-coded STS locations on sxm image
        if(self.stsPanel.imprint or self.stsPanel.active):
            for i in self.stsPanel.stsPos:                                      # Do a for loop to colour each plot.
                self.ax.plot(i[0],i[1],'.',markersize=12)
            
            customSTSPos = self.stsPanel.customSTSPos
            if(len(customSTSPos)):
                self.ax.plot(customSTSPos[0],customSTSPos[1],'x',markersize=12)
                
    def customSTSBind(self):
        if(self.bound): return
        self.bound = True
        self.lcMarkSTSBind     = self.canvas.get_tk_widget().bind('<Button-1>', self._setMarkSTS)
        self.rcMarkSTSBind     = self.canvas.get_tk_widget().bind('<Button-3>', self._cancelMarkSTS)
        self.motionMarkSTSBind = self.canvas.get_tk_widget().bind('<Motion>',   self._motionMarkSTS)
        
        self.updateHelpLabel("Place the STS marker on the SXM image\n"
                             + "Left click to place\mRight click to cancel")
        
    def _customSTSUnbind(self):
        self.canvas.get_tk_widget().unbind('<Button-1>', self.lcMarkSTSBind)
        self.canvas.get_tk_widget().unbind('<Button-3>', self.rcMarkSTSBind)
        self.canvas.get_tk_widget().unbind('<Motion>',   self.motionMarkSTSBind)
        self.bound = False
        
        self.updateHelpLabel("")
        
    def _cancelMarkSTS(self,event):
        self._customSTSUnbind()
        self.stsPanel.cancelMarkSTS()
        self.update(upd=[0])
        
    def _setMarkSTS(self,event):
        size = self.fig.get_size_inches()*self.fig.dpi                          # size in pixels
        x = event.x
        y = size[1] - event.y
        X = super()._getX(x)
        Y = super()._getY(y)
        
        stsPos = [X,Y]
        self.stsPanel.setMarkSTS(stsPos,True)
        
        self._customSTSUnbind()
        self.update(upd=[0])
    
    def _motionMarkSTS(self,event):
        size = self.fig.get_size_inches()*self.fig.dpi                          # size in pixels
        x = event.x
        y = size[1] - event.y
        X = super()._getX(x)
        Y = super()._getY(y)
        
        stsPos = [X,Y]
        self.stsPanel.setMarkSTS(stsPos)
        self.update(upd=[0])
    ###########################################################################
    # Placing Molecules
    ###########################################################################
    def placeMolecule(self,filename):
        if(not filename):
            filename = self._browseFile()
        self.curMol = filename
        self._placeMoleculeBind()
        
    def _placeMoleculeBind(self):
        if(self.bound): return
        self.bound = True
        self.canvas.get_tk_widget().focus_set()
        self.lcMolBind       = self.canvas.get_tk_widget().bind('<Button-1>', self._setMolecule)
        self.rcMolBind       = self.canvas.get_tk_widget().bind('<Button-3>', self._cancelMolecule)
        self.motionMolBind   = self.canvas.get_tk_widget().bind('<Motion>',   self._moveMolecule)
        self.rotMolUpBind    = self.canvas.get_tk_widget().bind('<Up>',       self._rotUpMolecule)
        self.rotMolDownBind  = self.canvas.get_tk_widget().bind('<Down>',     self._rotDownMolecule)
        self.rotMolRightBind = self.canvas.get_tk_widget().bind('<Right>',    self._rotRightMolecule)
        self.rotMolLeftBind  = self.canvas.get_tk_widget().bind('<Left>',     self._rotLeftMolecule)
        
        self.updateHelpLabel("Place atoms on the image\n"
                             + "LR arrows course rotate around z, UD arrows fine rotate around z\n"
                             + "Shift active: LR arrows rotate around x, UD arrows rotate around y\n"
                             + "Left click to confirm, right click to cancel")
        
    def _placeMoleculeUnbind(self):
        self.canvas.get_tk_widget().unbind('<Button-1>',    self.lcMolBind)
        self.canvas.get_tk_widget().unbind('<Button-3>',    self.rcMolBind)
        self.canvas.get_tk_widget().unbind('<Motion>',      self.motionMolBind)
        self.canvas.get_tk_widget().unbind('<Up>',          self.rotMolUpBind)
        self.canvas.get_tk_widget().unbind('<Down>',        self.rotMolDownBind)
        self.canvas.get_tk_widget().unbind('<Right>',       self.rotMolRightBind)
        self.canvas.get_tk_widget().unbind('<Left>',        self.rotMolLeftBind)
        self.curMol = []
        self.curAtoms = []
        self.bound = False
        
        self.updateHelpLabel("")
    
    def _rotUpMolecule(self,event):
        if(self.shiftL):
            self.curRotX += self.fineRot*2
        else:
            self.curRot += self.fineRot
        self._moveMolecule(event)
        self.update(upd=[0])
    
    def _rotDownMolecule(self,event):
        if(self.shiftL):
            self.curRotX -= self.fineRot*2
        else:
            self.curRot -= self.fineRot
        self._moveMolecule(event)
        self.update(upd=[0])
            
    def _rotRightMolecule(self,event):
        if(self.shiftL):
            self.curRotY -= self.fineRot*2
        else:
            self.curRot -= self.courseRot
        self._moveMolecule(event)
        self.update(upd=[0])
    
    def _rotLeftMolecule(self,event):
        if(self.shiftL):
            self.curRotY += self.fineRot*2
        else:
            self.curRot += self.courseRot
        self._moveMolecule(event)
        self.update(upd=[0])
            
    def _setMolecule(self,event):
        size = self.fig.get_size_inches()*self.fig.dpi                          # size in pixels
        x = event.x
        y = size[1] - event.y
        X = super()._getX(x)
        Y = super()._getY(y)
        
        self.molFiles.append(self.curMol)                                       # Store this molecule's filename on click
        self.molPos.append([X,Y,0])                                             # Store this molecule's position on click
        self.molRot.append(self.curRot)                                         # Store this molecule's rotation on click
        self.molRotX.append(self.curRotX)                                       # Store this molecule's rotation on click
        self.molRotY.append(self.curRotY)                                       # Store this molecule's rotation on click
        
        atoms = read(self.curMol)
        atoms.positions *= 1e-10                                                # Convert to angstroms
        if(self.curRot != 0):                                                   # Doesn't like rotating zero for some reason
            atoms.rotate(self.curRot, 'z')                                      # Rotate the molecule
        if(self.curRotX != 0):                                                  # Doesn't like rotating zero for some reason
            atoms.rotate(self.curRotX, 'x')                                     # Rotate the molecule
        if(self.curRotY != 0):                                                  # Doesn't like rotating zero for some reason
            atoms.rotate(self.curRotY, 'y')                                     # Rotate the molecule
        atoms.positions += np.array([X,Y,0])                                    # Molecule goes where the click happened
        
        self.atoms.append(atoms)
        
        self._placeMoleculeUnbind()
        
        self.update(upd=[0])
    
    def _cancelMolecule(self,event):
        self._placeMoleculeUnbind()
        self.update(upd=[0])
        
    def _moveMolecule(self,event):
        size = self.fig.get_size_inches()*self.fig.dpi                          # size in pixels
        x = event.x
        y = size[1] - event.y
        X = super()._getX(x)
        Y = super()._getY(y)
        
        self.curAtoms = read(self.curMol)
        if(self.curRot != 0):                                                   # Doesn't like rotating zero for some reason
            self.curAtoms.rotate(self.curRot, 'z')                              # Rotate the molecule
        if(self.curRotX != 0):                                                  # Doesn't like rotating zero for some reason
            self.curAtoms.rotate(self.curRotX, 'x')                             # Rotate the molecule
        if(self.curRotY != 0):                                                  # Doesn't like rotating zero for some reason
            self.curAtoms.rotate(self.curRotY, 'y')                             # Rotate the molecule
        self.curAtoms.positions *= 1e-10                                        # Convert to angstroms
        self.curAtoms.positions += np.array([X,Y,0])                            # Molecule goes where the click happened
        
        self.update(upd=[0])
        
    def _drawMolecules(self):                                                   # Draws all the placed/set molecules on the sxm image
        if(len(self.atoms)):
            for m in self.atoms:
                plot_atoms(m, self.ax, radii=0.5e-10)                           # Plot this molecule
                
        if(self.curMol):
            plot_atoms(self.curAtoms, self.ax, radii=0.5e-10)                   # Plot this molecule
    
    def _loadMolecules(self):
        self.atoms = []
        for i in range(len(self.molFiles)):                                     # Draw all the saved molecules
            atoms = read(self.molFiles[i])                                      # Read in xyz file
            atoms.positions *= 1e-10                                            # Convert to angstroms
            if(self.molRot[i] != 0):                                            # Doesn't like rotating zero for some reason
                atoms.rotate(self.molRot[i], 'z')                               # Rotate the molecule
            if(self.molRotX[i] != 0):                                           # Doesn't like rotating zero for some reason
                atoms.rotate(self.molRotX[i], 'x')                              # Rotate the molecule
            if(self.molRotY[i] != 0):                                           # Doesn't like rotating zero for some reason
                atoms.rotate(self.molRotY[i], 'y')                              # Rotate the molecule
            atoms.positions += np.array(self.molPos[i])                         # Molecule goes where the click happened
            
            self.atoms.append(atoms)
    
    def undoMolecule(self):
        self.atoms    = self.atoms[:-1]
        self.molFiles = self.molFiles[:-1]
        self.molRot   = self.molRot[:-1]
        self.molRotX  = self.molRotX[:-1]
        self.molRotY  = self.molRotY[:-1]
        self.molPos   = self.molPos[:-1]
        self.update(upd=[0])
        
    ###########################################################################
    # STS Grid Box
    ###########################################################################
    def _drawGridBox(self):
        if(self.gridPanel.imprint or self.gridPanel.active):
            if(not self.gridPanel.gridData): return
            bbox = self.gridPanel.bbox
            ox,oy = bbox['origin']
            angle = bbox['angle']
            r = patches.Rectangle(*bbox['bbox'],linewidth=2,edgecolor='black',facecolor='none')
            transform  = matplotlib.transforms.Affine2D().rotate_deg_around(ox,oy,angle)
            transform += self.ax.transData
            r.set_transform(transform)
            self.ax.add_patch(r)
            
    ###########################################################################
    # Tilt Methods
    ###########################################################################
    def tilt(self):                                                             # This will eventually take tiltFactor as a user input
        if(self.tiltActive):
            self._setTilt()
        else:
            self._tiltBind()
            
    def _tiltPlane(self):
        # Calculate plane through three points. Taken from...
        # https://kitchingroup.cheme.cmu.edu/blog/2015/01/18/Equation-of-a-plane-through-three-points/
        _range = abs(self.vmax - self.vmin)*self.lxy[0]
        p1 = self.curTilt[0]; p2 = self.curTilt[1]; p3 = self.curTilt[2];       # Three points that define the plane
        v1 = p3 - p1;   v2 = p2 - p1;                                           # Vectors in plane
        cp = np.cross(v1,v2)                                                    # Cross product
        a,b,c = cp                                                              # Coefs
        d = np.dot(cp,p3)
        
        x = np.linspace(0, len(self.im[0]),len(self.im[0]))
        y = np.linspace(0, len(self.im),len(self.im))
        X,Y = np.meshgrid(x,y)
        Z = np.array((d - a * X - b * Y) / c)                                   # Plane defined by points curTilt
        return Z*_range                                                         # Scale Factor
    
    def _tiltBind(self):
        if(self.bound): return
        self.bound = True
        self.canvas.get_tk_widget().focus_set()
        self.lcTiltBind     = self.canvas.get_tk_widget().bind('<Button-1>',self._setTilt)     # Left click calls self._setTilt
        self.rcTiltBind     = self.canvas.get_tk_widget().bind('<Button-3>',self._cancelTilt)  # Right click bind
        self.upTiltBind     = self.canvas.get_tk_widget().bind('<Up>',      self._upTilt)      # Up arrow key bind
        self.downTiltBind   = self.canvas.get_tk_widget().bind('<Down>',    self._downTilt)    # Down arrow key bind
        self.rightTiltBind  = self.canvas.get_tk_widget().bind('<Right>',   self._rightTilt)   # Right arrow key bind
        self.leftTiltBind   = self.canvas.get_tk_widget().bind('<Left>',    self._leftTilt)    # Left arrow key bind
        self.tiltActive = True                                                  
        self.btn['Tilt'].configure(bg='Red')
        
        self.updateHelpLabel("Use the arrow keys to tilt the image. Recommend openning the Profiles Panel\n"
                             + "Shift inactive: fine tilt\nShift active: course tilt\n"
                             + "Left click to confirm, right click to cancel")
        
    def _tiltUnbind(self):
        self.canvas.get_tk_widget().unbind('<Button-1>',self.lcTiltBind)
        self.canvas.get_tk_widget().unbind('<Button-3>',self.rcTiltBind)
        self.canvas.get_tk_widget().unbind('<Up>',      self.upTiltBind)
        self.canvas.get_tk_widget().unbind('<Down>',    self.downTiltBind)
        self.canvas.get_tk_widget().unbind('<Right>',   self.rightTiltBind)
        self.canvas.get_tk_widget().unbind('<Left>',    self.leftTiltBind)
        self.tiltActive = False
        self.btn['Tilt'].configure(bg='SystemButtonFace')
        self.bound = False
        
        self.updateHelpLabel("")
    
    def _upTilt(self,event=[]):
        self.curTilt[0][2] -= self.tiltFactor
        self.update(upd=[0,1])
        
    def _downTilt(self,event=[]):
        self.curTilt[0][2] += self.tiltFactor
        self.update(upd=[0,1])
    
    def _rightTilt(self,event=[]):
        self.curTilt[1][2] += self.tiltFactor
        self.update(upd=[0,1])
        
    def _leftTilt(self,event=[]):
        self.curTilt[1][2] -= self.tiltFactor
        self.update(upd=[0,1])
    
    def _setTilt(self,event=[]):
        self.tiltPlane = np.copy(self.curTilt)
        self._tiltUnbind()
        self.update()
        
    def _cancelTilt(self,event=[]):
        self.curTilt = np.copy(self.tiltPlane)
        self._tiltUnbind()
        self.update()
    ###########################################################################
    # Plane fit
    ###########################################################################
    def planeFit(self):
        self.planeFitBind()
        
    def planeFitBind(self):
        if(self.bound): return
        self.bound = True
        self.lcPFBind       = self.canvas.get_tk_widget().bind('<Button-1>', self.setPlaneFitArea)
        self.rcPFBind       = self.canvas.get_tk_widget().bind('<Button-3>', self.cancelPlaneFit)
        self.motionPFBind   = self.canvas.get_tk_widget().bind('<Motion>',   self.placePlaneFitArea)
        self.curPlaneFitArea = [[0,0],[0,0]]
        self.planeFitCursor = 0
        
        self.updateHelpLabel("Click at two locations on the image to highlight region to plane fit.")
        
    def planeFitUnbind(self):
        self.canvas.get_tk_widget().unbind('<Button-1>',    self.lcPFBind)
        self.canvas.get_tk_widget().unbind('<Button-3>',    self.rcPFBind)
        self.canvas.get_tk_widget().unbind('<Motion>',      self.motionPFBind)
        self.planeFitCursor = -1
        self.bound = False
        
        self.updateHelpLabel("")
        
    def placePlaneFitArea(self,event):
        size = self.fig.get_size_inches()*self.fig.dpi                          # size in pixels
        x = event.x
        y = size[1] - event.y
        X = super()._getX(x)
        Y = super()._getY(y)
        
        self.curPlaneFitArea[self.planeFitCursor] = np.array([X,Y])
        self.update(upd=[0])
    
    def setPlaneFitArea(self,event):
        self.planeFitCursor += 1
        if(self.planeFitCursor == 2):
            self.planeFitUnbind()
            self.planeFitArea = self.curPlaneFitArea
            self.tiltPlane = np.array([[1/2, 0, 0],[0,1,0],[1,1,0]]);
            
            topx = np.min([self.curPlaneFitArea[0][0],self.curPlaneFitArea[1][0]])
            topy = np.max([self.curPlaneFitArea[0][1],self.curPlaneFitArea[1][1]])
            topy = self.lxy[1] - topy
            topCorner = ((np.array([topx,topy])/self.lxy)*np.flip(self.im.shape)).astype(int)
            
            botx = np.max([self.curPlaneFitArea[0][0],self.curPlaneFitArea[1][0]])
            boty = np.min([self.curPlaneFitArea[0][1],self.curPlaneFitArea[1][1]])
            boty = self.lxy[1] - boty
            botCorner = ((np.array([botx,boty])/self.lxy)*np.flip(self.im.shape)).astype(int)
            
            self.planeFitBox = [topCorner,botCorner]
            
            self.curTilt   = np.array([[1/2, 0, 0],[0,1,0],[1,1,0]]);
            self.tiltPlane = np.array([[1/2, 0, 0],[0,1,0],[1,1,0]]);
            
            self.planeFitActive = True
        self.update()
    
    def _planeFit(self,im):
        if(not self.planeFitActive): return im
        topCorner = self.planeFitBox[0]
        botCorner = self.planeFitBox[1]
        im = napfit.plane_fit_2d(im,region=[topCorner,botCorner])
        # self.vmin, self.vmax = napfit.filter_sigma(self.im)                     # cmap saturation. 3 sigma by default
        return im
    
    def cancelPlaneFit(self,event):
        self.planeFitUnbind()
        self.update()
    
    def _drawFitArea(self):
        if(self.planeFitCursor == -1): return
        self.ax.plot(self.curPlaneFitArea[0][0],self.curPlaneFitArea[0][1],'x',markersize=10,color='white')
        if(self.planeFitCursor == 1):
            w = abs(self.curPlaneFitArea[0][0] - self.curPlaneFitArea[1][0])
            h = abs(self.curPlaneFitArea[0][1] - self.curPlaneFitArea[1][1])
            topx = np.min([self.curPlaneFitArea[0][0],self.curPlaneFitArea[1][0]])
            topy = np.min([self.curPlaneFitArea[0][1],self.curPlaneFitArea[1][1]])
            topCorner = (topx,topy)
            
            rect = patches.Rectangle(topCorner, w, h, linewidth=2, edgecolor='white', facecolor='red',alpha=0.3)
        
            self.ax.add_patch(rect)
            self.ax.plot(self.curPlaneFitArea[1][0],self.curPlaneFitArea[1][1],'x',markersize=10,color='white')
    ###########################################################################
    # Global alternate functions
    ###########################################################################
    def toggleShiftL(self,event=[]):
        self.shiftL = not self.shiftL
        
        # self.tiltFactor = 1e-14*(19*self.shiftL + 1)                            # Here, shift toggles the tiltFactor
        self.tiltFactor = (900*self.shiftL + 100)                               # Here, shift toggles the tiltFactor
        
        btnColour = ['SystemButtonFace','red']                                  # Highlight the shift button red when shift is active
        self.btn['Shift'].configure(bg = btnColour[self.shiftL])
    ###########################################################################
    # Save
    ###########################################################################
    def _exportPNG(self,dpi=0):
        if(not dpi): dpi = self.dpi
        initialfile = self.sxm.header['scan_file'].rsplit('\\',1)[1].rsplit('.')[0]
        
        path = filedialog.asksaveasfilename(title="Save as",initialfile=initialfile + '.png')
        if path == "":
            return None
        if(not path.endswith('.png')): path += '.png'
        self.fig.savefig(path,format='png',dpi=dpi)
                
    def _save(self):
        saveString  = self._buildSaveString()
        saveString += self.linePanel.buildSaveString()
        saveString += self.fftPanel.buildSaveString()
        
        default = self.sxm.header['scan_file'].rsplit('\\',1)[1].rsplit('.')[0] + '.g80'
        path = filedialog.asksaveasfilename(title="Save as",initialfile=default)
        if(not path.endswith('.g80')): path += '.g80'
        with open(path,'w') as f:
            f.write(saveString)
        
    def _buildSaveString(self):
        saveString = "#MainPanel\n"                                             # Header
        
        saveString += self.filename + "\n"
        
        tilt = []
        for tiltPlane in self.tiltPlane:
            tilt.append(','.join("{:.17f}".format(t) for t in tiltPlane))
        saveString += '|'.join(i for i in tilt) + "\n"                          # Line 3: 3 points defining tilt plane
        
        saveString += ("{:.20f}".format(self.vmin) + "," +                      # Line 4: min/max for colour map
                       "{:.20f}".format(self.vmax) + "\n")
        
        saveString += str(self.cmap) + "\n"                                     # Line 5: Colour map index
        
        molPos   = []
        for i in range(len(self.molFiles)):
            molPos.append(','.join(str(pos) for pos in self.molPos[i]))
        saveString += ('|'.join(i for i in self.molFiles) + "\n")               # Line 6: List of molecules to place (.xyz filename)
        saveString += ('|'.join(i for i in molPos) + "\n")                      # Line 7: Position of each molecule
        saveString += ('|'.join(i for i in [str(int)                            # Line 8: Rotation of each molecule in z
                                     for int in self.molRot]) + "\n")
        saveString += ('|'.join(i for i in [str(int)                            # Line 9: Rotation of each molecule in x
                                     for int in self.molRotX]) + "\n")
        saveString += ('|'.join(i for i in [str(int)                            # Line 10: Rotation of each molecule in y
                                     for int in self.molRotY]) + "\n")
                
        return saveString
    ###########################################################################
    # Load
    ###########################################################################
    def _load(self,g80File=[]):
        if(not g80File):
            g80File = self._browseFile()
        if(not g80File.endswith(".g80")):                                       # Needs to be .g80 file
            print("Expecting .g80 file")
            return
        
        self._loadFromFile(g80File)
        self.linePanel.loadFromFile(g80File)
        self.fftPanel.loadFromFile(g80File)
    
    def _loadFromFile(self,g80File):
        headerFound = False
        with open(g80File, 'r') as f:
            line = "begin"
            while(not headerFound and line):
                line = f.readline()[:-1]
                if(line == "#MainPanel"): headerFound = True
            if(not headerFound): print("Could not find header"); return
            
            filename = f.readline()[:-1]
            
            self.tiltPlane = []
            tilt = f.readline()[:-1].split('|')                                 # Line 3: 3 points defining tilt plane
            for i in tilt:
                self.tiltPlane.append([float(t) for t in i.split(',')])
            self.tiltPlane = np.array(self.tiltPlane)
            self.curTilt   = np.copy(self.tiltPlane)
            
            cmapSat = f.readline()[:-1].split(',')                              # Line 4: Colour map saturation (num STDs around median)
            self.vmin = float(cmapSat[0])
            self.vmax = float(cmapSat[1])
            
            cmap = f.readline()[:-1]                                            # Line 5: colour map index
            self.cmap = int(cmap) - 1                                           # Subtract 1 from cmap
            super()._cmap()                                                     # Then cycle back to it using this function which updates cmap button text
            
            self.molFiles = []; self.molPos = []; self.molRot = []
            molFiles = f.readline()[:-1]                                        # Line 6: List of molecules to place (.xyz filename)
            molPos   = f.readline()[:-1].split('|')                             # Line 7: Position of each molecule
            molRot   = f.readline()[:-1].split('|')                             # Line 8: Rotation of each molecule in z
            molRotX  = f.readline()[:-1].split('|')                             # Line 9: Rotation of each molecule in x
            molRotY  = f.readline()[:-1].split('|')                             # Line 10: Rotation of each molecule in y
            if(molFiles):
                self.molFiles = molFiles.split('|')
                for i in molPos:
                    self.molPos.append([float(pos) for pos in i.split(',')])
                self.molRot  = [float(rot) for rot in molRot]
                self.molRotX = [float(rot) for rot in molRotX]
                self.molRotY = [float(rot) for rot in molRotY]
            
            self._loadMolecules()                                               # Stores the molecules in an array of atoms. Otherwise need to repeatedly load in xyz files every canvas update
            
            self.loadSXM(filename=filename,load=True)                           # Reload
    ###########################################################################
    # Misc
    ###########################################################################
    def _toggleCaption(self):
        self.plotCaption = not self.plotCaption
        self.update(upd=[0])
        
    def _channel(self,getChannels=False):
        channels = list(self.sxm.signals.keys())
        if(getChannels): return channels
        idx = channels.index(self.plotChannel) + 1
        if(idx == len(channels)): idx = 0
        self.plotChannel = channels[idx]
        self.btn['Channel'].configure(text=self.plotChannel)
        rawim = np.array(self.sxm.signals[self.plotChannel]['forward'])         # Raw sxm image. Take the forward scan by convention. 
        rawim = np.nan_to_num(rawim,nan=np.nanmin(rawim))                       # nan to zero
        self.im = rawim - rawim.min()                                           # Set min value to zero
        
        way  = self.sxm.header['scan_dir']                                      # Scan direction (up/down)
        if way == 'up':                                                         # Invert image according to scan direction
            self.im = np.flipud(self.im)
        
        self.update()
        self.vmin = np.min(self.finalim)
        self.vmax = np.max(self.finalim)
        self.update(upd=[0])
        
    def _flipScan(self):
        self.im = np.flipud(self.im)
        self.update()
    
    def _rotateScan(self):
        self.im = np.rot90(self.im)
        self.update()
        
    def quit(self):
        self.master.destroy()                                                   # Close the Tkinter GUI
        os._exit(00)                                                            # Easier to restart kernal than to figure out Tkinter memory leaks (think this is a problem with spyder, not tkinter)