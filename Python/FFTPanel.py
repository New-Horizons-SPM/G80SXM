# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 18:09:36 2022

@author: jced0001
"""

from Panel import Panel
import tkinter as tk
import numpy as np
class FFTPanel(Panel):
    fftZoom = 1;
    fftLabel = []; 
    fftDataActive = False
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width, height, dpi, mainPanel):
        super().__init__(master, width, height, dpi, mainPanel=mainPanel)
        self.canvas.get_tk_widget().bind("<MouseWheel>", self._zoom)            # Used for zooming on FFT canvas
        self.buttons()
    ###########################################################################
    # Panel
    ###########################################################################
    def buttons(self):
        self.btn = {
            "cmap":  tk.Button(self.master, text="viridis",     command=super()._cmap), # Button to cycle through colour maps
            "Place": tk.Button(self.master, text="Data Point",  command=self._fftDataPointBind),
            "Undo":  tk.Button(self.master, text="Undo",        command=self._undo),
            "Reset": tk.Button(self.master, text="Reset",       command=self._reset),
            "Inset": tk.Button(self.master, text="Inset",       command=super().addInset),
            "Close": tk.Button(self.master, text="Close",       command=self.destroy)
            }
    ###########################################################################
    # Update and Plotting
    ###########################################################################
    def update(self):
        if(not self.mainPanel.init): return
        
        self.ax.cla()                                                           # Clear the axis
        self._takeFFT()
        self._placeFftLabel()
        
        self.ax.set_position([0, 0, 1, 1])                                      # Make it take up the whole canvas
        self.ax.axis('off')                                                     # Hide the axis
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()                                                      # Redraw the canvas with the updated figure
    
    def _takeFFT(self):
        im  = np.copy(self.mainPanel.finalim)
        dxy = np.copy(self.mainPanel.dxy)
        fx = self.xunit/(2*dxy[0])                                              # xlim in Fourier domain
        fy = self.xunit/(2*dxy[1])                                              # ylim in Fourier domain
        
        self.fftExtent = [-fx,fx,-fy,fy]                                        # Extent in Fourier domain
        fft = np.log(np.abs(np.fft.fftshift(np.fft.fft2(im)))**2)               # Take the log of |FFT|^2 and shift low freqs to origin
        
        cmap = self.cmaps[self.cmap][1]
        self.ax.imshow(fft,extent=self.fftExtent,cmap=cmap())                   # Show the image
        self.ax.set_xlim((-fx/self.fftZoom,fx/self.fftZoom))                    # Zoom in (fftZoom set by _canvas2Motion)
        self.ax.set_ylim((-fy/self.fftZoom,fy/self.fftZoom))
    
    def _placeFftLabel(self):                                                   # Place labels at selected points (magnitude of wavelength)
        for i in self.fftLabel:                                                 # fftLabel contains [1/lambda_x,1/lambda_y]
            wavelength = 1/np.sqrt(sum(i**2))
            self.ax.plot(i[0], i[1],'k',marker='.',markersize=5)
            self.ax.annotate("{:.3f} nm".format(wavelength), xy=(i[0], i[1]),fontsize=8,color='black')
            
        if(self.fftDataActive):                                                 # Also take care of the point we're currently labeling with the mouse
            wavelength = 1/np.sqrt(sum(self.fftMotionLabel**2))
            self.ax.plot(self.fftMotionLabel[0], self.fftMotionLabel[1],'k',marker='.',markersize=5)
            self.ax.annotate("{:.3f}".format(wavelength), xy=(self.fftMotionLabel[0], self.fftMotionLabel[1]),fontsize=8,color='black')
    ###########################################################################
    # Placing Labels
    ###########################################################################
    def _fftDataPointBind(self):
        if(self.fftDataActive):
            self._fftDataPointUnbind()
            return
        self.motionBind    = self.canvas.get_tk_widget().bind('<Motion>',   self._motionDataPoint)
        self.lcBind        = self.canvas.get_tk_widget().bind('<Button-1>', self._setDataPoint)
        self.rcBind        = self.canvas.get_tk_widget().bind('<Button-3>', self._cancelDataPoint)
        self.fftDataActive = True
        
    def _fftDataPointUnbind(self):
        self.canvas.get_tk_widget().unbind('<Motion>', self.motionBind)
        self.canvas.get_tk_widget().unbind('<Button-1>', self.lcBind)
        self.canvas.get_tk_widget().unbind('<Button-3>', self.rcBind)
        self.fftDataActive = False
    
    def _motionDataPoint(self,event):
        size = self.fig.get_size_inches()*self.fig.dpi                          # size in pixels
        x = event.x
        y = size[1] - event.y
        X = super()._getX(x)
        Y = super()._getY(y)
        
        self.fftMotionLabel = np.array([X,Y])
        self.update()
        
    def _setDataPoint(self,event):
        size = self.fig.get_size_inches()*self.fig.dpi                          # size in pixels
        x = event.x
        y = size[1] - event.y
        X = super()._getX(x)
        Y = super()._getY(y)
        
        self.fftLabel.append(np.array([X,Y]))
        self._fftDataPointUnbind()
        self.update()
        
    def _cancelDataPoint(self,event):
        self._fftDataPointUnbind()
        self.update()
    ###########################################################################
    # Misc Button Functions
    ###########################################################################
    def _undo(self):
        self.fftLabel = self.fftLabel[:-1]
        self.update()
        
    def _zoom(self,event):
        self.fftZoom += self.fftZoom*event.delta/1200
        if(self.fftZoom < 1):  self.fftZoom = 1
        if(self.fftZoom > 10): self.fftZoom = 10
        self.update()
        
    def _reset(self):
        self.fftLabel = []
        self.update()
    ###########################################################################
    # Save
    ###########################################################################
    def buildSaveString(self):
        saveString = "#FFTPanel\n"                                              # Line 1: Header
        
        saveString += str(self.cmap) + "\n"                                     # Line 2: Colour map
        
        saveString += str(self.fftZoom) + "\n"                                  # Line 3: Zoom
        
        saveString += str(len(self.fftLabel)) + "\n"                            # Line 4: Number of labels to follow
        for label in self.fftLabel:
            saveString += ','.join("{:.5f}".format(pos)                         # Line 5...: FFT Labels
                             for pos in label) + "\n"
        
        return saveString
    ###########################################################################
    # Load
    ###########################################################################
    def loadFromFile(self,g80File):
        headerFound = False
        with open(g80File, 'r') as f:
            line = "begin"
            while(not headerFound and line):
                line = f.readline()[:-1]
                if(line == "#FFTPanel"): headerFound = True
            if(not headerFound): print("Missing #FFTPanel"); return
            
            cmap = f.readline()[:-1]                                            # Line 2: colour map index
            self.cmap = int(cmap) - 1                                           # Subtract 1 from cmap
            super()._cmap()                                                     # Then cycle back to it using this function which updates cmap button text
            
            self.fftZoom = float(f.readline()[:-1])                             # Line 3: Plot Mode
            
            numLabels = int(f.readline()[:-1])                                  # Line 4: Number of labels to follow
            
            self.fftLabel = []                                                  # Line 5...: label
            for l in range(numLabels):
                label = np.array(f.readline()[:-1].split(','))
                self.fftLabel.append(label.astype(float))