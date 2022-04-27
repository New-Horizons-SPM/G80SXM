# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 18:09:36 2022

@author: jced0001
"""

from Panel import Panel
import tkinter as tk
class FilterPanel(Panel):
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
            "Close": tk.Button(self.master, text="Close",   command=self.destroy)
            }
    ###########################################################################
    # Update and Plotting (WIP)
    ###########################################################################
    def update(self):
        if(not self.mainPanel.init): return
        
        self.ax.cla()                                                           # Clear the axis
        
        self.ax.set_position([0, 0, 1, 1])                                      # Make it take up the whole canvas
        self.ax.axis('off')                                                     # Hide the axis
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()                                                      # Redraw the canvas with the updated figure
    