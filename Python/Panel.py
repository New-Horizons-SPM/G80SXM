# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 18:56:27 2022

@author: jced0001
"""

import tkinter as tk
import customtkinter as ctk
from tkinter import filedialog
import numpy as np
import math
# import matplotlib
# matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import colors
from matplotlib_scalebar.scalebar import ScaleBar

class Panel():
    pos    = 0                                                                  # This panel's column position. init to zero. gets added to the RHS of the window when created
    active = False; cmap  = 0; shift = False                                    # Use alternate values - bound to the left shift key
    zunit = 1e-12; xunit = 1e-9                                                 # convert xy units to nm and z units to pm
    helpText = ""
    imprint  = False
    def __init__(self, master, width, height, dpi, mainPanel=None):
        self.master = master
        self.width  = width
        self.height = height
        self.dpi    = dpi
        self.mainPanel = mainPanel
        
        # Set up canvas
        self.canvas = FigureCanvasTkAgg(master=master)
        self.canvas.get_tk_widget().configure(width=width, height=height)
        
        # Figure
        self.fig = plt.figure(figsize=(width/dpi,height/dpi),dpi=dpi)
        self.ax  = self.fig.add_subplot(111)
        
        # Global bindings - Bind these keys to all panels
        self.canvas.get_tk_widget().bind('<Shift_L>',self.toggleShiftL)
        
        # Misc init
        self._initCmaps()
    
    def addPlotCaption(self, plotCaption, pos, props=None, fontsize=16, ax=None):
        if(not props):
            props = dict(boxstyle='round',facecolor='white',alpha=0.5)          # Box properties. alpha = transparency
        if(not ax): ax = self.ax
        ax.text(pos[0],pos[1],plotCaption,bbox=props,fontsize=fontsize)
    
    def addPlotScalebar(self, fontsize=20, ax=None):
        # Add a scalebar to the sxm image
        font_properties = {"size": fontsize}
        scalebar = ScaleBar(dx=1, units='m', length_fraction=.5,
        font_properties=font_properties, frameon=False, location=4, color='w')
        
        if(not ax): ax = self.ax
        ax.add_artist(scalebar)
        
    def addButtons(self):
        numCols = 4                                                             # We're displaying 4 buttons per row per panel
        row = 5; col = 0                                                        # Row=5 because rows 0-4 are taken up by the canvas
        self.master.rowconfigure(index=row,minsize=40)
        for btn in self.btn.values():                                           # Now loop through each button and put it in the correct position. 4 per row
            btn.grid(row=row,column=col+self.pos)
            btn.configure(width=120,height=35)
            col += 1
            if(col == numCols):
                col  = 0
                row += 1
                self.master.rowconfigure(index=row,minsize=40)
    def removeButtons(self):                                                    # Remove all the buttons after destroying the canvas
        for btn in self.btn.values():
            btn.grid_forget()

    def special(self):
        pass
    def removeSpecial(self):
        pass
    
    def _helpLabel(self):
        self.helpLabel = tk.Label(self.master,text="",justify=tk.LEFT)
        self.helpLabel.grid(row=12,column=self.pos,columnspan=4,rowspan=4)
    
    def removeHelpLabel(self):
        self.helpLabel.grid_forget()
    
    def updateHelpLabel(self,helpText):
        self.helpText = helpText
        self.helpLabel.configure(text=helpText)
        
    def create(self):                                                           # Displays the panel to the right of the previous one
        if(self.mainPanel):
            if(not self.mainPanel.init): return
        if(self.active): return                                                 # Do nothing if the panel is already active
        self.pos = self.master.grid_size()[0]                                   # Position of the end of the last panel
        self.canvas.get_tk_widget().grid(row=0,column=self.pos, rowspan=4,columnspan=4) # Put this panel after the last one (left to right)
        self.addButtons()                                                       # Display the buttons
        self.special()
        self._helpLabel()
        self.active = True
        self.update()
        if(self.mainPanel): self.mainPanel.update()
    
    def destroy(self):                                                          # Hide this canvas, it's panel is not active anymore
        self.canvas.get_tk_widget().grid_forget()
        self.removeButtons()                                                    # Also hide the buttons
        self.removeSpecial()
        self.removeHelpLabel()
        self.active = False
        self.mainPanel.update()
    
    def addInset(self):
        self.mainPanel.addInset(self.fig)
        self.mainPanel.update()
    
    def _imprint(self):
        self.imprint = not self.imprint
        self.btn['Imprint'].configure(bg=['SystemButtonFace','red'][self.imprint])
        self.mainPanel.update()
        
    def toggleShiftL(self,event=[]):
        self.mainPanel.toggleShiftL()
    
    def _initCmaps(self):
        self.cmaps = {0 : ['viridis',lambda : "viridis"],                       # Colour map for all images. More at https://matplotlib.org/stable/gallery/color/colormap_reference.html
                      1 : ['plasma', lambda : "plasma"],
                      2 : ['inferno',lambda : "inferno"],
                      3 : ['magma',  lambda : "magma"],
                      4 : ['cividis',lambda : "cividis"],
                      5 : ['flame',  lambda : self.customCmap(cmap='flame')],
                      6 : ['Greys',  lambda : "Greys_r"],
                      7 : ['Blues',  lambda : "Blues_r"],
                      8 : ['Purples',lambda : "Purples_r"]
                      }
        
    def _cmap(self):
        self.cmap += 1
        if(self.cmap == len(self.cmaps)):
            self.cmap = 0
        self.btn['cmap'].config(text = self.cmaps[self.cmap][0])
        self.update()
        
    def customCmap(self,cmap):
        nodes=[0, 0.2, 0.5, 0.8, 1]
        col = [[0, 0, 0, 255], [7, 0, 220, 255], [236, 0, 134, 255], [246, 246, 0, 255], [255, 255, 255, 255]]
        col = np.asarray(col) / 256
        return colors.LinearSegmentedColormap.from_list(cmap, list(zip(nodes, col)))
        
    def _getX(self, x):                                                         # Takes in x pixel position and returns the corresponding value on the x axis
        size = self.fig.get_size_inches()*self.fig.dpi # size in pixels
        pos = np.array(self.ax.get_position())
        lim = np.array(self.ax.get_xlim())
        dx = (lim[1] - lim[0])/((pos[1][0] - pos[0][0])*size[0])
        x = x - pos[0][0]*size[0]
        x = x*dx + lim[0]
        if(x < lim[0]): x = lim[0]
        if(x > lim[1]): x = lim[1]
        return x
    
    def _getY(self, y):                                                         # Takes in y pixel position and returns the corresponding value on the y axis
        size = self.fig.get_size_inches()*self.fig.dpi # size in pixels
        pos = np.array(self.ax.get_position())
        lim = np.array(self.ax.get_ylim())
        dy = (lim[1] - lim[0])/((pos[1][1] - pos[0][1])*size[1])
        y = y - pos[0][1]*size[1]
        y = y*dy + lim[0]
        if(y < lim[0]): y = lim[0]
        if(y > lim[1]): y = lim[1]
        return y
    
    def rotate(self,origin, point, angle):
        """
        Taken from:
        https://stackoverflow.com/questions/34372480/rotate-point-about-another-point-in-degrees-python
        Rotate a point counterclockwise by a given angle around a given origin.
    
        The angle should be given in radians.
        """
        ox, oy = origin
        px, py = point
    
        qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
        qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
        return qx, qy
    
    def _browseFile(self):
        path = filedialog.askopenfilename(title='Select File')                  # shows dialog box and return the path\filename
        return path
    
    def _browseFolder(self):
        path = filedialog.askdirectory(title='Select Folder')                   # shows dialog box and return the path
        return path