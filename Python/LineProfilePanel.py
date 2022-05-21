# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 18:03:13 2022

@author: jced0001
"""

from Panel import Panel
import tkinter as tk
import numpy as np
import math
class LineProfilePanel(Panel):
    plotModes = ["XY", "P1P2"]; plotMode = 0                                    # XY plots 1D profiles through a crosshair at cPos[0].
    cPos = [np.array([[0.5,0.5],[0.75,0.75]])]                                  # P1P2 plots a 1D profile through the line segment from cPos[0] to cPos[1]
    activeCursor = np.array([-1,0])                                             # [0] -1=not placing any cursors atm. 0=Placing P1 (used in XY and P1P2 modes) 1=Placing P2 (used in only P2 mode when drawing a line b/w P1 and P2)
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
            "Mode"          : tk.Button(self.master, text="Mode: xy",   command=self.toggleMode),
            "Cursor 1"      : tk.Button(self.master, text="Cursor 1",   command=lambda:self.cursor(0)),
            "Cursor 2"      : tk.Button(self.master, text="Cursor 2",   command=lambda:self.cursor(1)),
            "Next Cursor"   : tk.Button(self.master, text="Next Cursor",command=self.nextCursor),
            "Add Cursor"    : tk.Button(self.master, text="Add Cursor", command=self.addCursor),
            "Rem Cursor"    : tk.Button(self.master, text="Rem Cursor", command=self.remCursor),
            "Inset"         : tk.Button(self.master, text="Inset",      command=super().addInset),
            "Imprint"       : tk.Button(self.master, text="Imprint",    command=super()._imprint),
            "Close"         : tk.Button(self.master, text="Close",      command=self.destroy)
            }
    
    ###########################################################################
    # Update and Plotting
    ###########################################################################
    def update(self):
        if(not self.mainPanel.init): return
        
        self.ax.cla()                                                           # Clear the axis
        if(not self.plotMode):  self.plot1D();
        else:                   self.plotP1P2()
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()                                                      # Redraw the canvas with the updated figure
    
    def plot1D(self):                                                           # 1D x and y plots through current cursor pos
        im  = np.copy(self.mainPanel.finalim)
        dxy = np.copy(self.mainPanel.dxy)
        sx = len(im[0]); sy = len(im)                                           # size of the image in pixels
        
        idx = self.activeCursor[1]
        zx = im[int(sy*(1-self.cPos[idx][0][1])),:]/self.zunit                               # Raw 1D horizontal cut
        zy = im[:,int(sx*self.cPos[idx][0][0])]/self.zunit                                   # Raw 1D vertical cut
        offset =  0*1.05*np.max(zx)                                             # Offset the 1D cuts so easier to see. (currently set to zero)
        
        xx = np.linspace(0,sx*dxy[0],sx)/self.xunit                             # Pixel size in x and y might be different (lines vs pixels)
        yy = np.linspace(0,sy*dxy[1],sy)/self.xunit                             # so different x axis for horizontal and vertical cuts
        
        # Plotting
        self.ax.plot(xx,zx,c='b'); self.ax.plot(yy,zy + offset,c='r')           #
        self.ax.set_xlabel('x, y (nm)'); self.ax.set_ylabel('z (pm)')
        mn = np.min([np.min(zx),np.min(zy) + offset]);                          # Minimum y-axis value
        mx = np.max([np.max(zx),np.max(zy) + offset]);                          # Maximum y-axis value
        wt = 0.05*(mx - mn)                                                     # Add a little either side to look nice
        self.ax.set_ylim((mn - wt,mx + wt))                                     # Set the ylim
        self.ax.grid()                                                          # Put some grid lines on the plot
        
        self.ax.set_position([0.13, 0.1, 0.83, 0.83])                           # Leave room for axis labels and title
        # self.ax.set_title("Line Profiles")                                      # Slap a title on there
    
    def plotP1P2(self):
        im  = np.copy(self.mainPanel.finalim)
        dxy = np.copy(self.mainPanel.dxy)
        sx = len(im[0]); sy = len(im)                                           # size of the image in pixels
        
        self.segInfo = []
        for cPos in self.cPos:
            p1 = cPos[0]*sx
            p2 = cPos[1]*sy
        
            temp = []
            if(p1[0] > p2[0]): temp = p2; p2 = p1; p1 = temp
            
            m  = (p1[1] - p2[1])/(p1[0] - p2[0])
            theta = math.atan(m)
            dz = abs(1/math.cos(theta))
            c  = p1[1] - m*p1[0]
            Ix = np.linspace(p1[0],p2[0],int(dz*(p2[0] - p1[0])))
            Iy = sy - np.array(m*Ix + c).astype(int)
            Ix = Ix.astype(int)
            zx = im[Iy,Ix]/self.zunit                                               # Raw 1D cut through two points
            
            xx = np.linspace(0,(abs(p2[0] - p1[0]))*dz,len(Ix))*dxy[0]/self.xunit
            
            p1 = p1*dxy
            p2 = p2*dxy 
            length = math.sqrt((p1[1] - p2[1])**2 + (p1[0] - p2[0])**2)/self.xunit
            self.segInfo.append(np.array([length,180*theta/math.pi-self.mainPanel.scanAngle]))
            
            self.ax.plot(xx,zx);
        self.ax.set_xlabel('Position (nm)'); self.ax.set_ylabel('z (pm)')
        self.ax.grid()
    ###########################################################################
    # Placing Cursors
    ###########################################################################
    def cursor(self,c):
        if(self.activeCursor[0] > -1): return                                   # If we're already placing a cursor, don't let the button do anything
        if(c > 0 and self.plotModes[self.plotMode] == "XY"): return             # Can't place cursor 1 if we're in xy mode..
        self.activeCursor[0] = c                                                # If we made it this far, we're placing cursor c
        self.mainPanel.cursorBind()                                             # Bind the mouse on main panel
    
    def setCursor(self,cPos,finalSet = False):
        self.cPos[self.activeCursor[1]][self.activeCursor[0]] = cPos
        if(finalSet): self.activeCursor[0] = -1
    ###########################################################################
    # Multiple Cursors
    ###########################################################################
    def nextCursor(self):
        numCursors = len(self.cPos)
        self.activeCursor[1] += 1
        if(self.activeCursor[1] == numCursors): self.activeCursor[1] = 0
        c = self.mainPanel.mplibColours[self.activeCursor[1]]                   # Get the default matplotlib colour for this line so it matches the colour on profile panel                                    
        self.btn['Next Cursor'].configure(bg=c)
        self.update()
        self.mainPanel.update()
        
    def addCursor(self):
        self.cPos.append(np.array([[0.5,0.5],[0.75,0.75]]))                     # Just add this as initial cursor positions
        self.activeCursor[1] = len(self.cPos) - 1                               # Auto select this new cursor set
        c = self.mainPanel.mplibColours[self.activeCursor[1]]                   # Get the default matplotlib colour for this line so it matches the colour on profile panel                                    
        self.btn['Next Cursor'].configure(bg=c)
        self.update()
        self.mainPanel.update()
    
    def remCursor(self):
        if(len(self.cPos) > 1):                                                 # Keep at least one cursor there always
            del self.cPos[self.activeCursor[1]]                                 # Pop this cursor from the list
            if(self.activeCursor[1] > len(self.cPos) - 1):                      # In case the active cursor was the last one
                self.activeCursor[1] = len(self.cPos) - 1
        c = self.mainPanel.mplibColours[self.activeCursor[1]]                   # Get the default matplotlib colour for this line so it matches the colour on profile panel                                    
        self.btn['Next Cursor'].configure(bg=c)
        self.update()
        self.mainPanel.update()
    ###########################################################################
    # Misc Button Functions
    ###########################################################################
    def toggleMode(self):
        self.plotMode += 1
        if(self.plotMode == len(self.plotModes)): self.plotMode = 0
        self.btn["Mode"].configure(text="Mode: " + self.plotModes[self.plotMode])
        self.update()
        self.mainPanel.update()
    ###########################################################################
    # Save
    ###########################################################################
    def buildSaveString(self):
        saveString = "#LineProfilePanel\n"                                      # Line 1: Header
        
        saveString += str(self.plotMode)  + "\n"                                # Line 2: Plot Mode
        saveString += str(len(self.cPos)) + "\n"                                # Line 3: Number of cursors to follow (might have more than two cursors in future)
        for cursor in self.cPos:
            saveString += ','.join("{:.5f}".format(pos)                         # Line 4...: cursor position
                             for pos in cursor) + "\n"
        
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
                if(line == "#LineProfilePanel"): headerFound = True
            if(not headerFound): print("Missing #LineProfilePanel"); return
            
            self.plotMode = int(f.readline()[:-1])                              # Line 2: Plot Mode
            
            numCursor = int(f.readline()[:-1])                                  # Line 3: Number of cursors to follow
            
            self.cPos = np.empty((numCursor,2))                                 # Line 4...: cursor position
            for c in range(numCursor):
                cPos = np.array(f.readline()[:-1].split(','))
                self.cPos[c] = cPos.astype(float)
        
        self.mainPanel.update()