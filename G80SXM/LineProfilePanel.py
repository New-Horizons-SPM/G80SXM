# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 18:03:13 2022

@author: jced0001
"""

from Panel import Panel
import customtkinter as ctk
import numpy as np
import math
class LineProfilePanel(Panel):
    plotModes = ["XY", "P1P2"]; plotMode = 0                                    # XY plots 1D profiles through a crosshair at cPos[0].
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
        self.cPos = [np.array([[0.5,0.5],[0.75,0.75]])]                         # P1P2 plots a 1D profile through the line segment from cPos[0] to cPos[1]
        self.activeCursor = np.array([-1,0])                                    # [0] -1=not placing any cursors atm. 0=Placing P1 (used in XY and P1P2 modes) 1=Placing P2 (used in only P2 mode when drawing a line b/w P1 and P2)
        self.segInfo  = [[-1,-1,1,0]]                                           # [0]: segment length. [1]: segment angle. [2]: showInfo. [3]: segment line colour.
        self.fitLocations = [[]]                                                # List of locations along the 1D line segment (in P1P2 mode) to fit linear lines between
        self.fitProfileActive = False                                           # Flag for when we're currently fitting a line segment (in P1P2 mode)
    ###########################################################################
    # Panel
    ###########################################################################
    def buttons(self):
        self.btn = {
            "Mode"          : ctk.CTkButton(self.master, text="Mode: xy",       command=self.toggleMode),
            "Add Cursor"    : ctk.CTkButton(self.master, text="Add Line",       command=self.addCursor),
            "Rem Cursor"    : ctk.CTkButton(self.master, text="Remove Line",    command=self.remCursor),
            "Next Cursor"   : ctk.CTkButton(self.master, text="Next Cursor",    command=self.nextCursor),
            "Cursor 1"      : ctk.CTkButton(self.master, text="P1 Location",    command=lambda:self.cursor(0)),
            "Cursor 2"      : ctk.CTkButton(self.master, text="P2 Location",    command=lambda:self.cursor(1)),
            "Info"          : ctk.CTkComboBox(self.master,values=["Line Info"], command=self.labelInfo),
            "LineColour"    : ctk.CTkButton(self.master, text="Line Colour",    command=self.changeLineColour),
            "Fit Steps"     : ctk.CTkButton(self.master, text="Fit Steps",      command=self.fitSteps),
            "Inset"         : ctk.CTkButton(self.master, text="Inset",          command=super().addInset),
            "Imprint"       : ctk.CTkButton(self.master, text="Imprint",        command=super()._imprint),
            "Close"         : ctk.CTkButton(self.master, text="Close",          command=self.destroy)
            }
        
        infoValues = ["Line Info","Show/Hide","Move (Coming)"]
        self.btn['Info'].configure(values=infoValues,fg_color=['#3B8ED0', '#1F6AA5'])
        
    def buttonHelp(self):
        helpStr = "Toggle between XY and Point-toPoint modes.\nXY mode plots horizontal and vertical cuts through the sxm image at a point defined by 'P1'.\nPoint-to-Point mode plots a line segment defined by 'P1' and 'P2'"
        self.btn['Mode'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Add a new line to the sxm plot (Point-to-Point Mode Only)"
        self.btn['Add Cursor'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Remove the currently selected line"
        self.btn['Rem Cursor'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Cycle to the next line to modify its properties (Point-to-Point mode only)"
        self.btn['Next Cursor'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Reposition the lines in XY mode or P1 in Point-to-Point mode by clicking on the main figure"
        self.btn['Cursor 1'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Reposition P2 in Point-to-Point mode by clicking on the main figure"
        self.btn['Cursor 2'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Show/Hide the line segment info or Move it by clicking on the main figure (Point-to-Point mode only)"
        self.btn['Info'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Change the colour of the selected line (Point-toPoint mode only)"
        self.btn['LineColour'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Fit step edges to the 1D line cuts by placing a cursor either side of each step edge.\nLeft click to place a cursor. Right click to end"
        self.btn['Fit Steps'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Add the above plot as an inset on the main figure. Double click a location in the main figure to repoisition the inset and use the scroll wheel to change its size"
        self.btn['Inset'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Imprint the overlay drawn by this panel on the main figure so it persits after closing this panel"
        self.btn['Imprint'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Close this panel"
        self.btn['Close'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
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
        zx = im[int(sy*(1-self.cPos[idx][0][1])),:]/self.zunit                  # Raw 1D horizontal cut
        zy = im[:,int(sx*self.cPos[idx][0][0])]/self.zunit                      # Raw 1D vertical cut
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
        
        for idx,cPos in enumerate(self.cPos):
            p1 = cPos[0]*sx
            p2 = cPos[1]*sy
        
            temp = []
            if(p1[0] > p2[0]): temp = p2; p2 = p1; p1 = temp
            
            m  = (p1[1] - p2[1])/(p1[0] - p2[0])
            theta = math.atan(m)
            dz = abs(1/math.cos(theta))
            c  = p1[1] - m*p1[0]
            Ix = np.linspace(p1[0],p2[0],int(dz*(p2[0] - p1[0])))
            Iy = sy - np.array(m*Ix + c).astype(int) - 1
            Ix = Ix.astype(int)
            zx = im[Iy,Ix]/self.zunit                                               # Raw 1D cut through two points
            
            xx = np.linspace(0,(abs(p2[0] - p1[0]))*dz,len(Ix))*dxy[0]/self.xunit
            
            p1 = p1*dxy
            p2 = p2*dxy 
            length = math.sqrt((p1[1] - p2[1])**2 + (p1[0] - p2[0])**2)/self.xunit
            
            self.segInfo[idx][0] = length
            self.segInfo[idx][1] = 180*theta/math.pi-self.mainPanel.scanAngle
            col = self.mainPanel.mplibColours[self.segInfo[idx][3]]
            
            # If we're placing the fit lines, leave 1D data (zx) as is and plot fit lines on top
            if(idx == self.activeCursor[1] and self.fitProfileActive):
                self.ax.axvline(x=self.motionFitX,c=col)
                
                for X in self.fitLocations[idx]:
                    self.ax.axvline(x=X,c=col)
            
            # If we have already placed fitLocations, fit the data (zx) with linear line segments
            hasFit = len(self.fitLocations[idx])
            if(hasFit and not self.fitProfileActive):
                zx = self.fitzx(idx,xx,zx)
            
            # If we're currently placing fitLocations for another segment, don't show this one
            if(self.fitProfileActive and idx != self.activeCursor[1]):
                continue
            
            if(col == 'white'): continue                                        # Don't plot if white. Looks weird on the grid
        
            self.ax.plot(xx,zx,c=col)
            
        self.ax.set_xlabel('Position (nm)'); self.ax.set_ylabel('z (pm)')
        self.ax.grid()
    
    def fitzx(self,idx,xx,zx):
        fs = 10                                                                 # Annotation font size
        xSteps = []; xPos = []                                                  # xStep edge heights and xpos for annotations
        xFit = self.fitLocations[idx]
        if(xFit):                                                               # Take care of the first marker
            idx = xx < xFit[0]                                                  # Index of all the points to the left of marker 0
            m,c = np.polyfit(xx[idx],zx[idx],1)                                 # Fit all these points with a straight line
            zx[idx] = m*xx[idx] + c                                             # Equation of the straight line
            
            lastOne = len(idx) - list(idx[::-1]).index(1) - 1
            xSteps.append(np.average(zx[lastOne]));                             # Base of the step to the right
            xPos.append(xx[lastOne])                                            # Position of the step edge (for annotating)
            
            # self.ax[1].annotate("{:.0f}".format(zx[idx.tolist().index(1)]), xy=(xx[idx.tolist().index(1)], zx[idx.tolist().index(1)] + 10),fontsize=fs,color='green')
            
        for x in range(2,len(xFit),2):                                          # Take care of all steps in the middle
            idx = (xx > xFit[x-1]) & (xx < xFit[x])                             # Step is between this fit point and the previous one
            m,c = np.polyfit(xx[idx],zx[idx],1)                                 # Linear fit
            zx[idx] = m*xx[idx] + c                                             # Equation of the line
            if(sum(idx)):                                                       # If there are any points between them figure out the step height
                firstOne = list(idx).index(1)
                xSteps.append(np.average(zx[firstOne]));                        # Height of the step to the left
                lastOne = len(idx) - list(idx[::-1]).index(1) - 1
                xSteps.append(np.average(zx[lastOne]));                         # Base of the step to the Right
                xPos.append(xx[lastOne])                                        # The location of the step along x axis... used for annotating
                
                # self.ax[1].annotate("{:.0f}".format(zx[idx.tolist().index(1)]), xy=(xx[idx.tolist().index(1)], zx[idx.tolist().index(1)] + 10),fontsize=fs,color='green')
        
        if(xFit):                                                               # Take care of the last step
            idx = xx > xFit[-1]                                                 # Which is to the right of the last marker/fit point
            m,c = np.polyfit(xx[idx],zx[idx],1)
            zx[idx] = m*xx[idx] + c
            
            firstOne = list(idx).index(1)
            xSteps.append(np.average(zx[firstOne]));                            # Height of the step to the left
            
            # self.ax[1].annotate("{:.0f}".format(zx[idx.tolist().index(1)]), xy=(xx[idx.tolist().index(1)], zx[idx.tolist().index(1)] + 10),fontsize=fs,color='green')
        
        for i in range(len(xPos)):                                              # For each step edge...
            idx=i*2
            dz = xSteps[idx+1]-xSteps[idx]                                      # Calculate step height
            yPos = xSteps[idx] + dz/2                                           # Go to 
            
            self.ax.annotate("{:.0f} pm".format(abs(dz)), xy=(xPos[i], yPos),fontsize=fs,color='black')
        
        return zx
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
        if(not self.plotModes[self.plotMode] == "P1P2"): return
        numCursors = len(self.cPos)
        self.activeCursor[1] += 1
        if(self.activeCursor[1] == numCursors): self.activeCursor[1] = 0
        cidx = self.segInfo[self.activeCursor[1]][3]
        c = self.mainPanel.mplibColours[cidx]                                   # Colour of the next cursor
        self.btn['Next Cursor'].configure(fg_color=c)
        self.mainPanel.update(upd=[0,1])
        
    def addCursor(self):
        if(not self.plotModes[self.plotMode] == "P1P2"): return
        self.cPos.append(np.array([[0.5,0.5],[0.75,0.75]]))                     # Just add this as initial cursor positions
        self.activeCursor[1] = len(self.cPos) - 1                               # Auto select this new cursor set
        c = self.mainPanel.mplibColours[self.activeCursor[1]]                   # Next colour in the list
        self.btn['Next Cursor'].configure(fg_color=c)
        self.segInfo.append([-1,-1,1,self.activeCursor[1]])                     # Default to c colour and show segment info = 1
        self.fitLocations.append([])
        self.mainPanel.update(upd=[0,1])
    
    def remCursor(self):
        if(not self.plotModes[self.plotMode] == "P1P2"): return
        if(len(self.cPos) > 1):                                                 # Keep at least one cursor there always
            del self.cPos[self.activeCursor[1]]                                 # Remove this cursor from the list
            del self.segInfo[self.activeCursor[1]]
            del self.fitLocations[self.activeCursor[1]]
            if(self.activeCursor[1] > len(self.cPos) - 1):                      # In case the active cursor was the last one
                self.activeCursor[1] = len(self.cPos) - 1
                
        cidx = self.segInfo[self.activeCursor[1]][3]
        c = self.mainPanel.mplibColours[cidx]                                   # Colour of the next available cursor
        self.btn['Next Cursor'].configure(bg=c)
        self.mainPanel.update(upd=[0,1])
    
    def changeLineColour(self):
        if(not self.plotModes[self.plotMode] == "P1P2"): return
        numColours = len(self.mainPanel.mplibColours)
        self.segInfo[self.activeCursor[1]][3] += 1                              # Next colour in the list
        if(self.segInfo[self.activeCursor[1]][3] == numColours):
            self.segInfo[self.activeCursor[1]][3] = 0                           # Cycle back to the first colour
        
        cidx = self.segInfo[self.activeCursor[1]][3]
        c = self.mainPanel.mplibColours[cidx]
        self.btn['Next Cursor'].configure(fg_color=c)
        self.mainPanel.update(upd=[0,1])
        
    ###########################################################################
    # Move Label Info
    ###########################################################################
    def toggleShowInfo(self):
        self.segInfo[self.activeCursor[1]][2] = not self.segInfo[self.activeCursor[1]][2]
        self.mainPanel.update()
    
    def moveLabel(self):
        pass
    ###########################################################################
    # Fitting Step Edges
    ###########################################################################
    def fitSteps(self):
        if(self.plotModes[self.plotMode] == "XY"): return
        self.fitLocations[self.activeCursor[1]] = []
        self.fitBind()
        
    def fitBind(self):
        self.placeFitBind  = self.canvas.get_tk_widget().bind('<Button-1>', self.placeFit)
        self.setFitBind    = self.canvas.get_tk_widget().bind('<Button-3>', self.setFit)
        self.motionFitBind = self.canvas.get_tk_widget().bind('<Motion>', self.motionFit)
        self.fitProfileActive  = True;
        
    def fitUnbind(self):
        self.canvas.get_tk_widget().unbind('<Button-1>', self.placeFitBind)
        self.canvas.get_tk_widget().unbind('<Button-3>', self.setFitBind)
        self.canvas.get_tk_widget().unbind('<Motion>', self.motionFitBind)
        self.fitProfileActive  = False;
        
    def motionFit(self,event):
        x = event.x
        X = super()._getX(x)
        self.motionFitX = X
        self.update()
    
    def placeFit(self,event):
        x = event.x
        X = super()._getX(x)
        
        if(not self.fitLocations[self.activeCursor[1]]):                        # This is hacky but seems to be a bug e.g...
            self.fitLocations[self.activeCursor[1]] = [X]                       # a = [[]]*3; a[1].append(1); print(a)
        else:                                                                   # yields [[1], [1], [1]],
            self.fitLocations[self.activeCursor[1]].append(X)                   # not [[], [1], []] as expected
        self.update()
        
    def setFit(self,event):
        self.fitUnbind()
        self.update()
    ###########################################################################
    # Misc Button Functions
    ###########################################################################
    def toggleMode(self):
        self.plotMode += 1
        if(self.plotMode == len(self.plotModes)): self.plotMode = 0
        self.btn["Mode"].configure(text="Mode: " + self.plotModes[self.plotMode])
        self.mainPanel.update(upd=[0,1])
        
    def labelInfo(self,option):
        if(not self.plotModes[self.plotMode] == "P1P2"):
            self.btn['Info'].set("Line Info")
            return
        if(option == "Show/Hide"): self.toggleShowInfo()
        if(option == "Move"): self.moveLabel()
        self.btn['Info'].set("Line Info")
        
    ###########################################################################
    # Save
    ###########################################################################
    def buildSaveDict(self):
        saveDict = {}
        
        saveDict['plotMode'] = self.plotMode
        saveDict['cPos']     = self.cPos
        saveDict['segInfo']  = self.segInfo
        saveDict['fitLocations']  = self.fitLocations
        saveDict['imprint']  = self.imprint
        
        return saveDict
    
    ###########################################################################
    # Load
    ###########################################################################
    def loadFromDict(self,loadDict):
        for key,value in loadDict.items():
            setattr(self,key,value)