# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 11:23:21 2022

@author: jced0001
"""

from Panel import Panel
import tkinter as tk
import numpy as np

from lmfit import Model, Parameters, fit_report
from lmfit.models import GaussianModel, ConstantModel

class FitPanel(Panel):
    curveIdx = -1
    componentIdx = -1
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width, height, dpi, mainPanel):
        super().__init__(master, width, height, dpi, mainPanel=mainPanel)
        self.buttons()
        self.reset()
        # self.fermiDiracForm()
        
    ###########################################################################
    # Panel
    ###########################################################################
    def buttons(self):
        self.btn = {
            "Next":     tk.Button(self.master, text="Next Curve", command=self.nextCurve),
            "Add":      tk.Menubutton(self.master, text="Add",  relief=tk.RAISED),
            "Edit":     tk.Menubutton(self.master, text="Edit", relief=tk.RAISED),
            # "Scale":    tk.Button(self.master, text="Linear",     command=self._scale),
            # "Undo":     tk.Button(self.master, text="Undo Last",  command=self._undo),
            # "Reset":    tk.Button(self.master, text="Reset",      command=self._reset),
            # "Imprint":  tk.Button(self.master, text="Imprint",    command=super()._imprint),
            "Close":    tk.Button(self.master, text="Close",      command=self.destroy)
            }
        
        menu = tk.Menu(self.btn['Add'],tearoff=0)
        menu.add_command(label="Reference",   command=self.addReference)
        menu.add_command(label="Gaussian",    command=self.addGaussian)
        menu.add_command(label="Fermi-Dirac", command=self.addFermiDirac)
        self.btn['Add']['menu'] = menu
        
    def special(self):
        params = []
        params.append(['A','Amin','Amax'])
        params.append(['x0','x0min','x0max'])
        params.append(['T'])
        
        idr = 0; row = 7
        self.fermiEntries = []
        self.fermiLabels  = []
        for idr,p in enumerate(params):
            for idp,param in enumerate(p):
                self.fermiLabels.append([tk.Label(self.master, text=param),row+idr,self.pos + 2*idp])
                self.fermiEntries.append([tk.Entry(self.master),row+idr,self.pos+2*idp+1])
                
        idr += 1
        self.fermiBtn = []
        self.fermiBtn.append([tk.Button(self.master, text="submit", command=self.submitFermiDiracForm),row+idr,self.pos + 1])
        self.fermiBtn.append([tk.Button(self.master, text="cancel", command=self.cancelFermiDiracForm),row+idr,self.pos + 2])
        self.fermiBtn.append([tk.Button(self.master, text="remove", command=self.removeFermiDiracForm),row+idr,self.pos + 3])
        
    def removeSpecial(self):
        self.hideFermiDiracForm()
        
    ###########################################################################
    # Update and Plotting
    ###########################################################################
    def update(self):
        if(not self.mainPanel.init): return
        if(not self.active): return
        
        self.updateHelpLabel("Use the 'Add' button to add fitting components\n"
                        +    "Then use the 'Fit Params' button to adjust the parameters for each component")
        
        self.ax.cla()                                                           # Clear the axis
        self.plotSTS()                                                          # Loops through .dat files, takes the derivative of IV curves and plot dI/dV
        self.ax.set_position([0.13, 0.1, 0.83, 0.83])                           # Leave room for axis labels and title
        
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()                                                      # Redraw the canvas with the updated figure
    
    def plotSTS(self):
        self.validateCurveIdx()
        if(self.curveIdx < 0):
            self.updateHelpLabel("Add at least one spectrum to the STSPanel to start fitting")
            return
        
        datFiles = self.mainPanel.stsPanel.datFile.copy()
        df = datFiles[self.curveIdx]
        
        self.curve = self.mainPanel.stsPanel.getDIDV(df)
        V, didv = self.curve
        
        c = self.mainPanel.mplibColours[self.curveIdx + 1]
        self.ax.plot(V,didv,linewidth=1.3,c=c)
       
            
        self.ax.set_xlabel("Bias (V)")
        self.ax.set_ylabel("dI/dV (arb)")
        self.ax.set_title("Point Spectroscopy")
    
    ###########################################################################
    # Curve Fitting
    ###########################################################################
    def addReference(self):
        print("Add Reference")
        
    def addGaussian(self):
        print("Add Gaussian")
        
    def addFermiDirac(self):
        if(not "FermiDirac" in self.fitDict):
            self.fitDict["FermiDirac"] = []
        self.showFermiDiracForm()
    
    def showFermiDiracForm(self,index=-1):
        for l in self.fermiLabels:
            l[0].grid(row=l[1],column=l[2])
            
        for e in self.fermiEntries:
            e[0].grid(row=e[1],column=e[2])
        
        for b in self.fermiBtn:
            b[0].grid(row=b[1],column=b[2])
        
    def hideFermiDiracForm(self):
        for l in self.fermiLabels:
            l[0].grid_forget()
            
        for e in self.fermiEntries:
            e[0].grid_forget()
        
        for b in self.fermiBtn:
            b[0].grid_forget()
            
    def submitFermiDiracForm(self,index=-1):
        params = []
        for e in self.fermiEntries:
            try:
                params.append(float(e[0].get()))
            except:
                self.updateHelpLabel("Error in form: All values must be numeric.")
                return
        
        if(index == -1):
            self.fitDict["FermiDirac"].append(params)
        else:
            self.fitDict["FermiDirac"][index] = params
        
        self.hideFermiDiracForm()
        self.update()
        self.updateEditButton()
    
    def updateEditButton(self):
        menu = tk.Menu(self.btn['Edit'], tearoff=0)
        if("FermiDirac" in self.fitDict):
            for idx,f in enumerate(self.fitDict["FermiDirac"]):
                menu.add_command(label="Fermi-Dirac " + str(idx),command=lambda: self.editFermiDirac(idx))
        
        self.btn['Edit']['menu'] = menu
        
    def editFermiDirac(self,index):
        self.showFermiDiracForm(index)
        
    def cancelFermiDiracForm(self):
        self.hideFermiDiracForm()
        
    def reset(self):
        self.pars = Parameters()
        self.fitDict = {}
    
    def fit(self):
        fermiDirac = self.fitDict["FermiDirac"]
        for idx,fd in enumerate(fermiDirac):
            xx = self.curve[0]
            yy = self.curve[1]
            tol = 0.05*abs(np.max(xx) - np.min(xx))
            A  = fd[0]; Amin  = A  - tol; Amax  = A  + tol
            x0 = fd[1]; x0min = x0 - tol; x0max = x0 + tol
            kT = 8.617e-5*fd[2]
            self.pars.add('FD' + idx + '_A',  value=A,  min=Amin,  max=Amax)
            self.pars.add('FD' + idx + '_x0', value=x0, min=x0min, max=x0max)
            self.pars.add('FD' + idx + '_T',  value=-kT, vary=False)
    ###########################################################################
    # Custom Fitting Curves (not in lmfit)
    ###########################################################################
    def fermiDirac(x, A, x0, T):
        """
        Fermi-Dirac function

        Parameters
        ----------
        x  : x-axis
        A  : Amplitude
        x0 : Onset
        T : Temperature (K)

        Returns
        -------
        Fermi-Dirac Function

        """
        return A/(np.exp((x-x0)/(-8.617e-5*T)) + 1)
    
    ###########################################################################
    # Misc Button Functions
    ###########################################################################
    def nextCurve(self):
        self.curveIdx += 1
        numCurves = len(self.mainPanel.stsPanel.datFile)
        if(self.curveIdx >= numCurves): self.curveIdx = 0
        if(not numCurves): self.curveIdx = -1
        self.update()
    
    def validateCurveIdx(self):
        numCurves = len(self.mainPanel.stsPanel.datFile)
        if(not numCurves): self.curveIdx = -1
        if(self.curveIdx >= numCurves): self.curveIdx = numCurves - 1
        if(self.curveIdx < 0 and numCurves): self.curveIdx = 0