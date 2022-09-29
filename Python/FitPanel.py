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
        self.plotSTS()                                                          # Plots the sts curve selected from the sts panel
        self.plotFit()                                                          # 
        self.ax.set_position([0.13, 0.1, 0.83, 0.83])                           # Leave room for axis labels and title
        
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()                                                      # Redraw the canvas with the updated figure
        print(self.fitDict)
    
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
    
    def plotFit(self):
        result = self.fit()
        if(not result): return
        print(fit_report(result))
        c = self.mainPanel.mplibColours[self.curveIdx + 1]
        self.ax.plot(self.curve[0], result.init_fit, '--', label='init fit', c=c)
        self.ax.plot(self.curve[0], result.best_fit, '-', label='best fit', c=c)
        
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
    
    def showFermiDiracForm(self):
        for l in self.fermiLabels:
            l[0].grid(row=l[1],column=l[2])
            
        for idx,e in enumerate(self.fermiEntries):
            e[0].grid(row=e[1],column=e[2])
            if(self.componentIdx < 0): continue
            print(">>",self.componentIdx)
            e[0].delete(0,tk.END)
            e[0].insert(0,self.fitDict["FermiDirac"][self.componentIdx][idx])
        
        for b in self.fermiBtn:
            if(b[0]['text'] == "remove" and self.componentIdx == -1): continue
            b[0].grid(row=b[1],column=b[2])
        
    def hideFermiDiracForm(self):
        for l in self.fermiLabels:
            l[0].grid_forget()
            
        for e in self.fermiEntries:
            e[0].grid_forget()
        
        for b in self.fermiBtn:
            b[0].grid_forget()
            
    def submitFermiDiracForm(self):
        params = []
        for e in self.fermiEntries:
            try:
                params.append(np.float32(e[0].get()))
            except:
                self.updateHelpLabel("Error in form: All values must be numeric.")
                return
        
        if(self.componentIdx == -1):
            self.fitDict["FermiDirac"].append(params)
        else:
            self.fitDict["FermiDirac"][self.componentIdx] = params
        
        self.hideFermiDiracForm()
        self.update()
        self.updateEditButton()
        self.componentIdx = -1
    
    def updateEditButton(self):
        menu = tk.Menu(self.btn['Edit'], tearoff=0)
        if("FermiDirac" in self.fitDict):
            for idx,f in enumerate(self.fitDict["FermiDirac"]):
                menu.add_command(label="Fermi-Dirac " + str(idx),command=lambda x=idx: self.editFermiDirac(x)) # Can't bind loop variable to lambda in Python. Need to do it this way (indexes[idx])
        
        self.btn['Edit']['menu'] = menu
        
    def editFermiDirac(self,index):
        self.componentIdx = index
        self.showFermiDiracForm()
        
    def cancelFermiDiracForm(self):
        self.hideFermiDiracForm()
        self.componentIdx = -1
    
    def removeFermiDiracForm(self):
        del self.fitDict["FermiDirac"][self.componentIdx]
        self.componentIdx = -1
        self.updateEditButton()
        self.hideFermiDiracForm()
        self.update()
        
    def reset(self):
        self.fitDict = {}
        self.update()
    
    def fit(self):
        model = 0
        pars = Parameters()
        xx = self.curve[0]
        yy = self.curve[1]
        if("FermiDirac" in self.fitDict):
            fermiEdge = []
            fermiDirac = self.fitDict["FermiDirac"]
            for idx,fd in enumerate(fermiDirac):
                A  = fd[0]; Amin  = fd[1]; Amax  = fd[2]
                x0 = fd[3]; x0min = fd[4]; x0max = fd[5]
                kT = 8.617e-5*fd[6]
                
                print(x0,x0min,x0max)
                pars.add('FD' + str(idx) + '_A',  value=A,  min=Amin,  max=Amax)
                pars.add('FD' + str(idx) + '_x0', value=x0, min=x0min, max=x0max)
                pars.add('FD' + str(idx) + '_T',  value=-kT, vary=False)
                pars.add('FD' + str(idx) + '_c',  value=A,  min=0,     max=2*A)
                
                fermiEdge.append(Model(self.fermiDiracFunc,prefix='FD' + str(idx) + '_'))
                if(not model): model = fermiEdge[-1]
                else:          model += fermiEdge[-1]
        
        if("Gaussian" in self.fitDict):
            pass
        
        if("Reference" in self.fitDict):
            pass
        
        if(model == 0): return 0
        
        model.eval(pars,x=xx)
        return model.fit(yy,pars,x=xx)
        
    ###########################################################################
    # Custom Fitting Curves (not in lmfit)
    ###########################################################################
    def fermiDiracFunc(self,x, A, x0, T, c):
        """
        Fermi-Dirac function

        Parameters
        ----------
        x  : x-axis
        A  : Amplitude
        x0 : Onset
        T  : Temperature (K)

        Returns
        -------
        Fermi-Dirac Function

        """
        # step = x*0
        # step[x < x0] = 0
        # step[x >=x0] = A
        # step += c
        # return step
        # return A*np.heaviside(x-x0,0.5) + c
        # return -A/(np.exp((x-x0)/(-8.617e-5*T)) + 1) + c
        return A/(1+np.exp((x0-x)/(-8.617e-5*T))) + c
    
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