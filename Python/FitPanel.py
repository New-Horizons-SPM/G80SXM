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
    forms = {}
    formActive = False
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width, height, dpi, mainPanel):
        super().__init__(master, width, height, dpi, mainPanel=mainPanel)
        self.buttons()
        self.reset()
        
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
        menu.add_command(label="Reference",   command=lambda: self.addFitCurve("Reference"))
        menu.add_command(label="Gaussian",    command=lambda: self.addFitCurve("Gaussian"))
        menu.add_command(label="Fermi-Dirac", command=lambda: self.addFitCurve("FermiDirac"))
        self.btn['Add']['menu'] = menu
        
    def special(self):
        # Reference Form
        params = []
        params.append(['A','Amin','Amax'])
        params.append(['c','cmin','cmax'])
        self.buildForm(name="Reference", params=params)
        
        # Fermi-Dirac Form
        params = []
        params.append(['A','Amin','Amax'])
        params.append(['x0','x0min','x0max'])
        params.append(['T'])
        self.buildForm(name="FermiDirac", params=params)
        
        # Gaussian Form
        params = []
        params.append(['A','Amin','Amax'])
        params.append(['x0','x0min','x0max'])
        params.append(['sigma','sigmaMin','sigmaMax'])
        params.append(['c','cmin','cmax'])
        self.buildForm(name="Gaussian", params=params)
        
    def removeSpecial(self):
        self.hideFermiDiracForm()
    
    def buildForm(self,name,params):
        self.forms[name] = {"labels"  : [],
                            "entries" : [],
                            "buttons" : []}
        idr = 0; row = 7
        for idr,p in enumerate(params):
            for idp,param in enumerate(p):
                self.forms[name]['labels'].append([tk.Label(self.master, text=param),row+idr,self.pos + 2*idp])
                self.forms[name]['entries'].append([tk.Entry(self.master),row+idr,self.pos+2*idp+1])
        
        idr += 1
        self.forms[name]['buttons'] = []
        self.forms[name]['buttons'].append([tk.Button(self.master, text="submit", command=lambda n=name: self.submitForm(n)),row+idr,self.pos + 1])
        self.forms[name]['buttons'].append([tk.Button(self.master, text="cancel", command=lambda n=name: self.cancelForm(n)),row+idr,self.pos + 2])
        self.forms[name]['buttons'].append([tk.Button(self.master, text="remove", command=lambda n=name: self.removeForm(n)),row+idr,self.pos + 3])
        
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
        self.ax.plot(V,didv,linewidth=1.5,c=c)
       
        self.ax.set_xlabel("Bias (V)")
        self.ax.set_ylabel("dI/dV (arb)")
        self.ax.set_title("Point Spectroscopy")
    
    def plotFit(self):
        result = self.fit()
        if(not result): return
        print(fit_report(result))
        c = self.mainPanel.mplibColours[self.curveIdx + 1]
        self.ax.plot(self.curve[0], result.init_fit, '--', label='init fit', c=c, linewidth=0.75, alpha=0.25)
        self.ax.plot(self.curve[0], result.best_fit, '--', label='best fit', c=c, linewidth=1,    alpha=0.8)
        
    ###########################################################################
    # Curve Fitting...
    ###########################################################################
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
                
                pars.add('FD' + str(idx) + '_A',  value=A,  min=Amin,  max=Amax)
                pars.add('FD' + str(idx) + '_x0', value=x0, min=x0min, max=x0max)
                pars.add('FD' + str(idx) + '_T',  value=-kT, vary=False)
                pars.add('FD' + str(idx) + '_c',  value=A,  min=0,     max=2*A)
                
                fermiEdge.append(Model(self.fermiDiracFunc,prefix='FD' + str(idx) + '_'))
                if(not model): model = fermiEdge[-1]
                else:          model += fermiEdge[-1]
        
        if("Gaussian" in self.fitDict):
            gaussian = self.fitDict["Gaussian"]
            for idx,g in enumerate(gaussian):
                A  = g[0]; Amin  = g[1]; Amax  = g[2]
                x0 = g[3]; x0min = g[4]; x0max = g[5]
                sigma = g[6]; sigmaMin = g[7]; sigmaMax = g[8]
                c  = g[9]; cmin  = g[10]; cmax = g[11]
                
                pars.add('GAUSS' + str(idx) + '_center',    value=x0,    min=x0min,    max=x0max)
                pars.add('GAUSS' + str(idx) + '_amplitude', value=A,     min=Amin,     max=Amax)
                pars.add('GAUSS' + str(idx) + '_sigma',     value=sigma, min=sigmaMin, max=sigmaMax)
                pars.add('GAUSS' + str(idx) + '_c',         value=c,     min=cmin,     max=cmax)
                
                if(not model): model  = GaussianModel(prefix='GAUSS' + str(idx) + '_')
                else:          model += GaussianModel(prefix='GAUSS' + str(idx) + '_')
                
                model += ConstantModel(prefix='GAUSS' + str(idx) + '_')
        
        if("Reference" in self.fitDict):
            reference = self.fitDict["Reference"]
            for idx,ref in enumerate(reference):
                A = ref[0]; Amin = ref[1]; Amax = ref[2]
                c = ref[3]; cmin = ref[4]; cmax = ref[5]
                
                pars.add('REF' + str(idx) + '_A', value=A, min=Amin, max=Amax)
                pars.add('REF' + str(idx) + '_c', value=c, min=cmin, max=cmax)
                
                if(not model): model  = Model(self.referenceFunc,prefix='REF' + str(idx) + '_')
                else:          model += Model(self.referenceFunc,prefix='REF' + str(idx) + '_')
            
        if(model == 0): return 0
        
        model.eval(pars,x=xx)
        return model.fit(yy,pars,x=xx)
    
    ###########################################################################
    # Form Actions (Show, Submit, Cancel, Remove, etc.)
    ###########################################################################
    def showForm(self,name):
        if(self.formActive): return
        self.formActive = True
        
        for l in self.forms[name]['labels']:
            l[0].grid(row=l[1],column=l[2])
            
        for idx,e in enumerate(self.forms[name]['entries']):
            e[0].grid(row=e[1],column=e[2])
            if(self.componentIdx < 0): continue
            e[0].delete(0,tk.END)
            e[0].insert(0,self.fitDict[name][self.componentIdx][idx])
        
        for b in self.forms[name]['buttons']:
            if(b[0]['text'] == "remove" and self.componentIdx == -1): continue  # Only show the 'remove' button if we're editing a selection
            b[0].grid(row=b[1],column=b[2])
            
    def hideForm(self,name):
        self.formActive = False
        
        for l in self.forms[name]['labels']:
            l[0].grid_forget()
            
        for e in self.forms[name]['entries']:
            e[0].grid_forget()
        
        for b in self.forms[name]['buttons']:
            b[0].grid_forget()
            
    def submitForm(self,name):
        params = []
        for e in self.forms[name]['entries']:
            try:
                params.append(np.float32(e[0].get()))
            except:
                self.updateHelpLabel("Error in form: All values must be numeric.")
                return
        
        if(self.componentIdx == -1):
            self.fitDict[name].append(params)
        else:
            self.fitDict[name][self.componentIdx] = params
        
        self.hideForm(name)
        self.update()
        self.updateEditButton()
        self.componentIdx = -1
        
    def cancelForm(self,name):
        self.hideForm(name=name)
        self.componentIdx = -1
        
    def removeForm(self,name):
        del self.fitDict[name][self.componentIdx]
        self.componentIdx = -1
        self.updateEditButton()
        self.hideForm(name=name)
        self.update()
        
    def editForm(self,name,index):
        self.componentIdx = index
        self.showForm(name=name)
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
    
    def referenceFunc(self,x,A,c):
        reference = self.mainPanel.stsPanel.getReferenceForCurve(x=x)
        return A*reference + c
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
        
    def addFitCurve(self,name):
        if(not name in self.fitDict):
            self.fitDict[name] = []
        self.showForm(name=name)
        
    def updateEditButton(self):
        menu = tk.Menu(self.btn['Edit'], tearoff=0)
        for key in self.fitDict.keys():
            for idx,p in enumerate(self.fitDict[key]):
                menu.add_command(label=key + " " + str(idx),command=lambda k=key,x=idx: self.editForm(k,x)) # Can't bind loop variable to lambda in Python. Need to do it this way (indexes[idx])
                
        # if("Reference" in self.fitDict):
        #     for idx,f in enumerate(self.fitDict["Reference"]):
        #         menu.add_command(label="Reference " + str(idx),command=lambda x=idx: self.editForm("Reference",x)) # Can't bind loop variable to lambda in Python. Need to do it this way (indexes[idx])
                
        # if("FermiDirac" in self.fitDict):
        #     for idx,f in enumerate(self.fitDict["FermiDirac"]):
        #         menu.add_command(label="Fermi-Dirac " + str(idx),command=lambda x=idx: self.editForm("FermiDirac",x)) # Can't bind loop variable to lambda in Python. Need to do it this way (indexes[idx])
        
        # if("Gaussian" in self.fitDict):
        #     for idx,f in enumerate(self.fitDict["Gaussian"]):
        #         menu.add_command(label="Gaussian " + str(idx),command=lambda x=idx: self.editForm("Gaussian",x)) # Can't bind loop variable to lambda in Python. Need to do it this way (indexes[idx])
                
        self.btn['Edit']['menu'] = menu
        
    def reset(self):
        self.fitDict = {}
        self.update()