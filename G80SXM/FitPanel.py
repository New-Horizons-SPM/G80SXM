# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 11:23:21 2022

@author: jced0001
"""

from .Panel import Panel
import tkinter as tk
import customtkinter as ctk
import numpy as np

from lmfit import Model, Parameters, fit_report
from lmfit.models import GaussianModel, ConstantModel

class FitPanel(Panel):
    curveIdx = -1
    componentIdx = -1
    forms = {}
    formActive = False
    curve = []
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width, height, dpi, mainPanel):
        super().__init__(master, width, height, dpi, mainPanel=mainPanel,length=8,btnSize=2)
        self.buttons()
        self.reset()
        
    ###########################################################################
    # Panel
    ###########################################################################
    def buttons(self):
        self.btn = {
            "Next":     ctk.CTkButton(self.master, text="Next Curve", command=self.nextCurve),
            "Add":      ctk.CTkComboBox(self.master,values=["Add Fit"], command=self.addFitCurve),
            "Edit":     ctk.CTkComboBox(self.master,values=["Edit Fit"],command=self.editForm),
            "Inset":    ctk.CTkButton(self.master, text="Inset",      command=super().addInset),
            "Reset":    ctk.CTkButton(self.master, text="Reset",      command=self.reset),
            "Close":    ctk.CTkButton(self.master, text="Close",      command=self.destroy)
            }
        
        addValues=["Add Fit","Reference","Exponential","Gaussian","Fermi-Dirac","Point-Spectrum"]
        self.btn['Add'].configure(values=addValues,variable="Add Fit")
    
    def buttonHelp(self):
        helpStr = "Cycle through curves on the STS Panel"
        self.btn['Next'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Add a new component to the fit"
        self.btn['Add'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Edit one of the existing fit components"
        self.btn['Edit'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Add the above plot as an inset on the main figure. Double click a location in the main figure to repoisition the inset and use the scroll wheel to change its size"
        self.btn['Inset'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Start from scratch"
        self.btn['Reset'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Close this panel"
        self.btn['Close'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
    def special(self):
        # Reference Form
        params = []
        params.append(['A','Amin','Amax'])
        params.append(['c','cmin','cmax'])
        self.buildForm(name="Reference", params=params)
        
        # Exponential Form
        params = []
        params.append(['A','Amin','Amax'])
        params.append(['c','cmin','cmax'])
        params.append(['P','Pmin','Pmax'])
        self.buildForm(name="Exponential", params=params)
        
        # Fermi-Dirac Form
        params = []
        params.append(['A','Amin','Amax'])
        params.append(['x0','x0min','x0max'])
        params.append(['T'])
        params.append(['Tb','Tbmin','Tbmax'])
        self.buildForm(name="Fermi-Dirac", params=params)
        
        # Gaussian Form
        params = []
        params.append(['A','Amin','Amax'])
        params.append(['x0','x0min','x0max'])
        params.append(['sigma','sigmaMin','sigmaMax'])
        params.append(['c','cmin','cmax'])
        self.buildForm(name="Gaussian", params=params)
        
        # Point Spectrum Form
        params = []
        params.append(['A','Amin','Amax'])
        params.append(['c','cmin','cmax'])
        self.buildForm(name="Point-Spectrum", params=params)
        
    def removeSpecial(self):
        self.hideForm()
    
    def buildForm(self,name,params):
        self.forms[name] = {"labels"  : [],
                            "entries" : [],
                            "buttons" : []}
        idr = 0; row = 7
        for idr,p in enumerate(params):
            for idp,param in enumerate(p):
                self.forms[name]['labels'].append([ctk.CTkLabel(self.master, text=param),row+idr,self.pos + 2*idp])
                self.forms[name]['entries'].append([ctk.CTkEntry(self.master),row+idr,self.pos+2*idp+1])
        
        idr += 1
        self.forms[name]['buttons'] = []
        self.forms[name]['buttons'].append([ctk.CTkButton(self.master, text="submit", command=lambda n=name: self.submitForm(n)),row+idr,self.pos])
        self.forms[name]['buttons'].append([ctk.CTkButton(self.master, text="cancel", command=lambda n=name: self.cancelForm(n)),row+idr,self.pos])
        self.forms[name]['buttons'].append([ctk.CTkButton(self.master, text="remove", command=lambda n=name: self.removeForm(n)),row+idr,self.pos])
        
    ###########################################################################
    # Update and Plotting
    ###########################################################################
    def update(self):
        if(not self.mainPanel.init): return
        if(not self.active): return
        
        self.ax.cla()                                                           # Clear the axis
        self.plotSTS()                                                          # Plots the sts curve selected from the sts panel
        self.plotFit()                                                          # 
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
        self.ax.plot(V,didv,linewidth=1.5,c=c)
       
        self.ax.set_xlabel("Bias (V)")
        self.ax.set_ylabel("dI/dV (arb)")
        self.ax.set_title("Curve Fitting")
    
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
        if(not self.curve): return
        xx = self.curve[0]
        yy = self.curve[1]
        if("Fermi-Dirac" in self.fitDict):
            fermiEdge = []
            fermiDirac = self.fitDict["Fermi-Dirac"]
            for idx,fd in enumerate(fermiDirac):
                A  = fd[0]; Amin  = fd[1]; Amax  = fd[2]
                x0 = fd[3]; x0min = fd[4]; x0max = fd[5]
                kT = 8.617e-5*fd[6]
                Tb = fd[7]; Tbmin = fd[8]; Tbmax = fd[9]
                print(fd)
                pars.add('FD' + str(idx) + '_A',  value=A,  min=Amin,  max=Amax)
                pars.add('FD' + str(idx) + '_x0', value=x0, min=x0min, max=x0max)
                pars.add('FD' + str(idx) + '_T',  value=-kT, vary=False)
                pars.add('FD' + str(idx) + '_Tb', value=Tb, min=Tbmin, max=Tbmax)
                pars.add('FD' + str(idx) + '_c',  value=0,  min=-1,     max=1)
                
                fermiEdge.append(Model(self.fermiDiracFunc,prefix='FD' + str(idx) + '_'))
                if(not model): model = fermiEdge[-1]
                else:          model += fermiEdge[-1]
        
        if("Exponential" in self.fitDict):
            exponential = self.fitDict["Exponential"]
            for idx,e in enumerate(exponential):
                A  = e[0]; Amin  = e[1]; Amax = e[2]
                c  = e[3]; cmin  = e[4]; cmax = e[5]
                P  = e[6]; Pmin  = e[7]; Pmax = e[8]
                
                pars.add('EXP' + str(idx) + '_A', value=A, min=Amin, max=Amax)
                pars.add('EXP' + str(idx) + '_c', value=c, min=cmin, max=cmax)
                pars.add('EXP' + str(idx) + '_P', value=P, min=Pmin, max=Pmax)
                
                if(not model): model  = Model(self.exponentialFunc,prefix='EXP' + str(idx) + '_')
                else:          model += Model(self.exponentialFunc,prefix='EXP' + str(idx) + '_')
                
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
        
        if("Point-Spectrum" in self.fitDict):
            pointSpec = self.fitDict["Point-Spectrum"]
            for idx,ps in enumerate(pointSpec):
                A = ps[0]; Amin = ps[1]; Amax = ps[2]
                c = ps[3]; cmin = ps[4]; cmax = ps[5]
                
                pars.add('PS' + str(idx) + '_A', value=A, min=Amin, max=Amax)
                pars.add('PS' + str(idx) + '_c', value=c, min=cmin, max=cmax)
                pars.add('PS' + str(idx) + '_datFile', value=idx, vary=False)
                
                if(not model): model  = Model(self.pointSpecFunc,prefix='PS' + str(idx) + '_')
                else:          model += Model(self.pointSpecFunc,prefix='PS' + str(idx) + '_')
                
        if(model == 0): return 0
        
        model.eval(pars,x=xx)
        return model.fit(yy,pars,x=xx)
    
    ###########################################################################
    # Form Actions (Show, Submit, Cancel, Remove, etc.)
    ###########################################################################
    def showForm(self,name):
        if(self.formActive): self.hideForm()
        self.formActive = True
        
        for l in self.forms[name]['labels']:
            l[0].grid(row=l[1],column=l[2],columnspan=1)
            l[0].configure(width=int(self.width/self.length),height=27)
            
        for idx,e in enumerate(self.forms[name]['entries']):
            e[0].grid(row=e[1],column=e[2],columnspan=1)
            e[0].configure(width=int(self.width/self.length),height=27)
            if(self.componentIdx < 0): continue
            e[0].delete(0,tk.END)
            e[0].insert(0,self.fitDict[name][self.componentIdx][idx])
        
        for idx,b in enumerate(self.forms[name]['buttons']):
            if(b[0].cget('text') == "remove" and self.componentIdx == -1): continue  # Only show the 'remove' button if we're editing a selection
            b[0].grid(row=b[1],column=b[2] + 2*idx,columnspan=2)
            b[0].configure(width=int(self.width/self.length),height=27)
            
    def hideForm(self,name=""):
        self.formActive = False
        names = [name]
        if(not names[0]): names = self.forms.keys()
        
        for name in names:
            for l in self.forms[name]['labels']:
                l[0].grid_forget()
                
            for e in self.forms[name]['entries']:
                e[0].grid_forget()
            
            for b in self.forms[name]['buttons']:
                b[0].grid_forget()
        
        self.btn['Add'].set("Add Fit")
        self.btn['Edit'].set("Edit Fit")
            
    def submitForm(self,name):
        filename = ""
        if(name == "Point-Spectrum"):
            filename = self._browseFile()
            if(not filename.endswith(".dat")):
                print("Excpecting .dat file")
                return
        
        params = []
        for e in self.forms[name]['entries']:
            try:
                params.append(float(e[0].get()))
            except:
                self.updateHelpLabel("Error in form: All values must be numeric.")
                return
        
        if(filename): params.append(filename)
        
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
        
    def editForm(self,name):
        if(name == "Edit Fit"): return
        name,index = name.split(" ")
        self.componentIdx = int(index)
        self.showForm(name=name)
        
    ###########################################################################
    # Custom Fitting Curves (not in lmfit)
    ###########################################################################
    def fermiDiracFunc(self,x, A, x0, T, Tb, c):
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
        return A/(1+np.exp((x0-x)/(-8.617e-5*T*Tb))) + c
    
    def referenceFunc(self,x,A,c):
        reference = self.mainPanel.stsPanel.getReferenceForCurve(x=x)
        return A*reference + c
    
    def pointSpecFunc(self,x,A,c,datFile):
        datFile = self.fitDict['Point-Spectrum'][datFile][6]
        ps = self.mainPanel.stsPanel.getDIDV(datFile=datFile)
        PS = self.mainPanel.stsPanel.getReferenceForCurve(x,reference=ps)
        return A*PS + c
    
    def exponentialFunc(self,x,A,c,P):
        return A*np.exp(x*P) + c
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
        if(name == "Add Fit"): return
        if(not name in self.fitDict):
            self.fitDict[name] = []
        self.showForm(name=name)
        
    def updateEditButton(self):
        editValues = ["Edit Fit"]
        for key in self.fitDict.keys():
            for idx,p in enumerate(self.fitDict[key]):
                editValues.append(key + " " + str(idx))
        self.btn['Edit'].configure(values=editValues,variable="Edit Fit")
        
    def reset(self):
        self.fitDict = {}
        self.update()