# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 17:12:16 2022

@author: jced0001
"""

import tkinter as tk
import customtkinter as ctk
from MainPanel import MainPanel as mp

ctk.set_appearance_mode("Dark")  # Modes: system (default), light, dark
ctk.set_default_color_theme("blue")  # Themes: blue (default), dark-blue, green

# master    = tk.Tk()

class App(ctk.CTk):
    WIDTH = 512
    HEIGHT = 512
    def __init__(self):
        super().__init__()
        self.title("G80SXM")
        self.geometry(f"{App.WIDTH}x{App.HEIGHT}")
        self.protocol("WM_DELETE_WINDOW", self.on_closing)  # call .on_closing() when app gets closed
        
        dpi       = self.winfo_fpixels('1i')
        mainPanel = mp(self, width=512, height=512, dpi=dpi)

    def on_closing(self, event=0):
        self.destroy()
app = App()
app.mainloop()