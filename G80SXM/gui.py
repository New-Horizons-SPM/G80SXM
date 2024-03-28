# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 17:12:16 2022

@author: jced0001
"""

import customtkinter as ctk
from .MainPanel import MainPanel as mp

ctk.set_appearance_mode("Dark")                                                 # Modes: system (default), light, dark
ctk.set_default_color_theme("blue")                                             # Themes: blue (default), dark-blue, green

class App(ctk.CTk):
    WIDTH = 512
    HEIGHT = 512
    def __init__(self):
        super().__init__()
        self.title("G80SXM")
        self.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        dpi = self.winfo_fpixels('1i')
        self.mainPanel = mp(self, width=self.WIDTH, height=self.HEIGHT, dpi=dpi)
        self.geometry("%dx%d" % (self.WIDTH, 800))
    def on_closing(self, event=0):
        self.mainPanel.quit()

def main():
    app = App()
    app.mainloop()