# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 17:12:16 2022

@author: jced0001
"""

import tkinter as tk
from MainPanel import MainPanel as mp

master    = tk.Tk()
dpi       = master.winfo_fpixels('1i')
mainPanel = mp(master, width=512, height=512, dpi=dpi)

master.mainloop()