# G80SXM
G80SXM is a tool written in python to quickly analyse and create figures for Nanonis SPM data (STM/STS/nc-AFM).
I am actively developing this tool so if you have any suggestions or find any bugs, feel free to create an issue

# Installation
clone the repository and run 'pip install .' from the root directory

To run: python gui.py

# Usage
There are a number of panels in the gui, each with a specific function. See below figures for some examples. Hovering over buttons in the GUI will also provide help

![nc-AFM_with_STS](./Documentation/nc-AFM_with_STS_and_molecules.png)
* Easily overlay atoms with correct scaling from custom .xyz files
* Easily add, position, and style inset figures
* Automatically or manually locate markers on the main figure that correspond to point spectra locations

![STM_with_Grid_and_STS](./Documentation/TOPO_with_Grid_and_AveragedSTS.png)
* Visualise grid data from .3ds files
* Easily pull out individual curves from a grid
* Ability to average curves from several pixels within a grid

![STM_with_STS_and_CurveFitting](./Documentation/TOPO_with_STS_and_CurveFitting.png)
* Quickly fit point spectra with
  * Reference curves/STS curves
  * Gaussian
  * Fermi-Dirac
  * More
