#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

version = '0.1.0'

setup(
    name='G80SXM',
    version=version,
    author='Julian Ceddia',
    author_email='jdceddia@gmail.com',
    description='Data analysis tool for nanonis .sxm, .dat, and .3ds files.',
    long_description=long_description,
    url='https://github.com/New-Horizons-SPM/G80SXM',
    project_urls = {
        "Bug Tracker": "https://github.com/New-Horizons-SPM/G80SXM/issues"
    },
    license='MIT',
    packages=find_packages(),
    install_requires=['numpy', 'matplotlib', 'customtkinter','matplotlib_scalebar','tk','nanonispy','scipy','ase','lmfit'],
)
