#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst


from setuptools import setup, find_packages


with open("README.rst", 'r') as f:
    long_description = f.read()

setup(
    name='EDIBLES',
    version='0.0.1',
    description='Software for EDIBLES data analysis and plotting',
    license="BSD",
    long_description=long_description,
    author='Jan Cami',
    author_email='jcami@uwo.ca',
    packages=find_packages(),  # same as name
    install_requires=[
        'numpy',
        'matplotlib',
        'pandas',
        'astropy',
        'scipy',
        'lmfit'
    ],  # external packages as dependencies
)
