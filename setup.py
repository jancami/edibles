#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst


from setuptools import setup, find_packages
import sys

with open("README.rst", 'r') as f:
    long_description = f.read()

min_version = (3, 8)

if sys.version_info < min_version:
    error = """
EDIBLES does not support Python {0}.{1}.
Python {2}.{3} or above is required. Check your Python version like so:

python3 --version

This may be due to an out-of-date pip. Make sure you have pip >= 9.0.1.
Upgrade pip like so:

pip install --upgrade pip
""".format(*(sys.version_info[:2] + min_version))
    sys.exit(error)

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
        'lmfit',
        'specutils'
    ],  # external packages as dependencies
)
