# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='edibles',
    version='0.1.0',
    description='Software for EDIBLES data analysis and plotting.',
    long_description=readme,
    author='Jan Cami',
    author_email='jcami@uwo.ca',
    url='https://github.com/jancami/edibles',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)