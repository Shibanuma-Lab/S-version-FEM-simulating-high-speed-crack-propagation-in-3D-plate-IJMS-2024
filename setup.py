#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
from setuptools import setup, find_packages

def read_requirements():
    """Parse requirements from requirements.txt."""
    reqs_path = os.path.join('.', 'requirements.txt')
    with open(reqs_path, 'r') as f:
        requirements = [line.rstrip() for line in f]
    return requirements

__author__ = 'FuruhashiFumihito'

setup(
    name='arrest',  
    version='0.0',  
    description='This code is the simulation model of arrest',  
    author='Fumihito Furuhashi', 
    author_email='furuhashi-fumihito418@g.ecc.u-tokyo.ac.jp',  
    url='https://github.com/FuruhashiFumihito/arrest', 
    classifiers=[ 
        'Intended Audience :: Science/Research',
        'License :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3 :: Only',
    ],
    packages=['SourceDirectory'], 
    include_package_data=True,  
    license='MIT License', 
    install_requires=read_requirements(),
)