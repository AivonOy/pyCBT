# -*- coding: utf-8 -*-
"""
Created on Wed May 15 14:54:55 2013

@author: Leif Roschier
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May 15 14:29:28 2013

@author: Leif Roschier
"""
import pyximport
pyximport.install(setup_args={"script_args":["--compiler=mingw32"]}, reload_support=True)
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("CBT_lib", ["CBT_lib.pyx"])]
)