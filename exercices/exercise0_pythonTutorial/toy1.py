#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 13:19:44 2021

@author: ami
"""

import numpy as np

rho = 1000
L   = 10
d   = 0.05
nu  = 1e-6
V   = 10

Re = V*d/nu
f  = (-1.8*np.log10(6.9/Re))**-2
dp = 0.5*rho*V**2*L/d*f

print("Re =", Re)
print("f  =", f)
print("dp =", dp, "Pa")




