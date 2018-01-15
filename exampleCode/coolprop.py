#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 16:08:36 2018

@author: ami
"""

from CoolProp.CoolProp import PropsSI

# Density
rho = PropsSI("D","T", 300, "P", 1e5, "air")
print(rho)

