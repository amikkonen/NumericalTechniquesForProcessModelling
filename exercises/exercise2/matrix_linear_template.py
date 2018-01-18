#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A course example of steady one-dimensional heat conduction in a fin.


Created on Thu Jan  18 9:49 2018

@author: Antti Mikkonen, a.mikkonen@iki.fi
"""
import scipy as sp
from scipy import linalg
from matplotlib import pyplot as plt

##############################################################################

def main():
    # Source vector
    b = sp.array([1,1,1])
    # Coefficient matrix
    A = sp.array([[1,0,0],
                  [0,1,0],
                  [0,0,1]
                 ]
                )
    # Solution
    x = sp.linalg.solve(A,b)


    print("A\n",A)
    print("b\n",b)
    print("x\n",x)


##############################################################################

if __name__ == "__main__":
    print("START")
    main()
    print("END")
