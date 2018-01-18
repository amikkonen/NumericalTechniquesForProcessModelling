#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A simple example of linear algebra.

Created on Thu Jan  18 9:49 2018

@author: Antti Mikkonen, a.mikkonen@iki.fi
"""
import scipy as sp
from scipy import linalg
from matplotlib import pyplot as plt

##############################################################################

def main():
    # Source vector
    b = sp.array([4,3,1])
    # Coefficient matrix
    A = sp.array([[2,8,0],
                  [0,3,3],
                  [7,0,8]
                 ]
                )
    # Solution
    x = sp.linalg.solve(A,b)

    # Post
    print("A\n",A)
    print("b\n",b)
    print("x\n",x)


##############################################################################

if __name__ == "__main__":
    print("START")
    main()
    print("END")
