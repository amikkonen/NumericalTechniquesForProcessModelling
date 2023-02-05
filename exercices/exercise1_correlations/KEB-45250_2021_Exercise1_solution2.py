# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 23:48:11 2017

@author: ami
"""
import numpy as np
from numpy import linalg

A = np.array([[3,2,0],
              [1,-1,0],
              [0,5,1]])
b = np.array([2,4,-1])

x = linalg.solve(A,b)
print("A", A)
print("b", b)
print("x", x)
    
A = np.array([[1,3,4],
              [7,2,3],
              [7,9,2],
                ])    
b = np.array([3,5,6])
x = linalg.solve(A,b)
print("A", A)
print("b", b)
print("x", x)
    
    


