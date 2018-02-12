#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 14:18:02 2018

@author: ami
"""

import scipy as sp

import scipy.sparse as sparse
from scipy.sparse import linalg as splinalg

class SparseMatrixA(object):
    def __init__(self):
        self.row = []
        self.col = []
        self.val = []

    def add(self, row, col, val):
        self.row.append(row)
        self.col.append(col)
        self.val.append(val)

    def finalize(self):
        return sparse.csr_matrix(sparse.coo_matrix(
                               (self.val, [self.row, self.col])
                              ))
    
    def todense(self):
        return self.finalize().todense()
    
    def solve(self,b, A=None):
        if A is None:
            A = self.finalize()
        return splinalg.spsolve(A, b)
    
    
if __name__ == "__main__":
    import heatConduction1DConstantTemperatureBoundariesNoSources
    
    
    def solver_1d(n, kt, Ac, L, T_A, T_B):
    
        # Cell size
        dx = L / n
        # Coefficients
        aW = kt/dx*Ac
        aE = kt/dx*Ac
        
        # Source vector
        b = sp.zeros(n)
        # Coefficient matrix
#        A = sp.zeros((n,n))
        A = SparseMatrixA()
    
        # Add aW to matrix A
        for k in range(1,n):
#            A[k,k]   += aW
#            A[k,k-1] += -aW
            A.add(k,k, aW)
            A.add(k,k-1, -aW)

    #    print("A with aW\n", A)
        
        # Add aE to matrix A
        for k in range(n-1):
#            A[k,k]   += aE
#            A[k,k+1] += -aE
            A.add(k,k,aE)
            A.add(k,k+1, -aE)
    #    print("A with aE\n", A)
        
        # Boundary A on the left
#        A[0,0] += kt/(dx/2)*Ac
        A.add(0,0,kt/(dx/2)*Ac)
        b[0]   += kt/(dx/2)*Ac*T_A
        
        # Boundary B on the right
#        A[-1,-1] += kt/(dx/2)*Ac
        A.add(n-1,n-1, kt/(dx/2)*Ac)
        b[-1]    += kt/(dx/2)*Ac*T_B
        
        # Solution
#        T = sp.linalg.solve(A,b)
        T = A.solve(b)
        
        return T
    
    print("START")
    # Number of control volumes
    n = 5
    # Thermal conductivity (W/mK)
    kt = 1000
    # Cross-sectional area (m2)
    Ac = 10e-3
    # Lenght (m)
    L = 0.5
    # Boundary temperatures (C)
    T_A = 100
    T_B = 500
    
    Tfull, dxFull = heatConduction1DConstantTemperatureBoundariesNoSources.main(n, kt, Ac, L, 
                                                                T_A, T_B)
    
    Tsparse = solver_1d(n, kt, Ac, L, T_A, T_B)
    
    assert all(sp.isclose(Tfull, Tsparse))
    print("OK")
    print("END")
