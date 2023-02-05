import numpy as np
from numpy import linalg


A = np.array([[3.,2,0],
              [1,-1,0],
              [0,5,1],
                ])
b = np.array([2,4,-1])


# A = np.array([[3,2,0],
#               [1,-1,0],
#               [0,5,1]])
# b = np.array([2,4,-1])

x = linalg.solve(A,b)
print("A", A)
print("b", b)
print("x", x)
    

