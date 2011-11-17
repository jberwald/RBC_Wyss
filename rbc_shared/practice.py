import numpy as np
from pylab import *

A=np.arange(9)
#print A

A.resize([3,3])
#print A

B=np.arange(3)
B=np.matrix(B)

C=A*B.T
#print C

#print "B.T",B.T

#print "B",B

# for i in range(3):
#     print A[i]

# for i in range(3):
#     print A[:,i]

def func( mat1,mat2 ):
    """comment"""
    C=A*B.T
    return C

D=func(A,B)
print D

def plot_result( mat1,mat2 ):
    r=func(mat1,mat2)
    plot(r)
    show()

plot_result(A,B)

    

    
