'''
Eduard Larrañaga
Computational Astrophysics 
2020

Cramers Rule
'''
import numpy as np

def CramersRule(A, b):
    '''
    ------------------------------------------
    CramersRule(A,b)
    ------------------------------------------
    Returns the solution fo the Linear System
    A x = b
    where
    A: nxn matrix
    b: n vector
    
    Arguments:
    A: numpy array of size nxn
    b: numpy array of size n
    ------------------------------------------
    '''
    detA = np.linalg.det(A)
    n = len(b)
    # Create an empty vector for the solution
    x = np.zeros(n)
    # Create an empty matrix for the Cramer's rule
    Aj = np.zeros_like(A)
    
    #Main loop of Cramer's Rule
    for j in range(n):
        Aj[:] = A
        for i in range(n):
            Aj[i,j] = b[i]
        detAj = np.linalg.det(Aj)
        x[j] = detAj/detA
    return x