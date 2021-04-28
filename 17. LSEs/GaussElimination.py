'''
Eduard Larrañaga
Computational Astrophysics 
2020

Gauss Elimination
'''
def GaussElim(A,b):
    '''
    ------------------------------------------
    GaussElim(A,b)
    ------------------------------------------
    Returns an upper-diagonal Linear System
    the solution fo the Linear System
    A x = b
    where
    A: nxn upper diagonal matrix
    b: n vector
    
    Arguments:
    A: upper-diagonal numpy array of size nxn
    b: numpy array of size n
    ------------------------------------------
    '''
    n = len(b)
    # Check that the pivots are not zero
    for k in range(n):
        if(A[k,k]==0):
            print(f'Pivot A[{k+1:d}, {k+1:d}] is zero')
            return None, None
        
    # Main Loop of the Gauss Elimination
    for i in range(n-1):
        for j in range(i+1,n):
            C=A[j,i]/A[i,i]
            for k in range(n):
                A[j,k] = A[j,k] - A[i,k]*C
            b[j] = b[j] - C*b[i]
    return A, b