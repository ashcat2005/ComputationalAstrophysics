'''
Eduard Larrañaga
Computational Astrophysics 
2020

Lagrange Interpolation Method
'''

import numpy as np

#Lagrange Coefficients
def L(x, xi, j):
	'''
	------------------------------------------
	L(x, xi, j)
	------------------------------------------
    Returns the Lagrange coefficient for the 
    interpolation evaluated at points x
    Receives as arguments:
    x : array of points where the interpolated
    polynomial will be evaluated
    xi : array of N data points 
    j : index of the coefficient to be 
    calculated
	------------------------------------------
	'''
	# Number of points
	N = len(xi) 

	prod = 1
	for k in range(N):
		if (k != j):
			prod = prod * (x - xi[k])/(xi[j] - xi[k])
	return prod





# Interpolated Polynomial
def p(x, xi, fi):
	'''
	------------------------------------------
    p(x, xi, fi)
    ------------------------------------------
    Returns the values of the Lagrange 
    interpolated polynomial in a set of points
    defined by x
    x : array of points where the interpolated
    polynomial will be evaluated
    xi : array of N data points points
    fi : values of the function to be 
    interpolated
    ------------------------------------------
	'''
	# Number of points
	N = len(xi)

	summ = 0
	for j in range(N):
		summ = summ + fi[j]*L(x, xi, j)
	return summ
