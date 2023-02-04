'''
Eduard Larrañaga
Computational Astrophysics 
2020

Hermite Interpolation Method
'''

import numpy as np

#Hermite Coefficients
def psi0(z):
	'''
	------------------------------------------
	psi0(z)
	------------------------------------------
	Returns the Hermite coefficients Psi_0
	for the interpolation
	Receives as arguments: z
	------------------------------------------
	'''
	psi_0 = 2*z**3 - 3*z**2 + 1
	return psi_0

def psi1(z):
	'''
	------------------------------------------
	psi1(z)
	------------------------------------------
	Returns the Hermite coefficients Psi_1 for 
	the interpolation
	Receives as arguments: z
	------------------------------------------
	'''
	psi_1 = z**3 - 2*z**2 + z
	return psi_1


# Interpolated Polynomial
def H3(x, xi, fi, dfidx):
	'''
	------------------------------------------
    H3(x, xi, fi, dfidx)
    ------------------------------------------
    Returns the values of the Cubic Hermite 
    interpolated polynomial in a set of points
    defined by x
    x : array of points where the interpolated
    polynomial will be evaluated
    xi : array of 2 data points 
    fi : array of values of the function at xi
    dfidx : array of values of the derivative 
    of the function at xi
    ------------------------------------------
	'''
	# variable z in the interpolation
	z = (x - xi[0])/(xi[1] - x[0])
	
	h1 = psi0(z) * fi[0]
	h2 = psi0(1-z)*fi[1]
	h3 = psi1(z)*(xi[1] - xi[0])*dfidx[0]
	h4 = psi1(1-z)*(xi[1] - xi[0])*dfidx[1]
	H =  h1 + h2 + h3 - h4
	return H
