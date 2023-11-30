'''
Plot of the analytic gravitatinal potential for a
uniform sphere inside and outside. The surface of the 
sphere is located at the radius R.
All values are considered in cgs untis
'''
import numpy as np
import matplotlib.pyplot as plt

# Constants of the problem
G = 6.67E-8 				# Newtonian Grav. Constant (cgs units)
rho = 1. 					# Uniform density of the sphere (gr/cm^3)
R = 1E9 					# Radius of the sphere (cm)
M = (4/3)*np.pi*rho*R**3 	# Total Mass (gr)

def potential_in(r):
	# Interior Potential
	Phi = (2.*np.pi*G/3)*(r**2 - 3*R**2)
	return Phi

def potential_out(r): 
	# Exterior Potential
	Phi = -G*M/r
	return Phi

# Variables for plotting
r_range_in = np.arange(0,1E9,1E5)
Phi_in = potential_in(r_range_in)
r_range_out = np.arange(1E9, 1E10, 1E5)
Phi_out = potential_out(r_range_out)


# Plot
plt.plot(r_range_in, Phi_in)
plt.plot(r_range_out, Phi_out)
plt.xlabel(r'$r$')
plt.ylabel(r'$\phi$')
plt.savefig('GravitationalPotential.pdf')
plt.show()
