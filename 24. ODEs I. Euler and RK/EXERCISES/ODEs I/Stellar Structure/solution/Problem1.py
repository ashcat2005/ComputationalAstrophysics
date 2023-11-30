'''
Computational Astrophysics 
2019

Problem 1
Stellar Evolution
'''
import numpy as np
#import scipy as sp


# global constants (cgs units)
G = 6.67e-8
Msun  = 1.99e33

# EOS parameters for white dwarfs
Gamma = 4.0/3.0
K = 1.244e15*0.5**Gamma


#######################################
# function definitions
def ODE(radius, p, m): 
    # ODE RHS

    # invert EOS to get density (needed in rhs)
    rho = (p/K)**(1/Gamma)

    rhs = np.zeros(2)
    if(radius > 1.0e-10):
        rhs[0] = -G*m*rho/radius**2
        rhs[1] = 4*np.pi*rho*radius**2
    else:
        rhs[0] = 0.0
        rhs[1] = 0.0

    return rhs

def FEuler(radius, dr, p, m):
    # Forward-Euler Integrator

    new = np.zeros(2)
    old = np.zeros(2)
    old[0] = p
    old[1] = m

    # forward Euler integrator
    # for this, must call RHS routine
    new = old + dr*ODE(radius, p, m)
    
    # assign outputs
    pnew = new[0]
    mnew = new[1]
    
    return (pnew,mnew)

#######################################

# Set up grid
npoints = 1000
radmax = 2.0e8 # 2000 km
radius = np.linspace(0, radmax, npoints)
dr = radius[1]-radius[0]

# Set up variables
press = np.zeros(npoints)
rho   = np.zeros(npoints)
mass  = np.zeros(npoints)

# Set up central values (initial values)
rho[0]   = 1.0e10 #g/cm^3
press[0] = K * rho[0]**Gamma
mass[0]  = 0.0

# Set up termination criterion (surface pressure)
press_min = 1.0e-10 * press[0] 

nsurf = 0
for n in range(npoints-1):
    
    (press[n+1],mass[n+1]) = FEuler(radius[n],
                                    dr,
                                    press[n],
                                    mass[n])
    # check for termination criterion
    if(press[n+1] < press_min and nsurf==0):
        nsurf = n

    if(n+1 > nsurf and nsurf > 0):
        press[n+1] = press[nsurf]
        rho[n+1]   = rho[nsurf]
        mass[n+1]  = mass[nsurf]

    # invert the EOS to get density
    rho[n+1] = (press[n+1]/K)**(1/Gamma)

print(radius[nsurf]/1.0e5)
print(mass[nsurf]/Msun)

print(rho[nsurf])

dataout = np.stack((radius,rho), axis=1)
np.savetxt('rho.data', dataout)

# Plot the density function
import matplotlib.pyplot as plt 


plt.figure()
plt.plot(radius/1.0e5,rho)
plt.xlabel(r'$r$')
plt.ylabel(r'$\rho(r)$')
plt.show()
