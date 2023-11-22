import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from numpy import sqrt, pi, exp, sum, hstack
import sys
import time

"""
Structure of a star with SPH
"""

def W(x, y, z, h):
	"""
    3D Gausssian Smoothing kernel
	x : vector/matrix of x positions
	y : vector/matrix of y positions
	z : vector/matrix of z positions
	h : smoothing length
	w : evaluated smoothing function
	"""
	r = sqrt(x**2 + y**2 + z**2)
	w = (1. / (h*sqrt(pi)))**3 * exp(-(r/h)**2)
	return w
	
	
def gradW(x, y, z, h):
    """
    Gradient of the 3D Gausssian Smoothing kernel
    x : vector/matrix of x positions
	y : vector/matrix of y positions
	z : vector/matrix of z positions
	h : smoothing length
    wx, wy, wz : evaluated gradient
    """
    r = sqrt(x**2 + y**2 + z**2)
    dw = -2 * exp(-(r/h)**2) / (h**5 * (pi)**(3/2))
    dWx = dw*x
    dWy = dw*y
    dWz = dw*z
    return dWx, dWy, dWz
	
	
def PairwiseSeparations(ri, rj):
	"""
	Get pairwise separations between 2 sets of coordinates
	ri : M x 3 matrix of positions
	rj : N x 3 matrix of positions
	dx, dy, dz : M x N matrices of separations
	"""
	
	M = ri.shape[0]
	N = rj.shape[0]

	# positions ri = (x,y,z)
	rix = ri[:,0].reshape((M,1))
	riy = ri[:,1].reshape((M,1))
	riz = ri[:,2].reshape((M,1))

	# other set of points positions rj = (x,y,z)
	rjx = rj[:,0].reshape((1,N))
	rjy = rj[:,1].reshape((1,N))
	rjz = rj[:,2].reshape((1,N))
	
	# matrices that store all pairwise particle separations: r_i - r_j
	dx = rix - rjx
	dy = riy - rjy
	dz = riz - rjz
	
	return dx, dy, dz
	

def Density( r, pos, m, h ):
	"""
	Get Density at sampling loctions from SPH particle distribution
	r   : M x 3 matrix of sampling locations
	pos : N x 3 matrix of SPH particle positions
	m   : particle mass
	h   : smoothing length
	rho : M x 1 vector of accelerations
	"""
	M = r.shape[0]
	dx, dy, dz = PairwiseSeparations( r, pos )
	rho = sum( m * W(dx, dy, dz, h), 1 ).reshape((M,1))
	return rho
	
	
def Pressure(rho, k, n):
	"""
	Equation of State
	rho : vector of densities
	k   : equation of state constant
	n   : polytropic index
	P   : pressure
	"""
	P = k * rho**(1+1/n)
	return P
	

def Acceleration( pos, vel, m, h, k, n, lmbda, nu ):
	"""
	Calculate the acceleration on each SPH particle
	[ Euler Equation ]
	pos   : N x 3 matrix of positions
	vel   : N x 3 matrix of velocities
	m     : particle mass
	h     : smoothing length
	k     : equation of state constant
	n     : polytropic index
	lmbda : external force constant
	nu    : viscosity
	a     : N x 3 matrix of accelerations
	"""
	N = pos.shape[0]
	
	# Calculate densities at the position of the particles
	rho = Density( pos, pos, m, h )
	
	# Get the pressures
	P = Pressure(rho, k, n)
	
	# Get pairwise distances and gradients
	dx, dy, dz = PairwiseSeparations( pos, pos )
	
	dWx, dWy, dWz = gradW( dx, dy, dz, h )
	
	# Add Pressure contribution to accelerations
	ax = - sum( m * ( P/rho**2 + P.T/rho.T**2  ) * dWx, 1).reshape((N,1))
	ay = - sum( m * ( P/rho**2 + P.T/rho.T**2  ) * dWy, 1).reshape((N,1))
	az = - sum( m * ( P/rho**2 + P.T/rho.T**2  ) * dWz, 1).reshape((N,1))
	
	# pack together the acceleration components
	a = hstack((ax,ay,az))

	# Add gravtiational force
	a += gravForce(lmbda, pos)
	
	# Add viscosity force
	a += viscosForce(nu,vel)
	
	return a
	

def gravForce(lmbda, r):
	'''
	Simple Gravitational Force Model
	representing a constant gravitational potential
	'''
	return -lmbda*r


def viscosForce(nu, v):
	'''
	Viscosity force
	'''
	return -nu*v


def RealTimePlot( data, rho_data, density_data, rho_analytic ):
	fig = plt.figure()
	grid = plt.GridSpec(3, 1, wspace=0.0, hspace=0.3)
	ax1 = plt.subplot(grid[0:2,0])
	ax2 = plt.subplot(grid[2,0])
	rlin = np.linspace(0,1,100)

	for i in range(data.shape[0]):
		plt.sca(ax1)
		plt.cla()
		cval = np.minimum((rho_data[i]-3)/3,1).flatten()
		plt.scatter(data[i,:,0],data[i,:,1], c=cval, cmap=plt.cm.autumn, s=10, alpha=0.5)
		
		ax1.set(xlim=(-1.4, 1.4), ylim=(-1.2, 1.2))
		ax1.set_aspect('equal', 'box')
		ax1.set_xticks([-1,0,1])
		ax1.set_yticks([-1,0,1])
		ax1.set_facecolor('black')
		ax1.set_facecolor((.1,.1,.1))
			
		plt.sca(ax2)
		plt.cla()
		ax2.set(xlim=(0, 1), ylim=(0, 3))
		ax2.set_aspect(0.1)
		plt.plot(rlin, rho_analytic, color='gray', linewidth=2)
		plt.plot(rlin, density_data[i], color='blue')
		plt.xlabel('radius')
		plt.ylabel('density')

		plt.pause(0.00001)
	
	print()
	print('Done!')
	# Save figure
	#plt.savefig('sph.png',dpi=240)
	plt.show()




def main():
	""" N-body simulation """
	
	# Simulation parameters
	N     = 1000    # Number of particles
	t     = 0      # current time of the simulation
	tEnd  = 12     # time at which simulation ends
	dt    = 0.04   # timestep
	M     = 2      # star mass
	R     = 0.75   # star radius
	h     = 0.1    # smoothing length
	k     = 0.1    # equation of state constant
	n     = 1      # polytropic index
	nu    = 2      # damping
	
	# Generate Initial Conditions
	np.random.seed(413)            # set the random number generator seed
	
	lmbda = 2*k*(1+n)*np.pi**(-3/(2*n)) * (M*gamma(5/2+n)/R**3/gamma(1+n))**(1/n) / R**2  # ~ 2.01
	m     = M/N                    # single particle mass
	pos   = np.random.randn(N,3)   # randomly selected positions and velocities
	vel   = np.zeros(pos.shape)
	#print('lambda = ', lmbda)
	#print('m = ', m)

	# calculate initial gravitational accelerations
	acc = Acceleration( pos, vel, m, h, k, n, lmbda, nu )

	# number of timesteps
	Nt = int(np.ceil(tEnd/dt))
	
	rr = np.zeros((100,3))
	rlin = np.linspace(0,1,100)
	rr[:,0] = rlin
	rho_analytic = lmbda/(4*k) * (R**2 - rlin**2)
	
	data = np.zeros([Nt, N, 3])
	rho_data = np.zeros([Nt, N])
	density_data = np.zeros([Nt, 100])

	# Simulation Main Loop
	for i in range(Nt):
		# (1/2) kick
		vel += acc * dt/2
		# drift
		pos += vel * dt
		# update accelerations
		acc = Acceleration( pos, vel, m, h, k, n, lmbda, nu )
		# (1/2) kick
		vel += acc * dt/2
		
		# update time
		t += dt

		rho = Density( pos, pos, m, h )
		data[i,:,:] = pos
		rho_data[i] = rho.reshape(N,)

		density_data[i,:] = Density( rr, pos, m, h ).reshape((100,))
		sys.stdout.write("\rStep # %d" %i)
		sys.stdout.flush()
	    
	return data, rho_data, density_data, rho_analytic


	
  
if __name__== "__main__":
	print('\n Starting SPH Simulation \n')
	start_time = time.time()
	data, rho_data, density_data, rho_analytic = main()
	print()
	print("--- Execution Time : %s seconds ---" % (time.time() - start_time))
	np.save('data/data.npy', data)
	np.save('data/rho_data.npy', rho_data)
	np.save('data/density_data.npy', density_data)
	np.save('data/rho_analytic.npy', rho_analytic)
	RealTimePlot(data, rho_data, density_data, rho_analytic)
