import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

"""
Structure of a star with SPH
"""

def W(x, y, z, h):
	"""
    3D Gausssian Smoothing kernel
	x     is a vector/matrix of x positions
	y     is a vector/matrix of y positions
	z     is a vector/matrix of z positions
	h     is the smoothing length
	w     is the evaluated smoothing function
	"""
	r = np.sqrt(x**2 + y**2 + z**2)
	w = (1.0 / (h*np.sqrt(np.pi)))**3 * np.exp( -r**2 / h**2)
	return w
	
	
def gradW(x, y, z, h):
    """
    Gradient of the 3D Gausssian Smoothing kernel
    x     is a vector/matrix of x positions
    y     is a vector/matrix of y positions
    z     is a vector/matrix of z positions
    h     is the smoothing length
    wx, wy, wz     is the evaluated gradient
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    n = -2 * np.exp( -r**2 / h**2) / (h**5 * (np.pi)**(3/2))
    dWx = n * x
    dWy = n * y
    dWz = n * z
    return dWx, dWy, dWz
	
	
def PairwiseSeparations(ri, rj):
	"""
	Get pairwise desprations between 2 sets of coordinates
	ri    is an M x 3 matrix of positions
	rj    is an N x 3 matrix of positions
	dx, dy, dz   are M x N matrices of separations
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
	r     is an M x 3 matrix of sampling locations
	pos   is an N x 3 matrix of SPH particle positions
	m     is the particle mass
	h     is the smoothing length
	rho   is M x 1 vector of accelerations
	"""
	M = r.shape[0]
	dx, dy, dz = PairwiseSeparations( r, pos )
	rho = np.sum( m * W(dx, dy, dz, h), 1 ).reshape((M,1))
	return rho
	
	
def Pressure(rho, k, n):
	"""
	Equation of State
	rho   vector of densities
	k     equation of state constant
	n     polytropic index
	P     pressure
	"""
	P = k * rho**(1+1/n)
	return P
	

def Acceleration( pos, vel, m, h, k, n, lmbda, nu ):
	"""
	Calculate the acceleration on each SPH particle
	pos   is an N x 3 matrix of positions
	vel   is an N x 3 matrix of velocities
	m     is the particle mass
	h     is the smoothing length
	k     equation of state constant
	n     polytropic index
	lmbda external force constant
	nu    viscosity
	a     is N x 3 matrix of accelerations
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
	ax = - np.sum( m * ( P/rho**2 + P.T/rho.T**2  ) * dWx, 1).reshape((N,1))
	ay = - np.sum( m * ( P/rho**2 + P.T/rho.T**2  ) * dWy, 1).reshape((N,1))
	az = - np.sum( m * ( P/rho**2 + P.T/rho.T**2  ) * dWz, 1).reshape((N,1))
	
	# pack together the acceleration components
	a = np.hstack((ax,ay,az))
	
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



def main():
	""" N-body simulation """
	
	# Simulation parameters
	N         = 1000    # Number of particles
	t         = 0      # current time of the simulation
	tEnd      = 12     # time at which simulation ends
	dt        = 0.04   # timestep
	M         = 2      # star mass
	R         = 0.75   # star radius
	h         = 0.1    # smoothing length
	k         = 0.1    # equation of state constant
	n         = 1      # polytropic index
	nu        = 1      # damping
	plotRealTime = True # switch on for plotting as the simulation goes along
	
	# Generate Initial Conditions
	np.random.seed(413)            # set the random number generator seed
	
	lmbda = 2*k*(1+n)*np.pi**(-3/(2*n)) * (M*gamma(5/2+n)/R**3/gamma(1+n))**(1/n) / R**2  # ~ 2.01
	m     = M/N                    # single particle mass
	pos   = np.random.randn(N,3)   # randomly selected positions and velocities
	vel   = np.zeros(pos.shape)
	
	# calculate initial gravitational accelerations
	acc = Acceleration( pos, vel, m, h, k, n, lmbda, nu )
	
	# number of timesteps
	Nt = int(np.ceil(tEnd/dt))
	
	# prep figure
	fig = plt.figure()
	grid = plt.GridSpec(3, 1, wspace=0.0, hspace=0.3)
	ax1 = plt.subplot(grid[0:2,0])
	ax2 = plt.subplot(grid[2,0])
	rr = np.zeros((100,3))
	rlin = np.linspace(0,1,100)
	rr[:,0] = rlin
	rho_analytic = lmbda/(4*k) * (R**2 - rlin**2)
	
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
		
		# get density for plottiny
		rho = Density( pos, pos, m, h )
		
		# plot in real time
		if plotRealTime or (i == Nt-1):
			plt.sca(ax1)
			plt.cla()
			cval = np.minimum((rho-3)/3,1).flatten()
			plt.scatter(pos[:,0],pos[:,1], c=cval, cmap=plt.cm.autumn, s=10, alpha=0.5)
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
			rho_radial = Density( rr, pos, m, h )
			plt.plot(rlin, rho_radial, color='blue')
			plt.pause(0.0001)
	    
	
	
	# add labels/legend
	plt.sca(ax2)
	plt.xlabel('radius')
	plt.ylabel('density')

	print('Done!')
	# Save figure
	plt.savefig('sph.png',dpi=240)
	plt.show()
	    
	return 0
	


  
if __name__== "__main__":
  main()
