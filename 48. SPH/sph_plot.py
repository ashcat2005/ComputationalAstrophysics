import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from numpy import sqrt, pi, exp


"""
Animation of the Structure of a star with SPH from the data
"""

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
		#rho_radial = Density( rr, pos, m, h )
		plt.plot(rlin, density_data[i], color='blue')
		plt.xlabel('radius')
		plt.ylabel('density')

		plt.pause(0.00001)
	
	print()
	print('Done!')
	# Save figure
	#plt.savefig('sph.png',dpi=240)
	plt.show()



  
if __name__== "__main__":
	data = np.load('data.npz')
	rho_data = np.load('rho_data.npz')
	density_data = np.load('density_data.npz')
	rho_analytic = np.load('rho_analytic.npz')
	RealTimePlot(data, rho_data, density_data, rho_analytic)
