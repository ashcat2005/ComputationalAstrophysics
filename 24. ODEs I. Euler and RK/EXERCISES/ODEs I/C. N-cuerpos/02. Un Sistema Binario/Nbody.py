import numpy as np
# Computation time
import time

# Relevant Constants in the units
# (years, AU, Solar_masses)
# Newtonian Gravitational Constant
G = 4.*np.pi**2

# Conversion factors
m_Sun = 1.98855e30 # Solar mass in kg
au_in_meters = 1.49598261e11 # 1 au in meters
year_in_seconds = 3600*24*365 # 1 year in seconds

def a(t0, qn):
    '''
    ------------------------------------------
    a(t,qn) 
    ------------------------------------------
    ODEs system for the motion of N-particles 
    ------------------------------------------
    Arguments:
    t0: time parameter (not necessary for 
        the Newtonian problem)
    qn: numpy array with the position data for
        the N-particles at time n
        qn[0] = particle 1
        qn[1] = particle 2
        etc.
        qn[0] = [x0, y0, z0]
    mass: masses of the particles
    ------------------------------------------
    Returns: Acceleration of the N-particles
             at time n
    ------------------------------------------
    Dependences : NumPy
    ------------------------------------------
    '''
    f = np.zeros([N,3])
    
    # Gravitational Force
    for i in range(0,N):
        Deltaxyz = qn[i,0:3] - qn[:,0:3]
        # Distance between particles
        r = np.sqrt(np.sum(Deltaxyz*Deltaxyz,axis=1))
        # To avoid divivision by zero in the self-force terms
        # The terms vanish due to the Delta term in the numerator! 
        r[i] = 1 
        f[i,0] = -G*np.sum(Deltaxyz[:,0]*mass/(r**3))
        f[i,1] = -G*np.sum(Deltaxyz[:,1]*mass/(r**3))
        f[i,2] = -G*np.sum(Deltaxyz[:,2]*mass/(r**3))
    return f


def velVerlet(a, t_grid, q0, v0):
    '''
    ------------------------------------------
    velVerlet(h, t0, q0)
    Velocity-Verlet method for solving 
    a system of ODEs.
    ------------------------------------------
    Arguments:
    a: function defining the RHS of ODEs
    t_grid: array with the times at which the 
            solution will be computed
    q0: numpy array with the initial values of
        the positions
    v0: numpy array with the initial values of 
        the velocities
    ------------------------------------------
    Dependences: NumPy
    ------------------------------------------
    '''
    n = len(t_grid) # Number of steps
    dt = t_grid[1] - t_grid[0]
    # Initial condition
    q = np.zeros([n,N,3])
    v = np.zeros([n,N,3])
    energy = np.zeros(n)
    q[0] = np.transpose(q0)
    v[0] = np.transpose(v0)
    energy[0] = TotalEnergy(q[0],v[0])

    for i in range(1,n):
        v_half = v[i-1] + a(t_grid[i-1], q[i-1])*dt/2
        q[i] = q[i-1] + v_half*dt
        v[i] = v_half + a(t_grid[i], q[i])*dt/2
        energy[i] = TotalEnergy(q[i],v[i])
    return q,v,energy


def TotalEnergy(qn,vn):
    '''
    Total Energy calculation
    '''
    (x,y,z) = qn.transpose()
    (vx,vy,vz) = vn.transpose()
    v2 = vx**2+vy**2+vz**2
    Ekin = 0.5*np.sum(mass*v2)
    Egrav = 0.
    for i in range(0,N):
        deltax = x[i] - x
        deltay = y[i] - y
        deltaz = z[i] - z
        r = np.sqrt(deltax**2 + deltay**2 + deltaz**2)
        # To avoid divivision by zero and put the value of the term to zero!
        r[i] = 1e300 
        Egrav += - 0.5*G*mass[i]*np.sum(mass/r)
    return Ekin + Egrav



def NBodyplot(Q, save_opt=False):
    '''
    ------------------------------------------
    NBodyplot(Q)
    ------------------------------------------
    Plots the trajectories of the N-bodies in
    a 3D graph
    ------------------------------------------
    Arguments:
    Q: Array with the information of positions 
       for the N-body system
    ------------------------------------------
    Dependences: matplotlib
    ------------------------------------------
    '''
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d as Axes3D
    # Limits for the plot
    boundary = max(abs(np.amax(Q[:,:,0:3])),abs(np.amin(Q[:,:,0:3])))*1.1
    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(projection='3d')
    for i in range(len(mass)):
    	ax.plot(Q[:,i,0],Q[:,i,1],Q[:,i,2], label=f'Body {i+1:.0f}')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim3d(-boundary, boundary)
    ax.set_ylim3d(-boundary, boundary)
    ax.set_zlim3d(-boundary, boundary)
    ax.legend()
    if save_opt==True:
    	plt.savefig('NBody-output.jpg')
    else:
    	plt.show()



###################################################



# Read the initial data
#initial_data_file = "BinarySystemData.asc"
initial_data_file = "2BinarySystemData.asc"

(x,y,z,vx,vy,vz,mass) = np.loadtxt(initial_data_file, unpack = True)

# Convert from SI units to (years, AU, Solar_Mass) units
Q0 = np.array([x,y,z])/au_in_meters
V0 = np.array([vx,vy,vz])*year_in_seconds/au_in_meters
mass = mass/m_Sun # Global Variable

# Number of particles
N = len(mass) # Global Variable

# Creation of the time grid (in years)
t0 = 0.
tf = 5.
n = 40000 # Number of steps in the grid
t_grid = np.linspace(t0, tf, n)
dt = (tf - t0)/n # Constant stepsize


start = time.time()
Q,V,energy = velVerlet(a, t_grid, Q0, V0)  
end = time.time()

print('\n\nEl tiempo de computo fue:', end - start)

energychange = (energy[n-1]-energy[0])/energy[0]
print(f'\n\nEl cambio relativo de energía es {energychange:.5E}% con un tamaño de paso dt = {dt:.1E}\n\n')

NBodyplot(Q)