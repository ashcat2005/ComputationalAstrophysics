import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from IPython.display import HTML


# In all objects, we will consider the first index as the x-component 
# and the second index as the y-component!


def gaussian2D(X, Y, x0=2, y0=.0):
    sigma2 = 0.5
    z = np.exp(-((X.T-x0)**2 + (Y.T-y0)**2)/(2*sigma2))
    return z

def initVel(X, Y):
    vx0 = 0.01*np.exp(-Y.T**2/0.25)
    vy0 = -0.001*Y.T
    return np.array([vx0, vy0])





# FTCS Method
def FTCS(init_profile, init_vel, tgrid, xgrid, ygrid):
    print('FTCS Method')
    dt = tgrid[1] - tgrid[0]
    dx = xgrid[1] - xgrid[0]
    dy = ygrid[1] - ygrid[0]
    print('dt = ', dt)
    print('dx = ', dx)
    print('dy = ', dy)

    rho = np.zeros([len(tgrid), len(xgrid), len(ygrid)])
    vx = np.zeros([len(tgrid), len(xgrid), len(ygrid)])
    vy = np.zeros([len(tgrid), len(xgrid), len(ygrid)])
    rho[0] = init_profile # initial condition
    vx[0] = init_vel[0]
    vy[0] = init_vel[1]
    for n in range(len(tgrid)-1):
        for i in range(1,len(xgrid)-1):
            for j in range(1,len(ygrid)-1):
                rho[n+1,i,j] = rho[n,i,j] \
                               - (dt/(2*dx))*vx[n,i,j]*(rho[n,i+1,j] - rho[n,i-1,j]) \
                               - (dt/(2*dy))*vy[n,i,j]*(rho[n,i,j+1] - rho[n,i,j-1])
                vx[n+1,i,j] = vx[n,i,j]  # No evolution of the velocity field
                vy[n+1,i,j] = vy[n,i,j]  # No evolution of the velocity field
        # Outflow boundary conditions (up and down)
        rho[n+1,:,0] = rho[n+1,:,1]
        rho[n+1,:,-1] = rho[n+1,:,-2]
        vx[n+1,:,0] = vx[n+1,:,1]
        vx[n+1,:,-1] = vx[n+1,:,-2]
        vy[n+1,:,0] = vy[n+1,:,1]
        vy[n+1,:,-1] = vy[n+1,:,-2]
        # Periodic boudnary conditions (left and right)
        rho[n+1,0,:] = rho[n+1,-2,:]
        rho[n+1,-1,:] = rho[n+1,-2,:]
        vx[n+1,0,:] = vx[n+1,-2,:]
        vx[n+1,-1,:] = vx[n+1,-2,:]
        vy[n+1,0,:] = vy[n+1,-2,:]
        vy[n+1,-1,:] = vy[n+1,-2,:]
    return rho, vx, vy



# Spatial grids
xgrid = np.linspace(0,10,100) 
ygrid = np.linspace(-1,1,20) 
X, Y = np.meshgrid(xgrid, ygrid)

# Initial Conditions
rho0 = gaussian2D(X,Y)
v0 = initVel(X,Y)




n = 2000 
tgrid = np.linspace(0,800,n)
rho, vx, vy = FTCS(rho0, v0, tgrid, xgrid, ygrid)



fig,ax = plt.subplots(2,1)
fstep = 5
nframe = 0
while nframe < n:
    ax[0].set_title('Density   Frame :'+str(nframe))
    ax[0].set_xlim(0,10)
    ax[0].set_ylim(-1,1)
    ax[0].set_xlabel(r'$x$')
    ax[0].set_ylabel(r'$y$')
    ax[0].imshow(rho[nframe].T, origin='lower', extent=[0,10,-1,1])
    #ax[0].streamplot(X,Y,vx[nframe].T,vy[nframe].T, density=.7, 
    #                linewidth=0.1, arrowsize=0.5, color='white')
    
    #ax[1].set_title('Vorticity')
    #ax[1].imshow((-vx[nframe].T)*Y + X*(vy[nframe].T), 
    #           origin='lower', extent=[0,10,-1,1], cmap='bwr')
    ax[1].set_title('Velocity')
    ax[1].imshow(np.log(np.sqrt((vx[nframe].T)**2+(vy[nframe].T)**2)), 
               origin='lower', extent=[0,10,-1,1], cmap='Reds')
    ax[1].set_xlabel(r'$x$')
    ax[1].set_ylabel(r'$y$')
    plt.pause(0.0001)
    ax[0].cla()
    ax[1].cla()
    nframe += fstep

print("Done!")
plt.show()




