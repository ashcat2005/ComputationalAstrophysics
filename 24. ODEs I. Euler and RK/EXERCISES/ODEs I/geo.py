
import numpy as np
import matplotlib.pyplot as plt


# Runge-Kutta 4 Algorithm
def RK4(ODE, tau0, q0, tauf, n):
    '''
    RK4
    '''
    dtau = (tauf - tau0)/(n-1)
    q = np.zeros([n,len(q0)+1])
    q[0,0] = tau0
    q[0,1:] = q0

    for i in range(1,n):
        q[i,0] = q[i-1,0] + dtau
        k1 = dtau*ODE(q[i-1,1:])
        k2 = dtau*ODE(q[i-1,1:] + k1/2)
        k3 = dtau*ODE(q[i-1,1:] + k2/2)
        k4 = dtau*ODE(q[i-1,1:] + k3)
        q[i,1:] = q[i-1,1:] + (k1 + 2*k2 + 2*k3 + k4)/6
    
    return q



def geodesics(x):
    '''
    This function contains the geodesic equations in the Schwarzschild metric.
    '''

    # Coordinates and velocity components
    t = x[0]
    r = x[1]
    theta = x[2]
    phi = x[3]
    omega = x[4]
    rho = x[5]
    mu = x[6]
    nu = x[7]

    # Geodesics differential equations 
    dtdtau = omega
    drdtau = rho
    dthdtau = mu
    dphidtau = nu
    
    domegadtau = -2/(r*(r-2))*omega*rho
    drhodtau =  rho**2/(r*(r-2)) + (r-2)*(mu**2 + np.sin(theta)**2*nu**2 - omega**2/r**3)
    dmutau = -2*rho*mu/r + np.sin(theta)*np.cos(theta)*nu**2
    dnudtau = -2*np.cos(theta)*mu*nu/np.sin(theta) -2*rho*nu/r
    
    dxdtau = np.array([dtdtau, drdtau, dthdtau, dphidtau, 
              domegadtau, drhodtau, dmutau, dnudtau])
    return dxdtau



# Malla para integración
tau0 = 0.
tauf = 2000.
n = 50000

# Condiciones iniciales
t0 = 0.
r0 = 15.
theta0 = np.pi/2
phi0 = 0.

L = 5. # Momento angular
E = np.sqrt((1-2/r0)*(1+L**2/r0**2)) # Energia con el Potencial efectivo

omega0 = E*r0/(r0-2)
mu0 = 0.
nu0 = L/(r0**2*np.sin(theta0)**2)

#A1 =(1 - 2/r0)*(-1 + (1-2/r0)*omega0**2 - r0**2*mu0**2 - r0**2*np.sin(theta0)**2*nu0**2) 
A1 = (1 - 2/r0)*(-1 + r0*E**2/(r0-2) - r0**2*mu0**2 - L**2/(r0**2*np.sin(theta0)**2))
rho0 = np.sqrt(A1)

q0 = np.array([t0, r0, theta0, phi0, omega0, rho0, mu0, nu0])


# Llamado a RK4
Q = RK4(geodesics, tau0, q0, tauf, n)





x = Q[:,2]*np.cos(Q[:,4])
y = Q[:,2]*np.sin(Q[:,4])

plt.figure()
plt.plot(x,y)
plt.xlabel('x')
plt.ylabel('y')
plt.show()






