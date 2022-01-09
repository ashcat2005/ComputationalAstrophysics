'''
Eduard Larrañaga
Computational Astrophysics 
2020

Runge-Kutta 4
'''

def RK4(ODE, h, t0, q0):
    '''
    ------------------------------------------
    RK4(h, t0, q0)
    ------------------------------------------
    4th Order Runge-Kutta method for solving 
    a system of ODEs.
    Arguments:
    ODE: function defining the system of ODEs
    h: stepsize for the iteration
    t0: independent parameter initial value
    q0: numpy array with the initial values of
        the functions in the ODEs system
    ------------------------------------------
    '''
    k1 = h*ODE(t0, q0)
    k2 = h*ODE(t0 + h/2, q0 + k1/2)
    k3 = h*ODE(t0 + h/2, q0 + k2/2)
    k4 = h*ODE(t0 + h, q0 + k3)
    q1 = q0 + (k1 + 2*k2 + 2*k3 + k4)/6
    return q1