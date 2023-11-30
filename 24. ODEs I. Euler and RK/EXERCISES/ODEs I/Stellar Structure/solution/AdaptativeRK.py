'''
Eduard Larrañaga
Computational Astrophysics 
2020

Runge-Kutta 3/4 with and adaptative sizestep
'''
import numpy as np

def RK34(ODE, h, t0, q0):
    '''
    ------------------------------------------
    RK34(h, t0, q0)
    ------------------------------------------
    Runge-Kutta method of order 3/4 for solving 
    a ODEs system using an adaptative sizestep.
    Arguments:
    h: stepsize for the iteration
    t0: independent parameter initial value
    q0: numpy array with the initial values of
        the functions in the ODEs system
    ------------------------------------------
    '''
    k1 = h*ODE(t0, q0)
    k2 = h*ODE(t0 + h/2, q0 + k1/2)
    k3 = h*ODE(t0 + 3*h/4, q0 + 3*k2/4)
    q1 = q0 + (2*k1 + 3*k2 + 4*k3)/9
    k4 = h*ODE(t0 + h, q1)
    q1s = q0 + (7*k1 + 6*k2 + 8*k3 +3*k4)/24
    return q1, q1s


def ARK(ODE, h, t0, q0):
    '''
    ------------------------------------------
    ARK(h, t0, q0)
    ------------------------------------------
    Adaptative Runge-Kutta method for solving 
    a ODEs system.
    Arguments:
    h: stepsize for the iteration
    t0: independent parameter initial value
    q0: numpy array with the initial values of
        the functions in the ODEs system
    
    The error tolerance is given by the 
    parameter epsilon and the fudge factor is
    given by the parameter S.
    ------------------------------------------
    '''
    S = 0.999
    epsilon = 1.E-11 
    q1, q1s = RK34(ODE, h, t0, q0)
    q_error = q1[1] - q1s[1]
    Delta = np.abs(q_error)/epsilon
    h = h*S/(Delta**(1/4))
    while Delta>1:
        q1, q1s = RK34(ODE, h, t0, q0)
        q_error = q1[1] - q1s[1]
        Delta = np.abs(q_error)/epsilon
        h = h*S/(Delta**(1/4))
    return h, q1