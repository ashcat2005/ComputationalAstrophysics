import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time

L = 1.0  
kappa = 0.1
lambd = 2.01
nu = 0.9

Nx = 601  
x = np.linspace(-L, L, Nx)
dx = x[1] - x[0]
T_final = 10.0  
Nt = 5000      
dt = T_final / Nt

rho_initial = np.ones(Nx)
v_initial = np.zeros(Nx)   
p_initial = kappa*(rho_initial)**2

rho_history = [rho_initial.copy()]
v_history = [v_initial.copy()]
p_history = [p_initial.copy()]

def ddx(f):
    d = np.zeros_like(f)
    d[1:-1] = (f[2:] - f[:-2])/(2*dx)
    d[0] = (f[1] - f[0])/dx
    d[-1] = (f[-1] - f[-2])/dx
    return d

def lax_friedrichs(rho, v):
    """Método de Lax-Friedrichs."""
    F = rho * v
    alpha = dx / (2 * dt) 
    fluxes = np.zeros(Nx - 1)
    fluxes[:] = 0.5 * (F[:-1] + F[1:]) - alpha * (rho[1:] - rho[:-1])
    return fluxes

def dv_dt(v, rho):
    drho_dx = ddx(rho)
    return -2 * kappa * drho_dx - lambd * x - nu * v

def drho_dt(rho, v):
    fluxes = lax_friedrichs(rho, v)
    drho_dt = np.zeros_like(rho)
    drho_dt[1:-1] = -(fluxes[1:] - fluxes[:-1]) / dx
    drho_dt[0] = -(fluxes[0] - 0) / dx   #Con flujo distinto de cero
    drho_dt[-1] = -(0 - fluxes[-1]) / dx  #Con flujo distinto de cero
    # drho_dt[0] = 0  # Sin flujo en los bordes
    # drho_dt[-1] = 0 # Sin flujo en los bordes
    return drho_dt


def dRho_dt(v, rho):
    dv_dx = ddx(v)
    return - rho*dv_dx

rho_current = rho_initial.copy()
v_current = v_initial.copy()
p_current = p_initial.copy()

start_time = time.time()

# El loop principal es un RK4

for i in range(Nt):
    t_n = i * dt

    k1_v = dv_dt(v_current, rho_current)
    k2_v = dv_dt(v_current + 0.5 * dt * k1_v, rho_current)
    k3_v = dv_dt(v_current + 0.5 * dt * k2_v, rho_current)
    k4_v = dv_dt(v_current + dt * k3_v, rho_current)

    v_next = v_current + (dt / 6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v)

    #k1_r = dRho_dt(v_current, rho_current)
    #k2_r = dRho_dt(v_current, rho_current + 0.5 * dt * k1_r)
    #k3_r = dRho_dt(v_current, rho_current + 0.5 * dt * k2_r)
    #k4_r = dRho_dt(v_current, rho_current + dt * k3_r)

    #rho_next = rho_current + (dt / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r)

    drhodt = drho_dt(rho_current, v_current)
    rho_next = rho_current + dt * drhodt #Paso euler

    p_next = kappa*(rho_next)**2

    rho_current = rho_next
    v_current = v_next
    p_current = p_next

    rho_history.append(rho_current.copy())
    v_history.append(v_current.copy())
    p_history.append(p_current.copy())

    if i % 100 == 0 or i == Nt - 1:
        elapsed = time.time() - start_time
        frac   = (i + 1) / Nt
        eta    = elapsed * (1 - frac) / frac if frac > 0 else float('nan')
        print(f'\rProgreso: {frac*100:5.1f}%', end='')

# Animación

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7, 7), sharex = True)
line1, = ax1.plot(x, rho_history[0], label=r'$\rho(x, t)$')
line2, = ax2.plot(x, v_history[0], label=r'$v(x, t)$', color='orange')
line3, = ax3.plot(x, p_history[0])

ax1.set_xlim([-L-0.25, L+0.25])
ax1.set_ylim([-0.5, np.max(rho_history)+0.5]) 
ax1.set_xlabel('x')
ax1.set_ylabel(r'$\rho$')
ax1.set_title(r'Densidad')
ax1.grid(True)

ax2.set_xlim([-L-0.25, L+0.25])
ax2.set_ylim([-2.2, 2.2]) 
ax2.set_xlabel('x')
ax2.set_ylabel('v')
ax2.set_title(r'Velocidad')
ax2.grid(True)

ax3.set_xlim([-L-0.25, L+0.25])
ax3.set_ylim([-0.5, np.max(p_history)+0.5]) 
ax3.set_xlabel('x')
ax3.set_ylabel('P')
ax3.set_title(r'Presión')
ax3.grid(True)

fig.subplots_adjust(hspace=0.5)

time_text = ax1.text(0.02, 0.75, '', transform=ax1.transAxes)

def update(frame):
    current_time = frame * dt
    line1.set_ydata(rho_history[frame])
    line2.set_ydata(v_history[frame])
    line3.set_ydata(p_history[frame])
    time_text.set_text(f'Tiempo: {current_time:.3f} s')
    return line1, line2, line3, time_text

frame_indices = range(0, len(rho_history), 10)
ani = FuncAnimation(fig, update, frames=frame_indices, blit=True, interval=20)
# ani.save('Toy_Star.gif', writer='pillow', fps=15, dpi=80)


#plt.tight_layout()
plt.show()
