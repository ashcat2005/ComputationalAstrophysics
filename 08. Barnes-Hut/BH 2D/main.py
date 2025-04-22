from common import *
from numpy import concatenate
import time


######### GENERAL SETTINGS ###############################

# Mass of the N bodies.
max_mass = 100. # Solar masses

# Number of time-iterations executed by the program.
n = 2000

# Frequency at which .PNG images are written.
img_step = 10

# Name of the generated video
video_name = 'E1.mp4'

 
######### MODEL OPTIONS ##################################

# Model 1: SMBH + Keplerian velocity
# Supermassive Central Black Hole data
BHM = 4.e6 # [Solar masses] # Black Hole Mass
center = array([0.5, 0.5]) #[kpc] # Location of the SBH
BHmomentum = array([0.,0.]) # Momentum of the SBH
N = 500 # Number of bodies
R = 30. #[kpc]   # Initial radius of the distribution


# Model 2: Exponential Surface Density
Mtotal = 1.e4 #[Solar masses] # Total mass of the galaxy
center = array([0.5, 0.5]) #[kpc]
R = 30. #[kpc] # Galactic Radius
Rd = 0.1*R #[kpc] # Length Scale


# Model 3: Galaxy Model
# Supermassive Central Black Hole data
BHM = 4.e6 # [Solar masses] # Black Hole Mass
center = array([0.5, 0.5]) #[kpc] # Location of the SBH
BHmomentum = array([0.,0.]) # Momentum of the SBH
Mtotal = 5.e4 #[Solar masses] # Total mass of the galaxy
R = 30. #[kpc] # Galactic Radius
Rd = 0.1*R #[kpc] # Length Scale


######### MAIN PROGRAM ###################################

start = time.time()
# For Model 1 inizialization
#bodies = keplerDisk(N, max_mass, BHM, center, BHmomentum, R)

# For Model 2 inizialization
#bodies = exponentialDisk(Mtotal, max_mass, center, R, Rd)

# For Model 3 inizialization
bodies = galaxy(BHM, center, BHmomentum, Mtotal, max_mass, R, Rd)


velocity_plot(bodies, center, saveimage=True, plt_name='initial.png')
evolve(bodies, n, center, R, img_step)
end = time.time()

total_time = end - start
print(f'\nEl tiempo de cómputo para {len(bodies):.0f} cuerpos fue de {total_time:.2f} \n')
velocity_plot(bodies, center, saveimage=True, plt_name='final.png')
create_video(video_name = video_name)

