from common import *

######### MAIN PROGRAM ########################################################

# Number of bodies (may be smaller according to the distribution chosen).
N = 100

# Mass of the N bodies.
max_mass = 50. # Solar masses

# Supermassive Central Black Hole data
BHM = 1.e6 # Solar masses
BHposition = array([.5, .4, .5]) # Location of the SBH
BHmomentum = array([0.,0.,0.]) # Momentum of the SBH

#Parameters of the galaxy plane orientation 
beta=0      #Inclination
alpha=0     #Angle in the plain x,y

# Initial radius of the distribution
ini_radius = 10 #kpc

# Number of time-iterations executed by the program.
n = 100000

# Frequency at which .PNG images are written.
img_step = 250

# Folder to save the images
image_folder = 'images/'

# Name of the generated video
video_name = '200bodies50000.mp4'

bodies = system_init(N, max_mass, BHM, BHposition, BHmomentum, ini_radius, alpha, beta)
print('Total number of bodies: ', len(bodies))
evolve(bodies, n, BHposition, ini_radius, img_step, image_folder, video_name)
create_video(image_folder, video_name)