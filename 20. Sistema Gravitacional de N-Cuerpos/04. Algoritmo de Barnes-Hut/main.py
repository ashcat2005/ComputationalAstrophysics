from common import *

######### MAIN PROGRAM ########################################################

# Number of bodies (the actual number is smaller, because all bodies
# outside the initial radius are removed).
N = 200

# Mass of the N bodies.
max_mass = 50. # Solar masses

# Supermassive Central Black Hole data
BHM = 4.e6 # Solar masses
center = array([0.5, 0.5, 0.5]) # Location of the SBH
BHmomentum = array([0.,0.,0.]) # Momentum of the SBH

# Initial radius of the distribution
ini_radius = 10. #kpc

# Number of time-iterations executed by the program.
n = 50000

# Frequency at which .PNG images are written.
img_step = 250

# Folder to save the images
image_folder = 'images/'

# Name of the generated video
video_name = '200bodies50000.mp4'



bodies = system_init(N, max_mass, BHM, center, BHmomentum, ini_radius)
print('Total number of bodies: ', len(bodies))
evolve(bodies, n, center, ini_radius, img_step, image_folder, video_name)
create_video(image_folder, video_name)