# Barnes-Hut Algorithm

from copy import deepcopy
from numpy import array, ones, empty, random, sqrt, exp, pi, sin, cos, arctan, zeros
from numpy.linalg import norm
import scipy.integrate as integrate
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D


##### Simulation Parameters #########################################################

# Gravitational constant in units of kpc^3 M_sun^-1 Gyr-2
G = 4.4985022e-6

# Discrete time step.
dt = 1.e-3 # Gyr

# Theta-criterion of the Barnes-Hut algorithm.
theta = 0.3

#Radius on the image
scale_factor=.025 #---> 0.25/ini_radius

#####################################################################################

class Node:
    '''
    A node object will represent a body (if node.child is None)
    or an abstract node of the octant-tree if it has node.child attributes.
    '''
    def __init__(self, m, position, momentum):
        '''
        Creates a child-less node using the arguments
        .mass : scalar
        .position : NumPy array  with the coordinates [x,y,z]
        .momentum : NumPy array  with the components [px,py,pz]
        '''
        self.m = m
        self.m_pos = m * position
        self.momentum = momentum
        self.child = None
        self.force = None
    
    def position(self):
        '''
        Returns the physical coordinates of the node.
        '''
        return self.m_pos / self.m
        
    def reset_location(self):
        '''
        Resets the position of the node to the 0th-order octant.
        The size of the octant is reset to the value 1.0
        '''
        self.size = 1.0
        # The relative position inside the 0th-order octant is equal
        # to the current physical position.
        self.relative_position = self.position().copy()
        
    def place_into_octant(self):
        '''
        Places the node into next order octant.
        Returns the octant number according to the labels defined in the
        documentation.
        '''
        # The next order octant will have half the size of the current octant
        self.size = 0.5 * self.size
        return self.subdivide(2) + 2*self.subdivide(1) + 4*self.subdivide(0)

    def subdivide(self, i):
        '''
        Places the node node into the next order octant along the direction i
        and recalculates the relative_position of the node inside this octant.
        '''
        self.relative_position[i] *= 2.0
        if self.relative_position[i] < 1.0:
            octant = 0
        else:
            octant = 1
            self.relative_position[i] -= 1.0
        return octant

def add(body, node):
    '''
    Defines the octo-tree by introducing a body and locating it
    according to three conditions (see documentation for details).
    Returns the updated node containing the body.
    '''
    smallest_quadrant = 1.e-4 # Lower limit for the size of the octants
    # Case 1. If node does not contain a body, the body is put in here
    new_node = body if node is None else None
    
    if node is not None and node.size > smallest_quadrant:
        # Case 3. If node is an external node, then the new body can not
        # be put in there. We have to verify if it has .child attribute
        if node.child is None:
            new_node = deepcopy(node)
            # Subdivide the node creating 8 children
            new_node.child = [None for i in range(8)]
            # Place the body in the appropiate octant
            quadrant = node.place_into_octant()
            new_node.child[quadrant] = node
        # Case 2. If node is an internal node, it already has .child attribute
        else:
            new_node = node

        # For cases 2 and 3, it is needed to update the mass and the position
        # of the node
        new_node.m += body.m
        new_node.m_pos += body.m_pos
        # Add the new body into the appropriate octant.
        octant = body.place_into_octant()
        new_node.child[octant] = add(body, new_node.child[octant])
    return new_node

def distance_between(node1, node2):
    '''
    Returns the distance between node1 and node2.
    '''
    #print(node1.position())
    #print(node2.position())
    #print(norm(node1.position() - node2.position()))
    #print("")
    return norm(node1.position() - node2.position())/scale_factor

def gravitational_force(node1, node2):
    '''
    Returns the gravitational force that node1 exerts on node2.
    A short distance cutoff is introduced in order to avoid numerical
    divergences in the gravitational force.
    '''
    cutoff_dist = 2.e-4
    d = distance_between(node1, node2)
    if d < cutoff_dist:
        #print('Collision!')
        #Returns no Force!
        return array([0., 0., 0.])
        
    else:
        # Gravitational force
        return G*node1.m*node2.m*(node1.position() - node2.position())/d**3/scale_factor
    
def force_on(body, node, theta):
    '''
    # Barnes-Hut algorithm: usage of the quad-tree. This function computes
    # the net force on a body exerted by all bodies in node "node".
    # Note how the code is shorter and more expressive than the human-language
    # description of the algorithm.
    '''
    # 1. If the current node is an external node,
    #    calculate the force exerted by the current node on b.
    if node.child is None:
        return gravitational_force(node,body)#node.force_on(body)

    # 2. Otherwise, calculate the ratio s/d. If s/d < θ, treat this internal
    #    node as a single body, and calculate the force it exerts on body b.
    if node.size < scale_factor*distance_between(node,body)*theta:#node.distance(body) * theta:
        return gravitational_force(node,body)#node.force_on(body)
        
    # 3. Otherwise, run the procedure recursively on each child.
    return sum(force_on(body, c, theta) for c in node.child if c is not None)

def verlet(bodies, root, theta, dt):
    '''
    Verlet method for time evolution.
    '''
    for body in bodies:
        body.force = force_on(body, root, theta)
        body.momentum += body.force*dt
        body.m_pos += scale_factor*body.momentum*dt
                  
def uniform_galaxy(N, max_mass, BHM, center, ini_radius):
    '''
    Randomly generate the system of N particles.
    Returns
    - Masses
    - Positions
    - Momentum (Based on a Keplerian velocity)
    '''
    # We will generate K=2*N random particles from which we will chose
    # only N-1 bodies for the system
    K = 2*N
    random.seed(2)
    masses = empty([N-1])
    positions = empty([N-1,3])
    momenta = empty([N-1,3])
    # Random masses between 1 solar mass and max_mass solar masses
    mass = random.random(K)*(max_mass-1.) + 1.
    # x-, y- and z- positions are initialized inside a square with
    # a side of length = 2*ini_radius.
    posx = random.random(K) *2.*ini_radius-ini_radius
    posy = random.random(K) *2.*ini_radius-ini_radius
    posz = zeros(K)
    i=0
    j=0
    #Loop until complete the random N-1 bodies or use the K generated bodies
    while i<K and j<N-1:
        position = array([posx[i],posy[i],posz[i]])
        norm_r = norm(position)
        if norm_r < ini_radius:
            masses[j] = mass[i]
            positions[j] = scale_factor*position+center
            # We use the projection of the Keplerina velocity to define the momentum
            Kep_v = sqrt(G*BHM/norm_r) # Keplerian velocity
            momenta[j] = mass[i]*Kep_v*array([-position[1], position[0], 0.])/norm_r
            j+=1
        i+=1
    return masses, positions, momenta

def func(x,Distribution,Point): 
    """
    Equation that follows the point of the wanted distribution that matches the 
    random one of a uniform distribution
    """
    return integrate.quad(Distribution,0,x)[0]-Point

def spiral_galaxy(N, max_mass, BHM, center, ini_radius, alpha, beta):
    '''
    Version 1.1:
    Use a radial distrubution of masses proportional to the brightness surface
    distributation to create a plain Bulb and Disk resembling an spiral galaxy
    Returns
    - Masses
    - Positions
    - Momentum (Based on a Keplerian velocity)
    '''
    # We will generate N random particles 
    positions = empty([N,3])
    momenta = empty([N,3])
    # Random masses between 1 solar mass and max_mass solar masses
    masses = random.random(N)*(max_mass-1.) + 1.
    # r- initializec acording to a distribution, gamma- initialized ramdomly
    #Parameters of the model of density
    initial_density=.1
    const_bulb=.3
    const_disc=.8
    bulb_radius=0.2     #according to the initial radius 
    bulb_hight=3/20     #according to the initial radius 
    #Model of density normalized
    f1 = lambda x: initial_density*exp(-x**(1/4)/const_bulb)        #Bulb
    f2 = lambda x: f1(bulb_radius)*exp(-(x-bulb_radius)/const_disc) #Disc
    f = lambda x:  f1(x) if x<bulb_radius else f2(x)                #Piecewise 
    norm = integrate.quad(f,0,1)[0]                                  
    uf=lambda x: f(x)/norm                                          #Density function with integral=1
    #Uniform distribution to get random points as normal
    random.seed(413)
    Uniform=random.random(N)
    #Random angle generation
    gamma=random.random(N)*2*pi
    #Empty array for the points mapped from the uniform distribution
    Map=empty(N)
    #Noise in the perpendicular direction in order to create the bulb and the width of the disc
    random.seed(1)
    w_noise=(random.random(N)*2-1) 
    width_disc=1/20*ini_radius
    Rad=empty(N)
    phi=empty(N)
    a=(bulb_radius**2)/(1-(bulb_hight/width_disc)**2)*ini_radius**2
    for i in range(N):
        #Calls the function that maps the ramdom points to the wanted distribution for the radius 
        Map[i]=fsolve(func,.5,args=(uf,Uniform[i]))*ini_radius
        if Map[i]<0:
            Map[i]=ini_radius*Uniform[i]
        #Correction of the radius and angle according to the width of the disc
        if Map[i]>0.2*ini_radius:
            Width=w_noise[i]*width_disc
        else:
            Width=w_noise[i]*bulb_hight*sqrt(1-(Rad[i]**2)/a)
        Rad[i]=sqrt(Map[i]**2+Width**2)
        phi[i]=arctan(Width/Map[i])
        #Change to cartesian coordinates
        positions[i][0] = Rad[i]*(cos(gamma[i])*cos(alpha)+
                                  sin(gamma[i])*cos(beta+phi[i])*sin(alpha)) + center[0]
        positions[i][1] = Rad[i]*(sin(gamma[i])*cos(beta+phi[i])*cos(alpha)-
                                  cos(gamma[i])*sin(alpha))+ center[1]
        positions[i][2] = Rad[i]*sin(gamma[i])*sin(beta+phi[i]) + center[2]
        # We use the Keplerina velocity to define the momentum in the plain of the disc 
        Kep_v = sqrt(G*BHM/Map[i]) # Keplerian velocity
        vec_mom=array([-Map[i]*(sin(gamma[i])*cos(alpha)-cos(gamma[i])*cos(beta)*sin(alpha)),
                       Map[i]*(cos(gamma[i])*cos(beta)*cos(alpha)+sin(gamma[i])*sin(alpha)), 
                       Map[i]*cos(gamma[i])*sin(beta)])/Map[i]  #unit vector for direction
        momenta[i] = masses[i]*Kep_v*vec_mom
    return masses, positions, momenta

def second_galaxy(bodies,ini_radius):
    """
    Creation of a secundary galaxy to simulate the merge of two galaxies
    """
    #BlackHole in the center of the secondary galaxy
    BHM_2 = 1.e5                            # Solar masses
    center_2 = array([.8,.8,0])         # Location of the SBH
    BHmomentum_2 = array([-BHM_2*.1,0,0.])  # Momentum of the SBH
    #Number of bodies
    N_2=30
    #Mass of the N bodies
    max_mass_2= 50.                     # Solar masses
    # Initial radius of the distribution
    ini_radius_2 = .5*ini_radius        #kpc
    masses_2, positions_2, momenta_2 = uniform_galaxy(N_2, max_mass_2, BHM_2, center_2, 
                                                     ini_radius_2)
    bodies.append(Node(BHM_2, position=center_2, momentum=BHmomentum_2))
  
    for i in range(len(masses_2)-1):
        momenta_2[i]=momenta_2[i]+BHmomentum_2/BHM_2*masses_2[i]
        bodies.append(Node(masses_2[i], positions_2[i], momenta_2[i]))
    return bodies

def system_init(N, max_mass, BHM, center, BHmomentum, ini_radius, alpha, beta):
    '''
    This function initialize the N-body system by randomly defining
    the position vectors fo the bodies and creating the corresponding
    objects of the Node class
    '''
    bodies = []
    bodies.append(Node(BHM, position=center, momentum=BHmomentum))
        
    masses, positions, momenta = uniform_galaxy(N, max_mass, BHM, center, ini_radius)
    #masses, positions, momenta = spiral_galaxy(N, max_mass, BHM, center, ini_radius, alpha, beta)
    for i in range(N-1):
       bodies.append(Node(masses[i], positions[i], momenta[i])) 
    bodies=second_galaxy(bodies, ini_radius)
    return bodies

def evolve(bodies, n, center, ini_radius, img_step, image_folder='images/', video_name='my_video.mp4'):
    '''
    This function evolves the system in time using the Verlet algorithm and the Barnes-Hut octo-tree
    '''
    # Principal loop over time iterations.
    for i in range(n+1):
        # The octo-tree is recomputed at each iteration.
        root = None
        for body in bodies:
            body.reset_location()
            root = add(body, root) 
        # Evolution using the integration method
        verlet(bodies, root, theta, dt)
        # Write the image files
        if i%img_step==0:
            print("Writing image at time {0}".format(i))
            plot_bodies(bodies, i//img_step,image_folder)

def plot_bodies(bodies, i, image_folder='images/'):
    '''
    Writes an image file with the current position of the bodies
    '''
    plt.rcParams['grid.color'] = 'dimgray'
    #plt.rcParams['axes.edgecolor'] = 'dimgray'
    #plt.rcParams['axes.labelcolor'] = 'dimgray'
    fig = plt.figure(figsize=(10,10), facecolor='black')
    ax = plt.gcf().add_subplot(111, projection='3d')
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_zlim([0,1])
    #ax.set_proj_type('ortho')
    ax.set_facecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('dimgray')
    ax.yaxis.pane.set_edgecolor('dimgray')
    ax.zaxis.pane.set_edgecolor('dimgray')
    ax.view_init(90, -90)
    #ax.xaxis.label.set_color('dimgray')
    #ax.yaxis.label.set_color('dimgray')
    #ax.zaxis.label.set_color('dimgray')
    #ax.set_xlabel('x')
    #ax.set_ylabel('y')
    #ax.set_zlabel('z')
    #ax.grid(False)
    for body in bodies:
        pos = body.position()
        if body.m>100.:
            ax.scatter(pos[0], pos[1], pos[2], marker='.', color='lightcyan')
        else:
            ax.scatter(pos[0], pos[1], pos[2], marker='.', color='red')
    print(" ")
    plt.gcf().savefig(image_folder+'bodies3D_{0:06}.png'.format(i))
    plt.close()

def create_video(image_folder='images/', video_name='my_video.mp4'):
    '''
    Creates a .mp4 video using the stored files images
    '''
    from os import listdir
    #from os import environ
    #environ["IMAGEIO_FFMPEG_EXE"] = "/Users/beyonder88/Downloads/ffmpeg"
    import moviepy.video.io.ImageSequenceClip
    fps = 15
    image_files = [image_folder+img for img in sorted(listdir(image_folder)) if img.endswith(".png")]
    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
    clip.write_videofile(video_name)

def create_avi_video(image_folder='images/', video_name = 'video.avi'):
    '''
    Creates a .avi video using the stored files images
    '''
    import cv2
    from os import listdir
    from os.path import join
    images = [img for img in listdir(image_folder) if img.endswith(".png")]
    frame = cv2.imread(join(image_folder, images[0]))
    height, width, layers = frame.shape
    video = cv2.VideoWriter(video_name, 0, 1, (width,height))
    for image in images:
        video.write(cv2.imread(join(image_folder, image)))
    cv2.destroyAllWindows()
    video.release()    