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
dt = 1.e-2 # Gyr

# Theta-criterion of the Barnes-Hut algorithm.
theta = 0.3

#Radius on the image
scale_factor=.05 

#####################################################################################

class Node:
    '''---------------------------------------------------------------------
    A node object will represent a body (if node.child is None)
    or an abstract node of the octant-tree if it has node.child attributes.
    ----------------------------------------------------------------------'''
    def __init__(self, m, position, momentum):
        '''-------------------------------------------------------
        Creates a child-less node using the arguments
        ----------------------------------------------------------
        .mass     : scalar
        .position : NumPy array  with the coordinates [x,y,z]
        .momentum : NumPy array  with the components [px,py,pz]
        -------------------------------------------------------'''
        self.m = m
        self.m_pos = m * position #Mass times position. Usefull in calculating center of mass
        self.momentum = momentum
        self.child = None
        self.force = None         #Force on the node    
    
    def position(self):
        '''------------------------------------------
        Returns the physical coordinates of the node.
        -------------------------------------------'''
        return self.m_pos / self.m
        
    def reset_location(self):
        '''-----------------------------------------------------
        Resets the position of the node to the 0th-order octant.
        The size of the octant is reset to the value 1.0
        ------------------------------------------------------'''
        self.size = 1.0
        # The relative position inside the 0th-order octant is equal
        # to the current physical position
        self.relative_position = self.position().copy()
        
    def place_into_octant(self):
        '''-------------------------------------------------------------
        Places the node into next order octant.
        Returns the octant number according to the labels defined in the
        documentation.
        --------------------------------------------------------------'''
        # The next order octant will have half the size of the current octant
        self.size = 0.5 * self.size
        return self.subdivide(2) + 2*self.subdivide(1) + 4*self.subdivide(0)

    def subdivide(self, i):
        '''-------------------------------------------------------------------
        Places the node node into the next order octant along the direction i
        and recalculates the relative_position of the node inside this octant.
        --------------------------------------------------------------------'''
        self.relative_position[i] *= 2.0
        if self.relative_position[i] < 1.0:
            octant = 0
        else:
            octant = 1
            self.relative_position[i] -= 1.0
        return octant

def add(body, node):
    '''------------------------------------------------------------
    Defines the octo-tree by introducing a body and locating it
    according to three conditions (see documentation for details).
    Returns the updated node containing the body.
    ------------------------------------------------------------'''
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
    '''--------------------------------------------------------
    Returns the distance between node1 and node2. (Scaled down)
    ---------------------------------------------------------'''
    return norm(node1.position() - node2.position())

def gravitational_force(node1, node2):
    '''--------------------------------------------------------------
    Returns the gravitational force that node1 exerts on node2.
    A short distance cutoff is introduced in order to avoid numerical
    divergences in the gravitational force.
    ---------------------------------------------------------------'''
    cutoff_dist = 2.e-4
    d = distance_between(node1, node2)
    if d < cutoff_dist:
        return array([0., 0., 0.])
    else:
        return G*node1.m*node2.m*(node1.position() - node2.position())/d**3*scale_factor**2
    
def force_on(body, node, theta):
    '''-----------------------------------------------------------------------
    # Barnes-Hut algorithm: usage of the quad-tree. This function computes
    # the net force on a body exerted by all bodies in node "node".
    # Note how the code is shorter and more expressive than the human-language
    # description of the algorithm.
    ------------------------------------------------------------------------'''
    # 1. If the current node is an external node,
    #    calculate the force exerted by the current node on b.
    if node.child is None:
        return gravitational_force(node,body)#node.force_on(body)

    # 2. Otherwise, calculate the ratio s/d. If s/d < θ, treat this internal
    #    node as a single body, and calculate the force it exerts on body b.
    if node.size < distance_between(node,body)*theta:#node.distance(body) * theta:
        return gravitational_force(node,body)#node.force_on(body)
        
    # 3. Otherwise, run the procedure recursively on each child.
    return sum(force_on(body, c, theta) for c in node.child if c is not None)

def step(bodies, root, theta, dt):
    '''-----------------------------
    Euler method for time evolution.
    ------------------------------'''
    for body in bodies:
        body.force = force_on(body, root, theta)
        body.momentum += body.force*dt
        body.m_pos += scale_factor*body.momentum*dt

def PEFRL(bodies, root, theta, dt):
    '''-----------------------------
    PEFRL method for time evolution.
    ------------------------------'''
    epsilon=0.1786178958448091
    lannnda=-0.2123418310626054
    xsi=-0.06626458266981849
    for body in bodies:
        body.m_pos += epsilon*body.momentum*dt*scale_factor
        body.momentum += .5*(1-2*lannnda)*dt*force_on(body, root, theta)
        body.m_pos += xsi*body.momentum*dt*scale_factor
        body.momentum += lannnda*dt*force_on(body, root, theta)
        body.m_pos += (1-2*(xsi*epsilon))*body.momentum*dt*scale_factor
        body.momentum += lannnda*dt*force_on(body, root, theta)
        body.m_pos += xsi*body.momentum*dt*scale_factor
        body.momentum += .5*(1-2*lannnda)*dt*force_on(body, root, theta)
        body.m_pos += epsilon*body.momentum*dt*scale_factor
        
def func(x,Distribution,Point): 
    """------------------------------------------------------------------------
    Equation that follows the point of the wanted distribution that matches the 
    random one of a uniform distribution
    ---------------------------------------------------------------------------
       x            : Random variable in the wanted distribution (unkonwn)
       Distribution : Wanted distribution
       Point        : Random variable in the uniform distribution
    ------------------------------------------------------------------------"""
    return integrate.quad(Distribution,0,x)[0]-Point

def spiral_galaxy(N, max_mass, BHM, center, ini_radius, alpha, beta):
    '''-----------------------------------------------------------------------
    Use a radial distrubution of masses proportional to the brightness surface
    distributation to create a plain Bulb and Disk resembling an spiral galaxy
    --------------------------------------------------------------------------
       N            : Number of particles
       max_mass     : Biggest mass of the stars in the system 
       BHM          : Black Hole's mass
       center       : Black Hole position
       ini_radius   : Galaxy radius
       alpha        : Angle in the x,y plane
       beta         : Inclination
    ------------------------------------------------------------------------'''
    N -= 1
    random.seed(10)
    # Generates N random particles 
    positions = empty([N,3])
    momenta = empty([N,3])
    # Random masses varies between 1 solar mass and max_mass solar masses
    masses = random.random(N)*(max_mass-1.) + 1.
    #Parameters of the model of density of starts
    initial_density=.1
    const_bulb=.3
    const_disc=.8
    bulb_radius=0.2
    #Model of density normalized
    f1 = lambda x: initial_density*exp(-x**(1/4)/const_bulb)        #Bulge
    f2 = lambda x: f1(bulb_radius)*exp(-(x-bulb_radius)/const_disc) #Disc
    f = lambda x:  f1(x) if x<bulb_radius else f2(x)                #Piecewise 
    norm = integrate.quad(f,0,1)[0]                                  
    uf=lambda x: f(x)/norm                                          #Density function with integral=1
    #Random angle generation
    gamma = random.random(N)*2*pi
    #Random width
    width = .05*ini_radius                                          #Half of with in relation to the radius of the galaxy
    gross  = random.random(N)*2*width-width
    temp = beta
    #Uniform distribution to get random points
    Uniform = random.random(N)
    #Empty array for the points mapped from the uniform distribution
    Map=zeros(N)   
    for i in range(N):
        #Calls the function that maps the ramdom points to the wanted distribution for the radius 
        Map[i]=fsolve(func,0,args=(uf,Uniform[i]))*ini_radius
        #Creates an elipsoid in the region of the bulge
        if Map[i] < bulb_radius*ini_radius:
            a = 0.18*ini_radius
            bulg_countour = a*sqrt(1-(Map[i]/(bulb_radius*ini_radius))**2)
            gross[i] = random.random(1)*2*bulg_countour-bulg_countour
        #Adjustment for width
        beta += arctan(gross[i]/Map[i])
        Map[i] = sqrt(Map[i]**2+gross[i]**2)
        #Change to cartesian coordinates
        positions[i][0] = scale_factor*Map[i]*(cos(gamma[i])*cos(alpha)+
                                   sin(gamma[i])*cos(beta)*sin(alpha)) + center[0]
        positions[i][1] = scale_factor*Map[i]*(sin(gamma[i])*cos(beta)*cos(alpha)-
                                   cos(gamma[i])*sin(alpha))+ center[1]
        positions[i][2] = scale_factor*Map[i]*sin(gamma[i])*sin(beta) + center[2]
        # Keplerina velocity in the plain of the disc 
        Kep_v = sqrt(G*BHM/Map[i])
        vec_vel=array([-Map[i]*(sin(gamma[i])*cos(alpha)-cos(gamma[i])*cos(beta)*sin(alpha)),
                       Map[i]*(cos(gamma[i])*cos(beta)*cos(alpha)+sin(gamma[i])*sin(alpha)), 
                       Map[i]*cos(gamma[i])*sin(beta)])/Map[i]
        momenta[i][0] = masses[i]*Kep_v*vec_vel[0]
        momenta[i][1] = masses[i]*Kep_v*vec_vel[1]
        momenta[i][2] = masses[i]*Kep_v*vec_vel[2]
        beta = temp
    return masses, positions, momenta

def system_init(N, max_mass, BHM, center, BHmomentum, ini_radius, alpha, beta):
    '''-------------------------------------------------------------------
    Initializes the N-body system by defining the position and momentum
    of the bodies and creating the corresponding objects of the Node class
    --------------------------------------------------------------------'''
    #Defines initial conditions
    bodies = []
    bodies.append(Node(BHM, position=center, momentum=BHmomentum))   
    masses, positions, momenta = spiral_galaxy(N, max_mass, BHM, center, ini_radius, alpha, beta)
    #Creates nodes
    for i in range(N-1):
       bodies.append(Node(masses[i], positions[i], momenta[i]))
    return bodies

def evolve(bodies, n, center, ini_radius, img_step, image_folder='09. Barnes-Hut/BH 3D/main/Images/', video_name='my_video.mp4'):
    '''---------------------------------------------------------------------------------------------
    This function evolves the system in time using the Euler algorithm and the Barnes-Hut octo-tree
    ----------------------------------------------------------------------------------------------'''
    # Principal loop over time iterations.
    for i in range(n+1):
        # The octo-tree is recomputed at each iteration.
        root = None
        for body in bodies:
            body.reset_location()
            root = add(body, root) 
        # Evolution using the integration method
        step(bodies, root, theta, dt)
        #PEFRL(bodies, root, theta, dt)
        # Write the image files
        if i%img_step==0:
            print("Writing image at time {0}".format(i))
            plot_bodies(bodies, i//img_step,image_folder)

def plot_bodies(bodies, i, image_folder='09. Barnes-Hut/BH 3D/main/Images/'):
    '''---------------------------------------------------------
    Writes an image file with the current position of the bodies
    ---------------------------------------------------------'''
    plt.rcParams['grid.color'] = 'dimgray'
    fig = plt.figure(figsize=(10,10), facecolor='black')
    ax = plt.gcf().add_subplot(111, projection='3d')
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_zlim([0,1])
    ax.set_facecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('dimgray')
    ax.yaxis.pane.set_edgecolor('dimgray')
    ax.zaxis.pane.set_edgecolor('dimgray')
    #ax.view_init(90, -90)
    for body in bodies:
        pos = body.position()
        if body.m>100.:
            ax.scatter(pos[0], pos[1], pos[2], marker='.', color='lightcyan')
        else:
            ax.scatter(pos[0], pos[1], pos[2], marker='.', color='darkorchid')
    print(" ")
    plt.gcf().savefig(image_folder+'bodies3D_{0:06}.png'.format(i))
    plt.close()

def create_video(image_folder='09. Barnes-Hut/BH 3D/main/Images/', video_name='09. Barnes-Hut/BH 3D/main/vid/my_videoN.mp4'):
    '''-----------------------------------------------
    Creates a .mp4 video using the stored files images
    -----------------------------------------------'''
    from os import listdir
    import moviepy.video.io.ImageSequenceClip
    fps = 15
    image_files = [image_folder+img for img in sorted(listdir(image_folder)) if img.endswith(".png")]
    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
    clip.write_videofile(video_name) 
