# 2-D Barnes-Hut Algorithm

from copy import deepcopy
from numpy import array, ones, empty, random, sqrt, pi, exp, linspace
from numpy import sin, cos
import matplotlib.pyplot as plt


##### Simulation Parameters #########################################################

# Gravitational constant in units of kpc^3 M_sun^-1 Gyr-2
G = 4.4985022e-6

# Discrete time step.
dt = 5.e-3 # Gyr

# Theta-criterion of the Barnes-Hut algorithm.
theta = 0.3

#####################################################################################

class Node:
    '''
    A node object will represent a body (if node.child is None)
    or an abstract node of the quad-tree if it has node.child attributes.
    '''
    def __init__(self, m, position, momentum, R):
        '''
        Creates a child-less node using the arguments
        .mass : scalar mass of the node
        .position : NumPy array  with the coordinates [x,y]
        .momentum : NumPy array  with the components [px,py]
        .R = initial radius of the distribution
        '''
        self.m = m
        self.m_pos = m * position
        self.momentum = momentum
        self.R = R
        self.child = None
    
    def position(self):
        '''
        Returns the physical coordinates of the node.
        '''
        return self.m_pos / self.m
        
    def reset_location(self):
        '''
        Resets the position of the node to the 0th-order quadrant.
        The size of the quadrant is reset to the value 2*R.
        '''
        self.size = 2*self.R
        # The relative position inside the 0th-order quadrat is equal
        # to the current physical position.
        self.relative_position = self.position().copy()
        
    def place_into_quadrant(self):
        '''
        Places the node into next order quadrant.
        Returns the quadrant number according to the labels defined in the
        documentation.
        '''
        # The next order quadrant will have half the side of the current quadrant
        self.size = 0.5 * self.size
        return self.subdivide(1) + 2*self.subdivide(0)

    def subdivide(self, i):
        '''
        Places the node node into the next order quadrant along the direction i
        and re-calculates the relative_position of the node inside this quadrant.
        '''
        self.relative_position[i] *= 2.0
        if self.relative_position[i] < 1.0:
            quadrant = 0
        else:
            quadrant = 1
            self.relative_position[i] -= 1.0
        return quadrant



def add(body, node):
    '''
    Defines the quad-tree by introducing a body and locating it
    according to three conditions (see documentation for details).
    Returns the updated node containing the body.
    '''
    smallest_quadrant = 1.e-4 # Lower limit for the side-size of the quadrants

    # Case 1. If node does not contain a body, the body is put in here
    new_node = body if node is None else None
    
    if node is not None and node.size > smallest_quadrant:
        # Case 3. If node is an external node, then the new body can not
        # be put in there. We have to verify if it has .child attribute
        if node.child is None:
            new_node = deepcopy(node)
            # Subdivide the node creating 4 children
            new_node.child = [None for i in range(4)]
            # Place the body in the appropiate quadrant
            quadrant = node.place_into_quadrant()
            new_node.child[quadrant] = node
        # Case 2. If node is an internal node, it already has .child attribute
        else:
            new_node = node

        # For cases 2 and 3, it is needed to update the mass and the position
        # of the node
        new_node.m += body.m
        new_node.m_pos += body.m_pos
        # Add the new body into the appropriate quadrant.
        quadrant = body.place_into_quadrant()
        new_node.child[quadrant] = add(body, new_node.child[quadrant])
    return new_node


def distance_between(node1, node2):
    '''
    Returns the distance between node1 and node2.
    '''
    d12 = node1.position() - node2.position()
    return sqrt(d12.dot(d12))


def gravitational_force(node1, node2):
    '''
    Returns the gravitational force that node1 exerts on node2.
    A short distance cutoff is introduced in order to avoid numerical
    divergences in the gravitational force.
    '''
    cutoff_dist = 1.e-3 #6  # [kpc]
    d12 = node1.position() - node2.position()
    d = sqrt(d12.dot(d12))
    if d < cutoff_dist:
        # Returns no Force to prevent divergences!
        return array([0., 0.])
    else:
        # Gravitational force
        return G*node1.m*node2.m*(d12)/d**3


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
        return gravitational_force(node,body)

    # 2. Otherwise, calculate the ratio s/d. If s/d < θ, treat this internal
    #    node as a single body, and calculate the force it exerts on body b.
    if node.size < distance_between(node,body)*theta:
        return gravitational_force(node,body)

    # 3. Otherwise, run the procedure recursively on each child.
    return sum(force_on(body, c, theta) for c in node.child if c is not None)


def verlet(bodies, root, theta, dt):
    '''
    Velocity-Verlet method for time evolution.
    '''
    for body in bodies:
        force = force_on(body, root, theta)
        body.momentum += 0.5*force*dt
        body.m_pos += body.momentum*dt
        body.momentum += 0.5*force_on(body, root, theta)*dt

        
def keplerDisk(N, max_mass, BHM, center, BHmomentum, R):
    '''
    This function initialize the N-body system by randomly defining
    the position vectors fo the bodies and creating the corresponding
    objects of the Node class
    '''
    bodies = []
    bodies.append(Node(BHM, position=center, momentum=BHmomentum, R=R))

    K = 2*N
    random.seed(413)
    # Random masses between 1 solar mass and max_mass solar masses
    mass = random.random(K)*(max_mass-1.) + 1.
    # x- and y- positions are initialized inside a square with
    # a side of length = 2*R.
    posx = random.random(K)*2*R + center[0] - R
    posy = random.random(K)*2*R + center[1] - R
    i=0
    j=0
    #Loop until complete the random N-1 bodies or use the K generated bodies
    while i<K and j<N-1:
        position = array([posx[i],posy[i]])
        r = position - center
        norm_r = sqrt(r.dot(r))
        if norm_r < R:
            # We use the projection of the Keplerian velocity to define the momentum
            Kep_v = sqrt(G*BHM/norm_r) # Keplerian velocity
            v_mod = 1. # Modification of the Keplerian velocity
            momentum = mass[i]*Kep_v*v_mod*array([-r[1], r[0]])/norm_r
            bodies.append(Node(mass[i], position, momentum, R))
            j+=1
        i+=1

    print(f'\nNúmero Total de Cuerpos creados: {len(bodies):.0f}\n')
    return bodies


def exponentialDisk(Mtotal, max_mass, center, R, Rd):
    '''
    This function initialize the N-body system using a
    exponential surface density model. 
    The mass and  position of the particles are
    randomly created 
    '''
    from scipy.special import iv,kv
    bodies = []
    
    M = Mtotal
    Sigma0 = M/(2*pi*Rd**2*(1-exp(-R/Rd)*(1+R/Rd)))
    R_range = linspace(0,R,30)
    for i in range(1,len(R_range)):
        m_in = 2*pi*Sigma0*Rd**2*(1-exp(-R_range[i-1]/Rd)*(1+R_range[i-1]/Rd))
        m_out = 2*pi*Sigma0*Rd**2*(1-exp(-R_range[i]/Rd)*(1+R_range[i]/Rd))
        m = m_out - m_in
        while m>0:
            mass = random.random()*(max_mass-1.) + 1.
            r = random.random()*(R_range[i]-R_range[i-1]) + R_range[i-1]
            th = random.random()*2*pi
            posx = r*cos(th) + center[0]
            posy = r*sin(th) + center[1]
            position = array([posx,posy])
            y = r/(2*Rd)
            v = sqrt(4*pi*Sigma0*G*Rd*(y**2)*(iv(0,y)*kv(0,y) - iv(1,y)*kv(1,y)))
            momentum = mass*v*array([-sin(th), cos(th)])
            bodies.append(Node(mass, position, momentum, R))
            m -= mass

    print(f'\nNúmero Total de Cuerpos creados: {len(bodies):.0f}\n')
    return bodies

def galaxy(BHM, center, BHmomentum, Mtotal, max_mass, R, Rd):
    '''
    This function initialize the N-body system using a
    galaxy model with a surface density model. 
    The mass and position of the particles are
    randomly created 
    '''
    from scipy.special import iv,kv
    bodies = []
    bodies.append(Node(BHM, position=center, momentum=BHmomentum, R=R))
    
    M = Mtotal
    Sigma0 = M/(2*pi*Rd**2*(1-exp(-R/Rd)*(1+R/Rd)))
    R_range = linspace(0,R,30)
    for i in range(1,len(R_range)):
        m_in = 2*pi*Sigma0*Rd**2*(1-exp(-R_range[i-1]/Rd)*(1+R_range[i-1]/Rd))
        m_out = 2*pi*Sigma0*Rd**2*(1-exp(-R_range[i]/Rd)*(1+R_range[i]/Rd))
        m = m_out - m_in
        while m > 0:
            mass = random.random()*(max_mass-1.) + 1.
            r = random.random()*(R_range[i]-R_range[i-1]) + R_range[i-1]
            th = random.random()*2*pi
            posx = r*cos(th) + center[0]
            posy = r*sin(th) + center[1]
            position = array([posx,posy])
            r = position - center
            norm_r = sqrt(r.dot(r))
            y = r/(2*Rd)
            # We use the projection of the Keplerian velocity to define the momentum
            Kep_v = sqrt(G*BHM/norm_r) # Keplerian velocity
            v_mod = 1. # Modification of the Keplerian velocity
            momentum = mass*Kep_v*v_mod*array([-r[1], r[0]])/norm_r
            #v = sqrt(4*pi*Sigma0*G*Rd*(y**2)*(iv(0,y)*kv(0,y) - iv(1,y)*kv(1,y)))
            #momentum = mass*v*array([-sin(th), cos(th)])
            bodies.append(Node(mass, position, momentum, R))
            m -= mass

    print(f'\nNúmero Total de Cuerpos creados: {len(bodies):.0f}\n')
    return bodies

def evolve(bodies, n, center, R, img_step):
    '''
    This function evolves the system in time using the Verlet algorithm 
    and the Barnes-Hut quad-tree
    '''
    # Limits for the axes in the plot
    axis_limit = 1.1*R
    lim_inf = [center[0]-axis_limit, center[1]-axis_limit]
    lim_sup = [center[0]+axis_limit, center[1]+axis_limit]
    # Principal loop over time iterations.
    for i in range(n+1):
        # The quad-tree is recomputed at each iteration.
        root = None
        for body in bodies:
            body.reset_location()
            root = add(body, root)
        # Evolution using the Verlet method
        verlet(bodies, root, theta, dt)
        # Write the image files
        if i%img_step==0:
            print("Creando imagen en el paso {0}".format(i))
            plot_bodies(bodies, i//img_step, lim_inf, lim_sup)


def plot_bodies(bodies, i, lim_inf, lim_sup, image_folder='09. Barnes-Hut/BH 2D/imagesBH/'):
    '''
    Writes an image file with the current position of the bodies
    '''
    fig,ax = plt.subplots(figsize=(10,10), facecolor='black')
    ax.set_xlim([lim_inf[0], lim_sup[0]])
    ax.set_ylim([lim_inf[1], lim_sup[1]])
    ax.set_facecolor('black')
    ax.scatter(bodies[0].position()[0], bodies[0].position()[1], marker='.', color='lightcyan')
    for body in bodies[1:]:
        pos = body.position()
        ax.scatter(pos[0], pos[1], marker='.', color='lightcyan', s=body.m*0.06+0.94)
    plt.savefig(image_folder+'bodies_{0:06}.png'.format(i))
    plt.close()




def create_video(im_folder='09. Barnes-Hut/BH 2D/imagesBH/', v_folder='09. Barnes-Hut/BH 2D/videos/', video_name='my_video.mp4'):
    '''
    Creates a .mp4 video using the stored files images
    '''
    from os import listdir
    import moviepy.video.io.ImageSequenceClip
    fps = 15
    image_files = [im_folder+img for img in sorted(listdir(im_folder)) if img.endswith(".png")]
    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
    clip.write_videofile(v_folder+video_name)



def create_avi_video(im_folder='09. Barnes-Hut/BH 2D/imagesBH/', video_name = 'video.avi'):
    '''
    Creates a .avi video using the stored files images
    '''
    import cv2
    from os import listdir
    from os.path import join
    images = [img for img in listdir(im_folder) if img.endswith(".png")]
    frame = cv2.imread(join(im_folder, images[0]))
    height, width, layers = frame.shape
    video = cv2.VideoWriter(video_name, 0, 1, (width,height))
    for image in images:
        video.write(cv2.imread(join(im_folder, image)))
    cv2.destroyAllWindows()
    video.release()


def velocity_plot(bodies, center, saveimage=False, plt_folder='09. Barnes-Hut/BH 2D/plots/', plt_name='plot.png'):
    plt.figure(figsize=(10,7))
    plt.xlabel(r'$r$')
    plt.ylabel(r'$\left| v \right| $')
    for body in bodies:
        pos = body.position()
        vel = body.momentum/body.m
        r = sqrt((pos[0]-center[0])**2 + (pos[1]-center[1])**2)
        v = sqrt(vel[0]**2 + vel[1]**2)
        plt.scatter(r, v, color='black', s=1)
    if saveimage:
        plt.savefig(plt_folder+plt_name)
    else:
        plt.show()




if __name__=="__main__":
    '''
    Example of a randomly generated N-body system to be evolved using
    the Barnes-Hut Algorithm and the Verlet method
    '''
    ######### GENERAL SETTINGS ###############################
    import time

    # Mass of the N bodies.
    max_mass = 50. # Solar masses

    # Number of time-iterations executed by the program.
    n = 2600

    # Frequency at which .PNG images are written.
    img_step = 50

    # Name of the generated video
    video_name = 'N-bodies.mp4'


    ######### MODEL OPTIONS ##################################

    # Model 1: SMBH + Keplerian velocity
    # Supermassive Central Black Hole data
    BHM = 4.e6 # [Solar masses] # Black Hole Mass
    center = array([0.5, 0.5]) #[kpc] # Location of the SBH
    BHmomentum = array([0.,0.]) # Momentum of the SBH

    N = 2000 # Number of bodies
    R = 30. #[kpc]   # Initial radius of the distribution

    # Model 2: Exponential Surface Density
    Mtotal = 5.e4 #[Solar masses] # Total mass of the galaxy
    R = 30. #[kpc] # Galactic Radius
    Rd = 0.1*R #[kpc] # Length Scale


    ######### MAIN PROGRAM ###################################

    start = time.time()
    # For Model 1 inizialization
    #bodies = system_init(N, max_mass, BHM, center, BHmomentum, R)

    # For Model 2 inizialization
    bodies = exponentialDisk(Mtotal, max_mass, center, R, Rd)
    velocity_plot(bodies, center, saveimage=True, plt_name='initial.png')
    evolve(bodies, n, center, R, img_step)
    end = time.time()

    total_time = end - start
    print(f'\nEl tiempo de cómputo para {len(bodies):.0f} cuerpos fue de {total_time:.2f} \n')
    velocity_plot(bodies, center, saveimage=True, plt_name='final.png')
    create_video(video_name = video_name)


