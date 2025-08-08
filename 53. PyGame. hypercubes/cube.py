'''
Rotating 3D Cube
'''

import arcade
from numpy import sin, cos, array, dot, copy

####################
SCREEN_WIDTH = 800
SCREEN_HEIGHT = 600
SCREEN_TITTLE = "Rotating 3D Cube"

PROJECTION = 'perspective'
cube_scale = 100
####################

def Rx(theta):
	return array([[1,0,0],
			        [0,cos(theta), sin(theta)],
			        [0,-sin(theta), cos(theta)]])

def Ry(theta):
	return array([[cos(theta), 0, -sin(theta)],
			         [0, 1, 0],
			         [sin(theta), 0, cos(theta)]])

def Rz(theta):
	return array([[cos(theta), -sin(theta), 0],
			         [sin(theta), cos(theta), 0],
			         [0, 0, 1]])

def persp_proj(z):
	'''
	Perspective projection
	'''
	global cube_scale 
	distance = 4
	cube_scale = 400
	return array([[1/(distance-z),0,0],
			      [0,1/(distance-z),0],
			      [0,0,0]])

def orth_proj():
	'''
	Orthograph projection
	'''
	return array([[1,0,0],
			      [0,1,0],
			      [0,0,0]])

def draw_edge(p1,p2):
	'''
	Draws a line between points p1 and p2
	'''
	arcade.draw_line(p1[0], p1[1], p2[0], p2[1], arcade.color.WHITE, 3)




class cube(arcade.Window):
	'''
	cube class
	'''
	def __init__(self):
		'''
		initializer
		'''
		super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITTLE)

		arcade.set_background_color(arcade.color.BLACK)

	def setup(self):
		'''
		Setup the all  the variables 
		'''
		# Define the vertices of the cube
		self.points = array([array([ 1., 1., 1.]),
		               array([ 1.,-1., 1.]),
		               array([-1.,-1., 1.]),
		               array([-1., 1., 1.]),
		               array([ 1., 1.,-1.]),
		               array([ 1.,-1.,-1.]),
		               array([-1.,-1.,-1.]),
		               array([-1., 1.,-1.])])

		self.rot_points = copy(self.points)
		
		self.edges = [[0,1], [1,2], [2,3], [3,0]]
		
		self.theta = 0.

	def on_draw(self):
		'''
		Draw the cube (vertices and edges)
		'''		
		arcade.start_render()
		self.clear()

		for k in range(4):
			i,j = self.edges[k]
			draw_edge(self.rot_points[i], self.rot_points[j])
			draw_edge(self.rot_points[i+4], self.rot_points[j+4])
			draw_edge(self.rot_points[i], self.rot_points[i+4])
			
		for point in self.rot_points:
			x,y = point[0], point[1]
			arcade.draw_circle_filled(x, y, 5, arcade.color.RED)

	def on_update(self, delta_time):
		'''
		Update the rotation of the cube
		'''

		for k in range(len(self.points)):
			point = self.points[k]
			#point = dot(Rz(self.theta),point) # z-rotation
			#point = dot(Rx(self.theta),point) # x-rotation
			point = dot(Ry(self.theta),point) # y-rotation
			if PROJECTION == 'perspective':
				proj_point = dot(persp_proj(point[2]),point) # perspective projection 
			else:
				proj_point = dot(orth_proj(),point) # orthograph projection  
			x = proj_point[0]*cube_scale + SCREEN_WIDTH//2
			y = proj_point[1]*cube_scale + SCREEN_HEIGHT//2
			z = proj_point[2]*cube_scale
			self.rot_points[k] = array([x,y,z])

		self.theta += 0.01



def main():
	'''
	Main function
	'''
	window = cube()
	window.setup()
	arcade.run()

if __name__ == "__main__":
	main()

