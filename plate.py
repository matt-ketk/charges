import numpy as np

from conductor import Conductor
from constants import Constants

class Plate(Conductor):
	vertex1 = tuple()  # 3D coordinate of one vertex
	vertex2 = tuple()  # 3D coordinate of the opposing vertex
	vertices = []  # list of all 8 vertices
	shape = tuple()  # tuple of positive floats indicating dimensions
	faces = []  # list of all 6 faces on the plate stored as tuples of the rectangle's corners
	thicknessDimension = 0  # 0, 1, or 2 indicating which dimension of the plate is its thickness (the index of the smallest value in shape)
	thickness = 0

	def __init__(self, vertex1, vertex2, resistivity=Constants.COPPER_RESISTIVITY):
		"""
		input the 3D coordinate numpy arrays of 2 opposing vertices of the plate
		requires that the plate's edges lie parallel to the coordinate system
		generates all the corners of the prism
		"""
		self.vertex1 = vertex1
		self.vertex2 = vertex2
		self.shape = np.absolute(vertex2 - vertex1)
		self.thicknessDimension = np.argmin(self.shape)
		self.thickness = self.shape[self.thicknessDimension]
		self.vertices = []
		for x in (vertex1[0], vertex2[0]):
			for y in (vertex1[1], vertex2[1]):
				for z in (vertex1[2], vertex2[2]):
					self.vertices.append(np.array((x, y, z)))
		self.faces = []
		self.faces.append(tuple(v for v in self.vertices[:4]))
		self.faces.append(tuple(v for v in self.vertices[4:]))
		self.faces.append(tuple(self.vertices[i] for i in range(len(self.vertices)) if i % 4 < 2))
		self.faces.append(tuple(self.vertices[i] for i in range(len(self.vertices)) if i % 4 >= 2))
		self.faces.append(tuple(v for v in self.vertices[::2]))
		self.faces.append(tuple(v for v in self.vertices[1::2]))

		super(Plate, self).__init__(resistivity=resistivity)

	@staticmethod
	def rectangleCollision(prevPos, pos, rCorners, velocity, dampingFactor=1):
		"""
		rectangle must be parallel to xy, xz, or yz for now
		:param prevPos: position of a particle at the previous iteration
		:param pos: position of that particle at the current iteration
		:param rCorners: a tuple of the corner points of a rectangle in 3D space to check for collisions
		:param velocity: the initial velocity vector of the particle
		:param dampingFactor: the factor by which the particle's velocity decreases upon impact
		:return: the collision position of the particle and its new velocity. if no collision returns None
		"""
		# find orientation of the rectangle
		for i in range(3):
			if all(rCorners[0][i] == rCorners[j][i] for j in range(4)):
				z = i
				break
		if z != 0:
			x = 0
			if z != 1:
				y = 1
			else:
				y = 2
		else:
			x = 1
			y = 2
		planeZ = rCorners[0][z]
		xBounds = (min(c[x] for c in rCorners), max(c[x] for c in rCorners))
		yBounds = (min(c[y] for c in rCorners), max(c[y] for c in rCorners))

		# ###---find coordinate of collision---###
		m = (prevPos[y] - pos[y]) / (prevPos[x] - pos[x])
		# calculate with x and y swapped if the slope is infinite so now the slope is 0
		if np.isinf(m):
			x, y = y, x
			m = (prevPos[y] - pos[y]) / (prevPos[x] - pos[x])
		if Plate.inInterval(planeZ, (prevPos[z], pos[z])):
			# this tells you how far along the path between the 2 points the collision occurred (assuming infinite plane)
			zRatio = abs(prevPos[z] - planeZ) / abs(prevPos[z] - pos[z])
			if np.isinf(zRatio):
				return None
			xCol = zRatio * (pos[x] - prevPos[x]) + prevPos[x]
		else:
			return None

		yCol = m * (xCol - prevPos[x]) + prevPos[y]
		if Plate.inInterval(xCol, xBounds) and Plate.inInterval(yCol, yBounds):
			colPos = np.zeros(3)
			colPos[x] = xCol
			colPos[y] = yCol
			colPos[z] = planeZ
			velocity[z] *= -1
			velocity *= dampingFactor
			return colPos, velocity
		return None

	@staticmethod
	def inInterval(x, bounds, inclusive=True):
		"""
		checks if x is between bound1 and bound2
		:param x: value to check position of
		:param bounds: tuple of the 2 bounds
		:param inclusive: whether or not both bounds are inclusive
		:return: boolean
		"""
		if inclusive:
			return min(bounds[0], bounds[1]) <= x <= max(bounds[0], bounds[1])
		return min(bounds[0], bounds[1]) < x < max(bounds[0], bounds[1])

