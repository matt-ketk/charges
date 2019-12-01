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
		self.faces.append(tuple(self.vertices[i] for i in [0, 1, 3, 2]))
		self.faces.append(tuple(self.vertices[i] for i in [4, 5, 7, 6]))
		self.faces.append(tuple(self.vertices[i] for i in [0, 1, 5, 4]))
		self.faces.append(tuple(self.vertices[i] for i in [2, 3, 7, 6]))
		self.faces.append(tuple(self.vertices[i] for i in [0, 2, 6, 4]))
		self.faces.append(tuple(self.vertices[i] for i in [1, 3, 7, 5]))

		super(Plate, self).__init__(resistivity=resistivity)

	def checkCollision(self, prevPos, pos, velocity, dampingFactor=1):
		"""
		checks to see if the specified particle collided with any of the faces of the plate
		:param prevPos: previous position vector of the particle
		:param pos: current position of the particle
		:param velocity: current velocity vector of the particle
		:param dampingFactor: the factor by which the velocity decreases upon impact
		:return: the position of the collision, the particle's new velocity vector, the corners of face the particle
		collided with. If there was no collision, returns None
		"""
		collisions = []
		for face in self.faces:
			col = Plate.rectangleCollision(prevPos, pos, face, velocity, dampingFactor=dampingFactor)
			if col is not None:
				collisions.append((face, col))

		if len(collisions) == 0:
			return None

		col = min(collisions, key=lambda x: abs(x[1][0][0] - prevPos[0]))  # the actual collision is the first one
		return col[1][0], col[1][1], col[0]

	@staticmethod
	def rectangleCollision(prevPos, pos, rCorners, velocity, dampingFactor=1):
		"""
		checks to see if a particle collided with the defined rectangle and gives you that particles updated info
		rectangle must be parallel to xy, xz, or yz for now
		:param prevPos: position of a particle at the previous iteration
		:param pos: position of that particle at the current iteration
		:param rCorners: a tuple of the corner points of a rectangle in 3D space to check for collisions
		stored such that ADJACENT CORNERS IN THE LIST ARE ADJACENT IN THE RECTANGLE
		:param velocity: the initial velocity vector of the particle
		:param dampingFactor: the factor by which the particle's velocity decreases upon impact
		:return: the collision position of the particle and its new velocity. if no collision returns None
		"""
		# find basis for the plane that the rectangle is parallel to
		basis = list()
		basis.append(rCorners[1] - rCorners[0])
		basis.append(rCorners[2] - rCorners[1])
		# so that the new space obeys right hand rule and can be represented like cartesian with z point out of the page
		basis.append(np.cross(basis[0], basis[1]))
		basis = np.array(basis, dtype=float)
		for i in range(len(basis)):
			basis[i] /= np.linalg.norm(basis[i])
		changeBasisM = basis.transpose()
		invChangeBasisM = np.linalg.inv(changeBasisM)
		newRectCorners = tuple(invChangeBasisM @ corner for corner in rCorners)
		planeZ = newRectCorners[0][2]

		nPrevPos = invChangeBasisM @ prevPos
		nPos = invChangeBasisM @ pos

		xBounds = (min(c[0] for c in newRectCorners), max(c[0] for c in newRectCorners))
		yBounds = (min(c[1] for c in newRectCorners), max(c[1] for c in newRectCorners))

		# ###---find coordinate of collision---###
		m = (nPrevPos[1] - nPos[1]) / (nPrevPos[0] - nPos[0])
		# calculate with x and y swapped if the slope is infinite so now the slope is 0
		if np.isinf(m):
			permuteXY = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
			changeBasisM = permuteXY @ changeBasisM
			invChangeBasisM = invChangeBasisM @ permuteXY
			newRectCorners = tuple(invChangeBasisM @ corner for corner in rCorners)
			nPrevPos = invChangeBasisM @ prevPos
			nPos = invChangeBasisM @ pos
			xBounds = (min(c[0] for c in newRectCorners), max(c[0] for c in newRectCorners))
			yBounds = (min(c[1] for c in newRectCorners), max(c[1] for c in newRectCorners))

			m = (nPrevPos[1] - nPos[1]) / (nPrevPos[0] - nPos[0])
		if Plate.inInterval(planeZ, (nPrevPos[2], nPos[2])):
			# this tells you how far along the path between the 2 points the collision occurred (assuming infinite plane)
			zRatio = abs(nPrevPos[2] - planeZ) / abs(nPrevPos[2] - nPos[2])
			if np.isinf(zRatio):
				return None
			xCol = zRatio * (nPos[0] - nPrevPos[0]) + nPrevPos[0]
		else:
			return None

		yCol = m * (xCol - nPrevPos[0]) + nPrevPos[1]
		if Plate.inInterval(xCol, xBounds) and Plate.inInterval(yCol, yBounds):
			colPos = np.zeros(3)
			colPos[0] = xCol
			colPos[1] = yCol
			colPos[2] = planeZ
			originalBasisColPos = changeBasisM @ colPos

			newVelocity = dampingFactor * invChangeBasisM @ velocity
			newVelocity[2] *= -1
			originalBasisNewVel = changeBasisM @ newVelocity
			return originalBasisColPos, originalBasisNewVel
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
