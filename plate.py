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
	def rectangleCollision(prevPos, pos, rCorners, velocity, dampeningFactor=1):
		"""
		checks to see if a particle collided with the defined rectangle and gives you that particle's updated info
		:param prevPos: position of a particle at the previous iteration
		:param pos: position of that particle at the current iteration
		:param rCorners: a tuple of the corner points of a rectangle in 3D space to check for collisions
		stored such that ADJACENT CORNERS IN THE LIST ARE ADJACENT IN THE RECTANGLE
		:param velocity: the initial velocity vector of the particle
		:param dampeningFactor: the factor by which the particle's velocity decreases upon impact
		:return: the collision position of the particle and its new velocity. if no collision returns None
		"""
		basis = list()
		basis.append(rCorners[1] - rCorners[0])
		basis.append(rCorners[2] - rCorners[1])
		normal = np.cross(basis[0], basis[1])
		normal = normal / np.linalg.norm(normal)
		basis.append(normal)

		# find the value of t that satisfies both:
		# 	colPos = prevPos + velocity*t
		# 	normal dot (colPos - pointOnPlane) = 0
		# where colPos represents the collision point
		# all of this is still in the original basis
		pointOnPlane = rCorners[0]  # arbitrary
		nDotV = np.dot(normal, velocity)
		if nDotV == 0:
			return None
		t = np.dot(normal, pointOnPlane - prevPos) / nDotV
		colPos = prevPos + velocity * t

		if pos[0] - prevPos[0] != 0:
			if not Plate.inInterval(colPos[0], (pos[0], prevPos[0])):
				return None
		elif pos[1] - prevPos[1] != 0:
			if not Plate.inInterval(colPos[1], (pos[1], prevPos[1])):
				return None
		else:
			if not Plate.inInterval(colPos[2], (pos[2], prevPos[2])):
				return None

		basis = np.array(basis, dtype=float)
		changeBasisM = basis.transpose()
		invChangeBasisM = np.linalg.inv(changeBasisM)
		newRectCorners = tuple(invChangeBasisM @ corner for corner in rCorners)
		newColPos = invChangeBasisM @ colPos

		xBounds = (min(c[0] for c in newRectCorners), max(c[0] for c in newRectCorners))
		yBounds = (min(c[1] for c in newRectCorners), max(c[1] for c in newRectCorners))

		if Plate.inInterval(newColPos[0], xBounds) and Plate.inInterval(newColPos[1], yBounds):
			newVelocity = dampeningFactor * (velocity - 2 * normal * np.dot(normal, velocity))
			return colPos, newVelocity
		return None

	@staticmethod
	def inInterval(x, bounds, inclusive=False):
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


def main():
	prevPos = np.array([1, 3, 1])
	pos = np.array([1, 3, 5])
	v = pos - prevPos
	print("velocity:", v)
	rCorners = np.array([[1, 1, 1], [3, 3, 3], [1, 5, 5], [-1, 3, 3]])
	print("rectangle corners:", rCorners)
	print("collision:", Plate.rectangleCollision(prevPos, pos, rCorners, v))


if __name__ == "__main__":
	main()
