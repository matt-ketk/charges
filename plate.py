import numpy as np

from conductor import Conductor
from constants import Constants

class Plate(Conductor):
	vertex1 = tuple()  # 3D coordinate of one vertex
	vertex2 = tuple()  # 3D coordinate of the opposing vertex
	vertices = list()  # list of all 8 vertices
	shape = tuple()  # tuple of positive floats indicating dimensions
	thicknessDimension = 0  # 0, 1, or 2 indicating which dimension of the plate is its thickness (the index of the smallest value in shape)
	thickness = 0

	def __init__(self, vertex1, vertex2, resistivity=Constants.COPPER_RESISTIVITY):
		# input the 3D coordinate numpy arrays of 2 opposing vertices of the plate
		# generates all the corners of the prism
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

		super(Plate, self).__init__(resistivity=resistivity)

