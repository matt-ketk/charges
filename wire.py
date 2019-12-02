import numpy as np

from constants import Constants
from conductor import conductor

class Wire(conductor):

	def __init__(start, lengthV, r, resistivity=Constants.COPPER_RESISTIVITY):
		"""
		start: length-3 numpy array representing the starting coordinates of the cylindrical wire
		lengthV: length-3 numpy array (vector) that represents the length of the wire. 2 of the entries should be zero to make the cylinder orthogonal to the basis
		r: radius of wire
		"""
		if len(i for i in lengthV if i != 0) != 1:
			print("Warning: non-orthogonal wires are not yet supported")
		self.start = start
		self.lengthV = lengthV
		self.end = start + lengthV
		self.r = r

		super(Wire, self).__init__(resistivity=resistivity)

	def checkCollision(self, pos, prevPos):
		"""
		checks for collision of a particle with the cylinder given its position (as a numpy array) at 2 successive iterations
		returns: the position of the collision and the particle's new velocity vector. If there was no colision returns None.
		"""
	
		

