import numpy as np

from constants import Constants
from conductor import Conductor

class Wire(Conductor):
	def __init__(self, start, lengthV, r, resistivity=Constants.COPPER_RESISTIVITY):
		"""
		start: length-3 numpy array representing the starting coordinates of the cylindrical wire
		lengthV: length-3 numpy array (vector) that represents the length of the wire. 2 of the entries should be zero to make the cylinder orthogonal to the basis
		r: radius of wire
		"""
		super(Wire, self).__init__(resistivity=resistivity)

		if len([i for i in lengthV if i != 0]) != 1:
			print("Warning: non-orthogonal wires are not yet supported")
		self.start = start
		self.lengthV = lengthV
		self.end = start + lengthV
		self.r = r


	def checkCollision(self, pos, prevPos):
		"""
		checks for collision of a particle with the cylinder given its position (as a numpy array) at 2 successive iterations
		returns: the position of the collision and the particle's new velocity vector. If there was no colision returns None.
		"""

		# find orientation of wire
		z = np.argmax(np.abs(self.lengthV))
		if z != 0:
			x = 0
			if z != 1:
				y = 1
			else:
				y = 2
		else:
			x = 1
			y = 2

		m = (prevPos[y] - pos[y]) / (prevPos[x] - pos[x])
		A = -m * prevPos[x] + pos[y] - self.start[y]
		xCol = Wire.quad(1 + m**2, 2*(A*m + self.start[x]), A**2 + self.start[x]**2 - self.r**2)
		
		if len(xCol) == 0:
			return None
		if len(xCol) == 1:
			xCol = xCol[0]
			withinSegment = min(prevPos[x], pos[x]) <= xCol <= max(prevPos[x], pos[x])
		else:
			withinSegment0 = min(prevPos[x], pos[x]) <= xCol[0] <= max(prevPos[x], pos[x])
			withinSegment1 = min(prevPos[x], pos[x]) <= xCol[1] <= max(prevPos[x], pos[x])
			if withinSegment0:
				xCol = xCol[0]
			if withinSegment1:
				xCol = xCol[1]
		
		yCol = m*(xCol - prevPos[x]) + prevPos[y]
		
		mz = (prevPos[z] - pos[z]) / (prevPos[x] - pos[x])
		zCol = m*(xCol - prevPos[x]) + prevPos[y]
		withinHeight = zCol - self.start[z] <= self.lengthV[z]

		if not withinHeight or not withinSegment:
			return None
	  pos = np.zeros(3)
		pos[z] = zCol
		pos[x] = xCol
		pos[y] = yCol
		return pos, np.zeros(3)
	
	@staticmethod
	def quad(a, b, c):
		"""
		calculates real roots of quadratic with coefficients a, b, c
		returns: tuple of roots
		"""
		disc = b**2 - 4*a*c
		if disc < 0:
			return tuple()
		if disc == 0:
			return (-b / (2*a),)
		return ((-b + disc**0.5) / (2*a), (-b - disc**0.5) / (2*a))
