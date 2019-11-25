from plate import Plate
from constants import Constants

class Capacitor:
	plate1 = None
	plate2 = None
	thickness = 0
	resistivity = 0
	
	def __init__(self, plate, gap, direction=1, dielectric=None):
		"""
		plate: the plate to start from. a copy will be made of plate offset from it
		gap > 0: the distance between the plates (in meters)
		direction: 1 or -1: direction is assumed to be along the shortest dimension of the plate
		"""
		self.plate1 = plate
		self.resistivity = plate.resistivity
		offset = direction * (gap + plate.thickness)
		self.plate2 = Plate(
			plate.vertex1 + offset,
			plate.vertex2 + offset,
			resistivity=self.resistivity
		)