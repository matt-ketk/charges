from constants import Constants

class Conductor:

	def __init__(self, resistivity=Constants.COPPER_RESISTIVITY):
		self.resistivity = resistivity