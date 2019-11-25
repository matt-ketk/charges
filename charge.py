from constants import Constants

class Charge:
  def __init__(self, charge):
    self.charge = charge * Constants.E
    self.mass = Constants.MASS_ELECTRON