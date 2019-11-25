from constants import Constants
import numpy as np

class Charge:
  def __init__(self, charge, position = np.zeros(3), velocity = np.zeros(3)):
    self.charge = charge * Constants.E
    self.mass = Constants.MASS_ELECTRON
    self.position = position
    self.velocity = velocity

  
