from constants import Constants
import numpy as np

class Charge:
  massList = np.array()
  posList = np.ndarray((, 3))
  posXList = np.array()
  posYList = np.array()
  posZList = np.array()

  chargeList = np.array()
  
  def __init__(self, charge, position = np.zeros(3), velocity = np.zeros(3)):
    '''
    charge: charge of the particle in Coulombs
    position: 3D numpy array for position vector
    velocity: 3D numpy array for velocity vector
    '''
    self.charge = charge * Constants.E
    self.mass = Constants.MASS_ELECTRON
    self.position = position
    self.prevPosition = np.zeros(3)
    self.velocity = velocity

    # every time charge obj is init'ed, add the info to a class-wide list
    np.stack()
    np.append(self.posXList, self.position[0])
    np.append(self.posYList, self.position[1])    np.append(self.posZList, self.position[2])
    np.append(self.massList, self.mass)


    
    



  
