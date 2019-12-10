import numpy as np
from constants import Constants

def forces(coords, charges):
  '''
  coords: a 2D numpy array with columns X, Y, Z and rows for each particle
    [ [x1, y1, z1],
      [x2, y2, z2],
      ...,
      xN, yN, zN] ]
  charges: a numpy array containing the charges of each particle

  returns a 2D numpy (float) array of size N where the ith row contains the vector sum of the other particles' electrostatic forces [Fx, Fy, Fz] on the ith particle
  '''

  N = coords.shape[0] # number of particles

  # find the distance vectors between each pair of particles
	# np.array([coords,] * N, dtype = float) is equivalent to np.repeat(coords, N, axis=0)
  coordsRow = np.transpose(np.array([coords,] * N, dtype = float), axes = [2, 1, 0])
  coordsCol = np.transpose(coordsRow, axes = [0, 2, 1])
  distanceVectors = np.array(coordsRow - coordsCol) # distance vectors in x, y, and z for each particle to each other particle

  # find the distances (magnitude only) between each pair of particles
  distances = np.sum(distanceVectors ** 2, axis = 0, keepdims = True) ** 0.5

  # find the product of the charges of each pair of particles (3D array mirroring the distanceVectors array)
  chargesCol = np.transpose([charges])

  # chargeProducts -> 1 x n x n 
  chargeProducts = np.matmul(chargesCol, np.array([charges]))[np.newaxis, :, :]

  #print(distanceVectors.shape, chargeProducts.shape, distances.shape)

  # find forces between each pair of particles (3D array)
  forcesSplit = -1 * Constants.K * distanceVectors * np.divide(chargeProducts, distances ** 3 + 1E-20) # added 1E-12 to avoid divide by zero error 

  # find forces on each particle from by summing the forces from all the other particles from forcesSplit (2D array)
  forces = np.sum(forcesSplit, axis = 2)

  return -forces

def deltaPosition(dt, forcesList, velList, massList):
  '''
  dt: the timestep (which will be determined by pygame)
  forcesList: a 2D (shape Nx3) numpy array of forces on each particle, the ith row contains [Fx, Fy, Fz] for the ith particle
  velList: a 2D (shape Nx3) numpy array of the original velocities of each particle, [vx, vy, vz]
  massList: a numpy array of length N of the masses of each particle 

  returns a 2D (shape Nx3) of the change in positions of each particle 
  '''
  # make Nx1 array of masses
  # EDIT: Deleted the transpose, so shape is now 1 x N
  massList = np.array([massList])

  # calculate acclerations from F=ma
  accelerations = np.divide(forcesList, massList)

  # calculate the change in position by using deltaPos = 0.5at^2 + vt
  deltaPositions = 0.5 * dt**2 * accelerations + velList * dt
  return deltaPositions
  
def deltaVelocity(dt, forcesList, massList):
  '''
  dt: the timestep (which will be determined by pygame)
  forcesList: a 2D (shape Nx3) numpy array of forces on each particle, the ith row contains [Fx, Fy, Fz] for the ith particle
  massList: a numpy array of length N of the masses of each particle 

  returns a 2D (shape Nx3) of the change in positions of each particle 
  '''
  # make Nx1 array of masses
  # EDIT: Deleted the transpose, so shape is now 1 x N
  massList = np.array([massList])

  # calculate acceleration from F=ma
  accelerations = np.divide(forcesList, massList)

  # calculate change in velocity
  deltaVelocities = dt * accelerations
  return deltaVelocities  