import numpy as np

from wire import Wire
from plate import Plate
from charge import Charge
from constants import Constants
from plate import Plate


class LatticeIon:

    def __init__(self, center, radius, charge, mass):
        '''
        center: numpy array of length 2
        '''

        self.center = center
        self.radius = radius
        self.charge = charge
        self.mass = mass

    def checkCollision(self, prevPos, pos, velocity, dampeningFactor=.9):
        """
        checks for collision of a particle with this ion given its position (as a numpy array) at 2 successive iterations
        prevPos, pos, velocity: numpy arrays for all charges
        returns: the collision positions and the particles' new velocity vectors. If there was no collision returns the original position and velocity vector.
        """
        # TODO: make the input an array of positions/velocities and use numpy operations rather than a for loop
        
        a = np.sum(np.square(velocity), axis=1)
        b = 2 * np.sum(velocity * (prevPos - self.center), axis=1)
        c = np.sum(np.square(prevPos - self.center), axis=1) - self.radius**2

        t0,  t1 = LatticeIon.quad(a,b,c)
        print (t0)
        print (t1)
        tMax = np.linalg.norm(pos - prevPos, axis=1) / np.linalg.norm(velocity, axis=1)
        print (tMax)

        t0In = np.where((t0 > 0) * (tMax > t0), t0, 0)
        t1In = np.where((t1 > 0) * (tMax > t1), t1, 0)

        print (t0In)

        t = np.where((t0 > 0) * (t1 > 0), np.minimum(t0In, t1In), np.nan)
         
        print (t)

        colPos = np.where(np.isfinite(t), prevPos + np.array(list(map(np.multiply, velocity, t))), pos)  # TODO: figure out reshaping

        return colPos, reflectVector(pos, velocity, colPos, dampeningFactor=dampeningFactor) 


    def reflectVector(self, position, velocity, collisionPoint, dampeningFactor=1): 
        d = collisionPoint - self.center

        # array of normal vectors from collision point to center
        norm = np.linalg.norm(d)

        # for each particle, check if there is a collision point and if so, take normal vector; if not, take np.nan
        n = np.where(collisionPoint != position, d / norm, np.nan)  
        print (n)
        return np.where(np.isfinite(n), dampeningFactor * (velocity - 2 * n * map(np.dot, n, velocity)), velocity)
   
    @staticmethod
    def quad(a, b, c):
        """
        calculates real roots of quadratic with coefficients a, b, c
        returns: tuple of np arrays of roots
        """
        disc = b ** 2 - 4 * a * c

        t1 = np.divide(-b + disc ** 0.5, 2 * a, out=np.zeros(len(disc)), where=disc>0)
        t2 = np.divide(-b - disc ** 0.5, 2 * a, out=np.zeros(len(disc)), where=disc>0)
        return (t1, t2)


    @staticmethod
    def isInWire(point, wire, orientation=None):
        if orientation is None:
            x, y, z = wire.getXYZ()
        else:
            x, y, z = orientation
        inLength = (min(wire.start[z], wire.end[z]) + Constants.COPPER_ION_RADIUS) \
                    <= point[z] \
                    <= (max(wire.start[z], wire.end[z]) - Constants.COPPER_ION_RADIUS)

        inCircle = ((point[x] - wire.start[x]) ** 2 + (point[y] - wire.start[y]) ** 2) ** 0.5 + Constants.COPPER_ION_RADIUS < wire.r
        return inLength and inCircle

    @staticmethod
    def isInPlate(point, plate):
        # this does not work for non-orthogonal plates :(
        x, y, z = (0, 1, 2)
        inX = (plate.vertex1[x] <= point[x] <= plate.vertex2[x]) or (plate.vertex1[x] >= point[x] >= plate.vertex2[x])
        inY = (plate.vertex1[y] <= point[y] <= plate.vertex2[y]) or (plate.vertex1[y] >= point[y] >= plate.vertex2[y])
        inZ = (plate.vertex1[z] <= point[z] <= plate.vertex2[z]) or (plate.vertex1[z] >= point[z] >= plate.vertex2[z])
        return inX and inY and inZ

    @staticmethod
    def generateLatticePoints(cond, radius = Constants.COPPER_ION_RADIUS, charge = Constants.E, mass = Constants.COPPER_MASS, offset = np.array([0.,0.,0.])):
        if type(cond) is Wire:
            latticePoints = set()
            centers = []
            x, y, z = cond.getXYZ()
            for k in np.arange(cond.start[z], cond.end[z], 2 * LatticeIon.COPD):
                for j in np.arange(cond.start[y] - cond.r, cond.start[y] + cond.r, 2 * LatticeIon.COPD):
                    for i in np.arange(cond.start[x] - cond.r, cond.start[x] + cond.r, 2 * LatticeIon.COPD):
                        center = [0,0,0]
                        center[x] += i + LatticeIon.COPD
                        center[y] += j + LatticeIon.COPD
                        center[z] += k + LatticeIon.COPD
                        centers.append(np.array(center))
            print (len(centers))
            print (centers)
            for center in centers:
                for relPoint in LatticeIon.COPPER_LATTICE_UNIT:
                    p = center + relPoint
                    if LatticeIon.isInWire(p, cond, orientation=(x, y, z)):
                        point = tuple(((center + relPoint) / LatticeIon.COPD * 100).astype(int))
                        if point in latticePoints:
                            print (point)
                        latticePoints.add(point)
                    
            return [LatticeIon(np.array(p) * LatticeIon.COPD / 100 + offset, radius, charge, mass, 1) for p in latticePoints]
        if type(cond) is Plate:
            latticePoints = set()
            centers = []
            x, y, z = (0, 1, 2)
            for k in np.arange(cond.vertex1[z], cond.vertex2[z], 2 * LatticeIon.COPD * np.sign(cond.vertex2[z] - cond.vertex1[z])):
                for j in np.arange(cond.vertex1[y], cond.vertex2[y], 2 * LatticeIon.COPD * np.sign(cond.vertex2[y] - cond.vertex1[y])):
                    for i in np.arange(cond.vertex1[x], cond.vertex2[x], 2 * LatticeIon.COPD * np.sign(cond.vertex2[x] - cond.vertex1[x])):
                        center = np.copy(offset)
                        center[x] += i + LatticeIon.COPD
                        center[y] += j + LatticeIon.COPD
                        center[z] += k + LatticeIon.COPD
                        centers.append(np.array(center))
            print (centers)
            print (len(centers))
            for center in centers:
                for relPoint in LatticeIon.COPPER_LATTICE_UNIT:
                    p = center + relPoint
                    if LatticeIon.isInPlate(p, cond):
                        latticePoints.add(tuple(((center + relPoint) / LatticeIon.COPD).astype(int)))
            return [LatticeIon(np.array(p) * LatticeIon.COPD, radius, charge, mass, 1) for p in latticePoints]
        else:
            raise NotImplementedError("This kind of conductor is not yet supported")
    

    def compile(self):
        return (self.center, self.charge, self.mass, 1)

def main():
    """c = Charge(1.0, np.zeros(3), np.array([1, 1, 1]))
    n = LatticeIon(np.array([57, 56, 58]), 5, 0.9)
    print("point of collision:", n.collisionPoint(c))
    print("normal vector:", n.collisionNormal(c))
    print("reflected velocity vector:", n.reflectVector(c))
    """

    # start = np.array([0, 0, 0])
    # lengthV = np.array([3E-9, 0, 0])
    # c = Wire(start, lengthV, 5E-10)
    # print(LatticeIon.generateLatticePoints(c))
    prevPos = np.array([[0, 0], [0.5, 0], [-0.5, 0]])
    pos = np.array([[0,5], [0.5, 5], [-0.5, 5]])
    v = pos - prevPos

    lat = LatticeIon(np.array([0,3]), 1, 1, 1)
    print(lat.checkCollision(prevPos, pos, v))

if __name__ == "__main__":
    main()
