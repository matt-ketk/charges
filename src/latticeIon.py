import numpy as np

# from wire import Wire
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

    def checkCollision(self, prevPos, pos, velocity, dt, dampeningFactor=1):
        """
        checks for collision of a particle with this ion given its position (as a numpy array) at 2 successive iterations
        prevPos, pos, velocity: numpy arrays for all charges
        returns: the collision positions and the particles' new velocity vectors. If there was no collision returns the original position and velocity vector.
        """

        a = np.linalg.norm(velocity)**2
        b = 2 * (velocity[0] * (prevPos - self.center)[0] + velocity[1] * (prevPos - self.center)[1])
        c = np.linalg.norm(prevPos - self.center)**2 - self.radius**2

        ans = LatticeIon.quad(a, b, c)
        if not ans:
            return pos, velocity

        t0, t1 = ans
        tMax = dt

        tFiltered = set()
        if 0 < t0 < tMax:
            tFiltered.add(t0)
        if 0 < t1 < tMax:
            tFiltered.add(t1)

        if not tFiltered:
            return pos, velocity

        t = min(tFiltered)

        colPos = prevPos + velocity * t

        adjColPos = self.center + 1.01 * (colPos - self.center)
        return adjColPos, self.reflectVector(pos, velocity, colPos, dampeningFactor=dampeningFactor)


    def reflectVector(self, position, velocity, collisionPoint, dampeningFactor=1): 
        # normal vector from collision point to center
        d = collisionPoint - self.center
        norm = np.linalg.norm(d)
        n = d / norm
        # print("normal", n)
        return dampeningFactor * (velocity - 2 * n * np.dot(n, velocity))

    @staticmethod
    def quad(a, b, c):
        """
        calculates real roots of quadratic with coefficients a, b, c
        returns: tuple of roots
        """
        disc = b ** 2 - 4 * a * c

        if disc > 0:
            return (-b + disc ** 0.5)/(2 * a), (-b - disc ** 0.5)/(2 * a)
        else:
            return None

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

    @staticmethod
    def notEquals(a, b):
        return not all(a == b)

    @staticmethod
    def convMap(func, a, b):
        return np.array(list(map(func, a, b)))

    def drawIon(self, env, color=(255, 0, 0)):
        # points = [[[self.radius * np.cos(t) + self.center[0], self.radius * np.sin(t) + self.center[1]] for t in np.linspace(0, 2 * np.pi, 15)]]
        env.drawParticle(self.center, color, int(self.radius * env.zoom))

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
    prevPos = np.array([[0, 0], [2, 0], [-0.5, 0]])
    pos = np.array([[0, 0.5], [2, 5], [-0.5, 5]])
    v = pos - prevPos

    lat = LatticeIon(np.array([0, 3]), 1, 1, 1)
    print("final", lat.checkCollision(prevPos, pos, v))

if __name__ == "__main__":
    main()
