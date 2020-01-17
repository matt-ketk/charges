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
        # TODO: make the input an array of positions/velocities and use numpy operations rather than a for loop
        
        a = np.sum(np.square(velocity), axis=1)
        b = 2 * np.sum(velocity * (prevPos - self.center), axis=1)
        c = np.sum(np.square(prevPos - self.center), axis=1) - self.radius**2

        t0,  t1 = LatticeIon.quad(a,b,c)
        # print("t0", t0)
        # print("t1", t1)
        tMax = np.ones(len(pos)) * dt
        # print("tmax", tMax)

        t0filtered = np.where((t0 > 0) * (tMax > t0), t0, 2 * dt)
        t1filtered = np.where((t1 > 0) * (tMax > t1), t1, 2 * dt)
        #
        # print("t0 filtered", t0filtered)
        # print("t1 filtered", t1filtered)

        t = np.where((t0filtered != 2 * dt) + (t1filtered != 2 * dt), np.minimum(t0filtered, t1filtered), np.nan)  #TODO: fix logic
         
        # print("t", t)

        colPos = np.where(np.repeat(np.isfinite(t), 2).reshape((len(t), 2)), prevPos + LatticeIon.convMap(np.multiply, velocity, t), pos)  # TODO: figure out reshaping

        return colPos, self.reflectVector(pos, velocity, colPos, dampeningFactor=dampeningFactor)


    def reflectVector(self, position, velocity, collisionPoint, dampeningFactor=1): 
        d = collisionPoint - self.center

        # array of normal vectors from collision point to center
        norm = np.linalg.norm(d, axis=1)

        # for each particle, check if there is a collision point and if so, take normal vector; if not, take np.nan
        n = np.where(np.repeat(LatticeIon.convMap(LatticeIon.notEquals, collisionPoint, position), 2).reshape((len(position), 2)), d / norm.reshape((len(d), 1)), np.nan)
        # print("normal", n)
        return np.where(np.isfinite(n), dampeningFactor * (velocity - 2 * LatticeIon.convMap(np.multiply, n, np.sum(np.multiply(n, velocity), axis=1))), velocity)  # TODO fix this :((

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

    @staticmethod
    def notEquals(a, b):
        return not all(a == b)

    @staticmethod
    def convMap(func, a, b):
        return np.array(list(map(func, a, b)))

    def drawIon(self, env, color=(255, 0, 0)):
        # points = [[[self.radius * np.cos(t) + self.center[0], self.radius * np.sin(t) + self.center[1]] for t in np.linspace(0, 2 * np.pi, 15)]]
        env.drawParticle(self.center, color, int(self.radius * env.zoom * 1 * 1.2))

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
