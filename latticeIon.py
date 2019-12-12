import numpy as np

from wire import Wire
from plate import Plate
from charge import Charge
from constants import Constants


class LatticeIon:
    COPD = Constants.COPPER_CUBE_SIZE / 2
    COPPER_LATTICE_UNIT = np.array(
        [
             [ COPD,  COPD,  COPD],
             [ COPD,  COPD, -COPD],
             [ COPD, -COPD, -COPD],
             [ COPD, -COPD,  COPD],
             [-COPD, -COPD,  COPD],
             [-COPD,  COPD,  COPD],
             [-COPD,  COPD, -COPD],
             [-COPD, -COPD, -COPD],
             [ COPD,     0,     0],
             [-COPD,     0,     0],
             [    0,  COPD,     0],
             [    0, -COPD,     0],
             [    0,     0,  COPD],
             [    0,     0, -COPD]
         ]
    )

    def __init__(self, center, radius, charge, dampeningFactor):
        '''
        center: numpy array of length 3 
        '''

        self.center = center
        self.radius = radius
        self.charge = charge
        self.dampeningFactor = dampeningFactor

    def reflectVector(self, particle):
        n = self.collisionNormal(particle)
        v = particle.velocity
        return self.dampeningFactor * (v - 2 * n * np.dot(n, v))

    def collisionNormal(self, particle):
        # returns a unit vector that is a surface normal at the point of collision

        v = self.collisionPoint(particle) - self.center
        return v / np.linalg.norm(v)

    def collisionPoint(self, particle):
        unitV = particle.velocity / np.linalg.norm(particle.velocity)
        origin = particle.position
        distance = 0.0

        # this is where it gets tricky... so basically,

        # between a line and a sphere, either one, two, or no intersections are formed

        # put simply, the calculations look kind of like when solving a quadratic equation, and they are...
        # ad^2 + bd + c = 0

        d = origin - self.center  # this is for readability

        # a = np.linalg.norm(unitV) ** 2 # a just becomes 1
        b = 2 * (np.dot(unitV, d))
        c = np.dot(d, d) - self.radius ** 2

        # in the case that the quadratic formula outputs two possible results, I CURRENTLY DO NOT KNOW HOW TO DETERMINE WHICH ONE IS WHICH. There is the entry point and there is the exit point (if you know what I mean) imma leave it here fo now...

        distance = -(b / 2) - np.sqrt((b / 2) ** 2 - c)
        # distance = -(b / 2) - np.sqrt(b ** 2 - c) # the alternative solution

        return (origin + distance * unitV)

    def willCollide(self, particle):
        # essentially boolean output whether particle is gonna hit given its trajectory OF THAT TICK (as in its current position + velocity * dt)
        # has the same meat and bones of the method above, but now it outputs a boolean based on the quadratic discriminant formula, b**2 -4ac >= 0.
        unitV = particle.velocity / np.linalg.norm(particle.velocity)

        d = particle.position - self.center

        a = np.linalg.norm(unitV) ** 2
        b = 2 * (np.dot(unitV, d))
        c = np.dot(d, d) - self.radius ** 2

        return (b ** 2 - 4 * (a * c) >= 0)

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
        x, y, z = (0, 1, 2)
        inX = (plate.vertex1[x] <= point[x] <= plate.vertex2[x]) or (plate.vertex1[x] >= point[x] >= plate.vertex2[x])
        inY = (plate.vertex1[y] <= point[y] <= plate.vertex2[y]) or (plate.vertex1[y] >= point[y] >= plate.vertex2[y])
        inZ = (plate.vertex1[z] <= point[z] <= plate.vertex2[z]) or (plate.vertex1[z] >= point[z] >= plate.vertex2[z])
        return inX and inY and inZ

    @staticmethod
    def generateLatticePoints(cond):
        if type(cond) is Wire:
            latticePoints = set()
            centers = []
            x, y, z = cond.getXYZ()
            for k in np.arange(cond.start[z], cond.end[z], 2 * LatticeIon.COPD):
                for j in np.arange(cond.start[y] - cond.r, cond.start[y] + cond.r, 2 * LatticeIon.COPD):
                    for i in np.arange(cond.start[x] - cond.r, cond.start[x] + cond.r, 2 * LatticeIon.COPD):
                        center = [0, 0, 0]
                        center[x] = i + LatticeIon.COPD
                        center[y] = j + LatticeIon.COPD
                        center[z] = k + LatticeIon.COPD
                        centers.append(np.array(center))

            for center in centers:
                for relPoint in LatticeIon.COPPER_LATTICE_UNIT:
                    p = center + relPoint
                    if LatticeIon.isInWire(p, cond, orientation=(x, y, z)):
                        latticePoints.add(tuple(((center + relPoint) / LatticeIon.COPD).astype(int)))
            return [LatticeIon(np.array(p) * LatticeIon.COPD, Constants.COPPER_ION_RADIUS, Constants.E, 1) for p in latticePoints]
        if type(cond) is Plate:
            latticePoints = set()
            centers = []
            x, y, z = (0, 1, 2)
            for k in np.arange(cond.vertex1[z], cond.vertex2[z], 2 * LatticeIon.COPD * np.sign(cond.vertex2[z] - cond.vertex1[z])):
                for j in np.arange(cond.vertex1[y], cond.vertex2[y], 2 * LatticeIon.COPD * np.sign(cond.vertex2[y] - cond.vertex1[y])):
                    for i in np.arange(cond.vertex1[x], cond.vertex2[x], 2 * LatticeIon.COPD * np.sign(cond.vertex2[x] - cond.vertex1[x])):
                        center = [0, 0, 0]
                        center[x] = i + LatticeIon.COPD
                        center[y] = j + LatticeIon.COPD
                        center[z] = k + LatticeIon.COPD
                        centers.append(np.array(center))
            print (centers)
            for center in centers:
                for relPoint in LatticeIon.COPPER_LATTICE_UNIT:
                    p = center + relPoint
                    if LatticeIon.isInPlate(p, cond):
                        latticePoints.add(tuple(((center + relPoint) / LatticeIon.COPD).astype(int)))
            return [LatticeIon(np.array(p) * LatticeIon.COPD, Constants.COPPER_ION_RADIUS, Constants.E, 1) for p in latticePoints]
        else:
            raise NotImplementedError("This kind of conductor is not yet supported")
    

    def compile(self):
        return (self.center, self.charge, Constants.COPPER_MASS, 1)

def main():
    """c = Charge(1.0, np.zeros(3), np.array([1, 1, 1]))
    n = LatticeIon(np.array([57, 56, 58]), 5, 0.9)
    print("point of collision:", n.collisionPoint(c))
    print("normal vector:", n.collisionNormal(c))
    print("reflected velocity vector:", n.reflectVector(c))
    """

    start = np.array([0, 0, 0])
    lengthV = np.array([3E-9, 0, 0])
    c = Wire(start, lengthV, 5E-10)
    print(LatticeIon.generateLatticePoints(c))


if __name__ == "__main__":
    main()
