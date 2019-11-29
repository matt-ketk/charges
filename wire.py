import numpy as np

from constants import Constants
from conductor import Conductor
from plate import Plate

class Wire(Conductor):
    def __init__(self, start, lengthV, r, resistivity=Constants.COPPER_RESISTIVITY):
        """
        start: length-3 numpy array representing the starting coordinates of the cylindrical wire
        lengthV: length-3 numpy array (vector) that represents the length of the wire. 2 of the entries should be zero to make the cylinder orthogonal to the basis
        r: radius of wire
        """
        super(Wire, self).__init__(resistivity=resistivity)

        if len([i for i in lengthV if i != 0]) != 1:
            print("Warning: non-orthogonal wires are not yet supported")
        self.start = start
        self.lengthV = lengthV
        self.end = start + lengthV
        self.r = r

    def checkCollision(self, prevPos, pos):
        """
        checks for collision of a particle with the cylinder given its position (as a numpy array) at 2 successive iterations
        returns: the position of the collision and the particle's new velocity vector. If there was no colision returns None.
        """

        # if the point didn't move ignore it
        if all(pos == prevPos):
            return None

        # find orientation of wire
        z = np.argmax(np.abs(self.lengthV))
        if z != 0:
            x = 0
            if z != 1:
                y = 1
            else:
                y = 2
        else:
            x = 1
            y = 2

        # ###---find x coordinate of collision (xCol) if its within the circle---###
        m = (prevPos[y] - pos[y]) / (prevPos[x] - pos[x])
        # calculate with x and y swapped if the slope is infinite so now the slope is 0
        if np.isinf(m):
            x, y = y, x
            m = (prevPos[y] - pos[y]) / (prevPos[x] - pos[x])

        A = -m * prevPos[x] + prevPos[y] - self.start[y]
        mz = (prevPos[z] - pos[z]) / (prevPos[x] - pos[x])
        if np.isinf(mz):
            return None

        xCol = Wire.quad(1 + m ** 2, 2 * (A * m + self.start[x]), A ** 2 + self.start[x] ** 2 - self.r ** 2)

        withinSegment = False
        if len(xCol) == 0:
            return None
        if len(xCol) == 1:
            xCol = xCol[0]
            withinSegment = Plate.inInterval(xCol, (prevPos[x], pos[x]))
        elif Plate.inInterval(xCol[0], (prevPos[x], pos[x])):
            xCol = xCol[0]
            withinSegment = True
        elif Plate.inInterval(xCol[1], (prevPos[x], pos[x])):
            xCol = xCol[1]
            withinSegment = True

        if not withinSegment:
            return None

        # find y coordinate of collision (yCol)
        yCol = m * (xCol - prevPos[x]) + prevPos[y]

        # find z coordinate of collision (zCol) and check if it's within the length vector of the cylinder
        zCol = mz * (xCol - prevPos[x]) + prevPos[z]
        withinHeight = min(0, self.lengthV[z]) <= zCol <= max(0, self.lengthV[z])

        if not withinHeight:
            return None

        colPos = np.zeros(3)
        colPos[z] = zCol
        colPos[x] = xCol
        colPos[y] = yCol
        return colPos, np.zeros(3)


    @staticmethod
    def quad(a, b, c):
        """
        calculates real roots of quadratic with coefficients a, b, c
        returns: tuple of roots
        """
        disc = b ** 2 - 4 * a * c
        if disc < 0:
            return tuple()
        if disc == 0:
            return -b / (2 * a),
        return (-b + disc ** 0.5) / (2 * a), (-b - disc ** 0.5) / (2 * a)
