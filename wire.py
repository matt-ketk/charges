import numpy as np

from constants import Constants
from conductor import Conductor
from plate import Plate

from pyquaternion import Quaternion as Q


class Wire(Conductor):
    def __init__(self, start, lengthV, r, resistivity=Constants.COPPER_RESISTIVITY):
        """
        start: length-3 numpy array representing the starting coordinates of the cylindrical wire
        lengthV: length-3 numpy array (vector) that represents the length of the wire.
        r: radius of wire
        """
        super(Wire, self).__init__(resistivity=resistivity)

        self.start = start
        self.lengthV = lengthV
        self.end = start + lengthV
        self.r = r
        
        self.coords = []
        self.charges = []
        self.masses = []
        self.stationary = []
        

    def getXYZ(self):
        """
        :return: (x, y, z), where z is the axis aligned with the length of the cylinder
        """
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
        return x, y, z

    def checkCollision(self, prevPos, pos, vel, dampeningFactor=1):
        """
        checks for collision of a particle with the cylinder given its position (as a numpy array) at 2 successive iterations
        :param prevPos: particle's previous position
        :param pos: particle's current position
        :param vel: particle's velocity
        :param dampeningFactor: factor by which the collision should slow down the velocity of the particle
        :return: position of the collision and a 2D numpy array of particle's new velocity vector
        returns None if no collision
        """
        # checks for collision of a particle with the cylinder given its position (as a numpy array) at 2 successive iterations
        # :param prevPos: 2D numpy array where each row represents a particle's previous position
        # :param pos: 2D numpy array where each row represents a particle's current position
        # :param vel: 2D numpy array where each row represents a particle's velocity
        # :param dampeningFactor: factor by which the collisions should slow down the velocity of the particle
        # :return: 2D numpy array of positions of the collision (pos[i] if no collision)
        # and a 2D numpy array of each particle's new velocity vectors (vel[i] if no collision)

        # if the point didn't move ignore it
        if all(pos == prevPos):
            return None

        basis = list()
        # generate some random vector that is not in the same direction as lengthV
        # randMatrix = np.array([[1, 1, 0],
        #                        [0, 1, 2],
        #                        [-1, -1, 2]])
        # randVec = randMatrix @ self.lengthV
        randVec = np.array([0, 0, 1])
        basis.append(-np.cross(self.lengthV, randVec))
        basis[0] /= np.linalg.norm(basis[0])
        basis.append(np.cross(self.lengthV, basis[0]))
        basis[1] /= np.linalg.norm(basis[1])
        basis.append(self.lengthV / np.linalg.norm(self.lengthV))
        basis = np.array(basis, dtype=float)
        invChangeBasisM = basis.transpose()
        changeBasisM = np.linalg.inv(invChangeBasisM)

        nPrevPos = changeBasisM @ prevPos
        nPos = changeBasisM @ pos
        nStart = changeBasisM @ self.start
        nEnd = changeBasisM @ self.end
        nVel = changeBasisM @ vel

        tempColPos = np.zeros(3)
        # ###---find x coordinate of collision (xCol) if its within the circle---###
        m = (nPrevPos[1] - nPos[1]) / (nPrevPos[0] - nPos[0])
        # calculate with x and y swapped if the slope is infinite so now the slope is 0

        a = -m * nPrevPos[0] + nPrevPos[1] - nStart[1]
        mz = (nPrevPos[2] - nPos[2]) / (nPrevPos[0] - nPos[0])
        # if np.isinf(mz):
        #     return None

        xColTemp = Wire.quad(1 + m ** 2, 2 * (a + nStart[0]), a ** 2 + nStart[0] ** 2 - self.r ** 2)

        withinSegment = False
        if len(xColTemp) == 0 or len(xColTemp) == 1:
            return None
        else:
            xColTemp = [xc for xc in xColTemp if Plate.inInterval(xc, (nPrevPos[0], nPos[0]))]
            if len(xColTemp) != 0:
                withinSegment = True
                tempColPos[0] = min(xColTemp, key=lambda p: abs(p - nPrevPos[0]))

        if not withinSegment:
            return None

        # find y coordinate of collision (yCol)
        tempColPos[1] = m * (tempColPos[0] - nPrevPos[0]) + nPrevPos[1]

        # find z coordinate of collision (zCol) and check if it's within the length vector of the cylinder
        tempColPos[2] = mz * (tempColPos[0] - nPrevPos[0]) + nPrevPos[2]
        withinHeight = min(nStart[2], nEnd[2]) <= tempColPos[2] <= max(nStart[2], nEnd[2])

        if not withinHeight:
            return None

        # find normal vector to the cylinder to calculate the new velocity
        normal = np.zeros(3)
        normal[0] = tempColPos[0] - nStart[0]
        normal[1] = tempColPos[1] - nStart[1]
        normal[2] = 0
        normal /= np.linalg.norm(normal)

        tempVel = dampeningFactor * (nVel - 2 * normal * np.dot(normal, nVel))

        outputVel = invChangeBasisM @ tempVel
        colPos = invChangeBasisM @ tempColPos

        return colPos, outputVel
    
    def compile(self, list):
        self.coords = []
        self.charges = []
        self.masses = []
        self.stationary = []
        
        for entry in list:
            self.coords.append(list[0])
            self.charges.append(list[1])
            self.masses.append(list[2])
            self.stationary.append(list[3])
        return self.coords, self.charges, self.masses, self.coords

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

    def draw(self, environment, color, thickness = 1):
        l = self.lengthV
        ln = np.sqrt(np.sum(l ** 2))

        # Form the rotation axis by getting the cross product between the z-axis vector and the normalized length vector
        rotationAxis = np.cross(np.array([0, 0, 1]), l / ln)


        # get the magnitude of rotationAxis ( sin(theta) )
        # normalize rotationAxis
        rotationMag = np.sqrt(np.sum(rotationAxis ** 2))
        
        angle = np.math.asin(rotationMag)

        # Create a quaternion
        if l[2] < 0:
            angle = np.pi/2 - angle
        if rotationMag == 0:
            angle = 0
            rotationAxis_ = np.array([0,0,1])
        else:
            rotationAxis_ = rotationAxis / rotationMag
        rotationQ = Q(axis = rotationAxis_, angle = angle)

        # Draw two circles centered around the point and rotate them by the rotation quaternion
        points = [[], []]
        for i in range (20):
            p = rotationQ.rotate(np.array([np.sin(i/20 * 2 * 3.1415) * self.r, np.cos(i/20 * 2 * 3.1415) * self.r, 0]))
            points[0].append(p + self.start)

            p = rotationQ.rotate(np.array([np.sin(i/20 * 2 * 3.1415) * self.r, np.cos(i/20 * 2 * 3.1415) * self.r, ln]))
            points[1].append(p + self.start)
        
        environment.drawMesh(points, color)

        tpoints = environment.transformPoint(np.array([self.start, self.end]))

        if (abs(tpoints[1][0] - tpoints[0][0]) < 0.01):
            offset = np.array([self.r, 0]) * environment.zoom 
        elif (abs(tpoints[1][1] - tpoints[0][1]) < 0.01):
            offset = np.array([0, self.r]) * environment.zoom 

        else:
            slope = (tpoints[1][1] - tpoints[0][1]) / (tpoints[1][0] - tpoints[0][0])
            perp = - 1 / slope

            norm = np.sqrt(1 + 1 / slope ** 2)

            offset = np.array([self.r / norm, perp * self.r / norm]) * environment.zoom 

        l1 = np.array([tpoints[0] + offset, tpoints[1] + offset]).astype(int)
        l2 = np.array([tpoints[0] - offset, tpoints[1] - offset]).astype(int)
        
        points = [[l1[0], l1[1]], [l2[0], l2[1]]]
        environment.drawMeshRaw(points, color)
