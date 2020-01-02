import numpy as np

from constants import Constants
from conductor import Conductor
from plate import Plate

from pyquaternion import Quaternion as Q


class Wire_old(Conductor):
    def __init__(self, start, lengthV, r, resistivity=Constants.COPPER_RESISTIVITY):
        """
        start: length-3 numpy array representing the starting coordinates of the cylindrical wire
        lengthV: length-3 numpy array (vector) that represents the length of the wire. 2 of the entries should be zero to make the cylinder orthogonal to the basis
        r: radius of wire
        """
        super(Wire, self).__init__(resistivity=resistivity)

        if len([i for i in lengthV if i != 0]) != 1:
            raise NotImplementedError("non-orthogonal wires are not yet supported")
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
        returns: the position of the collision and the particle's new velocity vector. If there was no colision returns None.
        """
        #TODO: make the input an array of positions/velocities and use numpy operations rather than a for loop

        # if the point didn't move ignore it
        if all(pos == prevPos):
            return None

        # find orientation of wire
        x, y, z = self.getXYZ()

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

        xColTemp = Wire.quad(1 + m ** 2, 2 * (A * m + self.start[x]), A ** 2 + self.start[x] ** 2 - self.r ** 2)

        withinSegment = False
        if len(xColTemp) == 0:
            return None
        if len(xColTemp) == 1:
            return None
        else:
            xColTemp = [xc for xc in xColTemp if Plate.inInterval(xc, (prevPos[x], pos[x]))]
            if len(xColTemp) != 0:
                withinSegment = True
                xCol = min(xColTemp, key=lambda p: abs(p - prevPos[0]))

        if not withinSegment:
            return None

        # find y coordinate of collision (yCol)
        yCol = m * (xCol - prevPos[x]) + prevPos[y]

        # find z coordinate of collision (zCol) and check if it's within the length vector of the cylinder
        zCol = mz * (xCol - prevPos[x]) + prevPos[z]
        withinHeight = min(self.start[z], self.end[z]) <= zCol <= max(self.start[z], self.end[z])
        

        if not withinHeight:
            return None

        colPos = np.zeros(3)
        colPos[z] = zCol
        colPos[x] = xCol
        colPos[y] = yCol

        # find normal vector to the cylinder to calculate the new velocity
        normal = np.zeros(3)
        normal[x] = xCol - self.start[x]
        normal[y] = yCol - self.start[y]
        normal[z] = 0
        normal /= self.r

        #print(vel, np.sum(normal ** 2 ))

        newVel = dampeningFactor * (vel - 2 * normal * np.dot(normal, vel))
        return colPos, newVel
    
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


    
