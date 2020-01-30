import numpy as np

from constants import Constants
from conductor import Conductor
from plate import Plate

class Wire(Conductor):
    def __init__(self, start, lengthV, r, resistivity=Constants.COPPER_RESISTIVITY):
        '''
        start: length-2 np array representing start of the wire
        lengthV: length-2 np array representing length (and direction) of wire
        r: radius of wire
        '''

        super(Wire, self).__init__(resistivity=resistivity)

        self.start = start
        self.lengthV = lengthV
        self.r = r
        self.end = start + lengthV
        
        length = np.linalg.norm(lengthV)
        horizontal_end = np.array([length, 0])
        angle = np.arccos(lengthV[0] / length)
        print('angle {}'.format(angle))
        print('start: {} end: {}'.format(self.start, self.end))
        pre_rotated = np.array([
            -np.array([0, self.r]),
            np.array([0, self.r]),
            horizontal_end + np.array([0, self.r]),
            horizontal_end - np.array([0, self.r]) 
        ])

        self.corners = pre_rotated @ self.rotMatrix(angle)
        print(pre_rotated)
        print(self.rotMatrix(angle))
        print(self.corners)
    # vectorized version, only works for axis-aligned wires

    def rotMatrix(self, theta):
        c = np.cos(theta)
        s = np.sin(theta)
        return np.array([
            [c, -s],
            [s, c]
        ])

    def checkCollision(self, prevPos, pos, vel, dampeningFactor=1):
        collisionPositionsA, collisionT_A = self.lineIntersection(
            prevPos,
            pos,
            np.array(self.corners[1]),
            np.array(self.corners[2])
        )

        collisionPositionsB, collisionT_B = self.lineIntersection(
            prevPos,
            pos,
            np.array(self.corners[0]),
            np.array(self.corners[3])
        )

        colPos = None

        validA, validB = not (collisionT_A < 0 or collisionT_A > 1), not (collisionT_B < 0 or collisionT_B > 1)
        
        print(collisionT_A, collisionT_B)

        if not(validA) and not(validB):
            return False

        if (validA and validB):
            if validA < validB:
                colPos = collisionPositionsA
            else:
                colPos = collisionPositionsB
        
        elif (validA):
            colPos = collisionPositionsA
        elif validB:
            colPos = collisionPositionsB

        vectorWallA = self.corners[2] - self.corners[1]
      
        normal = np.array([1, -vectorWallA[0] / vectorWallA[1]])
        normal /= np.linalg.norm(normal)

        newVel = self.reflect(vel, normal) 

        return colPos, newVel

    def reflect(self, v0, normal):
        return v0 - 2 * np.dot(v0, normal) * normal

    def lineIntersection(self, s1, e1, s2, e2):
        # represent line 1 as y = mx + b
        v1 = (e1 - s1).astype(float)
        v2 = (e2 - s2).astype(float)

        a = s2[1] - s1[1] - (v1[1] * s2[0] - v1[1] * s1[0]) / v1[0]
        b = (v1[1] * v2[0] - v2[1] * v1[0]) / (v1[0])

        t = (a / b)
        s = (s2[0] + t * v2[0] - s1[0]) / v1[0]

        print(np.array([s2 +v2 * t]).shape)
        return s2 + v2 * t, s

    def draw(self, environment, color=(255, 255, 0), thickness=1):
        if self.lengthV[0] == 0:
            points = [
                self.start + np.array([self.r, 0]),
                self.start - np.array([self.r, 0]),
                self.end - np.array([self.r, 0]),
                self.end + np.array([self.r, 0])
            ]
            
        if self.lengthV[1] == 0:
            points = [
                self.start + np.array([0, self.r]),
                self.start - np.array([0, self.r]),
                self.end - np.array([0, self.r]),
                self.end + np.array([0, self.r])
            ]
        environment.drawMesh([self.corners], color)