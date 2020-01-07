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
        if all(pos == prevPos):
            return None

        
    def intersection(l0, l1):
        xdiff = (l0[0][0] - l0[1][0], l1[0][0] - l1[1][0])
        ydiff = (l0[0][1] - l0[1][1], l1[0][1] - l1[1][1])

        def det(a, b):
            return a[0] * b[1] - a[1] * b[0]

        div = det(xdiff, ydiff)
        if div == 0:
            return None
        d = np.array([det(*l0), det(*l1)])
        x = det(d, xdiff) / div
        y = det(d, ydiff) / div
        return np.array([x, y])
                
    def checkCollision(self, prevPos, pos, vel, dampeningFactor=1):
        # if nothing's moving, don't bother
        if all(pos == prevPos):
            return None
        # if length is pointing in y-direction,
        
        if self.lengthV[0] == 0:
            wall0 = self.start[0] - self.r
            wall1 = self.start[0] + self.r
            if Plate.inInterval(wall0,(prevPos[0], pos[0])): 
                collisionY = (vel[1] / vel[0]) * (wall0 - prevPos[0]) + prevPos[1]
                newVel = np.array([-vel[0],vel[1]])
                return np.array([wall0, collisionY]), newVel 
            elif Plate.inInterval(wall1,(prevPos[0], pos[0])):
                collisionY = (vel[1] / vel[0]) * (wall1 - prevPos[0]) + prevPos[1]
                newVel = np.array([-vel[0],vel[1]])
                return np.array([wall1, collisionY]), newVel
        if self.lengthV[1] == 0:
            wall0 = self.start[1] - self.r
            wall1 = self.start[1] + self.r
            if Plate.inInterval(wall0,(prevPos[1], pos[1])): 
                collisionX = (vel[0] / vel[1]) * (wall0 - prevPos[1]) + prevPos[0]
                newVel = np.array([vel[0],-vel[1]])
                return np.array([collisionX, wall0]), newVel 
            elif Plate.inInterval(wall1,(prevPos[1], pos[1])):
                collisionX = (vel[0] / vel[1]) * (wall1 - prevPos[1]) + prevPos[0]
                newVel = np.array(vel[0],-vel[1]
                )
                return np.array([collisionX, wall1]), newVel
        return None 

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