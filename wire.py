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
        if np.array_equal(pos, prevPos):
            return None
        # finds collision points for one wall of the wire, return None if not
        wallACollisions = self.lineIntersection(
            self.corners[1],
            self.corners[2],
            prevPos,
            pos
        )
        # ...and the other
        wallBCollisions = self.lineIntersection(
            self.corners[0],
            self.corners[3],
            prevPos,
            pos
        )
        # get the distances of collisions to each wall for comparison
        print(wallACollisions)
        colDistanceA = np.linalg.norm(
            wallACollisions,
            axis=1
        )
        
        colDistanceB = np.linalg.norm(
            wallBCollisions,
            axis=1
        )
        collisionDistances = np.hstack((
            colDistanceA,
            colDistanceB
        ))

        wallCollisions = np.dstack((
            wallACollisions,
            wallBCollisions
        ))
        collidedIndex = np.argmin(collisionDistances, axis=1)
        collisionPositions = wallCollisions[collidedIndex]

        vectorWallA = self.corners[2] - self.corners[1]
        vectorWallB = self.corners[3] - self.corners[0]

        normalWall = np.array([
            [1, -vectorWallA[0] / vectorWallA[1]],
            [1, -vectorWallB[0] / vectorWallB[1]]
        ])

        return collisionPositions, dampeningFactor * self.reflect(vel, normalWall[collidedIndex])         

    def reflect(self, v0, normal):
        return v0 - 2 * np.dot(v0, normal) * normal

    def lineIntersection(self, end0, end1, prevPos, pos):
        x = np.array([
            np.ones(pos.shape[0]) * end0[0],
            np.ones_like(pos.shape[0]) * end1[0],
            prevPos[:,0],
            pos[:,0]
        ])

        y = np.array([
            np.ones(pos.shape[0]) * end0[1],
            np.ones(pos.shape[0]) * end1[1],
            prevPos[:,1],
            pos[:,1]
        ])

        xOut = ((x[0] * y[1] - y[0] * x[1]) * (x[2] - x[3]) - (x[0] - x[1]) * (x[2] * y[3] - y[2] * y[3])) / ((x[0] - x[1]) * (y[2] - y[3]) - (y[0] - y[1]) * (x[2] - x[3]))        
        yOut = ((x[0] * y[1] - y[0] * x[1]) * (y[2] - y[3]) - (y[0] - y[1]) * (x[2] * y[3] - y[2] * x[3])) / ((x[0] - x[1]) * (y[2] - y[3]) - (y[0] - y[1]) * (x[2] - x[3]))

        return np.vstack((xOut, yOut))
    '''
    def lineIntersection(self, end0, end1, prevPos, pos):
        # print(end0, end1)
        l0 = np.vstack((end0, end1))
        l1 = np.dstack((prevPos, pos))
        diff = np.array([
            [l0[0][0] - l0[1][0], l1[0][0] - l1[1][0]],
            [l0[0][1] - l0[1][1], l1[0][1] - l1[1][1]]
        ])
        xDiff = np.array([
            np.ones_like(pos) * (l0[0,0] - l1[1,0]),
            l1[:,0,0] - l1[:,1,0]
        ])
        yDiff = np.array([
            np.ones_like(pos) * (l0[0,1] - l0[1,1]),
            l1[:,0,1] - l1[:,1,1]
        ])

        div = np.linalg.det(diff, axis=1)
        if div == 0: # when div == 0, there's no intersection thus no collision
            return None
        d = np.array([np.linalg.det(l0), np.linalg.det(l1)])
        x = np.linalg.det(np.vstack((d, diff[0])))
        y = np.linalg.det(np.vstack((d, diff[1])))

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
    '''

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