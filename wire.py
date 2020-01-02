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
        
        environment.drawMesh([points], color)