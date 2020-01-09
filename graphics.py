import pygame
from pygame.locals import * 

import numpy as np

class Environment:
    def __init__(self, surface):
        self.surface = surface 

        self.xRotate = 0
        self.yRotate = 0

        self.zoom = 10
        self.offset = np.array([0, 0])  # TODO: I'm not sure if the offset works exactly as I think it does

        self.perspectiveMatrix = self.createPerspective(self.xRotate, self.yRotate)
    
    def createPerspective(self, x, y):
        self.xRotate = x
        self.yRotate = y

        cx = np.cos(y_rotate)
        sx = np.sin(y_rotate)

        cy = np.cos(y_rotate)
        sy = np.sin(y_rotate)

        yTransform = np.array([[cx , 0, sx, 0],
                               [0  , 1, 0 , 0],
                               [-sx, 0, cx, 0],
                               [0  , 0, 0,  1]])
        xTransform = np.array([[1, 0 , 0 , 0],
                               [0, cy,-sy, 0],
                               [0, sy, cy, 0],
                               [0, 0 , 0 , 1]])
        
        #transform = np.matmul(yTransform, xTransform)
        transform = xTransform

        return transform

    def changePerspective(self, dx, dy):
        self.xRotate += dx
        self.yRotate += dy
    
        self.perspectiveMatrix = self.createPerspective(self.xRotate, self.yRotate)

    def transformPoint(self, position):
        # Transform the position with the perspective matrix first
        # position -> 1 x 3, perspectiveMatrix -> 3 x3
        position_ = np.matmul(position, self.perspectiveMatrix)

        # Offset, zoom and cast 
        position_ = (position_ * self.zoom + self.offset).astype(int)

        return position_
    
    def drawParticle(self, position, color, radius = 4):
        """Draws a particle on screen, accounting for perspective, panning and zoom.

        position: must be numpy array, the x, y and z coordinates of the particle
        color: a 3-tuple of RGB values
        radius: an integer
        """

        position_ = self.transformPoint(position)

        # We only display the projection of the particle onto the XY plane
        pygame.gfxdraw.aacircle(self.surface, position_[0], position_[1], radius, color)
    
    def drawMesh(self, polygons, color):
        """Draws an arbitrary mesh, this will probably be used for drawing cubes/rectangular prisms/plates"""

        for p in polygons:
            # Each polygon has some number of points, if it's just 2 then draw a line
            if len(p) == 2:
                  p1, p2 = self.transformPoint(p)
                  pygame.gfxdraw.aaline(self.surface, p1[0], p1[1], p2[0], p2[1], color)
                
            else:
                  p_ = self.transformPoint(p)
                  pygame.gfxdraw.polygon(self.surface, p_, color)
     
        
        
    

    
    
