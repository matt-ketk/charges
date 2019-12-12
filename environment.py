import pygame
from pygame import gfxdraw
from pygame.locals import * 

import numpy as np

class Environment:
    def __init__(self, surface, surfaceSize):
        self.surface = surface 
        self.surfaceSize = surfaceSize 

        self.xRotate = 0
        self.yRotate = 0

        self.zoom = 10
        self.offset = np.array([surfaceSize[0] // 2, surfaceSize[0] // 2])  # TODO: I'm not sure if the offset works exactly as I think it does

        self.perspectiveMatrix = self.createPerspective(self.xRotate, self.yRotate)
    
    def createPerspective(self, x, y):
        self.xRotate = x
        self.yRotate = y

        cx = np.cos(self.xRotate)
        sx = np.sin(self.xRotate)

        cy = np.cos(self.yRotate)
        sy = np.sin(self.yRotate)

        yTransform = np.array([[cx , 0, sx, 0],
                               [0  , 1, 0 , 0],
                               [-sx, 0, cx, 0],
                               [0  , 0, 0,  1]])
        xTransform = np.array([[1, 0 , 0 , 0],
                               [0, cy,-sy, 0],
                               [0, sy, cy, 0],
                               [0, 0 , 0 , 1]])
        
        transform = np.dot(yTransform, xTransform)
        
        return transform

    def changePerspective(self, dx, dy):
        self.xRotate += dx
        self.yRotate += dy
    
        self.perspectiveMatrix = self.createPerspective(self.xRotate, self.yRotate)

    def transformPoint(self, position):
        # Transform the position with the perspective matrix first
        # position -> 1 x 3, perspectiveMatrix -> 3 x3
        if len(position.shape) == 1:
            position = position[np.newaxis, :]

        position = np.insert(position, 3, 1, axis = 1)     # Add the one to the end of each row
        
        position_ = np.matmul(position, self.perspectiveMatrix)

        # Offset, zoom and cast 
        position_ = (position_[:, :2] * self.zoom + self.offset).astype(int) 
        #print(position_)

        return position_
    
    def drawParticle(self, position, color, radius = 1):
        """Draws a particle on screen, accounting for perspective, panning and zoom.

        position: must be numpy array, the x, y and z coordinates of the particle
        color: a 3-tuple of RGB values
        radius: an integer
        """

        position_ = self.transformPoint(position)[0]

        # We only display the projection of the particle onto the XY plane
        #print (position_[0], position_[1], radius)
        if radius == 0:
            self.surface.set_at((position_[0], position_[1]), color)
        else:
            pygame.gfxdraw.aacircle(self.surface, position_[0], position_[1], radius, color)
            
    
    def drawMesh(self, polygons, color, filled = False):
        """Draws an arbitrary mesh, this will probably be used for drawing cubes/rectangular prisms/plates"""

        for p in polygons:
            # Each polygon has some number of points, if it's just 2 then draw a line
            if len(p) == 2:
                  p1, p2 = self.transformPoint(np.array(p))
                  
                  pygame.gfxdraw.line(self.surface, p1[0], p1[1], p2[0], p2[1], color)
                
            else:
                  p_ = self.transformPoint(np.array(p))
                  
                  if filled:
                    pygame.gfxdraw.filled_polygon(self.surface, p_, color)

                  else:
                    pygame.gfxdraw.polygon(self.surface, p_, color)
    
    def drawObject(self, obj, color, thickness = 1):
        obj.draw(self, color, thickness)
     
    def drawMeshRaw(self, polygons, color):
        for p in polygons:
            # Each polygon has some number of points, if it's just 2 then draw a line
            if len(p) == 2:
                  p1, p2 = p
                  
                  pygame.gfxdraw.line(self.surface, p1[0], p1[1], p2[0], p2[1], color)
                
            else:
                  p_ = p
                  
                  if filled:
                    pygame.gfxdraw.filled_polygon(self.surface, p_, color)

                  else:
                    pygame.gfxdraw.polygon(self.surface, p_, color)
     
        
        
    

    
    
