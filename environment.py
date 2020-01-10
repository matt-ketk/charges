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

    def transformPoint(self, position):

        return (position * self.zoom + self.offset).astype(int)
    
    def drawParticle(self, position, color, radius=1):
        """Draws a particle on screen, accounting for perspective, panning and zoom.

        position: must be numpy array, the x, y and z coordinates of the particle
        color: a 3-tuple of RGB values
        radius: an integer
        """

        position_ = self.transformPoint(position)
        
        if not self.surface.get_rect().collidepoint(position_):
            return

        # We only display the projection of the particle onto the XY plane
        if radius == 0:
            self.surface.set_at((position_[0], position_[1]), color)
        else:
            pygame.gfxdraw.aacircle(self.surface, position_[0], position_[1], radius, color)
            
    
    def drawMesh(self, polygons, color, filled = False):
        """Draws an arbitrary mesh, this will probably be used for drawing rectangles/circles"""

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
    
    def drawObject(self, obj, color, thickness=1):
        obj.draw(self, color, thickness)
     
    def drawMeshRaw(self, polygons, color, filled=False):
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
