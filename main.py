import pygame
from pygame.locals import *

import numpy as np

from environment import Environment

def main():
  pygame.init()
  screen = pygame.display.set_mode((1000, 1000))

  env = Environment(screen, (1000, 1000))

  # Relevant variables
  events = {
    "mouseDown" : False,
    "prevMouse" : None
  }

  # Houeskeeping for the event loop
  done = False
  while not done:
        for event in pygame.event.get():
              if event.type == pygame.QUIT:
                  pygame.quit()
                  done = True
        
              elif event.type == pygame.MOUSEBUTTONDOWN:
                  events['mouseDown'] = True
              elif event.type == pygame.MOUSEBUTTONUP:
                  events['mouseDown'] = False
              elif event.type == pygame.MOUSEMOTION:
                  position = np.array(pygame.mouse.get_pos())
                  
                  # Get rel and shift perspective
                  rel = event.rel

                  if events['mouseDown']:
                    dy, dx = rel[1] * 3.1415/180 / 5, rel[0] * 3.1415/180 / 5
                    env.changePerspective(dx, -dy)
          
        screen.fill((255, 255, 255))
      
        env.drawCylinder(np.array([0, 0, 0]), np.array([0, 0, 20]), 5, (255, 100, 100))
        pygame.display.flip()                

if __name__ == "__main__":
  main()