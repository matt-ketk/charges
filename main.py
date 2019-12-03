import pygame
from pygame.locals import *

import numpy as np

from environment import Environment
from charge import Charge
from p2p import forces, deltaPosition, deltaVelocity

import debugtools as dbt

def main():
  lCharges = [Charge(np.random.choice([-1, 1]) * 1E-8, np.random.random(3))for i in range (100)]

  n = 100

  coords = np.random.random((n, 3))
  charges = np.random.choice([-1, 1], n) * 1E-8
  masses = np.ones(n) * 1

  vel = np.zeros((n, 3))

  dt = 0.1

  f = forces(coords, charges)
  print(masses)
  vel = vel + deltaVelocity(dt, f, masses)
  coords = coords + deltaPosition(dt, f, vel, masses)

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

        for c in range (n):
          pos = coords[c]
          color = (255, 100, 100)
          if charges[c] < 0:
                color = (100, 100, 255)
          
          env.drawParticle(pos, color, radius = 5)

        pygame.display.flip()                

if __name__ == "__main__":
  main()