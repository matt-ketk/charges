import pygame
from pygame.locals import *

import numpy as np

from environment import Environment
from charge import Charge
from p2p import forces, deltaPosition, deltaVelocity
from constants import Constants

import debugtools as dbt

def main():
  lCharges = [Charge(np.random.choice([-1, 1]) * 1E-8, np.random.random(3))for i in range (100)]

  n = 120
  coords = []
  charges = []
  masses = []
  stationary = []
  for i in range (10):
        for j in range (10):
              
              coords.append([-1 + i * 0.2, -1 + j * 0.2, 0])
              charges.append(1 * 1E-8)
              masses.append(1E-3)
              stationary.append([1])
    
  for i in range (20):
        coords.append(np.random.uniform(-1, 1, 3) * 2)
        charges.append(-1 * 1E-8)
        masses.append(1E-3)
        stationary.append([0])

  coords, charges, masses, stationary = np.array(coords), np.array(charges), np.array(masses), np.array(stationary)

  vel = np.zeros((n, 3)) 
  vel[0] = np.array([0, 0, -0.1]).astype(float)

  vel2 = np.zeros((n, 3))
  
  dt = 0.1

  pygame.init()
  screen = pygame.display.set_mode((1000, 1000))

  env = Environment(screen, (1000, 1000))
  env2 = Environment(screen, (1000, 1000))
  env.offset = np.array([400, 500])
  env2.offset = np.array([600, 500])
  env.changePerspective(-0.05, 0)
  env2.changePerspective(0.05, 0)
  env2.zoom = 40
  env.zoom = 40

  # Relevant variables
  events = {
    "mouseDown" : False,
    "prevMouse" : None
  }
  screen.fill((255, 255, 255))

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
                    env2.changePerspective(dx, -dy)
          
        screen.fill((255, 255, 255))

        f = forces(coords, charges)
        
        vel = vel + deltaVelocity(dt, f, masses).T * (1 - stationary)
        coords = coords + deltaPosition(dt, f, vel.T, masses).T * (1  -stationary)

        # also a magnetic field in the positive z-axis, [0, 0, 1]
        

        """
        ts = 0.1
        accs =[ ]
        for c in range (n):
          acc = 0
          for c2 in range (n):
            vec = coords2[c] - coords2[c2]
            env.drawMesh([[coords2[c2], coords2[c]]], (100, 250, 100))
            dist = np.sqrt(np.sum(vec ** 2))

            if dist > 0.0001:
              pf = 1 * Constants.K * charges[c] * charges[c2] / (dist ** 3 + 1E-4) * vec
            
              acc += pf / masses[c]

            vel2[c] += acc * ts
            accs.append(acc)
    
        for c in range (n):
          #pass
          coords2[c] += vel2[c] * ts + 0.5 * ts ** 2 * accs[c]
        """
            
        for c in range (n):
          pos = coords[c]
          color = (255, 100, 100)
          r = 3
          if charges[c] < 0:
               color = (100, 100, 255)
               r = 1

               if pos[-1] > 0:
                     color = (200, 200, 255)
          
          env.drawParticle(pos, color, radius = r)
          env2.drawParticle(pos, color, radius = r)

          pygame.draw.rect(screen, (255, 0, 0), (400 - 10, 650, 20, 15), 2)
          pygame.draw.rect(screen, (255, 0, 0), (400 - 10 + 200, 650, 20, 15))
          #env.drawParticle(pos2, color2, radius = 6)

        pygame.display.flip()                

if __name__ == "__main__":
  main()