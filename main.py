import pygame
from pygame.locals import *

import numpy as np

from environment import Environment
from charge import Charge
from p2p import forces, deltaPosition, deltaVelocity
from constants import Constants
from wire import Wire

import debugtools as dbt

def main():
  np.random.seed(10)

  coords = []
  charges = []
  masses = []
  stationary = []

  w = Wire (np.array([0,0,-25]), np.array([0,0,50]), 0.5)
  
  # protons
  for i in range (2):
        for j in range (2):
            for k in range (25):
                coords.append([-0.3 + i * 0.6,-0.3 + j * 0.6, k * 0.08])
                charges.append(1 * 1E-8)
                masses.append(1E-3)
                stationary.append([1])
  
  # electrons
  for i in range (20):
        coords.append(np.random.uniform(-0.2, 0.2, 3) * 2 + np.array([0,0,0.6]))
        charges.append(-1 * 1E-8)
        masses.append(1E-3)
        stationary.append([0])

  coords, prevCoords, charges, masses, stationary = np.array(coords), np.array(coords), np.array(charges), np.array(masses), np.array(stationary)
  n = charges.size
  vel = np.zeros((n, 3)) 
  vel[0] = np.array([0.01, 0.02, 0]).astype(float)

  vel2 = np.zeros((n, 3))
  
  dt = 0.1

  pygame.init()
  screenSize = (800, 800)
  screen = pygame.display.set_mode(screenSize)
  env = Environment(screen, screenSize)
  env.changePerspective(0, 0)
  env.zoom = 100

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

        screen.fill((255, 255, 255))

        f = forces(coords, charges)
        
        prevCoords = coords
        vel = vel + deltaVelocity(dt, f, masses).T * (1 - stationary)
        coords = coords + deltaPosition(dt, f, vel.T, masses).T * (1  -stationary)

        for i in range (n):
          result = w.checkCollision(prevCoords[i], coords[i], vel[i])
          if (coords[i][0] ** 2 +coords[i][1]**2 > w.r**2) and not result:
              w.checkCollision(prevCoords[i], coords[i], vel[i])
          if result:
            colPos, newVel = result
            vel[i] = newVel

            #prevCoords[i] = coords[i]
            coords[i] = colPos + 0.005 * newVel

          collision = bool(result)

          
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
               if collision:
                 color = (50, 200, 50)
               r = 1
          
          
          env.drawParticle(pos, color, radius = r)

          #env.drawParticle(pos2, color2, radius = 6)

        env.drawCylinder(np.array([0,0,-25]), np.array([0,0,25]), 0.5, color = (255, 100, 100))

        pygame.display.flip()                

if __name__ == "__main__":
  main()