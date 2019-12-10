import pygame
from pygame.locals import *

import numpy as np

from environment import Environment
from charge import Charge
from p2p import forces, deltaPosition, deltaVelocity
from constants import Constants
from wire import Wire
from latticeIon import LatticeIon

import debugtools as dbt

def main():
  #np.random.seed(10)

  coords = []
  charges = []
  masses = []
  stationary = []
  w = Wire (np.array([0,0,-20E-10]), np.array([0,0,4E-9]), 4E-10)
  
  lattice = LatticeIon.generateLatticePoints(w)

  for ion in lattice:
    center, charge, mass, stat = ion.compile()
    coords.append(center)
    charges.append(Constants.E)
    masses.append(Constants.COPPER_MASS)
    stationary.append([stat])
  
  
  
  # protons
  #for i in range (2):
  #      for j in range (2):
  #          for k in range (2):
  #              coords.append([-9E-10 * i, -9E-10 * j, k * 9E-10 - 4E-10])
  #              charges.append(Constants.E)
  #              masses.append(Constants.COPPER_MASS)
  #              stationary.append([1])
  
  
  # electrons
  for i in range (20):
        coords.append(np.random.uniform(-3E-10, 3E-10, 3)  + np.array([0,0,np.random.uniform(-5E-10, 5E-10)]))
        charges.append(-Constants.E)
        masses.append(Constants.MASS_ELECTRON)
        stationary.append([0])
  
  #print (coords)
  coords, prevCoords, charges, masses, stationary = np.array(coords), np.array(coords), np.array(charges), np.array(masses), np.array(stationary)
  
  print ("PROTONS")
  print (coords[:-1])
  print ("\n\n\n")
  print ("ELECTRONS")
  print (coords[-1:])
   
  n = charges.size
  print (charges)
  print (n)
  vel = np.zeros((n, 3)) 
  #vel = vel + np.random.uniform(-10E-10, 10E-10, 3)
  #vel = vel + 1E-4
  #vel[0] = np.array([0.01, 0.02, 0]).astype(float)

  vel2 = np.zeros((n, 3))
  
  dt = 1E-17

  pygame.init()
  screenSize = (800, 800)
  screen = pygame.display.set_mode(screenSize)
  env = Environment(screen, screenSize)
  env.changePerspective(0, 1.57)
  env.zoom = 10E10

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

        f = np.array([[0], [0], [-5e-10]]) + forces(coords, charges) 
   
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
            coords[i] = colPos + dt * newVel * 2
            print(newVel, colPos, coords[i])

          collision = bool(result)

         
        # also a magnetic field in the positive z-axis, [0, 0, 1]
        #print ("\n")

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
          
          
          env.drawParticle(pos, color, radius = r)

          #env.drawParticle(pos2, color2, radius = 6)

        env.drawCylinder(w.start, w.end, w.r, color = (255, 100, 100))

        pygame.display.flip()                

if __name__ == "__main__":
  main()