import pygame
from pygame.locals import *

import numpy as np

from environment import Environment
from charge import Charge
from p2p import forces, deltaPosition, deltaVelocity
from constants import Constants

from wire import Wire
# from plate import Plate
# from latticeIon import LatticeIon

import debugtools as dbt

def main():
    coords = []
    charges = []
    masses = []
    stationary = []

    wireLength = 1E-8
    wireRadius = .0625E-8

    wire0 = Wire(np.array([0,-wireLength/2]), np.array([0,wireLength]), wireRadius)

    # protons
    for i in range (2):
        for j in range (2):
            coords.append([-9E-10 * i, -9E-10 * j])
            charges.append(Constants.E)
            masses.append(Constants.COPPER_MASS)
            stationary.append([1])
    

    # electrons
    for i in range(4):
        coords.append(np.random.uniform(-1E-9, 1E-9, 2))
        charges.append(-Constants.E)
        masses.append(Constants.MASS_ELECTRON)
        stationary.append([0])

    #print (coords)
    
    coords, prevCoords, charges, masses, stationary = np.array(coords), np.array(coords), np.array(charges), np.array(masses), np.array(stationary)
    n = charges.size
    print(n)
    vel = np.zeros((n, 2))

    dt = 1E-13

    pygame.init()
    screenSize = (800, 800)
    screen = pygame.display.set_mode(screenSize)
    env = Environment(screen, screenSize)
    env.zoom = 1E11

    # Relevant variables
    events = {
        "mouseDown" : False,
        "prevMouse" : None
    }
    screen.fill((255, 255, 255))

    # Houeskeeping for the event loop
    done = False
    paused = False

    while not done:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                done = True
            elif event.type == pygame.MOUSEBUTTONDOWN:
                events['mouseDown'] = True
                if event.button == 4:
                    env.zoom *= 0.8
                if event.button == 5:
                    env.zoom *= 1.2
            elif event.type == pygame.MOUSEBUTTONUP:
                events['mouseDown'] = False

            elif event.type == pygame.KEYDOWN and event.key == pygame.K_SPACE:
                paused = not paused
        if not done:
            if not paused:
                f = forces(coords, charges)
                prevCoords = coords
                vel = vel + deltaVelocity(dt, f, masses).T * (1 - stationary)
                coords = coords + deltaPosition(dt, f, vel.T, masses).T * (1  -stationary)

                # collision detection for wire
                for i in range (n):
                    if stationary[i][0]:
                        continue
                    result = wire0.checkCollision(prevCoords[i], coords[i], vel[i])
                    # for ion in latticeIons:
                    #     result = ion.checkCollision(prevCoords[i], coords[i], vel[i])
                    #     if result:
                    #         break
                    if result:
                        colPos, newVel = result
                        vel[i] = newVel
                
                        coords[i] = colPos + 0.01 * dt * newVel
                
                    collision = bool(result)
            '''
            # collision detection for plate
            for i in range (n):
                result = p1.checkCollision(prevCoords[i], coords[i], vel[i])
                if result:
                    colPos, newVel = result
                    vel[i] = newVel
    
                    coords[i] = colPos + 0.01 * dt * newVel
    
                collision = bool(result)
            
            # also a magnetic field in the positive z-axis, [0, 0, 1]
            '''
            screen.fill((255, 255, 255))
            for c in range (n):
                pos = coords[c]
                color = (255, 100, 100)
                r = 8
                if charges[c] < 0:
                    color = (100, 100, 255)
                    r = 1

                env.drawParticle(pos, color, radius = r)

            #env.drawObject(w, color = (255, 100, 100))
            #env.drawCylinder(w.start, w.end, w.r, color = (255, 100, 100))
            wire0.draw(env, color=(255, 100, 100))

            pygame.display.flip()

if __name__ == "__main__":
    main()