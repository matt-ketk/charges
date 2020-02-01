import pygame
from pygame.locals import *

import numpy as np

from environment import Environment
from charge import Charge
from p2p import forces, deltaPosition, deltaVelocity
from constants import Constants

from wire import Wire
# from plate import Plate
from latticeIon import LatticeIon

import matplotlib.pyplot as plt
import pylab

import debugtools as dbt

def main():
    coords = []
    charges = []
    masses = []
    stationary = []

    coordsAll = []
    chargesAll = []
    massesAll = []
    stationaryAll = []

    wireLength = 1E-7
    wireRadius = .125E-8

    wire0 = Wire(np.array([-wireLength / 2, 0]), np.array([wireLength, 0]), wireRadius)

    latticeIons = []
    # protons
    xNum = 25
    yNum = 3
    for i in range(0, xNum):
        for j in range(0, yNum):
            latticeIons.append(LatticeIon(np.array([-12E-10 * i + 10E-9, 6E-10 * j - 6E-10]), Constants.COPPER_ION_RADIUS, Constants.E, Constants.COPPER_MASS))

    for ion in latticeIons:
        coords.append(ion.center)
        charges.append(ion.charge)
        masses.append(ion.mass)
        stationary.append([1])
    

    # electrons
    for i in range(0, xNum):
        for j in range(0, yNum):
            coords.append([-12E-10 * i + 10E-9 + np.random.uniform(-1E-9, 1E-9), 6E-10 * j - 4E-10])
            charges.append(-Constants.E)
            masses.append(Constants.MASS_ELECTRON)
            stationary.append([0])


    coords, prevCoords, charges, masses, stationary = np.array(coords), np.array(coords), np.array(charges), np.array(masses), np.array(stationary)
    n = charges.size
    vel = np.zeros((n, 2))

    dt = 5E-12

    currents = []

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
    iterNum = 0
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

                f = forces(coords, charges) + np.array([5E-19, 0]).repeat(len(charges)).reshape((2, len(charges)))
                prevCoords = coords
                prevVel = vel
                vel = vel + deltaVelocity(dt, f, masses).T * (1 - stationary)
                coords = coords + deltaPosition(dt, f, vel.T, masses).T * (1 - stationary)

                if iterNum % 10 == 0:
                    currents.append(0)
                for k in range(len(latticeIons), len(charges)):
                    collision = wire0.checkCollision(prevCoords[k], coords[k], vel[k])
                    # if prevVel[k][0] < 0 and vel[k][0] > 0:
                    #     currents[-1] += 1
                    # if prevVel[k][0] > 0 and vel[k][0] < 0:
                    #     currents[-1] -= 1
                    currents[-1] += vel[k][0]
                    if collision:
                        wire0.checkCollision(prevCoords[k], coords[k], vel[k])
                        coords[k], vel[k] = collision
                        coords[k] += vel[k] * dt



                for ion in latticeIons:
                    for i in range(len(latticeIons), len(coords)):
                        coords[i], vel[i] = ion.checkCollision(prevCoords[i], coords[i], vel[i], dt, dampeningFactor=0.9)
                # collision detection for wire
                # print('prev coords:', prevCoords)
                # print('coords:', coords)
                # print('vel:', vel)
                # collisionPositions, vel = wire0.checkCollision(prevCoords, coords, vel)
                #coords = collisionPositions + 0.01 * dt * vel
                

                #print('done with loop')
                #done = True

                # collision detection for wire
                #print('prev coords:', prevCoords)
                #print('coords:', coords)
                #print('vel:', vel)
                
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
            for c in range(n):
                if charges[c] < 0:
                    pos = coords[c]
                    color = (100, 100, 255)
                    r = 1
                    env.drawParticle(pos, color, radius=1)

            wire0.draw(env, color = (255, 100, 100))

            for ion in latticeIons:
                ion.drawIon(env)

            #env.drawObject(w, color = (255, 100, 100))
            #env.drawCylinder(w.start, w.end, w.r, color = (255, 100, 100))
            wire0.draw(env, color=(50, 255, 150))
            # wire1.draw(env, color=(255,100,100))
            # wire2.draw(env, color=(255,100,100))
            pygame.display.flip()
            iterNum += 1
    currents = currents[:-1]
    print(currents)
    with open("currents.txt", "w") as f:
        f.write(",".join(str(current for current in currents)))
    iters = np.array(list(range(len(currents))))
    m, b = pylab.polyfit(iters, currents, 1)
    plt.plot(iters, currents)
    plt.plot(iters, m*iters + b)
    plt.show()

if __name__ == "__main__":
    main()