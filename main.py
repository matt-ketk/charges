import pygame
from pygame.locals import *

import numpy as np

from environment import Environment
from charge import Charge
from p2p import forces, deltaPosition, deltaVelocity
from constants import Constants
from wire import Wire
from plate import Plate
from latticeIon import LatticeIon

import debugtools as dbt

def main():
    coords = []
    charges = []
    masses = []
    stationary = []

    wireLength = 3E-8
    wireRadius = 3E-8
    
    w = Wire (np.array([0,0,-wireLength / 2]), np.array([0,0,wireLength]), wireRadius)
    latticeIons = LatticeIon.generateLatticePoints(w)
    
    print ("protons:", len(latticeIons))
    '''
    p1 = Plate (np.array([-1E-7, -1E-7, 0.25E-7]), np.array([0E-7, 0E-7, 0.5E-7]))
    lattice1 = LatticeIon.generateLatticePoints(p1)
    print (len(lattice1))

    print ("\n")
    
    p2 = Plate(np.array([-1E-7, -1E-7, 1.5E-7]), np.array([0E-7, 0E-7, 1.75E-7]))
    lattice2 = LatticeIon.generateLatticePoints(p2)
    print (len(lattice2))

    print ("\n")

    p3 = Plate (np.array([-1E-7, -1E-7, 2.75E-7]), np.array([0E-7, 0E-7, 3E-7]))
    lattice3 = LatticeIon.generateLatticePoints(p3)
    print (len(lattice3))
    '''
    electrons = LatticeIon.generateLatticePoints(w, charge = -Constants.E, mass = Constants.MASS_ELECTRON, offset = np.array([.1E-7, .1E-7, .1E-7]))
    print ("electrons:", len(electrons))
    
    
    for ion in latticeIons:
        center, charge, mass, stat = ion.compile()
        coords.append(center)
        charges.append(charge)
        masses.append(mass)
        stationary.append([stat])
    
    for ion in electrons:
        center, charge, mass, stat = ion.compile()
        coords.append(center)
        charges.append(charge)
        masses.append(mass)
        stationary.append([0])
    '''
    for ion in lattice2:
        center, charge, mass, stat = ion.compile()
        coords.append(center)
        charges.append(charge)
        masses.append(mass)
        stationary.append([stat])

    for ion in lattice3:
        center, charge, mass, stat = ion.compile()
        coords.append(center)
        charges.append(-charge)
        masses.append(mass)
        stationary.append([stat])
    
    
    # protons
    for i in range (2):
        for j in range (2):
            for k in range (5):
                coords.append([-9E-10 * i, -9E-10 * j, k * 9E-10 - 4E-10])
                charges.append(Constants.E)
                masses.append(Constants.COPPER_MASS)
                stationary.append([1])
    

    # electrons
    for i in range (1):
        coords.append(np.random.uniform(-.5E-7, -.5E-7, 3)  + np.array([0,0,np.random.uniform(2.75E-7, 2.75E-7)]))
        charges.append(-Constants.E)
        masses.append(Constants.MASS_ELECTRON)
        stationary.append([0])
    '''
    #print (coords)
    
    coords, prevCoords, charges, masses, stationary = np.array(coords), np.array(coords), np.array(charges), np.array(masses), np.array(stationary)
    n = charges.size
    print (n)
    vel = np.array([-1E-8, -1E-8,-1E-8]) #np.zeros((n, 3)) 

    dt = 5E-16

    pygame.init()
    screenSize = (800, 800)
    screen = pygame.display.set_mode(screenSize)
    env = Environment(screen, screenSize)
    env.changePerspective(0.8, 0.8)
    env.zoom = 15E8

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
            elif event.type == pygame.MOUSEMOTION:
                position = np.array(pygame.mouse.get_pos())
                    
                # Get rel and shift perspective
                rel = event.rel

                if events['mouseDown']:
                    dy, dx = rel[1] * 3.1415/180 / 5, rel[0] * 3.1415/180 / 5
                    env.changePerspective(dx, -dy)
            elif event.type == pygame.KEYDOWN and event.key == pygame.K_SPACE:
                paused = not paused

        if not paused:
            f = forces(coords, charges)
            prevCoords = coords
            vel = vel + deltaVelocity(dt, f, masses).T * (1 - stationary)
            coords = coords + deltaPosition(dt, f, vel.T, masses).T * (1  -stationary)

            # collision detection for wire  
            for i in range (n):
                if stationary[i][0]:
                    continue
                result = w.checkCollision(prevCoords[i], coords[i], vel[i])
                for ion in latticeIons:
                    result = ion.checkCollision(prevCoords[i], coords[i], vel[i])
                    if result:
                        break
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
            r = 11
            if charges[c] < 0:
                color = (100, 100, 255)
                r = 1
            
            env.drawParticle(pos, color, radius = r)

        #env.drawObject(w, color = (255, 100, 100))
        #env.drawCylinder(w.start, w.end, w.r, color = (255, 100, 100))
        w.draw(env, color = (255, 100, 100))

        pygame.display.flip()                

if __name__ == "__main__":
    main()