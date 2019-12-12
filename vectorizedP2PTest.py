import numpy as np

import pygame 
from pygame.locals import *

pygame.init()

screen = pygame.display.set_mode((1000, 1000))

charges = np.array([[-1, 1, 1, 1, -1, -1]]).astype(float)
masses = np.array([[1, 1, 1, 1, 1, 1]]).astype(float).T
positions = np.random.random([charges.shape[1], 3]) * 2 - 1
velocities = np.zeros([charges.shape[1], 3])

x, y, z = positions[:, 0:1], positions[:, 1:2], positions[:, 2:3]
vx, vy, vz = velocities[:, 0:1], velocities[:, 1:2], velocities[:, 2:3]

print(x.shape)

def step(charges, masses, x, y, z, vx, vy, vz, t = 0.01):
    xt = x.T
    yt = y.T
    zt = z.T

    Dx = np.matmul(x, xt) + x**2 + xt ** 2
    Dy = np.matmul(y, xt) + y**2 + yt ** 2
    Dz = np.matmul(z, xt) + z**2 + zt ** 2

    D = np.sqrt(Dx + Dy + Dz)

    Q = np.matmul(charges.T, charges)

    Fx = np.sum(
            np.multiply(np.multiply(1 / D**3, Q), Dx),
            axis = -1, keepdims = True)
    Fy = np.sum(
            np.multiply(np.multiply(1 / D**3, Q), Dy),
            axis = -1, keepdims = True)
    Fz = np.sum(
            np.multiply(np.multiply(1 / D**3, Q), Dz),
            axis = -1, keepdims = True)
    
    Ax, Ay, Az = Fx/masses, Fy/masses, Fz/masses
    print(Ax.shape, Fx.shape, masses.shape)
    print(vx.shape)

    vx += Ax * t
    vy += Ay * t
    vz += Az * t

step(charges, masses, x, y, z, vx, vy, vz)

print(vx)


