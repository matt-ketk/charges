import numpy as np

from wire import Wire
from plate import Plate

start = np.array([0, 0, 0])
lengthV = np.array([10, 0, 0])
w = Wire(start, lengthV, 0.5)

# print(Wire.quad(1, -1, -1))
prevPos = np.array([0, 1, 1])
pos = np.array([4, 1, 1])
# print(w.checkCollision(prevPos, pos))


p = Plate(np.ones(3), 3 * np.ones(3))
v = pos - prevPos
# print(p.checkCollision(prevPos, pos, v))
rCorners = np.array([[0, 0, 0], [2, 2, 2], [0, 4, 4], [-2, 2, 2]])
# print(rCorners)
# print(v)
print(Plate.rectangleCollision(prevPos, pos, rCorners, v))
