import numpy as np
from charge import Charge


class Nucleus:
  def __init__(self, center, radius, dampeningFactor):
    self.center = center
    self.radius = radius
    self.dampeningFactor = dampeningFactor

  def reflectVector(self, particle):
    n = self.collisionNormal(particle)
    v = particle.velocity
    return self.dampeningFactor * (v - 2 * n * np.dot(n, v))
  # returns a unit vector that is a surface normal at the point of collision
  def collisionNormal(self, particle):
    v = self.collisionPoint(particle) - self.center
    return v / np.linalg.norm(v)

  def collisionPoint(self, particle):
    unitV = particle.velocity / np.linalg.norm(particle.velocity)
    origin = particle.position
    distance = 0.0

    # this is where it gets tricky... so basically,

    # between a line and a sphere, either oneh, two, or no intersections are formed

    # put simply, the calculations look kind of like when solving a quadratic equation, and they are...
    # ad^2 + bd + c = 0

    d = origin - self.center # this is for readability

    # a = np.linalg.norm(unitV) ** 2 # a just becomes 1
    b = 2 * (np.dot(unitV, d))
    c = np.dot(d, d) - self.radius ** 2

    # in the case that the quadratic formula outputs two possible results, I CURRENTLY DO NOT KNOW HOW TO DETERMINE WHICH ONE IS WHICH. There is the entry point and there is the exit point (if you know what I mean) imma leave it here fo now...

    distance = -(b / 2) - np.sqrt((b / 2) ** 2 - c)
    # distance = -(b / 2) - np.sqrt(b ** 2 - c) # the alternative solution

    return (origin + distance * unitV)

  # essentially boolean output whether particle is gonna hit given its trajectory OF THAT TICK (as in its current position + velocity * dt)
  # has the same meat and bones of the method above, but now it outputs a boolean based on the quadratic discriminant formula, b**2 -4ac >= 0.
  
  def willCollide(self, particle):
    unitV = particle.velocity / np.linalg.norm(particle.velocity)

    d = particle.position - self.center 

    a = np.linalg.norm(unitV) ** 2
    b = 2 * (np.dot(unitV, d))
    c = np.dot(d, d) - self.radius ** 2

    return (b ** 2 - 4 * (a * c) >= 0)
      
def main():
  p = Charge (1, np.array([0,0,0]), np.array([1,1,1]))
  n = Nucleus(np.array([4,4,4]), 1, 0.5)

  print (n.collisionPoint(p))
  

if __name__ == "__main__":
    main()
