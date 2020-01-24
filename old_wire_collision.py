   def checkCollision(self, prevPos, pos, vel, dampeningFactor=1):
        print('pos', pos)
        print('prevpos', prevPos)
        if np.array_equal(pos, prevPos):
            return None
        # finds collision points for one wall of the wire, return None if not
        wallACollisions, colDistanceA = self.lineIntersection(
            prevPos,
            pos,
            np.array([self.corners[1]]),
            np.array([self.corners[2]])
        )
        # ...and the other
        wallBCollisions, colDistanceB = self.lineIntersection(
            prevPos,
            pos,
            np.array([self.corners[0]]),
            np.array([self.corners[3]])
        )
        # get the distances of collisions to each wall for comparison
        collisionDistances = np.hstack((
            colDistanceA,
            colDistanceB
        ))
        print('coldista:', colDistanceA)
        print('coldistb:', colDistanceB)
        wallCollisions = np.dstack((
            wallACollisions,
            wallBCollisions
        ))


        collidedT = np.where((collisionDistances > 0) & (collisionDistances < 1), collisionDistances, 2)        
        print('collided index', collidedT)
        print('wall collisions shape:', wallCollisions.shape)
        collisionPositions = np.choose(np.argmin(collidedT, axis=-1), np.transpose(wallCollisions, axes=(1,0,2)))
        
        '''
        colDistanceA = np.linalg.norm(
            wallACollisions,
            axis=1
        )
        
        colDistanceB = np.linalg.norm(
            wallBCollisions,
            axis=1
        )
        collisionDistances = np.vstack((
            colDistanceA,
            colDistanceB
        ))

        wallCollisions = np.dstack((
            wallACollisions,
            wallBCollisions
        ))
        collidedIndex = np.argmin(collisionDistances, axis=1)
        '''
        # collisionPositions = wallCollisions[collidedIndex]
        
        vectorWallA = self.corners[2] - self.corners[1]
        vectorWallB = self.corners[3] - self.corners[0]

        normalWall = np.array([
            [1, -vectorWallA[0] / vectorWallA[1]],
            [1, -vectorWallB[0] / vectorWallB[1]]
        ])
        print('collision positions', collisionPositions)
        return collisionPositions, dampeningFactor * self.reflect(vel, np.array([1,-vectorWallA[0] / vectorWallA[1]]))      

'''
    def lineIntersection(self, end0, end1, prevPos, pos):
        x = np.array([
            np.ones(pos.shape[0]) * end0[0],
            np.ones_like(pos.shape[0]) * end1[0],
            prevPos[:,0],
            pos[:,0]
        ])

        y = np.array([
            np.ones(pos.shape[0]) * end0[1],
            np.ones(pos.shape[0]) * end1[1],
            prevPos[:,1],
            pos[:,1]
        ])

        xOut = ((x[0] * y[1] - y[0] * x[1]) * (x[2] - x[3]) - (x[0] - x[1]) * (x[2] * y[3] - y[2] * y[3])) / ((x[0] - x[1]) * (y[2] - y[3]) - (y[0] - y[1]) * (x[2] - x[3]))        
        yOut = ((x[0] * y[1] - y[0] * x[1]) * (y[2] - y[3]) - (y[0] - y[1]) * (x[2] * y[3] - y[2] * x[3])) / ((x[0] - x[1]) * (y[2] - y[3]) - (y[0] - y[1]) * (x[2] - x[3]))


        return np.vstack((xOut, yOut))
    '''
    '''
    def lineIntersection(self, end0, end1, prevPos, pos):
        # print(end0, end1)
        l0 = np.vstack((end0, end1))
        l1 = np.dstack((prevPos, pos))
        diff = np.array([
            [l0[0][0] - l0[1][0], l1[0][0] - l1[1][0]],
            [l0[0][1] - l0[1][1], l1[0][1] - l1[1][1]]
        ])
        xDiff = np.array([
            np.ones_like(pos) * (l0[0,0] - l1[1,0]),
            l1[:,0,0] - l1[:,1,0]
        ])
        yDiff = np.array([
            np.ones_like(pos) * (l0[0,1] - l0[1,1]),
            l1[:,0,1] - l1[:,1,1]
        ])

        div = np.linalg.det(diff, axis=1)
        if div == 0: # when div == 0, there's no intersection thus no collision
            return None
        d = np.array([np.linalg.det(l0), np.linalg.det(l1)])
        x = np.linalg.det(np.vstack((d, diff[0])))
        y = np.linalg.det(np.vstack((d, diff[1])))

        return np.array([x, y])

    def checkCollision(self, prevPos, pos, vel, dampeningFactor=1):
        # if nothing's moving, don't bother
        if all(pos == prevPos):
            return None
        # if length is pointing in y-direction,
        
        if self.lengthV[0] == 0:
            wall0 = self.start[0] - self.r
            wall1 = self.start[0] + self.r
            if Plate.inInterval(wall0,(prevPos[0], pos[0])): 
                collisionY = (vel[1] / vel[0]) * (wall0 - prevPos[0]) + prevPos[1]
                newVel = np.array([-vel[0],vel[1]])
                return np.array([wall0, collisionY]), newVel 
            elif Plate.inInterval(wall1,(prevPos[0], pos[0])):
                collisionY = (vel[1] / vel[0]) * (wall1 - prevPos[0]) + prevPos[1]
                newVel = np.array([-vel[0],vel[1]])
                return np.array([wall1, collisionY]), newVel
        if self.lengthV[1] == 0:
            wall0 = self.start[1] - self.r
            wall1 = self.start[1] + self.r
            if Plate.inInterval(wall0,(prevPos[1], pos[1])): 
                collisionX = (vel[0] / vel[1]) * (wall0 - prevPos[1]) + prevPos[0]
                newVel = np.array([vel[0],-vel[1]])
                return np.array([collisionX, wall0]), newVel 
            elif Plate.inInterval(wall1,(prevPos[1], pos[1])):
                collisionX = (vel[0] / vel[1]) * (wall1 - prevPos[1]) + prevPos[0]
                newVel = np.array(vel[0],-vel[1]
                )
                return np.array([collisionX, wall1]), newVel
        return None 
    '''