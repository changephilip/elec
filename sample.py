from definition import *
from utils import *

class SA:
    "simulating annealing"
    frameStack=[]
    shiftStack=[]
    rotateStack=[]
    energyStack=[]
    gStack=[]
    def __init__(self,l:hFrame,r:hFrame):
        self.startFrame=l
        self.r=r
        self.backEnergy = calEU(l,r)
        #self.backEnergy = 65536.0
        self.energyStack.append(self.backEnergy)
        self.gStack.append(l.g_mean)
        self.backFrame=l
        self.currentFrame=l
        
        self.shift=numpy.array([1,1,1])
        self.rotate=numpy.array([0,0,0])

        
        self.currentEnergy=0.0

        self.T=300
        self.acceptRate=0

    def isTooNear(self):
        from scipy.spatial import distance as D
        for latom in self.currentFrame.coord:
            for ratom in self.r.coord:
                if D.euclidean(latom,ratom) < nearDistance:
                    return True
        return False

    def isNear(self):
        lBox = self.currentFrame.getBox()
        rBox = self.r.getBox()

        from scipy.spatial import distance as D
        directDistance = D.euclidean(self.currentFrame.g_mean,self.r.g_mean)

        if directDistance < (lBox + rBox + nearDistance):
            return False
        else:
            return self.isTooNear()

    def accept(self):
        self.currentEnergy = calEU(self.currentFrame,self.r)
        if self.currentEnergy > self.backEnergy:
            p=(- self.currentEnergy + self.backEnergy)/(self.backEnergy) 
            #print(-p*20)
            if numpy.exp(-p*30) > numpy.random.random():
                self.acceptRate+=1
                return True
            else:
                return False
        else:
            return True

    def shiftFrame(frame: hFrame,vector):
        frame.g_mean = frame.g_mean + vector

    
    def update(self):
        if self.isTooNear():
            if (self.accept()):
                self.backFrame=self.currentFrame
                self.backEnergy=self.currentEnergy
            
                self.shiftStack.append(self.shift)
                self.rotateStack.append(self.rotate)
                self.gStack.append(self.backFrame.g_mean)
                self.energyStack.append(self.currentEnergy)
                from copy import deepcopy
                self.frameStack.append(deepcopy(self.currentFrame))
    
    def stop(self):
        return True

    def move(self):
        #spin
        self.theta= smallTheta(1)
        self.currentFrame.spin(self.theta)
        
        #shift
        self.shift= smallVector(1)
        self.currentFrame.shift(self.shift)
        self.update()
        return True

    def testSA(self,N):
        for i in range(0,N):
            self.move()

        #print(self.gStack)
        #print(self.energyStack)
        print(len(self.energyStack))
        print(len(self.gStack))
        for j in range(0,len(self.energyStack)):
            print(self.energyStack[j])
            print(self.gStack[j])
        #print(len(self.energyStack))
        
        return self.frameStack
