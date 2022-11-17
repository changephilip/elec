#from multiprocessing.sharedctypes import RawArray
from platform import processor
import cProfile
import prody
import numpy
import pandas
import math
import scipy

#from sympy import N

# constant 
e0=78.4
fpi=4.0*math.pi
#theta=5.0
distance=2.0
nProcess=4

converseDistance=20.0
nearDistance=8.0

class Ligand:
    #relative coordinate system to geometry center
    def __init__(self,pqr):
        atoms = prody.parsePQR(pqr)

        resname= atoms.getResnames()
        atomType=atoms.getNames()

        x= atoms.getCoords().transpose()[0]
        y= atoms.getCoords().transpose()[1]
        z= atoms.getCoords().transpose()[2]
        charge=atoms.getCharges()
        radius=atoms.getRadii()

        pd=pandas.DataFrame({
            "atom":atomType,
            "res":resname,
            "x":x,
            "y":y,
            "z":z,
            "q":charge,
            "r":radius
        })
        pd["isH"]=pd["atom"].apply(lambda x: not (x=="Hlp" or x =="Hp" or x=="Hh" or x=="Hsh"))
        g_mean=numpy.array([x.mean(),y.mean(),z.mean()])
        coord=pandas.DataFrame({
            "x":x-x.mean(),
            "y":y-y.mean(),
            "z":z-z.mean()
        })

        hpd=pd[pd["isH"]==True]
        h_coord=numpy.array(pandas.DataFrame({
            "x":hpd.x,
            "y":hpd.y,
            "z":hpd.z
        }))
        hypd=pd[pd["isH"]==False]
        hy_coord=numpy.array(pandas.DataFrame({
            "x":hypd.x,
            "y":hypd.y,
            "z":hypd.z
        }))

        h_charge=numpy.array(pandas.DataFrame({
            "q":pd[pd["isH"]==True].q
        }))

        hy_charge=numpy.array(pandas.DataFrame({
            "q":pd[pd["isH"]==False].q
        }))
        h_len = len(h_coord)
        hpqr=prody.atomic.atomgroup.AtomGroup()
        hpqr.setResnums(atoms.getResnums()[0:h_len])
        hpqr.setNames(atomType[0:h_len])
        hpqr.setResnames(resname[0:h_len])
        hpqr.setCoords(h_coord)
        #hpqr.setCharges(h_charge)

        self.pqr=pd
        self.g_mean=g_mean
        self.coord=coord
        self.hpd=hpd
        self.h_coord=h_coord
        self.hypd=hypd
        self.hy_coord=hy_coord
        self.h_charge=h_charge
        self.hy_charge=hy_charge
        self.atoms = atoms
        self.hpqr=hpqr
        #return Ligand(pd,g_mean,coord,hpd,h_coord,hypd,hy_coord,h_charge,hy_charge)


class hFrame:
    def __init__(self,g_mean,coord,charge):
        self.g_mean=numpy.copy(g_mean)
        self.coord=numpy.copy(coord)
        self.charge=numpy.copy(charge)

    def spin_(self,theta):
        from scipy.spatial.transform import Rotation as R
        r= R.from_euler('xyz',theta,degrees=True)
        return r.apply(self.coord)

    def spin(self,theta):
        from scipy.spatial.transform import Rotation as R
        r= R.from_euler('xyz',theta,degrees=True)
        self.coord=r.apply(self.coord)

    def shift_(self,vector):
        return self.g_mean + vector 

    def shift(self,vector):
        self.g_mean = self.g_mean + vector

    def rotateFixPoint_(self,point,eularTheta):
        newCoord = self.coord +self.g_mean - point

        newMean = eularTheta.apply(self.g_mean - point) + point

        newCoord = eularTheta.apply(newCoord) +point - self.g_mean
        return newMean,newCoord

    def rotateFixPoint(self,point,eularTheta):
        newCoord = self.coord +self.g_mean - point

        newMean = eularTheta.apply(self.g_mean - point) + point

        newCoord = eularTheta.apply(newCoord) + point - newMean

        self.coord = newCoord
        self.g_mean = newMean


    def randomRotate(self):
        return 0

    def getBox(self):
        RBox=[]
        from scipy.spatial import distance as D
        for item in self.coord:
            RBox.append(D.euclidean(item,self.g_mean))
        return max(RBox)

def writePQR(pqr: Ligand,frame: hFrame,name):
    pqr.hpqr.setCoords(frame.coord + frame.g_mean)
    prody.writePQR(name+".pqr",pqr.hpqr)
    return True

def writePQRSerial(pqr:Ligand, frameList,name):
    i=0
    for frame in frameList:
        writePQR(pqr,frame,name+"_"+str(i))
        i=i+1
    return True
    
def dielectric(point1,point2):
    "return e*r*r 1/(e*r^2)"
    #print(point1)
    #print(point2)
    from scipy.spatial import distance as D
    r= D.euclidean(point1,point2)
    A= -8.5525
    B= e0 - A
    k= 7.7839
    l= 0.003627
    e= A + (B/(1+k*math.exp(-1.0*l*B*r)))
    return 1/(e*r*r)

def dielectricP(point1,point2,value,i,n1,n2):
    value=dielectric(point1[n1],point2[n2])
    
from functools import partial
def dielectric_wrapper(results,set1,set2,N):
    from multiprocessing import Pool
    p=Pool(nProcess)
    origin_shape=results.shape
    
    results=range(0,63)
    for i in range(0,len(N)):
        p.apply_async(dielectricP,args=(set1.h_coord,set2.h_coord,results[i],i,N[i][0],N[i][1]))
        
    
    p.close()
    p.join()


def smallTheta(N):
    return numpy.random.uniform(-60.0,60.0,[3])
    
def smallVector(N):
    return numpy.random.uniform(-5.0,5.0,[3])

def calEU(ligand1: hFrame,ligand2: hFrame):
    U=0.0
    aux1_q1q2=numpy.outer(ligand1.charge,ligand2.charge)

    N1=numpy.arange(len(ligand1.coord))
    N2=numpy.arange(len(ligand2.coord))
    N=[]
  
    for i in N1:
      for j in N2:
        N.append((i,j))
   
    aux2_r=numpy.ones((len(ligand1.coord),len(ligand2.coord)))
    #print(N1.shape)
    #print(N2.shape)

    
    for item in N:
        aux2_r[item[0]][item[1]]=dielectric(ligand1.coord[item[0]]+ligand1.g_mean,ligand2.coord[item[1]]+ligand2.g_mean)
    
    "q1*q2/(r^2)"
    aux3_q1q2_r2=aux2_r*aux1_q1q2
    U=numpy.sum(aux3_q1q2_r2)

    return -U/fpi

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
        self.energyStack.append(self.backEnergy)

        self.backFrame=l
        self.currentFrame=l
        
        self.shift=numpy.array([1,1,1])
        self.rotate=numpy.array([0,0,0])

        self.backEnergy=0.0
        self.currentEnergy=0.0

        self.T=300
        self.acceptRate=0.0

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
            #print(p)
            if numpy.exp(-p) > numpy.random.random():
                return True
            else:
                return False
        else:
            return True

    def shiftFrame(frame: hFrame,vector):
        frame.g_mean = frame.g_mean + vector

    
    def update(self):
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
        #print(self.rotateStack)
        #print(self.shiftStack)
        print(self.energyStack)
        print(len(self.energyStack))
        
        return self.frameStack
    

def initTwoFrames(ligand: hFrame, receptor: hFrame):
    "asure there is enough distance between two frames, init without spin"
    lBox=ligand.getBox()
    rBox=receptor.getBox()

    from scipy.spatial import distance as D
    directDistance = D.euclidean(ligand.g_mean,receptor.g_mean)

    patchBox=0.0
    if converseDistance < min(lBox,rBox):
        patchBox+= (lBox + rBox)*4.0
    else:
        patchBox+=converseDistance

    #adjust for test
    if directDistance < (lBox + rBox + converseDistance):
        directVector =  (ligand.g_mean - receptor.g_mean)/(directDistance) * (patchBox)
        ligand.shift(directVector)
    
    from scipy.spatial.transform import Rotation as R
    rotate=R.random()
    ligand.rotateFixPoint(receptor.g_mean,rotate)
    
        
    return False 



def run(ligand,receptor):
    
    l=Ligand(ligand)
    r=Ligand(receptor)
    lframe=hFrame(l.g_mean,l.h_coord,l.h_charge)
    rframe=hFrame(r.g_mean,r.h_coord,r.h_charge)
    
    initTwoFrames(lframe,rframe)
    
    writePQR(l,lframe,"init.pqr")

    sa=SA(lframe,rframe)
    fstack=sa.testSA(1)
    #print(fstack)
    writePQRSerial(l,fstack,"move")

    #print(calEU(lframe,rframe))
    return 0


import sys
if __name__=="__main__":
    cProfile.run('run(sys.argv[1],sys.argv[2])')