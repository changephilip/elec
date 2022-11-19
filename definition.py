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

converseDistance=10.0
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
            "x":hpd.x - hpd.x.mean(),
            "y":hpd.y - hpd.y.mean(),
            "z":hpd.z - hpd.z.mean()
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
        self.h_len = h_len
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

        newCoord = eularTheta.apply(newCoord) +point - newMean
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


    
def initOrthoGroup(N):
    #from scipy.stats import special_ortho_group
    from scipy.spatial.transform import Rotation as R
    return R.random(N)

def initTwoFrames(ligand: hFrame, receptor: hFrame):
    "asure there is enough distance between two frames, init without spin"
    lBox=ligand.getBox()
    rBox=receptor.getBox()

    from scipy.spatial import distance as D
    directDistance = D.euclidean(ligand.g_mean,receptor.g_mean)

    patchBox=0.0
    if converseDistance < min(lBox,rBox):
        patchBox+= converseDistance
    else:
        patchBox+=converseDistance

    #adjust for test
    if directDistance < (lBox + rBox + converseDistance):
        directVector =  (ligand.g_mean - receptor.g_mean) + patchBox
        ligand.shift(directVector)
    #print(D.euclidean(ligand.g_mean,receptor.g_mean))
    from scipy.spatial.transform import Rotation as R
    #rotate=R.random()
    #ligand.rotateFixPoint(receptor.g_mean,rotate)
    #ligand.coord= ligand.coord - ligand.g_mean
    return False 

def initOrthoFrames(ligand: hFrame, receptor: hFrame):
    lBox = ligand.getBox()
    rBox = receptor.getBox()

    from scipy.spatial import distance as D
    directDistance = D.euclidean(ligand.g_mean, receptor.g_mean)

    patchBox = 0.0
    if directDistance < lBox + rBox + converseDistance:
        patchBox += converseDistance
        directVector = (ligand.g_mean - receptor.g_mean) + patchBox
        ligand.shift(directVector)

    from scipy.spatial.transform import Rotation as R
    orthoFrames=[]
    orthoGroup=initOrthoGroup(64)

    for r in orthoGroup:
        new_mean,new_coord = ligand.rotateFixPoint_(receptor.g_mean,r)
        orthoFrames.append(hFrame(new_mean,new_coord,ligand.charge))
    return orthoFrames

