from importPqr import *

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
    
def writeEnsemble(pqr:Ligand, frameList,name):
    e = prody.PDBEnsemble(pqr.atoms[0:pqr.h_len])
    
    e.delCoordset(0)
    for frame in frameList:
        #print(frame.coord+frame.g_mean)
        e.addCoordset(frame.coord+frame.g_mean)
    
    #prody.writePDB(name[0:]+".pdb",e)
    prody.writeDCD(name+".dcd",e)

def dielectric(point1,point2):
    "return e*r*r 1/(e*r^2)"
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
    return numpy.random.uniform(-3.0,3.0,[3])
    
def smallVector(N):
    return numpy.random.uniform(-1.0,1.0,[3])


from interface import lib
#numpy.float_ = numpy.float32

def calEU(ligand1: hFrame,ligand2: hFrame):
    U=0.0
    aux1_q1q2=numpy.outer(ligand1.charge,ligand2.charge)

    N1=numpy.arange(len(ligand1.coord))
    N2=numpy.arange(len(ligand2.coord))
  
    
    aux2_r=numpy.ones((len(ligand1.coord),len(ligand2.coord))).astype('float32')

    
    lib.calc_wrap(ligand1.coord.astype('float32'),ligand2.coord.astype('float32'),ligand1.g_mean.astype('float32'),ligand2.g_mean.astype('float32'),len(N1),len(N2),aux2_r)      

    aux2_r = aux2_r.astype('float64')
    
    "q1*q2/(r^2)"
    aux3_q1q2_r2=aux2_r*aux1_q1q2
    U=numpy.sum(aux3_q1q2_r2)

    return -U/fpi
