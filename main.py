from sample import *
from definition import *
from utils import *

def runOneStart(ligand,receptor,N,name):
    
    l=Ligand(ligand)
    r=Ligand(receptor)
    lframe=hFrame(l.g_mean,l.h_coord,l.h_charge)
    rframe=hFrame(r.g_mean,r.h_coord,r.h_charge)
    
    print(calEU(lframe,rframe))
    initTwoFrames(lframe,rframe)
    writePQR(l,lframe,"init")
    writePQR(r,rframe,"r")

    sa=SA(lframe,rframe)
    fstack=sa.testSA(N)

    #writePQRSerial(l,fstack,"out/move")
    writeEnsemble(l,fstack,name)
    from scipy.spatial import distance as D
    for item in fstack:
        print(D.euclidean(item.g_mean,l.g_mean))
    return 0

def runOrthoStart(ligand,receptor, N ,name):
    l=Ligand(ligand)
    r=Ligand(receptor)
    lframe=hFrame(l.g_mean,l.h_coord,l.h_charge)
    rframe=hFrame(r.g_mean,r.h_coord,r.h_charge)
    
    print(calEU(lframe,rframe))
    startFrames=initOrthoFrames(lframe,rframe)
    
    writePQRSerial(l,startFrames,"out/start")

    #writePQRSerial(l,fstack,"out/move")
    #writeEnsemble(l,fstack,name)
    
    return 0

import sys
if __name__=="__main__":
    #cProfile.run('run(sys.argv[1],sys.argv[2],int(sys.argv[3]))')
    #runOneStart(sys.argv[1],sys.argv[2],int(sys.argv[3]),sys.argv[4])
    runOrthoStart(sys.argv[1],sys.argv[2],int(sys.argv[3]),sys.argv[4])
