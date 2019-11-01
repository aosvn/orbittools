import math
import numpy as np
from scipy.optimize import brentq
import sympy as symp 
from orbittools import *

# geom.py is a python adapation of Mark Tapley's OrbitTools.nb section C "In-Plane Geometry, Energy, Speed, Momentum"
# Adapted by Amanda Zangari from 8/4/2017 to 8/24/2017

def orbPerSP(mu,smAxis):
    return 2 * symp.pi * symp.sqrt(smAxis**3/mu)

def totOrbTest(alt):
    gmEarth=398600.433
    rEarth=6378.137
    totOrbitsAll=86164*12/orbPer(gmEarth,alt+rEarth)
    totOrbits=(totOrbitsAll-np.round(totOrbitsAll))
    return totOrbits


def main():
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import matplotlib as mpl

#test rOfTrueAnom
    
    
    npts = 100
    
    eccrange=np.linspace(0.0,3.0,npts)
    plotecc = (np.outer(eccrange,np.ones(npts))).flatten()
    anomrange = np.linspace(-np.pi/2,np.pi/2,npts)
    plotanom = (np.outer(np.ones(npts),anomrange)).flatten()
    
    plotsma = []
    
    for x in range(len(plotecc)):
        if plotecc[x] != 1.:
            plotsma.append(1./(1.-plotecc[x]**2.))
        else:
            plotsma.append(1.)
            
    rota = []
                        
    for x in range(len(plotecc)):
        rota.append(rOfTrueAnom(plotsma[x],plotecc[x],plotanom[x]))
                
#print rota
                

                
#no test for trueAnomOfR
                
#test norm
    print('Norm of [4,3] is',norm([4,3]))
                
                
#Test for totEnergy and totEnerFmRV
                
    gmTest=398600.433
    rTest=6378.137
    rOrbit=np.linspace(1,50,npts)
                
    tEFRV = []
                
                
    for x in range(len(rOrbit)):
        tEFRV.append(totEnerFmRV(gmTest,[rOrbit[x]*rTest,0,0],[0,np.sqrt(gmTest/(rOrbit[x]*rTest)),0]))
                    
    totE = []
    
    for x in range(len(rOrbit)):
        totE.append(totEnergy(gmTest,rOrbit[x]*rTest))
                        
    ax2=plt.subplot(111)
                        
    ax2.plot(rOrbit,tEFRV,'b')
    ax2.plot(rOrbit,totE,'r',linestyle='dashed')
    plt.show()
                        
#test passed, objects are identical
                        
#test for orbPer
                        
    gmTest=398600.433
    rTest=6378.137+200.0
                        
    print('Low Earth Orbit, Period should be 5309.64 seconds',orbPer(gmTest,rTest))
                        

#Test for orbPer

# orbPer had to be 
    print('Solve for Geosynchonous alitude, should be 42164.1 km)')
    symp.var('testSMAvar')
    print(symp.solve(orbPerSP(gmTest,testSMAvar)-86164.0,testSMAvar))



    daysToCheck=12
    alt=np.linspace(500,600,500)
    numDay=np.linspace(1,daysToCheck,daysToCheck)


    gmEarth=gmTest
    rEarth=rTest-200.0

#    pdb.set_trace()

    ax3=plt.subplot(111)

    cmap = mpl.cm.autumn

    for x in range(len(numDay)):
        totOrbitsAll=86164*numDay[x]/orbPer(gmEarth,alt+rEarth)
        totOrbits=np.absolute(totOrbitsAll-np.round(totOrbitsAll))
 #       totOrbits=(totOrbitsAll-np.round(totOrbitsAll))
        ax3.plot(alt,totOrbits,color=cmap(x/float(len(numDay))))

    plt.show()

                     
#    var('alts')
    thisday=12

    root=brentq(totOrbTest,520,530)
    print(root)

    gmTest= 398600.433
    rTest=6378.137
    print('Test of vCirc: vCirc(gmTest,200+rtest)')
    print(vCirc(gmTest,200.+rTest))

    rCirc=np.array(range(50))+1
    
    vCirctest=vCirc(gmTest,rCirc*rTest)

    ax2=plt.subplot(111)
    ax2.plot(rCirc,vCirctest)
    plt.show()


    plotRadius = np.array(range(12,200))/10.0*rTest
    vOR=vOfRadius(gmTest,12*rTest,plotRadius)

    ax2=plt.subplot(111)
    ax2.plot(plotRadius,vOR)
    plt.show()


    gmTest= 398600.433
    rTest=6378.137

    trueAn = np.pi*3.0*(np.array(range(201))/100.0-1.0)
    
    radVOTA = []
    tanVOTA = []

    for i in range(len(trueAn)):
        vOTA=(velOfTrueAnom(gmTest,2*rTest,0.3,trueAn[i]))
        radVOTA.append(vOTA[0])
        tanVOTA.append(vOTA[1])
#pdb.set_trace()


    ax2=plt.subplot(111)
    ax2.plot(trueAn,radVOTA,color='r')
    ax2.plot(trueAn,tanVOTA,color='k')
    plt.show()

#now test eccentri = 1

    trueAn = np.pi*(np.array(range(201))/100.0-1)

    radVOTA = []
    tanVOTA = []

    for i in range(len(trueAn)):
        vOTA=(velOfTrueAnom(gmTest,2*rTest,1.0,trueAn[i]))
        radVOTA.append(vOTA[0])
        tanVOTA.append(vOTA[1])
#pdb.set_trace()


    ax2=plt.subplot(111)
    ax2.plot(trueAn,radVOTA,color='r')
    ax2.plot(trueAn,tanVOTA,color='k')
    plt.show()
    


    hypAsympt = np.pi - (np.pi - 2.0 * np.arcsin(0.5))/2.0
    trueAn = hypAsympt*(np.array(range(201))/100.0-1)

    radVOTA = []
    tanVOTA = []

    for i in range(len(trueAn)):
        vOTA=(velOfTrueAnom(gmTest,-2.0*rTest,2.0,trueAn[i]))
        radVOTA.append(vOTA[0])
        tanVOTA.append(vOTA[1])
#pdb.set_trace()


    ax2=plt.subplot(111)
    ax2.plot(trueAn,radVOTA,color='r')
    ax2.plot(trueAn,tanVOTA,color='k')
    plt.show()

#test periApseRad and apoApseRad
    eccPlot=np.array(range(76))/100.0

    pAR = []
    aAR = []

    for i in range(len(eccPlot)):
        pAR.append(periApseRad(gmTest,4*rTest,angMomentum(gmTest,4*rTest,eccPlot[i])))
        aAR.append(apoApseRad(gmTest,4*rTest,angMomentum(gmTest,4*rTest,eccPlot[i])))

    ax2=plt.subplot(111)
    ax2.plot(eccPlot,pAR,color='r')
    ax2.plot(eccPlot,aAR,color='k')
    plt.show()
    
    #no test for mean motion but it is used in later functions with success.
    
if __name__ == "__main__":
    main()
                            
