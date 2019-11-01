import numpy as np
from scipy import optimize as opt #for tests
from orbittools import *

def vosmin(b):
    outBoundFmJup = np.array([ -6.98853123, -17.08305122 ,  0.9735782 ] )
    inBoundToJup = np.array([-11.34809306, -14.5795811,    0.38830128])

    rJupiter = 71492.0
    gmJupiter = 126712767.863       

    vos=np.linalg.norm(np.array(vOutSwingby(gmJupiter,inBoundToJup,rJupiter*b))-outBoundFmJup)
    return vos


def vosminvenus(b):
    venusVinfIn =  np.array([ 1.08553078,  4.19934229, -5.46907035])
    venusVinfOut = np.array([ 5.30065096, -1.48975192, -4.29060155])

    gmVenus    = 324858.599
    rVenus = 6051.8 

    vos=np.linalg.norm(np.array(vOutSwingby(gmVenus,venusVinfIn,rVenus*b))-venusVinfOut)
    return vos

def main():
    from mpl_toolkits.mplot3d import Axes3D
    #import dttmfunc as dttm
    import matplotlib.pyplot as plt
    import matplotlib as mpl


    varlist = origConstDict()
    secPerDay = varlist['secPerDay']
    
    print('Test of planetSwingby from Earth to Venus to Mercury')
    print('Results should be [-74.4484,-0.718316]')

    swingday = julDay(2020,9,10,0,0,0)

    print(np.array(planetSwingby("EMB",swingday-167.0,"Venus",swingday,"Mercury",swingday+84))/np.array([np.pi/180,1.0]))

    gmVenus    = 324858.599
    rVenus = 6051.8 


    print("Test of turnAngleOfPeriApse at Venus for above maneover")
    print("result should be 77.2214")
    print("turnAngleOfPeriApse[gmVenus,5.55107,rVenus+300.0]/Degree")
    print(turnAngleOfPeriApse(gmVenus,5.55107,rVenus+300.0)/np.pi*180.0)



    print("Test of vOutSwingby")
    
    testNHlaunch = julDay(2006,1,19,18,0,0)
    testNHswingby = julDay(2007,2,28,18,0,0)
    testNHencounter = julDay(2015,7,14,12,0,0)
    
    swingPlot = np.linspace(testNHswingby-2,testNHswingby+2,48)

    #'''
    # AMZ: This section commented out because it takes a long time to run.
    psPlot = []

    for x in range(len(swingPlot)):
        print(x)
        #print(planetSwingby("EMB",testNHlaunch,"Jupiter",swingPlot[x],"Pluto",testNHencounter)[1])
        psPlot.append(planetSwingby("EMB",testNHlaunch,"Jupiter",swingPlot[x],"Pluto",testNHencounter)[1])
    

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(swingPlot,psPlot)
    plt.show()

        
    #'''

    print('inBoundToJup = gaussPlanetary("EMB",testNHlaunch,"Jupiter",testNHswingby,0.][1]')
    print('result should be [-11.3481,-14.5796,0.388301]')
    inBoundToJup = (gaussPlanetary("EMB",testNHlaunch,"Jupiter",testNHswingby,0.))[1]
    print(inBoundToJup)


    print('outBoundToJup = gaussPlanetary("Jupiter",testNHswingby,"Pluto",testNHencounter,0.)[1]')
    print('results should be [-6.98853,-17.0831,0.973578]')
    outBoundFmJup = (gaussPlanetary("Jupiter",testNHswingby,"Pluto",testNHencounter,0.))[0]
    print(outBoundFmJup)

    print('np.linalg.norm(outBoundFmJup, should be 18.4829')
    print(np.linalg.norm(outBoundFmJup))
    
    
    print('np.linalg,norm(inBoundToJup, should be 18.4796')
    print(np.linalg.norm(inBoundToJup))

    rJupiter = 71492.0
    gmJupiter = 126712767.863       

    print('vOutSwingby(gmJupiter,inBoundToJup,np.array([3,40])*rJupiter)')
    print('answer should be [-7.26862,-16.9745,0.727204]')
    print(vOutSwingby(gmJupiter,inBoundToJup,np.array([3,40])*rJupiter))

    branges = np.array([[3.0,6.0],[36.0,39.0]])

    res = opt.minimize(vosmin,np.array([4.5,37.0]),method = 'TNC',tol = 1e-9)

    npts = 50

    b1range = np.linspace(3,6,npts)
    b2range = np.linspace(36,39,npts)

    plotb1 = (np.outer(b1range,np.ones(npts))).flatten()
    plotb2 = (np.outer(np.ones(npts),b2range)).flatten()    

    plotb = []

    for x in range(len(plotb1)):
        plotb.append(np.linalg.norm(vOutSwingby(gmJupiter,inBoundToJup,rJupiter*np.array([plotb1[x],plotb2[x]]))-np.array(outBoundFmJup)))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(plotb1,plotb2,plotb,'o')
    plt.show()

    
    print(res.x)
    print(res)


    swingday = julDay(2020,9,6,0,0,0) #note that this is different from above

    venusVinfIn = (gaussPlanetary("EMB", (swingday - 145),"Venus",swingday,0.))[1]
    venusVinfOut = (gaussPlanetary("Venus", swingday,"Mercury",swingday + 86.98,0.))[0]
    
    print('venusVinfIn, should be {1.08553,4.19934,-5.46907}',venusVinfIn)
    print('venusVinfOut, should be {5.29827,-1.49782,-4.38526}',venusVinfOut)

    vDiff = np.linalg.norm(venusVinfIn)-np.linalg.norm(venusVinfOut)
    print("vDiff, should be -0.0586391",vDiff)

    (turnAngleVen,vDiff) = np.array(planetSwingby("EMB",swingday - 145,"Venus",swingday,"Mercury",swingday+86.98))*np.array([180.0/np.pi,1.0])

    print('turnAngle,vDiff, should be 61.4991, -0.0586391')
    print(turnAngleVen,vDiff)

    print("turnAngleOfPeriApse(gmVenus,np.linalg.norm(venusVinfIn),(rVenus+300))*180/np.pi")
    print("should be 61.6103")
    print(turnAngleOfPeriApse(gmVenus,np.linalg.norm(venusVinfIn),(rVenus+300))*180/np.pi)

    print('missDist, should be 321.216')
    print(missDist(gmVenus,np.linalg.norm(venusVinfIn),turnAngleVen*np.pi/180)-rVenus)
    venusAtSwing=(posVelEcl("Venus",swingday))

    kmPerAU = 149597870.700
    print("venusAtSwing, should be (-23.9656,25.4872,1.73272)")
    venusAtSwing=((posVelEcl("Venus",swingday))[1]*kmPerAU/(secPerDay*36525))[0]
    print(venusAtSwing)


    bradrange = np.linspace(1,5,npts)
    bthetarange = np.linspace(-1.0*np.pi,1.0*np.pi,npts)

    plotbrad = (np.outer(bradrange,np.ones(npts))).flatten()
    plotbtheta = (np.outer(np.ones(npts),bthetarange)).flatten()

    #'''
    vosoutx = []
    vosouty = []
    vosoutz = []

    for x in range(len(plotbrad)):
        vosout = np.array(vOutSwingby(gmVenus,venusVinfIn,rVenus*np.array([plotbrad[x]*np.cos(plotbtheta[x]),plotbrad[x]*np.sin(plotbtheta[x])])))+np.array(venusAtSwing)
        vosoutx.append(vosout[0])
        vosouty.append(vosout[1])
        vosoutz.append(vosout[2])
    

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(vosoutx,vosouty,vosoutz,'o')
    plt.show()


    
    smackPlot = np.linspace(swingday+86.38,swingday+87.58,npts)

    gtoutx =[]
    gtouty =[]
    gtoutz =[]

    for x in range(len(smackPlot)):
        print(x)
        gtout = np.array((gaussPlanetary("Venus",swingday,"Mercury",smackPlot[x],0.))[0])+np.array(venusAtSwing)

        gtoutx.append(gtout[0])
        gtouty.append(gtout[1])
        gtoutz.append(gtout[2])
        

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(gtoutx,gtouty,gtoutz)
    plt.show()


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(gtoutx,gtouty,gtoutz)
    ax.scatter(vosoutx,vosouty,vosoutz,'o')
    plt.show()



    smackPlot2 = np.linspace(86.95,86.951,npts)+swingday

    psPlot = []
    for x in range(len(smackPlot2)):
        psPlot.append((planetSwingby("EMB",swingday - 145,"Venus",swingday,"Mercury",smackPlot2[x]))[1])
        print(x,smackPlot2[x],psPlot[x])


    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(smackPlot2-swingday,psPlot) #modfied so real results can be read off plot
    plt.show()
        
    #'''
    trueSmack = swingday + 86.9507

    venusVinfOut = (gaussPlanetary("Venus", swingday,"Mercury",trueSmack,0.))[0]
    print('venusVinfOut = (gaussPlanetary("Venus", swingday,"Mercury",trueSmack,0.))[0]')
    print('should be [5.30065,-1.48975,-4.2906]')
    print(venusVinfOut)
 


    branges = np.array([[-1,-0.5],[-1.5,-1]])
    res = opt.minimize(vosminvenus,np.array([-.5,-1.5]),method = 'TNC',tol = 1e-9)


    print(res.x)
    print(res)

    bTarg = np.array(res.x)

    periAOB = periApseOfB(gmVenus,np.linalg.norm(venusVinfIn),rVenus*bTarg)-rVenus
    print('periApseOfB(gmVenus,np.linalg.norm(venusVinfIn),rVenus*bTarg)-rVenus')
    print('should be 248.537')
    print(periAOB)


if __name__ == "__main__":
    main()

