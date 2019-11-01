# Translation of Mark Tapley's OrbitTools.nb, functions only
# Module written by Amanda Zangari, September 26, 2017
# Requires Python 3!!!!!!

#import other necessary python libraries
import copy
import math
import matplotlib
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pdb
import pickle
import re
import requests,os,sqlite3,datetime
#import spiceypy as spice
import time


from html.parser import HTMLParser
from scipy import optimize as opt



#GAUSS PLANETARY HAS HARD CODED CONSTANTS, AS DOES TESTCOMPAREMAT (a test, probably, test function removed totOrbPer

#Section A: Constants

#Note. first function, which defines global variables is original.  Other functions were added for more up-to-date testing.

def origConst():
   #Lists of Planet masses, semi-major axes of orbits, Planetary radii, and various useful constants and conversion factors. #Fundamental Units
#AU in km is a defined (exact) value, meaning semi-major axis of Earth's orbit is not exactly 1 AU.
#From http://neo.jpl.nasa.gov/glossary/au.html

    global kmPerAU
    kmPerAU = 149597870.700

#Speed of light in km/s is an exact value
#From http://physics.nist.gov/cuu/Constants/index.html

    global cLight
    cLight = 299792.458

#Calculation of seconds/century, convenient conversion factor for use with planetary ephems in (I).

    global secPerCy
    secPerCy = 86400.0 * 36525.0

#Seconds per day, 86400
    global secPerDay
    secPerDay = 86400.0

#"Gravitational constant in km^3 / (kg s^2)
#from http://physics.nist.gov/cuu/Constants/index.html

    global gravConst
    gravConst = 6.67384e-20

#Planetary Constants
    
#Planetary masses
    
#Planetary Gravitational Constants from DE405 (according to old notebook) except Pluto, revised value from 
#    http://ssd.jpl.nasa.gov/?planet_phys_par
#All units are km^3/s^2
    
    global gmSun      
    gmSun              = 132712440017.987 
    global gmMercury   
    gmMercury  = 22032.080
    global gmVenus     
    gmVenus    = 324858.599
    global gmEarth     
    gmEarth    = 398600.433
    global gmMars      
    gmMars     = 42828.314
    global gmJupiter    
    gmJupiter  = 126712767.863
    global gmSaturn    
    gmSaturn   = 37940626.063
    global gmUranus    
    gmUranus   = 5794549.007
    global gmNeptune   
    gmNeptune  = 6836534.064
    global gmPluto     
    gmPluto     =873.663
    

#Planetary orbit semi-major axes from wikipedia
    
#Not Consistent with AU or DE405 values! Do not Match planetary Ephem in (I) below! Cite Standish, get 2020 values, and refer to I for plots of variation.
    
    global smaMercury
    smaMercury = 57909100.0  
    global smaEarth
    smaEarth = 149597887.5 
    global smaVenus 
    smaVenus = 108208930.
    global smaMars
    smaMars = 1.523679 * smaEarth
    global smaJupiter
    smaJupiter = 5.204267 * smaEarth
    
#Planetary radius values from wikipedia
#Need to complete table
#Units are km


    global rMercury
    rMercury = 2439.7 
    global rVenus 
    rVenus = 6051.8 
    global rEarth 
    rEarth = 6378.137 
    global rJupiter 
    rJupiter = 71492.0 


#Planetary Oblateness
#Value for Earth from http://en.wikipedia.org/wiki/Nodal_precession.

    global j2Mercury
    j2Mercury = 0.0
    global j2Earth
    j2Earth = 0.00108262668


#moon constants

    global gmLuna
    gmLuna= 4902.800013
    global gmGany 
    gmGany = 0.025*gmEarth

    global rGany
    rGany = 2631.2
#Propellant/Rocket System constants
#Rough estimates; STAR solids based on STAR catalog
#Values in seconds (i.e. N-s over N weight of propellant on Earth)

# monopro
    global iSPHydrazpellant
    iSPHydraz = 215 
# STAR 
    global iSPSTARcatalog 
    iSPSTAR = 295
 
    return


def origConstDict():
    smaEarth = 149597887.5
    gmEarth = 398600.433

    constdict = {'kmPerAU': 149597870.700,
                 'cLight': 299792.458,
                 'secPerCy': 86400.0 * 36525.0,
                 'secPerDay': 86400.0,
                 'gravConst': 6.67384e-20,
                 'gmSun': 132712440017.987 ,
                 'gmMercury': 22032.080,
                 'gmVenus': 324858.599,
                 'gmEarth': 398600.433,
                 'gmMars': 42828.314,
                 'gmJupiter': 126712767.863,
                 'gmSaturn': 37940626.063,
                 'gmUranus': 5794549.007,
                 'gmNeptune': 6836534.064,
                 'gmPluto':873.663,
                 'smaMercury': 57909100.0  ,
                 'smaEarth': 149597887.5 ,
                 'smaVenus': 108208930.,
                 'smaMars': 1.523679 * smaEarth,
                 'smaJupiter': 5.204267 * smaEarth,
                 'rMercury': 2439.7 ,
                 'rVenus': 6051.8 ,
                 'rEarth': 6378.137 ,
                 'rJupiter': 71492.0 ,
                 'rSaturn': 60268.0,
                 'rUranus': 25559.0,
                 'rNeptune':24764.0,
                 'j2Mercury': 0.0,
                 'j2Earth': 0.00108262668,
                 'gmLuna': 4902.800013,
                 'gmGany': 0.025*gmEarth,
                 'rGany': 2631.2,
                 'iSPHydraz': 215 ,
                 'iSPSTAR': 295}
 
    return constdict

def newConstDict():
    smaEarth = 149597887.5
    gmEarth = 398600.433

    constdict = {'kmPerAU': 149597870.700,
                 'cLight': 299792.458,
                 'secPerCy': 86400.0 * 36525.0,
                 'secPerDay': 86400.0,
                 'gravConst': 6.67408e-20,
                 'gmSun': 1.3271244004193938E+11,
                 'gmMercurySys': 2.2031780000000021E+04,
                 'gmVenusSys': 3.2485859200000006E+05,
                 'gmEarthSys': 4.0350323550225981E+05,
                 'gmMarsSys': 4.2828375214000022E+04,
                 'gmJupiterSys': 1.2671276480000021E+08,
                 'gmSaturnSys': 3.7940585200000003E+07,
                 'gmUranusSys': 5.7945486000000080E+06,
                 'gmNeptuneSys': 6.8365271005800236E+06,
                 'gmPlutoSys': 9.7700000000000068E+02,
                 'gmMercury': 2.2031780000000021E+04,
                 'gmVenus': 3.2485859200000006E+05,
                 'gmEarth': 3.9860043543609598E+05,
                 'gmMars': 4.282837362069909E+04,
                 'gmJupiter': 1.266865349218008E+08,
                 'gmSaturn': 3.793120749865224E+07,
                 'gmUranus': 5.793951322279009E+06,
                 'gmNeptune': 6.835099502439672E+06,
                 'gmPluto': 8.696138177608748E+02,
                 'gmLuna': 4.9028000661637961E+03,
                 'gmPhobos': 7.087546066894452E-04,
                 'gmDeimos': 9.615569648120313E-05,
                 'gmIo': 5.959916033410404E+03,
                 'gmEuropa': 3.202738774922892E+03,
                 'gmGanymede': 9.887834453334144E+03,
                 'gmCallisto': 7.179289361397270E+03,
                 'gmAmalthea': 1.378480571202615E-01,
                 'gmMimas': 2.503522884661795E+00,
                 'gmEnceladus': 7.211292085479989E+00,
                 'gmTethys': 4.121117207701302E+01,
                 'gmDione': 7.311635322923193E+01,
                 'gmRhea': 1.539422045545342E+02,
                 'gmTitan': 8.978138845307376E+03,
                 'gmHyperion': 3.718791714191668E-01,
                 'gmIapetus': 1.205134781724041E+02,
                 'gmPhoebe': 5.531110414633374E-01,
                 'gmJanus': 1.266231296945636E-01,
                 'gmEpimetheus': 3.513977490568457E-02,
                 'gmAtlas': 3.759718886965353E-04,
                 'gmPrometheus': 1.066368426666134E-02,
                 'gmPandora': 9.103768311054300E-03,
                 'gmAriel': 8.346344431770477E+01,
                 'gmUmbriel': 8.509338094489388E+01,
                 'gmTitania': 2.269437003741248E+02,
                 'gmOberon': 2.053234302535623E+02,
                 'gmMiranda': 4.319516899232100E+00,
                 'gmTriton': 1.427598140725034E+03,
                 'gmCharon': 1.058799888601881E+02,
                 'gmNix': 3.048175648169760E-03,
                 'gmHydra': 3.211039206155255E-03,
                 'gmStyx': 1.110040850536676E-03,
                 'gmCeres': 6.3130000000000003E+01,
                 'gmPallas': 1.3730000000000000E+01,
                 'gmJuno': 1.8200000000000001E+00,
                 'gmVesta': 1.7289999999999999E+01,
                 'gmHebe': 9.3000000000000005E-01,
                 'gmIris': 8.5999999999999999E-01,
                 'gmHygiea': 5.7800000000000002E+00,
                 'gmEunomia': 2.1000000000000001E+00,
                 'gmPsyche': 1.8100000000000001E+00,
                 'gmAmphitrite': 8.5999999999999999E-01,
                 'gmEuropaAsteroid': 1.5900000000000001E+00,
                 'gmCybele': 9.1000000000000003E-01,
                 'gmSylvia': 9.8999999999999999E-01,
                 'gmThisbe': 1.0200000000000000E+00,
                 'gmEros':  4.463E-4,
                 'gmDavida': 2.2599999999999998E+00,
                 'gmInteramnia': 2.1899999999999999E+00,
                 'rMercury': 2439.7,
                 'rVenus': 6051.8,
                 'rEarth': 6378.137,
                 'rMars':3396.19,
                 'rJupiter': 71492.0,
                 'rSaturn': 60268.0,
                 'rUranus': 25559.0,
                 'rNeptune':24764.0,
                 'rPluto':1188.0,
                 'j2Mercury': 0.0,
                 'j2Earth': 0.00108262668,
                 'rGany': 2631.2,
                 'iSPHydraz': 215,
                 'iSPSTAR': 295}
 
    return constdict


def fundamentalUnits():

    # from https://cneos.jpl.nasa.gov/glossary/au.html (new link found by Amanda)
    global kmPerAU
    kmPerAU = 149597870.700

    #Speed of light in km/s is an exact value
    #from http://physics.nist.gov/cuu/Constants/index.html 
    # more specfically https://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=C
    global cLight
    cLight = 299792.458

    #Calculation of seconds per century, convenient converstion factor for use with planetary ephems i (I)
    global secPerCy 
    secPerCy = 86400 * 36525
    
    #Seconds per year
    #https://www.grc.nasa.gov/www/k-12/Numbers/Math/Mathematical_Thinking/calendar_calculations.htm
    secperyear = 365*24*3600.0+5*3600.0+48*60.0+46.0

    #Seconds per day, 86400
    global secPerDay
    secPerDay = 86400.0
    
    #Gravitational constant in km^3/(kg s^2)
    #from http://physics.nist.gov/cuu/Constants/index.html
    # more specifcally https://physics.nist.gov/cgi-bin/cuu/Value?bg|search_for=gravitational+constant
    # Note, value is updated from orbittools.nb
    # difference likely due because code was started before "2014 CODATA recommended values" came out
    global gravConst
    gravConst = 6.67384e-20 #(note this was in original, useful for testing errors?)
    #gravConst = 6.67408e-20

def planetaryConstants():
    # all units in km^3/s^2
    # Planetary constants from https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc
    # ^This is a new source from DE405, last updated in 2014
    # Improved to both use SPICE numbers and common names for planets

    # Planetary constants and moon constants have been merged from original file
    #we hope to upgrage code to just directly use SPICE and not these hardcoded values

    global BODY1_GM
    BODY1_GM       = ( 2.2031780000000021E+04 )
    global BODY2_GM
    BODY2_GM       = ( 3.2485859200000006E+05 )
    global     BODY3_GM
    BODY3_GM       = ( 4.0350323550225981E+05 )
    global     BODY4_GM
    BODY4_GM       = ( 4.2828375214000022E+04 )
    global     BODY5_GM
    BODY5_GM       = ( 1.2671276480000021E+08 )
    global     BODY6_GM
    BODY6_GM       = ( 3.7940585200000003E+07 )
    global     BODY7_GM
    BODY7_GM       = ( 5.7945486000000080E+06 )
    global     BODY8_GM
    BODY8_GM       = ( 6.8365271005800236E+06 )
    global     BODY9_GM
    BODY9_GM       = ( 9.7700000000000068E+02 )
    global     BODY10_GM
    BODY10_GM      = ( 1.3271244004193938E+11 )

    global      BODY199_GM
    BODY199_GM     = ( 2.2031780000000021E+04 )
    global      BODY299_GM
    BODY299_GM     = ( 3.2485859200000006E+05 )
    
    global     BODY399_GM
    BODY399_GM     = ( 3.9860043543609598E+05 )
    global     BODY499_GM
    BODY499_GM     = ( 4.282837362069909E+04  )
    global     BODY599_GM
    BODY599_GM     = ( 1.266865349218008E+08  )
    global     BODY699_GM
    BODY699_GM     = ( 3.793120749865224E+07  )
    global     BODY799_GM
    BODY799_GM     = ( 5.793951322279009E+06  )
    global     BODY899_GM
    BODY899_GM     = ( 6.835099502439672E+06  )
    global     BODY999_GM
    BODY999_GM     = ( 8.696138177608748E+02  )
    
    global     BODY301_GM
    BODY301_GM     = ( 4.9028000661637961E+03 )
    
    global     BODY401_GM
    BODY401_GM     = ( 7.087546066894452E-04 )
    global     BODY402_GM
    BODY402_GM     = ( 9.615569648120313E-05 )
    
    global     BODY501_GM
    BODY501_GM     = ( 5.959916033410404E+03 )
    global     BODY502_GM
    BODY502_GM     = ( 3.202738774922892E+03 )
    global     BODY503_GM
    BODY503_GM     = ( 9.887834453334144E+03 )
    global     BODY504_GM
    BODY504_GM     = ( 7.179289361397270E+03 )
    global     BODY505_GM
    BODY505_GM     = ( 1.378480571202615E-01 )
    
    global     BODY601_GM
    BODY601_GM     = ( 2.503522884661795E+00 )
    global     BODY602_GM
    BODY602_GM     = ( 7.211292085479989E+00 )
    global     BODY603_GM
    BODY603_GM     = ( 4.121117207701302E+01 )
    global     BODY604_GM
    BODY604_GM     = ( 7.311635322923193E+01 )
    global     BODY605_GM
    BODY605_GM     = ( 1.539422045545342E+02 )
    global     BODY606_GM
    BODY606_GM     = ( 8.978138845307376E+03 )
    global     BODY607_GM
    BODY607_GM     = ( 3.718791714191668E-01 )
    global     BODY608_GM
    BODY608_GM     = ( 1.205134781724041E+02 )
    global     BODY609_GM
    BODY609_GM     = ( 5.531110414633374E-01 )
    global     BODY610_GM
    BODY610_GM     = ( 1.266231296945636E-01 )
    global     BODY611_GM
    BODY611_GM     = ( 3.513977490568457E-02 )
    global     BODY615_GM
    BODY615_GM     = ( 3.759718886965353E-04 )
    global     BODY616_GM
    BODY616_GM     = ( 1.066368426666134E-02 )
    global     BODY617_GM
    BODY617_GM     = ( 9.103768311054300E-03 )
    
    global     BODY701_GM
    BODY701_GM     = ( 8.346344431770477E+01 )
    global     BODY702_GM
    BODY702_GM     = ( 8.509338094489388E+01 )
    global     BODY703_GM
    BODY703_GM     = ( 2.269437003741248E+02 )
    global     BODY704_GM
    BODY704_GM     = ( 2.053234302535623E+02 )
    global     BODY705_GM
    BODY705_GM     = ( 4.319516899232100E+00 )
    
    global BODY801_GM
    BODY801_GM     = ( 1.427598140725034E+03 )
    
    global     BODY901_GM
    BODY901_GM     = ( 1.058799888601881E+02 )
    global     BODY902_GM
    BODY902_GM     = ( 3.048175648169760E-03 )
    global     BODY903_GM
    BODY903_GM     = ( 3.211039206155255E-03 )
    global     BODY904_GM 
    BODY904_GM     = ( 1.110040850536676E-03 )
    
    global     BODY2000001_GM
    BODY2000001_GM = ( 6.3130000000000003E+01 )
    global     BODY2000002_GM
    BODY2000002_GM = ( 1.3730000000000000E+01 )
    global     BODY2000003_GM
    BODY2000003_GM = ( 1.8200000000000001E+00 )
    global     BODY2000004_GM
    BODY2000004_GM = ( 1.7289999999999999E+01 )
    global     BODY2000006_GM
    BODY2000006_GM = ( 9.3000000000000005E-01 )
    global     BODY2000007_GM
    BODY2000007_GM = ( 8.5999999999999999E-01 )
    global     BODY2000010_GM
    BODY2000010_GM = ( 5.7800000000000002E+00 )
    global     BODY2000015_GM
    BODY2000015_GM = ( 2.1000000000000001E+00 )
    global     BODY2000016_GM
    BODY2000016_GM = ( 1.8100000000000001E+00 )
    global     BODY2000029_GM
    BODY2000029_GM = ( 8.5999999999999999E-01 )
    global     BODY2000052_GM
    BODY2000052_GM = ( 1.5900000000000001E+00 )
    global     BODY2000065_GM
    BODY2000065_GM = ( 9.1000000000000003E-01 )
    global     BODY2000087_GM
    BODY2000087_GM = ( 9.8999999999999999E-01 )
    global     BODY2000088_GM
    BODY2000088_GM = ( 1.0200000000000000E+00 )
    global     BODY2000433_GM
    BODY2000433_GM = (  4.463E-4 )
    global     BODY2000511_GM
    BODY2000511_GM = ( 2.2599999999999998E+00 )
    global     BODY2000704_GM
    BODY2000704_GM = ( 2.1899999999999999E+00 )



    global gmMercurySys
    gmMercurySys       = ( 2.2031780000000021E+04 )
    global gmVenusSys
    gmVenusSys       = ( 3.2485859200000006E+05 )
    global     gmEarthSys
    gmEarthSys       = ( 4.0350323550225981E+05 )
    global     gmMarsSys
    gmMarsSys       = ( 4.2828375214000022E+04 )
    global     gmJupiterSys
    gmJupiterSys       = ( 1.2671276480000021E+08 )
    global     gmSaturnSys
    gmSaturnSys       = ( 3.7940585200000003E+07 )
    global     gmUranusSys
    gmUranusSys       = ( 5.7945486000000080E+06 )
    global     gmNeptuneSys
    gmNeptuneSys       = ( 6.8365271005800236E+06 )
    global     gmPlutoSys
    gmPlutoSys      = ( 9.7700000000000068E+02 )
    global    gmSun
    gmSun      = ( 1.3271244004193938E+11 )

    global      gmMercury
    gmMercury     = ( 2.2031780000000021E+04 )
    global      gmVenus
    gmVenus     = ( 3.2485859200000006E+05 )
    
    global     gmEarth
    gmEarth     = ( 3.9860043543609598E+05 )
    global     gmMars
    gmMars     = ( 4.282837362069909E+04  )
    global     gmJupiter
    gmJupiter     = ( 1.266865349218008E+08  )
    global    gmSaturn
    gmSaturn     = ( 3.793120749865224E+07  )
    global     gmUranus
    gmUranus     = ( 5.793951322279009E+06  )
    global     gmNeptune
    gmNeptune     = ( 6.835099502439672E+06  )
    global     gmPluto
    gmPluto     = ( 8.696138177608748E+02  )
    
    global     gmLuna
    gmLuna     = ( 4.9028000661637961E+03 )
    global     gmMoon
    gmMoon = gmLuna
    
    global     gmPhobos
    gmPhobos     = ( 7.087546066894452E-04 )
    global     gmDeimos
    gmDeimos     = ( 9.615569648120313E-05 )
    
    global     gmIo
    gmIo     = ( 5.959916033410404E+03 )
    global     gmEuropa
    gmEuropa     = ( 3.202738774922892E+03 )
    global     gmGanymede
    gmGanymede     = ( 9.887834453334144E+03 )
    global     gmCallisto
    gmCallisto     = ( 7.179289361397270E+03 )
    global     gmAmalthea
    gmAmalthea     = ( 1.378480571202615E-01 )
    
    global     gmMimas
    gmMimas     = ( 2.503522884661795E+00 )
    global     gmEnceladus
    gmEnceladus     = ( 7.211292085479989E+00 )
    global     gmTethys
    gmTethys     = ( 4.121117207701302E+01 )
    global     gmDione
    gmDione     = ( 7.311635322923193E+01 )
    global     gmRhea
    gmRhea     = ( 1.539422045545342E+02 )
    global     gmTitan
    gmTitan     = ( 8.978138845307376E+03 )
    global     gmHyperion
    gmHyperion     = ( 3.718791714191668E-01 )
    global     gmIapetus
    gmIapetus     = ( 1.205134781724041E+02 )
    global     gmPhoebe
    gmPhoebe     = ( 5.531110414633374E-01 )
    global     gmJanus
    gmJanus     = ( 1.266231296945636E-01 )
    global     gmEpimetheus_GM
    gmEpimetheus     = ( 3.513977490568457E-02 )
    global     gmAtlas
    gmAtlas     = ( 3.759718886965353E-04 )
    global     gmPrometheus
    gmPrometheus     = ( 1.066368426666134E-02 )
    global     gmPandora
    gmPandora     = ( 9.103768311054300E-03 )
    
    global     gmAriel
    gmAriel     = ( 8.346344431770477E+01 )
    global     gmUmbriel
    gmUmbriel     = ( 8.509338094489388E+01 )
    global     gmTitania
    gmTitania     = ( 2.269437003741248E+02 )
    global     gmOberon
    gmOberon     = ( 2.053234302535623E+02 )
    global     gmMiranda
    gmMiranda     = ( 4.319516899232100E+00 )
    
    global gmTriton
    gmTitan     = ( 1.427598140725034E+03 )
    
    global     gmCharon
    gmCharon     = ( 1.058799888601881E+02 )
    global     gmNix
    gmNix     = ( 3.048175648169760E-03 )
    global     gmHydra
    gmHydra     = ( 3.211039206155255E-03 )
    global     gmStyx
    gmStyx     = ( 1.110040850536676E-03 )
    
    global     gmCeres
    gmCeres = ( 6.3130000000000003E+01 )
    global     gmPallas
    gmPallas = ( 1.3730000000000000E+01 )
    global     gmJuno
    gmJuno = ( 1.8200000000000001E+00 )
    global     gmVesta
    gmVesta = ( 1.7289999999999999E+01 )
    global     gmHebe
    gmHebe = ( 9.3000000000000005E-01 )
    global     gmIris
    gmIris = ( 8.5999999999999999E-01 )
    global     gmHygiea
    gmHygiea = ( 5.7800000000000002E+00 )
    global     gmEunomia
    gmEunomia = ( 2.1000000000000001E+00 )
    global     gmPsyche
    gmPsyche = ( 1.8100000000000001E+00 )
    global     gmAmphitrite
    gmAmphitrite = ( 8.5999999999999999E-01 )
    global     gmEuropaAsteroid
    gmEuropaAsteroid = ( 1.5900000000000001E+00 )
    global     gmCybele
    gmCybele = ( 9.1000000000000003E-01 )
    global     gmSylvia
    gmSylvia = ( 9.8999999999999999E-01 )
    global     gmThisbe
    gmThisbe = ( 1.0200000000000000E+00 )
    global     gmEros
    gmEros = (  4.463E-4 )
    global     gmDavida
    gmDavida = ( 2.2599999999999998E+00 )
    global     gmInteramnia
    gmInteramnia = ( 2.1899999999999999E+00 )


    #planet radii.  Source pck00010.tpc

    global BODY10_RADII      
    BODY10_RADII      = (696000., 696000., 696000.  )
    global rSun      
    rSun      =  696000.
    global BODY199_RADII     
    BODY199_RADII     = ( 2439.7, 2439.7, 2439.7 )
    global rMercury     
    rMercury     = 2439.7
    global BODY299_RADII     
    BODY299_RADII     = ( 6051.8, 6051.8, 6051.8 )
    global rVenus     
    rVenus     =  6051.8
    global BODY399_RADII     
    BODY399_RADII     = ( 6378.1366, 6378.1366, 6356.7519 )
    global rEarth 
    rEarth = 6378.1366
    global BODY499_RADII       
    BODY499_RADII       = ( 3396.19, 3396.19, 3376.20 )
    global rMars      
    rMars      = 3396.19 
    global BODY599_RADII     
    BODY599_RADII     = ( 71492, 71492, 66854 )
    global rJupiter 
    rJupiter = 71492.0
    global BODY699_RADII     
    BODY699_RADII     = ( 60268,60268,54364 )
    global rSaturn
    rSaturn= 60268
    global BODY799_RADII     
    BODY799_RADII     = ( 25559, 25559, 24973 )
    global rUranus     
    rUranus     = ( 25559, 25559, 24973 )
    global BODY899_RADII     
    BODY899_RADII     = ( 24764,   24764,  24341 )
    global rNeptune     
    rNeptune     =  24764
    global BODY999_RADII     
    BODY999_RADII     = ( 1195,   1195,   1195 )
    global rPluto     
    rPluto     =  1195
    #Note Pluto's radius from NH is actually around 1187.  For reasons of matching this up with PCK10, it will remain for now, but may be changed in the future

    #holdover from old code
    global BODY503_RADII     
    BODY503_RADII     = ( 2631.2,  2631.2,    2631.2  )
    global rGany     
    rGany     = 2631.2
    
    #Planetary Oblateness, source geophysical.ker
    global j2Mercury
    j2Mercury=0
    
    global BODY399_J2 
    BODY399_J2 =    1.082616E-3
    global BODY399_J3 
    BODY399_J3 =   -2.53881E-6
    global BODY399_J4 
    BODY399_J4 =   -1.65597E-6
    global j2Earth 
    j2Earth =    1.082616E-3
    
    #Planetary orbit semi-major axes, units are in KM
    #source https://ssd.jpl.nasa.gov/txt/p_elem_t2.txt

    kmPerAU = 149597870.700 #also defined previously as global

    global smMercury 
    smMercury = 0.38709843 * kmPerAU
    global smVenus 
    smVenus = 0.72332102 * kmPerAU
    global smMars  
    smMars  = 1.52371243 * kmPerAU
    global smJupiter 
    smJupiter =  5.20248019 * kmPerAU
    global smSaturn  
    smSaturn  =  9.54149883 * kmPerAU
    global smUranus  
    smUranus  = 19.18797948 * kmPerAU
    global smNeptune 
    smNeptune = 30.06952752 * kmPerAU
    global smPluto  
    smPluto  =  39.48686035 * kmPerAU
    
def propellantConstants():
    #AMZ note: don't know where to find sources for these
    iSPHydraz = 215 #monoprollant
    iSPTAR = 295 #STAR CATLOG


def initConst():
    fundamentalUnits()
    planetaryConstants()
    propellantConstants()

#Section B: Dates and Times

#Julian day from Gregorian calendar per wikipedia/Julian_day
#Units are days (real number, including fractional part, around 2.4 million for contemporary days

def julDay(yr,mo,day,hr,minutes,seconds):
    fundamentalUnits()
    a = math.floor((14.0-mo)/12.0)
    y = yr+4800.0-a
    jd = day+math.floor((2.0+153.0*(mo+12.0*a-3.0))/5.0)+365.0*y+math.floor(y/4.0)-math.floor(y/100.0)+math.floor(y/400.0)-32045.0+(hr-12.0)/24.0+minutes/1440.0+seconds/secPerDay
    return jd

#def getET(yr,mo,day,hr,minutes,seconds):
#    #ET NAME IS BAD, NOT RETURNING ET.  NEED ZERO PADDING.
#    isocstring = str(yr)+'-'+str(mo)+'-'+str(day)+'T'+str(hr)+':'+str(minutes)+':'+str(seconds)
#    return jd
    


#Based on tests of old MMA code, output should be (and is) a real number
def centT(datein,dateformat = 'J'):
    #modified to handle other dates formats, default is usage as done in the notebook.
    if dateformat == 'S':
        #SPICE string that str2et will take  If you have a JD string, e.g. 'JD 2451545.0', use this
        #centout = (spice.str2et(datein))/36525.0/3600.0/24.0
        print('Error: SPICE not included, str2et not supported at this time. Please reformat your date.')
    elif dateformat == 'ET':
        centout = (datein)/36525.0/3600.0/24.0
    else:
        centout = (datein - 2451545.0)/36525.0
    return centout

#Modified Julian day from Gregorian calendar per wikipedia/Julian_day and http://tycho.usno.navy.mil/mjd
#Units are days (real number, including fractional part; roughly 5 digits of day number for contemporary days)
def modJulDay(yr,mo,day,hr,minutes,seconds):
    modJD = julDay(yr,mo,day,hr,minutes,seconds)-2400000.5
    return modJD

#Units are days (real number, including fractional part; roughly 5 digits of Modified Julian day number for contemporary days)
def toModJulDay(julDay):
    mjdout = julDay-2.4000005e6
    return mjdout

#Units are days (real number, including fractional part; roughly 7 digits of Julian day number for contemporary days)
def toJulDay(modJulDay):
    jdout = modJulDay+2.4000005e6
    return jdout
    
#Section C: In-Plane Geometry, Energy, Speed, Momentum

def rOfTrueAnom(smAxis,eccentri,trueAnom):
    if eccentri == 1:
        semiLatusRectum = smAxis
    else:
        semiLatusRectum = smAxis * (1-eccentri**2)
    if eccentri * np.cos(trueAnom) > -1.0:
        rOTA=semiLatusRectum / (1+eccentri * np.cos(trueAnom))
    else:
        rOTA=float('inf')
    return rOTA

def trueAnomOfR(smAxis,eccentri,radius,radialV):
    #note from Mark, put constrains on R
    tAOR = np.acos((smAxis*(1-eccentri**2)-radius)/(radius*eccentri))*np.sign(radialV)
    return tAOR

def norm(valIn):
    #shortcut for norm.  Does not do order or other tricks.
    return np.linalg.norm(valIn)

def totEnerFmRV(muCentral,radVector,velVector):
    #not sure why vectors grouped as a list in original
    return np.dot(velVector,velVector)/2-muCentral/norm(radVector) 

def totEnergy(muCentral,smAxis):
    return -muCentral/(2*smAxis)

def smaOfEnergy(muCentral,orbEnergy):
    return -muCentral/(2*orbEnergy)

def eccentri(radAp,radPeri):
    return (radApo - radPeri) / (radApo + radPeri)

def orbPer(mu,smAxis):
    return 2 * np.pi * np.sqrt(smAxis**3/mu)

def vCirc(muCentral,radCirc):
    vCirc = np.sqrt(muCentral/radCirc)
    return vCirc

def vOfRadius(muCentral,smaOrbit,currRadius):
    vOfRadius=np.sqrt(muCentral * (2 / currRadius - 1 / smaOrbit))
    return vOfRadius

def angMomentum(muCentral,smAxis,eccentri):
    #units aggre with muCentralBody and smAxis, eg km^3/s2^2 and km give km^2/s per unit of mass
    #NOTE: if eccentri = 1.0 exactly, orbit is parabolic and mAxis is assumed to give semi-latus rectum, not SMA.
    #NOTE: for hyperbolic orbit, smAxis must be negative and eccentricity must be > 1: result will still be positive.

    if eccentri == 1:
        aM=np.sqrt(muCentral*smAxis)
    else:
        aM=np.sqrt(muCentral*smAxis*(1-eccentri ** 2))
    return aM

def velOfTrueAnom(muCent,smAxis,eccentri,trueAnom):
    modAnom = ((trueAnom+np.pi) % (2*np.pi))-np.pi
    if eccentri == 1:
        radiusNow = smAxis / (1.0 + np.cos(modAnom))
        totSpeedSq = (2 * muCent) / radiusNow
    else:
        radiusNow = rOfTrueAnom(smAxis,eccentri,modAnom)
        totSpeedSq = muCent * (2 / radiusNow - 1 / smAxis)

    tangSpeed = angMomentum(muCent,smAxis,eccentri) / radiusNow

    radSpeed = np.sign(modAnom) * np.sqrt(totSpeedSq-tangSpeed**2)

    return [tangSpeed,radSpeed]

def periApseRad(muCent,smaCent,angMom):
    pAR=smaCent*(1-np.sqrt(1-angMom**2/(smaCent*muCent)))
    return pAR

def apoApseRad(muCent,smaCent,angMom):
    aAR=smaCent*(1+np.sqrt(1-angMom**2/(smaCent*muCent)))
    return aAR

def meanMotion(muCent,smaOrbit):
    mM = np.sqrt(muCent/smaOrbit)/smaOrbit
    return mM

#Section D:  "Geometry in 3D Earth-Fixed Frames"
#note, there seeem to be no test functions....

# These comments are true for all functions

# geo Long - radians of longitude with respect to the Vernal equinox (inertial frame J2000 equatorial)
# vernal Long = radians of longitude with respect to the Greenwiuch meridan (longitude 0 degrees)
# time in seconds since Greenwich Meridian((long 0 degress was aligned with Vernal equinox)

# x toward vernal equinox
# y pi/2 radians forward
# z towards north pole


def geoLongitude(vernalLong,time):
    gL = ((venalLong - time * 2 * np.pi) / 86164.0 ) % (2 * np.pi)
    return gL

def vernLongitude(geoLong,time):
    gL = ((geoLong + time * 2 * np.pi) / 86164.0 ) % (2 * np.pi)
    return gL

def xyzPos(r,lat,vernLon):
    x = r * np.cos(lat) * np.cos(vernLon)
    y = r * np.cos(lat) * np.sin(vernLon)
    z = r * np.sin(lat) 
    return [x,y,z]

def xyzPosFmGeo(r,lat,geoLong,time):
    #note input is not  vector + time
    x = r * np.cos(lat) * np.cos(vernLongitude(geoLong,time))
    y = r * np.cos(lat) * np.sin(vernLongitude(geoLong,time))
    z = r * np.sin(lat) 
    return [x,y,z]

def vernalRLatLon(xin,yin,zin):
    #note input is not a vector 
    #not sure what the deal is with the if statement.  Must be a rounding thing.  No likey.  -AMZ
    vernalR = np.sqrt(xin2+yin**2+zin**2)
    lat = np.arcsin(zin / np.sqrt(xin**2+yin**2+zin**2))
    if (xin**2 + yin**2) > 0.00001:
        lon = np.arctan(yin,xin)
    else:
        lon = 0.0
    return [vernalR,lat,lon]


def geoRLatLon(xin,yin,zin,time):
    #note input is not a vector 
    #not sure what the deal is with the if statement.  Must be a rounding thing.  No likey.  -AMZ
    vernalR = np.sqrt(xin2+yin**2+zin**2)
    lat = np.arcsin(zin / np.sqrt(xin**2+yin**2+zin**2))
    if (xin**2 + yin**2) > 0.00001:
        lontemp = np.arctan(yin,xin)
    else:
        lontemp = 0.0
        
    lon = (lontemp - (2 * np.pi * time / 86164)) % (2 * np.pi)
    return [vernalR,lat,lon]

#Section E: Relative Motion Near a Reference Circular Orbit

def normMatrixHCW(muCent,smaOrbit):
    # Matrix to normalize {R,V} to distance and rate units.
    nMHCW=np.asmatrix(np.diagflat([np.ones(3)*1.0/smaOrbit,np.ones(3)*1.0/np.sqrt(muCent/smaOrbit)]))
    return nMHCW

def stmHCW(tAnomal):
    stmH=np.asmatrix([[4.0 - 3.0*np.cos(tAnomal), 0.0,0.0, np.sin(tAnomal), 2.0*(1.0-np.cos(tAnomal)), 0.0],
          [-6.0*(tAnomal-np.sin(tAnomal)), 1.0, 0.0, 2.0*(np.cos(tAnomal)-1.0), 4.0*np.sin(tAnomal)-3.0*tAnomal, 0.0],
          [0.0,0.0,np.cos(tAnomal), 0.0, 0.0, np.sin(tAnomal)],
          [3.0*np.sin(tAnomal), 0.0,0.0,np.cos(tAnomal), 2.0*np.sin(tAnomal), 0.0],
          [-6.0*(1.0-np.cos(tAnomal)), 0.0, 0.0, -2.0*np.sin(tAnomal), 4.0*np.cos(tAnomal) -3.0, 0.0],
          [0.0,0.0, -np.sin(tAnomal), 0.0,0.0, np.cos(tAnomal) ]])
    return stmH

def finalStateHCW(muCent,smaRef,tAnomRef,initialState):
    fSHCW=np.dot(np.linalg.inv(normMatrixHCW(muCent,smaRef)),np.dot(stmHCW(tAnomRef),np.dot(normMatrixHCW(muCent,smaRef),initialState)))
    return fSHCW

def crossRanges(muCent,smaRef,tSinceEpo,rvArray):
    #still needs to throw an error if only one sat is in the array
    numSat=rvArray.shape[0]
    satDist=[]

    for i in range(numSat):
        fSHCWout=finalStateHCW(muCent,smaRef,tSinceEpo,np.transpose(np.asmatrix(rvArray[i])))
        satDist.append(np.array([fSHCWout.item(0),fSHCWout.item(1),fSHCWout.item(3)]))

        xRange=np.array(np.linalg.norm(satDist-satDist[0],axis=1))
    for i in range(1,numSat):
        xRange=np.vstack([xRange,np.array(np.linalg.norm(satDist-satDist[i],axis=1))])
    return xRange

#Section F:  Hohmann and Zero-Energy Transfers and Capture at Perigree

def dVInnerHohmann(muCent,smaInner,smaOuter):
    dvIH = vOfRadius(muCent,(smaInner+smaOuter)/2.0,smaInner)-vCirc(muCent,smaInner)
    return dvIH

def dVOuterHohmann(muCent,smaInner,smaOuter):
    dvOH = vCirc(muCent,smaOuter)-vOfRadius(muCent,(smaInner+smaOuter)/2.0,smaOuter)
    return dvOH

def dVHohmann(muCent,smaInner,smaOuter):
    return dVInnerHohmann(muCent,smaInner,smaOuter)+dVOuterHohmann(muCent,smaInner,smaOuter)

def dVZeroEnergy(muCent,smaInner,smaOuter):
    return (np.sqrt(2.0)-1.0)*(vCirc(muCent,smaInner)+vCirc(muCent,smaOuter))

def dVOrbitInsert(mu,vInf,periR,apoR):
    dVOI = np.sqrt(2.0*(((vInf**2.0)/2.0)+(mu/periR)))-vOfRadius(mu,((periR+apoR)/2.0),periR)
    return dVOI

#Section G: Bi-Elliptic In-plane transfers
def nuCPfunc(nuOffset,aone,eccone,atwo,ecctwo,periOff):
    return aone*(1.0-eccone**2.0)/(1.0+eccone*np.cos(nuOffset))-atwo*(1.0-ecctwo**2.0)/(1.0+ecctwo*np.cos(nuOffset-periOff)) 

def nuCrossPoint(aone,eccone,atwo,ecctwo,periOff):
    nCP=opt.root(nuCPfunc,np.pi/20.0,args=(aone,eccone,atwo,ecctwo,periOff))#,options = {'maxiter': 60}) #didn't like this, hm.
    return nCP

def crossAngle(mu,totEner,angMom,radNow):
    cA=np.arccos((angMom/radNow)/(vOfRadius(mu,smaOfEnergy(mu,totEner),radNow)))
    return cA

def energyDep(mu,aPlan,dV,angle):
    eD = (-mu/aPlan) + ((vCirc(mu,aPlan)+dV*np.cos(angle))**2.0 + (dV*np.sin(angle))**2.0)/2.0
    return eD

def angMomDep(mu,aPlan,dV,angle):
    aMD = (vCirc(mu,aPlan)+ dV*np.cos(angle))*aPlan
    return aMD

def vInfArriv(mu,energyDep,angMomDep,aArrive):
    ca=crossAngle(mu,energyDep,angMomDep,aArrive)
    vOR=vOfRadius(mu,smaOfEnergy(mu,energyDep),aArrive)
    vIA=norm(vOR*np.array([np.cos(ca),np.sin(ca)])-np.array([vCirc(mu,aArrive),0.0]))
    return vIA

#Section H: Oblateness Perturbations

def ascNodeRate(j2,rCentral,muCentral,smAxis,inclination,eccentri):
    aNR = -1.5*meanMotion(muCentral,smAxis)*j2*(rCentral/smAxis)**2.0*np.cos(inclination)*(1.0-eccentri**2.0)**(-2.0)
    return aNR

def periApseRate(j2,rCentral,muCentral,smAxis,inclination,eccentri):
    pAR = 0.75 * meanMotion(muCentral,smAxis)*j2*(rCentral/smAxis)**2.0*(4.0-5.0*(np.sin(inclination)**2.0))*(1.0-eccentri**2.0)**(-2.0)
    return pAR

#Section I: Planetary Ephemerides

def eltRateTable(planet):
    #this makes python users want to cry but there are no dictionaries in MMA.
    eltRateTbl={"Mercury":[[0.38709927, 0.20563593, 7.00497902, 252.2503235, 77.45779628, 48.33076593],
                           [3.7e-7, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081]], 
                "Venus":[[0.72333566, 0.00677672, 3.39467605, 181.9790995, 131.60246718,76.67984255],
                         [3.9e-6, -0.00004107, -0.0007889, 58517.81538729, 0.00268329,-0.27769418]],
                "EMB":[[1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0],
                           [5.62e-6, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0]], 
                "Mars":[[1.52371034, 0.0933941, 1.84969142, -4.55343205, -23.94362959,49.55953891], 
                        [0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343]], 
                "Jupiter":[[5.202887, 0.04838624, 1.30439695, 34.39644051, 14.72847983,100.47390909],
                           [-0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668,0.20469106]],
                "Saturn":[[9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831,113.66242448],
                          [-0.0012506, -0.00050991, 0.00193609, 1222.49362201, -0.41897216,-0.28867794]],
                "Uranus":[[19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.9542763,74.01692503],
                          [-0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281,0.04240589]],
                "Neptune":[[30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227,131.78422574],
                           [0.00026291, 0.00005105, 0.00035372, 218.45945325, -0.32241464,-0.00508664]], 
                "Pluto":[[39.48211675, 0.2488273, 17.14001206, 238.92903833, 224.06891629,110.30393684], 
                         [-0.00031596, 0.0000517, 0.00004818, 145.20780515, -0.04062942,-0.01183482]]}
    return eltRateTbl[planet]


def smaT(planet,julianD,dateformat='J'):
    eRT=((eltRateTable(planet))[0])[0]+centT(julianD,dateformat=dateformat)*((eltRateTable(planet))[1])[0]
    return eRT

def eccT(planet,julianD,dateformat='J'):
    eRT=((eltRateTable(planet))[0])[1]+centT(julianD,dateformat=dateformat)*((eltRateTable(planet))[1])[1]
    return eRT


def inclinT(planet,julianD,dateformat='J'):
    eRT=((eltRateTable(planet))[0])[2]+centT(julianD,dateformat=dateformat)*((eltRateTable(planet))[1])[2]
    return eRT

def mnLnT(planet,julianD,dateformat='J'):
    eRT=((eltRateTable(planet))[0])[3]+centT(julianD,dateformat=dateformat)*((eltRateTable(planet))[1])[3]
    return eRT

def lnPeriT(planet,julianD,dateformat='J'):
    eRT=((eltRateTable(planet))[0])[4]+centT(julianD,dateformat=dateformat)*((eltRateTable(planet))[1])[4]
    return eRT

def lnAscNT(planet,julianD,dateformat='J'):
    eRT=((eltRateTable(planet))[0])[5]+centT(julianD,dateformat=dateformat)*((eltRateTable(planet))[1])[5]
    return eRT

def argPeriHel(planet,julianD,dateformat='J'):
    return lnPeriT(planet,julianD,dateformat=dateformat)-lnAscNT(planet,julianD,dateformat=dateformat)

def meanAnom(planet,julianD,dateformat='J'):
    return ((mnLnT(planet,julianD,dateformat=dateformat)-lnPeriT(planet,julianD,dateformat=dateformat)+180.0) % 360.0)-180.0

def eccAnomVelPE(planet,julianD,dateformat='J'):
    eccen = eccT(planet, julianD,dateformat=dateformat)
    meanAn = meanAnom(planet,julianD,dateformat=dateformat)*np.pi/180.0
    eccAn = meanAn + eccen * np.sin(meanAn)
    delE =1.0
    iterCount = 0
    while np.absolute(delE) > (1e-6*np.pi/180.0) and iterCount < 100:
        delE = (meanAn - (eccAn -  eccen* np.sin(eccAn)))/(1.0 - eccen*np.cos(eccAn))
        eccAn = eccAn + delE
        iterCount = iterCount + 1
        if iterCount > 90:
            print('eccAnomVelPE itercount, eccAn:',iterCount,', ',eccAn)
    eccAnDeg = eccAn*180.0/np.pi
    eccAnDEdt = ((eltRateTable(planet))[1])[3]/(1.0-eccen*np.cos(eccAn))
    return [eccAnDeg,eccAnDEdt]

def xyVelPQ(planet,julianD,dateformat='J'):
    eccen = eccT(planet, julianD,dateformat=dateformat)
    eccenAnomNRate = eccAnomVelPE(planet,julianD,dateformat=dateformat)
    eANR0=eccenAnomNRate[0]
    eANR1=eccenAnomNRate[1]
    xyVPQ=smaT(planet,julianD,dateformat=dateformat)*np.array([[np.cos(eANR0*np.pi/180.0)-eccen,
                                                                np.sin(eANR0*np.pi/180.0)*np.sqrt(1.0-eccen**2.0)],
                                                               [-np.sin(eANR0*np.pi/180.0)*eANR1*np.pi/180.0,
                                                                 np.cos(eANR0*np.pi/180.0)*np.sqrt(1.0-eccen**2.0)*eANR1*np.pi/180.0]],
                                                              dtype='float64')
    return xyVPQ


def posVelEcl(planet,julianD,dateformat='J'):
    colw = np.cos(argPeriHel(planet,julianD,dateformat=dateformat)*np.pi/180.0)
    silw = np.sin(argPeriHel(planet,julianD,dateformat=dateformat)*np.pi/180.0)
    coOm = np.cos(lnAscNT(planet,julianD,dateformat=dateformat)*np.pi/180.0)
    siOm = np.sin(lnAscNT(planet,julianD,dateformat=dateformat)*np.pi/180.0)
    coIn = np.cos(inclinT(planet,julianD,dateformat=dateformat)*np.pi/180.0)
    siIn = np.sin(inclinT(planet,julianD,dateformat=dateformat)*np.pi/180.0)
    xyVxVyPq = xyVelPQ(planet,julianD,dateformat=dateformat)
    xyVxVyPq1 = np.transpose(np.asmatrix(xyVxVyPq[0]))
    xyVxVyPq2 = np.transpose(np.asmatrix(xyVxVyPq[1]))
    
    rotMat = np.asmatrix([[(colw * coOm - silw * siOm * coIn ), (-silw * coOm - colw * siOm * coIn)],
              [(colw * siOm + silw * coOm * coIn ),(-silw * siOm + colw * coOm * coIn) ],
              [(silw * siIn), (colw * siIn)]])
    pV1=np.array(np.transpose(np.dot(rotMat,xyVxVyPq1)))
    pV2=np.array(np.transpose(np.dot(rotMat,xyVxVyPq2)))
    return [pV1,pV2]

#Section J: Gauss Transfers

def sSFunc(zZ):
    if zZ > 1.0:
         sSout = (np.sqrt(zZ)-np.sin(np.sqrt(zZ))) / np.sqrt(zZ**3.0)
    elif zZ < -1.0:
        sSout = (np.sinh(np.sqrt(-1.0*zZ))-np.sqrt(-1.0*zZ))/np.sqrt((-1.0*zZ)**3.0)
    else:
        sSout = (1.0/math.factorial(3) - (1.0*zZ)/math.factorial(5) + (zZ**2.0)/math.factorial(7) 
                 - (zZ**3.0)/math.factorial(9) + (zZ**4.0)/math.factorial(11) - 
                 (zZ**5.0)/math.factorial(13) + (zZ**6.0)/math.factorial(15) - (zZ**7.0)/math.factorial(17))

    return sSout


def cCFunc(zZ):
    if zZ > 1.0:
        cCout = (1.0-np.cos(np.sqrt(zZ)))/zZ
    elif zZ < -1.0:
        cCout = (1.0-np.cosh(np.sqrt(-zZ)))/zZ
    else:
        cCout = (1.0/math.factorial(2) - 1.0*zZ/math.factorial(4) + zZ**2.0/math.factorial(6) -
                 zZ**3.0/math.factorial(8) + zZ**4.0/math.factorial(10) - zZ**5.0/math.factorial(12)
                 + zZ**6.0/math.factorial(14) - zZ**7.0/math.factorial(16))
    return cCout

def sSPrime(zZ):
    if np.absolute(zZ) > 9.0:
        sSP = (cCFunc(zZ) - 3.0 * sSFunc(zZ))/(2.0*zZ)
    else:
        sSP = (-1.0/math.factorial(5)+2.0*zZ/math.factorial(7)-3.0*(zZ**2.0/math.factorial(9))
                +4.0*(zZ**3.0/math.factorial(11))-5.0*(zZ**4.0/math.factorial(13))
                +6.0*(zZ**5.0/math.factorial(15))-7.0*(zZ**6.0/math.factorial(17))
                +8.0*(zZ**7.0/math.factorial(19))-9.0*(zZ**8.0/math.factorial(21))
                +10.0*(zZ**9.0/math.factorial(23))-11.0*(zZ**10.0/math.factorial(25))
                +12.0*(zZ**11.0/math.factorial(27))-13.0*(zZ**12.0/math.factorial(29)))
    return sSP

def cCPrime(zZ):
    if np.absolute(zZ) > 5.0:
        cCP = (1.0 - zZ*sSFunc(zZ) - 2.0*cCFunc(zZ))/(2.0*zZ)
    else:
        cCP = (-1.0/math.factorial(4)
                +2.0*zZ**1.0/math.factorial(6) -3.0*zZ**2.0/math.factorial(8)
                +4.0*zZ**3.0/math.factorial(10)-5.0*zZ**4.0/math.factorial(12)
                +6.0*zZ**5.0/math.factorial(14)-7.0*zZ**6.0/math.factorial(16)
                +8.0*zZ**7.0/math.factorial(18)-9.0*zZ**8.0/math.factorial(20)
                +10.0*zZ**9.0/math.factorial(22)-11.0*zZ**10.0/math.factorial(24))
    return cCP




def zLimCalcFunc(zLimCalcX,rLnch,rArrv,aConst):
    return aConst * (1 - zLimCalcX * sSFunc(zLimCalcX) ) / np.sqrt(cCFunc(zLimCalcX)) - (rLnch + rArrv) 



def gaussTGuessMin(delMeanAnSq,rDepart,rArrive,delNu,aConst,mu,debugnow):
    return (gaussTGuess(delMeanAnSq,rDepart,rArrive,delNu,aConst,mu,debugnow))[0]

def gaussTGuess(delMeanAnSq,rDepart,rArrive,delNu,aConst,mu,debugnow):
    sS = sSFunc(delMeanAnSq)
    cC = cCFunc(delMeanAnSq)
    yY = np.max([0.0,rDepart + rArrive - aConst *((1.0-delMeanAnSq * sS)/np.sqrt(cC))])

    if debugnow == 1:
        print('gaussTGuess in: delMeanAnSq: {:}, rDepart: {:}, rArrive: {:}, delNu: {:}, aConst: {:}, mu: {:}'.format(delMeanAnSq,rDepart,rArrive,delNu,aConst,mu))

    '''
    if cC < 0 or yY < 0:
        print('zZ: ',delMeanAnSq,'cC: ',cC,' yY: ',yY,' sS: ',sS,' aconst: ',aConst)
        if cC < 14:
            return
        if aConst > 6.0e7:
            return
    '''

    xX = np.sqrt(yY / cC)
    gtgout = [((xX**3.0) * sS + aConst * np.sqrt(yY)) / np.sqrt(mu),
              ((xX**3.0) * (sSPrime(delMeanAnSq) - 3.0 * sS * cCPrime(delMeanAnSq) / ( 2.0 * cC )) + 
               (aConst / 8.0) * (3.0 * sS * np.sqrt(yY) / cC + aConst / xX) ) / np.sqrt(mu),yY]
    return gtgout

def gaussTrajec(mu,rVecLnch,rVecArrv,tXfer,goPrograde,numLoops):
    # note that this function has a lot of comments in the original mathematica about test cases.
    #print("gaussTraject: rLaunch, rArriv, tXfer: ", rVecLnch, rVecArrv, tXfer) 
    
    # From BMW Section 5.3, page 234, Universal variables formulation.
    # mu: Mass times gravitation constant for the center of attraction (normally km^3/s^2) 
    # rVecLnch: Position vector for start of trajectory; units agree with mu (normally km) 
    # rVecArrv: Position vector for end of trajectory; units agree with mu (normally km)
    # tXfer: Time between leaving  rVecLnch and arriving at rVecArrv; units agree with mu (normally s)
    # goPrograde: True or False flag which should force right-hand orbit around system Z vector, which is minimum energy for most solar-system orbits 
    # numLoops: Number of phasing loops (complete orbits) between start and end of trajectory; normally 0
    # Output is two 3-vectors of velocity  relative to center of system. 
    #     First vector is velocity departing rVecLnch. Units agree with mu, rVecLnch, rVecArriv, and tXfer on input (normally km/s)
    #     Second vector is velocity arriving at rVecArrv. Units agree with mu, rVecLnch, rVecArriv, and tXfer on input (normally km/s) 
    # NOTES:
    #     If tXfer is specified as too short for a multi-phasing-orbit transfer, gaussTraject will print a warning and return with Null results. This is checked using FindMinimum, which is sometimes flaky.
    #     In some cases Newton search will project a value for zZ (square of mean anomaly change) which is too negative, and would lead to imaginary results. The If[(zZ < zLimLow).... block prevents this, by placing the new guess for zZ 20% of the way from the lower limit to the current guess for zZ. 
    #     First If[] block uses goPrograde to decide whether to select the trajectory with initial V closest to departing planet V, suitable for low-energy trajectories (as opposed to low-time trajectories).
    #     zZ is limited to stay within the correct (2 Pi)^2 interval to match the desired number of Loops. 
    #     "deBug" T/F parameter if set to True will print out input and intermediate values during execution to allow insight into calculations. Normally leave this set to False.

    machinePrecision = 16.0
    zGuard = 1.0*10**(11.0-machinePrecision) # not sure what this is.    Precision for float 64 is about 2e-16,-1e-16
    rLnch = np.linalg.norm(rVecLnch)
    rArrv = np.linalg.norm(rVecArrv)
    rVecLnch = rVecLnch
    rVecArrv = rVecArrv
    iterCount = 0

    veccross = np.cross(rVecLnch,rVecArrv)

    if goPrograde and veccross[2] < 0.0:
        #print("delNu forced to go long way to go Prograde")
        delNu = 2.0*np.pi * (1.0+numLoops) - np.arccos(np.dot(rVecLnch,rVecArrv)/rLnch/rArrv)
    else:
        delNu = 2.0*np.pi * (numLoops) + np.arccos(np.dot(rVecLnch,rVecArrv)/rLnch/rArrv)

    #print('delNu: ',delNu)

    aConst = np.sign(np.sin(delNu)) * np.sqrt(rLnch*rArrv*(1.0+np.cos(delNu)))

    if numLoops > 0:
        zLimLow = (2.0*np.pi*numLoops)**2.0+zGuard
    else:
        zLimLow = -8.0*np.pi**2

    zLimHigh = (((numLoops+1.0)*2*np.pi)**2.0) - zGuard

    #these three lines save hours and stop the program from breaking.
    if delNu < np.pi / 2.0:
        zLimCalc = opt.fsolve(zLimCalcFunc,[0.00001],args =(rLnch,rArrv,aConst))
        zLimLow = zLimCalc[0]

    zZ = delNu ** 2.0

    precisionGoal = tXfer * 10.0**(5 - machinePrecision)

    start1=time.time()

    tXtest = []

    npts = 100000

    zsearch =  np.linspace(zLimLow,zLimHigh,npts)

    debugnow = 0

    # we should optimimize, but this function can't take an array (yet)
    '''
    for x in range(len(zsearch)):
        gtgout = gaussTGuess(zsearch[x],rLnch,rArrv,delNu,aConst,mu,debugnow)
        tXtest.append(gtgout[0])
       
    mingtg = np.argmin(tXtest)
    
    end1=time.time()



    ax2=plt.subplot(111)
    ax2.plot(zsearch,tXtest)
    plt.show()


    mintXfer = tXtest[mingtg]
    '''    
    

    #print('From brute force searching',mingtg,zsearch[mingtg],mintXfer)
    #print("--- %s seconds ---" % (end1-start1))

    start2=time.time()

    zsearch =  np.linspace(zLimLow,zLimHigh,npts)

    #gtgout = gaussTGuess(zsearch[0],rLnch,rArrv,delNu,aConst,mu,debugnow)
    #print(type(gtgout),gtgout)

    mintXfer2 = opt.minimize_scalar(lambda zsea: gaussTGuessMin(zsea,rLnch,rArrv,delNu,aConst,mu,debugnow),bounds = (zLimLow,zLimHigh),method='bounded')
    #,args=([rLnch,rArrv,delNu,aConst,mu,debugnow])
    #tol = precisionGoal
    end2=time.time()
    #print('From actual opt fn',mintXfer2)
    #print("--- %s seconds ---" % (end2-start2))


    mintXfer = mintXfer2.fun

    if numLoops > 0 and tXfer < mintXfer:
        print('WARNING: in gaussTraj, specified transfer time too low for multi-orbit solution')
        print('   Launch Position: ',rVecLnch)
        print('   Arrival Position: ',rVecArrv)
        print('   Transfer Time: ',tXfer)
        return None

    tFunctions = np.array([2.0*tXfer,1.0,1.0])

    while np.absolute(tFunctions[0] - tXfer) > precisionGoal and iterCount < 100:
        iterCount = iterCount + 1 
        debugnow =0
        tFunctions = gaussTGuess(zZ,rLnch,rArrv,delNu,aConst,mu,debugnow)
        #print("pre-update zz: ",zZ,", tFunctions[0]: ",tFunctions[0],", tXfer: ",tXfer," , tFunctions[1]: ",tFunctions[1])
        zZold = zZ;
        zZ = zZold - (tFunctions[0] - tXfer)/tFunctions[1]
        #print("Post-update pre-limit zz ", zZ,", tdiff:",(tFunctions[0] - tXfer))
        if zZ < zLimLow:
            zZ = zLimLow + (zZold - zLimLow) * 0.2
        #zZ = np.max([zZ,zLimLow])
        zZ = np.min([zZ,zLimHigh])

    gG = aConst*np.sqrt(tFunctions[2]/mu)

    #pdb.set_trace()
    outvec = [(rVecArrv - (1.0 - tFunctions[2]/rLnch)*rVecLnch)/gG,((1 - tFunctions[2]/rArrv)*rVecArrv - rVecLnch)/gG]
    return outvec


def gaussPlanetary(fromPlanet,fromJD,toPlanet,toJD,phasingOrbits):
# Initial and Final Velocity vectors to go from R1 to R2 in given time
# Input quantities and units:
#     fromPlanet: Planet to depart from. String variable agreeing with Planets listed in Part I Planetary Ephemerides. 
#         Those routines are used to derive planet position from Julian Date of departure.
#     fromJD: Julian Date of departure. May be generated by julDay[yr_,mo_,day_,hr_,minutes_,seconds_] (see Part B Dates and Times) 
#     toPlanet: Planet to arrive at. String variable, as above.
#     toJD:  Julian Date of arrival, as above.
#     phasingOrbits: number of complete orbits (2Pi radians of true anomaly) between departure and arrival. Normally 0.

#     First vector is velocity departing fromPlanet, planetocentric (km/s). 
#     Second vector is velocity arriving at toPlanet, planetocentric (km/s).
# NOTES: 
#     To generate heliocentric velocity it is necessary to add planet's heliocentric velocity to output velocity vectors; 
#     use posVelEcl[planet_, julianD_] from Section I to generate that velocity and be sure to convert from AU/Century to km/s.
#     gaussPlanetary automatically assumes prograde choice is preferable.

    gmSun  = 132712440017.987 #find a better way for this, in case constant is different
    kmPerAU = 149597870.700
    
    pvFrom = posVelEcl(fromPlanet,fromJD)
    pvTo = posVelEcl(toPlanet,toJD)

    gP1 = np.array(gaussTrajec(gmSun,pvFrom[0].flatten()*kmPerAU,pvTo[0].flatten()*kmPerAU,(toJD-fromJD)*secPerDay,1,phasingOrbits))
    gP2 = np.array([pvFrom[1].flatten()*kmPerAU/secPerCy,pvTo[1].flatten()*kmPerAU/secPerCy])

#    print('gP1',gP1)

    gP = gP1-gP2
    return gP


#Section K: General Position and Velocity
# rnvFmElts[mu_, orbElts_, tNow_] - Given Keplerian elements, generate R and V vectors as a function of time.
# Inputs and suggested units (note other units work but must be consistent across all input parameters; output will appear in same units): 
#    Mass times gravitation constants for the center of attraction (km^3/s^2)
#    Keplerian orbit elements
#    {
#       a   semimajor axis , (km)  NOTE: this argument is assumed to be p, semi-latus rectum, if e = 1.0, denoting a parabolic orbit)
#       e   eccentricity, (NOTE: e=1.0 denotes a parabolic orbit and changes interpretation of a)
#             NOTE also; if a < 0, e must be > 1, denoting hyperbolic orbit.
#       i   inclination (wrt equator) (radians),
#       Om  longitude of asc. node (rel. to vernal equinox) (radians),
#       w   argument of perigee (wrt ascending node) (radians),
#       tPaps   time of periapsis passage (s)
#    }
#    tNow    current time (seconds)
# Result: 
#    {
#       {r} Position vector at time tNow in IJK (vernal equinox, equator, north polar axis) (km)
#       {v}     Velocity vector at time tNow in IJK (vernal equinox, equator, north polar axis) (km/s)
#    }
def rnvFmElts(mu,orbElts,tNow):
    machinePrecision = 16.0
    a = orbElts[0]
    e = orbElts[1]
    i = orbElts[2]
    Om = orbElts[3]
    w = orbElts[4]
    tPaps = orbElts[5]
    iterCount = 0

    if (a > 0.0):
        tThisOrbit = (tNow-tPaps) % orbPer(mu,a)
    else:
        tThisOrbit = tNow-tPaps #hyperbolic case


    if (e < 1.0):
        meanAn = tThisOrbit * np.sqrt(mu/(a**3.0))
        eccAnom = meanAn + e * np.sin(meanAn)
        
        while (np.absolute(meanAn - (eccAnom - e * np.sin(eccAnom))) > (10**(5.0-machinePrecision))) and iterCount < 100:
            iterCount = iterCount + 1
            delE = (meanAn - (eccAnom - e *np.sin(eccAnom))) / (1.0 - e*np.cos(eccAnom))
            eccAnom = eccAnom + delE
        eccAnom = ((eccAnom +np.pi) % (2.0*np.pi))-np.pi
        nu = np.arccos((e-np.cos(eccAnom))/(e*np.cos(eccAnom)-1.0))*np.sign(eccAnom)
        p = a*(1.0-e**2.0) 
    elif (e == 1.0):
        meanAn = tThisOrbit *2.0 * np.sqrt(mu)
        eccAnom = (3.0*meanAn*(1.0-a))**(1.0/3.0)
        while(np.absolute(meanAn - (a*eccAnom+eccAnom**3.0/3.0)) >  (10**(5.0-machinePrecision))) and iterCount < 100:
            iterCount = iterCount + 1
            delE = (meanAn - (a*eccAnom+eccAnom**3.0/3.0))/ (eccAnom**2.0 + a);
            eccAnom = eccAnom + delE
        nu = 2*np.arctan(eccAnom/np.sqrt(a))
        p = a # Added by AMZ, was missing before.  But this work won't in equations!
    else:  # e > 1.0 -> hyperbolic
        meanAn = tThisOrbit*np.sqrt(mu/(-a**3.0))
        eccAnom = -1.0*(meanAn-e*np.sinh(meanAn))
        while(np.absolute(meanAn - (e*np.sinh(eccAnom) - eccAnom)) > (10**(5.0-machinePrecision))) and iterCount < 100:
            iterCount = iterCount + 1
            delE = (meanAn - (e*np.sinh(eccAnom)-eccAnom))/(e*np.cosh(eccAnom)-1.0)
            if np.absolute(delE) > np.pi/2.0:
                delE = np.pi/2.0 * (delE / np.linalg.norm(delE)) 
                # (* Ad-hoc limit to prevent delE runaway *)
            eccAnom = eccAnom + delE
        nu = np.arccos((e-np.cosh(eccAnom))/(e*np.cosh(eccAnom)-1.0))*np.sign(eccAnom)
        p = a*(1.0-e**2.0)

    r = p / (1.0+e*np.cos(nu))

    posVecInPQ = np.asmatrix([[r*np.cos(nu)],[r*np.sin(nu)],[0.0]]) 
    velVecInPQ = np.sqrt(mu/p)*np.asmatrix([[-1.0*np.sin(nu)],[e+np.cos(nu)],[0.0]]) #e = 1 case will not work

    RotMatPQijk = np.asmatrix([[np.cos(Om)*np.cos(w) - np.sin(Om)*np.sin(w)*np.cos(i) , 
                                - np.cos(Om)*np.sin(w) - np.sin(Om)*np.cos(w)*np.cos(i), 
                                np.sin(Om)*np.sin(i) ], 
                               [np.sin(Om)*np.cos(w) + np.cos(Om)*np.sin(w)*np.cos(i), 
                                - np.sin(Om)*np.sin(w) + np.cos(Om)*np.cos(w)*np.cos(i), 
                                - np.cos(Om)*np.sin(i) ], 
                               [np.sin(w)*np.sin(i), np.cos(w)*np.sin(i),np.cos(i) ] ])

    #print(RotMatPQijk)
    #print(posVecInPQ)
    #print(velVecInPQ)

    posrot = np.matmul(RotMatPQijk,posVecInPQ)
    velrot = np.matmul(RotMatPQijk,velVecInPQ)

    posVelInIJK = [posrot,velrot]

    return posVelInIJK


def eltsFmRV(muCent,rrVec,vvVec,tEpo):
    # Input Units:
    # muCent:     km^3/s^2 for mass of central body times gravitational constant (or other suitable units agreeing with rrVec, vvVec, and t) 
    # rrVec:      position vector at epoch time in km (or matching units)
    # vvVec:      Velocity vector at epoch time in km/s (or matching units)
    # tEpo:       Epoch time in seconds (or matching units)

    #Note it is very important that rrVec and vvVec come in as np.arrays of 3 elments long.
    #print('rrVec',rrVec,'vvVec',vvVec)
    hhVec=np.cross(rrVec,vvVec)

    #print('hhVec',hhVec)
    nnVec=np.cross(np.array([0,0,1]),hhVec)
    #print('nnVec',nnVec)
    rscal = np.linalg.norm(rrVec)

    eccVec=(1.0/muCent)*((np.dot(vvVec,vvVec)-muCent/rscal)*rrVec-np.dot(rrVec,vvVec)*vvVec)
    #print('eccVec',eccVec)

    #print('test here: ',np.dot(rrVec,vvVec))
    ecc = np.linalg.norm(eccVec)
    if (np.dot(rrVec,vvVec) > 0.0):
        trueAnom = np.arccos(np.dot(eccVec,rrVec)/ecc/rscal)

    else:
        trueAnom = 2.0*np.pi-np.arccos(np.dot(eccVec,rrVec)/ecc/rscal)

    if ecc < 1.0:
        if np.dot(rrVec,vvVec) > 0.0:
            eccAnom = np.arccos((ecc+np.cos(trueAnom))/(1.0+ecc*np.cos(trueAnom)))
        else:
            eccAnom = -1.0*np.arccos((ecc+np.cos(trueAnom))/(1.0+ecc*np.cos(trueAnom)))
    else:
        #print('True Anom',trueAnom)
        if trueAnom < np.pi:
            eccAnom = np.arccosh((ecc+np.cos(trueAnom))/(1.0+ecc*np.cos(trueAnom)))
        else:
            eccAnom = -1.0*np.arccosh((ecc+np.cos(trueAnom))/(1.0+ecc*np.cos(trueAnom)))

    sma = (np.dot(hhVec,hhVec)/muCent)/(1.0-(np.dot(eccVec,eccVec)))


    incl=np.arccos(np.dot(hhVec,np.array([0.0,0.0,1.0]))/np.linalg.norm(hhVec))
   

 
    if np.dot(nnVec,np.array([0.0,1.0,0.0])) > 0.0:
        Om = np.arccos(np.dot(nnVec,np.array([1.0,0.0,0.0]))/np.linalg.norm(nnVec))
    else:
        Om = -1.0 * np.arccos(np.dot(nnVec,np.array([1.0,0.0,0.0]))/np.linalg.norm(nnVec))

    if np.dot(eccVec,np.array([0.0,0.0,1.0])) > 0.0:
        w = np.arccos(np.dot(nnVec,eccVec)/np.linalg.norm(nnVec)/np.linalg.norm(eccVec))
    else:
        w = -1.0 * np.arccos(np.dot(nnVec,eccVec)/np.linalg.norm(nnVec)/np.linalg.norm(eccVec))
        
    if (ecc < 1.0):
        tPaps = tEpo - np.sqrt(sma**3.0/muCent)*(eccAnom-ecc*np.sin(eccAnom))
    else:
        tPaps = tEpo - np.sqrt(-1.0*sma**3.0/muCent)*(ecc*np.sinh(eccAnom)-eccAnom)
    
    return [sma,ecc,incl,Om,w,tPaps]


#note testComparMat has been moved to orbittoolsk.py because it was a testing function

def tDiffKepler(muCent,alpha,x,rrVec,vvVec,r0,tNow,tEpo):
    tDKout = ((cCFunc(alpha*x**2.0)*x**2.0*np.dot(rrVec,vvVec)/np.sqrt(muCent) +
               (1.0-r0*alpha)*x**3.0*sSFunc(alpha*x**2.0)+
               r0*x)/np.sqrt(muCent)-(tNow-tEpo))
    return tDKout

def dTdXKepler(muCent,alpha,x,rrVec,vvVec,r0):
    dTdXKout = (x**2.0*cCFunc(alpha*x**2.0)+
                np.dot(rrVec,vvVec)*x*(1.0-(alpha*x**2.0)*sSFunc(alpha*x**2.0))/np.sqrt(muCent)+
                r0*(1.0-(alpha*x**2.0)*cCFunc(alpha*x**2.0)))/np.sqrt(muCent)
    return dTdXKout


def rvKeplerOfT(muCent,rrVec,vvVec,tEpo,tNow):
    machinePrecision = 16.0
    debug = 0
    iterCount = 0
    #Position and velocity as a function of rrVec, vvVect and T
    r0 = np.linalg.norm(rrVec)
    alpha = (2.0*muCent/r0-np.dot(vvVec,vvVec))/muCent
    if debug == 1:
        print('rvK Entry: r0: ',r0,' v0:',np.linalg.norm(vvVec))
        print('r.v, muCent, alpha: ',np.dot(rrVec,vvVec),' muCent: ',muCent,' alpha: ',alpha)

    perEllips = (2 * np.pi/np.sqrt(muCent*alpha**3.0))

    if alpha > 0.0 and np.absolute(tEpo - tNow) > (perEllips/2.0):
        tCalc = tEpo + (((tNow-tEpo)+perEllips/2.0) % perEllips) - perEllips/2.0
        if debug == 1:
            print('rvKeptlerOfT: Changed tNow to less than one half orbit from tEpo.')
            print('  Old tNow = ',tNow,' New tNow = ', tCalc)
    else:
        tCalc = copy.deepcopy(tNow)

    if tCalc == tEpo:
        x = 1.0
    elif alpha > 0.0:
        x = np.sqrt(muCent)*alpha*(tCalc-tEpo)
    elif alpha == 0.0:
        x = np.sqrt(muCent)*(tCalc-tEpo)/np.linalg.norm(rrVec)
    elif alpha < 0.0:
        xRadic = np.log((-2.0*muCent*alpha*(tCalc-tEpo))/(np.dot(rrVec,vvVec)+np.sign(tCalc-tEpo)*np.sqrt(-1.0*muCent/alpha)*(1.0-r0*alpha)))
        if debug == 1:
            print('xRadic: ',xRadic)
        if xRadic > 0.5:
            x = np.sign(tCalc-tEpo)*np.sqrt(-1.0/alpha)*xRadic
            if debug ==1:
                print('  Hyperbolic large angle xguess: ',x)
        else:
            x = np.sqrt(np.linalg.norm(rrVec))*np.linalg.norm(vvVec)*(tCalc-tEpo)/np.linalg.norm(rrVec)
            if debug ==1:
                print('  Hyperbolic small angle xguess: ',x)
    else:
        x = np.sqrt(np.linalg.norm(rrVec))*np.linalg.norm(vvVec)*(tCalc-tEpo)/np.linalg.norm(rrVec)

    if debug == 1:
        defX = np.sqrt(np.linalg.norm(rrVec))*np.linalg.norm(vvVec)*(tCalc-tEpo)/np.linalg.norm(rrVec)
        print('Default x: ',defX)
        
    if x == 0:
        x = 1 #don't get this.
        print('x made to be 1 because it was zero')

    maxPrecis = 10.0**(4.0+np.ceil(np.log10(np.absolute(x*dTdXKepler(muCent,alpha,x,rrVec,vvVec,r0))))-machinePrecision)

    tDiff = tDiffKepler(muCent,alpha,x,rrVec,vvVec,r0,tCalc,tEpo)

    if debug == 1:
        print('tCalc-tEpo: ',tCalc-tEpo)
        print('initial guess for x: ',x)
        print('max precision: ',maxPrecis)
        print('tDiff init: ',tDiff)



    while (np.absolute(tDiff) > maxPrecis) and iterCount < 100:
        iterCount = iterCount + 1
        dtdxkep = dTdXKepler(muCent,alpha,x,rrVec,vvVec,r0)
        x = x - tDiff / dtdxkep
        tDiff = tDiffKepler(muCent,alpha,x,rrVec,vvVec,r0,tCalc,tEpo)
        if debug == 1:
            print('x, tDiff, slope: ',x,' , ',tDiff,' , ',dtdxkep)
        
    if debug == 1:
        print('Final: x, tDiff, slope: ',x,' , ',tDiff,' , ',dtdxkep)
        
    fF = 1.0 - x**2.0 * cCFunc(alpha*x**2.0)/r0
    gG = (tCalc - tEpo)  - x**3.0*sSFunc(alpha*x**2.0)/np.sqrt(muCent)
    if debug == 1:
        print('fF ',fF,'gG ',gG)
        print('rrVec ',rrVec,'vvVec ',vvVec)
    rVecNow = fF * rrVec + gG * vvVec
    rrNow = np.linalg.norm(rVecNow)
    fFDot = np.sqrt(muCent)*x*(alpha*x**2.0*sSFunc(alpha*x**2.0)-1.0)/(r0*rrNow)
    gGDot = 1.0 - x**2.0*cCFunc(alpha*x**2.0)/rrNow
    vvVecNow = fFDot*rrVec + gGDot*vvVec
    return [rVecNow,vvVecNow]
                      
#Section L: Delta V and Rocket Equation
def dVRocket(iSp,mRatio):
    dvr = iSp *0.00981*np.log(mRatio)
    return dvr

def massRatioRocket(deltaV,iSp):
    mrr = np.exp(deltaV/iSp/0.00981)
    return mrr

#Section M: Swingby Calculations

def missDist(muPlanet,vInf,turnDelta):
    md = (muPlanet/vInf**2.0)*((1.0/np.sin(turnDelta/2.0))-1.0)
    return md

def planetSwingby(origPlanet,origJD,swingPlanet,swingJD,tgtPlanet,tgtJD):
    #note use of gauss transfers function
    velIn = (gaussPlanetary(origPlanet,origJD,swingPlanet,swingJD,0.))[1]
    velOut = (gaussPlanetary(swingPlanet,swingJD,tgtPlanet,tgtJD,0.))[0]
    print(velIn)
    print(velOut)
    ps1 = np.arccos(np.dot(velIn,velOut)/np.linalg.norm(velIn)/np.linalg.norm(velOut))
    ps2 = np.linalg.norm(velIn)-np.linalg.norm(velOut)
    return [ps1,ps2]


def periApseOfB(muPlanet,vInf,bVec):
    paob = (-1.0*muPlanet+np.sqrt(muPlanet**2.0 + np.dot(bVec,bVec)*vInf**4.0))/vInf**2.0
    return paob

def bMagOfPeriApse(muPlanet,vInf,rPeri):
    bmofpa = np.sqrt(((rPeri*vInf**2.0+muPlanet)**2.0-muPlanet**2.0)/vInf**4.0)
    return bmofpa


def turnAngleOfPeriApse(muPlanet,vInf,rPeri):
    taopa = 2.0 * np.arcsin(1.0/(1.0+rPeri*(vInf ** 2.0/muPlanet)))
    return taopa


def impactParameter(muPlanet,vInf,radPlanet):
    ip = (radPlanet/vInf)*np.sqrt(vInf**2.0 + 2 * muPlanet/radPlanet)
    return ip


def rPeriApseOfTurnAngle(muPlanet,vInf,turnAngle):
    rpaota = muPlanet/vInf**2.0*(1.0/np.sin(turnAngle/2.0)-1.0)
    return rpaota
                                 
                                 
def vOutSwingby(muPlanet,vInfIn,bVec):
    vMag = np.linalg.norm(vInfIn)
    deflAng = turnAngleOfPeriApse(muPlanet,vMag,periApseOfB(muPlanet,vMag,bVec))
    ttVec = np.cross(vInfIn,np.array([0.,0.,1.]))/np.linalg.norm(np.cross(vInfIn,np.array([0.,0.,1.])))
    rrVec = np.cross(vInfIn,ttVec)/np.linalg.norm(np.cross(vInfIn,ttVec))
    bUnitInIJK = (bVec[0]*rrVec+bVec[1]*ttVec)/np.linalg.norm(bVec)
    vos = vInfIn * np.cos(deflAng)-vMag*bUnitInIJK * np.sin(deflAng)
    return vos




def main():
    pass

if __name__ == "__main__":
    main()

