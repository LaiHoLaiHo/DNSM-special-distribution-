"""xyz here is define by the orbital plane
x:   M2 -> M1 direction
y: orbital velocity direction orbital plane 
z: perpendicular to orbital plane 
according to 
ORBITAL CHARACTERISTICS OF BINARY SYSTEMS AFTER ASYMMETRIC 
SUPERNOVA EXPLOSIONS VASSILIKI KALOGERA """


import math as m
from scipy.integrate import quad
c = 299792458  #meter
G = 6.674*10**-11 #m3 /kg /s2
SM = 1.989*10**30 #kg  #solar mass
SR = 6.957*10**8 #meter #solar radii
from scipy.integrate import quad
def FVbr(MHe,M2,Ai):  #relative velo binary
    Vbr = (G*(MHe+M2)/Ai)**0.5
    return Vbr
def FPDMv(casi,v): #Probability distribution of velo megnitude
    P = (2/m.pi)**0.5/casi**3*v**2*m.exp(-v**2/2/casi**2)
    return P
def FVkx(Vk,theta, phi):
    #this is not the system velovity direction
    #it is the orbital xyz
    #20200720 found it is wrong about the direction
    #Vky is the vr dir
    Vkx = Vk*m.sin(theta)*m.cos(phi)
    return Vkx
def FVky(Vk,theta, phi):
    #this is not the system velovity direction
    #it is the orbital xyz
    #Vky = Vk*m.sin(theta)*m.sin(phi)
    Vky = Vk*m.cos(theta)
    return Vky
def FVkz(Vk,theta, phi):
    #this is not the system velovity direction
    #it is the orbital xyz
    #Vkz = Vk*m.cos(theta)
    Vkz = Vk*m.sin(theta)*m.sin(phi)
    return Vkz
def Fef(Vkz,Vky,Vr,Ai,MNS,M2,Af):
    ef = (1-(Vkz**2+(Vky+Vr)**2)*Ai**2/(G*(MNS+M2)*Af))**0.5
    return ef
def FVsysx(Vkx,Vky,Vkz,MNS,M2):
    #this XYZ is define by the orbital plane
    Vsysx = MNS/(MNS+M2)*Vkx  
    return Vsysx
def FVsysy(Vkx,Vky,Vkz,MNS,M2,M1,Vr):
    #this XYZ is define by the orbital plane
    Vsysy = MNS/(MNS+M2)*Vky-(M1-MNS)/(M1+M2)*(M2)/(MNS+M2)*Vr 
    return Vsysy
def FVsysz(Vkx,Vky,Vkz,MNS,M2):
    #this XYZ is define by the orbital plane
    Vsysz = MNS/(MNS+M2)*Vkz
    return Vsysz


#################################    
def FAf(MNS,M2,Ai,Vk,Vr,Vky):
    Af = G*(MNS+M2)/(2*G*(MNS+M2)/Ai-Vk**2-Vr**2-2*Vky*Vr)
    return Af
def FbG(MNS,M2):
    bG = 64/5*(G**3*MNS*M2*(MNS+M2))/c**5
    return bG

def FifHe(x): #be integ-ed fun
    return x**(29/19)*(1+121/304*x**2)**(1181/2299)/(1-x**2)**1.5
def FHe(ef):  #MT merger time
    inte, err = quad(FifHe,0,ef)
    He = (1-ef**2)**4*(ef**(-48/19))*(1+121/304*ef**2)**(-3840/2299)*inte
    return He

def FMT(bG, Af,He):  #Function Merger time  #He is the H(e) in the merger time formula
    MT = Af**4/bG*12/19*He
    return MT


#####################################
def FSiSM(ip):  #input solar mass converter
    SiSM = ip*SM
    return SiSM
def FSiSR(ip):  #input solar radii converter
    SiSR = ip*SR
    return SiSR
def FMwPD(x,a):    #maxwell prob distri func   #maximan is at 2**0.5*a
    pd = (2/m.pi)**0.5*x**2*m.exp(-x**2/(2*a**2))/a**3 #probability distribution
    return pd
def FBeta(MHe,M2,MNS):
    beta = (MNS+M2)/(MHe+M2)
    return beta
def FStGyr(sec):    #Function second to Gyr
    Gyr = sec/60/60/24/365.422/10**9
    return Gyr

