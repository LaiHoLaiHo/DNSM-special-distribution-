#!/usr/bin/python
#1.5 py 
#to calculation the vsys
#need to read par.csv
#it will output a file contain vsys megnitude and a mergertime file
import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random
import math 
from allfuc import FVky
from allfuc import FVkx
from allfuc import FVkz
from allfuc import FVkz
from allfuc import *
import pandas as pd
from scipy.integrate import quad
#read parameter
Par = pd.read_csv("/hetghome/jordan/miop/par.csv")
for i in Par:
    #print(i)
    #print(Par[i][0])
    if Par[i][0] != Par['lt'][0]:
        Par[i][0] = float(Par[i][0])

ft = Par.ft[0]
ft = 2.0
vs = Par.vs[0]     #v distribution range  for single meg. vk only
dlmv = Par.dlmv[0]   #maximan medium vsys meg
dlmvg = Par.dlmvg[0]
sn = Par.sn[0]  
#number of the grid velo not 
# particulary vsys or vk
ss = Par.ss[0]  #steps of image
r0 = Par.r0[0]
v0 = Par.v0[0]
lt = Par["lt"][0]
lt = "Aug26_22:16:35"
MNS = Par.MNS[0]
vg = Par["vvgg"][0]
M2 = Par.M2[0]
Ai = Par.Ai[0]
MHe = Par.MHe[0]
M1 = Par.MHe[0]
SiM1 = FSiSM(M1)
Vku = Par.Vku[0]
SiMHe = FSiSM(MHe)
SiM2 = FSiSM(M2)
SiMNS = FSiSM(MNS)
SiAi = FSiSR(Ai)
#read file  #############################
IRZ = pd.read_csv("/hetghome/jordan/miop/"+lt+"/IRZ.csv")
R = IRZ.initialR
Z = IRZ.initialz
AMT = []
ASVSYS = []
AVkL = [] #All kick V label
for v in range(int(dlmvg)): #this v is now represent maxwellian a
    #here is to put the all data inside the list
    PHI = []  #the coor in binary orbital
    THE = []  #the coor in binary orbital
    VSYS = []
    GyMt = []#Merger time in Gyr
    CTHE = [] #cos theta
    MT = [] # merger time
    ecce = []  #eccentriity
    apa = [] #alpha Af/Ai 
    Bond = [] #if it survives
    #test region begain ###############
    VKX = []
    VKY = []
    VKZ = []
    #test region end ###############
    print(v)
    v = v+1  #just let the dispersion of Maxwellian distribution is not zero
    VeRe = [] 
    SVk= []  #survived  Vk
    AVk = []
    for k in range((10)):  
        if  len(VeRe)<len(R): 
            VeRe = [] # velocity pool , with maxwell velo inside, #Non-D inside
            RPDD = []
            TMD = []
            for j in range(int(len(10*R)*(k+1))):
                RV = random.uniform(0.00000001,5*v) #Random v in non-D unit
                RPD = random.uniform(0,3*0.2*(v**-1))#Random Prob Distri
                if RPD < FPDMv(v,RV):
                    VeRe.append(RV) #Non-D inside
                    RPDD.append(RPD)
                    TMD.append(FPDMv(v,RV))
    plt.scatter(VeRe,TMD, s = 0.5)
    plt.scatter(VeRe,RPDD,s = 0.5)
    plt.savefig("/hetghome/jordan/miop/"+lt+"/"+"VTMD"+str(v)+".png")
    MV = Vku*v
    AVkL.append(str(MV)+'km $\mathregular{s^{-1}}$')
    for i in range(int(len(R))):
        #may be 100 is better per grid
        #Vr is the relative velocity in binary system
        #vbr velo binary relative
        SiVbr = FVbr(SiMHe,SiM2,SiAi)    #in si unit m/s  #binary relative velocity
        kmVbr = SiVbr/1000.0  #Vbr in km/s
        Vk = Vku*VeRe[i]
        MV = Vku*v
        #Vk = Vku*v*vg
        phi = random.uniform(0,2* math.pi)
        theta = math.acos(1-2*random.uniform(0, 1))
        PHI.append(phi)
        THE.append(theta)
        CTHE.append(math.cos(theta))
        Vkx = FVkx(Vk,theta,phi)
        Vky = FVky(Vk,theta,phi)
        #Vky = Vk*m.cos(theta)
        Vkz = FVkz(Vk,theta,phi)
        #Vkz = Vk*m.sin(theta)*m.sin(phi)
        SiVk = Vk*1000
        Vsysx = FVsysx(Vkx,Vky,Vkz,SiMNS,SiM2) #using the function in the allfuc.py
        Vsysz = FVsysz(Vkx,Vky,Vkz,SiMNS,SiM2)
        Vsysy = FVsysy(Vkx,Vky,Vkz,SiMNS,SiM2,SiM1,kmVbr)
        Vsys = (Vsysx**2+Vsysy**2+Vsysz**2)**0.5
        SiVky = Vky*1000 #from km/s to m/s  in SI unit
        SiVkz = Vkz*1000
        VKY.append(Vky)
        VKZ.append(Vkz)
        AVk.append(Vk)
        #print(Vky)
        #print(Vkz)
        SiAf = FAf(SiMNS,SiM2,SiAi,SiVk,SiVbr,SiVky)#  Af must larger than zero!!!!!!!
        if SiAf>0:
            #bG is beta G as a parameter to calculate the merger time
            SibG = FbG(SiMNS,SiM2)
            ef=Fef(SiVkz,SiVky,SiVbr,SiAi,SiMNS,SiM2,SiAf)
            ecce.append(ef)
            apa.append(SiAf/SiAi)
            He = FHe(ef)
            SiMt = FMT(SibG,SiAf,He)  #He no unit it is only integral e
            MT.append(SiMt)
            GyMt.append(FStGyr(SiMt))
            Bond.append(0)
            VKY.append(Vky)
            VKZ.append(Vkz)
            VSYS.append(Vsys)
            SVk.append(Vk)
        else:
            MT.append(np.nan)
            ecce.append(np.nan)
            apa.append(np.nan)
            VKY.append(np.nan)
            VKZ.append(np.nan)
            GyMt.append(np.nan)
            VSYS.append(np.nan)
            SVk.append(np.nan)
            Bond.append(-1)
    # test begain #######################
    #not necessary in the simulation
    #here is the plotting reagion to check that if the vsys is generated in a right way
    plt.scatter(VKY,VKZ,s = 0.03)
    plt.title("Vky-Vkz", loc = "right")
    plt.xlabel("Vky")
    plt.ylabel("Vkz")
    #plt.xlim(0, 40)
    #plt.ylim(-10, 10)
    #plt.xscale('log')
    plt.savefig("/hetghome/jordan/miop/"+lt+"/"+"MVky-z"+str(int(MV))+".png")
    plt.close("all")
    ####################
    plt.scatter(THE,PHI,s = 0.03)
    plt.title("THE-PHI", loc = "right")
    plt.xlabel("THE")
    plt.ylabel("PHI")
    #plt.xlim(0, 40)
    #plt.ylim(-10, 10)
    #plt.xscale('log')
    plt.savefig("/hetghome/jordan/miop/"+lt+"/"+"MTHE-PHI"+str(int(MV))+".png")
    plt.close("all")
    ###########################
    plt.scatter(apa,ecce,s = 0.03)
    plt.title("apa-ecce", loc = "right")
    plt.xlabel("apa")
    plt.ylabel("ecce")
    plt.xlim(10**-2, 10**2)
    plt.ylim(-0.1, 1.0)
    plt.xscale('log')
    plt.savefig("/hetghome/jordan/miop/"+lt+"/"+"MAE"+str(int(MV))+".png")
    plt.close("all")
    ########################
    DDDD = []
    for l in range(len(GyMt)):
        if GyMt[l] >0:
            DDDD.append(GyMt[l])
    SVSYS = []
    for l in range(len(VSYS)):
        if VSYS[l]>-1:
            SVSYS.append(VSYS[l])    
    plt.scatter(SVSYS,DDDD,s = 0.3) #Survive V SYS
    plt.title("Vsys vs Merger time", loc = "right")
    plt.xlabel("Vsys(km/s)")
    plt.ylabel("Merger time(Gyr)")
    #plt.xlim(10**-2, 10**2)
    plt.ylim(10**-3, 10**2)
    plt.yscale('log')
    plt.savefig("/hetghome/jordan/miop/"+lt+"/"+"MVsysMt"+str(int(MV))+".png")
    plt.close("all")  
    ############################
    """DDDD = []
    for l in range(len(GyMt)):
        if GyMt[l] >0:
            DDDD.append(GyMt[l])"""
    EEEE = []
    for i in range(len(DDDD)):
        if DDDD[i]>0:
            EEEE.append(float(m.log10(DDDD[i])))
    AMT.append(EEEE)
    plt.hist(EEEE,log = True,bins = 40)
    plt.title("Distribution of Merger time", loc = "right")
    plt.xlabel("Merger time(Gyr)(log scale)")
    plt.ylabel("Number of DNS")
    #plt.xlim(10**-2, 10**2)
    #plt.ylim(-0.1, 1.0)
    #plt.xscale('log')
    plt.savefig("/hetghome/jordan/miop/"+lt+"/"+"MtH"+str(int(MV))+".png")
    plt.close("all") 
    #######################
    AAAA = []
    for l in range(len(SVk)):
        if SVk[l]>-1:
            AAAA.append(SVk[l])
    plt.hist(AAAA,log = True,bins = 40)
    plt.title("Distribution of kick velocity(Survived)", loc = "right")
    plt.xlabel("Kick Velocity(km/s)")
    plt.ylabel("Number of DNS")
    #plt.xlim(10**-2, 10**2)
    #plt.ylim(-0.1, 1.0)
    plt.yscale('linear')
    plt.savefig("/hetghome/jordan/miop/"+lt+"/"+"MVkH"+str(int(MV))+".png")
    plt.close("all") 
    #################################
    plt.hist(AVk,log = True,bins = 40)
    plt.title("Distribution of All kick velocity", loc = "right")
    plt.xlabel("Kick Velocity(km $\mathregular{s^{-1}}$)")
    plt.ylabel("Number of DNS")
    #plt.xlim(10**-2, 10**2)
    #plt.ylim(-0.1, 1.0)#
    plt.yscale('linear')
    plt.savefig("/hetghome/jordan/miop/"+lt+"/"+"MAVkH"+str(int(MV))+".png")
    plt.close("all") 
    ##################################
    plt.scatter(VeRe,RPDD,s = 0.3) #Survive V SYS
    plt.title("Test Maxwell", loc = "right")
    plt.xlabel("Vsys(Non-D)")
    plt.ylabel("Probability distribution")
    #plt.xlim(10**-2, 10**2)
    #plt.ylim(10**-3, 10**2)
    plt.yscale('linear')
    plt.savefig("/hetghome/jordan/miop/"+lt+"/"+"MDT"+str(int(MV))+".png")
    plt.close("all")  
    ############################
    ASVSYS.append(SVSYS)
    plt.hist(SVSYS,log = True,bins = 40)
    plt.title("Distribution of Systemetic V", loc = "right")
    plt.xlabel("Systemic Velocity")
    plt.ylabel("Number of DNS")
    #plt.xlim(10**-2, 10**2)
    #plt.ylim(-0.1, 1.0)
    #plt.xscale('log')
    plt.savefig("/hetghome/jordan/miop/"+lt+"/"+"MVsysH"+str(int(MV))+".png")
    plt.close("all")   
    # test end #######################
    print(SiAi)
    print(SiM2)
    print(SiMHe)
    print('vbr='+str(kmVbr))
    print("Vsysx="+str(Vsysx))
    IVSYS = pd.DataFrame(VSYS,columns = ['VsysMeg_kms'])
    IVSYS.insert(1,'Bond',Bond)
    IVSYS.insert(1,'alpha',apa)
    IVSYS.insert(1,'ecce',ecce)
    MMT = pd.DataFrame(MT,columns = ['MT_sec'])
    MMT.to_csv("/hetghome/jordan/miop/"+str(lt)+"/MT_Mvk~"+str(MV)+".csv",index = False)
    FA = pd.DataFrame(CTHE,columns = ['cos_theta'])
    FA.insert(1,"phi",PHI)
    FA.to_csv("/hetghome/jordan/miop/"+str(lt)+"/angles_Mvk~"+str(MV)+".csv",index = False)
    #FA.csv is not used in the following code
    #initial system velo
    IVSYS.to_csv("/hetghome/jordan/miop/"+str(lt)+"/IVSYS_MVk~"+str(MV)+".csv",index = False)
    print(v)
    print("fin")

plt.hist([AMT[0],AMT[1],AMT[2]],log = True,bins = 40,label = AVkL)#,stacked=True)
plt.title("Distribution of Merger time", loc = "right")
plt.xlabel("Merger time(Gyr)(log scale)")
plt.ylabel("Number of DNS")
plt.legend()
plt.savefig("/hetghome/jordan/miop/"+lt+"/"+"AMtH"+".png",dpi = 1200)
plt.close("all")

plt.hist([ASVSYS[0],ASVSYS[1],ASVSYS[2]],log = True,bins = 40,label = AVkL)#,stacked=True)
plt.title("Distribution of Systemic Velocity", loc = "right")
plt.xlabel("Systemic Velocity(km/s)")
plt.ylabel("Number of DNS")
plt.legend()
plt.savefig("/hetghome/jordan/miop/"+lt+"/"+"ASVSYSH"+".png",dpi = 1200)
plt.close("all")
plt.clf()

plt.scatter(ASVSYS[0],AMT[0],label = AVkL[0],s = 0.1)#,c = 'b')
plt.scatter(ASVSYS[1],AMT[1],label = AVkL[1],s = 0.1)
plt.scatter(ASVSYS[2],AMT[2],label = AVkL[2],s = 0.1)
plt.title("Systemic Velocity vs Merger time", loc = "right")
plt.xlabel("Vsys(km/s)")
plt.ylabel("Merger time(Gyr)(log scale)")
#plt.xscale('log')
#plt.yscale('log')
plt.legend()
plt.savefig("/hetghome/jordan/miop/"+lt+"/"+"ASVSMT"+".png",dpi = 1200)
plt.close("all")

