#!/usr/bin/python
#generate initial position 
#1111111111111111111111
"""generate position of samples"""
#and it will out put the IRZ.csv for the initial position
print("begin")
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
import galpy.potential
from galpy.orbit import Orbit
import numpy
from galpy.potential import KeplerPotential
from galpy.util import bovy_conversion
from galpy.potential import PowerSphericalPotentialwCutoff
from galpy.potential import MiyamotoNagaiPotential
from galpy.potential import NFWPotential
from astropy import units
from galpy.potential import evaluatePotentials
from galpy.potential import plotRotcurve
import math
from galpy.potential import plotPotentials
import random
from galpy.potential import evaluateDensities
from galpy.potential import MWPotential2014
import csv
import pandas as pd
import imageio
import math
import os
import time 
#write time in par file
localtime = time . asctime ( time . localtime ( time . time ( ) ) ) 
tc = localtime.split( )
lt = tc[1]+tc[2]+"_"+tc[3]
os.makedirs("/hetghome/jordan/miop/"+lt)

Par = pd.read_csv("/hetghome/jordan/miop/par.csv")

Par["lt"][0] = lt
Par.to_csv("/hetghome/jordan/miop/par.csv",index = False)


"""
mp = []
test = 30
MP = pd.DataFrame(mp, columns = ['test'])

for i in range(len(Par)):
    if Par[i] = Par['lt']:
        mp = []  #just a mid product
        mp.append(float(Par[i][0]))
"""
#re-save the par file to encode the date and time in the par file
Par.to_csv("/hetghome/jordan/miop/"+lt+"/par.csv",index = False)
Par = pd.read_csv("/hetghome/jordan/miop/"+lt+"/par.csv")
if 'HOME' not in os.environ:
    import pwd
    userhome = pwd.getpwuid(os.getuid()).pw_dir
else:
    userhome = os.environ['HOME']
"""
#############################################################
#
#galpy origional
bpg= PowerSphericalPotentialwCutoff(alpha=1.8,rc=1.9/8.,normalize=0.05)
mpg= MiyamotoNagaiPotential(a=3./8.,b=0.28/8.,normalize=.6)
npg= NFWPotential(a=16/8.,normalize=.35)
MWPotential2014g= [bpg,mpg,npg]
MWPotential2014gwBH= [MWPotential2014g,KeplerPotential(amp=4*10**6./bovy_conversion.mass_in_msol(220.,8.))]
#generate a New born DNS
# limit unit kpc  density 3 solarmass/ pc3
################################################################
"""

print(lt)
print(Par['lt'][0])
print("read parameter")
for i in Par:
    #print(i)
    #print(Par[i][0])
    if Par[i][0] != Par['lt'][0]:
        Par[i][0] = float(Par[i][0])
rr = Par.rr[0]
zr = Par.zr[0]
dr = Par.dr[0]
sn = Par.sn[0]
r0 = Par.r0[0]
v0 = Par.v0[0]

#print(type(sn))

#print("dr = " + str(dr))
"""rr = 30
zr = 5
dr = 3
sn = 1000"""
##################################
#start to generate initial RZ by Monte Carlo 
R = []
Z = []
D = []
RR = []
for i in range(int(sn)):
    r = random.uniform(0, 1)*rr
    #RR.append(r)
    z = random.uniform(-1, 1)*zr
    d = random.random()*dr
    td = evaluateDensities(MWPotential2014,r/r0,z/r0)
    if d < td:
    #bpg.dens(r/8,z/8)+mpg.dens(r/8,z/8)+npg.dens(r/8,z/8):
        R.append(r)
        Z.append(z)
IRZ = pd.DataFrame(R, columns = ['initialR'])
IRZ.insert(1, 'initialz',Z)
IRZ.to_csv("/hetghome/jordan/miop/"+lt+"/IRZ.csv",index = False)
print("IRZ out put")
#plt.hist(R, bins =  numpy.linspace(0,30,100))
#plt.show()
"""plt.hist(RR, bins =  numpy.linspace(-31,31 , 63))
plt.title("Random Number dirtribution  loop")
plt.show()"""
print(len(R))
"""
plt.hist(Z, bins =  numpy.linspace(-zr-2,zr+2,(zr+2)*4+1))
plt.title("Initial z position distribution")
plt.xlabel("z (kpc)")
plt.ylabel("DNS Number")
plt.savefig("/hetghome/jordan/iop/"+lt+"/ZID.png")
plt.clf()"""
print("fin")
print(lt)
#make new direction to put all the followig files into this path
newPath = shutil.copy('/hetghome/jordan/miop/par.csv', '/hetghome/jordan/miop/'+lt+'/')
