# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 15:11:10 2021

@author: jim
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 08:26:54 2020

@author: BK
"""

import cmath
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mplot
import scipy.io
from utils import *
import pandas as pd

# matlab 參數匯入
mat = scipy.io.loadmat('ar1.mat')

f = np.logspace(1.301,4.301,1000)
w = 2*np.pi*f 
j = cmath.sqrt(-1)
s=w*j

rho0 = 1.29                                       
c = 343           
k=w/c                           
mu = 1.86e-5                                                                    
Bi = 8/(3*np.pi)
Br = 1/2                           
r = 0.1                                    
Pref = 20e-6

# TS Paranmeters
eg = 0.217                                       
Bl = 1.46                                      
AD = 5.54e-4 
aD = math.sqrt(AD/np.pi)                    
Re = 43.39
Le = 0.09e-3                           
R2 = 3.92
L2 = 0.05e-3
Rm = 0.01
Cm = 15.88e-3
Mm = 0.04e-3          

# acoustic mechanism sizes
Aarh1 = 3.5e-6                                  
aarh1 = math.sqrt(Aarh1/np.pi)                                 
larh1 = 0.1e-3                            
Aarh2 = 1.018e-4                              
aarh2 = math.sqrt(Aarh2/np.pi)                                
larh2 = 1.5e-3                             
Varc  = 1.99e-6                            
larc  = 2e-3

# mechanic domain 
Rafrad = (Bi**2/Br)*((rho0*c)/(AD))                   
Mafrad = (Bi)*((rho0*aD)/AD)

# acoustic domain                                                                                                                                                                                                                                                                                 «á         
Carc = Varc/(rho0*c**2)                     
Rarh1 = (8*mu*larh1)/(np.pi*aarh1**4)                     
Marh1 = (4*rho0*larh1)/(3*np.pi*aarh1**2)        
Rap = 0
Map = 0
Nah = 9
Rarrad1 = (Bi**2/Br)*((rho0*c)/Aarh1)                  
Marrad1 = (Bi)*((rho0*aarh1)/Aarh1)
Rarrad2 = (Bi**2/Br)*((rho0*c)/Aarh2)                  
Marrad2 = (Bi)*((rho0*aarh2)/Aarh2)    


df=pd.read_excel("SPLunmesh.xlsx",sheet_name=4)
Hz_memsure=df.iloc[:,0]
voltage=df.iloc[:,1]

df=pd.read_excel("impedancecurvenew.xlsx",sheet_name=4)
Hz_measure_2=df.iloc[:,0]
impedance_measurement=df.iloc[:,1]

ecm_temp=[]
sub_ecm_temp=[]
mag=np.zeros([len(f),1],dtype=complex)
volocity=np.zeros([len(f),1],dtype=complex)
impedance=np.zeros([len(f),1],dtype=complex)

dB=np.zeros([len(f),1],dtype=complex)
Zafrad = par(Rafrad,s*Mafrad)                   
Zeb = Re+s*Le+par(s*L2,R2)           
Zm = Rm+s*Mm+(1/(s*Cm))                 
Zarc = 1/(s*Carc)                  
Zarh1 = Rarh1+Marh1                   
Zarh2 = j*((rho0*c)/Aarh2)*np.tan((k*larh2)/2)            
Zarh2c = ((rho0*c)/Aarh2)/(j*np.sin(k*larh2))           
Zap = (Rap+s*Map)/Nah     
Zarrad1 = par(Rarrad1,s*Marrad1)       
Zarrad2 = par(Rarrad2,s*Marrad2)      
Z1=(Bl**2/Zeb)+Zm+((Zafrad+Zarc)*(AD**2))
Z2=(Zarh1+Zarrad1+Zap+Zarh2+Zarh2+Zarrad2)*(AD**2)
Z3=(Zarc+Zap+Zarh2+Zarh2c)*(AD**2);
Z4=(Zarh2+Zarh2c+Zarrad2)*(AD**2);
Z13=-Zarc*(AD**2)
Z23=-(Zap+Zarh2)*(AD**2)
Z24=-(Zarh2+Zarrad2)*(AD**2)
Z34=-Zarh2c*(AD**2)

A = Zarh2+Zarrad2
B = Zarh2c
C = Zap+Zarh2
D = Zarh1+Zarrad1
E = Zarc
F = Zafrad
A1 = par(A,B)
A2 = A1+C
A3 = par(A2,D)
A4 = par(A3,E)
A5 = F+A4
A6 = (AD**2)*(A5)
Ztotal = np.abs(Bl**2/(A6+Zm)+Zeb);

for i in range(len(f)):

    ecm_temp=circuit_list2mat_diagonal([Z1[i],Z2[i],Z3[i],Z4[i]],anti_diagonal=0)
    sub_ecm_temp=([0,Z13[i],0],[0,Z23[i],Z24[i]],[0,0,Z34[i]])
    sub_ecm_temp=np.asarray(sub_ecm_temp)
    ecm_temp=circuit_mat_subelement(ecm_temp,sub_ecm_temp)
    impedance[i],_=(circuit_impedence_cal(ecm_temp,1,1,0))
    volocity[i]=impedance[i]*((eg*Bl)/Zeb[i])
    mag=np.abs((s*rho0*AD*volocity[i])/(2*np.pi*r))
    dB[i]=20*cmath.log10(mag[i]/Pref)

# fig1, ax1 = plt.subplots(figsize=(20,12))
# plt.grid()
# plt.title("ECM of electret condenser microphone")
# plt.semilogx(f,dB)
# plt.semilogx(f,Ztotal)
# plt.semilogx(Hz_memsure,voltage)
# ax1.set_xticks([20,100,1000,3000,4000,5000,6000,10000,20000])
# ax1.get_xaxis().set_major_formatter(mplot.ticker.ScalarFormatter())
# plt.savefig("ecm.png",dpi=900)
# plt.show()

fig, host=plt.subplots(figsize=(20,12))
par1=host.twinx()
host.set_xlabel("frequency")
host.set_ylabel("magnititude")
par1.set_ylabel("Ohm")
color1=plt.cm.viridis(0)
color2=plt.cm.viridis(0.5)
p1,=host.plot(f,dB,label="frequency response",color=color1)
p2,=par1.plot(f,Ztotal,label="impedance curve",color=color2)
p3,=host.plot(Hz_memsure,voltage,'--',label="measurement frequency response",color=color1)
p4,=par1.plot(Hz_measure_2,impedance_measurement,'--',label="measurement impedance curve",color=color2)
host.yaxis.label.set_color(p1.get_color())
par1.yaxis.label.set_color(p2.get_color())
host.set_xscale("log")
host.set_xticks([20,100,1000,3000,4000,5000,6000,10000,20000])
par1.set_yticks([50,100,150,200,250,300,350,400])
host.get_xaxis().set_major_formatter(mplot.ticker.ScalarFormatter())
lns = [p1, p2, p4]
host.legend(handles=lns, loc='best')
plt.draw()
plt.savefig("ecm.png",dpi=900)
plt.show()
