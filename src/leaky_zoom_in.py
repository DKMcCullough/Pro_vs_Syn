#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name:   leaky_contour.py 

location: /Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Competitions/Pro_vs_Syn/src

Goal: 
    Compete Pro and Syn on one nutrient N in leaky HOOH detox 
    Syn has worse k1 value to compensate for HOOH detox energy/size costs
    Pro has kdam caused by HOOH 
    HOOH detoxed innately (deltah) or by Syn (phi)
get coexistance in range of H and N supply rates? 

@author: dkm

"""




import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import *
from scipy.integrate import odeint
import sys

##################################################3
# parameter and variable Set UP 
#############################################


step = 0.001
ndays = 2500
mtimes = np.linspace(0,ndays,int(ndays/step))
SNs = np.linspace(0, 1000, num = 20)
Shs = np.linspace(0, 1000, num =20)
Z = np.zeros((int(SNs.shape[0]),int(Shs.shape[0])),float)


#parameters
Qnp = 1#(9.4e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertilison 
Qns = 1#(20.0e-15*(1/(14.0))*1e+9) 

k1p =  0.00002     #Pro alpha
k1s =  0.00001      #Syn alpha 
k2 =  0.88    #Vmax    shared for P and S here
ksp = k2/k1p
kss = k2/k1s
dp = 0.2   #pro delta
ds =  0.2   #syn delta
kdam = 0.005   #hooh mediated damage rate of Pro  
deltah = 0.002       #decay rate of HOOH via Syn 
phi = 0.02    #0007  #detoxification-based decay of HOOH via Syn in this case
rho =  0.002

#params = [ksp,kss,k2,dp,ds,kdam,deltah,phi,rho]



#empty arrays to be populated by odeint when calling leak function
P = np.array([])
S  = np.array([])
N = np.array([])
H = np.array([])
y = [P,S,N,H]

#initial values to be used for odeint start 
P0 = 1e4
S0 = 1e4
N0 = 1.0e5        #nM 
H0 = 1     #nM
inits = (P0,S0,N0,H0)

## ODE functions  to be solved with odeint 


def leak(y,t,params):
    ksp,kss,k2,dp,ds,kdam,deltah,phi,rho,SN,Sh = params[0], params[1], params[2], params[3],params[4], params[5], params[6], params[7],params[8],params[9], params[-1]
    P,S,N,H = y[0],y[1],y[2],y[3]
    dPdt = (k2 * N /( (ksp) + N) )*P*Qnp - (dp *P) - kdam*H*P
    dSdt =(k2 * N /( (kss) + N))*S*Qns - (ds *S) #- kdams*H*S      
    dNdt =  SN - ((k2 * N /( (ksp) + N) )*P* Qnp) - ((k2 * N /( (kss) + N))*S* Qns) - rho*N    
    dHdt = Sh - deltah*H  -phi*S*H  #phi being S cell-specific detox rate
    return [dPdt,dSdt,dNdt,dHdt]



##############################
#Contour and model calculations 
##############################

for (i,SN) in zip(range(SNs.shape[0]),SNs):
    for (j,Sh) in zip(range(Shs.shape[0]),Shs):
        params = [ksp,kss,k2,dp,ds,kdam,deltah,phi,rho, SN, Sh]
        leaky  = odeint(leak, inits, mtimes, args = (params,))
        if (i == 17) and (j == 3):
            Ps = leaky[:,0]
            Ss = leaky[:,1]
            Ns = leaky[:,2]
            Hs = leaky[:,3]


#####################################

#  Graphing dynamic model 

#####################################



fig4,  (ax1,ax2) = plt.subplots(2, 1, sharex=True,figsize=(9,5))
fig4.suptitle('Leaky Zoom in')

ax1.plot(mtimes, np.clip(Ps,1,np.max(Ps)) , linewidth = 3, color = 'g', label = 'Pro')
ax1.plot(mtimes, np.clip(Ss,1,np.max(Ss)), linewidth = 3, color = 'orange', label = 'Syn')

ax1.set(ylabel='cells per ml')

ax2.set_ylabel('Nutrient (per ml)')
ax2.plot(mtimes, np.clip(Ns,1,np.max(Ns)),linewidth = 3, color = 'purple', label = "Nutrient")

fig2, ax3 = plt.subplots()
ax3.plot(mtimes, np.clip(Hs,1,np.max(Hs)),linewidth = 3, color = 'red', label = "HOOH")

ax3.set(xlabel='Time (days)', ylabel='HOOH per ml')


ax1.semilogy()
ax2.semilogy()


###### coesistance equilibrium (star equations ) ###########

#Coexist
Nstar = (kss*ds)/(k2-ds)
Hstar = (((k2*Nstar)/(Nstar + ksp))-(dp))*(1/kdam)
Sstar = (Sh - deltah*Hstar)/(phi*Hstar)
Pstar = ((SN-rho*Nstar)*(Nstar + ksp))/(k2*Nstar*Qnp)



#Pwin 
Nstarp = ((ksp*dp )+(ksp*kdam))/((k2*Qnp) - dp - kdam)
Pstarp = (SN - rho*Nstarp)*((Nstarp + ksp)/((k2*Nstarp)*Qnp))
Hstarp = Sh/(deltah+phi*Pstar)  #do we need toassume H must be 0 for P to win?????

#Swin 
Nstars = (ds*kss)/((k2*Qns)-ds)
Sstars = (SN - rho*Nstars)*(((Nstars + kss)/(k2*Nstars*Qns)))
Hstars = Sh/(deltah)


Nstarph = ((ksp*dp )+(ksp*kdam*Hstar))/((k2*Qnp) - dp - (kdam*Hstar))
vHline = ((deltah)/(Pstar*kdam)*((Nstarp+ksp)/(k2*Nstarp*Pstar*Qnp)+(dp*Pstar)))



##### graphing stars equations ############### 

ax1.axhline(Sstar,color = 'brown', linestyle = "-.",label = 'S*')
ax1.axhline(Pstar,color = 'green', linestyle = ":",label = 'P*')
#ax2.axhline(Nstar,color = 'purple', linestyle = "-.",label = 'Nstar')

ax2.axhline(Nstars,color = 'purple', linestyle = "-.",label = 'N*s')
ax2.axhline(Nstarp,color = 'magenta', linestyle = ":",label = 'N*p')
ax3.axhline(Hstar,color = 'red', linestyle = "-.",label = 'Hstar')

ax1.legend(loc = 'lower center')
ax2.legend(loc = 'lower center')
ax3.legend(loc = 'lower center')

#fig.savefig('../figures/leaky_zoom_in_auto',dpi=300)



print('') #printing blank line 
print('*** Nstars? or Sh cut off?  ***')
print('*** Done ***')