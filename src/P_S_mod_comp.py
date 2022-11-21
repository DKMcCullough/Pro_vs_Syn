#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 11:50:19 2022
name: P_S_mod_comp 

location: /Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Competitions/Pro_vs_Syn/src

Goals: 
    make Pro and Syn modular competition file for leaky and non-leaky HOOH detox via Syn in N-limited growth
    create 2 parameter sets for leaky or non-leaky Pro vs Syn Competition 
    differnt Nstar and Hstars depending on who wins and what equilibrium is being/
        reached by what means? 



@author: dkm

"""




import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import *
from scipy.integrate import odeint



##################################################3
# parameter and variable Set UP 
#############################################


step = 0.001
ndays = 1000
mtimes = np.linspace(0,ndays,int(ndays/step))
SNs = np.linspace(0, 3000, num = 15)
Shs = np.linspace(0, 3000, num =15)
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
        Psc = leaky[:,0]
        Ssc = leaky[:,1]
        Nsc = leaky[:,2]
        Hsc = leaky[:,3]
        #print(Z[i,j])
        if (i == 1) and (j == 0):
            Ps = leaky[:,0]
            Ss = leaky[:,1]
            Ns = leaky[:,2]
            Hs = leaky[:,3]
        Ssc_av = np.mean(Ssc[-200:])
        Psc_av = np.mean(Psc[-200:])
        ratio = (Ssc_av/( Ssc_av+ Psc_av))
        Z[i,j] = ratio
        if np.all([g <= 1e-3 for g in Psc[-200:]]) and np.all([h >= 10 for h in Ssc[-200:]]) : 
            Z[i,j] = 1
        if np.all([g >= 10 for g in Psc[-200:]]) and np.all([h <= 1e-3 for h in Ssc[-200:]]) : 
            Z[i,j] = 0
        if np.all([g <= 1e-3 for g in Psc[-200:]]) and np.all([h <= 1e-3 for h in Ssc[-200:]]) : 
            Z[i,j] = -1
        #elif Z[i,j] > 10:
            #print(i,j)
            #print(Ssc[-1],Psc[-1])
            #sys.exit()'''



#####################################

#  Graphing dynamic model 

#####################################


#fig,ax1 = plt.subplots()
fig, (ax1, ax2,ax3) = plt.subplots(3,1, sharex=True, figsize=(9,5))
fig.suptitle('Growth Competition Projections')
plt.subplots_adjust(wspace = 0.5, top = 0.9,bottom = 0.1)


ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'Pro ') #'k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn ') #'k1 =' + str(k1s))
ax1.set(xlabel='Time (days)', ylabel='cells per ml')
#ax1.set_ylim(bottom = -20)


ax2.plot(mtimes, Ns, linewidth = 3, color = 'purple', label = "Nutrient Concentration ")
ax2.set(xlabel='Time (days)', ylabel='Nutrient [ ]')

ax3.plot(mtimes, Hs,  linewidth = 3, color = 'red', label = "HOOH concentration ")
ax3.set(xlabel='Time (days)', ylabel='HOOH [ ]')

ax1.semilogy()
ax2.semilogy()
ax3.semilogy()

#ax1.legend(loc = 'lower right')




#plt.show()


##########################################################

# Calculated analytical solutions at equilibrium

##########################################################


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

ax1.axhline(Sstar,color = 'orange', linestyle = "-.",label = 'Sstar')
ax1.axhline(Pstar,color = 'green', linestyle = "-.",label = 'Pstar')
ax2.axhline(Nstar,color = 'purple', linestyle = "-.",label = 'Nstar')
ax3.axhline(Hstar,color = 'red', linestyle = "-.",label = 'Hstar')

ax1.legend(loc = 'best')
ax2.legend(loc = 'best')
ax3.legend(loc = 'best')

fig.savefig('../figures/leaky_calcs_auto',dpi=300)



#######################################
# Graphing Cotour plots from 
######################################

fig3,(ax1) = plt.subplots(sharex = True, sharey = True, figsize = (8,5))
fig3.suptitle('Leaky HOOH Contour')

grid = ax1.pcolormesh( Shs, SNs, np.where(Z == -1, np.nan, Z), vmin=0, vmax=np.max(Z), cmap = 'summer', shading= 'auto'  )  #'gouraud'
ax1.set(xlabel='SupplyHOOH')
ax1.set(ylabel='Supply nutrient')

#np.max(Z)
ax1.axhline((rho*Nstars),color = 'purple', linestyle = "-.",label = 'SN cutoff for S growth?')
ax1.axhline((rho*Nstarp),color = 'magenta', linestyle = "-.",label = 'SN cutoff for P growth?')

ax1.axvline((vHline),color = 'c', linestyle = "-.",label = 'H cutoff?')
#ax1.axvline((Hstar*(deltah+(phi*Pstar))),color = 'b', linestyle = "-.",label = 'H cutoff?')

#fig3.legend()
fig3.colorbar(grid, cmap= 'summer',label = 'S / S+P')
#plt.colorbar.set_label('S/P+S')

fig3.savefig('../figures/no_leak_contour_f3_auto',dpi=300)

print('') #printing blank line 
print('*** Nstars? or Sh cut off?  ***')
print('*** Done ***')