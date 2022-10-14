#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name:   leaky_hooh_detox_calcs.py 

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


##################################################3
# parameter and variable Set UP 
#############################################


step = 0.01
ndays = 300
mtimes = np.linspace(0,ndays,int(ndays/step))
Shs = np.linspace(0, 200, num=10)
SNs = np.linspace(0, 1000, num = 10)
Z = np.zeros((int(SNs.shape[0]),int(Shs.shape[0])),int)


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
kdam = 0.05   #hooh mediated damage rate of Pro  
deltah = 0.002       #decay rate of HOOH via Syn 
phi = 0.0002    #0007  #detoxification-based decay of HOOH via Syn in this case
rho =  0.002

params = [ksp,kss,k2,dp,ds,kdam,deltah,phi,rho]



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


 ###############################

 # ODE models broken down 

 ###################################

'''
 plankton = (max growth rate)(nutient concentration)/michelis-mention constant) + (nutrinet concentration) (Population size) 
 nutrients = (Supply of nutrient) - (max growth rate)(nutient concentration)/michelis-mention constant) + (nutrinet concentration)


 k1 (P or S) = alpha value aka. nutrient uptake rate(per day rate)
 k1 differes between S and P in this model and is used to as drawback for S having catalase detox
 k2 = vmax   aka. usage maximum of nutrient(per day rate)
 P = Prochlorococcus abundance (cells/ml)
 S = Synechococcys abundance (cells/ml)
 N = nutrient concetration (nM conentration per ml)
 H = HOOH concentration 


'''
##############################
#Contour and model calculations 
##############################




for (i,SN) in zip(range(SNs.shape[0]),SNs):
    for (j,Sh) in zip(range(Shs.shape[0]),Shs):
        params = [k1p,k1s,k2,dp,ds,kdam,deltah,rho, SN, Sh]
        leaky  = odeint(leak, inits, mtimes, args = (params,))
        Psc = leaky[:,0]
        Ssc = leaky[:,1]
        Nsc = leaky[:,2]
        Hsc = leaky[:,3]
        if (i == 8) and (j == 1):
            #Hsupply = Shs[j]
            #Nsupply = SNs[i]
            Ps = leaky[:,0]
            Ss = leaky[:,1]
            Ns = leaky[:,2]
            Hs = leaky[:,3]
        Z[i,j] = Ssc[-1]/(Psc[-1]+Ssc[-1])
        if np.all([g <= 1e-3 for g in Psc[-10:]]) and np.all([h <= 1e-3 for h in Ssc[-10:]]) : 
            Z[i,j] = -1
        print('Ssc = '+ str(Ssc[-1]) + ' and Psc = '+ str(Psc[-1]))
        print(Z[i,j])



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
###### with QN ##################

Nstar = (kss*ds)/(k2-ds)

Hstar = (((k2*Nstar)/(Nstar + ksp))-(dp))*(1/kdam)


Sstar = (Sh - deltah*Hstar)/(phi*Hstar)


Pstar = (SN - (((k2*Nstar)/(Nstar + kss))*Sstar*Qns) - (rho*Nstar))*((Nstar + ksp)/(k2*Nstar*Qnp))





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


grid = ax1.pcolormesh( Shs, SNs, np.where(Z == -1, np.nan, Z), vmin=0, vmax=np.max(Z), cmap = 'summer', shading= 'auto'  )  #'gouraud'
ax1.set(xlabel='Supply hooh')
ax1.set(ylabel='Supply nutrient')

#ax1.axvline((deltah*Hstar),color = 'purple', linestyle = "-.",label = 'H cutoff??')
#ax1.axhline(((rho*Nstar)+(((k2*Nstar)/(Nstar-kss))*Sstar*Qns)),color = 'magenta', linestyle = "-.",label = 'Sn cut off on S?')

plt.legend()
fig3.colorbar(grid, cmap= 'summer',label = 'S / S+P')
#plt.colorbar.set_label('S/P+S')

fig3.savefig('../figures/no_leak_contour_f3_auto',dpi=300)

print('*** SStar =/= dsdt solutions ***')
print('*** Done ***')