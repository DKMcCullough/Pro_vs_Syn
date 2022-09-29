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





################################

# model 

###################################

'''
plankton = (max growth rate)(nutient concentration)/michelis-mention constant) + (nutrinet concentration) (Population size) 
nutrients = (Supply of nutrient) - (max growth rate)(nutient concentration)/michelis-mention constant) + (nutrinet concentration)


k1 (p or s) = alpha (per day rate)
k2 = vmax (per day rate)
P = Prochlorococcus abundance (cells/ml)
S = Synechococcys abundance (cells/ml)
N = nutrient concetration (nM conentration per ml)
H = HOOH concentration 


'''
#####################
# Set UP
######################


step = 0.01
ndays = 300
mtimes = np.linspace(0,ndays,int(ndays/step))
Shs = np.linspace(0, 200, num=8)
SNs = np.linspace(0, 1000, num = 8)
Z = np.zeros((int(SNs.shape[0]),int(Shs.shape[0])),int)

#initial values 
P0 = 1e4
S0 = 1e4
N0 = 1.0e5        #nM 
H0 = 1     #nM
inits = (P0,S0,N0,H0)

#parameters

k1p =  0.00002     #Pro alpha
k1s =  0.00001      #Syn alpha 
k2 =  0.88    #Vmax    shared for P and S here
dp = 0.2   #pro delta
ds =  0.2   #syn delta
kdam = 0.05   #hooh mediated damage rate of Pro  
deltah = 0.002       #decay rate of HOOH via Syn 
phi = 0.0002    #0007  #detoxification-based decay of HOOH via Syn in this case
rho =  0.002
Nsupply = SNs[7] #max(SNs)   
Hsupply = Shs[7] #max(Shs)

params = [k1p,k1s,k2,dp,ds,kdam,deltah,phi,rho,Nsupply,Hsupply]

#Nsupply = N0
#Hsupply = H0
Qnp = (9.4e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertilison 
Qns = (20.0e-15*(1/(14.0))*1e+9) 
ks1 = k2/k1p
ks2 = k2/k1s

#empty arrays 
P = np.array([])
S  = np.array([])
N = np.array([])
H = np.array([])
y = [P,S,N,H]


#function set up for ode int


def leak(y,t,params):
    k1p,k1s,k2,dp,ds,kdam,deltah,phi,rho,Nsupply,Hsupply = params[0], params[1], params[2], params[3],params[4], params[5], params[6], params[7],params[8],params[9], params[-1]
    P,S,N,H = y[0],y[1],y[2],y[3]
    dPdt = (k2 * N /( (ks1) + N) )*P*Qnp - (dp *P) - kdam*H*P
    dSdt =(k2 * N /( (ks2) + N))*S*Qns - (ds *S) #- kdams*H*S      
    dNdt =  Nsupply - ((k2 * N /( (ks1) + N) )*P* Qnp) - ((k2 * N /( (ks2) + N))*S* Qns) - rho*N    
    dHdt = Hsupply - deltah*H  -phi*S*H  #phi *S*H with phi being Scell-specific detox rate
    #if t<0.05:
        #print(muP,muS)
    return [dPdt,dSdt,dNdt,dHdt]

#solve ODEs via odeint
leaky  = odeint(leak, inits, mtimes, args = (params,))

#redefine where P and N are in returned matrix from ode int
Ps = leaky[:,0]
Ss = leaky[:,1]
Ns = leaky[:,2]
Hs = leaky[:,3]
 #need to set a limit on y lower axis bound bc cutting it off at graph just leaves a blank graph

#####################################

#  Graphing

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

'''
# Models
    #muP = (k2 * N /( (k2/k1p) + N) )
    #muS = (k2 * N /( (k2/k1s) + N) )
    dPdt = P * muP - (dp *P) - kdamp*H*P
    dSdt = S * muS - (ds *S) - kdams*H*S      
    dNdt =  Nsupply - (muP * P * Qnp) - (muS * S * Qns) -rho*N)
    dHdt = S_HOOH - deltah*H  - phi*S*H 
    
#Analytical solutions at equilibrium (dx/dt equations above set to 0 and solved)

#Equilibrium  (5 cases: both die, Syn dies after enough H is gone, coexistance, and Pro dies in 2 different equs due to radical/signs )
'''

###### coesistance equilibrium without Qns ###########

'''
###matmatica
#Nstar = - (ds*k2)/((k1s)*( ds-k2))
#Hstar = ((k1p*ds*(dp-k2)*k2 + (k1s*dp*k2*(-ds + k2)) / (kdam*((k1s*k2*(ds-k2)-k1p*ds*k2)))))
#Sstar = -(k1s*(deltah*dp + Hsupply*kdam)*k2*(ds-k2) + k1p*ds*(deltaH*dp + Hsupply*kdam-deltah*k2)*k2) / ((k1s*dp*k2*(ds-k2)+k1p*ds*(-dp+k2)*k2)*phi)
#Pstar = ((k1p*k2*(ds-k2)-k1p*ds*k2)*(k1p*((ds*ds))*(deltah*dp + Hsupply*kdam-deltah*k2)((k2*k2*k2)) + ((k1s*k1sk1s))*k2*((ds-k2)*(ds-k2))*(Hsupply*kdam*k2 + dp*(deltah*k2 + Nsuppply*phi))-((k1s*k1s))*ds*(ds-k2)*(k2)*(-rho*dp*k2*phi + k1p*(Hsupply*kdam*k2 +dp*(deltah*k2 + Nsupply*phi)-k2*(deltah*k2 + Nsupply*phi))) + k1s*ds*(k2*k2)*(k2*(Hsupply*kdam*k2 + ds*(-Hsupply*kdam + rho*k1p*phi))-dp*(-deltah*k2*k2 + ds*(deltah*k2 + rho*k1p*phi))))) / (k1p*(k1s*k1s)*ds*k2*k2* (-ds+k2) *(k1s*dp*k2*(ds-k2) + k1p*ds*(-dp+k2)*k2)*phi)


'''

###### with QN ##################

#Qn = (9.4e-15*(1/14.0)*1e+9)  #Nitrogen Quota for Pro 
#9.4*10^-15 g N per cell   #Nitrogen fg quota of Pro cells from Bertillison et al 2003
#14 g per Mol N     #N g to M converstion from periodic table
#10^9 ng per g      #conversion from grams to n grams


Nstar = (ks2*ds)/(k2-ds)

Hstar = (((k2*Nstar)/(Nstar + ks1))-(dp))*(1/kdam)


Sstar = (Hsupply - deltah*Hstar)/(phi*Hstar)


Pstar = (Nsupply - (((k2*Nstar)/(Nstar + ks2))*Sstar*Qns) - (rho*Nstar))*((Nstar + ks1)/(k2*Nstar*Qnp))


##### graphing stars 

ax1.axhline(Sstar,color = 'orange', linestyle = "-.",label = 'Sstar')
ax1.axhline(Pstar,color = 'green', linestyle = "-.",label = 'Pstar')
ax2.axhline(Nstar,color = 'purple', linestyle = "-.",label = 'Nstar')
ax3.axhline(Hstar,color = 'red', linestyle = "-.",label = 'Hstar')

ax1.legend(loc = 'best')
ax2.legend(loc = 'best')
ax3.legend(loc = 'best')

fig.savefig('../figures/leaky_calcs_auto',dpi=300)



##############################
#Contour 
##############################



#print(Z)


for (i,SN) in zip(range(SNs.shape[0]),SNs):
    for (j,Sh) in zip(range(Shs.shape[0]),Shs):
        params = [k1p,k1s,k2,dp,ds,kdam,deltah,rho, SN, Sh]
        leaky  = odeint(leak, inits, mtimes, args = (params,))
        Psc = leaky[:,0]
        Ssc = leaky[:,1]
        Nsc = leaky[:,2]
        Hsc = leaky[:,3]
        #print([i,j])
        if np.all(i <= 0.001 for i in Psc[-10:]): 
            Psc[-1] =0.0000000000000000001 
        if Ssc[-1] <0:
            Z[i,j] = 4
        Z[i,j] = Ssc[-1]/(Psc[-1]+Ssc[-1])
        print('Ssc = '+ str(Ssc[-1]) + ' and Psc = '+ str(Psc[-1]))
        print(Z[i,j])
        #if Ps[-1] <=10 and Ss[-1] <=10:
            #Z[i,j] = 4
            #print('dead')


fig3,(ax1) = plt.subplots(sharex = True, sharey = True, figsize = (8,5))


grid = ax1.pcolormesh( Shs, SNs, np.where(Z == 4, np.nan, Z), vmin=np.min(Z), vmax=np.max(Z), cmap = 'summer', shading= 'auto'  )  #'gouraud'
ax1.set(xlabel='Supply hooh')
ax1.set(ylabel='Supply nutrient')

fig3.colorbar(grid, cmap= 'summer',label = 'S / S+P')
#plt.colorbar.set_label('S/P+S')

fig3.savefig('../figures/no_leak_contour_f3_auto',dpi=300)

print('*** SStar =/= dsdt solutions ***')
print('*** Done ***')