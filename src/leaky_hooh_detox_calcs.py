#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name:   leaky_hooh_detox_calcs.py 

location: /Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Competitions/Pro_vs_Syn/src

Goal: 
    Compete Pro and Syn on one nutrient N
    Pro to have higher affinity for N (as well as usage perhaps) than Syn
    Pro dies from HOOH and SYn detoxes only for itselfc

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

muP = (k2 * N /( (k2/k1p) + N) )
muS = (k2 * N /( (k2/k1s) + N) )
dPdt = P * muP 
dSdt = S * muS  
dNdt = supply - (muP * P) - (muS *S)


'''
#####################
# Set UP
######################


step = 0.01
ndays = 800
mtimes = np.linspace(0,ndays,int(ndays/step))

#initial values 
P0 = 1e4
S0 = 1e4
N0 = 1.0e5        #nM 
H0 = 1485     #nM
inits = (P0,S0,N0,H0)

#parameters

k1p =  0.00002     #Pro alpha
k1s =  0.00001      #Syn alpha 
k2 =  0.88    #Vmax    shared for P and S here
dp = 0.2   #pro delta
ds =  0.2   #syn delta
kdamp = 0.005   #hooh mediated damage rate of Pro 
kdams = 0.0004   #hooh mediated damage rate of Syn 
deltah = 0.002       #decay rate of HOOH via Syn 
phi = 0.0002    #0007  #detoxification-based decay of HOOH via Syn in this case
params = [k1p,k1s,k2,dp,ds,kdamp,kdams,delta,phi]
#params = list(zip(k1s,k2s,kdams))

Nsupply = N0
S_HOOH = H0


#empty arrays 
P = np.array([])
S  = np.array([])
N = np.array([])
H = np.array([])
y = [P,S,N,H]


#function set up for ode int


def comp(y,t,params):
    k1p,k1s,k2,dp,ds,kdamp,kdams,deltah,phi = params[0], params[1], params[2],params[3],params[4], params[5],params[6], params[7],params[8]
    P,S,N,H = y[0],y[1],y[2],y[3]
    Qnp = (9.6e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertillison? 
    Qns = (20.0e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Syn from Bertillison? 
    muP = (k2 * N /( (k2/k1p) + N) )
    muS = (k2 * N /( (k2/k1s) + N) )
    dPdt = P * muP - (dp *P) - kdamp*H*P
    dSdt = S * muS - (ds *S) - kdams*H*S      
    dNdt =  Nsupply - (muP * P * Qnp) - (muS * S * Qns)    
    dHdt = S_HOOH - deltah*H  -phi*S*H  #phi *S*H with phi being Scell-specific detox rate
    if t<0.05:
        print(muP,muS)
    return [dPdt,dSdt,dNdt,dHdt]

#solve ODEs via odeint
competition  = odeint(comp, inits, mtimes, args = (params,))

#redefine where P and N are in returned matrix from ode int
Ps = competition[:,0]
Ss = competition[:,1]
Ns = competition[:,2]
Hs = competition[:,3]
 #need to set a limit on y lower axis bound bc cutting it off at graph just leaves a blank graph

#####################################

#  Graphing

#####################################


#fig,ax1 = plt.subplots()
fig, (ax1, ax2,ax3) = plt.subplots(1,3,figsize=(9,4.5))
fig.suptitle('Growth Competition Projections')
plt.subplots_adjust(wspace = 0.5, top = 0.85,bottom = 0.1)


ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'Pro k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn k1 =' + str(k1s))
ax1.set(xlabel='Time (days)', ylabel='cells per ml')
#ax1.set_ylim(bottom = -20)


ax2.plot(mtimes, Ns, label = "Nutrient Concentration over time")
ax2.set(xlabel='Time (days)', ylabel='Nutrient concentration')

ax3.plot(mtimes, Hs, label = "HOOH concentration ")
ax3.set(xlabel='Time (days)', ylabel='HOOH concentration')

ax1.semilogy()
ax2.semilogy()
ax3.semilogy()

#ax1.legend(loc = 'lower right')

'''
pnames = ('k1p','k1s','k2','dp','ds','kdamp','kdams','delta','phi')
subplot =fig.add_subplot()
subplot.table([params], colLabels=(pnames) ,loc='lower center')

'''


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
    N --> Nsupply/rho
    P --> 0
    S --> 0
    H --> Hsupply/deltah
    
    #or
    
    N --> ((deltah*dp)+(Hsupply*kdamp))/deltah*muP
    P --> - ((deltah*rho*dp)+(Hsupply*rho*kdamp)-(Nsupply*deltah*muP))/((deltah*dp*muP)+(Hsupply*kdamp*muP))
    S --> 0
    H --> Hsupply/deltah
    
    #or
    
    N --> ((ds*kdamp-dp*kdams)/(-kdams*muP +kdamp*muS))
    P --> ((((ds^2)*kdamp*muP*(deltah*muS - rho*phi))+ (ds*(((-dp)*(kdams*muP + kdamp*muS)*(deltah*muS - rho*phi)) -((kdams*muP-kdamp*muS)*((-Hsupply)*kdamp*muS + Nsupply*muP*phi)))) + ((dp*muS)*(((-Hsupply)*(kdams^2)*muP)-(Nsupply*kdamp*muS*phi)+kdams*(Hsupply*kdamp*muS + Nsupply*muP*phi + dp*(deltah_muP - rho*phi)))))/((ds*kdamp-dp*kdams)*muP*(ds*muP-dp*muS)*phi))
    S --> (((-deltah*ds*muP)-(Hsupply*kdams*muP)+(deltah*muP + Hsupply*kdamp)*muS))/((dS*muP - dp*muS)*phi))
    H --> ((ds*muP-dp*muS)/(-kdams*muP +kdamp*muS))
    
    #or
    
    N --> 
    P --> 
    S --> 
    H --> 
    
    #or
    
    N --> 
    P --> 
    S --> 
    H --> 
    
'''



