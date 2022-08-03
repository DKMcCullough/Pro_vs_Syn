#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name:   no_leak_calcs.py 

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
ndays = 200
mtimes = np.linspace(0,ndays,int(ndays/step))

#initial values 
P0 = 1e4
S0 = 1e4
N0 = 1.0e4        #nM 
H0 = 1500     #nm
inits = (P0,S0,N0,H0)

#parameters

k1p =  0.00002     #Pro alpha
k1s =  0.00001      #Syn alpha 
k2 =  0.88    #Vmax    shared for P and S here
dp = 0.2   #pro delta
ds =  0.2   #syn delta
kdam = 0.5   #hooh mediated damage rate of Pro 
deltah = 0.2       #decay rate of HOOH via Syn 
params = [k1p,k1s,k2,dp,ds,kdam,delta]
#params = list(zip(k1s,k2s,kdams))
Qnp = (9.6e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertillison? 
Qns = (20.0e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Syn from Bertillison? 
muP = (k2 * N /( (k2/k1p) + N) )
muS = (k2 * N /( (k2/k1s) + N) )


#empty arrays 
P = np.array([])
S  = np.array([])
N = np.array([])
H = np.array([])
y = [P,S,N,H]


#function set up for ode int


def comp(y,t,params):
    k1p,k1s,k2,dp,ds,kdam,deltah = params[0], params[1], params[2],params[3],params[4], params[5],params[6]
    P,S,N,H = y[0],y[1],y[2],y[3]
    Nsupply = N0
    S_HOOH = H0
    Qnp = (9.6e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertillison? 
    Qns = (20.0e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Syn from Bertillison? 
    muP = (k2 * N /( (k2/k1p) + N) )
    muS = (k2 * N /( (k2/k1s) + N) )
    dPdt = P * muP - (dp *P) - kdam*H
    dSdt = S * muS - (ds *S) 
    dNdt = Nsupply - (muP * P * Qnp) - (muS * S * Qns)
    dHdt = S_HOOH - deltah*H
    return [dPdt,dSdt,dNdt,dHdt]

#solve ODEs via odeint
competition  = odeint(comp, inits, mtimes, args = (params,))

#redefine where P and N are in returned matrix from ode int
Ps = competition[:,0]
Ss = competition[:,1]
Ns = competition[:,2]
Hs = competition[:,3]



#####################################

#  Graphing

#####################################



fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Growth Competition Projections')
plt.subplots_adjust(wspace = 0.3, top = 0.85)


ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'Pro k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn k1 =' + str(k1s))
ax1.set(xlabel='Time (days)', ylabel='cells per ml')

ax2.plot(mtimes, Ns, label = "Nutrient Concentration over time")
ax2.set(xlabel='Time (days)', ylabel='Nutrient concentration')

ax1.semilogy()
ax2.semilogy()

ax1.legend(loc = 'lower right')

#plt.show()




##########################################################

# Calculated analytical solutions at equilibrium

##########################################################

'''
# Models
    #muP = (k2 * N /( (k2/k1p) + N) )
    #muS = (k2 * N /( (k2/k1s) + N) )
    dPdt = P * muP - (dp *P) - kdam*H 
    dSdt = S * muS - (ds *S) 
    dNdt = Nsupply - (muP * P * Qnp) - (muS * S * Qns)-rho*N
    dHdt = S_HOOH - deltah*H
    
    
#Analytical solutions at equilibrium (dx/dt equations above set to 0 and solved)

#Equilibrium  (3 cases: both P and S die, P wins while S->0, S wins while P->0))
    N --> Nsupply/rho
    P --> 0
    S --> 0
    H --> Hsupply/deltah
    
    #or
    
    N --> (deltah*dp + Hsupply*kdam)/deltah*muP
    P --> ((deltah*rho*dp)+(Hsupply*rho*kdam)-(Nsupply*deltah*muP))/((deltah*dp*muP) + (Hsupply*kdam*muP))
    S --> 0
    H --> Hsupply/deltah
    
    #or
    
    N --> ds/muS
    P --> 0
    S --> (Nsupply/ds) - (rho/muS)
    H --> Hsupply/deltah
'''

















