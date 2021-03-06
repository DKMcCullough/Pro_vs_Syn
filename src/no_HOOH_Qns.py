#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name:   no_leak.py 

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
inits = (P0,S0,N0)

#parameters

k1p =  0.000003     #Pro alpha
k1s =  0.000002      #Syn alpha 
k2 =  0.88    #Vmax    shared for P and S here
dp = 0.2   #pro delta
ds =  0.2   #syn delta
params = [k1p,k1s,k2,dp,ds]
#params = list(zip(k1s,k2s,kdams))



#empty arrays 
P = np.array([])
S  = np.array([])
N = np.array([])
y = [P,S,N]


#function set up for ode int


def comp(y,t,params):
    k1p,k1s,k2,dp,ds = params[0], params[1], params[2],params[3],params[4]
    P,S,N = y[0],y[1],y[2]
    Nsupply = N0
    Qnp = (9.6e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertillison? 
    Qns = (20.0e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Syn from Bertillison? 
    muP = (k2 * N /( (k2/k1p) + N) )
    muS = (k2 * N /( (k2/k1s) + N) )
    dPdt = P * muP - (dp *P)
    dSdt = S * muS - (ds *S)
    dNdt = Nsupply - (muP * P * Qnp) - (muS * S * Qns)
    print(dNdt)
    return [dPdt,dSdt,dNdt]

#solve ODEs via odeint
competition  = odeint(comp, inits, mtimes, args = (params,))

#redefine where P and N are in returned matrix from ode int
Ps = competition[:,0]
Ss = competition[:,1]
Ns = competition[:,2]



#####################################

#  Graphing

#####################################



fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Growth Competition Projections')
plt.subplots_adjust(wspace = 0.3, top = 0.85)


ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'Pro') #label = 'Pro k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn') #label = 'Syn k1 =' + str(k1s))
ax1.set(xlabel='Time (days)', ylabel='cells per ml')

ax2.plot(mtimes, Ns, label = "Nutrient Concentration over time")
ax2.set(xlabel='Time (days)', ylabel='Nutrient concentration')

ax1.semilogy()
ax2.semilogy()

ax1.legend(loc = 'lower right')

plt.show()








