#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name:   no_HOOH_Qns_calcs.py 

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
rho =  0.002                  #innate N loss from system 
params = [k1p,k1s,k2,dp,ds,rho]
#params = list(zip(k1s,k2s,kdams))
#muP = (k2 * N /( (k2/k1p) + N) )
#muS = (k2 * N /( (k2/k1s) + N) )


#empty arrays 
P = np.array([])
S  = np.array([])
N = np.array([])
y = [P,S,N]


#function set up for ode int


def comp(y,t,params):
    k1p,k1s,k2,dp,ds,rho = params[0], params[1], params[2],params[3],params[4],params[5]
    P,S,N = y[0],y[1],y[2]
    Nsupply = N0
    Qnp = (9.6e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertillison? 
    Qns = (20.0e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Syn from Bertillison? 
    muP = (k2 * N /( (k2/k1p) + N) )
    muS = (k2 * N /( (k2/k1s) + N) )
    dPdt = P * muP - (dp *P)
    dSdt = S * muS - (ds *S)
    dNdt = Nsupply - (muP * P * Qnp) - (muS * S * Qns)-rho*N
    #print(dNdt)
    return [dPdt,dSdt,dNdt]

#solve ODEs via odeint
competition  = odeint(comp, inits, mtimes, args = (params,))

#redefine where P and  N are in returned matrix from ode int
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
#ax2.plot(mtimes, y = (N0/rho), linestyle = ';', color = 'c', label = "N equilibrium solution")

ax1.semilogy()
ax2.semilogy()

ax1.legend(loc = 'lower right')
ax2.legend(loc = 'lower right')

#plt.show()



##########################################################

# Calculated analytical solutions at equilibrium

##########################################################
'''
# Models
    #muP = (k2 * N /( (k2/k1p) + N) )
    #muS = (k2 * N /( (k2/k1s) + N) )
    dPdt = P * muP - (dp *P)
    dSdt = S * muS - (ds *S)
    dNdt = Nsupply - (muP * P * Qnp) - (muS * S * Qns) -rho*N)            #ignore Qns for equilibrium calculation in matmatica [']]']

#Equilibrium  (3 cases: both P and S die, P wins while S->0, S wins while P->0) 
    N --> Nsupply/rho
    P --> 0
    S --> 0 
    
    #or
    
    N --> dp/muP
    P --> Nsupply/dp - rho/muP
    S -->  0
    
    #or
    
    N --> ds/muS
    P --> 0
    S --> Nsupply/ds - rho/muS
'''

###################################################

#Getting values for analytical solution printing

################################################
    #dPdt = P * muP - (dp *P)
    #dSdt = S * muS - (ds *S)
    
    
#have dsdt,dPdt, and dNdt, but we need muS and muP for all time points
muPs = np.array([])
muSs = np.array([])

'''
# no winner
#N = N0/rho
ax2.hlines(y=N0/rho,xmin=0, xmax=200, linestyle = ':', linewidth=2, color='c', label = 'steady sol N')
'''
# P wins
for p in Ps:
    muP = (p/P)+dp
    np.append(muPs, muP)
    #print(muPs)


ax1.hlines(y = 10**10, xmin=0, xmax=200, linestyle = ':', linewidth=2, color='g', label = 'steady sol P')
ax2.hlines(y = 10**5, xmin=0, xmax=200, linestyle = ':', linewidth=2, color='g', label = 'steady sol N')

#N0/dp - rho/muP
#dp/muP
'''
# S wins
for s in Ss:
    muS = (s/S)+ds
    np.append(muSs, muS)



ax1.hlines(y =N0/ds - rho/muS , xmin=0, xmax=200, linestyle = ':', linewidth=2, color='r', label = 'steady sol S')
ax2.hlines(y =ds/muS , xmin=0, xmax=200, linestyle = ':', linewidth=2, color='c', label = 'steady sol N')
'''


plt.show()

print('no worky')
3


'''
#current issues :
    
    length of single solution wont graph against whole mtimes array
    need N,P,S values for analytical solution calculation but only have dxdt for each'
    muPs array empty so not calculating correctly from odeint output (Ps array) -_-
    
'''