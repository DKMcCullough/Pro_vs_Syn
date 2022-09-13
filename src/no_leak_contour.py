#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name:   no_leak_contour.py 

location: /Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Competitions/Pro_vs_Syn/src

Goal: 
Compete Pro and Syn on one nutrient N with HOOH only effecting Pro 
create contour graph to show Hsupply and N supply rates on x and y axis 
print out color of green (pro wins) or yellow (syn wins) along arrays of H and N supply 

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
H0 = 1490     #nm
inits = (P0,S0,N0,H0)

#parameters

k1p =  0.00002     #Pro alpha
k1s =  0.00001      #Syn alpha 
k2 =  0.88    #Vmax    shared for P and S here
dp = 0.2   #pro delta
ds =  0.2   #syn delta
kdam = 0.5   #hooh mediated damage rate of Pro 
deltah = 0.2       #decay rate of HOOH via Syn 
rho =  0.002                  #innate N loss from system 

params = [k1p,k1s,k2,dp,ds,kdam,deltah,rho]


#empty arrays 
P = np.array([])
S  = np.array([])
N = np.array([])
H = np.array([])
y = [P,S,N,H]


#function set up for ode int


def comp(y,t,params):
    k1p,k1s,k2,dp,ds,kdam,deltah,rho = params[0], params[1], params[2],params[3],params[4], params[5],params[6],params[7]
    P,S,N,H = y[0],y[1],y[2],y[3]
    Nsupply = N0
    S_HOOH = H0
    Qnp = (9.6e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertillison? 
    Qns = (20.0e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Syn from Bertillison? 
    muP = (k2 * N /( (k2/k1p) + N) )
    muS = (k2 * N /( (k2/k1s) + N) )
    dPdt = P * muP - (dp *P) - kdam*H*P
    dSdt = S * muS - (ds *S) 
    dNdt = Nsupply - (muP * P * Qnp) - (muS * S * Qns)-rho*N
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

fig1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True,figsize=(12,7))
fig1.suptitle('Non_leaky HOOH Growth Projections')

#fig2, (ax1, ax2) = plt.subplots(1, 2)
#fig2.suptitle('Non_leaky HOOH Growth Projections')
#plt.subplots_adjust(wspace = 0.3, top = 0.85)

ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'Pro')#' k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn')#' k1 =' + str(k1s))
ax1.set(ylabel='cells per ml')

ax2.plot(mtimes, Ns,linewidth = 3, color = 'purple', label = "Nutrient")
ax2.set(ylabel='Nutrient per ml')

#ax3 = fig.add_subplot(2,1,2)
ax3.plot(mtimes, Hs,linewidth = 3, color = 'red', label = "HOOH")
ax3.set(xlabel='Time (days)', ylabel='HOOH per ml')


ax1.semilogy()
ax2.semilogy()
ax3.semilogy()

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

'''

#mumaxS or P = Qnp*k2 or Qns*k2 respectively

#Swins
Qnp = (9.6e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertillison?  #giving N in micro units? 
Qns = (20.0e-15*(1/(14.0))*1e+9)

Hsupply = H0
Nsupply = N0

Nstar = -(((deltah*dp)+(Hsupply*kdam))*k2)/((k1p)*(deltah*dp)+(Hsupply*kdam)-(deltah*k2))

Sstar = -((k1s)*((Nsupply*k1s)*(ds-k2)+(rho*ds*k2))/((k2)*(ds*k2)+((k1s*k1s)*(-ds + k2))))

Hstar = Hsupply/deltah

#graphing equilibrium values 

ax1.axhline(Sstar,color = 'orange', linestyle = "-.",label = 'Sstar')
ax2.axhline(Nstar,color = 'purple', linestyle = "-.",label = 'Nstar')
ax3.axhline(Hstar,color = 'red', linestyle = "-.",label = 'Hstar')


ax1.legend(loc = 'best')
ax2.legend(loc = 'best')
ax3.legend(loc = 'best')

plt.show()

fig1.savefig('../figures/no_leak_contour_auto',dpi=300)



#NSTAR NOT showing up
#Sstar value off by 10*

##############################
#Contour 
##############################

Shs = np.linspace(0, 800, num=50)
SNs = np.linspace(1000, 100000, num = 50)
Shs = Shs.astype(int)
SNs = SNs.astype(int)

Z = np.random.rand(50,50)
'''
for (i, SN) in zip(range(len(SNs),SNs)): 
                   for (j,Sh) in zip(range(len(Shs),Shs)):
                                     Z[i,j] = Sr*Sh

'''

fig3,(ax1) = plt.subplots(sharex = True, sharey = True, figsize = (8,5))

ax1.pcolor( Z, vmin=np.min(Z), vmax=np.max(Z), cmap = 'summer', shading='auto')






