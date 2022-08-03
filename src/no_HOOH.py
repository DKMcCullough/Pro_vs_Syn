#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name:   no_HOOH.py 

location: /Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Competitions/Pro_vs_Syn/src

Goal: 
    Compete Pro and Syn on one nutrient N
    Pro to have higher affinity for N (as well as usage perhaps) than Syn

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
#P_star0 = 1e4
#N_star0 = 1e4
inits = (P0,S0,N0) #P_star0,N_star0)

#parameters

k1p =  0.000003     #Pro alpha
k1s =  0.000002      #Syn alpha 
k2 =  0.88    #Vmax    shared for P and S here
dp = 0.2   #pro delta
ds =  0.2   #syn delta
params = [k1p,k1s,k2,dp,ds]
#params = list(zip(k1s,k2s,kdams))
#muP = (k2 * N /( (k2/k1p) + N) )
#muS = (k2 * N /( (k2/k1s) + N) )


#empty arrays 
P = np.array([])
S  = np.array([])
N = np.array([])
#P_star = np.array([])
#N_star = np.array([])
y = [P,S,N]    #P_star,N_star]
muPs = np.array([])


#function set up for ode int


def comp(y,t,params):
    k1p,k1s,k2,dp,ds = params[0], params[1], params[2],params[3],params[4]
    P,S,N = y[0],y[1],y[2] #,y[3],y[4]
    Nsupply = N0
    Qn = (9.6e-15*(1/(14.0))*1e+9)   #use Qn of Pro for both 
    muP = (k2 * N /( (k2/k1p) + N) )
    muS = (k2 * N /( (k2/k1s) + N) )
    dPdt = P * muP - (dp *P)
    dSdt = S * muS - (ds *S)
    dNdt = Nsupply - (muP * P * Qn) - (muS * S * Qn)
    #P_star = Nsupply/muP
    #N_star = dp/muP
    np.append(muPs,muP)
   
    return [dPdt,dSdt,dNdt] #N_star,P_star]

#solve ODEs via odeint
competition  = odeint(comp, inits, mtimes, args = (params,))

#redefine where P and N are in returned matrix from ode int
Ps = competition[:,0]
Ss = competition[:,1]
Ns = competition[:,2]
#N_stars = competition[:,3]
#P_stars = competition[:,4]



#####################################

#  Graphing

#####################################



fig, (ax1, ax2) = plt.subplots(1, 2,figsize=[7,5])
fig.suptitle('Growth Projections abscent HOOH')
plt.subplots_adjust(wspace = 0.3, top = 0.9)


ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'Pro') #label = 'Pro k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn') #label = 'Syn k1 =' + str(k1s))
ax1.set(xlabel='Time (days)', ylabel='cells  L$^{-1}$')

ax2.plot(mtimes, Ns, label = "Nutrient Concentration over time")
ax2.set(xlabel='Time (days)', ylabel='Nutrient concentration')

ax1.semilogy()
ax2.semilogy()

ax1.legend(loc = 'lower right')
ax1.set_ylim(bottom = 0.1, top =10000000000000)
#plt.show()



##########################################################

# Calculated analytical solutions at equilibrium

##########################################################
'''
# Models
    muP = (k2 * N /( (k2/k1p) + N) )
    muS = (k2 * N /( (k2/k1s) + N) )
    dPdt = P * muP - (dp *P)
    dSdt = S * muS - (ds *S)
    dNdt = Nsupply - (muP * P * Qnp) - (muS * S * Qns)   #ignore Qns for equilibrium calculation in matmatica [']]']

#Equilibrium  (2 cases: P wins while S->0, S wins while P->0) 
    N --> dp/muP
    P --> Nsupply/muP
    S --> 0 
    
    #or
    
    N --> ds/muS
    P --> 0
    S --> Nsupply/muS 
    
'''


'''

for P, S, N in zip(Ps, Ss, Ns): 
    if p-1 < p: # & s-1 > s: 
        N = dp/muP
        P = Nsupply/muP
        S = 0
    else:
        N = ds/muP
        P = 0
        S = Nsupply/muS
'''

#ax1.plot(mtimes, N0/muPs , linestyle = ':', color = 'g', label = 'P* projection')
#ax2.plot(mtimes, N_stars , linestyle = ':', color = 'b', label = 'P* projection')
# print out dotted P* and N* lines with starting params onto time dependantcompetition projections graph
#plt.show()

fig.savefig('../figures/no_hooh_auto',dpi=300)

