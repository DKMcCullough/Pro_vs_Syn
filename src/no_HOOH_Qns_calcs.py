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
ndays = 300
mtimes = np.linspace(0,ndays,int(ndays/step))
Shs = np.linspace(0, 200, num =10)
SNs = np.linspace(0, 1000, num = 10)
Z = np.zeros((int(SNs.shape[0]),int(Shs.shape[0])),int)

#initial values 
P0 = 1e4
S0 = 1e4
N0 = 1.0e4        #nM 
inits = (P0,S0,N0)


#initial values 
P0 = 1e4
S0 = 1e4
N0 = 1e5        #nM 
H0 = 1    #nm
inits = (P0,S0,N0,H0)



#parameters
Qnp = 1#(9.4e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertillison?  #giving N in micro units? 
Qns = 1#(20.0e-15*(1/(14.0))*1e+9)

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


def noH(y,t,params):
    k1p,k1s,k2,dp,ds,rho = params[0], params[1], params[2],params[3],params[4],params[5]
    P,S,N = y[0],y[1],y[2]
    Nsupply = N0
    Qnp = (9.4e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertillison?  #giving N in micro units? 
    Qns = (20.0e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Syn from Bertillison? 
    muP = (k2 * N /( (k2/k1p) + N) )
    muS = (k2 * N /( (k2/k1s) + N) )
    dPdt = P * muP - (dp *P)
    dSdt = S * muS - (ds *S)
    dNdt = Nsupply - (muP * P * Qnp) - (muS * S * Qns)-rho*N
    #print(dNdt)
    return [dPdt,dSdt,dNdt]

#solve ODEs via odeint
noHooh  = odeint(noH, inits, mtimes, args = (params,))

#redefine where P and  N are in returned matrix from ode int
Ps = noHooh[:,0]
Ss = noHooh[:,1]
Ns = noHooh[:,2]



#####################################

#  Graphing

#####################################



fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Growth Competition Projections')
plt.subplots_adjust(wspace = 0.3, top = 0.85)


ax1.plot(mtimes, np.clip(Ps,0.001,10e13) , linewidth = 3, color = 'g', label = 'Pro') #label = 'Pro k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn') #label = 'Syn k1 =' + str(k1s))
ax1.set(xlabel='Time (days)', ylabel='cells per ml')


ax2.plot(mtimes, Ns, linewidth = 3, color = 'purple', label = "Nutrient")
ax2.set(xlabel='Time (days)', ylabel='Nutrient concentration')
#ax2.plot(mtimes, y = (N0/rho), linestyle = ';', color = 'c', label = "N equilibrium solution")

ax1.semilogy()
ax2.semilogy()




##########################################################

# Calculated analytical solutions at equilibrium

##########################################################


#Pwins
Qnp = (9.4e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertillison?  #giving N in micro units? 
Qns = (20.0e-15*(1/(14.0))*1e+9)
    
Nstar = (dp*(k2/k1p))/(k2-dp)
Pstar = (N0-rho*Nstar)/(Qnp*(k2*Nstar)/(Nstar + (k2/k1p)))

ax1.axhline(Pstar, color = 'green', linestyle = "-.",label = 'Pstar')
ax2.axhline(Nstar, color = 'purple', linestyle = "-.",label = 'Nstar')


'''
#M solutions retry
# Models
    #muP = (k2 * N /( (k2/k1p) + N) )
    #muS = (k2 * N /( (k2/k1s) + N) )
    dPdt = P * muP - (dp *P)
    dSdt = S * muS - (ds *S)
    dNdt = Nsupply - (muP * P * Qnp) - (muS * S * Qns) -rho*N)            #ignore Qns for equilibrium calculation in matmatica [']]']

#Equilibrium  (3 cases: both P and S die, P wins while S->0, S wins while P->0) 
    N --> 
    P --> 0
    S --> 0 
    
    #or
    
    N --> (dp*(k2/k1p))/(k2-dp)
    P --> (N0-rho*Nstar)/(Qnp*(k2*Nstar)/(Nstar + (k2/k1p)))
    S -->  0
    
    #or
    
    N --> 
    P --> 0
    S --> 
'''
ax1.legend(loc = 'lower right')
ax2.legend(loc = 'lower right')


#plt.xlim([0, 50])
plt.ylim([10,10e8])

plt.show()

fig.savefig('../figures/no_HOOH_Qns_calcs_auto',dpi=300)


for (i,SN) in zip(range(SNs.shape[0]),SNs):
    for (j,Sh) in zip(range(Shs.shape[0]),Shs):
        params = [ksp,kss,k2,dp,ds,kdam,deltah,rho, SN, Sh]
        nonleaky  = odeint(nleak, inits, mtimes, args = (params,))
        Psc = nonleaky[:,0]
        Ssc = nonleaky[:,1]
        Nsc = nonleaky[:,2]
        Hcs = nonleaky[:,3]
        if (i == 0) and (j == 3):
            Hsupply = Shs[j]
            Nsupply = SNs[i]
            Ps = nonleaky[:,0]
            Ss = nonleaky[:,1]
            Ns = nonleaky[:,2]
            Hs = nonleaky[:,3]
            #need to shoose star values here for dynamic print (N, P or S, H)
        Z[i,j] = Ssc[-1]/(Psc[-1]+Ssc[-1])
        if Z[i,j] < 0.5:
            print('SN = '+ str(i),', Sh = ' + str(j))
            print('Psc[-1] = '+ str(Psc[-1]), ', Ssc[-1] = '+ str(Ssc[-1]))
        #if ([Psc[g] < Psc[g-1] for g in Psc[:]]) and ([Ssc[h] < Ssc[h-1] for h in Ssc[:]]) : 
            #Z[i,j] = -1
        if np.all([g <= 1e-3 for g in Psc[-10:]]) and np.all([h <= 1e-3 for h in Ssc[-10:]]) : 
            Z[i,j] = -1


