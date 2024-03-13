#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
name:   leaky_contour.py 

location: /Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Competitions/Pro_vs_Syn/src


@author: dkm
'''




from functions import * 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import *
from scipy.integrate import odeint
import sys

##################################################3
# parameter and variable Set UP 
#############################################


step = 0.001
ndays = 50
mtimes = np.linspace(0,ndays,int(ndays/step))


P = np.array([])
S  = np.array([])
N = np.array([])
H = np.array([])

y = [P,S,N,H]

#initial values to be used for odeint start 

P0 = 1e4
S0 = 1e4
N0 = 1.0e1        #nM 
H0 = 1     #nM

inits = (P0,S0,N0,H0)


#####################################

#  P dominated dynamic model 

#####################################
#params for P to Win 
Sh = 10
SN = 100000
params = [ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh]

#get equilibria with function 
[Nstar, Pstar, Sstar, Hstar] = Coexist(params)[0,1,2,3]

#run model 

competition  = odeint(leak, inits, mtimes, args = (params,))


#grab values to graph 
Ps = competition[:,0]
Ss = competition[:,1]
Ns = competition[:,2]
Hs = competition[:,3]

#graaph dynamics where P dominates 
fig1, (ax1, ax2,ax3) = plt.subplots(3,1, sharex=True, figsize=(9,5),dpi = 300)
fig1.suptitle('Growth Competition Projections in Sn'+ str(SN)+' Sh '+str(Sh))
plt.subplots_adjust(wspace = 0.5, top = 0.9,bottom = 0.1)
ax1.set(xlabel='Time (days)', ylabel='cells per ml')
ax2.set(xlabel='Time (days)', ylabel='Nutrient [ ]')
ax3.set(xlabel='Time (days)', ylabel='HOOH [ ]')


ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'Pro ') #'k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn ') #'k1 =' + str(k1s))
ax2.plot(mtimes, Ns, linewidth = 3, color = 'purple', label = "Nutrient Concentration ")
ax3.plot(mtimes, Hs,  linewidth = 3, color = 'red', label = "HOOH concentration ")

##### graphing stars equations ############### 

ax1.axhline(Sstar,color = 'orangered', linestyle = "-.",label = 'S*')
ax1.axhline(Pstar,color = 'darkgreen', linestyle = ":",label = 'P*')
ax2.axhline(Nstar,color = 'purple', linestyle = "-.",label = 'N*')
ax3.axhline(Hstar,color = 'red', linestyle = "-.",label = 'H*')

ax1.semilogy()
ax2.semilogy()
ax3.semilogy()

fig1.savefig('../figures/leaky_dynamics_P',dpi=300)
plt.show()


#####################################

#  S dominated dynamic model 

#####################################
#params for S to Win 
Sh = 100
SN = 200
params = [ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh]

#get equilibria with function 
[Nstar, Pstar, Sstar, Hstar] = Coexist(params)[0,1,2,3]

#run model 

competition  = odeint(leak, inits, mtimes, args = (params,))


#grab values to graph 
Ps = competition[:,0]
Ss = competition[:,1]
Ns = competition[:,2]
Hs = competition[:,3]



#graaph dynamics where S dominates 
fig1, (ax1, ax2,ax3) = plt.subplots(3,1, sharex=True, figsize=(9,5),dpi = 300)
fig1.suptitle('Growth Competition Projections in Sn'+ str(SN)+' Sh '+str(Sh))
plt.subplots_adjust(wspace = 0.5, top = 0.9,bottom = 0.1)
ax1.set(xlabel='Time (days)', ylabel='cells per ml')
ax2.set(xlabel='Time (days)', ylabel='Nutrient [ ]')
ax3.set(xlabel='Time (days)', ylabel='HOOH [ ]')


ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'Pro ') #'k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn ') #'k1 =' + str(k1s))
ax2.plot(mtimes, Ns, linewidth = 3, color = 'purple', label = "Nutrient Concentration ")
ax3.plot(mtimes, Hs,  linewidth = 3, color = 'red', label = "HOOH concentration ")

##### graphing stars equations ############### 

ax1.axhline(Sstar,color = 'orangered', linestyle = "-.",label = 'S*')
ax1.axhline(Pstar,color = 'darkgreen', linestyle = ":",label = 'P*')
ax2.axhline(Nstar,color = 'purple', linestyle = "-.",label = 'Nstar')

ax2.axhline(Nstars,color = 'pink', linestyle = "-.",label = 'Ns*')
ax2.axhline(Nstarp,color = 'magenta', linestyle = ":",label = 'Np*')
ax3.axhline(Hstar,color = 'red', linestyle = "-.",label = 'H*')

ax1.semilogy()
ax2.semilogy()
ax3.semilogy()

fig1.savefig('../figures/leaky_dynamics_S',dpi=300)


##############################
#Contour and model calculations 
##############################


SNs = np.linspace(0, 500, num = 140)
Shs = np.linspace(0, 1000, num = 40)
Z = np.zeros((int(SNs.shape[0]),int(Shs.shape[0])),float)

for (i,SN) in zip(range(SNs.shape[0]),SNs):
    for (j,Sh) in zip(range(Shs.shape[0]),Shs):
        params = [ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho, SN, Sh]
        leaky  = odeint(leak, inits, mtimes, args = (params,))
        #get equilibria with function 
        [Nstar, Pstar, Sstar, Hstar] = Coexist(params)[0,1,2,3]
        Psc = leaky[:,0]
        Ssc = leaky[:,1]
        Nsc = leaky[:,2]
        Hsc = leaky[:,3]
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


#######################################
# Graphing Cotour plots from 
######################################

fig3,(ax1) = plt.subplots(sharex = True, sharey = True, dpi = 300)
fig3.suptitle('HOOH Contour with Phi')# + str(phi))
# figsize = (8,5)
grid = ax1.pcolormesh( Shs, SNs, np.where(Z == -1, np.nan, Z), vmin=0, vmax=np.max(Z), cmap = 'summer', shading= 'auto'  )  #'gouraud'
ax1.set(xlabel='Supply HOOH')
ax1.set(ylabel='Supply N')


#ax1.axhline((rho*Nstar),color = 'purple', linestyle = "-.",label = 'Nstar')
#ax1.axvline((vHline),color = 'c', linestyle = "-.",label = 'H ')
#ax1.axvline((Hstar*(deltah+(phi*Pstar))),color = 'b', linestyle = "-.",label = 'H cutoff?')

fig3.colorbar(grid, cmap= 'summer',label = 'S / S+P')

fig3.savefig('../figures/leaky_contour',dpi=300)

print('') 

print('*** Done ***')