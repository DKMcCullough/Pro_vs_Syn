#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
name:   leaky_contour.py 

location: /Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Competitions/Pro_vs_Syn/src


@author: dkm
'''




from functions_chap3 import * 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import *
from scipy.integrate import odeint
import sys

################################################################################
# parameter and variable Set UP 
##################################################################################


step = 0.001
ndays = 400
mtimes = np.linspace(0,ndays,int(ndays/step))


P = np.array([])
S  = np.array([])
N = np.array([])
H = np.array([])

y = [P,S,N,H]

#initial values to be used for odeint start 

P0 = 1e3
S0 = 1e3
N0 = 5.0e5        #nM 
H0 = 1     #nM

inits = (P0,S0,N0,H0)
'''

##########################################################################

#  P dominated dynamic model 

##########################################################################
#params for P to Win 
Sh = 0
SN  = 5e4

params = [k1p,k1s,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh]

#get equilibria with function 
Nstar, Pstar, Sstar, Hstar = Pwins(params)

#run model 

competition  = odeint(leak, inits, mtimes, args = (params,))


#grab values to graph 
Ps = competition[:,0]
Ss = competition[:,1]
Ns = competition[:,2]
Hs = competition[:,3]


#graaph dynamics where P dominates 
fig1, (ax1, ax2) = plt.subplots(2,1, sharex=True, figsize=(7,5),dpi = 300)
fig1.suptitle('Growth Competition Projections in Sn '+ str(SN)+' Sh '+str(Sh))
plt.subplots_adjust(right=0.95, wspace = 0.45, left = 0.10, hspace = 0.15, bottom = 0.10)
ax1.set(xlabel='Time (days)', ylabel='Cells per ml')
ax2.set(xlabel='Time (days)', ylabel='Nutrient per ml')


ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'Pro abundance') #'k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn abundance') #'k1 =' + str(k1s))
ax2.plot(mtimes, Ns, linewidth = 3, color = 'purple', label = "Nutrient Concentration ")

##### graphing stars equations ############### 

ax1.axhline(Sstar,color = 'orangered', linestyle = "-.",label = 'S*')
ax1.axhline(Pstar,color = 'darkgreen', linestyle = ":",label = 'P*')
ax2.axhline(Nstar,color = 'purple', linestyle = "-.",label = 'N*')


ax1.semilogy()
ax2.semilogy()
ax1.set_ylim(ymin = 100, ymax = 1e+8)
ax2.set_ylim(ymin = 1, ymax = 1e+5)



fig1.savefig('../figures/leaky_dynamics_P',dpi=300)
plt.show()


##########################################################################

#  S dominated dynamic model 

##########################################################################
#params for S to Win 
Sh = 1200
SN  = 7e5
params = [k1p,k1s,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh]

#get equilibria with function 
Nstar, Pstar, Sstar, Hstar = Swins(params)

#run model 

competition  = odeint(leak, inits, mtimes, args = (params,))


#grab values to graph 
Ps = competition[:,0]
Ss = competition[:,1]
Ns = competition[:,2]
Hs = competition[:,3]


#####################################
#graaph dynamics where S dominates 
#####################################
fig1, (ax1, ax2,ax3) = plt.subplots(3,1, sharex=True, figsize=(7,6),dpi = 300)
fig1.suptitle('Growth Competition Projections in SN '+ str(SN)+' Sh '+str(Sh))
plt.subplots_adjust(right=0.95, wspace = 0.45, left = 0.10, hspace = 0.15, bottom = 0.10)
ax1.set(xlabel='Time (days)', ylabel='cells per ml')
ax2.set(xlabel='Time (days)', ylabel='Nutrient per ml')
ax3.set(xlabel='Time (days)', ylabel='H$_2$O$_2$ per ml')


ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'Pro abundance') #'k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn abundance ') #'k1 =' + str(k1s))
ax2.plot(mtimes, Ns, linewidth = 3, color = 'purple', label = "Nutrient Concentration ")
ax3.plot(mtimes, Hs,  linewidth = 3, color = 'red', label = "H$_2$O$_2$ concentration ")

#####################################
##### graphing stars equations ############### 
#####################################
ax1.axhline(Sstar,color = 'orangered', linestyle = "-.",label = 'S*')
ax1.axhline(Pstar,color = 'darkgreen', linestyle = ":",label = 'P*')
ax2.axhline(Nstar,color = 'purple', linestyle = "-.",label = 'N*')
ax3.axhline(Hstar,color = 'red', linestyle = "-.",label = 'H*')

ax1.semilogy()
ax2.semilogy()
ax3.semilogy()

ax1.set_ylim(ymin = 100, ymax = 1e+8)
ax2.set_ylim(ymin = 1, ymax = 1e+5)
ax3.set_ylim(ymin = 10, ymax = 1000)

fig1.savefig('../figures/leaky_dynamics_S',dpi=300)


###############################################################################################################

#  dynamics for coexistance (gotten from contour) 

###############################################################################################################
#params for coexists
Sh = 800
SN = 8e6

params = [k1p,k1s,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh]


#run model 

competition  = odeint(leak, inits, mtimes, args = (params,))
#get equilibria with function 
Nstar, Pstar, Sstar, Hstar = Coexist(params)


#grab values to graph 
Ps = competition[:,0]
Ss = competition[:,1]
Ns = competition[:,2]
Hs = competition[:,3]



#graph coexistence dynamics 
fig1, (ax1, ax2,ax3) = plt.subplots(3,1, sharex=True, figsize=(7,6),dpi = 300)
fig1.suptitle('Growth Competition Projections in SN '+ str(SN)+' Sh '+str(Sh))
plt.subplots_adjust(right=0.95, wspace = 0.45, left = 0.10, hspace = 0.15, bottom = 0.10)
ax1.set(xlabel='Time (days)', ylabel='cells per ml')
ax2.set(xlabel='Time (days)', ylabel='Nutrient per ml')
ax3.set(xlabel='Time (days)', ylabel='H$_2$O$_2$ per ml')


ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'Pro abundance') #'k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn abundance') #'k1 =' + str(k1s))
ax2.plot(mtimes, Ns, linewidth = 3, color = 'purple', label = "Nutrient Concentration ")
ax3.plot(mtimes, Hs,  linewidth = 3, color = 'red', label = "H$_2$O$_2$ concentration ")

##### graphing stars equations ############### 

ax1.axhline(Sstar,color = 'orangered', linestyle = "-.",label = 'S*')
ax1.axhline(Pstar,color = 'darkgreen', linestyle = ":",label = 'P*')
ax2.axhline(Nstar,color = 'purple', linestyle = "-.",label = 'N*')
ax3.axhline(Hstar,color = 'red', linestyle = "-.",label = 'H*')

ax1.semilogy()
ax2.semilogy()
ax3.semilogy()

ax1.set_ylim(ymin = 100, ymax = 1e+8)
ax2.set_ylim(ymin = 1, ymax = 1e+5)
ax3.set_ylim(ymin = 10, ymax = 1000)

fig1.savefig('../figures/leaky_dynamics_coexist',dpi=300)

'''
########################################################################################################
#dynamics contour and model calculations 
########################################################################################################


SNs = np.linspace(1.0e4, 1.0e7, num = 100)
Shs = np.linspace(0, 2000, num = 100)
Z = np.zeros((int(SNs.shape[0]),int(Shs.shape[0])),float)
'''
#dynamics in non leaky  - only phio froom bacteria - none from sYn
for (i,SN) in zip(range(SNs.shape[0]),SNs):
    for (j,Sh) in zip(range(Shs.shape[0]),Shs):
        #phi = 0
        params = [k1p,k1s,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh]
        leaky  = odeint(leak, inits, mtimes, args = (params,))
        #get equilibria with function 
        Psc = leaky[:,0]
        Ssc = leaky[:,1]
        Nsc = leaky[:,2]
        Hsc = leaky[:,3]
        Ssc_av = np.mean(Ssc[-50:])
        Psc_av = np.mean(Psc[-50:])
        ratio = (Ssc_av/( Ssc_av+ Psc_av))
        Z[i,j] = ratio

'\u03C6'

#############################################
# Graphing non-leaky contour
#############################################

fig3,(ax1) = plt.subplots(sharex = True, sharey = True, dpi = 300)
fig3.suptitle('Contour from non -leaky dynamics')# + str(phi))
# figsize = (8,5)
grid = ax1.pcolormesh( Shs, SNs, np.where(Z == -1, np.nan, Z), vmin=0, vmax=np.max(Z), cmap = 'summer', shading= 'auto'  )  #'gouraud'
ax1.set(xlabel='Supply H$_2$O$_2$')
ax1.set(ylabel='Supply N')
fig3.colorbar(grid, cmap= 'summer',label = 'Syn / Syn + Pro')

fig3.savefig('../figures/nonleaky_dynamic_contour',dpi=300)


'''

###########################################################################
# dynamics contour with actual Syn WH8703 phi 
###########################################################################

#dynamics
for (i,SN) in zip(range(SNs.shape[0]),SNs):
    for (j,Sh) in zip(range(Shs.shape[0]),Shs):
        params = [k1p,k1s,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh]
        leaky  = odeint(leak, inits, mtimes, args = (params,))
        #get equilibria with function 
        Psc = leaky[:,0]
        Ssc = leaky[:,1]
        Nsc = leaky[:,2]
        Hsc = leaky[:,3]
        '''
        Nstar, Pstar, Sstar, Hstar = Coexist(params)
        if Pstar <=0:
            Nstar, Pstar, Sstar, Hstar = Swins(params)
        elif Sstar <=0:
            Nstar, Pstar, Sstar, Hstar = Pwins(params)
        '''
        Ssc_av = np.mean(Ssc[-50:])
        Psc_av = np.mean(Psc[-50:])
        ratio = (Ssc_av/( Ssc_av+ Psc_av))
        Z[i,j] = ratio



#############################################
# Graphing Cotnour plots from dynamics
#############################################

fig3,(ax1) = plt.subplots(sharex = True, sharey = True, dpi = 300)
fig3.suptitle('Contour from leaky dynamics')# + str(phi))
# figsize = (8,5)
grid = ax1.pcolormesh( Shs, SNs, np.where(Z == -1, np.nan, Z), vmin=0, vmax=np.max(Z), cmap = 'summer', shading= 'auto'  )  #'gouraud'
ax1.set(xlabel='Supply H$_2$O$_2$')
ax1.set(ylabel='Supply N')
fig3.colorbar(grid, cmap= 'summer',label = 'Syn / Syn + Pro')

fig3.savefig('../figures/leaky_dynamic_contour',dpi=300)



############################################################################################

#Equilibria contours

##############################################################################################

for (i,SN) in zip(range(SNs.shape[0]),SNs):
    for (j,Sh) in zip(range(Shs.shape[0]),Shs):
        params = [k1p,k1s,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh]
        #get equilibria with function 
        Nstar, Pstar, Sstar, Hstar = Coexist(params)
        if Pstar <=0:
            Nstar, Pstar, Sstar, Hstar = Swins(params)
        elif Sstar <=0:
            Nstar, Pstar, Sstar, Hstar = Pwins(params)
        ratio = (Sstar/( Sstar+ Pstar))
        Z[i,j] = ratio

# equilibria values 

#######################################
# Graphing Cotour plots from equilibria
######################################

fig4,(ax1) = plt.subplots(sharex = True, sharey = True, dpi = 300)
fig4.suptitle('Contour from Equilibria')# + str(phi))
# figsize = (8,5)
grid = ax1.pcolormesh( Shs, SNs, np.where(Z == -1, np.nan, Z), vmin=0, vmax=np.max(Z), cmap = 'summer', shading= 'auto'  )  #'gouraud'
ax1.set(xlabel='Supply H$_2$O$_2$')
ax1.set(ylabel='Supply N')
fig4.colorbar(grid, cmap= 'summer',label = 'Syn / Syn + Pro')


fig4.savefig('../figures/leaky_equilibria_contour',dpi=300)




print('') 

print('*** Done ***')