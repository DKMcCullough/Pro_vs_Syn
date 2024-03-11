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



#####################
# Set UP
######################


step = 0.001
ndays = 400 
mtimes = np.linspace(0,ndays,int(ndays/step))
Shs = np.linspace(0, 30000, num = 10)
SNs = np.linspace(0, 30000, num = 10)
Z = np.zeros((int(SNs.shape[0]),int(Shs.shape[0])),float)

#initial values 
P0 = 1e4
S0 = 1e4
N0 = 1e5        #nM 
H0 = 1    #nm
inits = (P0,S0,N0,H0)



#parameters
Qnp = 1#(9.4e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertillison?  #giving N in micro units? 
Qns = 1#(20.0e-15*(1/(14.0))*1e+9)

k1p =  0.00002     #Pro alpha
k1s =  0.00001      #Syn alpha 
k2 =  0.88    #Vmax    shared for P and S here
ksp = k2/k1p
kss = k2/k1s
dp = 0.2   #pro delta
ds =  0.2   #syn delta
kdam = 0.005   #hooh mediated damage rate of Pro 
deltah = 0.2       #decay rate of HOOH via Syn 
rho =  0.002                  #innate N loss from system 

#empty arrays 
P = np.array([])
S  = np.array([])
N = np.array([])
H = np.array([])
y = [P,S,N,H]


#function set up for ode int
#taken QNs and QNp out of S and P ODEs together

def nleak(y,t,params):
    ksp,kss,k2,dp,ds,kdam,deltah,rho,SN,Sh= params[0], params[1], params[2],params[3],params[4], params[5],params[6],params[7],params[8],params[9]
    P,S,N,H = y[0],y[1],y[2],y[3]
    dPdt = (k2 * N /( (ksp) + N) )*P - (dp *P) - kdam*H*P
    dSdt = (k2 * N /( (kss) + N) )*S - (ds *S) 
    dNdt = SN - ((k2 * N /( (ksp) + N) )* P * Qnp) - ((k2 * N /( (kss) + N) ) * S * Qns)-rho*N
    dHdt = Sh - deltah*H
    return [dPdt,dSdt,dNdt,dHdt]
'''
#solve ODEs via odeint
nonleaky  = odeint(nleak, inits, mtimes, args = (params,))
'''
##############################

#Contour 

##############################


for (i,SN) in zip(range(SNs.shape[0]),SNs):
    for (j,Sh) in zip(range(Shs.shape[0]),Shs):
        params = [ksp,kss,k2,dp,ds,kdam,deltah,rho, SN, Sh]
        nonleaky  = odeint(nleak, inits, mtimes, args = (params,))
        Psc = nonleaky[:,0]
        Ssc = nonleaky[:,1]
        Nsc = nonleaky[:,2]
        Hcs = nonleaky[:,3]
        if (i == 7) and (j == 7):
            Sh = Shs[j]
            SN = SNs[i]
            Ps = nonleaky[:,0]
            Ss = nonleaky[:,1]
            Ns = nonleaky[:,2]
            Hs = nonleaky[:,3]
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


##########################################################

# Calculated analytical solutions at equilibrium

#Analytical solutions at equilibrium (dx/dt equations above set to 0 and solved)
#Equilibrium  (3 cases: both P and S die, P wins while S->0, S wins while P->0))

##########################################################


##### S wins #########

Nstars = (ds*kss)/((k2*Qns)-ds)

#Hstar = Sh/deltah
Sstar = (SN - rho*Nstars)*(((Nstars + kss)/(k2*Nstars*Qns)))


##### P wins #########
Nstarp = ((ksp*dp )+(ksp*kdam))/((k2*Qnp) - dp - kdam)

Hstar = Sh/deltah

Pstar = (SN - rho*Nstarp)*((Nstarp + ksp)/((k2*Nstarp)*Qnp))

Nstarph = ((ksp*dp )+(ksp*kdam*Hstar))/((k2*Qnp) - dp - (kdam*Hstar))
vHline = ((deltah)/(Pstar*kdam)*((Nstarp+ksp)/(k2*Nstarp*Pstar*Qnp)+(dp*Pstar)))


#graphoinig zoom in 
fig4,  (ax1,ax2) = plt.subplots(2, 1, sharex=True,figsize=(9,5),dpi = 300)
fig4.suptitle('Non-leaky dynamics: no HOOH ')

ax1.plot(mtimes, np.clip(Ps,1,np.max(Ps)) , linewidth = 3, color = 'g', label = 'Pro')
ax1.plot(mtimes, np.clip(Ss,1,np.max(Ss)), linewidth = 3, color = 'orange', label = 'Syn')

ax1.set(ylabel='Cells per ml')

ax2.set_ylabel('Nutrient (per ml)')
ax2.plot(mtimes, np.clip(Ns,10,np.max(Ns)),linewidth = 3, color = 'purple', label = "Nutrient")


#ax3.plot(mtimes, np.clip(Hs,10,np.max(Hs)),linewidth = 3, color = 'red', label = "HOOH")

#ax3.set(xlabel='Time (days)', ylabel='HOOH per ml')

ax1.legend(loc = 'best')

ax1.semilogy()
ax2.semilogy()
#ax3.semilogy()




######################################

#        graphing equilibrium values 

######################################

ax1.axhline(Sstar,color = 'brown', linestyle = "-.",label = 'S*')
ax1.axhline(Pstar,color = 'g', linestyle = ":",label = 'P*')
ax2.axhline(Nstars,color = 'purple', linestyle = "-.",label = 'N*s')
ax2.axhline(Nstarp,color = 'magenta', linestyle = "-.",label = 'N*p')
#ax3.axhline(Hstar,color = 'red', linestyle = "-.",label = 'Hstar')


ax1.legend(loc = 'lower center')
ax2.legend(loc = 'lower center')
#ax3.legend(loc = 'best')






#####################################

#  Graphing

#####################################

fig1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True,figsize=(9,5),dpi = 300)
fig1.suptitle('Non_leaky HOOH')

ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'Pro')#' k1 =' + str(k1p))
ax1.plot(mtimes, Ss, linewidth = 3, color = 'orange', label = 'Syn')#' k1 =' + str(k1s))
ax1.set(ylabel='cells per ml')
#np.clip(Ss,10,10e10)
#np.clip(Ps,10,10e10)
ax2.plot(mtimes, Ns,linewidth = 3, color = 'purple', label = "Nutrient")
ax2.set(ylabel='Nutrient per ml')


ax3.plot(mtimes, Hs,linewidth = 3, color = 'red', label = "HOOH")

ax3.set(xlabel='Time (days)', ylabel='HOOH per ml')


ax1.semilogy()
ax2.semilogy()

######################################

#        graphing equilibrium values 

######################################

ax1.axhline(Sstar,color = 'orange', linestyle = "-.",label = 'S*')
ax1.axhline(Pstar,color = 'g', linestyle = "-.",label = 'P*')
ax2.axhline(Nstars,color = 'purple', linestyle = "-.",label = 'N*s')
ax2.axhline(Nstarp,color = 'magenta', linestyle = "-.",label = 'N*p')
ax3.axhline(Hstar,color = 'red', linestyle = "-.",label = 'H*')


ax1.legend(loc = 'best')
ax2.legend(loc = 'best')
ax3.legend(loc = 'best')

fig1.tight_layout()
plt.show()

fig1.savefig('../figures/no_leak_contour_f1_auto',dpi=300)

###############
#graphing contour
#####################

fig3,(ax1) = plt.subplots(sharex = True, sharey = True, figsize = (8,5),dpi = 300)
fig3.suptitle('Non_leaky Contour')
grid = ax1.pcolormesh( Shs,SNs, np.where(Z == -1, np.nan, Z), vmin=0, vmax=np.max(Z), cmap = 'summer', shading='auto')  #'summer_r is reversed color map shading
#np.where(Z == 17, np.nan, Z)

#ax1.axhline((rho*Nstars),color = 'purple', linestyle = "-.",label = 'SN cutoff for S growth?')
#ax1.axhline((rho*Nstarp),color = 'magenta', linestyle = "-.",label = 'SN cutoff for P growth?')
#ax1.axhline((rho*Nstarph),color = 'orange', linestyle = "-.",label = 'SN cutoff for P+H growth?')
#ax1.axvline((vHline),color = 'c', linestyle = "-.",label = 'H cutoff?')

#ax1.axhline(((rho*Nstar)+(((k2*Nstar)/(Nstar-kss))*Sstar*Qns)),color = 'magenta', linestyle = "-.",label = 'Sn cut off on S?')


ax1.set(xlabel='Supply hooh')
ax1.set(ylabel='Supply nutrient')

fig3.colorbar(grid, cmap= 'summer',label = 'S / (S+P)')
plt.legend()

fig3.savefig('../figures/no_leak_contour_f3_auto',dpi=300)






print('')
print('*** Sstar and Hstar on countour ***')
print('*** Done ***')