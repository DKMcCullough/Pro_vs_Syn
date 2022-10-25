#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name:   no_HOOH_contour.py 

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




#####################
# Set UP
######################


step = 0.01
ndays = 300
mtimes = np.linspace(0,ndays,int(ndays/step))
k2s = np.linspace(0, 1.0, num = 10)
SNs = np.linspace(0, 1000, num = 10)
Z = np.zeros((int(SNs.shape[0]),int(k2s.shape[0])),int)

#initial values 
P0 = 1e4
S0 = 1e4
N0 = 1e5
inits = (P0,S0,N0)



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
rho =  0.002                  #innate N loss from system 

#params = [ksp,kss,k2,dp,ds,rho, SN]

#empty arrays 
P = np.array([])
S  = np.array([])
N = np.array([])
y = [P,S,N]


#function set up for ode int

def noH(y,t,params):
    ksp,kss,k2,dp,ds,rho,SN = params[0], params[1], params[2],params[3],params[4], params[5],params[6]
    P,S,N = y[0],y[1],y[2]
    dPdt = (k2 * N /( (ksp) + N) )*P*Qnp - (dp *P) 
    dSdt = (k2 * N /( (kss) + N) )*S*Qns - (ds *S) 
    dNdt = SN - ((k2 * N /( (ksp) + N) )* P * Qnp) - ((k2 * N /( (kss) + N) ) * S * Qns)-rho*N
    return [dPdt,dSdt,dNdt]




##########################################################

# Calculated analytical solutions at equilibrium

##########################################################


#Pwins

Nstar = (dp*(k2/k1p))/(k2-dp)
Pstar = (N0-rho*Nstar)/(Qnp*(k2*Nstar)/(Nstar + (k2/k1p)))



 #Swins???



#plt.xlim([0, 50])
plt.ylim([10,10e8])

plt.show()

#fig.savefig('../figures/no_HOOH_contour',dpi=300)


for (i,SN) in zip(range(SNs.shape[0]),SNs):
    for (j,k2) in zip(range(k2s.shape[0]),k2s):
        params = [ksp,kss,k2,dp,ds,rho,SN]
        noHooh  = odeint(noH, inits, mtimes, args = (params,))
        Psc =  noHooh[:,0]
        Ssc =  noHooh[:,1]
        Nsc =  noHooh[:,2]
        if (i == 0) and (j == 3):
            k2 = k2s[j]
            Nsupply = SNs[i]
            Ps =  noHooh[:,0]
            Ss =  noHooh[:,1]
            Ns =  noHooh[:,2]
            #need to shoose star values here for dynamic print (N, P or S, H)
        Z[i,j] = Ssc[-1]/(Psc[-1]+Ssc[-1])

#####################################

#  Graphing

#####################################



fig1, (ax1, ax2) = plt.subplots(1, 2)
fig1.suptitle('Growth Competition Projections')
plt.subplots_adjust(wspace = 0.3, top = 0.85)


ax1.plot(mtimes, np.clip(Ps,0.001,10e13) , linewidth = 3, color = 'g', label = 'Pro') #label = 'Pro k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn') #label = 'Syn k1 =' + str(k1s))
ax1.set(xlabel='Time (days)', ylabel='cells per ml')


ax2.plot(mtimes, Ns, linewidth = 3, color = 'purple', label = "Nutrient")
ax2.set(xlabel='Time (days)', ylabel='Nutrient concentration')
#ax2.plot(mtimes, y = (N0/rho), linestyle = ';', color = 'c', label = "N equilibrium solution")

ax1.semilogy()
ax2.semilogy()






######################################

#        graphing equilibrium values 

######################################


ax1.axhline(Pstar, color = 'green', linestyle = "-.",label = 'Pstar')
ax2.axhline(Nstar, color = 'purple', linestyle = "-.",label = 'Nstar')
ax1.legend(loc = 'lower right')
ax2.legend(loc = 'lower right')



plt.show()

fig1.savefig('../figures/no_H_contour_f1_auto',dpi=300)


'''
#### graphing contour? 

fig3,(ax1) = plt.subplots(sharex = True, sharey = True, figsize = (8,5))


grid = ax1.pcolormesh( Shs,SBs, np.where(Z == -1, np.nan, Z), vmin=0, vmax=np.max(Z), cmap = 'summer', shading='auto')  #'summer_r is reversed color map shading
#np.where(Z == 17, np.nan, Z)
#is this just masking the values or actually being helpful to cleaner code? 

#grid2 = ax1.pcolormesh( Shs,SNs, np.where(Z == 4, 1, Z), cmap = 'Greys', shading='auto')

#ax1.axhline((rho*Nstar),color = 'purple', linestyle = "-.",label = 'SN cutoff for ccell growth?')
#ax1.axvline((deltah*Hstar),color = 'magenta', linestyle = "-.",label = 'H cutoff?')

#ax1.axhline(((rho*Nstar)+(((k2*Nstar)/(Nstar-kss))*Sstar*Qns)),color = 'magenta', linestyle = "-.",label = 'Sn cut off on S?')


ax1.set(xlabel='Supply hooh')
ax1.set(ylabel='Supply nutrient')

fig3.colorbar(grid, cmap= 'summer',label = 'S / (S+P)')
plt.legend()

fig3.savefig('../figures/no_leak_contour_f3_auto',dpi=300)

'''
print('*** Done ***')