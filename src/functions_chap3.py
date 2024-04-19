#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 13:57:23 2024

@author: dkm
"""

#Base parameter set  (leaky or no H are special cases) 

k1p =  0.02     #Pro alpha
k1s =  0.01    #Syn alpha 
k2p =  0.5    #Vmax     P
k2s =  0.5  #Vmax     S  #0.300748361
ksp = k2p/k1p
kss = k2s/k1s
dp = 0.2   #pro delta
ds =  0.2   #syn delta
kdam = 0.005555   #hooh mediated damage rate of Pro  
deltaH = 0.002       #decay rate of HOOH via Syn 
phi = 1.7e-6    #0007  #detoxification-based decay of HOOH via Syn in this case
rho =  0.002
Sh = 400
SN = 8000
Qnp = 1#(9.4e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertilison 
Qns = 1#(20.0e-15*(1/(14.0))*1e+9) 

B = 1e5   #this is ALLL bacteria, including EZ55 #EZ55 at 5e+3 can save bc of higher phi (0.0096 - DT's functions)
phib = 0.00004   #from Emily Landlot data and DT calcs
deltah = (deltaH + B*phib)

params = [k1p,k1s,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh]

##############
#make param for heterotroph detox value and heterotroph to be unchanging with current N and also to be callable into leaky function
#############

#functions for model 

def leak(y,t,params):
    k1p,k1s,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh = params[0], params[1], params[2], params[3],params[4],  params[5], params[6], params[7],params[8],params[9], params[10],params[11]
    ksp = k2p/k1p
    kss = k2s/k1s
    P,S,N,H = y[0],y[1],y[2],y[3]
    dPdt = (k2p * N /( (ksp) + N) )*P - (dp *P) - kdam*H*P
    dSdt =(k2s * N /( (kss) + N))*S - (ds *S)      
    dNdt =  SN - ((k2p * N /( (ksp) + N) )*P* Qnp) - ((k2s * N /( (kss) + N))*S* Qns) - rho*N    
    dHdt = Sh - deltah*H  - phi*S*H  
    return [dPdt,dSdt,dNdt,dHdt]

#0.0004   from EZ55 calculcations from E.L. 
#1e5 assumed bacterial steady state concentration

#################################################
#equilibrium solutions for leaky 
################################################
# no H needs own * EQS 


def death_all(params):
    Nstar = SN/(rho)
    Pstar = 0
    Sstar = 0
    Hstar = Sh/(deltah )
    return  Nstar, Pstar, Sstar, Hstar 



def Pwins (params): 
    k1p,k1s,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh = params[0], params[1], params[2], params[3],params[4],  params[5], params[6], params[7],params[8],params[9], params[10],params[11]
    #print(params)
    ksp = k2p/k1p
    kss = k2s/k1s
    Hstar = Sh/deltah
    Nstar =  ((ds +kdam*Hstar)*ksp)/(k2p - ((ds+kdam*Hstar)))
    Pstar = ((SN-rho*Nstar)*(Nstar+(ksp))) / (k2p*Nstar)
    Sstar = 0
   #do we need toassume H must be 0 for P to win?????
    return  Nstar, Pstar, Sstar, Hstar 




def Swins (params): 
    k1p,k1s,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh = params[0], params[1], params[2], params[3],params[4],  params[5], params[6], params[7],params[8],params[9], params[10],params[11]
    #print(params)
    Nstar = - (ds*k2s)/(k1s*(ds-k2s))
    Pstar = 0
    Hstar = Sh/(deltah)
    Sstar = ((SN-rho*Nstar)*(Nstar+(k2s/k1s))) / (k2s*Nstar)
    return  Nstar, Pstar, Sstar, Hstar 



def Coexist (params): 
    k1p,k1s,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh = params[0], params[1], params[2], params[3],params[4],  params[5], params[6], params[7],params[8],params[9], params[10],params[11]
   # print(params)
    ksp = k2p/k1p
    kss = k2s/k1s
    a1 = k1p*ds*(dp-k2p)*k2p 
    a2  = k1s*dp*k2p*(-1*(ds)+k2s)
    a3 = kdam*(k1s*k2p*(ds-k2s)-k1p*ds*k2s)
    Hstar = ((a1+a2)/a3)
    Nstar =  - (ds*k2s)/(k1s*(ds-k2s))
    Sstar = (Sh - deltah*Hstar) /( phi*Hstar)
    Pstar = ((SN - ((k2s*Nstar)/(Nstar + (kss)))*Sstar  - rho*Nstar)  /((k2p*Nstar)/ ((Nstar +( ksp)))))

    return  Nstar, Pstar, Sstar, Hstar 

