#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 13:57:23 2024

@author: dkm
"""

#Base parameter set  (leaky or no H are special cases) 


Qnp = 1#(9.4e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertilison 
Qns = 1#(20.0e-15*(1/(14.0))*1e+9) 

k1p =  0.02     #Pro alpha
k1s =  0.01    #Syn alpha 
k2p =  0.5    #Vmax     P
k2s =  0.5  #Vmax     S  #0.300748361
ksp = k2p/k1p
kss = k2s/k1s
dp = 0.2   #pro delta
ds =  0.2   #syn delta
kdam = 0.005   #hooh mediated damage rate of Pro  
deltah = 0.002       #decay rate of HOOH via Syn 
phi = 1.7e-6    #0007  #detoxification-based decay of HOOH via Syn in this case
rho =  0.002
Sh = 0
SN = 10000

B = 1e5
phib = 0.0004

params = [ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh,B,phib]

##############
#make param for heterotroph detox value and heterotroph to be unchanging with current N and also to be callable into leaky function
#############

#functions for model 

def leak(y,t,params):
    ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh,B,phib = params[0], params[1], params[2], params[3],params[4],  params[5], params[6], params[7],params[8],params[9], params[10],params[11], params[12], params[13]
    P,S,N,H = y[0],y[1],y[2],y[3]
    dPdt = (k2p * N /( (ksp) + N) )*P - (dp *P) - kdam*H*P
    dSdt =(k2s * N /( (kss) + N))*S - (ds *S) #- kdams*H*S      
    dNdt =  SN - ((k2p * N /( (ksp) + N) )*P* Qnp) - ((k2s * N /( (kss) + N))*S* Qns) - rho*N    
    dHdt = Sh - deltah*H  - phi*S*H  - (phib)*(B)*(H) #phi being S cell-specific detox rate
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
    Hstar = Sh/(deltah + ((B)*(phib)))
    return  Nstar, Pstar, Sstar, Hstar 



def Pwins (params): 
    Nstar = - (((k2p*(Sh*kdam + dp*(deltah+((B)*(phib))))/(k1p*(Sh*kdam + dp*(deltah*((B)*(phib)))-k2p*(deltah*(B)*(phib)))))))
    Pstar = (((deltah+(((B)*(phib)))*(rho*k2p*(Sh*kdam))+SN*k1p*(Sh*kdam +dp*(Sh + (B)*(phib) ) - k2p*(Sh + (B)*(phib))))/(k1p*(Sh*kdam + dp*(Sh + (B)*(phib)))*(Sh*kdam + dp*(Sh + (B)*(phib) ) - k2p*(Sh + (B)*(phib))))))
    Sstar = 0
    Hstar = Sh/(deltah+(B)*(phib))  #do we need toassume H must be 0 for P to win?????
    return  Nstar, Pstar, Sstar, Hstar 




def Swins (params): 
    Nstar = - (ds*k2s)/(k1s*(ds-k2s))
    Pstar = 0
    Hstar = ((Sh*k2s)*((ds*k2s)+((k1s)*k1s))*(-ds + k2s))/((-rho*k1s*ds*k1s*phi)+ ds*((k2s)*k2s)*(deltah+ (B*phib))*((SN*phi)+ (k2s)*(deltah+ (B*phib)) ))
    Sstar = -(k1s*(SN*k1s*(ds-k2s)+rho*ds*k1s))/(k2s*((ds*k2s)+((k1s)*(k1s)))*(-ds + k2s))
    return  Nstar, Pstar, Sstar, Hstar 



def Coexist (params): 
    Nstar = - (ds*k2s)/(k1s*(ds-k2s))
    
    Pstar = ((k1s*k2p*(ds-k2s) - (k1p*ds*k2s)*(k1p*ds*ds*k2s*k2s*k2s)*(Sh*kdam + dp*(deltah + (B*phib))- k1p*(deltah + (B*phib)))) + (k1s*k1s*k1s*k2p*((ds-k2s)*(ds-k2s))*(Sh*kdam*k2s + dp*(SN*phi + k2s*(deltah + B*phib))) + k1s*ds*k2s*k2s*(k2p*(Sh*kdam*k2p + ds*(-Sh*kdam + rho*k1p*phi)))) \
             - dp*(k2p*k2s*((deltah + (B*phib)))) + ds*((rho*k1p*phi + (k2p*(deltah + B*phib)))) - (k1s*k1s*ds*(ds - k2s)*k2s*((-rho*dp*k2p*phi)+(k1p*(Sh*kdam*k2p + dp*(SN*phi + k2s*(deltah*(B*phib))))) -  k2p*(SN*phi + k2s*(deltah + (B*phib))) ) / ((k1p*k1s*k1s*ds*k2p*k2s*(-ds+k2s)*(k1s*dp*k1p*(ds-k2s)+k1p*ds*(-dp+k2p))*k2s)*phi)))
    Hstar = ((k1p*dp*(dp-k2p)*k2s) + (k1s*dp*k1p)*(-dp+k2p) )/((kdam*(k1s*k2p*(dp-k1p) - k1p*ds*k2s)))
    Sstar = -(k1s*k1p*(ds-k2s)*(Sh*kdam + dp*(deltah + (B*phib)))) + (k1p*ds*k1s*(Sh*kdam + dp*(deltah + (B*phib)) - k2p*(deltah + (B*phib))))
    return  Nstar, Pstar, Sstar, Hstar 


def StarContour(Nstar, Pstar, Sstar, Hstar):
    #
    return
#N threshold 

#Nstarph = ((ksp*dp )+(ksp*kdam*Hstar))/((k2p*Qnp) - dp - (kdam*Hstar))

#h thresthold 
#vHline = ((deltah)/(Pstar*kdam)*((Nstarp+ksp)/(k2p*Nstarp*Pstar*Qnp)+(dp*Pstar)))

