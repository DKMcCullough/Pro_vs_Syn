#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 11:50:19 2022
name: P_S_mod_comp 

location: /Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Competitions/Pro_vs_Syn/src

Goals: 
    make Pro and Syn modular competition file for leaky and non-leaky HOOH detox via Syn in N-limited growth
    create 2 parameter sets for leaky or non-leaky Pro vs Syn Competition 
    differnt Nstar and Hstars depending on who wins and what equilibrium is being/
        reached by what means? 



@author: dkm

"""




import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import *
from scipy.integrate import odeint



