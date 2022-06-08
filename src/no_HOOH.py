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



'''


mu = (k2*N)/((k2/k1)+N)
dPdt = mu*P
dSdt = mu*S
dNdt = Supply - mu*P - muS

'''