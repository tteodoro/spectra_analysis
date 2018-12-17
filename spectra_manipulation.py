#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 10:26:32 2018

@author: tteodoro
"""
############## Modules ###############
import numpy as np


############# Functions ##############
def overlap_integrals (i, j):
    """Calculates the (self and) overlap integrals of spectra i and j."""
    
    if len(i) != len(j): 
        raise ValueError("Spectra i and j have different lengths.")
    else: 
        overlap = 0.
        for v in range(len(i)):
            overlap += i[v] * j[v]
    return overlap 


def SimIR (calc, obsv):
    """Given a spectral range (vmin to vmax) calculates the spectrum similarity
    for measuring the likeness between calculated and observed IR spectra."""

    Ico = overlap_integrals (calc, obsv)
    Icc = overlap_integrals (calc, calc)
    Ioo = overlap_integrals (obsv, obsv)
    
    return Ico / (Icc + Ioo - Ico)
    

############# Input #################
vmin = 1000.
vmax = 2500.    
    

############ External ##############
calc = np.loadtxt("calc.txt")
obsv = np.loadtxt("obsv.txt")

calc_freq, calc_ints = zip(*calc)
obsv_freq, obsv_ints = zip(*obsv)

calc_ints = [x / 10 for x in calc_ints]

selec_ints = [i for i in range(len(calc_ints)) \
              if (calc_freq[i]>=vmin and calc_freq[i]<=vmax)]

calc_selec_ints = calc_ints[selec_ints[0]:selec_ints[-1]]
obsv_selec_ints = obsv_ints[selec_ints[0]:selec_ints[-1]]


############ Running ###############
simil = SimIR(calc_selec_ints, obsv_selec_ints)



 

