#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: tteodoro
"""
############## Modules ###############
import numpy as np


############# Functions ##############
def sim_IR (calc, obsv):
    """Measures the likeness between two (e.g.,calculated and observed) 
    IR spectra.
    Eqn 1: J. Shen et al. / Spectrochimica Acta Part A 76 (2010) 418-422"""
    Ico = overlap_integrals (calc, obsv)
    Icc = overlap_integrals (calc, calc)
    Ioo = overlap_integrals (obsv, obsv)
    return Ico / (Icc + Ioo - Ico)


def overlap_integrals (i, j):
    """Calculates the (self and) overlap integrals of intensities i and j.
    Eqn 2: J. Shen et al. / Spectrochimica Acta Part A 76 (2010) 418-422"""
    if len(i) != len(j): 
        raise ValueError("Spectra i and j have different lengths.")
    else: 
        overlap = 0.
        for v in range(len(i)):
            overlap += i[v] * j[v]
    return overlap 


def broaden_spectrum(peaks, points, hw):
    """(Lorentzian) Broadening of a spectrum (peaks) for a given frequency
    range (points) given a certain half-width. Units are not considered here"""
    freqs, intns = zip(*peaks)    # Slice to take out the peak positions
    broad_intns  = np.zeros_like(points)
    for freq, ints in zip(freqs,intns):
        broad_intns += Lorentzian(points, freq, hw) * ints
    return broad_intns


def Lorentzian(x,x0,hw):
    """Lorentzian line shape function.
    math:: L(x) = \frac{1}{\pi} \alpha \frac{1}{(x-x_0)^2 + \alpha^2}"""
    return (1/np.pi) * hw / ((x - x0)**2 + hw**2)


def normalize(v):
    """Normalizes a given vector v"""
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm
    

############# Input #################
# If defalut (0 & 0), the range is determined by the frequency range in the
# observed spectrum
vmin = 0  #Beginning of frequency range
vmax = 0  #End of frequency range
# Path for the calculated and observed files
calc_path = "ir_BLYP_DZ.xy"
obsv_path = "ir_BP86_DZ.xy"  
hw = 4.     # Half-width
resol = 1.  # Resolution of broadened spectra in n_points/cm-1
    

############ Preparing data ##############
# Loads data
calc = np.loadtxt(calc_path)
obsv = np.loadtxt(obsv_path)
# Unzip data
calc_freq, calc_ints = zip(*calc)
obsv_freq, obsv_ints = zip(*obsv)

# Generates a new grid of frequencies based on observed data or input 
if not vmin and not vmax:
    broad_freqs = np.arange(obsv_freq[0],obsv_freq[-1],1/resol)
else:
    broad_freqs = np.arange(vmin,vmax,1/resol)
#Re-broadens both spectra for consistency
broad_cints = broaden_spectrum(calc, broad_freqs, hw)
broad_oints = broaden_spectrum(obsv, broad_freqs, hw)
#Normalize intensities
broad_cints = normalize(broad_cints)
broad_oints = normalize(broad_oints)


############ Running ###############
sim = sim_IR(broad_cints, broad_oints)


############ Output ################
print('S = {:>4.2f}'.format(sim))



 

