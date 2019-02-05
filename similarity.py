#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: tteodoro
"""
##################################   Modules   ################################
import numpy as np
import matplotlib.pyplot as plt


#################################    Classes   ################################
class Spectra:
    """Class to define a vibrational spectrum."""
    
    def __init__(self, file):
        """Initializes the class"""
        self.xy = np.loadtxt(file,comments=['C','#'])
        self.freq, self.ints = zip(*self.xy)
            
    def lorentzian (self, x, x0, w):
        """Lorentzian line shape function.
        math:: L(x) = \frac{1}{\pi} \alpha \frac{1}{(x-x_0)^2 + \alpha^2}"""
        alpha = w / 2
        return (1/np.pi) * alpha / ((x - x0)**2 + alpha**2)
    
    def gaussian(self, x, x0, w):
        """Gaussian line shape"""
        sigma = w / (2 * np.sqrt(2 * np.log(2)))
        return 1 / (np.sqrt(2*np.pi)*sigma) * np.exp( -((x-x0)/sigma)**2 / 2 )
    
    def line_shape (self, name):
        """
        Chooses between Lorentzian or Gaussian
        """
        name = name.lower()
        funcs = {
                'lorentzian' : self.lorentzian,
                'gaussian'   : self.gaussian
                }
        if name not in funcs:
            raise ValueError("Non-existent type of function.\
                             Must choose between 'Lorentzian' and 'Gaussian")
        return funcs[name]
    
    def broaden (self, freq_range, width=8.0, scaling=1.00, offset=1.00, \
                 shape = 'lorentzian', fix_width = True):
        """(Lorentzian) Broadening of a spectrum (peaks) for a given frequency
        range (points) given a certain half-width. Units are 
        not considered here"""
        new_ints = np.zeros_like(freq_range)
        shape_func = self.line_shape(shape)
        for freq, ints in zip(self.freq, self.ints):
            new_ints += shape_func(freq_range, (freq * scaling  +offset), \
                                   width) * ints
        return new_ints

################################## Functions ##################################
def simil_shen (calc, obsv):
    """Measures the likeness between two (e.g.,calculated and observed) 
    IR spectra.
    Eqn 1: J. Shen et al. / Spectrochimica Acta Part A 76 (2010) 418-422"""
    # Vectors are first normalized
    calc = normalize(calc)
    obsv = normalize(obsv)
    # Integrals are then calculated
    Ico = overlap (calc, obsv)
    Icc = overlap (calc, calc)
    Ioo = overlap (obsv, obsv)
    return Ico / (Icc + Ioo - abs(Ico))

def normalize (v):
    """Normalizes a given vector v"""
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm

def overlap (i, j):
    """Calculates the (self and) overlap integrals of intensities i and j.
    Eqn 2: J. Shen et al. / Spectrochimica Acta Part A 76 (2010) 418-422"""
    if len(i) != len(j): 
        raise ValueError("Spectra i and j have different lengths.")
    else: 
        overlap = 0.
        for v in range(len(i)):
            overlap += i[v] * j[v]
    return overlap 
    
    
################################   Input   ####################################
# Path for the experimental spectrum
obsv_path = "obsv.xy"
# Path to the calculated peaks (un-broadened)
calc_path = "calc.xy"

# Spectral settings  
width  = 25.   # Width
scale  =  1.01 # Shift in x
offset = 20.   # Offset in x
shape  = 'Lorentzian' # Shape function

vmin =  800.0   #Beginning of frequency range to zoom in/out
vmax = 1800.0   #End of frequency range to zoom in/out
    
############################   Preparing data   ###############################
# Loads data
calc = Spectra(calc_path)
obsv = Spectra(obsv_path)
    
#Broadens calculated data according to observed points (frequencies)
broad_calc_ints = calc.broaden(obsv.freq, width, scale, offset, shape)

##########################    Running     #####################################
glob_sim  = simil_shen(broad_calc_ints, obsv.ints)

range_freq = [i for i in range(len(obsv.freq)) if \
             (obsv.freq[i] >= vmin and obsv.freq[i] <= vmax)] 

local_sim = simil_shen(np.array(broad_calc_ints)[range_freq], 
                     np.array(obsv.ints)[range_freq])  


#############################    Output    ####################################
fig = plt.figure()

ax  = fig.add_subplot(111)
ax1 = fig.add_subplot(211) 
ax2 = fig.add_subplot(212) 

ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax.set_xlabel('Frequency (cm-1)')
ax.set_ylabel('Normalized intensities')
ax.set_yticklabels([])

ax1.plot(obsv.freq,normalize(obsv.ints),label=('Normal. Exp. Data'))
ax1.plot(obsv.freq,normalize(broad_calc_ints),label=('Normal. Calc. Data'))
ax1.text(obsv.freq[-1], max(normalize(obsv.ints))/2, \
         ' S = {:>4.2f}'.format(glob_sim))
ax1.set_yticklabels([])
ax1.set_xlim(obsv.freq[0],obsv.freq[-1])
ax1.tick_params(left=False)
ax1.legend()

ax2.plot(obsv.freq,normalize(obsv.ints),label=('Normal. Exp. Data'))
ax2.plot(obsv.freq,normalize(broad_calc_ints),label=('Normal. Calc. Data'))
ax2.text(vmax, max(normalize(broad_calc_ints))/2, \
         ' S = {:>4.2f}'.format(local_sim))
ax2.set_yticklabels([])
ax2.set_xlim(vmin, vmax)
ax2.tick_params(left=False)

plt.show()

 

