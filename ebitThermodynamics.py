#!/usr/bin/env python
from math import *
from IonizationEnergies import *
from rr import *
import copy
import sys
import time
import numpy as np

"""
This module contains many of the functions for calculating the thermodynamics inside of the EBIT.
"""

__EMASS__ = 5.11e5  # "Electron mass in eV"
__C__ = 3.0e10  # "Speed of light in cm/s"
__ECHG__ = 1.6e-19  # "Electron charge"
__VBOHR__ = 2.2e8  # "Bohr velocity in cm/s"
__AMU__ = 9.311e8  # "1 AMU in eV"
__TORR__ = 3.537e16  # "1 torr in cm-3"
__ALPHA__ = 7.2974e-3  # "fine-structure constant"
__CHEXCONST__ = 2.25e-16  # "Constant for use in charge-exchange "
__LCONV__ = 3.861e-11  # "length conversion factor to CGS"
__kB__ = 1.381e-23  # "Boltzmann constant"



# Charge exchange constants
__SALZBORNAK__ = 1.43E-12 # "Constants from Mueller & Salzborn, Sept 1977"
__SALZBORNALPHAK__ = 1.17
__SALZBORNBETAK__ = 2.76

__MAXCHARGE__ = 105
__MAXSPECIES__ = 1000



def calcCoulombLogarithm_ei(mySpecies, myEbitParams):
	""" The Coulomb logarithm for electron-ion collisions. 

	Semi-empirical formulae obtained from the 2019 Plasma Formulary from the US Naval Research Laboratory:

	https://www.nrl.navy.mil/ppd/content/nrl-plasma-formulary 

	------------

	This is required for Spitzer heating, wherein the ion beam is heated by Coulomb collisions with the electron beam.

	So, this simply returns a number in the distribution will shift when the population of ions changes.
	"""
	if Te > Ti*Me/Mi:
		if Te < 10*mySpecies.Z**2:
			return 23 - np.log((myEbitParams.population**1/2) * (mySpecies.Z) * (myEbitParams.electronTemperature**-3/2) )
		elif Te > 10*mySpecies.Z**2:
			return 24 - np.log((myEbitParams.population**1/2) * (myEbitParams.electronTemperature**-1))
		else:
			# If we come here there is a problem
			sys.exit(1)
	elif Te < Ti*Me/Mi:
		return 16 - np.log((mySpecies.population**1/2) * (mySpecies.popTemperature**-3/2) * (mySpecies.Z**2) * mu)
	else:
		# If we come here there is also a problem
		sys.exit(1)

def calcCoulombLogarithm_ii(mySpecies, species):
	""" The Coulomb logarithm for ion-ion collisions. 

	Semi-empirical formulae obtained from the 2019 Plasma Formulary from the US Naval Research Laboratory:

	https://www.nrl.navy.mil/ppd/content/nrl-plasma-formulary 

	------------

	This is required for the Fokker-Planck equation, which determines if an energetic ion will escape the trap potential.

	"""
	return









