#!/usr/bin/env python
from math import *
# from IonizationEnergies import *
# from rr import *
import copy
import sys
import time

"""
This module contains many of the functions for calculating energy dynamics of the charge state distribution. 
"""

__EMASS__     = 5.11e5  # "Electron mass in eV"
__C__         = 3.0e10  # "Speed of light in cm/s"
__ECHG__      = 1.6e-19  # "Electron charge"
__VBOHR__     = 2.2e8  # "Bohr velocity in cm/s"
__AMU__       = 9.311e8  # "1 AMU in eV"
__TORR__      = 3.537e16  # "1 torr in cm-3"
__ALPHA__     = 7.2974e-3  # "fine-structure constant"
__CHEXCONST__ = 2.25e-16  # "Constant for use in charge-exchange "
__LCONV__     = 3.861e-11  # "length conversion factor to CGS"
__KB__        = 1.381e-23  # "Boltzmann constant"
__EPSILON0__  = 8.8541878128e-14 # "Vaccuum permittivity in F/cm"
__EVCONV__    = 1.602176e-19 # Conversion for eV to J
__PMASS__     = 1.672e-27  # "Proton mass in kg"



# Charge exchange constants
__SALZBORNAK__ = 1.43E-12 # "Constants from Mueller & Salzborn, Sept 1977"
__SALZBORNALPHAK__ = 1.17
__SALZBORNBETAK__ = 2.76

__MAXCHARGE__ = 105
__MAXSPECIES__ = 1000


def coulombLogarithm_ei(Eb, qi):
	""" The Coulomb logarithm for electron-ion collisions. The natural log of the ratio of the maximum impact
	parameter to the minimum impact parameter. In a plasma the maximum parameter would usually be the Debye
	length, but in EBIT's it is the inner radius of the central drift tube.

	This is a dynamic variable which needs to be calculated within the calculateK() function!
	"""
	if qi==0:
		return 0
	r_dt = 7e-1 # 0.7 cm, 0.25 cm for neighbors
	b90 = __ECHG__**2*qi/(8*pi*__EPSILON0__*Eb*__EVCONV__) # 90 deg impact parameter. Convert eV to J.

	# should be on the order of 10
	# With argon and 7 keV beam, we get about 24--22
	return log(r_dt/b90)

def rateSpitzerHeating(mySpecies, myEbitParams):
	""" Landau-Spitzer heating of ion bunches by the electron beam. This calculates all of the static values, while the
	dynamic values are left to be calculated in the calculateK() function. Static values are called at the beginning
	of the simulation in the adaptiveRkStepper() loop. Look at the function calcRateMatrices().

	d_t(Ni*kTi) = f(e,i)*e^2*qi^2*n_e*ln(Lambda(e,i))/(6*pi*epsilon_0^2*m_i*v_e)

	n_e is electron number density in m^-3
	f(e,i) is unitless electron-ion overlap function
	v_e is speed of electrons

	Returns eV/s
	"""
	# number density of electrons, cm^-3
	numberDensity = myEbitParams.currentDensity / (__ECHG__*myEbitParams.electronVelocity)
	# mass of ion in kg
	mi = mySpecies.A * __PMASS__

	rates = [0]*(mySpecies.Z + 1)
	for i in range(0, mySpecies.Z+1):
		# the factor 1e4 ahead is conversion for cm^2 to m^2 to get the correct units.
		rates[i] = (1e4/__EVCONV__)*(__ECHG__*i*__ECHG__)**2 *numberDensity*coulombLogarithm_ei(myEbitParams.beamEnergy, i)/(6*pi*__EPSILON0__**2 * mi * myEbitParams.electronVelocity)
		# print("Coulomb Log: %s"%str(coulombLogarithm_ei(myEbitParams.beamEnergy, i)))
	return rates


# def calcCoulombLogarithm_ei(mySpecies, myEbitParams):
# 	""" The Coulomb logarithm for electron-ion collisions. 

# 	Semi-empirical formulae obtained from the 2019 Plasma Formulary from the US Naval Research Laboratory:

# 	https://www.nrl.navy.mil/ppd/content/nrl-plasma-formulary 

# 	------------

# 	This is required for Spitzer heating, wherein the ion beam is heated by Coulomb collisions with the electron beam.

# 	So, this simply returns a number in the distribution will shift when the population of ions changes.
# 	"""
# 	if Te > Ti*Me/Mi:
# 		if Te < 10*mySpecies.Z**2:
# 			return 23 - np.log((myEbitParams.population**1/2) * (mySpecies.Z) * (myEbitParams.electronTemperature**-3/2) )
# 		elif Te > 10*mySpecies.Z**2:
# 			return 24 - np.log((myEbitParams.population**1/2) * (myEbitParams.electronTemperature**-1))
# 		else:
# 			# If we come here there is a problem
# 			sys.exit(1)
# 	elif Te < Ti*Me/Mi:
# 		return 16 - np.log((mySpecies.population**1/2) * (mySpecies.popTemperature**-3/2) * (mySpecies.Z**2) * mu)
# 	else:
# 		# If we come here there is also a problem
# 		sys.exit(1)

def coulombLogarithm_ii(mySpecies, species):
	""" The Coulomb logarithm for ion-ion collisions. 

	Semi-empirical formulae obtained from the 2019 Plasma Formulary from the US Naval Research Laboratory:

	https://www.nrl.navy.mil/ppd/content/nrl-plasma-formulary 

	------------

	This is required for the Fokker-Planck equation, which determines if an energetic ion will escape the trap potential.

	"""
	return









