"""
This is simply a file to store all of the docstrings for the documentation for CBSim. I place it into a new file so that it doesn't cause clutter for the scripts
"""


def docs_physics():
	"""	\n 	Physics implementation in CBSim
	===============================

General
-------
The current implementation of CBSim accounts for the following ionization and recombination mechanisms:
	- electron impact ionization (EI)
	- radiative recombination (RR)
	- charge exchange (CX)

For a species in a specified charge state i, the rate equation is written as

	dNi/dt = + (EI rate of charge state i-i) - (EI rate of charge state i  )
		 + (RR rate of charge state i+1) - (RR rate of charge state i  )
		 + (CX rate of charge state i+1) - (CX rate of charge state i  )
		 - Resc

The final term in the equation is accounting for the rate at which ions can escape the trap. This escape occurs either radially or axially when the ion obtains enough kinetic energy to overcome the poitentials of the trap.


Electron Impact Ionization
--------------------------

EI rates are calculated using formulae of the form:
	Ri = Je/e * Ni * sigmai * f(e, i)

where sigmai is the cross-section and f(e,i) is an electron-ion overlap factor. Je is the current density of electrons.


Radiative Recombination
-----------------------

The RR rates can be calculated using a formula mirroring the EI rates:
	Ri = Je/e * Ni * sigmai * f(e, i)

Charge Exchange
---------------

**NOTE: we aren't entirely sure where this semi-empirical formula comes from because it was originally written by Renee K. It's possible we might change this calculation to a formula given in the work of Mueller and Salzborn, 1977.**

The CX rate is calculated as:
	Ri = q_i * sigma * rho_H2

where rho_H2 is the density of H2, q_i is the charge state of ion i, and sigma is the CX cross section. The CX cross section is calculated by:
	sigma = Ccx * log(15/Vbohr) * Vi_avg

where Ccx is a charge exchange constant = 2.25E-16, Vi_avg is the average ion velocity and Vbohr is the average ion velocity relative to the Bohr velocity. Vi_avg can be calculated with:
	Vi_avg = c * sqrt(8 * Ti / (pi * Mi))

where Ti is the temperature, and Mi is the mass of the ion in AMU. This is the average velocity of a Maxwellian distribution (holds for low density plasmas).

Ion Escape
----------

The ion escape rate is written as:
	Ri = -Ni * Vi * ( exp(-omegai)/omegai - sqrt(omegai)[erf(omegai) - 1] )

where omegai = qi * e * V / kb / Ti. V is the potential trap depth (axially or radially), and Ti is the temperature of the ions.


Beam Dynamics
-------------

As inputs, CBSim takes the radius of the electron beam and the amount of spatial overlap between the electron and ion clouds. This calculation is performed in a 2D manner.




=====================

Please refer to the 'timestepping' topic of docs for more details.
	"""
	return

def docs_parameters():
	"""		Input Parameters
		================

This is a more detailed description of the physics-related input parameters that the user can input through the command line or through the .cfg file.

beamEnergy
----------

UNITS = eV

The electron beam energy is given in units of eV. This simulation assumes a unform beam energy across the radial profile of the beam. Outside of the parameter 'beamRadius', the beam energy is zero.

A general rule of thumb is that the beam energy should be between 3-4x the ionization energy to optimize for the charge state that you want. This is simply a result of competition with recombination processes in the trap.

breedingTime
------------

UNITS = seconds

The full calculation time for the time stepping solver. 

probeEvery
----------

UNITS = seconds

For saving space, the output plot showing the charge state distribution only has as much resolution as is given by 'probeEvery'. Keep in mind that this value does not affect the time stepping size used for calculation population rates, it is purely for display purposes.

ionEbeamOverlap
---------------

UNITS = unitless

This is currently implemented as a value for the amount of spatial overlap between the electron beam and ion cloud. In a future implementation we will use this to calculate the overlap function, f(e, i).

beamCurrent
-----------

UNITS = amps

The total current of the electron beam. Is used with 'beamRadius' to calculate the current density and hence the ionization and recombination cross-sections.

beamRadius
----------

UNITS = meters

The radius of the electron beam. Current implementation is a hard cutoff for the electron continuum, not a tapered distribution. A top-hat distribution.

pressure
--------

UNITS = torr

The background pressure in the EBIT trap. This is used to determine the charge exchange rates using a background of H2 gas.

ionTemperature
--------------

UNiTS = Kelvin

The __initial__ temperature of the ion cloud. It is used to determine the charge exchange rates by estimating the average ion velocity of a Maxwellian distribution at this temperature.

A general rule of thumb for setting this value is...

	"""

	return

def docs_timestepping():
	"""
		Time Stepping in CBSim
		======================

	Time stepping is using a Runge-Kutte 4 method with an adaptive time step. Because interactions between various populations in the EBIT are accounted for, the overall adapted time step for a single step is limited by any one of the populations. The following illustrates the algorithm:

	blah blah blah

	"""