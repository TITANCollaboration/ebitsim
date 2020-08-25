"""
This is simply a file to store all of the docstrings for the documentation for CBSim. I place it into a new file so that it doesn't cause clutter for the scripts
"""


def docs_physics():
	"""	\n 	Physics implementation in CBSim
	===============================

General
-------

Following the 2005 paper by Fred Currell and Gerd Fussmann, we consider the following simplifications:
 
 1. The electron beam has a radial top hat profile. The radius prescribed in the configuration file is nearly the same as the Herrmann radius. Inside of this radius is an electron beam of uniform density and energy and outside of the radius is zero charge.

 2. For both the electron beam and the ion cloud, we assume that the axial distributiona are uniform along the length of the trap.

 3. (currently being implemented) The radial distribution of ions depends on charge state. They follow a Boltzmann distribution.

 The current implementation of CBSim accounts for the following ionization and recombination mechanisms:
	- electron impact ionization (EI)
	- radiative recombination (RR)
	- charge exchange (CX)

For a species in a specified charge state i, the rate equation is written as

	dNi/dt = + (EI rate of charge state i-i) - (EI rate of charge state i  )
		 + (RR rate of charge state i+1) - (RR rate of charge state i  )
		 + (CX rate of charge state i+1) - (CX rate of charge state i  )
		 - Resc

Where Ni is the total number of ions per length in the trap. It's a good idea to look this up. At TITAN we can inject about 10^6 ions per bunch, but we can load up to a total capacity of about 10^8? I'm not sure, but check the Thesis of Annika Lennarz and the section where she discusses the stacked injection scheme.

The final term in the equation is accounting for the rate at which ions can escape the trap. This escape occurs either radially or axially when the ion obtains enough kinetic energy to overcome the poitentials of the trap.


Electron Impact Ionization
--------------------------

EI rates are calculated using formulae of the form:
	Ri = Je/e * Ni * sigmai * f(e, i)

where sigmai is the cross-section and f(e,i) is an electron-ion overlap factor. Je is the current density of electrons.

sigmai is calculated using the Lotz formula.


Radiative Recombination
-----------------------

The RR rates can be calculated using a formula mirroring the EI rates:
	Ri = Je/e * Ni * sigmai * f(e, i)

The cross section is calculated using a time-reversed photonionization cross section.

Charge Exchange
---------------

The CX rate is calculated as:
    Ri = vi_avg * N0 * Ni * sigmai

where vi_avg is the average speed of the ion based on a Maxwellian speed distribution, N0 is the number density of the background gas, Ni is the number density of ions, and sigma is the cross section.

In this implementation the cross section is calculated by the semi-empirical formula of Mueller and Salzborn, published September 1977.

sigmai = Ak * Epgas^betak * qi^alphak   for k ranging from i=1 to 4.

This is the cross section for charge exchange from charge state i to charge state i-k. Epgas is the ionization potential for the background gas, and qi is the charge state of the ion. The constants Ak, betak, and alphak are given for each integer k. Sor far this has only been implemented for k=1.

Ion Escape
----------

NOT YET IMPLEMENTED

The ion escape rate is written as:
	Ri = -Ni * Vi * ( exp(-omegai)/omegai - sqrt(omegai)[erf(omegai) - 1] )

where omegai = qi * e * V / kb / Ti. V is the potential trap depth (axially or radially), and Ti is the temperature of the ions.


Geometry
--------

We only consider the trapping region in the EBIT




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

UNITS = Kelvin

The __initial__ temperature of the ion cloud. It is used to determine the charge exchange rates by estimating the average ion velocity of a Maxwellian distribution at this temperature.

A general rule of thumb for setting this value is...

populationPercent
-----------------

UNITS = fraction

A fraction of the total population given for the species. If we start with a single species and populationPercent=1.0, then 100% of the population is this species. The program will renormalize the inputs, therefore if we have two species, each with populationPercent=1.0, then they each garner 50% of the total population.

Please note that the current configutation is that the initial population is ALL singly charged ions (SCI). We might make this customizable in the future.

	"""

	return

def docs_timestepping():
	"""
		Time Stepping in CBSim
		======================

	Time stepping is using a Runge-Kutte 4 method with an adaptive time step. Because interactions between various populations in the EBIT are accounted for, the overall adapted time step for a single step is limited by any one of the populations. The following illustrates the algorithm:

	blah blah blah

	"""