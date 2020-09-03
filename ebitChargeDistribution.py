#!/usr/bin/env python
from math import *
from IonizationEnergies import *
from ebitThermodynamics import *
from rr import *
import copy
import sys
import time

__EMASS__ = 5.11e5  # "Electron mass in eV"
__C__ = 3.0e10  # "Speed of light in cm/s"
__ECHG__ = 1.6e-19  # "Electron charge"
__VBOHR__ = 2.2e8  # "Bohr velocity in cm/s"
__AMU__ = 9.311e8  # "1 AMU in eV"
__TORR__ = 3.537e16  # "1 torr in cm-3"
__blah__ = 24.1521e23
__ALPHA__ = 7.2974e-3  # "fine-structure constant"
__CHEXCONST__ = 2.25e-16  # "Constant for use in charge-exchange "
__LCONV__ = 3.861e-11  # "length conversion factor to CGS"
__kB__ = 1.381e-23  # "Boltzmann constant"
__Epgas__ = 7.0 # first ionization potential of H2 gas in eV (was 15.42593 eV, Oliver says 7 is correct)


#Random comment


# Charge exchange constants (for q --> q-1)
__SALZBORNAK__ = 1.43E-12 # "Constants from Mueller & Salzborn, Sept 1977"
__SALZBORNALPHAK__ = 1.17
__SALZBORNBETAK__ = -2.76

__MAXCHARGE__ = 105
__MAXSPECIES__ = 1000
INITILIZEDEVERYTHING = 0


class Species:
    """ The Species class is created for each charge state of each species in the EBIT. Each of these populations has their
    own rates for interactions with other populations in the EBIT. With time evolution, the rates cause the populations to
    exchange particles (2+ charge state can "feed" 1+ or 3+ by various interactions) or they can exchange heat. This heat exchange
    will lead to temperature changes which can cause the LOSS of ions. Therefore (once that is implemented) the number of ions
    in the trap will no longer be conserved.
    """
    def __init__(self,
                 Z,
                 A,
                 decaysTo=0.0,
                 betaHalfLife=0.0,
                 initSCIPop=1.0,
                 chargeStates=None,
                 halfLife=0.0,
                 populationNumber=0.0,
                 initSCITemp = None,
                 # The lines above need to stay in the same order
                 NkT=0.0,
                 qkT = None, 
                 population=None,
                 decayConstant=0.0,
                 ionizationRates=None,
                 rrRates=None,
                 chargeExchangeRates=None,
                 spitzerHeatingRates=None,
                 diff=0.0,
                 bestStepSize=0.0,
                 noTooBigSteps=0.0,
                 betaDecayDelta = 0.0,
                 truncationError = 0.0,
                 stepCounter = 0,
                 k1=[],
                 k2=[],
                 k3=[],
                 k4=[],
                 r1=[],
                 r2=[],
                 r3=[],
                 r4=[],
                 tmpPop=[],
                 tmpEnergy=[],
                 y1=[],
                 y12=[],
                 y22=[],
                 f1=[],
                 f12=[],
                 f22=[],
                 results=[]):
        self.Z = Z
        self.A = A
        self.decaysTo = decaysTo
        self.betaHalfLife = betaHalfLife
        self.halfLife = halfLife
        self.populationNumber = populationNumber
        self.betaDecayDelta = betaDecayDelta
        self.truncationError = truncationError
        self.stepCounter = stepCounter
        self.initSCIPop = initSCIPop
        self.initSCITemp = initSCITemp
        self.decayConstant = decayConstant
        self.chargeStates = chargeStates
        self.diff = diff
        self.bestStepSize = bestStepSize
        self.noTooBigSteps = noTooBigSteps
        self.results = results

# I think I need to move chargeStates up to here.. could concievibly do that for a bunch
# of the stuff in chargepopparams due to maybe they should be following the species around


class RkStepParams:
    def __init__(self, minCharge=5e-5, tStep=1e-6, desiredAccuracyPerChargeState=1e-11, desiredAccuracy=0):
        self.minCharge = minCharge
        self.tStep = tStep
        self.desiredAccuracyPerChargeState = desiredAccuracyPerChargeState
        self.desiredAccuracy = desiredAccuracy

class EbitParams:
    def __init__(self,
                 breedingTime=0.1,
                 probeEvery=0.1,
                 ionEbeamOverlap=1.0,
                 beamEnergy=3000.0,
                 beamCurrent=0.1,
                 beamRadius=200.0e-4,
                 pressure=1e-12,
                 ionTemperature=100.0,
                 electronVelocity=0.0,
                 ignoreBetaDecay=1,
                 currentDensity=0.0,
                 rkParams=None,
                 decayConstants=[]):  # I'm throwing decay constants in here so I can have those ready for the beta decay stuff without having to pass all the species objects
        self.breedingTime = breedingTime
        self.probeEvery = probeEvery
        self.ionEbeamOverlap = ionEbeamOverlap # unitless
        self.beamEnergy = beamEnergy # eV
        self.beamCurrent = beamCurrent # A
        self.beamRadius = beamRadius # cm^2
        self.pressure = pressure # Torr
        self.ionTemperature = ionTemperature # will get rid of
        self.electronVelocity = __C__*sqrt(1-(beamEnergy/__EMASS__+1)**-2) # cm/s
        self.currentDensity = beamCurrent / (pi*beamRadius**2) # A/cm^2
        self.rkParams = rkParams
        self.decayConstants = decayConstants
        self.ignoreBetaDecay = ignoreBetaDecay


def createEmptyList(sizeOfArray):
    myEmptyList = [0.0 for i in range(sizeOfArray)]
    return myEmptyList

def createEmptyListOfList(sizeOfArray):
    val = [[0.0, 0.0] for i in range(sizeOfArray)]
    return val


def createDecayConstants(betaHalfLife):
    if betaHalfLife <= 0:
        decayConstant = 0.0
    else:
        decayConstant = log(2) / betaHalfLife
    return decayConstant


def createChargeExchangeRates(Z, A, pressure, ionTemperature):
    """ We are not sure of where this formula for CX cross sections is derived from
    It compares the velocity of the electron with the velocity of a 1s electron in 
    the hydrogen atom. Looks similar to Lotz formula.
    """
    chargeExchangeRates = [0] * (Z + 1)

    # Need to return a Z+1 array
    h2Density = pressure * __TORR__ # using IGL, determine number density of H2 from the prescribed pressure
    ionMassInAMU = A * __AMU__
    avgIonV = __C__ * sqrt(8.0 * (ionTemperature / (pi * ionMassInAMU)))

    avgIonVinBohr = avgIonV / __VBOHR__
    sigV = __CHEXCONST__ * log( 15.0 / avgIonVinBohr) * avgIonV

    for i in range(1, Z + 1):
        chargeExchangeRates[i] = i * sigV * h2Density
    return chargeExchangeRates


def createChargeExchangeRates_MS(Z, A, pressure, ionTemperature):
    """ An implementation using the work of Salzborn and Mueller from 1977 paper.
    """

    chargeExchangeRates = [0] * (Z + 1)

    h2Density = pressure * __TORR__
    
    print("Number density of background gas... %s cm^-3"%'{:.2e}'.format(h2Density) )

    # Epgas is the ionization potential of the background gas. Make this adjustable?

    # The full cross section includes qi^alphak, but this is added in the loop below.
    sigV = __SALZBORNAK__*__Epgas__**__SALZBORNBETAK__
    for i in range(1, Z+1):
        chargeExchangeRates[i] = h2Density*sigV*i**__SALZBORNALPHAK__
    return chargeExchangeRates


def createInteractionRates(Z, beamEnergy, currentDensity, crossSections):
    # Interaction rate formulae for EII and RR:
    #       Rate = (number density)(electron velocity)(cross section)(number of reaction ions)(overlap function)
    #
    # This returns the "static" parts of the interaction rate which can be calculated before the time stepping.
    #
    # From Lisp Code :
    #  "Calculate the rate for a specific interaction (ionization/rr/..)
    #  inside the trap given the BEAM-ENERGY in eV, CURRENT-DENSITY in
    #  A/cm2, and the CROSS-SECTIONS in cm-2 corresponding to the
    #  interaction. Returns array of size Z-ION+1 with rates."
    interactionRate = [0] * (Z + 1)

    electronVelocity = __C__*sqrt(1-(beamEnergy/__EMASS__+1)**-2)
    # electronVelocity = __C__ * sqrt(2 * (beamEnergy / __EMASS__))

    # electron number density
    electronRate = currentDensity / __ECHG__ / electronVelocity

    for i in range(0, Z + 1):
        # number density X velocity X cross section = rate
        interactionRate[i] = crossSections[i] * electronVelocity * electronRate

    return interactionRate


def createDefaultInteractionRates(mySpecies, myEbitParams, ebitParams, crossSecFunc, runChargeExchange=0):  # move away from chex=0 and move into object
    
    interactionRates = createEmptyList(mySpecies.Z + 2)
    if runChargeExchange == 0:
        myFuncValues = createInteractionRates(mySpecies.Z, myEbitParams.beamEnergy, myEbitParams.currentDensity, crossSecFunc(myEbitParams.beamEnergy, mySpecies.Z))
    else:

        myFuncValues = crossSecFunc(mySpecies.Z, mySpecies.A, ebitParams[0].pressure, ebitParams[0].ionTemperature)

    # Reorganizes output for return
    for r in range(0, len(myFuncValues)):
        interactionRates[r] = myFuncValues[r]

    return interactionRates

def calculateIonHeating(mySpecies, Z, spitzerHeatingRates, tstep):
    retArray = [0.0 for i in range(0, Z+1)]

    for qindex in range(0, Z+1):
        retArray[qindex] = tstep*spitzerHeatingRates[qindex]*mySpecies.population[qindex]
    return retArray

def calculateJ(mySpecies, species, tmpEnergy, Z, spitzerHeatingRates, retJ, Ei, kfactor, weighting, tstep):
    """ 
    This will be used later for RK4 of the energy balance equation. The J is then analogous to the K's
    calculated for the RK4 algorithm.
    """

    # we iterate through charge states and use this to keep track of deltas
    delta = 0.0

    # this is only necessary for the ion-ion interaction and will need some modification
    for qindex in range(0,Z+1):
        # weighting is either 0, 0.5 or 1.0, depending on kfactor
        # jfactor is 0, j1, j2, j3 for producing j1, j2, j3, and j4.
        # Ei is the energy we start with
        tmpEnergy = Ei[qindex] + (weighting*kfactor)

    # implemented for Spitzer heating. Will need to change for later versions. Will also need the overlap factor
    for qindex in range(0,Z+1):
        retJ[qindex] = tstep*spitzerHeatingRates[qindex]*mySpecies.population[qindex]

    # I will use this later for when I need to calculate heat transfer between charge states
    # for qindex in range(0, Z):
    #     delta = tstep * (spitzerHeatingRates)
    #     retJ[qindex] = delta - lastDelta
    #     lastDelta = delta
    # retJ[Z] = -lastDelta

    return retJ

def calculateK(ebitParams, mySpecies, species, tmpPop, Z, ionizationRates, chargeExchangeRates, rrRates, retK, p1, p2, addWFactor, tstep):
    # k's are the increments used in the Runge Kutta 4 iterative method


    

    # Lengths here are Z+1 because we have a neutral charge state to account for.
    mySpecies.betaDecayDelta = [0.0]*(Z+1)
    nonDecayLastDelta = 0.0

    # tmpPop is the population at the beginning, midpoint or end of the interval.
    # this is used to estimate slopes for calculating k1, k2, k3, k4.
    # indexes through each value of ion charge state: 0 to Z
    for zindex in range(0, Z+1):
        # print("Z is %s" %str(mySpecies.Z))
        # print("addWFactor is %s" %str(addWFactor))
        # print("zindex is %s" %str(zindex))
        tmpPop[zindex] = p1[zindex] + (p2[zindex] * addWFactor)

    # This only indexes to Z-1 because the Z index is filled separately
    for zindex in range(0, Z):
        # Calculate changes in charge state populations:
        # dNi = dt * ( Rei(i-1) - Rei(i) + Rrr(i+1) - Rrr(i) + Rcx(i+1) - Rcx(i) )
        #     where Ni is population of charge state i
        #     = ionization rate of i-1 minus ionization rate of i and recombination rate of i+1 minus recombination rate of i

        # For each value of charge state q, only the rates between q and q+1 are calculated (not between q-1 and q). This
        # value is retained and used for the next step to account for rates between q-1 and q.
        nonDecayDelta = tstep * (- (        ionizationRates[zindex] * tmpPop[zindex]     )
                                 + (chargeExchangeRates[zindex + 1] * tmpPop[zindex + 1] )
                                 + (            rrRates[zindex + 1] * tmpPop[zindex + 1] ) )

        if ebitParams.ignoreBetaDecay != 1:  # ignore if we don't have any to speed things up
            mySpecies.betaDecayDelta[zindex] = mySpecies.tmpPop[zindex]*mySpecies.decayConstant*tstep
            # mySpecies.betaDecayDelta[zindex] = betaDecay1(mySpecies, species, ebitParams, zindex, tstep)


        # charge changing reactions
        retK[zindex] = (nonDecayDelta - nonDecayLastDelta)
        nonDecayLastDelta = nonDecayDelta


    # Filling of Z index
    retK[Z] = (-nonDecayLastDelta)


    # should this be tmpPop instead of mySpecies.tmpPop?
    mySpecies.betaDecayDelta[Z] = tmpPop[Z]*mySpecies.decayConstant*tstep

    # The portion deals with cross-species rates of beta decay. Determines if neighboring species affect each other.
    mySpeciesIndex = species.index(mySpecies)
    if mySpeciesIndex != 0 and mySpecies.Z - species[mySpeciesIndex-1].Z == 1:
        # first condition: if mySpeciesIndex==0, then nothing is populating our species with beta decay
        # second condition: if the next lower index of species is one unit of Z away, it populates mySpecies

        # this is confusing. To the original retK, we subtract the beta decay rate of mySpecies and add the beta decay rate
        # of species[mySpeciesIndex-1]. Since betaDecayDelta of the species at Z-1 is shorter in length by 1 we must add a [0] at beginning
        retK = [*map(lambda x: x[0]-x[1]+x[2], zip(retK, mySpecies.betaDecayDelta, [0]+species[mySpeciesIndex-1].betaDecayDelta))]
    else:
        # else we just do the regular and subtract the beta decay rate.
        retK = [x[0]-x[1] for x in zip(retK, mySpecies.betaDecayDelta)]

    # This will not work because the next species will not see it. Also, I don't have a conditional in here for checking
    # that the next species is neighboring in Z. The solution is to make this rate a parameter of the class. The best thing
    # about such a solution is that we dont need betaDecayLastDelta... just look at the previous species.

    return retK


def calculateKR(ebitParams, mySpecies, species, tmpPop, Z, ionizationRates, chargeExchangeRates, rrRates, retK, retR, p1, p2, addWFactor, tstep):
    """ 
    A new function from calculateK() so that we can allow heat transfer which occurs which changing charge state populations.
    """

    ionMassIneV = mySpecies.A*__AMU__
    # Lengths here are Z+1 because we have a neutral charge state to account for.
    mySpecies.betaDecayDelta = [0.0]*(Z+1)
    nonDecayLastDelta = 0.0

    # tmpPop is the population at the beginning, midpoint or end of the interval.
    # this is used to estimate slopes for calculating k1, k2, k3, k4 or the r's (for heat transfer)
    # indexes through each value of ion charge state: q=0 to q=Z
    for qindex in range(0, Z+1):
        tmpPop[qindex]    = p1[qindex] + (p2[qindex] * addWFactor)
        # tmpEnergy[qindex] = E1[qindex] + (E2[qindex] * addWFactor)

    # This only indexes to Z-1 because the Z index is filled separately
    for zindex in range(0, Z):
        # Calculate changes in charge state populations:
        # dNi = dt * ( Rei(i-1) - Rei(i) + Rrr(i+1) - Rrr(i) + Rcx(i+1) - Rcx(i) )
        #     where Ni is population of charge state i
        #     = ionization rate of i-1 minus ionization rate of i and recombination rate of i+1 minus recombination rate of i

        avgIonV = __C__*sqrt(8.0*mySpecies.NkT[zindex]/(pi*ionMassIneV))

        # For each value of charge state q, only the rates between q and q+1 are calculated (not between q-1 and q). This
        # value is retained and used for the next step to account for rates between q-1 and q.
        nonDecayDelta = tstep * (- (        ionizationRates[zindex] * tmpPop[zindex]     )
                                 + (chargeExchangeRates[zindex + 1] * tmpPop[zindex + 1] * avgIonV )
                                 + (            rrRates[zindex + 1] * tmpPop[zindex + 1] ) )

        if ebitParams.ignoreBetaDecay != 1:  # ignore if we don't have any to speed things up

            # changed from mySpecies.tmpPop to tmpPop... shouldn't have any adverse effects.
            mySpecies.betaDecayDelta[zindex] = tmpPop[zindex]*mySpecies.decayConstant*tstep
            # mySpecies.betaDecayDelta[zindex] = betaDecay1(mySpecies, species, ebitParams, zindex, tstep)

        # charge changing reactions
        retK[zindex] = (nonDecayDelta - nonDecayLastDelta)
        nonDecayLastDelta = nonDecayDelta


    # Filling of Z index
    retK[Z] = (-nonDecayLastDelta)

    # A new attempt for populating retR. Trying to reduce occurance of rounding errors blowing up.
    # It is on purpose that this only indexes up to Z-1
    for qindex in range(0,Z):
        # retR = [ amount lost by ionization of q(i),
        #          amount gained by RR and CX of q(i+1)]
        retR[qindex] = [ tstep*ionizationRates[qindex]*tmpPop[qindex], 
                         tstep*(chargeExchangeRates[qindex+1]+rrRates[qindex+1])*tmpPop[qindex+1]]



    mySpecies.betaDecayDelta[Z] = tmpPop[Z]*mySpecies.decayConstant*tstep

    # The portion deals with cross-species rates of beta decay. Determines if neighboring species affect each other.
    mySpeciesIndex = species.index(mySpecies)
    if mySpeciesIndex != 0 and mySpecies.Z - species[mySpeciesIndex-1].Z == 1:
        # first condition: if mySpeciesIndex==0, then nothing is populating our species with beta decay
        # second condition: if the next lower index of species is one unit of Z away, it populates mySpecies

        # this is confusing. To the original retK, we subtract the beta decay rate of mySpecies and add the beta decay rate
        # of species[mySpeciesIndex-1]. Since betaDecayDelta of the species at Z-1 is shorter in length by 1 we must add a [0] at beginning
        retK = [*map(lambda x: x[0]-x[1]+x[2], zip(retK, mySpecies.betaDecayDelta, [0]+species[mySpeciesIndex-1].betaDecayDelta))]
    else:
        # else we just do the regular and subtract the beta decay rate.
        retK = [x[0]-x[1] for x in zip(retK, mySpecies.betaDecayDelta)]

    return retK, retR


def rkStep(ebitParams, mySpecies, species, tstep, populationAtT0, populationAtTtstep, energyAtT0, energyAtTtstep):
    # longer function param calls yes but it speeds it up calculateK by 10%...
    """
    Here we are returning k's and r's, but it might be easier to just return the r's because the population change can still be
    calulated from this. In addition, rounding errors wouldn't be a problem anymore.
    """
    # print("\nRunning an RK step... ")

    mySpecies.k1, mySpecies.r1 = calculateKR(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k1, mySpecies.r1, populationAtT0, mySpecies.tmpPop, 0.0, tstep)
    mySpecies.k2, mySpecies.r2 = calculateKR(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k2, mySpecies.r2, populationAtT0,     mySpecies.k1, 0.5, tstep)
    mySpecies.k3, mySpecies.r3 = calculateKR(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k3, mySpecies.r3, populationAtT0,     mySpecies.k2, 0.5, tstep)
    mySpecies.k4, mySpecies.r4 = calculateKR(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k4, mySpecies.r4, populationAtT0,     mySpecies.k3, 1.0, tstep)
    
    # mySpecies.k1 = calculateK(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k1, populationAtT0, mySpecies.tmpPop, 0.0, tstep)
    # mySpecies.k2 = calculateK(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k2, populationAtT0,     mySpecies.k1, 0.5, tstep)
    # mySpecies.k3 = calculateK(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k3, populationAtT0,     mySpecies.k2, 0.5, tstep)
    # mySpecies.k4 = calculateK(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k4, populationAtT0,     mySpecies.k3, 1.0, tstep)
   
    # print("k values for q=1:")
    # print("k1 %s"%mySpecies.k1[1])
    # print("k2 %s"%mySpecies.k2[1])
    # print("k3 %s"%mySpecies.k3[1])
    # print("k4 %s"%mySpecies.k4[1])

    # Updates the population of each charge state in the species.
    for qindex in range(0, mySpecies.Z + 1):
        # new energy value     = ( kT(q-1)(pop gained by q-1) - kT(q)(lost by q) + kT(q+1)(gained by q+1) ) / total change in population 
        # populationAtTtstep[qindex] = populationAtT0[qindex] + ((1 / 6) * (sum(mySpecies.r1[qindex]) + (2 * sum(mySpecies.r2[qindex]) + sum(mySpecies.r3[qindex]) )  + sum(mySpecies.r4[qindex]) ))
        populationAtTtstep[qindex] = populationAtT0[qindex] + ((1 / 6) * (mySpecies.k1[qindex] + (2 * (mySpecies.k2[qindex] + mySpecies.k3[qindex])) + mySpecies.k4[qindex]) )

    # New calculation of time stepped energy
    # deltaPop = [loss by q(i) from EI, gain by q(i) from CX or RR]
    for q in range(0, mySpecies.Z+1):
        deltaPop = [(mySpecies.r1[q][i] + (2 * (mySpecies.r2[q][i] + mySpecies.r3[q][i])) + mySpecies.r4[q][i])/6 for i in range(0,2)]
        # print("DeltaPop for q=%s"%q+": %s"%deltaPop)
        if q==0:
            try:
                #this one is with gain only...
                energyAtTtstep[q] = (energyAtT0[q]*(populationAtT0[q]) + energyAtT0[q+1]*deltaPop[1]) / (populationAtT0[q]+deltaPop[1])
                #this one is with gain and loss... caused problems
                # energyAtTtstep[q] = (energyAtT0[q]*(populationAtT0[q]-deltaPop[0]) + energyAtT0[q+1]*deltaPop[1]) / (populationAtTtstep[q])
            except ZeroDivisionError:
                energyAtTtstep[q] = energyAtT0[q]
        elif q==mySpecies.Z:
            lowerQ = [(mySpecies.r1[q-1][i] + (2 * (mySpecies.r2[q-1][i] + mySpecies.r3[q-1][i])) + mySpecies.r4[q-1][i])/6 for i in range(0,2)]
            try:
                #gain only
                energyAtTtstep[q] = (energyAtT0[q]*(populationAtT0[q]) + energyAtT0[q-1]*lowerQ[0]) / (populationAtT0[q]+deltaPop[1]+lowerQ[0])
                # gain and loss
                # energyAtTtstep[q] = (energyAtT0[q]*(populationAtT0[q]-deltaPop[0]) + energyAtT0[q-1]*lowerQ[0]) / (populationAtTtstep[q])
            except ZeroDivisionError:
                energyAtTtstep[q] = energyAtT0[q]
        else:
            lowerQ = [(mySpecies.r1[q-1][i] + (2 * (mySpecies.r2[q-1][i] + mySpecies.r3[q-1][i])) + mySpecies.r4[q-1][i])/6 for i in range(0,2)]
            # print("lowerQ: %s"%lowerQ)
            try:
                #gain
                # print("energyAtT0[q-1] = %s"%energyAtT0[q-1] + ", lowerQ[0] = %s"%lowerQ[0]+", populationAtT0[q]=%s"%populationAtT0[q]+", deltaPop[1]=%s"%deltaPop[1])
                energyAtTtstep[q] = (energyAtT0[q]*(populationAtT0[q]) + energyAtT0[q-1]*lowerQ[0] + energyAtT0[q+1]*deltaPop[1]) / (populationAtT0[q]+deltaPop[1]+lowerQ[0])
                #gain and loss
                # energyAtTtstep[q] = (energyAtT0[q]*(populationAtT0[q]-deltaPop[0]) + energyAtT0[q-1]*lowerQ[0] + energyAtT0[q+1]*deltaPop[1]) / (populationAtTtstep[q])
            except ZeroDivisionError:
                energyAtTtstep[q] = energyAtT0[q]

    
    # print("Initial pop: %s"%populationAtT0 + ",\nfinal pop:   %s"%populationAtTtstep)
    # print("Initial temp: %s"%energyAtT0 + ",\nfinal temp:   %s"%energyAtTtstep)
    return

def rkStepEnergy(ebitParams, mySpecies, species, tstep, energyAtT0):
    """ This is currently implemented as an RK step that doesn't get included in consideration for the adaptive time step...
    although that is a mistake.

    It take in energyAtT0, and then calculates how much heat was transferred to the charge state q(i) by three avenues:
        Heat given by q(i-1)
        Heat given by q(i+1)
        Heat lost by q(i)
    """
    return 


def probeFnAddPop(time, ebitParams, mySpecies, tstep):
    #  In the original lisp code this function was passed along to calcChargePopulations
    #  so for now I will do the same as apparently this is something someone might want
    #  to change to a different function... If someone knows the reason for this please
    #  let me know.
    for chargeIndex in range(0, len(mySpecies.chargeStates)):
        if len(mySpecies.results) < len(mySpecies.chargeStates):
            mySpecies.results.append([]) # Pre-allocate the size of .results
        mySpecies.results[chargeIndex].append([time, mySpecies.population[mySpecies.chargeStates[chargeIndex]], mySpecies.NkT[mySpecies.chargeStates[chargeIndex]], tstep, mySpecies.stepCounter])



def adaptiveRkStepper(species, ebitParams, probeFnAddPop):
    # This is an attempt at restructuring adaptiveRkStepper() so that it can perform beta decay much easier. -Zach
    # The essential difference is that I will not be looping through a species, rather performing the time advance
    # for ALL species in one iteration.
    # To accommodate for the adaptive stepping, I calculate the bestStepSize of each species and use the limiting
    # value for the next iteration. Therefore the rate of change of one population can limit the whole calculation.
    # It is slightly slower than the original algorithm and it allows for us to better implement interactions
    # between the different species in the EBIT.


    t = 0.0
    bestStepSize = 0 # This is an overall value that is used after the mySpecies.bestStepSize's are compared.
    nextPrint = 0.0
    step = ebitParams[0].rkParams.tStep
    for mySpecies in species:
        # mySpecies.noTooBigSteps is already defaulted to 0.0 upon instantiation
        mySpecies.bestStepSize = step
    timeToBreed = 0.0

    # f"{Decimal('40800000000.00000000000000'):.2E}"
    # '{:.2e}'.format(x)

    for myEbitParams in ebitParams:
        print("Simulating using beam energy %s eV" % myEbitParams.beamEnergy)
        print("Electron velocity ... %s c"%str(myEbitParams.electronVelocity/__C__))
        print("Electron number density ... %s cm^-3"%'{:.2e}'.format(myEbitParams.currentDensity/(myEbitParams.electronVelocity*__ECHG__)) )
        print("Current density ... %s A/cm^2"%'{:.2e}'.format(myEbitParams.currentDensity))

        print("Tracking species Z = %s" % ', '.join(map(lambda x: str(x.Z), species) ) )
        if myEbitParams.ignoreBetaDecay == 0:
            # the list below is incredibly ugly, but it just gets the Z's of species which have nonzero decay constants.
            betaDecayList = [species[i].Z for i,j in enumerate(map(lambda x: x.decayConstant, species)) if j!=0.0]
            print("Detected nonzero beta decay constants for Z = %s" %', '.join(map(str, betaDecayList )) )
            print("\n***Please note that the inclusion of beta decay physics significantly increases computation time!***")

        timeToBreed = timeToBreed + myEbitParams.breedingTime

        # Calculate the static portions of the interaction rates. New one for each new EBIT configuration
        for mySpecies in species:
            calcRateMatrices(mySpecies, myEbitParams, ebitParams)

            if mySpecies.initSCITemp != -1:
                print("initSCITemp detected... including energy dynamics")
                # calculates the static portion of the energy rates.... similar to what is done with calcRateMatrices()
                calcEnergyRates(mySpecies, myEbitParams, ebitParams)
                # initialize energy array with the total SCI energy
                mySpecies.NkT[1] = mySpecies.initSCITemp

        # Enter rk stepping loop
        while t <= timeToBreed:

            # Compute population change for all species in the EBIT.
            for mySpecies in species:

                # See Wikipedia page for "Adaptive stepsize" for an explanation of how this works
                # The last input parameters y1, y12 and y22 are actually populations which are being calculated.
                # y1 is calculated using the derivatives at t0 and evolving forward by 2*tstep
                # y12 is calculated using the derivatives at t0 and evolving forward by tstep
                # y22 is calculated using derivatives of y12 at t0 + tstep and evolves forward by tstep, resulting in a total of 2*tstep.
                # Therfore y22 is a more accurate estimation at t0 + 2*tstep than y1.

                # start = time.time()
                # print("\n-----------------Starting first RK step-----------------")
                #                                                   initial population,   extrap. pop,  initial energy, extrap. energy
                rkStep(ebitParams[0], mySpecies, species, 2 * step, mySpecies.population, mySpecies.y1 , mySpecies.NkT, mySpecies.f1)
                rkStep(ebitParams[0], mySpecies, species,     step, mySpecies.population, mySpecies.y12, mySpecies.NkT, mySpecies.f12)
                rkStep(ebitParams[0], mySpecies, species,     step, mySpecies.y12,        mySpecies.y22, mySpecies.f12, mySpecies.f22)


                # rkStep(ebitParams[0], mySpecies, species, 2 * step, mySpecies.population, mySpecies.y1 )
                # rkStep(ebitParams[0], mySpecies, species,     step, mySpecies.population, mySpecies.y12)
                # rkStep(ebitParams[0], mySpecies, species,     step, mySpecies.y12,        mySpecies.y22)
                # print("-----------------Finished with last RK step--------------\n")


                # end = time.time()
                # print("t = %s"%str(t)+", duration of Z=%s "%str(mySpecies.Z)+ ", 3 rkStep calculations: %s" %str(end-start))

                # Calculates for EACH species.
                # y22 is objectively a better estimate than y1.
                # If the sum of absolute differences between y1 and y22 is nonzero, modify bestStepSize to be smaller
                # than step. This means that step was too large.
                mySpecies.diff = sum([abs(x - y) for (x, y) in zip(mySpecies.y1, mySpecies.y22)])

                # Keep track of global truncation error
                mySpecies.truncationError += mySpecies.diff

            # If ANY of the diff's are >0, then we perform the bestStepSize calculation
            if any(i>0 for i in map(lambda x: x.diff, species)):
                # If y22 != y1, we determine the best step size.
                # print("wasn't the perfect step size... adjusting bestStepSize")
                for mySpecies in species:
                    # print("original step size: %s"%mySpecies.bestStepSize)
                    if mySpecies.diff != 0:
                        mySpecies.bestStepSize = step * ((ebitParams[0].rkParams.desiredAccuracy / mySpecies.diff) ** 0.2)
                        # print("new step size: %s"%mySpecies.bestStepSize)
            else:
                # We pretty much never enter this condition
                # If y1 == y22 for all species (unlikely), then we chose the perfect step size.
                # we commit y22 and continue to the next step
                t += 2*step

                for mySpecies in species:
                    mySpecies.stepCounter += 1
                    # print("initial temperatur is %s"%mySpecies.NkT)
                    # print("initial population is %s"%mySpecies.population)
                    # print("final   temperatur is %s"%mySpecies.f22)
                    # print("final   population is %s"%mySpecies.y22)
                    # print("\n")
                    if mySpecies.initSCITemp != -1:
                        mySpecies.NkT = copy.copy(mySpecies.f22)
                    # print(min([k for k in [abs(i-j) for i,j in zip(mySpecies.y22,mySpecies.population)] if k!=0.0] ))
                    mySpecies.population = copy.copy(mySpecies.y22)
                    

                # Continue means we start again at the beginning of the species loop!
                continue

            # If we get here, the first if statement was entered and bestStepSize has been adjusted.
            # If ALL of the bestStepSizes are greater than step, we commit y22, increment time, and adjust step
            if all(i >= step for i in map(lambda x: x.bestStepSize, species)):
                t += 2*step
                bestStepSize = min(map(lambda x: x.bestStepSize, species))
                for mySpecies in species:
                    mySpecies.stepCounter += 1
                    # print("initial temperatur is %s"%mySpecies.NkT)
                    # print("initial population is %s"%mySpecies.population)
                    # print("final   temperatur is %s"%mySpecies.f22)
                    # print("final   population is %s"%mySpecies.y22)
                    # print("\n")
                    # if min(mySpecies.f22) < 0.0:
                    #     print("min of mySpecies.f22 is %s"%min(mySpecies.f22)+", at q is %s"%mySpecies.f22.index(min(mySpecies.f22)))
                    #     sys.exit("Energy became negative... ")
                    if mySpecies.initSCITemp != -1:
                        mySpecies.NkT = copy.copy(mySpecies.f22)
                    # print(min([k for k in [abs(i-j) for i,j in zip(mySpecies.y22,mySpecies.population)] if k!=0.0]) )
                    mySpecies.population = copy.copy(mySpecies.y22)
            else:
                # If we get here, one of the mySpecies.bestStepSize values is smaller than step
                # and we need to re-perform this time loop with an adjusted step value.
                # We do NOT increment the time
                # print("Step rejected, it was too small...")
                bestStepSize = min(map(lambda x: x.bestStepSize, species))
                mySpecies.noTooBigSteps += 1

            # If we get here, one of the mySpecies.bestStepSize's was either too large or too small. In the statements above, depending on
            # largeness or smallness, we assigned the correct value to bestStepSize.
            step = 0.9 * bestStepSize

            # Probing the population of each species
            if t >= nextPrint:
                # print("Probing... %s"%t)
                nextPrint += ebitParams[0].probeEvery
                # print("noTooBigSteps: %s"%str(mySpecies.noTooBigSteps))
                for mySpecies in species:
                    probeFnAddPop(t, ebitParams[0], mySpecies, step)
        # print("Final step size: %s"%str(step))

    return

def calcEnergyRates(mySpecies, myEbitParams, ebitParams):

    # charge states is range(0, mySpecies.Z + 1)
    # For each charge state of the species, we need a rate calculation for Spitzer heating
    mySpecies.spitzerHeatingRates = rateSpitzerHeating(mySpecies, myEbitParams)


def calcRateMatrices(mySpecies, myEbitParams, ebitParams):
    print("Calculating rate matrices...")
    """ These here are what I would call the "static" portion of the rates. They are calculated before the time stepping occurs because
    the information they hold is not dynamic.
    """
    # myEbitParams.currentDensity   = (ebitParams[0].ionEbeamOverlap * ebitParams[0].beamCurrent) / (pi * (ebitParams[0].beamRadius ** 2))
    # removed above to calculate it during the class instantiation
    mySpecies.ionizationRates     = createDefaultInteractionRates(mySpecies, myEbitParams, ebitParams, createIonizationCrossSections   )
    mySpecies.rrRates             = createDefaultInteractionRates(mySpecies, myEbitParams, ebitParams,         createRRCrossSections   )
    mySpecies.chargeExchangeRates = createDefaultInteractionRates(mySpecies, myEbitParams, ebitParams,     createChargeExchangeRates_MS, 1)


def initEverything(species, ebitParams):
    global INITILIZEDEVERYTHING
    ebitParams[0].rkParams = RkStepParams()  # We're just going to shove this information into the first beam because why not?

    ebitParams[0].rkParams.desiredAccuracy = ebitParams[0].rkParams.desiredAccuracyPerChargeState / species[0].Z

    species.sort(key=lambda x: x.Z, reverse=True)  # Sort species list by Z in decending order
    largestZValue = species[0].Z  # We need the largest value for the size of the population array

    species.reverse() # For beta decay calculations, species must be ascending in Z

    # if species[0].betaHalfLife != 0.0:  # The largest Z value species is not allowed to have a beta decay
    #     raise ValueError("Can not handle having a beta decay for highest Z species")
        

    for mySpecies in species:
        # Initilize everything ...
        mySpecies.population = createEmptyList(largestZValue + 2) # I forget why I needed this...
        mySpecies.population[1] = mySpecies.initSCIPop
        mySpecies.decayConstant = createDecayConstants(mySpecies.betaHalfLife)

        if mySpecies.decayConstant > 0:  # If we don't need it we will use this to ignore the beta calculation function
            ebitParams[0].ignoreBetaDecay = 0

        ebitParams[0].decayConstants.append(mySpecies.decayConstant)

        # These arrays are dynamic and only hold instantaneous information for calculations of the population
        # of each charge state of a species.
        mySpecies.k1     = createEmptyList(mySpecies.Z + 2)
        mySpecies.k2     = createEmptyList(mySpecies.Z + 2)
        mySpecies.k3     = createEmptyList(mySpecies.Z + 2)
        mySpecies.k4     = createEmptyList(mySpecies.Z + 2)
        mySpecies.tmpPop = createEmptyList(mySpecies.Z + 2)
        mySpecies.tmpEnergy = createEmptyList(mySpecies.Z + 2)
        mySpecies.NkT    = createEmptyList(mySpecies.Z + 2)
        mySpecies.y1     = createEmptyList(mySpecies.Z + 2)
        mySpecies.y12    = createEmptyList(mySpecies.Z + 2)
        mySpecies.y22    = createEmptyList(mySpecies.Z + 2)
        mySpecies.f1     = createEmptyList(mySpecies.Z + 2)
        mySpecies.f12    = createEmptyList(mySpecies.Z + 2)
        mySpecies.f22    = createEmptyList(mySpecies.Z + 2)
        # these need to be list of list because we separate each k into contributions from q-1, q, and q+1
        mySpecies.r1     = createEmptyListOfList(mySpecies.Z + 2)
        mySpecies.r2     = createEmptyListOfList(mySpecies.Z + 2)
        mySpecies.r3     = createEmptyListOfList(mySpecies.Z + 2)
        mySpecies.r4     = createEmptyListOfList(mySpecies.Z + 2)

        mySpecies.results = []
    INITILIZEDEVERYTHING = 1


def calcChargePopulations(species, ebitParams, probeFnAddPop):
    if INITILIZEDEVERYTHING == 0:
        initEverything(species, ebitParams)

    start = time.time()
    adaptiveRkStepper(species, ebitParams, probeFnAddPop)  # this does all the heavy lifting
    end = time.time()
    for mySpecies in species:
        print("------------------------------------------------------- \n")
        print("Species diagnostics for Z = %s:" %str(mySpecies.Z))
        print("    Checking final EBIT population of species... %s" %str(sum([i[-1][1] for i in mySpecies.results]) ))
        print("    Total number of time steps used... %s" %str(mySpecies.stepCounter) )
        print("    Number of bad step guesses... %s" %str(mySpecies.noTooBigSteps))
        print("    Global truncation error... %s" %str(mySpecies.truncationError))
        print("------------------------------------------------------- \n")

    print("Clock time for adaptive RK4 algorithm: %s seconds" %str(end-start))
    print("Simulation finished...")
    return
