#!/usr/bin/env python
from math import *
from IonizationEnergies import *
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
INITILIZEDEVERYTHING = 0


class Species:
    def __init__(self,
                 Z,
                 A,
                 decaysTo=0.0,
                 betaHalfLife=0.0,
                 initSCIPop=1.0,
                 chargeStates=None,
                 halfLife=0.0,
                 populationNumber=0.0,
                 population=None,
                 decayConstant=0.0,
                 ionizationRates=None,
                 rrRates=None,
                 chargeExchangeRates=None,
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
                 tmpPop=[],
                 y1=[],
                 y12=[],
                 y22=[],
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
        self.decayConstant = decayConstant
        self.chargeStates = chargeStates
        self.diff = diff
        self.bestStepSize = bestStepSize
        self.noTooBigSteps = noTooBigSteps
        self.results = results

# I think I need to move chargeStates up to here.. could concievibly do that for a bunch
# of the stuff in chargepopparams due to maybe they should be following the species around


class RkStepParams:
    def __init__(self, minCharge=5e-5, tStep=1e-6, desiredAccuracyPerChargeState=1e-7, desiredAccuracy=0):
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
                 ignoreBetaDecay=1,
                 currentDensity=None,
                 rkParams=None,
                 decayConstants=[]):  # I'm throwing decay constants in here so I can have those ready for the beta decay stuff without having to pass all the species objects
        self.breedingTime = breedingTime
        self.probeEvery = probeEvery
        self.ionEbeamOverlap = ionEbeamOverlap
        self.beamEnergy = beamEnergy
        self.beamCurrent = beamCurrent
        self.beamRadius = beamRadius
        self.pressure = pressure
        self.ionTemperature = ionTemperature
        self.currentDensity = currentDensity
        self.rkParams = rkParams
        self.decayConstants = decayConstants
        self.ignoreBetaDecay = ignoreBetaDecay


def createEmptyList(sizeOfArray):
    myEmptyList = [0.0 for i in range(sizeOfArray)]
    return myEmptyList


def createDecayConstants(betaHalfLife):
    if betaHalfLife <= 0:
        decayConstant = 0.0
    else:
        decayConstant = log(2) / betaHalfLife
    return decayConstant


def createChargeExchangeRates(Z, A, pressure, ionTemperature):
    # We are not sure of where this formula for CX cross sections is derived from
    # possibly might change to formulation of Salzborn & Mueller 1977 work.
    chargeExchangeRates = [0] * (Z + 1)

    # Need to return a Z+1 array
    h2Density = pressure * __TORR__
    ionMassInAMU = A * __AMU__
    avgIonV = __C__ * sqrt(8.0 * (ionTemperature / (pi * ionMassInAMU)))
    avgIonVinBohr = avgIonV / __VBOHR__
    sigV = __CHEXCONST__ * log( 15.0 / avgIonVinBohr) * avgIonV

    for i in range(1, Z + 1):
        chargeExchangeRates[i] = i * sigV * h2Density
    return chargeExchangeRates


def createInteractionRates(Z, beamEnergy, currentDensity, crossSections):
    # From Lisp Code :
    #  "Calculate the rate for a specific interaction (ionization/rr/..)
    #  inside the trap given the BEAM-ENERGY in eV, CURRENT-DENSITY in
    #  A/cm2, and the CROSS-SECTIONS in cm-2 corresponding to the
    #  interaction. Returns array of size Z-ION+1 with rates."
    interactionRate = [0] * (Z + 1)

    electronVelocity = __C__ * sqrt(2 * (beamEnergy / __EMASS__))
    electronRate = currentDensity / __ECHG__ / electronVelocity

    for i in range(0, Z + 1):
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




def calculateK(ebitParams, mySpecies, species, tmpPop, Z, ionizationRates, chargeExchangeRates, rrRates,  retK, p1, p2, addWFactor, tstep):
    # K's are the increments used in the Runge Kutta 4 iterative method

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
#                                 + (    ionizationRates[zindex - 1] * tmpPop[zindex - 1] )
#                                 - (    chargeExchangeRates[zindex] * tmpPop[zindex]     )
#                                 - (                rrRates[zindex] * tmpPop[zindex]     )
                                 + (chargeExchangeRates[zindex + 1] * tmpPop[zindex + 1] )
                                 + (            rrRates[zindex + 1] * tmpPop[zindex + 1] ) )

        if ebitParams.ignoreBetaDecay != 1:  # ignore if we don't have any to speed things up
            mySpecies.betaDecayDelta[zindex] = mySpecies.tmpPop[zindex]*mySpecies.decayConstant*tstep
            # mySpecies.betaDecayDelta[zindex] = betaDecay1(mySpecies, species, ebitParams, zindex, tstep)

        retK[zindex] = (nonDecayDelta - nonDecayLastDelta)
        nonDecayLastDelta = nonDecayDelta

    # Filling of Z index
    retK[Z] = (-nonDecayLastDelta)
    mySpecies.betaDecayDelta[Z] = mySpecies.tmpPop[Z]*mySpecies.decayConstant*tstep


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


def rkStep(ebitParams, mySpecies, species, tstep, populationAtT0, populationAtTtstep):
    # longer function param calls yes but it speeds it up calculateK by 10%...
    mySpecies.k1 = calculateK(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k1, populationAtT0, mySpecies.tmpPop, 0.0, tstep)
    mySpecies.k2 = calculateK(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k2, populationAtT0, mySpecies.k1,     0.5, tstep)
    mySpecies.k3 = calculateK(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k3, populationAtT0, mySpecies.k2,     0.5, tstep)
    mySpecies.k4 = calculateK(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k4, populationAtT0, mySpecies.k3,     1.0, tstep)

    # Updates the population of each charge state in a species. From 0 to Z
    for zindex in range(0, mySpecies.Z + 1):
        populationAtTtstep[zindex] = populationAtT0[zindex] + ((1 / 6) * (mySpecies.k1[zindex] + (2 * (mySpecies.k2[zindex] + mySpecies.k3[zindex])) + mySpecies.k4[zindex]) )
    return



def probeFnAddPop(time, ebitParams, mySpecies):
    #  In the original lisp code this function was passed along to calcChargePopulations
    #  so for now I will do the same as apparently this is something someone might want
    #  to change to a different function... If someone knows the reason for this please
    #  let me know.
    for chargeIndex in range(0, len(mySpecies.chargeStates)):
        if len(mySpecies.results) < len(mySpecies.chargeStates):
            mySpecies.results.append([]) # Pre-allocate the size of .results
        mySpecies.results[chargeIndex].append([time, mySpecies.population[mySpecies.chargeStates[chargeIndex]]])



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

    for myEbitParams in ebitParams:
        print("Simulating using beam energy %s" % myEbitParams.beamEnergy)
        print("Tracking species Z = %s" % ', '.join(map(lambda x: str(x.Z), species) ) )
        if myEbitParams.ignoreBetaDecay == 0:
            # the list below is incredibly ugly, but it just gets the Z's of species which have nonzero decay constants.
            betaDecayList = [species[i].Z for i,j in enumerate(map(lambda x: x.decayConstant, species)) if j!=0.0]
            print("Detected nonzero beta decay constants for Z = %s" %', '.join(map(str, betaDecayList )) )
            print("\n***Please note that the inclusion of beta decay physics significantly increases computation time!***")

        timeToBreed = timeToBreed + myEbitParams.breedingTime
        for mySpecies in species:
            calcRateMatrices(mySpecies, myEbitParams, ebitParams)

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
                rkStep(ebitParams[0], mySpecies, species, 2 * step, mySpecies.population, mySpecies.y1 )
                rkStep(ebitParams[0], mySpecies, species,     step, mySpecies.population, mySpecies.y12)
                rkStep(ebitParams[0], mySpecies, species,     step, mySpecies.y12,        mySpecies.y22)
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
                for mySpecies in species:
                    if mySpecies.diff != 0:
                        mySpecies.bestStepSize = step * ((ebitParams[0].rkParams.desiredAccuracy / mySpecies.diff) ** 0.2)
            else:
                # If y1 == y22 for all species (unlikely), then we commit the y22 population of every species and continue without adjusting the step size.
                # But we DO increment the time... but I can't do it here because then time would change between
                # species.
                t += 2*step

                for mySpecies in species:
                    mySpecies.stepCounter += 1
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
                    mySpecies.population = copy.copy(mySpecies.y22)
            else:
                # If we get here, one of the mySpecies.bestStepSize values is smaller than step
                # and we need to re-perform this time loop with an adjusted step value.
                # We do NOT increment the time
                bestStepSize = min(map(lambda x: x.bestStepSize, species))
                mySpecies.noTooBigSteps += 1

            # If we get here, one of the mySpecies.bestStepSize's was either too large or too small. In the statements above, depending on
            # largeness or smallness, we assigned the correct value to bestStepSize.
            step = 0.9 * bestStepSize

            # Probing the population of each species
            if t >= nextPrint:
                nextPrint += ebitParams[0].probeEvery
                for mySpecies in species:
                    probeFnAddPop(t, ebitParams[0], mySpecies)

    return

def calcRateMatrices(mySpecies, myEbitParams, ebitParams):
        myEbitParams.currentDensity   = (ebitParams[0].ionEbeamOverlap * ebitParams[0].beamCurrent) / (pi * (ebitParams[0].beamRadius ** 2))
        mySpecies.ionizationRates     = createDefaultInteractionRates(mySpecies, myEbitParams, ebitParams, createIonizationCrossSections   )
        mySpecies.rrRates             = createDefaultInteractionRates(mySpecies, myEbitParams, ebitParams,         createRRCrossSections   )
        mySpecies.chargeExchangeRates = createDefaultInteractionRates(mySpecies, myEbitParams, ebitParams,     createChargeExchangeRates, 1)


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
        mySpecies.population = createEmptyList(largestZValue + 2)
        mySpecies.population[0] = mySpecies.initSCIPop
        mySpecies.decayConstant = createDecayConstants(mySpecies.betaHalfLife)

        if mySpecies.decayConstant > 0:  # If we don't need it we will use this to ignore the beta calculation function
            ebitParams[0].ignoreBetaDecay = 0

        ebitParams[0].decayConstants.append(mySpecies.decayConstant)

        # These arrays are dynamic and only hold instantaneous information for calculations of the population
        # of each charge state of a species.
        mySpecies.k1 = createEmptyList(mySpecies.Z + 2)
        mySpecies.k2 = createEmptyList(mySpecies.Z + 2)
        mySpecies.k3 = createEmptyList(mySpecies.Z + 2)
        mySpecies.k4 = createEmptyList(mySpecies.Z + 2)
        mySpecies.tmpPop = createEmptyList(mySpecies.Z + 2)
        mySpecies.y1 = createEmptyList(mySpecies.Z + 2)
        mySpecies.y12 = createEmptyList(mySpecies.Z + 2)
        mySpecies.y22 = createEmptyList(mySpecies.Z + 2)

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
        print("    Checking final EBIT population of species... %s" %str(sum([i[-1][-1] for i in mySpecies.results]) ))
        print("    Total number of time steps used... %s" %str(mySpecies.stepCounter) )
        print("    Number of bad step guesses... %s" %str(mySpecies.noTooBigSteps))
        print("    Global truncation error... %s" %str(mySpecies.truncationError))
        print("------------------------------------------------------- \n")

    print("Clock time for adaptive RK4 algorithm: %s seconds" %str(end-start))
    print("Simulation finished...")
    return
