#!/usr/bin/env python
from math import *
from IonizationEnergies import *
from rr import *
import copy
import sys

__EMASS__ = 5.11e5  # "Electron mass in eV"
__C__ = 3.0e10  # "Speed of light"
__ECHG__ = 1.6e-19  # "Electron charge"
__VBOHR__ = 2.2e8
__AMU__ = 9.311e8  # "1 AMU in eV"
__TORR__ = 3.537e16  # "1 torr in cm-3"
__ALPHA__ = 7.2974e-3  # "fine-structure constant"
__CHEXCONST__ = 2.25e-16  # "Constant for use in charge-exchange "
__LCONV__ = 3.861e-11  # "length conversion factor to CGS"

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
                 population=None,
                 decayConstant=0.0,
                 ionizationRates=None,
                 rrRates=None,
                 chargeExchangeRates=None,
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
        self.initSCIPop = initSCIPop
        self.decayConstant = decayConstant
        self.chargeStates = chargeStates
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


def createDefaultInteractionRates(mySpecies, ebitParams, crossSecFunc, runChargeExchange=0):  # move away from chex=0 and move into object

    interactionRates = createEmptyList(mySpecies.Z + 2)
    if runChargeExchange == 0:
        myFuncValues = createInteractionRates(mySpecies.Z, ebitParams.beamEnergy, ebitParams.currentDensity, crossSecFunc(ebitParams.beamEnergy, mySpecies.Z))
    else:
        myFuncValues = crossSecFunc(mySpecies.Z, mySpecies.A, ebitParams.pressure, ebitParams.ionTemperature)

    for r in range(0, len(myFuncValues)):
        interactionRates[r] = myFuncValues[r]

    return interactionRates


def betaDecay(mySpecies, species, ebitParams, zindex, tstep):
    mySpeciesIndex = species.index(mySpecies)
    if zindex >= 1:
        try:
            myval = ((-1 * ebitParams.decayConstants[mySpeciesIndex] * mySpecies.tmpPop[zindex])
                     + (ebitParams.decayConstants[mySpeciesIndex + 1] * species[mySpeciesIndex + 1].tmpPop[zindex]))
        except IndexError:   # *** PLEASE MAKE SURE THIS IS A DECENT ASSUMPTION OF WHAT I SHOULD DO... TEST ME! ***
            myval = (-1 * ebitParams.decayConstants[mySpeciesIndex] * mySpecies.tmpPop[zindex])
    else:
        myval = 0
    return tstep * myval


def calculateK(ebitParams, mySpecies, species, tmpPop, Z, ionizationRates, chargeExchangeRates, rrRates,  retK, p1, p2, addWFactor, tstep):
    betaDecayDelta = 0
    nonDecayLastDelta = 0.0

    for zindex in range(0, Z):
        tmpPop[zindex] = p1[zindex] + (p2[zindex] * addWFactor)

    for zindex in range(0, Z):
        # Calculate chargeChanges
        nonDecayDelta = tstep * ((-1 * ionizationRates[zindex] * tmpPop[zindex])
                                 + ((chargeExchangeRates[zindex + 1] + rrRates[zindex + 1]) * tmpPop[zindex + 1]))

        if ebitParams.ignoreBetaDecay != 1:  # ignore if we don't have any to speed things up
            betaDecayDelta = betaDecay(mySpecies, species, ebitParams, zindex, tstep)

        retK[zindex] = (nonDecayDelta - nonDecayLastDelta) + betaDecayDelta
        nonDecayLastDelta = nonDecayDelta
        retK[Z] = (-nonDecayLastDelta) + betaDecayDelta

    return retK


def rkStep(ebitParams, mySpecies, species, tstep, populationAtT0, populationAtTtstep):
    # longer function param calls yes but it speeds it up calculateK by 10%...
    mySpecies.k1 = calculateK(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k1, populationAtT0, mySpecies.tmpPop, 0.0, tstep)
    mySpecies.k2 = calculateK(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k2, populationAtT0, mySpecies.k1, 0.5, tstep)
    mySpecies.k3 = calculateK(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k3, populationAtT0, mySpecies.k2, 0.5, tstep)
    mySpecies.k4 = calculateK(ebitParams, mySpecies, species, mySpecies.tmpPop, mySpecies.Z, mySpecies.ionizationRates, mySpecies.chargeExchangeRates, mySpecies.rrRates, mySpecies.k4, populationAtT0, mySpecies.k3, 1.0, tstep)

    for zindex in range(0, mySpecies.Z):
        populationAtTtstep[zindex] = populationAtT0[zindex] + ((1 / 6) * (mySpecies.k1[zindex] + (2 * (mySpecies.k2[zindex] + mySpecies.k3[zindex])) + mySpecies.k4[zindex]))
    return


def probeFnAddPop(time, ebitParams, mySpecies):
    #  In the original lisp code this function was passed along to calcChargePopulations
    #  so for now I will do the same as apparently this is something someone might want
    #  to change to a different function... If someone knows the reason for this please
    #  let me know.
    for chargeIndex in range(0, len(mySpecies.chargeStates)):
        if len(mySpecies.results) < len(mySpecies.chargeStates):
            mySpecies.results.append([])
        mySpecies.results[chargeIndex].append([time, mySpecies.population[mySpecies.chargeStates[chargeIndex]]])


def adaptiveRkStepper(species, ebitParams, probeFnAddPop):

    for mySpecies in species:

        time = 0.0  #  We are going to need to adjust the time to pick up where it left off instead of going back here? Or do the loop here?
        # or or or ????
        nextPrint = 0.0
        noTooBigSteps = 0
        step = ebitParams[0].rkParams.tStep
        bestStepSize = step
        timeToBreed = 0.0
        # For loop for beam changes here?
        for myEbitParams in ebitParams:
            print("Simulating using beam energy %s" % myEbitParams.beamEnergy)
            timeToBreed = timeToBreed + myEbitParams.breedingTime
            calcRateMatrices(mySpecies, myEbitParams)

            print("Simulating Species : %s" % mySpecies.Z)
     #       sys.exit(0)
            while time <= timeToBreed:
                if time >= nextPrint:
                    nextPrint = nextPrint + ebitParams[0].probeEvery
                    probeFnAddPop(time, ebitParams[0], mySpecies)

                rkStep(ebitParams[0], mySpecies, species, 2 * step, mySpecies.population, mySpecies.y1)
                rkStep(ebitParams[0], mySpecies, species, step, mySpecies.population, mySpecies.y12)
                rkStep(ebitParams[0], mySpecies, species, step, mySpecies.y12, mySpecies.y22)

                diff = sum([abs(x - y) for (x, y) in zip(mySpecies.y1, mySpecies.y22)])  # This just takes a differences and sums all the diffs

                if diff > 0:
                    bestStepSize = step * ((ebitParams[0].rkParams.desiredAccuracy / diff) ** 0.2)
                else:
                    time = time + step
                    mySpecies.population = copy.copy(mySpecies.y22)
                    continue  # if we got here then we can skip to the top of the while statement
                if bestStepSize >= step:
                    time = time + step
                    mySpecies.population = copy.copy(mySpecies.y22)
                else:
                    noTooBigSteps = noTooBigSteps + 1

                step = 0.9 * bestStepSize
        # print(mySpecies.results)
    return


def calcRateMatrices(mySpecies, myEbitParams):
        myEbitParams.currentDensity = (myEbitParams.ionEbeamOverlap * myEbitParams.beamCurrent) / (pi * (myEbitParams.beamRadius ** 2))
        mySpecies.ionizationRates = createDefaultInteractionRates(mySpecies, myEbitParams, createIonizationCrossSections)
        mySpecies.rrRates = createDefaultInteractionRates(mySpecies, myEbitParams, createRRCrossSections)
        mySpecies.chargeExchangeRates = createDefaultInteractionRates(mySpecies, myEbitParams, createChargeExchangeRates, 1)


def initEverything(species, ebitParams):
    global INITILIZEDEVERYTHING
    ebitParams[0].rkParams = RkStepParams()  # We're just going to shove this information into the first beam because why not?

    ebitParams[0].rkParams.desiredAccuracy = ebitParams[0].rkParams.desiredAccuracyPerChargeState / species[0].Z

    species.sort(key=lambda x: x.Z, reverse=True)  # Sort species list by Z in decending order
    largestZValue = species[0].Z  # We need the largest value for the size of the population array

    if species[0].betaHalfLife != 0.0:  # The largest Z value species is not allowed to have a beta decay
        raise ValueError("Can not handle having a beta decay for highest Z species")

    for mySpecies in species:
        # Initilize everything ...
        mySpecies.population = createEmptyList(largestZValue + 2)
        mySpecies.population[1] = mySpecies.initSCIPop
        mySpecies.decayConstant = createDecayConstants(mySpecies.betaHalfLife)

        if mySpecies.decayConstant > 0:  # If we don't need it we will use this to ignore the beta calculation function
            ebitParams.ignoreBetaDecay = 0

        ebitParams[0].decayConstants.append(mySpecies.decayConstant)

        mySpecies.k1 = createEmptyList(species[0].Z + 2)
        mySpecies.k2 = createEmptyList(species[0].Z + 2)
        mySpecies.k3 = createEmptyList(species[0].Z + 2)
        mySpecies.k4 = createEmptyList(species[0].Z + 2)
        mySpecies.tmpPop = createEmptyList(species[0].Z + 2)
        mySpecies.y1 = createEmptyList(species[0].Z + 2)
        mySpecies.y12 = createEmptyList(species[0].Z + 2)
        mySpecies.y22 = createEmptyList(species[0].Z + 2)
        mySpecies.results = []
    INITILIZEDEVERYTHING = 1


def calcChargePopulations(species, ebitParams, probeFnAddPop):
    if INITILIZEDEVERYTHING == 0:
        initEverything(species, ebitParams)

    adaptiveRkStepper(species, ebitParams, probeFnAddPop)  # this does all the heavy lifting

    return
