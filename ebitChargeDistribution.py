#!/usr/bin/env python
from math import *
from IonizationEnergies import *
from rr import *


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


class Species:
    def __init__(self, Z, A, decaysTo=0.0, betaHalfLife=0.0, initSCIPop=1.0):
        self.Z = Z
        self.A = A
        self.decaysTo = decaysTo
        self.betaHalfLife = betaHalfLife
        self.initSCIPop = initSCIPop


class RkStepParams:
    def __init__(self, minCharge=5e-5, tStep=1e-6, desiredAccuracyPerChargeState=1e-4):
        self.minCharge = minCharge
        self.tStep = tStep
        self.desiredAccuracyPerChargeState = desiredAccuracyPerChargeState


def createEmptyListofLists(species):
    emptylist = [[0.0 for i in range(species[0].Z + 2)] for j in range(len(species) + 1)]  # Init multidimentional array w/ 0's
    return emptylist


def createDefaultPopulation(species):
# length of new array needs to be 1+# of species and 2+maxZ value
    populationList = createEmptyListofLists(species)
    for i in range(0, len(species)):
        populationList[i][1] = species[i].initSCIPop

    return populationList


def createDecayConstants(species):
    decayConstants = [0] * (len(species) + 1)

    if species[0].betaHalfLife != 0.0:
        raise ValueError("Can not handle having a beta decay for highest Z species")

    for i in range(0, len(species)):
        if species[i].betaHalfLife <= 0:
            decayConstants[i] = 0.0
        else:
            decayConstants[i] = log(2) / species[i].betaHalfLife
    return decayConstants

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


def createDefaultInteractionRates(species, beamEnergy, currentDensity, myfunc, pressure=0, ionTemperature=0, chex=0):
    interactionRates = createEmptyListofLists(species)
    for i in range(0, len(species)):

        if chex == 0:
            myFuncValues = createInteractionRates(species[i].Z, beamEnergy, currentDensity, myfunc(beamEnergy, species[i].Z))
        else:
            myFuncValues = myfunc(species[i].Z, species[i].A, pressure, ionTemperature)

        for r in range(0, len(myFuncValues)):
            interactionRates[i][r] = myFuncValues[r]

    return interactionRates


def calculateK():
    return


def rkStep():
    return


def adaptiveRkStepper(breedingTime, probeEvery, probeFn, population, species, decayConstants, ionizationRates, chexRates, rrRates, rkParams):
    #  Must create some arrays...
    k1 = k2 = k3 = k4 = tmp = y1 = y12 = y22 = createEmptyListofLists(species)
    time = 0.0
    nextPrint = 0.0
    noTooBigSteps = 0
    step = rkParams.tStep  # Not sure if there is method to this madness yet..
    bestStepSize = step  # Currently just translating directly will try to be more nuanced after another pass

    desiredAccuracy = rkParams.desiredAccuracyPerChargeState / species[0].Z

    return


def calcChargePopulations(species,
                          probeFn='something',
                          breedingTime=0.1,
                          probeEvery=0.1,
                          ionEbeamOverlap=1.0,
                          beamEnergy=3000.0,
                          beamCurrent=0.1,
                          beamRadius=200.0e-4,
                          pressure=1e-12,
                          ionTemperature=100.0,
                          currentDensity=None,
                          population=None,
                          decayConstants=None,
                          ionizationRates=None,
                          rrRates=None,
                          chexRates=None,
                          rkParams=None):

    if currentDensity is None:
        currentDensity = (ionEbeamOverlap * beamCurrent) / (pi * (beamRadius ** 2))

    species = []
    species.append(Species(4, 107, 0.0, 0.0, 7.0))
    species.append(Species(2, 107, 0.0, 1.0, 6.0))
    species.append(Species(6, 107, 0.0, 0.0, 5.0))
    species.sort(key=lambda x: x.Z, reverse=True)  # Sort by Z in decending order

    if population is None:
        population = createDefaultPopulation(species)

    if decayConstants is None:
        decayConstants = createDecayConstants(species)

    if ionizationRates is None:
        ionizationRates = createDefaultInteractionRates(species, beamEnergy, currentDensity, createIonizationCrossSections)

    if rrRates is None:
        rrRates = createDefaultInteractionRates(species, beamEnergy, currentDensity, createRRCrossSections)

    if chexRates is None:
        chexRates = createDefaultInteractionRates(species, beamEnergy, currentDensity, createChargeExchangeRates, pressure, ionTemperature, 1)

    rkParams = RkStepParams()

    adaptiveRkStepper(breedingTime, probeEvery, probeFn, population, species, decayConstants, ionizationRates, chexRates, rrRates, rkParams)