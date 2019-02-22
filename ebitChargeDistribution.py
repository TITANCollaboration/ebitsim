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


class ChargePopulationParams:
    def __init__(self,
                 species,
                 chargeStates,
                 probeFn,
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
                 rkParams=None,
                 results=[]):
        self.species = species
        self.chargeStates = chargeStates,
        self.probeFn = probeFn,
        self.breedingTime = breedingTime
        self.probeEvery = probeEvery
        self.ionEbeamOverlap = ionEbeamOverlap
        self.beamEnergy = beamEnergy
        self.beamCurrent = beamCurrent
        self.beamRadius = beamRadius
        self.pressure = pressure
        self.ionTemperature = ionTemperature
        self.currentDensity = currentDensity
        self.population = population
        self.decayConstants = decayConstants
        self.ionizationRates = ionizationRates
        self.rrRates = rrRates
        self.chexRates = chexRates
        self.rkParams = rkParams
        self.results = results

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


def calculateK(chargePopParams, retK, p1, p2, addWFactor):
    return


def rkStep(chargePopParams, tstep, populationT0, populationTtstep, k1, k2, k3, k4, tmp):

    return


def adaptiveRkStepper(chargePopParams):
    #  Must create some arrays...
    k1 = k2 = k3 = k4 = tmp = y1 = y12 = y22 = createEmptyListofLists(chargePopParams.species)
    time = 0.0
    nextPrint = 0.0
    noTooBigSteps = 0
    step = chargePopParams.rkParams.tStep  # Not sure if there is method to this madness yet..
    bestStepSize = step  # Currently just translating directly will try to be more nuanced after another pass

    desiredAccuracy = chargePopParams.rkParams.desiredAccuracyPerChargeState / chargePopParams.species[0].Z
    print(desiredAccuracy)
    while time <= chargePopParams.breedingTime:
        if time >= nextPrint:
            nextPrint = nextPrint + chargePopParams.probeEvery
            chargePopParams.probeFn(time, chargePopParams)
        rkStep(chargePopParams, 2 * step, chargePopParams.population, y1, k1, k2, k3, k4, tmp)
    return


def probeFnAddPop(time, chargePopParams):
    #  In the original lisp code this function was passed along to calcChargePopulations
    #  so for now I will do the same as apparently this is something someone might want
    #  to change to a different function... If someone knows the reasonf for this please
    #  let me know.
    newresult = []
    for charge in chargePopParams.chargeStates:
        newresult.append([time, chargePopParams.population[charge - 1]])
    chargePopParams.result.append(newresult)


def calcChargePopulations(chargePopParams):

    if chargePopParams.currentDensity is None:
        chargePopParams.currentDensity = (chargePopParams.ionEbeamOverlap * chargePopParams.beamCurrent) / (pi * (chargePopParams.beamRadius ** 2))

    chargePopParams.species.sort(key=lambda x: x.Z, reverse=True)  # Sort by Z in decending order

    if chargePopParams.population is None:
        chargePopParams.population = createDefaultPopulation(chargePopParams.species)

    if chargePopParams.decayConstants is None:
        chargePopParams.decayConstants = createDecayConstants(chargePopParams.species)

    if chargePopParams.ionizationRates is None:
        chargePopParams.ionizationRates = createDefaultInteractionRates(chargePopParams.species, chargePopParams.beamEnergy, chargePopParams.currentDensity, createIonizationCrossSections)

    if chargePopParams.rrRates is None:
        chargePopParams.rrRates = createDefaultInteractionRates(chargePopParams.species, chargePopParams.beamEnergy, chargePopParams.currentDensity, createRRCrossSections)

    if chargePopParams.chexRates is None:
        chargePopParams.chexRates = createDefaultInteractionRates(chargePopParams.species, chargePopParams.beamEnergy, chargePopParams.currentDensity, createChargeExchangeRates, chargePopParams.pressure, chargePopParams.ionTemperature, 1)

    chargePopParams.rkParams = RkStepParams()

    adaptiveRkStepper(chargePopParams)

    return


def main():
    chargeStates = [39, 40, 41, 42, 43]
    species = []
    species.append(Species(4, 107, 0.0, 0.0, 7.0))
    species.append(Species(2, 107, 0.0, 1.0, 6.0))
    species.append(Species(6, 107, 0.0, 0.0, 5.0))
    myparams = ChargePopulationParams(species, chargeStates, probeFnAddPop)
    calcChargePopulations(myparams)
    return


def testpopdata():
    result = [1.0773715045951986e-10,-2.5221903096086526e-10,2.3449007288242245e-10, \
        -1.2221336132501444e-10,4.026919164871968e-11,-9.343245253069252e-12, \
        1.4314856326808133e-12,-1.6629498403738547e-13,1.5050512947015083e-14, \
        -1.0786148204346698e-15,6.280613072141581e-17,-2.9130702329533736e-18, \
        1.1194883347241496e-19,-3.586879680201691e-21,9.650587100702459e-23, \
        -2.1886639635205885e-24,5.455806595921941e-26,1.684653698776325e-26, \
        1.3161935364266835e-26,5.531656225764603e-27,3.2398603494337613e-28, \
        -1.0338100473059566e-27,-2.8913558824647924e-27,-8.991390624069855e-27, \
        -1.2091860673384666e-26,2.0231038762218803e-26,7.085002901940573e-24, \
        1.2688700954569216e-21,1.8624840081805598e-19,2.210926329822218e-17, \
        2.116430358816347e-15,1.6355013404941251e-13,1.0084173983337306e-11, \
        4.937354709682791e-10,2.039748431910559e-8,6.670812700176424e-7, \
        1.6861321657554175e-5,3.2227812666130867e-4,0.004264150985831597e0, \
        0.03850582390819178e0,0.22895681488195369e0,0.5902408926186803e0, \
        0.1329114487663367e0,0.004781041407944411e0,0.0e0,0.0e0,0.0e0,0.0e0,0.0e0, \
        0.0e0,0.0e0,0.0e0,0.0e0]
    return result