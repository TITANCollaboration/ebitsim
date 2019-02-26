#!/usr/bin/env python
from math import *
from IonizationEnergies import *
from rr import *
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
                 chargeExchangeRates=None):
        self.Z = Z
        self.A = A
        self.decaysTo = decaysTo
        self.betaHalfLife = betaHalfLife
        self.initSCIPop = initSCIPop
        self.population = population
        self.decayConstant = decayConstant
        self.ionizationRates = ionizationRates
        self.rrRates = rrRates
        self.chargeExchangeRates = chargeExchangeRates
        self.chargeStates = chargeStates


# I think I need to move chargeStates up to here.. could concievibly do that for a bunch
# of the stuff in chargepopparams due to maybe they should be following the species around

class RkStepParams:
    def __init__(self, minCharge=5e-5, tStep=1e-6, desiredAccuracyPerChargeState=1e-4):
        self.minCharge = minCharge
        self.tStep = tStep
        self.desiredAccuracyPerChargeState = desiredAccuracyPerChargeState


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
                 currentDensity=None,
                 rkParams=None,
                 results=[]):
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
        self.results = results


def createEmptyListofLists(species):
    emptylist = [[0.0 for i in range(species[0].Z + 2)] for j in range(len(species) + 1)]  # Init multidimentional array w/ 0's
    return emptylist


def createEmptyList(sizeOfArray):
    myEmptyList = [0.0 for i in range(sizeOfArray)]
    return myEmptyList


#def createDefaultPopulation(species):
# length of new array needs to be 1+# of species and 2+maxZ value
    #populationList = createEmptyListofLists(species)
   # for i in range(0, len(species)):
       # populationList[i][1] = species[i].initSCIPop
#    emptylist = [0.0 for i in range(species[0].Z + 2)]
#    return populationList


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


#def createDefaultInteractionRates(myspecies, beamEnergy, currentDensity, myfunc, pressure=0, ionTemperature=0, chex=0):
def createDefaultInteractionRates(myspecies, ebitparams, probeFnAddPop, runChargeExchange=0):  # move away from chex=0 and move into object

    interactionRates = createEmptyList(myspecies.Z + 2)
    if runChargeExchange == 0:
        myFuncValues = createInteractionRates(myspecies.Z, ebitparams.beamEnergy, ebitparams.currentDensity, probeFnAddPop(ebitparams.beamEnergy, myspecies.Z))
    else:
#        myFuncValues = myfunc(species[i].Z, species[i].A, pressure, ionTemperature)
        myFuncValues = probeFnAddPop(myspecies.Z, myspecies.A, ebitparams.pressure, ebitparams.ionTemperature)

    for r in range(0, len(myFuncValues)):
        interactionRates[r] = myFuncValues[r]
    return interactionRates


def betaDecay(s, q, tstep, ionizationRates, tmpPop):
    if q >= 1:
        myval = (-1 * decayConstants[s] * tmpPop[s][q]) + ( decayConstants[s + 1] * tmpPop[s + 1][q])
    else:
        myval = 0
    return tstep * myval


def chargeChanges(s, q, tstep, decayConstants, chargeExchangeRates, ionizationRates, rrRates, tmpPop):
    myval = tstep * ((-1 * ionizationRates[s][q] * tmpPop[s][q]) + ((chargeExchangeRates[s][q + 1] + rrRates[s][q + 1]) * tmpPop[s][q + 1]))
    return myval


def calculateK(chargePopParams, retK, p1, p2, addWFactor, tstep, tmpPop):
    print(betaDecay(0, 21, tstep, chargePopParams.ionizationRates, tmpPop))
    sys.exit(0)
    return


def rkStep(chargePopParams, tstep, populationT0, populationTtstep, k1, k2, k3, k4, tmpPop):
    calculateK(chargePopParams, k1, populationT0, tmpPop, 0.0, tstep, tmpPop)
    return

def probeFnAddPop(time, chargePopParams):
    #  In the original lisp code this function was passed along to calcChargePopulations
    #  so for now I will do the same as apparently this is something someone might want
    #  to change to a different function... If someone knows the reasonf for this please
    #  let me know.
    newresult = []
    print(chargePopParams.chargeStates)
    for charge in chargePopParams.chargeStates:
        print(chargePopParams.population)
        newresult.append([time, chargePopParams.population[charge - 1]])
    chargePopParams.result.append(newresult)

def adaptiveRkStepper(chargePopParams, probeFnAddPop):
    #  Must create some arrays...
    k1 = k2 = k3 = k4 = tmpPop = y1 = y12 = y22 = createEmptyListofLists(chargePopParams.species)
    time = 0.0
    nextPrint = 0.0
    noTooBigSteps = 0
    step = chargePopParams.rkParams.tStep  # Not sure if there is method to this madness yet..
    bestStepSize = step  # Currently just translating directly will try to be more nuanced after another pass

    desiredAccuracy = chargePopParams.rkParams.desiredAccuracyPerChargeState / chargePopParams.species[0].Z
    while time <= chargePopParams.breedingTime:
        if time >= nextPrint:
            nextPrint = nextPrint + chargePopParams.probeEvery
            print(chargePopParams.probeFn)
            probeFnAddPop(time, chargePopParams)
        rkStep(chargePopParams, 2 * step, chargePopParams.population, y1, k1, k2, k3, k4, tmpPop)
    return


def calcChargePopulations(species, ebitparams, probeFnAddPop):
    if ebitparams.currentDensity is None:
        ebitparams.currentDensity = (ebitparams.ionEbeamOverlap * ebitparams.beamCurrent) / (pi * (ebitparams.beamRadius ** 2))

    species.sort(key=lambda x: x.Z, reverse=True)  # Sort species list by Z in decending order
    largestZValue = species[0].Z  # We need the largest value for the size of the population array

    if species[0].betaHalfLife != 0.0:  # The largest Z value species is not allowed to have a beta decay
        raise ValueError("Can not handle having a beta decay for highest Z species")

    for myspecies in species:

        myspecies.population = createEmptyList(largestZValue + 2)
        myspecies.decayConstant = createDecayConstants(myspecies.betaHalfLife)

#        species.ionizationRates = createDefaultInteractionRates(chargePopParams.species, chargePopParams.beamEnergy, chargePopParams.currentDensity, createIonizationCrossSections)
        myspecies.ionizationRates = createDefaultInteractionRates(myspecies, ebitparams, createIonizationCrossSections)

#        species.rrRates = createDefaultInteractionRates(chargePopParams.species, chargePopParams.beamEnergy, chargePopParams.currentDensity, createRRCrossSections)
        myspecies.rrRates = createDefaultInteractionRates(myspecies, ebitparams, createRRCrossSections)
        #print(myspecies.rrRates)

#        species.chargeExchangeRates = createDefaultInteractionRates(chargePopParams.species, chargePopParams.beamEnergy, chargePopParams.currentDensity, createChargeExchangeRates, chargePopParams.pressure, chargePopParams.ionTemperature, 1)
        myspecies.chargeExchangeRates = createDefaultInteractionRates(myspecies, ebitparams, createChargeExchangeRates, 1)
        print(myspecies.chargeExchangeRates)

        ebitparams.rkParams = RkStepParams()

#    adaptiveRkStepper(chargePopParams, probeFnAddPop)
        adaptiveRkStepper(myspecies, ebitparams, probeFnAddPop)

    return


def main():
    chargeStates = list([39, 40, 41, 42, 43])

    species = []

    # rewrite a bunch of this....
    species.append(Species(51, 121, 0.0, 0.0, 1.0, chargeStates))
    species.append(Species(22, 122, 0.0, 0.01, 1.0, chargeStates))
    species.append(Species(54, 123, 0.0, 0.0, 1.0, chargeStates))


#    species.append(Species(2, 107, 0.0, 1.0, 6.0))
#    species.append(Species(6, 107, 0.0, 0.0, 5.0))
    ebitparams = EbitParams(breedingTime=10.0, beamEnergy=7000.0, pressure=1e-10, beamCurrent=0.02, beamRadius=90e-4)
    calcChargePopulations(species, ebitparams, probeFnAddPop)
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