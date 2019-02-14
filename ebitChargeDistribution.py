#!/usr/bin/env python
from math import *

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
    Z = 0
    A = 0
    decaysTo = 0.0
    betaHalfLife = 0.0
    initSCIPop = 1.0

    def __init__(self, Z, A, decaysTo=0.0, betaHalfLife=0.0, initSCIPop=1.0):
        self.Z = Z
        self.A = A
        self.decaysTo = decaysTo
        self.betaHalfLife = betaHalfLife
        self.initSCIPop = initSCIPop


def createEmptyListofLists(species):
    emptylist = [[0.0 for i in range(species[0].Z + 2)] for j in range(len(species) + 1)]  # Init multidimentional array w/ 0's
    return emptylist


def createDefaultPopulation(species):
    species.sort(key=lambda x: x.Z, reverse=True)  # Sort by Z in decending order
# length of new array needs to be 1+# of species and 2+maxZ value
    populationList = createEmptyListofLists(species)
    for i in range(0, len(species)):
        populationList[i][1] = species[i].initSCIPop

    return populationList


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
                          decayConsants=None,
                          ionizationRates=None,
                          rrRates=None,
                          chexRates=None,
                          rkParams=None):

    if currentDensity is None:
        currentDensity = (ionEbeamOverlap * beamCurrent) / (pi * (beamRadius ** 2))
    mine = []

    mine.append(Species(4, 107, 0.0, 0.0, 7.0))
    mine.append(Species(2, 107, 0.0, 0.0, 6.0))
    mine.append(Species(6, 107, 0.0, 0.0, 5.0))
    createDefaultPopulation(mine)

#    if population is None:
#        population = createDefaultPopulation(species)
