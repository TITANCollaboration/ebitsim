#!/usr/bin/env python

import csv
from math import *

# this contains a bunch of misc functions for calculating ionization cross-sections in the EBIT.
# This is also directly translated from the somewhat original Lisp code.
__MINZ__ = 1
__MAXZ__ = 105


def readInCSV(filename):
    # Read in a CSV file and read the data as floats, returns a multidimentional array
    with open(filename, newline='') as f:
        csvarray = list(csv.reader(f, quoting=csv.QUOTE_NONNUMERIC))  # we want to make sure we are reading the data in as floats
    return csvarray


def reverseMyList(mylist):
    # Reverses the elements of a multidimentional list by the 2nd index and tacks on a 0 in front of each list
    # this is just to get it in the same format as the old Lisp code was in to make keeping track of index's easier
    # while translating.  It returns the new list of course
    newlist = [[0.0 for i in range(__MAXZ__ + 1)] for j in range(__MAXZ__ + 1)]  # Init multidimentional array w/ 0's

    for listIndex in range(0, len(mylist)):
        newlist[listIndex + 1] = [0.0] + list(reversed(mylist[listIndex][1:]))
    return newlist


def createIonizationEnergies(ionBindingEnergiesFromCSV):
    # Calculating the ionization energies from the binding energies and returns a multidimentional array
    ionizationEnergies = [[0.0 for i in range(__MAXZ__ + 1)] for j in range(__MAXZ__ + 1)]  # init multidimentional array w/ 0's

    for z in range(__MINZ__, __MAXZ__ + 1):
        for q in range(2, z + 1):

            # This round() was truncated to 2 decimal places. Not sure why, but I increased it and it didn't hurt anything.
            ionizationEnergies[z - 1][q - 1] = round((ionBindingEnergiesFromCSV[z - 1][z - q + 1] - ionBindingEnergiesFromCSV[z - 1][(z - q) + 2]), 4)

    return ionizationEnergies


def getCarlsonCorrectedBindEnergy(Z, q, e, ionizationEnergies, naturalBindingEnergies):
    # COMMENT STOLEN FROM OLD LISP CODE :"Given the proton number Z (1,105) and the charge state q (0,105),
    # returns the carlson corrected energy of a electron with number e for
    # the ion in charge state q. THIS NEXT PART IS NOT TESTED TO WORK : Be aware: If asked for a bogus charge
    # state, i.e. He10+, it will return 0d0."

    # I'm using this try/except so I don't have to pad the array with a bunch of
    # extra 0's
    try:
        mybindingenergy = -naturalBindingEnergies[Z][q]
    except IndexError:
        mybindingenergy = 0.0

    carlsonCorrection = naturalBindingEnergies[Z][e] + ionizationEnergies[Z][q] + mybindingenergy

    return carlsonCorrection


def lotzIonizationCrossSection(electronEnergy, Z, eToRemove, ionizationEnergies, naturalBindingEnergies):
    sigma = 0.0

    for j in range(eToRemove + 1, Z + 1):
        correctedBindingEnergy = getCarlsonCorrectedBindEnergy(Z, (eToRemove + 1), j, ionizationEnergies, naturalBindingEnergies)
        if (correctedBindingEnergy < electronEnergy) and (correctedBindingEnergy > 0):
            sigma = sigma + ((4.5e-14 * log(electronEnergy / correctedBindingEnergy)) / (electronEnergy * correctedBindingEnergy))
    return sigma

def createIonizationCrossSections(electronEnergy, Z):
    ionBindingEnergiesFromCSV = readInCSV('ionBindingEnergies.csv')  # convert from negatives to positives before maths
    naturalBindingEnergies = readInCSV('naturalBindingEnergiesList.csv')
    # Read function create-ionization-energies to get this...
    #ionizationEnergies = reverseMyList(ionizationEnergies)
    ionizationEnergies = createIonizationEnergies(ionBindingEnergiesFromCSV)
    naturalBindingEnergies = reverseMyList(naturalBindingEnergies)
    crossSections = [0] * (Z + 1)
    for q in range(0, Z + 1):
        crossSections[q] = lotzIonizationCrossSection(electronEnergy, Z, q, ionizationEnergies, naturalBindingEnergies)
    return crossSections
