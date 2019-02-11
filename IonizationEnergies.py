#!/usr/bin/env python

import csv
import sys

# work in progress...
# this contains a bunch of misc functions used for calculating ionization breeding
# in the ebit.
# This is also directly translated from the somewhat original Lisp code.
__MINZ__ = 1
__MAXZ__ = 105

def readInCSV(filename):
    with open(filename, newline='') as f:
        csvarray = list(csv.reader(f, quoting=csv.QUOTE_NONNUMERIC))
    return csvarray


def reverseMyList(mylist):
    #newlist = [0] * __MAXZ__
    newlist = [ [ 0.0 for i in range(__MAXZ__ + 1) ] for j in range(__MAXZ__ + 1) ]


    for listIndex in range(0, len(mylist)):
        newlist[listIndex + 1] = [0.0] + list(reversed(mylist[listIndex][1:]))
    return newlist


def createIonizationEnergies(ionBindingEnergiesFromCSV):
    ionizationEnergies = [ [ 0.0 for i in range(__MAXZ__ + 1) ] for j in range(__MAXZ__ + 1) ]

    for z in range(__MINZ__, __MAXZ__ + 1):
        for q in range(2, z + 1):

            ionizationEnergies[z - 1][q - 1] = round((ionBindingEnergiesFromCSV[z - 1][z - q + 1] - ionBindingEnergiesFromCSV[z - 1][(z - q) + 2]), 2)
    return ionizationEnergies

def getCarslonCorrectedBindEnergy(Z, q, e):
    # "Given the proton number Z (1,105) and the charge state q (0,105),
    # returns the carlson corrected energy of a electron with number e for
    # the ion in charge state q. Be aware: If asked for a bogus charge
    # state, i.e. He10+, it will return 0d0."
    ionBindingEnergiesFromCSV = readInCSV('ionBindingEnergies.csv') # convert from negatives to positives before maths
    naturalBindingEnergies = readInCSV('naturalBindingEnergiesList.csv')
    # Read function create-ionization-energies to get this...
    #ionizationEnergies = reverseMyList(ionizationEnergies)
    ionizationEnergies = createIonizationEnergies(ionBindingEnergiesFromCSV)
    naturalBindingEnergies = reverseMyList(naturalBindingEnergies)
    return
