#!/usr/bin/env python

import csv
#work in progress...
# this contains a bunch of misc functions used for calculating ionization breeding
# in the ebit.
# This is also directly translated from the somewhat original Lisp code.

def getCarslonCorrectedBindEnergy(Z, q, e):
    # "Given the proton number Z (1,105) and the charge state q (0,105),
    # returns the carlson corrected energy of a electron with number e for
    # the ion in charge state q. Be aware: If asked for a bogus charge
    # state, i.e. He10+, it will return 0d0."

    return
