#!/usr/bin/env python

from math import *
# These are set to the values they are because that is how the lisp
# code I translated this from did it...
__GSIZE__ = 101
__NMAX__ = 30


def g0():
    # I have little to no idea what this does nor did the person
    # who wrote the original lisp code I translated this from...
    # it's definitly a matrix though!
    # If you know what this does please let me know and I will owe
    # you a beer.
    global __GSIZE__
    global __NMAX__
    g0_list = [0.0] * __GSIZE__  # init list to all 0.0's
    g0_list[0] = sqrt(8 * pi) * exp(log(4) - 2)

    for n in range(1, __NMAX__):
        g0_list[n] = (g0_list[n-1] * 4 * n * exp((n + 1) * log((n + 1) / n))) \
            / sqrt((2 * (n + 1) - 1)) \
            / sqrt(2 * (n + 1) - 2) \
            / exp(2)
    return g0_list


def determineShell(q, z):

    nmin = 0
    lstart = 0
    num_electrons = z - q

    if q == z:
        nmin = 1
    else:
        while num_electrons > 0:
            num_electrons = num_electrons - 2
            nmin = nmin + 1

            if num_electrons < 0:
                break
            for l in range(1, nmin):
                num_electrons = num_electrons - ((4 * l) + 2)
                if num_electrons < 0:
                    lstart = l
                    break
    if num_electrons != 0:
        lastval = abs(num_electrons) / (((lstart * 2) + 1) * 2)
    else:
        lastval = 1
    return nmin, lstart, lastval


def createGGArray():
    # Create multidimention list that is 101 x 101... yup..
    gg = [[0] * __GSIZE__ for i in range(__GSIZE__)]
    return gg


def createGG(n, kappa, g0, gg):
    # Starts populating the gg matrix with stuff?
    f1 = 1 + (n * kappa) ** 2
    f = ((sqrt(1 / (1 - exp(-(2.0 * pi) / kappa))) * exp((2.0 * n) - (2.0 / kappa) * atan(n * kappa)) / f1 ** 2) * g0[n - 1])

    for s in range(1, n + 1):
        f = (f * sqrt(1+ (s * kappa) ** 2)) / f1
    gg[n - 1][n] = f

    if (n - 2) < 0:
        return gg

    gg[n - 2][n - 1] = 0.5 * sqrt(((2.0 * n) - 1) * f1) * f
    gg[n - 1][n - 2] = (0.5 / n) * sqrt(f1 / (1 + ((n - 1) * kappa) ** 2)) * f

    if (n - 3) < 0:
        return gg
    gg[n - 2][n - 3] = ((((n - 1) * f1) + 4) / 2 / n) * sqrt(((2 * n) - 1) / (1 + ((n - 2) * kappa) ** 2)) * gg[n - 1][n - 2]
    for l in range(n - 1, 1, -1):
        a = (((n ** 2) - (l ** 2)) * 4) + (((2 * l) - 1) * l * f1)
        b = -2 * n * sqrt((1 + ((1 + l) * kappa) ** 2) * ((n ** 2) - (l ** 2)))
        c = 2 * n * sqrt((1 + ((l * kappa) ** 2)) * ((n ** 2) - ((l - 1) ** 2)))
        d = ((4 * (n ** 2) ) - (4 * (l ** 2))) + (l * (1 + (2 * l)) * f1)
        e = -2 * n * sqrt(((n ** 2) - ((l + 1) ** 2)) * (1 + ((l * kappa) ** 2)))
        f = 2 * n * sqrt(((n ** 2) - (l ** 2)) * (1 + (((l - 1) * kappa) ** 2)))

        gg[l - 2][l - 1] = ((a * gg[l - 1][l]) + (b * gg[l][l + 1])) / c

        if l < (n - 1):
            gg[l - 1][l - 2] = ((d * gg[l][l - 1]) + (e * gg[l + 1][l])) / f

    return gg

def rrCrossSection(eKinetic, q, Z, gg=createGGArray(), g0=g0()):
    """ This looks to be a time reversed photoionization. So we determine the photonionization cross section and
    then multiply by some factor ahead to get the reversed cross section for RR.
    """
    ry = 13.605698
    alpha = 1 / 137.036
    a0 = 5.29177249e-9
    const = (4 * pi * alpha * (a0 ** 2)) / 3
    k = sqrt(eKinetic / ry)
    eta = (q + Z) / (2 * k)
    cross = 0.0

    nmin, lstart, factor = determineShell(q, Z)

    for n in range(nmin, __NMAX__ + 1):
        crossShell = 0
        if n == nmin:
            loop = lstart
        else:
            loop = 0
        for l in range(loop, n):
            if q == Z:
                zeff = Z
            else:
                zeff = (((q + Z) * 0.5) - (((Z - q) * 0.5 * (eta - 1)) / (eta + 1 + (3 * l)))) * exp(-(0.05 * (l - 1) ** 2))

            kappa = k / zeff
            gg = createGG(n, kappa, g0, gg)  # **Still need to figure out if this should be set to gg or a new var
            kk = (k * k) + ((zeff ** 2) / (n ** 2))

            tmpCross = (l + 1) * (gg[l][l + 1] ** 2) * ((1 + (n ** 2) * (kappa ** 2)))
            if l > 0:
                tmpCross = tmpCross +  ((l) * (gg[l][l - 1] ** 2) * ((1 + (n ** 2) * (kappa ** 2))))
            tmpCross = (tmpCross * n * n * const) / zeff / zeff
            crossSubShell = 0.5 * tmpCross * (((alpha * kk) / k) ** 2)
            if (n == nmin) and (l == lstart):
                crossSubShell = crossSubShell * factor

            if (l > loop) and (crossSubShell < 1e-50):
                break
            else:
                crossShell = crossShell + crossSubShell
            cross = cross + crossSubShell
    return(cross)

def createRRCrossSections(ElectronEnergy, Zion):
    crossSections = [0] * (Zion + 1)
    for q in range(1, Zion + 1):
        crossSections[q] = rrCrossSection(ElectronEnergy, q, Zion)

    return(crossSections)