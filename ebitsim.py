#!/usr/bin/env python
# Jon Ringuette - 03/2019

# I know in python one is supposed to use snake_case instead of camelCase but I kind of liked camelCase for this
# program better and I stuck with it throughout...  If someone is disturbed by this i'm happy to change it.
# other than that I have kept everything PEP8 complient for the most part except line length because we
# all have widescreen monitors these days..


import argparse
import configparser
import csv
import ebitChargeDistribution
from ebitChargeDistribution import probeFnAddPop

import platform
import sys
import os

# I suspect this class will grow with some additional options as they become needed
class OutputFormating:
    def __init__(self, outputFileName='output.png', outputType='matplotlib'):
        self.outputFileName = outputFileName
        self.outputType = outputType


def column(matrix, index):
    # Simple function to be able to pull out a single column of data from a list of lists, this lets us plot easily
    # as we can extract out our x,y values
    mycolumn = [i[index] for i in matrix]
    return mycolumn


def getElementAbv(z):
    with open('PeriodicTable.csv', newline='') as csvfile:
        elements = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in elements:
            if row[0] == str(z):
                abv = row[2]
    return abv


def plotSpeciesResults(species, ebitparams, outputConfig):
    #  Plot all species and charge state populations vs. time via matplotlib

    import matplotlib
    matplotlib.use("TKAgg")  # This is done as otherwise it goes wierd with a MacOS
    import matplotlib.pyplot as plt
    from cycler import cycler
    # Should figure out how to put tick marks on the top and right side of graphs, they're useful!
    plt.rcParams['legend.loc'] = 'best'
    plt.figure()

    for myspecies in species:
        for chargeStateResults in range(0, len(myspecies.results)):
            mylabel = getElementAbv(myspecies.Z) + str(myspecies.chargeStates[chargeStateResults]) + '+'
            default_cycler = (cycler(color=['k', 'g', 'b', 'y', 'c', 'm', 'k']) + cycler(linestyle=['-', '--', ':', '-.', '--', '-', ':']))
            plt.rc('lines', linewidth=2)
            plt.rc('grid', color='k', linestyle=':', linewidth=0.5)
            plt.rc('axes', grid=True, prop_cycle=default_cycler)
            plt.plot(column(myspecies.results[chargeStateResults], 0), column(myspecies.results[chargeStateResults], 1), label=mylabel)

    plt.ylabel('Population (%)')
    plt.xlabel('Breeding time (s)')
    plt.title("$I_e$ = %.2fA, $V_e$ = %ieV, $r_e$ = %.2Ecm, $P_H$ = %.2ET" % (ebitparams.beamCurrent, ebitparams.beamEnergy, ebitparams.beamRadius, ebitparams.pressure))
    plt.legend(framealpha=0.5)
    plt.savefig(outputConfig.outputFileName, dpi=300)

    return


def runSimulation(species, ebitparams, probeFnAddPop, outputConfig):
    #  Runs the actual simulation via ebitChargeDistribution.calcChargePopulations and then determines how to handle the output
    print("Running simulation! ....")
    ebitChargeDistribution.calcChargePopulations(species, ebitparams, probeFnAddPop)

    if outputConfig.outputType == 'matplotlib':
        print("Writing graph to : %s" % outputConfig.outputFileName)
        plotSpeciesResults(species, ebitparams, outputConfig)
    return


def getConfigEntry(config, heading, item, reqd=False, remove_spaces=True, default_val=''):
    #  Just a helper function to process config file lines, strip out white spaces and check if requred etc.
    if config.has_option(heading, item):
        if remove_spaces:
            config_item = config.get(heading, item).replace(" ", "")
        else:
            config_item = config.get(heading, item)
    elif reqd:
        print("The required config file setting \'%s\' under [%s] is missing") % (item, heading)
        sys.exit(1)
    else:
        config_item = default_val
    return config_item


def processConfigFile(configFileName):
    #
    # We will read the entire config file here and push it into a different classes and then call runSimulation()
    #
    config = configparser.RawConfigParser()
    baseDir = os.path.dirname(os.path.realpath(__file__))[:-3]

    if os.path.exists(configFileName):
        config.read(configFileName)
        ebitparams = ebitChargeDistribution.EbitParams()
        outputConfig = OutputFormating()

        configSections = config.sections()

        # Read output information
        outputConfig.outputType = getConfigEntry(config, 'Output', 'outputType', reqd=True, remove_spaces=True)
        outputConfig.outputFileName = getConfigEntry(config, 'Output', 'outputFileName', reqd=True, remove_spaces=True)

        # Read in Beam information
        ebitparams.beamEnergy = float(getConfigEntry(config, 'BeamAndTrap', 'beamEnergy', reqd=True, remove_spaces=True))
        ebitparams.breedingTime = float(getConfigEntry(config, 'BeamAndTrap', 'breedingTime', reqd=True, remove_spaces=True))
        ebitparams.probeEvery = float(getConfigEntry(config, 'BeamAndTrap', 'probeEvery', reqd=True, remove_spaces=True))
        ebitparams.ionEbeamOverlap = float(getConfigEntry(config, 'BeamAndTrap', 'ionEbeamOverlap', reqd=True, remove_spaces=True))
        ebitparams.beamCurrent = float(getConfigEntry(config, 'BeamAndTrap', 'beamCurrent', reqd=True, remove_spaces=True))
        ebitparams.beamRadius = float(getConfigEntry(config, 'BeamAndTrap', 'beamRadius', reqd=True, remove_spaces=True))
        ebitparams.pressure = float(getConfigEntry(config, 'BeamAndTrap', 'pressure', reqd=True, remove_spaces=True))
        ebitparams.ionTemperature = float(getConfigEntry(config, 'BeamAndTrap', 'ionTemperature', reqd=True, remove_spaces=True))

        # Read in Run section
        speciesList = tuple(getConfigEntry(config, 'Run', 'speciesList', reqd=True, remove_spaces=True).split(","))
        species = []

        for myspecies in speciesList:      # Collect all the species and parameters for each
            protons = int(getConfigEntry(config, myspecies, 'z', reqd=True, remove_spaces=True))
            nucleons = int(getConfigEntry(config, myspecies, 'nucleons', reqd=True, remove_spaces=True))
            population = float(getConfigEntry(config, myspecies, 'populationPercent', reqd=True, remove_spaces=True))
            chargeStates = list(map(int, getConfigEntry(config, myspecies, 'chargeStates', reqd=True, remove_spaces=True).split(",")))
            species.append(ebitChargeDistribution.Species(protons, nucleons, 0.0, 0.0, population, chargeStates))
    else:
        print("Config file does not appear to exist : %s" % configFileName)
        sys.exit(1)

#    ebitChargeDistribution.calcChargePopulations(species, ebitparams, probeFnAddPop)
#    plotSpeciesResults(species, ebitparams, outputFileName)
    runSimulation(species, ebitparams, probeFnAddPop, outputConfig)
    return 0


def processCommandLine(args):
    #  Process command line arguments and then call runSimulation()
    outputConfig = OutputFormating()
    species = []
    outputConfig.outputType = args.outputType
    outputConfig.outputFileName = args.outputFileName

    species.append(ebitChargeDistribution.Species(args.protons, args.nucleons, 0.0, 0.0, 1.0, args.chargeStates))
    ebitparams = ebitChargeDistribution.EbitParams(breedingTime=args.breedingTime, beamEnergy=args.beamEnergy, pressure=args.pressure, beamCurrent=args.beamCurrent, beamRadius=args.beamRadius, probeEvery=args.probeEvery)

    runSimulation(species, ebitparams, probeFnAddPop, outputConfig)

    return


def main():

    parser = argparse.ArgumentParser(description='EBIT Charge Breeding Simulation')

    parser.add_argument('--configFile', dest='configFile', required=False,
                        help="Specify the complete path to the config file, by default we'll use ebitsim.cfg")
    parser.add_argument('--outputType', dest='outputType', default='matplotlib', required=False,
                        help="Specify how to output the data, defaults to a matplotlib graph, csv will be available at some point")
    parser.add_argument('--outputFileName', dest='outputFileName', default='output.png', required=False,
                        help="Specify filename for csv or png file, please include .csv or .png extension")
    parser.add_argument('-z', dest='protons', type=int, default=0, required=False,
                        help="Specify the num of protons")
    parser.add_argument('-Z', dest='protons', type=int, default=0, required=False,  # there's probably a better way to handle this... sorry
                        help="Specify the num of protons")
    parser.add_argument('-a', dest='nucleons', type=int, required=False,
                        help="Specify the num of nucleons")
    parser.add_argument('--chargeStates', nargs='+', dest='chargeStates', type=int, required=False,
                        help="Specify charge states such as : --charge-states 32 33 34 35")
    parser.add_argument('--breedingTime', dest='breedingTime', default=11.0, type=float, required=False,
                        help="Specify breeding time (sec)")
    parser.add_argument('--probeEvery', dest='probeEvery', default=0.001, type=float, required=False, # You can probably crank this up if you want to or need to
                        help="Get results in this incriment (sec)")
    parser.add_argument('--beamEnergy', dest='beamEnergy', default=7000, type=float, required=False,
                        help="Beam energy (eV)")
    parser.add_argument('--beamCurrent', dest='beamCurrent', default=0.02, type=float, required=False,
                        help="Beam current (Amps)")
    parser.add_argument('--beamRadius', dest='beamRadius', default=90e-4, type=float, required=False,
                        help="Beam radius (cm)")
    parser.add_argument('--pressure', dest='pressure', default=1e-10, type=float, required=False,
                        help="ebit vacuum pressure (Torr / cm^3 I think)")
    parser.set_defaults(configFile="ebitsim.cfg")

    args, unknown = parser.parse_known_args()
    # sg.config_file = args.config_file

    if platform.python_implementation() == 'CPython':
        print("*** !!WARNING!! : While this can run via normal cPython it is highly recommended that you run it via pypy3 for a HUGE speedup ***")

    if args.protons != 0:
        processCommandLine(args)
    else:
        processConfigFile(args.configFile)

if __name__ == "__main__":
    main()
