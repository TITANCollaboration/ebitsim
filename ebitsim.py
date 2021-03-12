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
from ebitsim_docs import *
from ebitChargeDistribution import probeFnAddPop
from geant4MacroOutput import geant4MacroOutput
from commonUtils import getElementAbv, column

import platform
import sys
import os
import time

# I suspect this class will grow with some additional options as they become needed
class OutputFormating:
    def __init__(self, outputFileName='output.png',
                       outputType='matplotlib',
                       logOutput=0,
                       stepPlot = 0,
                       badGuessPlot = 0,
                       xmin=0,
                       xmax=0,
                       ymin=0,
                       ymax=0,
                       logx=0,
                       eventsPerTimeSlice=0,
                       subDivisionOfTime=0):
        self.outputFileName = outputFileName
        self.outputType = outputType
        self.logOutput = logOutput
        self.stepPlot = stepPlot
        self.badGuessPlot = badGuessPlot

        # matPlotLib stuff
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.logx = logx

        #geant4MacroOutput stuff
        self.eventsPerTimeSlice = eventsPerTimeSlice
        self.subDivisionOfTime = subDivisionOfTime


def plotSpeciesResults(species, ebitParams, outputConfig):
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
            default_cycler = (cycler(color=['b', 'k', 'g', 'y', 'c', 'm', 'r']) + cycler(linestyle=['--', '-', ':', '-.', '--', '-', ':']))
            plt.rc('lines', linewidth=2)
            plt.rc('grid', color='k', linestyle=':', linewidth=0.5)
            plt.rc('axes', grid=True, prop_cycle=default_cycler)
            plt.plot(column(myspecies.results[chargeStateResults], 0), column(myspecies.results[chargeStateResults], 1), label=mylabel)

            if outputConfig.xmin or outputConfig.xmax != 0:
                plt.xlim(outputConfig.xmin, outputConfig.xmax)
            if outputConfig.ymin or outputConfig.ymax != 0:
                plt.ylim(outputConfig.ymin, outputConfig.ymax)
            if outputConfig.logx == 1:
                plt.xscale('log')

    plt.ylabel('Population')
    plt.xlabel('Breeding time (s)')
    beamEnergies = ''

    for ebitIndex in range(0, len(ebitParams)):
        if ebitIndex != 0:
            beamEnergies = beamEnergies + '->'
        beamEnergies = beamEnergies + str(ebitParams[ebitIndex].beamEnergy)

    plt.title("$I_e$ = %.2fA, $V_e$ = (%s)eV, $r_e$ = %.2Ecm, $P_H$ = %.2ET" % (ebitParams[0].beamCurrent, beamEnergies, ebitParams[0].beamRadius, ebitParams[0].pressure))
    plt.legend(framealpha=0.5)
    plt.savefig(outputConfig.outputFileName, dpi=300)

    return

def plotSpeciesEnergies(species, ebitParams, outputConfig):
    import matplotlib
    # matplotlib.use("TKAgg")
    import matplotlib.pyplot as plt
    from cycler import cycler


    plt.rcParams['legend.loc'] = 'best'
    plt.figure()
    for mySpecies in species:
        for qResults in range(0, len(mySpecies.results)):
            mylabel = getElementAbv(mySpecies.Z) + str(mySpecies.chargeStates[qResults])

            plt.plot(column(mySpecies.results[qResults], 0), column(mySpecies.results[qResults], 2), label=mylabel)

            if outputConfig.logx == 1:
                plt.xscale('log')
    plt.ylabel("Energy [eV]")
    plt.xlabel("Breeding time [s]")
    plt.legend(framealpha=0.5)
    plt.savefig(outputConfig.outputFileName.split(".")[0]+"_energy.png", dpi=300)
    return


def writeCSVFile(species, ebitParams, outputConfig):
    newentry = []

    with open(outputConfig.outputFileName, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_NONE)

        for myspecies in species:
            for chargeStateResults in range(0, len(myspecies.results)):
                mylabel = getElementAbv(myspecies.Z) + str(myspecies.chargeStates[chargeStateResults])
                for myrow in range(len(myspecies.results[chargeStateResults])):
                    newentry = mylabel, myspecies.results[chargeStateResults][myrow][0], myspecies.results[chargeStateResults][myrow][1]
                    csvwriter.writerow(newentry)
    return


def writeRates(species, ebitParams, outputConfig):
    # An option to write out the rate matricies for diagnostic purposes.

    with open(outputConfig.outputFileName, 'w', newline='') as ratesfile:
        csvwriter = csv.writer(ratesfile, delimiter=',', quoting=csv.QUOTE_NONE)

        # Only showing one of the EBIT parameters, no beam energy scanning!
        csvwriter.writerow(['Breeding time: %s' % ebitParams[0].breedingTime])
        csvwriter.writerow(['Beam energy: %s' % ebitParams[0].beamEnergy])

        for myspecies in species:
            csvwriter.writerow(['Species = %s' % getElementAbv(myspecies.Z)])
            csvwriter.writerow(['ionization rates for q=0 to %s:' % str(len(myspecies.ionizationRates)-1)])
            for i, j in enumerate(myspecies.ionizationRates):
                csvwriter.writerow(['q = %s' %i] + [j])

            csvwriter.writerow(['radiative recombination rates:'])
            for i, j in enumerate(myspecies.rrRates):
                csvwriter.writerow(['q = %s' %i] + [j])

            csvwriter.writerow(['charge exchange rates:'])
            for i, j in enumerate(myspecies.chargeExchangeRates):
                csvwriter.writerow(['q = %s' %i] + [j])

        csvwriter.writerow(['END'])


def runSimulation(species, ebitParams, probeFnAddPop, outputConfig):
    #  Runs the actual simulation via ebitChargeDistribution.calcChargePopulations and then determines how to handle the output

    # Option to write to a log.log file in the same directory as the outputFileName
    if outputConfig.logOutput in ["True", "T", "Yes", "Y", 1]:
        print("Writing output to log file...")
        old_stdout = sys.stdout
        log_file = open(os.path.dirname(outputConfig.outputFileName)+os.sep+"log.log", "w")
        sys.stdout = log_file
        print(time.asctime(time.localtime()))

    print("Running simulation! ....")
    ebitChargeDistribution.calcChargePopulations(species, ebitParams, probeFnAddPop)

    if outputConfig.outputType == 'rates':
        print("Writing rates to csv: %s \n" % outputConfig.outputFileName)
        writeRates(species, ebitParams, outputConfig)
    if outputConfig.outputType == 'matplotlib':
        print("Writing charge state graph to : %s" % outputConfig.outputFileName)
        plotSpeciesResults(species, ebitParams, outputConfig)  # Think about fixing ebitParams later to deal with multiple beam energies..
        if species[0].initSCITemp != None: # if including energy dynamics, plot that too
            print("Writing energy graph to : %s" % (outputConfig.outputFileName.split(".")[0]+"_energy.png"))
            plotSpeciesEnergies(species, ebitParams, outputConfig)
    if outputConfig.outputType == 'csv':
        print("Writing csv to : %s \n" % outputConfig.outputFileName)
        writeCSVFile(species, ebitParams[0], outputConfig)
    if outputConfig.outputType == 'geant4Macro':
        print("Writing GEANT4 Macro to : %s \n" % outputConfig.outputFileName)
        geant4MacroOutput(species, ebitParams[0], outputConfig)

    if outputConfig.stepPlot != 0:
        print("Writing step size plot to %s \n" % os.path.dirname(outputConfig.outputFileName)+os.sep+outputConfig.stepPlot)
        plotStepSize(species, ebitParams, outputConfig)

    if outputConfig.logOutput in ["True", "T", "Yes", "Y", 1]:
        sys.stdout = old_stdout
        log_file.close
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
    # baseDir = os.path.dirname(os.path.realpath(__file__))[:-3]

    if os.path.exists(configFileName):
        config.read(configFileName)
        outputConfig = OutputFormating()

        # Read output information
        outputConfig.outputType = getConfigEntry(config, 'Output', 'outputType', reqd=True, remove_spaces=True)
        outputConfig.outputFileName = getConfigEntry(config, 'Output', 'outputFileName', reqd=True, remove_spaces=True)
        outputConfig.logOutput = getConfigEntry(config, 'Output', 'logOutput', reqd=False, remove_spaces=True, default_val=0)
        outputConfig.stepPlot = getConfigEntry(config, 'Output', 'stepPlot', reqd=False, remove_spaces=True, default_val=0)
        outputConfig.badGuessPlot = getConfigEntry(config, 'Output', 'badGuessPlot', reqd=False, remove_spaces=True, default_val=0)

        # matPlotLib stuff
        outputConfig.xmin = float(getConfigEntry(config, 'matPlotLib', 'graphXMinTime', reqd=False, remove_spaces=True, default_val=0))
        outputConfig.xmax = float(getConfigEntry(config, 'matPlotLib', 'graphXMaxTime', reqd=False, remove_spaces=True, default_val=0))
        outputConfig.ymin = float(getConfigEntry(config, 'matPlotLib', 'graphYMinPop', reqd=False, remove_spaces=True, default_val=0))
        outputConfig.ymax = float(getConfigEntry(config, 'matPlotLib', 'graphYMaxPop', reqd=False, remove_spaces=True, default_val=0))
        outputConfig.logx = float(getConfigEntry(config, 'matPlotLib', 'graphXScale', reqd=False, remove_spaces=True, default_val=0))

        # geant4MacroOutput stuff
        outputConfig.eventsPerTimeSlice = float(getConfigEntry(config, 'geant4MacroOutput', 'eventsPerTimeSlice', reqd=False, remove_spaces=True, default_val=0))
        outputConfig.subDivisionOfTime = float(getConfigEntry(config, 'geant4MacroOutput', 'subDivisionOfTime', reqd=False, remove_spaces=True, default_val=0))

        ebitParamsList = tuple(getConfigEntry(config, 'Run', 'beamList', reqd=True, remove_spaces=True).split(","))
        ebitParams = []
        for myebitParams in ebitParamsList:
            # Read in Beam information
            beamEnergy = float(getConfigEntry(config, myebitParams, 'beamEnergy', reqd=True, remove_spaces=True))
            breedingTime = float(getConfigEntry(config, myebitParams, 'breedingTime', reqd=True, remove_spaces=True))
            probeEvery = float(getConfigEntry(config, myebitParams, 'probeEvery', reqd=False, remove_spaces=True, default_val=0.0))
            ionEbeamOverlap = float(getConfigEntry(config, myebitParams, 'ionEbeamOverlap', reqd=False, remove_spaces=True, default_val=0.0))
            beamCurrent = float(getConfigEntry(config, myebitParams, 'beamCurrent', reqd=False, remove_spaces=True, default_val=0.0))
            beamRadius = float(getConfigEntry(config, myebitParams, 'beamRadius', reqd=False, remove_spaces=True, default_val=0.0))
            pressure = float(getConfigEntry(config, myebitParams, 'pressure', reqd=False, remove_spaces=True, default_val=0.0))
            ionTemperature = float(getConfigEntry(config, myebitParams, 'ionTemperature', reqd=False, remove_spaces=True, default_val=0.0))
            print("Beam Radius:", beamRadius)
            ebitParams.append(ebitChargeDistribution.EbitParams(breedingTime, probeEvery, ionEbeamOverlap, beamEnergy, beamCurrent, beamRadius, pressure, ionTemperature))

        # Gets the speciesList from under the [Run] headline
        speciesList = tuple(getConfigEntry(config, 'Run', 'speciesList', reqd=True, remove_spaces=True).split(","))
        species = []

        for myspecies in speciesList:
            # Collect all the species and parameters for each
            protons = int(getConfigEntry(config, myspecies, 'z', reqd=True, remove_spaces=True))
            nucleons = int(getConfigEntry(config, myspecies, 'nucleons', reqd=True, remove_spaces=True))
            population = float(getConfigEntry(config, myspecies, 'populationPercent', reqd=True, remove_spaces=True))
            chargeStates = list(map(int, getConfigEntry(config, myspecies, 'chargeStates', reqd=True, remove_spaces=True).split(",")))
            betaHalfLife = float(getConfigEntry(config, myspecies, 'betaHalfLife', reqd=False, remove_spaces=True, default_val=0.0))
            halfLife = float(getConfigEntry(config, myspecies, 'halfLife', reqd=False, remove_spaces=True, default_val=0.0))
            populationNumber = float(getConfigEntry(config, myspecies, 'populationNumber', reqd=False, remove_spaces=True, default_val=0.0))
            initSCITemp = float(getConfigEntry(config, myspecies, 'initSCITemp', reqd=False, remove_spaces=True, default_val=-1.0))
# NkT
            species.append(ebitChargeDistribution.Species(protons, nucleons, 0.0, betaHalfLife, population, chargeStates, halfLife, populationNumber, initSCITemp))
    else:
        print("Config file does not appear to exist : %s" % configFileName)
        sys.exit(1)

    if sum([*map(lambda x: x.initSCIPop, species)]) !=1.0:
        print("Population fractions do not sum to 1.0, please check the configuration file.")
        sys.exit(1)

    runSimulation(species, ebitParams, probeFnAddPop, outputConfig)
    return 0


def processCommandLine(args):
    #  Process command line arguments and then call runSimulation()
    outputConfig = OutputFormating()
    species = []
    ebitParams = []

    outputConfig.outputType = args.outputType
    outputConfig.outputFileName = args.outputFileName

    species.append(ebitChargeDistribution.Species(args.protons, args.nucleons, 0.0, 0.0, 1.0, args.chargeStates))
    ebitParams.append(ebitChargeDistribution.EbitParams(breedingTime=args.breedingTime, beamEnergy=args.beamEnergy, pressure=args.pressure, beamCurrent=args.beamCurrent, beamRadius=args.beamRadius, probeEvery=args.probeEvery))

    runSimulation(species, ebitParams, probeFnAddPop, outputConfig)

    return



def main():

    parser = argparse.ArgumentParser(description='EBIT Charge Breeding Simulation')

    parser.add_argument('-docs', dest='docs', required=False,
    					help='Obtain documentation about CBSim implementation. Options:\n  - physics\n  - timestepping\n')

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

    # Printing out documentation
    if args.docs:
    	try:
    		print (globals()['docs_'+args.docs].__doc__)
    	except KeyError:
    		print('Doc topic does not exist.')
    	sys.exit(1)


    if args.protons != 0:
        processCommandLine(args)
    else:
        processConfigFile(args.configFile)

if __name__ == "__main__":
    main()
