#!/usr/bin/env python
# Jon Ringuette - 03/2019

import argparse
import configparser
import ebitChargeDistribution
from ebitChargeDistribution import probeFnAddPop

import platform
import sys

def column(matrix, index):
    mycolumn = [i[index] for i in matrix]
    return mycolumn


def plotSpeciesResults(species, ebitparams):
    import matplotlib
    matplotlib.use("TKAgg")  # This is done as otherwise it goes wierd with a MacOS
    import matplotlib.pyplot as plt


    plt.rcParams['legend.loc'] = 'best'
    plt.ylabel('% Charge Breeding')
    plt.xlabel('Time after time')
    plt.title('Charge Breeding for stuff and things')
    for myspecies in species:
        for chargeStateResults in range(0, len(myspecies.results)):
            plt.plot(column(myspecies.results[chargeStateResults], 0), column(myspecies.results[chargeStateResults], 1), label=myspecies.chargeStates[chargeStateResults])
    plt.legend(framealpha=0.5)

    plt.savefig('mine.png', dpi=300)
    return


def buildAndRunClassesFromArgs(args, output='plot'):
    species = []

    # rewrite a bunch of this....
    species.append(ebitChargeDistribution.Species(args.protons, args.nucleons, 0.0, 0.0, 1.0, args.chargeStates))
    ebitparams = ebitChargeDistribution.EbitParams(breedingTime=args.breedingTime, beamEnergy=args.beamEnergy, pressure=args.pressure, beamCurrent=args.beamCurrent, beamRadius=args.beamRadius, probeEvery=args.probeEvery)
    ebitChargeDistribution.calcChargePopulations(species, ebitparams, probeFnAddPop)

    if output == 'plot':
        plotSpeciesResults(species, ebitparams)
    return


def main():

    parser = argparse.ArgumentParser(description='EBIT Charge Breeding Simulation')
#    parser.add_argument('--init', dest='init', action='store_true',
#                        help='Initialize the database if this is the first time running this')
#    parser.add_argument('--server', dest='server', action='store_true',
#                        help='Start a Still Task Server')
#    parser.add_argument('--client', dest='client', action='store_true',
#                        help='Start a Still Task Client')
    parser.add_argument('--config_file', dest='configFile', required=False,
                        help="Specify the complete path to the config file, by default we'll use ebitsim.cfg")
    parser.add_argument('-z', dest='protons', type=int, default=0, required=False,
                        help="Specify the num of protons")
    parser.add_argument('-a', dest='nucleons', type=int, required=False,
                        help="Specify the num of nucleons")
    parser.add_argument('--charge_states', nargs='+', dest='chargeStates', type=int, required=False,
                        help="Specify charge states such as : --charge-states 32 33 34 35")
    parser.add_argument('--breeding_time', dest='breedingTime', default=11.0, required=False,
                        help="Specify breeding time (sec)")
    parser.add_argument('--probe_every', dest='probeEvery', default=0.1, type=float, required=False,
                        help="Get results in this incriment (sec)")
    parser.add_argument('--beam_energy', dest='beamEnergy', default=7000, required=False,
                        help="Beam energy (eV)")
    parser.add_argument('--beam_current', dest='beamCurrent', default=0.02, required=False,
                        help="Beam current (Amps)")
    parser.add_argument('--beam_radius', dest='beamRadius', default=90e-4, required=False,
                        help="Beam radius (cm)")
    parser.add_argument('--pressure', dest='pressure', default=1e-10, required=False,
                        help="ebit vacuum pressure (Torr / cm^3 I think)")
    parser.set_defaults(configFile="ebitsim.cfg")

    args, unknown = parser.parse_known_args()
    # sg.config_file = args.config_file

    if platform.python_implementation() == 'CPython':
        print("*** !!WARNING!! : While this can run via normal cPython it is highly recommended that you run it via pypy3 for a HUGE speedup ***")

    if args.protons != 0:
        buildAndRunClassesFromArgs(args)
    else:
        print("we'll use a config instead...")

if __name__ == "__main__":
    main()