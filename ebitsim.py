#!/usr/bin/env python
# Jon Ringuette - 03/2019

import argparse
import configparser
import ebitChargeDistribution

def buildClassesFromArgs(args):
    ebitparams = EbitParams()

    return

def main():
  #  sg = SpawnerClass()
  #  workflow_objects = WorkFlow()

    # Probably accept config file location and maybe config file section as command line arguments
    # for the moment this is mostly just placeholder stuffs

    parser = argparse.ArgumentParser(description='EBIT Charge Breeding Simulation')
#    parser.add_argument('--init', dest='init', action='store_true',
#                        help='Initialize the database if this is the first time running this')
#    parser.add_argument('--server', dest='server', action='store_true',
#                        help='Start a Still Task Server')
#    parser.add_argument('--client', dest='client', action='store_true',
#                        help='Start a Still Task Client')
    parser.add_argument('--config_file', dest='config_file', required=False,
                        help="Specify the complete path to the config file, by default we'll use etc/still.cfg")
    parser.add_argument('-z', dest='protons', required=False,
                        help="Specify the num of protons")
    parser.add_argument('-a', dest='nucleons', required=False,
                        help="Specify the num of nucleons")
    parser.add_argument('--charge-states', nargs='+', dest='charge_states', type=int, required=False,
                        help="Specify charge states such as : --charge-states 32 33 34 35")
    parser.add_argument('--breeding_time', dest='breeding_time', required=False,
                        help="Specify breeding time (sec)")
    parser.add_argument('--beam_energy', dest='beam_energy', required=False,
                        help="Beam energy (eV)")
    parser.add_argument('--beam_current', dest='beam_current', required=False,
                        help="Beam current (Amps)")
    parser.add_argument('--beam_radius', dest='beam_radius', required=False,
                        help="Beam radius (um)")
    parser.add_argument('--pressure', dest='pressure', required=False,
                        help="ebit vacuum pressure (Torr / cm^3 I think)")
    parser.set_defaults(config_file="ebitsim.cfg")

    args, unknown = parser.parse_known_args()
    #sg.config_file = args.config_file

    print(args)
    buildClassesFromArgs(args)


if __name__ == "__main__":
    main()