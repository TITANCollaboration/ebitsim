# ebitsim
EBIT simulator for ion breeding

**REQUIRES Python3

I also highly recommend the use of PyPy3 with this.

-------------
Command line interface:
```
pypy3 ebitsim.py -h
usage: ebitsim.py [-h] [-docs DOCS] [--configFile CONFIGFILE]
                  [--outputType OUTPUTTYPE] [--outputFileName OUTPUTFILENAME]
                  [-z PROTONS] [-Z PROTONS] [-a NUCLEONS]
                  [--chargeStates CHARGESTATES [CHARGESTATES ...]]
                  [--breedingTime BREEDINGTIME] [--probeEvery PROBEEVERY]
                  [--beamEnergy BEAMENERGY] [--beamCurrent BEAMCURRENT]
                  [--beamRadius BEAMRADIUS] [--pressure PRESSURE]

EBIT Charge Breeding Simulation

optional arguments:
  -h, --help            show this help message and exit
  -docs DOCS            Obtain documentation about CBSim implementation.
                        Options: - physics - timestepping
  --configFile CONFIGFILE
                        Specify the complete path to the config file, by
                        default we'll use ebitsim.cfg
  --outputType OUTPUTTYPE
                        Specify how to output the data, defaults to a
                        matplotlib graph, csv will be available at some point
  --outputFileName OUTPUTFILENAME
                        Specify filename for csv or png file, please include
                        .csv or .png extension
  -z PROTONS            Specify the num of protons
  -Z PROTONS            Specify the num of protons
  -a NUCLEONS           Specify the num of nucleons
  --chargeStates CHARGESTATES [CHARGESTATES ...]
                        Specify charge states such as : --charge-states 32 33
                        34 35
  --breedingTime BREEDINGTIME
                        Specify breeding time (sec)
  --probeEvery PROBEEVERY
                        Get results in this incriment (sec)
  --beamEnergy BEAMENERGY
                        Beam energy (eV)
  --beamCurrent BEAMCURRENT
                        Beam current (Amps)
  --beamRadius BEAMRADIUS
                        Beam radius (cm)
  --pressure PRESSURE   ebit vacuum pressure (Torr / cm^3 I think)
```

-------------

Example usage :

Using command line arguments:
```
pypy3 ebitsim.py -z 51 -a 129 --chargeStates 46 47 48 49 50 51 --beamEnergy=25000.0 --breedingTime=12.0 --outputFile=BE25_I0.02_SB51_46_51.png
```
Using config file:
```
pypy3 ebitsim.py --configFile ebitsim.cfg
```
-------------

PyPy3 setup for using matplotlib:

```
pypy3 -m ensurepip

pypy3 -m pip install matplotlib
```
-------------




Grabbed periodic table of elements csv from : https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee  -  Thank you!