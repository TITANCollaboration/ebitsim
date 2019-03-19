# ebitsim
EBIT simulator for ion breeding


REQUIRES PYTHON3

I also highly recommend the use of PyPy3 with this.

-------------

Example usage :

Using command line arguments:

pypy3 ebitsim.py -z 51 -a 129 --chargeStates 46 47 48 49 50 51 --beamEnergy=25000.0 --breedingTime=12.0 --outputFile=BE25_I0.02_SB51_46_51.png

Using config file:

pypy3 ebitsim.py --configFile ebitsim.cfg

-------------
PyPy3 setup:

pypy3 -m ensurepip
pypy3 -m pip install matplotlib
-------------





Grabbed periodic table of elements csv from : https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee  -  Thank you!