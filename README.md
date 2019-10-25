# ebitsim
EBIT simulator for ion breeding

**REQUIRES Python3

I also highly recommend the use of PyPy3 with this.

-------------

Example usage :

Using command line arguments:

pypy3 ebitsim.py -z 51 -a 129 --chargeStates 46 47 48 49 50 51 --beamEnergy=25000.0 --breedingTime=12.0 --outputFile=BE25_I0.02_SB51_46_51.png

Using config file:

pypy3 ebitsim.py --configFile ebitsim.cfg

-------------

PyPy3 setup for using matplotlib:

pypy3 -m ensurepip

pypy3 -m pip install matplotlib

-------------

Windows 10 Pro 64bit instructions
 
Download Virtual Box from https://www.virtualbox.org/

Download Vagrant from https://www.vagrantup.com/

Restart PC (Installation may ask to reboot).

Download ZIP file from the link https://github.com/mineselectroweakgroup/ebitsim_vagrant

Unzip file to your preferred location, Open Command prompt and cd to this location.

Use command ‘vagrant up’; This will download all the files required.

Use command ‘vagrant ssh’ to login to virtual machine.




Grabbed periodic table of elements csv from : https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee  -  Thank you!