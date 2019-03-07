# ebitsim
EBIT simulator for ion breeding


REQUIRES PYTHON3

I also highly recommend the use of PyPy3 with this.

Example usage :

Using command line arguments:
pypy3 ebitsim.py -z 51 -a 129 --charge_states 39 40 41 42 43 --probe_every 0.001 --output_file somefile.png

Using config file:
pypy3 ebitsim.py --config_file ebitsim.cfg
