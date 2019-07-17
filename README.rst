Welcome to NumBAT!
--------------------

NumBAT, the Numerical Brillouin Analysis Tool, integrates electromagnetic and acoustic mode solvers to calculate the interactions of optical and acoustic waves in waveguides.


Origin
------

NumBAT was developed by Bjorn Sturmberg, Kokou Dossou, Chris Poulton and Michael Steel in a collaboration between Macquarie University and the University of Technology Sydney, as part of the Australian Research Council Discovery Project DP160101691.


Documentation
-------------

Is hosted on `ReadTheDocs <http://numbat-au.readthedocs.io/en/latest/>`_.


Compatability
-------------

NumBAT has been developed on Ubuntu 18.04 with the following package versions: Python 3.6, Numpy 1.16.2, Suitesparse 4.4.6, and Gmsh 3.0.6.

It has also been successfully installed by users on Debian, RedHat and on Windows 10 (installing Ubuntu after enabling the Windows Subystem for Linux) and with different versions of packages, but these installations have not been as thoroughly documented so may require user testing.

We also provide a `docker image <https://hub.docker.com/r/morblockdock/numbat>`_ that allows for easy cross platform operation. Information on how to use the docker image is contained in docker_notes.md file.


Installation is as simple as ::

    $ git clone https://github.com/bjornsturmberg/NumBAT.git 

Followed by ::

    $ sudo ./setup.sh

from the NumBAT/ directory.

For optimal results ::

	$ cp NumBAT/backend/NumBATstyle.mplstyle ~/.config/matplotlib/stylelib/

or replace plt.style.use('NumBATstyle') in NumBAT/backend/plotting.py with your own prefered matplotlib style file.