EDIBLES
-------

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

.. image:: https://travis-ci.org/jancami/edibles.svg?branch=master
    :target: https://travis-ci.org/jancami/edibles


.. image:: https://readthedocs.org/projects/edibles/badge/?version=latest
    :target: https://edibles.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status




Software for EDIBLES data analysis and plotting.


Installation
------------
Note: The following assumes a Linux system, and a bash shell. 
 
1.
The first step should be to create a folder somewhere on your computer that all of your python scripts will be kept (not just this project). As an example, you could choose /home/python. 

The next step is to add some environment variables to your system. This can be done by modifying the hidden file .bashrc in your home directory.

The syntax to create the example variable VARIABLE is as follows::

    VARIABLE='/path/to/folder'
    export VARIABLE

There are 4 environment variables required to install and use the EDIBLES package:

- PYTHONPATH
- EDIBLES_DATADIR
- EDIBLES_DATARELEASE
- EDIBLES_PYTHONDIR

The 3 beginning with EDIBLES_ are there to tell the edibles package where your data is, your data version, and the explicit path to the edibles project folder.

The PYTHONPATH variable represents a user defined path to the parent folder where python looks when trying to import python packages (in addition to the default path where pip would install packages). Ideally, all of your python coding should be done within this folder. The edibles python package will also reside within this folder.

once completed, your .bashrc file should look something like this::

    PYTHONPATH='/home/python'
    export PYTHONPATH

    EDIBLES_DATADIR='/data/DR4'
    export EDIBLES_DATADIR

    EDIBLES_DATARELEASE='DR4'
    export EDIBLES_DATARELEASE

    EDIBLES_PYTHONDIR='/home/python/edibles'
    export EDIBLES_PYTHONDIR

After saving, open a new terminal and type echo $PYTHONPATH to print the PYTHONPATH to the screen (you can also do this with the others if you like).
 
2.
Inside your newly created Python folder, type::
 
    git clone https://github.com/jancami/edibles.git

to clone the current repository, and then type::
 
    cd edibles
    git pull
    
to make sure it is updated to the latest version.
 
I would then ensure you have all the recommended packages – the only one that is likely to be missing is Sherpa, but if there are more then that’s alright.
 
Depending on whether you use pip or conda::

    pip install sherpa
 
or::

    conda install -c sherpa sherpa
 
After that, the easiest way to test to see if things work is to run the test script!

These can be found at /home/python/edibles/edibles/tests/tests/

Start with test_EdiblesSpectrum, then test_basic_fit, then test_advanced_fit. These tests also lay out a basic use case for the package so far, and can be replicated in your analysis scripts.


