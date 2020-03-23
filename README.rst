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
 
1.
The first step should be to create a folder somewhere on your computer that all of your python scripts will be kept (not just this project). As an example, you could choose /home/python. You will then want to permanently add this folder to your PYTHONPATH. This can be done by editing the /home/.bashrc or /home/.bash_profile files to include a line that says::

    PYTHONPATH="/home/python:$PYTHONPATH"
    export PYTHONPATH

(or whatever your folder is). Close and reopen the terminal and try echoing them (echo $PYTHONPATH). This is the first step and needs to be done for you to be able to import any handwritten/edited Python scripts. The PYTHONPATH is not usually set by default, so a blank line is what you would expect when you try echoing before adding these lines.  
 
2.
While you are working on the PYTHONPATH environment variable, you should also setup the EDIBLES_DATADIR, EDIBLES_DATARELEASE, and EDIBLES_PYTHONDIR environment variables. In total, you will probably end up adding::
 
    PYTHONPATH=' /home/python:$PYTHONPATH'
    export PYTHONPATH
 
    EDIBLES_DATADIR='/data/DR4_fits'
    export EDIBLES_DATADIR
 
    EDIBLES_DATARELEASE='DR4
    export EDIBLES_DATARELEASE
 
    EDIBLES_PYTHONDIR=' /home/python/edibles'
    export EDIBLES_PYTHONDIR
 
or something like this to your ~/.bashrc file. /data/DR4_fits is a folder where all my fits files are stored.
 
Once you have this part done, and you can echo the 4 environment variables, you can move on to actually installing the package.

Note: These environment variables can also be setup in an editor like sublime (there may also be others).
 
3.
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


