EDIBLES
-------

.. image:: https://travis-ci.org/jancami/edibles.svg?branch=master
    :target: https://travis-ci.org/jancami/edibles


.. image:: https://readthedocs.org/projects/edibles/badge/?version=latest
    :target: https://edibles.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status




Software for EDIBLES data analysis and plotting.



Installation from source
------------------------

The EDIBLES GitHub repository can only be `installed from source <https://edibles.readthedocs.io/en/latest/install.html>`_ for now.
 


to clone the current repository, and then type::
 
    cd edibles
    git pull

to make sure it is updated to the latest version.


the next step is to 'install' the package. There are two methods of doing this, and the one you choose should depend on how you will use this package. 

If you are using this package for the EDIBLES

I would then ensure you have all the recommended packages – the only one that is likely to be missing is Sherpa, but if there are more then that’s alright.
 
Depending on whether you use pip or conda::

    pip install sherpa
 
or::

    conda install -c sherpa sherpa
 
After that, the easiest way to test to see if things work is to run the test script!

These can be found at /home/python/edibles/edibles/tests/tests/

Start with test_EdiblesSpectrum, then test_basic_fit, then test_advanced_fit. These tests also lay out a basic use case for the package so far, and can be replicated in your analysis scripts.


