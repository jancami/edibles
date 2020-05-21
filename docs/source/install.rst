***************
Getting Started
***************

Requirements
============

EDIBLES has the following strict requirements:

* `Python <https://www.python.org/>`_ 

* `Numpy <https://numpy.org/>`_ 

* `Astropy <https://www.astropy.org>`_ 3.2 or later

* `Matplotlib <https://matplotlib.org/>`_

* `Pandas <https://pandas.pydata.org/>`_ 

* `Sherpa <https://sherpa.readthedocs.io/en/latest/>`_ 

EDIBLES also optionally depends on other packages for some features:

* `Scipy <https://www.scipy.org/>`_ 0.19 or later:  To power a variety of features in several
  modules (strongly recommended).

* `pytz <https://pypi.org/project/pytz/>`_ .

* `six <https://pypi.org/project/six/>`_ .

* `cycler <https://pypi.org/project/Cycler/>`_ .

* `pyparsing <https://pypi.org/project/pyparsing/>`_ .



EDIBLES depends on `pytest-astropy
<https://github.com/astropy/pytest-astropy>`_ (0.4 or later) to run
the test suite.

EDIBLES depends on `sphinx-astropy
<https://github.com/astropy/sphinx-astropy>`_ (0.4 or later) to build
the documentation.

Installation
============

Note: The following assumes a Linux system, and a bash shell. 


Basic Python Directory Setup
----------------------------


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

The 3 environment variables beginning with ``EDIBLES_`` are there to tell the edibles package where your data is, your data version, and the explicit path to the edibles project folder.

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

After saving, open a new terminal and type ``echo $PYTHONPATH`` to print the PYTHONPATH to the screen (you can also do this with the others if you like).
 

Installation of the EDIBLES Package
-----------------------------------

There are two methods to actually install the package, and the one you choose should depend on how you will use this package. 

For people that want to install and use the EDIBLES package code, you can install the package as you normally would using pip::

    git clone https://github.com/jancami/edibles.git
    cd edibles
    pip install .

This will install edibles alongside your other Python packages.

If you plan to develop the package, you should use the developer mode of pip installation::

    git clone https://github.com/jancami/edibles.git
    cd edibles
    pip install -e .


Your edibles Python package is now setup, and you are ready for some science!


Testing after installation
==========================

Testing is currently under development.

The easiest way to test your installed version of EDIBLES is running
correctly is to use the :func:`edibles.test.test_basic_fit.testBasicFit` function:

    >>> import edibles.test.test_basic_fit as test
    >>> test.testBasicFit()

Note that this may not work if you start Python from within the
EDIBLES source distribution directory.



.. _github: https://github.com/jancami/edibles
