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

Building and installing manually
--------------------------------

EDIBLES is being developed on `github`_.  The latest development
version of the EDIBLES source code can be retrieved using git::

    git clone https://github.com/jancami/edibles.git

Then to build and install EDIBLES, run::

    cd edibles
    python setup.py install


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
