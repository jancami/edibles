************
Installation
************

Requirements
============

EDIBLES has the following strict requirements:

* `Python <https://www.python.org/>`_ 

* `Numpy <https://numpy.org/>`_ 

* `Astropy`_ 3.2 or later

* `Matplotlib`_

* `Pandas`_ 

* `Sherpa`_ 

EDIBLES also optionally depends on other packages for some features:

* `Scipy <https://www.scipy.org/>`_ 0.19 or later:  To power a variety of features in several
  modules (strongly recommended).

* `matplotlib <https://matplotlib.org/>`_ 2.2 or later:  To power a
  variety of plotting features (e.g. plotting apertures).

* `pytz`_ .

* `six`_ .

* `cycler`_ .

* `pyparsing`_ .

* `scipy`_ .



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


Testing an installed EDIBLES
==============================

Testing is currently under development.

The easiest way to test your installed version of EDIBLES is running
correctly is to use the :func:`edibles.test.test_basic_fit.testBasicFit` function:

.. doctest-skip::

    >>> import edibles.test.test_basic_fit as test
    >>> test.testBasicFit()

Note that this may not work if you start Python from within the
EDIBLES source distribution directory.



.. _github: https://github.com/jancami/edibles
