# OpenMM Setup

This is an application for configuring and running simulations with [OpenMM](http://openmm.org).  It provides
a graphical interface for selecting input files, cleaning up PDB structures, and setting simulation options.
It can then either save a script for running the simulation later, or directly run the simulation itself.

## Installation

If you are using [Anaconda Python](https://www.continuum.io/downloads) or [Miniconda](http://conda.pydata.org/miniconda.html),
you can install OpenMM Setup and all its dependencies with the command:

    conda install -c conda-forge openmm-setup

## Installing From Source

Alternatively you can install from source code.  OpenMM Setup is a Python application.  It requires the
following software to be installed:

* [OpenMM](http://openmm.org)
* [PDBFixer](https://github.com/pandegroup/pdbfixer)
* [Flask](http://flask.pocoo.org)

If you are using conda you can install all of these with a single command:

    conda install -c conda-forge openmm pdbfixer flask

To install OpenMM Setup, clone this repository, then type

    python setup.py install

## Using OpenMM Setup

You can then execute it by typing

    openmm-setup

The user interface is displayed through a web browser, but it is still a single user desktop application,
not a web application.  It should automatically open a web browser displaying the user interface.  If for
any reason that does not happen, open a browser yourself and point it to the address displayed in the console
window (usually http://127.0.0.1:5000). If the port 5000 is already in use by another program, a different
port can be selected with the -p option, for example

    openmm-setup -p 5001

