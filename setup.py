"""OpenMM Setup: User interface for preparing and running OpenMM simulations

This is an application for configuring and running simulations with OpenMM. It
provides a graphical interface for selecting input files, cleaning up PDB
structures, and setting simulation options. It can then either save a script
for running the simulation later, or directly run the simulation itself.
"""
from __future__ import print_function
import os
import sys
from os.path import relpath, join
from setuptools import setup, find_packages
DOCLINES = __doc__.split("\n")

########################
__version__ = '1.0'
VERSION = __version__
ISRELEASED = False
########################
CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: End Users/Desktop
License :: OSI Approved :: MIT License
Programming Language :: Python
Programming Language :: Python :: 3
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Scientific/Engineering :: Chemistry
Topic :: Scientific/Engineering :: Physics
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""


def find_package_data():
    files = []
    for root, dirnames, filenames in os.walk('openmmsetup'):
        for fn in filenames:
            files.append(relpath(join(root, fn), 'openmmsetup'))
    return files

setup(
    name='OpenMM-Setup',
    author='Peter Eastman',
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    version=__version__,
    license='MIT',
    url='https://github.com/peastman/openmm-setup',
    platforms=['Linux', 'Mac OS-X', 'Unix', 'Windows'],
    classifiers=CLASSIFIERS.splitlines(),
    packages=find_packages(),
    package_data={'openmmsetup': find_package_data()},
    zip_safe=False,
    install_requires=['flask', 'openmm >= 7.0', 'pdbfixer'],
    entry_points={'console_scripts': ['openmm-setup = openmmsetup.openmmsetup:main']})

