#  _________________________________________________________________________
#
#  PyUtilib: A Python utility library.
#  Copyright (c) 2008 Sandia Corporation.
#  This software is distributed under the BSD License.
#  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#  the U.S. Government retains certain rights in this software.
#  _________________________________________________________________________
#

def configure(installer):
    installer.description = """This script manages the installation of SPOT."""
    installer.default_dirname = 'tevaspot'
    #
    # Add repositories
    #
    installer.add_repository('nose', pypi='nose')
    installer.add_repository('pyutilib', pyname='PyUtilib', root='https://software.sandia.gov/svn/public/pyutilib', dev=True)
    installer.add_repository('pyomo', pyname='Pyomo', root='https://software.sandia.gov/svn/public/pyomo', dev=True)
    installer.add_repository('tevaspot', pyname='TevaSpot', root='https://software.sandia.gov/svn/teva/spot/spot', trunk='/packages/tevaspot', tag='/packages/tevaspot', dev=True)
    #
    # Add cmd scripts, which need to be customized
    # for a local installation.
    #
    installer.add_dos_cmd('sp.cmd')
    installer.add_dos_cmd('teva-spot.cmd')
    #
    # Declare the default installation directories
    #
    installer.default_windir = 'C:\\tevaspot'
    installer.default_unixdir = './tevaspot'
    #
    # Add auxdirs
    #
    installer.add_auxdir('tevaspot', 'etc/mod', '/etc/mod')
    installer.add_auxdir('tevaspot', 'examples/simple', '/examples/simple')
    installer.add_auxdir('tevaspot', 'data/Net3', '/packages/test/data/Net3')
    installer.add_auxdir('tevaspot', 'data/test1', '/packages/test/data/test1')
    installer.add_auxdir('tevaspot', 'util', '/packages/tevaspot/util')
    installer.add_auxdir('tevaspot', 'doc', '/doc/pub')
    installer.add_auxdir('tevaspot', 'test', '/packages/test/pyunit')
    #
    # Define the contents of the README.txt file
    #
    installer.README="""
#
# Installation generated by the spot_install script.
#
# This directory is managed with virtualenv, which creates a
# virtual Python installation.  If the 'bin' directory is put in
# user's PATH environment, then the bin/python command can be used to
# employ SPOT without further installation.
#
# Directories:
#   admin      Administrative data for maintaining this distribution
#   bin        Scripts and executables
#   data       Test data
#   dist       Python packages that are not intended for development
#   doc        SPOT documentation and tutorials
#   examples   SPOT examples
#   include    Python header files
#   lib        Python libraries and installed packages
#   etc        Math programming model files used for optimization
#   src        Python packages whose source files can be
#              modified and used directly within this virtual Python
#              installation.
#   Scripts    Python bin directory (used on MS Windows)
#   test       Test directory
#   util       Utility scripts (including spot_install)
#
"""
    #
    # Return the modified installer
    #
    return installer
