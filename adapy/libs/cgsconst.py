#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
#  cgsconst.py
#
#  
#
#  Copyright 2012 Magnus Persson <http://vilhelm.nu>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  version 0.1 alpha
#
#
"""
CGS Constants for Astronomical Data Analysis

Module for various useful CGS constants.
CODATA values for most of them, i.e. from scipy.constants but converted
to CGS-units. With unit and description attributes.
Inspired by natconst.py by Jes Joergensen, and IDL's natconst.pro by
C. P. Dullemond.
"""

#-----------------------------------------------------------------------
#Top of the list TODO:
"""
To add: 

cm / light year


"""

"""
----[ Change log ]----


02 Oct 2012 - file created!



"""

########################################################################
# CONSTANTS OBJECT
#
class Constant(float):
    """
    Constant float type that can have a unit (unit)
    and a description (desc)

    Call : W = Constant(value, 'unit', desc = 'description')
    where 'unit' and desc, are optional.
    Then 'W' just returns the value assigned to it, and can
    be used in calculations as a normal float.
    e.g. W * 2e3 / 5.
    """

    def __new__(self, value, *args, **kwargs):
        # return without the *args and **kwargs
        # to make sure no error is raised in __new__
        return super(Constant, self).__new__(self, value)

    def __init__(self, value, *args, **kwargs):
        # store the different arguments into the class
        self.unit = args[0] if args else None
        self.desc = kwargs.pop('desc', None)


########################################################################
# CONSTANTS
# from "Physics Handbook for Science and Engineering" 2002
# and some  from script natconst.py by Jes Joergensen which
# some in turn are from C. P. Dullemond's IDL scripts
#

# Astronomy constants in cgs units

MEARTH = Constant(5.977e27, 'g',     desc = 'Mass of the Earth')
MMOON  = Constant(7.349e25, 'g',     desc = 'Mass of the Moon')
#AU   = Constant(1.496E13,   'cm / AU',    desc = 'Astronomical Unit') # old
AU   = Constant(1.495978707E13, 'cm / AU',    desc = 'Astronomical Unit') # new from IAU meeting int 2012
LSUN = Constant(3.8525e33,  'erg.s', desc = 'Solar luminosity')
RSUN = Constant(6.96e10,    'cm',    desc = 'Solar radius')
MSUN = Constant(1.99e33,    'g',     desc = 'Mass of the Sun')
PC   = Constant(3.08572e18, 'cm / parsec',    desc = 'Parsec')

# Physics constants in cgs units
from scipy import constants
CC = Constant(constants.c * 1e2,  'cm/s',         desc = 'Speed Of Light')
KK = Constant(constants.k * 1e7,  'erg/K',        desc = 'Bolzmann\'s constant')
HH = Constant(constants.h * 1e7,  'erg.s',        desc = 'Planck\'s constant')
GG = Constant(constants.G * 1e3,  'cm^3/(g.s^2)', desc = 'Gravitational constant')
MP = Constant(constants.m_p * 1e3,      'g?',     desc = 'Mass of proton')
ME = Constant(constants.m_e * 1e3,      'g?',     desc = 'Mass of electron')
try:
    SS = Constant(constants.Stefan_Boltzmann * 1e3, 'erg/cm^3/K^4', desc = 'Stefan-Boltzmann constant')
except (AttributeError):
    SS = Constant(constants.Stefan_Bolzmann * 1e3, 'erg/cm^3/K^4', desc = 'Stefan-Boltzmann constant')

#~ EE  = Constant(constants.,      '?',            desc = 'Unit charge')
#~ ST  = Constant(constants.,      'cm^2',         desc = 'Thomson cross-section')

# Alternatively (from natconst.py)
#~ CC  = Constant(2.9979e10,  'cm/s',         desc = 'Speed Of Light')
#~ KK  = Constant(1.3807e-16,  'erg/K',        desc = 'Bolzmann\'s constant')
#~ HH  = Constant(6.6262e-27,  'erg.s',        desc = 'Planck\'s constant')
#~ GG = Constant(6.672e-8,         'cm^3/(g.s^2)', desc = 'Grav. const.')
#~ MP  = Constant(1.6726e-24,      'g?',           desc = 'Mass of proton')
#~ ME  = Constant(9.1095e-28,      'g?',           desc = 'Mass of electron')
#~ SS  = Constant(5.6703e-5,       'erg/cm^3/K^4', desc = 'Stefan-Boltzmann constant')
EE  = Constant(4.8032e-10,      '?',            desc = 'Unit charge')
ST  = Constant(6.6524e-25,      'cm^2',         desc = 'Thomson cross-section')

AA  = Constant(7.5657e-15,      '?',            desc = '4 ss / CC')

#     Gas constants in cgs units
MUH2 = Constant(2.3000e0,   'amu?', desc='Mean molec weight H2+He+Metals')
#
#     Alternative units
#
EV   = Constant(1.6022e-12, 'erg',  desc = 'Electronvolt')
KEV  = Constant(1.6022e-9,  'erg',  desc = 'Kilo electronvolt')
MICR = Constant(1.e-4,      'cm',   desc = 'Micron')
KM   = Constant(1.e5,       'cm',   desc = 'Kilometer')
ANGS = Constant(1.e-8,      'cm',   desc = 'Angstroem')

#     Time units
YEAR = Constant(3.1536e7, 's', desc='Year')
HOUR = Constant(3.6000e3, 's', desc='Hour')
DAY  = Constant(8.64e4,   's', desc='Day')

