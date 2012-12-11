#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
#  errors.py
#
#  errors classes for adapy
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
#  version 0.1b

import logger as _logger

class ParError(Exception):
    def __init__(self, value, reason=None):
        self.value = value
    def __str__(self):
        s1 = 'Parameter(s) \"{0}\" is(are) malformed or missing.'.format(self.value)
        return stylify(s1, fg='r')

class FitsError(Exception):
    def __init__(self, value, reason=None):
        self.value = value
    def __str__(self):
        s1 = 'Error reading fits file: {0}'.format(self.value)
        return stylify(s1, fg='r')
