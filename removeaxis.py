#! /usr/bin/env python
#
# Simple script to remove one data axis that has a certain CTYPE
#
# Magnus Persson 31.03.2010.

import os
import sys
from optparse import OptionParser as op
import pyfits as pf


#version
ver = 1.0

desc="""Remove an axis (STOKES) from fitsfile (data and header)"""
usage = "usage: %prog [options] fitsfile"
epilog = """usage: %prog [options] fitsfile, \'fitsfile\' \n\
should contain the relative path to the file as well\n\n\
Removes the axis named STOKES (and the CRPIX\#, CRVAL\# etc for that axis), assumes this is the fastest axis,\n\
i.e. d = d[0] when removing it"""

parser = op(usage=usage,description=desc,epilog=epilog, version="%prog " + str(ver))

# the options permitted
parser.add_option("-l", "--list", action="store_true", dest="l", help="list the header, and do nothin", metavar="LIST")
parser.add_option("-N", "--Nax", dest="N", help="remove axis number n", metavar="NAX")
#parser.add_option("-v", "--value", dest="v", help="value to asign to keyword", metavar="VALUE")
# time to parse
(options, args) = parser.parse_args()

# if no filename is supplied
if len(args)!= 1:
    parser.error("Incorrect number of arguments")
    parser.print_usage()

# the options
listhdr = options.l
axisno = str(options.N)

# and the file...
filename = str(args[0])

# get the fitsfile
hdulist = pf.open(filename, mode='update')
hdr, d = hdulist[0].header, hdulist[0].data

# what if no option is given?
if listhdr==True:
    print hdr
    sys.exit(1)

if d.shape[0] == 1 and int(axisno) == 4:
    d = d[0]
    try:
        del hdr['CROTA'+axisno]
        del hdr['CDELT'+axisno]
        del hdr['CRPIX'+axisno]
        del hdr['CRVAL'+axisno]
        del hdr['CTYPE'+axisno]
        del hdr['NAXIS'+axisno]
        hdr['NAXIS'] = 3
    except KeyError:
        print 'Some error with deleting och changing hdr keywords/values, check the code'
        sys.exit()


hdulist.flush()




