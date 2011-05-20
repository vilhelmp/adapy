#! /usr/bin/env python
#
# Simple script to add a keyword and value to a fitsheader
#
# Magnus Persson 31.03.2010.
"""
TODO:
- what if one just wants to add a comment?

"""

import os
import sys
from optparse import OptionParser as op
import pyfits as pf

#version
ver = 1.0

desc="""Add a keyword and value to fitsheader"""
usage = "usage: %prog [options] fitsfile"
epilog = """usage: %prog [options] fitsfile, \'fitsfile\' \n\
should contain the relative path to the file as well"""

parser = op(usage=usage,description=desc,epilog=epilog, version="%prog " + str(ver))

# the options permitted
parser.add_option("-k", "--keyword", dest="k", help="keyword to add" , metavar="KEYWORD")
parser.add_option("-v", "--value", dest="v", help="value to asign to keyword", metavar="VALUE")
parser.add_option("-c", "--comment", dest="c", help="comment to add to keyword", metavar="COMMENT")
parser.add_option("-l", "--list", action="store_true", dest="l", help="list the header", metavar="LIST")
parser.add_option("-f", "--force", action="store_true", dest="f", help="force overwriting of keyword", metavar="FORCE")

# time to parse
(options, args) = parser.parse_args()

# if no filename is supplied
if len(args) != 1:
	parser.error("Incorrect number of arguments")
	parser.print_usage()
# the options
keyword = options.k
value = options.v
comment = options.c
listhdr = options.l
force = options.f

# and the file...
filename = str(args[0])

# get the fitsfile
hdulist = pf.open(filename, mode='update')
# and from that, the primary header
hdr = hdulist[0].header

# what if no option is given?
if keyword==None and value==None and comment==None and listhdr==None:
	print 'No options given, will just list header..'
	listhdr = True
if listhdr==True:
	print hdr
	sys.exit(1)

# NOW, assuming that some keyword,value where given..

# does the keyword we want to 
# add already exist?
test = hdr.has_key(str(keyword))
if test == 1 and force == None:
	# if it does exist and force is false
	ans = raw_input('The keyword '+str(keyword) +
		" exists (with value: "+str(value)+"\
) do you really want to continue and replace the value?(y/n)")
	if ans in ['y','Y','yes']:
		pass
	elif ans in ['n','N','no']:
		sys.exit(1)
elif test == 1 and force == True:
	pass

# either they answered 'y', it did not exist or the force-flag was set
try:
	value = float(value)
except (ValueError):
	value = str(value)

# update the fitsheader
hdr.update(str(keyword), value)

# save the changes back to the fits file
hdulist.flush()
