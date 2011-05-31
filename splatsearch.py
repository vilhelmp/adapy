#! /usr/bin/env python
#
# Simple script to search Splatalogue from cli
#
# Magnus Vilhelm Persson 19.05.2011.



"""
This script queries Splatalogue.net

Returns a data structure with the requested information from Splatalogue.

TODO : if there are more than 1 result page, get it to read them all \
        in to the data structure. 
TODO : page output if it is long

TODO : solve the parsing, the removal of tags...

"""



# import optparser
from urllib2 import urlopen
from ClientForm import ParseResponse
from BeautifulSoup import BeautifulSoup as bfs
from scipy import array, where, nan, arange

#import os
#import sys
from optparse import OptionParser as op
#import pyfits as pf
#from string import lower

#version
ver = 1.0

desc="""Search the splatalogue compilation"""
usage = "usage: %prog [options]"
epilog = """usage: %prog [options]"""

parser = op(usage=usage,description=desc,epilog=epilog, version="%prog " + str(ver))

# the options permitted
parser.add_option("-f", "--frequency", dest="f", help="frequency range given as f1,f2" , metavar="FREQUENCY")
#parser.add_option("-v", "--velocities", dest="v", help="between what velocities to plot", metavar="VELOCITIES")
#parser.add_option("-b", "--box", dest="b", help="zoom to a box X,Y size", metavar="BOX")
#parser.add_option("-l", "--list", action="store_true", dest="l", help="list the header", metavar="LIST")
#parser.add_option("-k", "--force", action="store_true", dest="k", help="force overwriting of keyword", metavar="FORCE")

# time to parse
(options, args) = parser.parse_args()


# get the form from the splatalogue search page
response = urlopen("http://www.cv.nrao.edu/php/splat/b.php")
forms = ParseResponse(response, backwards_compat=False)
response.close()
form = forms[0]

#now change the 'from' and 'to' controls with the frequency range requested
# if no filename is supplied
if options.f == None:
    print 'no options input, running example with f1,f2 = 203.406,203.408'
    parser.print_usage()
    form['from'] = str(203.406)
    form['to'] = str(203.408)
else:
    # the options
    f1,f2 = [x for x in (options.f).split(',')]
    form['from'] = f1
    form['to'] = f2


form['frequency_units'] = ['GHz']

# 'click' the form
#result = urlopen(form.click()).read() # this does NOT works

# workaround, add 'submit=Search' to the request url
result = urlopen(form.click_request_data()[0]+'&submit=Search').read() # this works
# parse the results

bfs.NESTABLE_TABLE_TAGS['td'].append('sub')
bfs.NESTABLE_TABLE_TAGS['td'].append('sup')
bfs.NESTABLE_TABLE_TAGS['td'].append('font')

soup = bfs(result)


def findreplace(isoup, x,p=0):
    r = isoup.findAll(x)
    if p:
        [i.replaceWith('('+i.renderContents()+')') for i in r]
    else:
        [i.replaceWith(i.renderContents()) for i in r]

# these work, but only in this order (!?)
findreplace(soup, 'font')
findreplace(soup, 'sub')
findreplace(soup, 'sup',p=1)
findreplace(soup, 'a')
findreplace(soup, 'br')

results_table = soup.find('table',"results")

data = [[col.renderContents() for col in row.findAll('td')] for row in results_table.findAll('tr')]

# put all the stuff in arrays
data[0][0] = 'N'
col_names = data[0]
data = data[1:]

N, chem, name, cfreq, mfreq, qns,cdmsjpl_I,lovas_I,EL,llist = array(data).transpose()


## this is a not so good coded part...
# running i.strip('(').strip(')').split()
# does not seemt to work either
cfreq[where(cfreq == '&nbsp;')] = 'nan nan'
mfreq[where(mfreq == '&nbsp;')] = 'nan nan'

cfreq = array([i.split(' ') for i in cfreq])
mfreq = array([i.split(' ') for i in mfreq])

cfreqerr = cfreq[:,1]
mfreqerr = mfreq[:,1]
cfreq = cfreq[:,0]
mfreq = mfreq[:,0]


mfreqerr = [i.strip('(').strip(')') for i in  mfreqerr]
cfreqerr = [i.strip('(').strip(')') for i in  cfreqerr]

cfreqerr = array(cfreqerr, dtype='float')
cfreq = array(cfreq, dtype='float')
mfreqerr = array(mfreqerr, dtype='float')
mfreq = array(mfreq, dtype='float')

print '{0:2} {1:10}\t{2:10}\t{3:9}\t{4:6} '.format('N', 'Form', 'Freq', 'FErr', 'List')
for i in arange(len(N)):
    print '{0:2} {1:10}\t{2:10}\t{3:9}\t{4:6} '.format(N[i], chem[i], cfreq[i], cfreqerr[i], llist[i])
