#! /usr/bin/env python
#
# Simple script search Splatalogue
#
# Magnus Vilhelm Persson 19.05.2011.



# import optparser
from urllib2 import urlopen
from ClientForm import ParseResponse

# parse the cmd line parameters, from, to,


# get the form from the splatalogue search page
response = urlopen("http://www.cv.nrao.edu/php/splat/b.php")
forms = ParseResponse(response, backwards_compat=False)
form = forms[0]

#now change the 'from' and 'to' controls
form['from'] = str(203.400)
form['to'] = str(203.410)
form['frequency_units'] = ['GHz']

# 'click' the form
result = urlopen(form.click()).read()
# workaround, add 'submit=Search' to the request url
result = urlopen(form.click_request_data()[0]+'&submit=Search').read()
# parse the results



## here I am cheating because the above did not work properly 
# what I am doing instead is just taking the search string, 
# and manipulating it, the control "submit=Search" is not included in the "click"

f = ['http://www.cv.nrao.edu/php/splat/c.php?calcIn=',
 'data_version=v2.0',
 'from=',
 'to=',
 'frequency_units=GHz',
 'energy_range_from=',
 'energy_range_to=',
 'lill=on',
 'tran=',
 'submit=Search',
 'no_atmospheric=no_atmospheric',
 'no_potential=no_potential',
 'no_probable=no_probable',
 'displayLovas=displayLovas',
 'displaySLAIM=displaySLAIM',
 'displayJPL=displayJPL',
 'displayCDMS=displayCDMS',
 'displayToyaMA=displayToyaMA',
 'displayOSU=displayOSU',
 'displayRecomb=displayRecomb',
 'displayLisa=displayLisa',
 'displayRFI=displayRFI',
 'ls1=ls1',
 'ls5=ls5',
 'el1=el1']

'&'.join(f)
result = urlopen(f).read()
