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

"""



# import optparser
from urllib2 import urlopen
from ClientForm import ParseResponse

# get the form from the splatalogue search page
response = urlopen("http://www.cv.nrao.edu/php/splat/b.php")
forms = ParseResponse(response, backwards_compat=False)
response.close()
form = forms[0]

#now change the 'from' and 'to' controls with the frequency range requested
form['from'] = str(203.406)
form['to'] = str(203.408)
form['frequency_units'] = ['GHz']

# 'click' the form
result = urlopen(form.click()).read() # this does NOT works

# workaround, add 'submit=Search' to the request url
result = urlopen(form.click_request_data()[0]+'&submit=Search').read() # this works
# parse the results


####### Alternative 1
from lxml import html

# breaks the <br> tag too, have to change the xpath string input parameter
table = html.fromstring(result)
# get everything directly
# change the xpath parameters. look at 
# http://www.stonehenge.com/merlyn/LinuxMag/col92.html for reference
findTR = '//table[@class="results"]/tr'
findTD = './td[position()=1]/a/@href | ./td[position()>1]/text() | self::node()[position()=1]/td/text()'
rows = [row.xpath(findTD) for row in table.xpath(findTR)]
for row in rows:
    print len(row)*'%7s\t' % tuple(row)
for row in table.xpath('//table[@class="results"]/tr'):
    column = row.xpath('./td[position()=1]/a/@href | ./td[position()>1]/text() | self::node()[position()=1]/td/text()')
    l = len(column)
    print l*'%7s\t' % tuple(column)


####### Alternative 2
# works well, code is really short
# need to remove <br> and &nbsp; <a href etc tags
# since &nbsp: means no data, replace with ND, Nan or similar
# inspired by http://blog.jgc.org/2009/11/parsing-html-in-python-with.html
from BeautifulSoup import BeautifulSoup as bfs

soup = bfs(result)
results_table = soup.find('table',"results")
data = [[col.renderContents() for col in row.findAll('td')] for row in results_table.findAll('tr')]
  
####### Alternative 3
# need file in parse cmd, I have result string
from xml.etree.ElementTree import ElementTree

doc = ElementTree().parse(result)

for t in doc.findall('.//table'):
    for tr in t.findall('./tr/')[1:]: # skip the header row
    tds = tr.findall('./td')
    print tds[0][0].attrib['href'], tds[1].text.strip(), tds[2].text.strip()


 




        



## here I am cheating because the above did not work properly 
# what I am doing instead is just taking the search string, 
# and manipulating it, the control "submit=Search" is not included in the "click"


#~ f = ['http://www.cv.nrao.edu/php/splat/c.php?calcIn=',
 #~ 'data_version=v2.0',
 #~ 'from=203.406',
 #~ 'to=203.408',
 #~ 'frequency_units=GHz',
 #~ 'energy_range_from=',
 #~ 'energy_range_to=',
 #~ 'lill=on',
 #~ 'tran=',
 #~ 'submit=Search',
 #~ 'no_atmospheric=no_atmospheric',
 #~ 'no_potential=no_potential',
 #~ 'no_probable=no_probable',
 #~ 'displayLovas=displayLovas',
 #~ 'displaySLAIM=displaySLAIM',
 #~ 'displayJPL=displayJPL',
 #~ 'displayCDMS=displayCDMS',
 #~ 'displayToyaMA=displayToyaMA',
 #~ 'displayOSU=displayOSU',
 #~ 'displayRecomb=displayRecomb',
 #~ 'displayLisa=displayLisa',
 #~ 'displayRFI=displayRFI',
 #~ 'ls1=ls1',
 #~ 'ls5=ls5',
 #~ 'el1=el1']

#~ f = '&'.join(f)
#~ result = urlopen(f).read() # this works

