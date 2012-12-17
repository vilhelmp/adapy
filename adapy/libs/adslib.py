#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
#  adslib.py
#
#  Module to search the ads
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
#  version 0.0.1a

"""
Script to search the NASA ADS directory

Need :  o scipy
        o mechanize module (standard in Python >2.6/2.7?)
        o urllib2 module (standard Python module, required(?) by mechanize)
        o beautiful soup/xml (xml standard in Python >2.6/2.7?)

"""

"""
ADSlib - Python Module to interact with NASA ADS
at

http://adswww.harvard.edu/

OR one of the mirrors

http://cdsads.u-strasbg.fr/
http://ukads.nottingham.ac.uk/
http://esoads.eso.org/
http://ads.ari.uni-heidelberg.de/
http://ads.inasan.ru/
http://ads.mao.kiev.ua/
http://ads.astro.puc.cl/
http://ads.nao.ac.jp/
http://ads.bao.ac.cn/
http://ads.iucaa.ernet.in/
http://ads.arsip.lipi.go.id/
http://saaoads.chpc.ac.za/
http://ads.on.br/




"""

"""
----[ Change log ]----

* 2012 Dec 15 
    Code cleanup. 

* 2012 Oct 29 
    Now only uses mechanize module(!) Yay.

* 2012 Oct 02
    File created.



"""

"""
NOTES:
# advanced search query
abstract_service.html

# quick search
index.html

"""


mirrors = [
        'http://adswww.harvard.edu/',
        'http://cdsads.u-strasbg.fr/',
        'http://ukads.nottingham.ac.uk/',
        'http://esoads.eso.org/',
        'http://ads.ari.uni-heidelberg.de/'
        'http://ads.inasan.ru/',
        'http://ads.nao.ac.jp/',
        'http://ads.iucaa.ernet.in/',
        'http://ads.arsip.lipi.go.id/',
        'http://saaoads.chpc.ac.za/',
        'http://ads.on.br/'
        ]

advanced_q = 'abstract_service.html'

def search(query, **kwargs):
    """
    query       :  Normal string to ADS
                   or dictionary for advanced search
    
    """
    
    
    ### test code to get it up and running
    
    # main access
    # TODO : either access via Z39.50 or via URLlib/mecahnise etc
    
    # wishlist
    # TODO : simple search
    # TODO : advanced search
    # TODO : browse
    
    
    import locale
    # this reads the environment and inits the right locale
    locale.setlocale(locale.LC_ALL, "")
    
    
    try:
        # the mechanize module exports urllib2 as well...
        import mechanize
        import urllib
    except (ImportError):
        print 'You need the \"mechanize\" and urllib module'
        ' for this script to work.'
    
    try:
        from BeautifulSoup import BeautifulSoup as bfs
    except (ImportError):
        print 'You need the BeautifulSoup module...'
    
    
    import scipy
    import sys
    
    #from string import lower, upper
    # search URL
    # http://adsabs.harvard.edu/cgi-bin/nph-basic_connect?qsearch=The+Search+String
    
    # to parse the search string from "The Search String" to "The+Search+String"
    # urllib.quote(url, safe=":/")
    
    ############################################
    ######## GET THE FORM
    
    #~  Ping to know which server to use.
    working_mirror = 0
    
    #~ while not got_reply:
       #~ try:
           #~ # try to get the form
           #~ response = mechanize.urlopen(mirrors[working_mirror] + types_q[search_type])
       #~ except mechanize.URLError:
           #~ # if we can't get it, try another mirror
           #~ if not i < len(mirrors):
               #~ break
           #~ else:
               #~ working_mirror += 1
           #~ pass
       #~ else:
           #~ got_reply = True
    #~ 
    #~ if not got_reply and working_mirror >= len(mirrors):
           #~ # TODO log output
           #~ sys.stderr.write('ERROR :  You have to be connected to the internet to access the NASA ADS database and it has to be online (on all mirrors).')
    #~ else:
            #~ # TODO log output
        #~ print ('got reply from : {0}'.format(mirrors[working_mirror]))
    
    
    
    
    #~  Then check if going for the advanced interface.
    #~ advanced = int((type(query) == type({}))
    if advanced:
        # ADVANCED QUERY 
        response = mechanize.urlopen(mirrors[working_mirror] + advanced_q)
        forms = mechanize.ParseResponse(response, backwards_compat=False)
        response.close()
        form = forms[0]
        #~ if arg.has_key('dbg_snd_form'): # for test purposes
        #~ return form
        #~ form['qsearch'] = '^Persson 2012'
        
        ######## SUBMIT FORM
        #~ clicked_form = form.click()
        
        #~ result = mechanize.urlopen(clicked_form)
        
        pass
        
    elif not advanced:
        # SIMPLE QUERY 
        baseurl = (mirrors[working_mirror] + 
        'cgi-bin/nph-basic_connect?qsearch=')
        
        result = mechanize.urlopen( urllib.quote(baseurl + query, safe = ":/=?^") )
        # test below
        data = urllib.urlencode({'qsearch' : '^Persson'})
        baseurl = (mirrors[working_mirror] + 
        'cgi-bin/nph-basic_connect?')
        f = urllib.urlopen(baseurl, data)
    ############################################
    ######## PARSE RESULTS
    
    page = result.readlines()
    result.close()
    
    # start parsing the results
    t = bfs(' '.join(page))
    tables = t.findAll('table')
    
    r = tables[1].findAll('td')[0]
    y = r.findAll('strong')[0].contents[0]
    nres = int(y)
    if nres<1:
        return 0
    
    # get table with results
    resulttable = tables[2]
    # get the rows of the table
    rows = resulttable.findAll('tr')
    # get each result entry per list item
    entries = [rows[i:i+3][1:] for i in scipy.arange(2,57,3)][:-1]

    ############################################
    ######## GET RESULTLIST

    ###### the problem with this is that web is in UNICODE, 
    # ie. Jørgensen, æ and åäö and ßü etc are represented by funny numbers and '\'
        
    #resultlist = [_Result(i) for i in entries]
    return _Resultlist(entries)


############################################
######## DEFINE RESULT(S) OBJECT


class _Resultlist:
    """
    Internal object to represent the result list
    """
    def __init__(self, entries):
        self.resultlist = [_Result(i) for i in entries]
    def sort(self,sortkey = 'author', reverse_bool = False):
        from operator import itemgetter, attrgetter
        #~ sorted(resultlist, key=attrgetter('author'), reverse=True)
        return sorted(self.resultlist, key=attrgetter(sortkey), reverse = reverse_bool)
    def __str__(self):
        printlist = []
        for i in self.resultlist[:-1]:
            printlist.append('Author : {0.author}\n' 
            'Title : {0.title}\n' 
            'Score : {0.ads_score}\n'.format(i))
        return '\n'.join(printlist)

class _Result:
    """
    Internal object to represent each result
    """
    def __init__(self, entry):
    #~ def __init__(self, author, 
                        #~ authors, 
                        #~ title, 
                        #~ score, 
                        #~ bibcode,
                        #~ pubdate,
                        #~ links):
        #~ self.author = author
        #~ self.authorlist = authors
        #~ self.title = title
        #~ self.score = score
        #~ self.bibcode = bibcode
        #~ self.pubdate = pubdate  # parse?
        #~ self.links = links      # dictionary of all the links
        #
        td_tags0 = entry[0].findAll('td')
        self.bibcode = td_tags0[1].findAll('input')[0]['value'].encode('UTF-8')
        self.url_abstract_page = td_tags0[1].findAll('a')[0]['href'].encode('UTF-8')
        self.ads_score = float(td_tags0[3].contents[0].encode('UTF-8'))
        self.rank = 100 - self.ads_score
        self.pubdate = td_tags0[4].contents[0].string.encode('UTF-8')
        self.pubday = self.pubdate[:2]
        self.pubyear = self.pubdate[3:]
        #
        self.links = dict()
        for link in td_tags0[5].findAll('a'):
            self.links[link.string.encode()] = link['href'].encode('UTF-8')
        #
        td_tags1 = entry[1].findAll('td')
        
        # second part of the result entry
        self.title = td_tags1[3].contents[0].string.encode('UTF-8')
        # still in unicode
        # TODO need to convert to normal UTF, not unicode
        authors = td_tags1[1].contents[0].encode('UTF-8').split(';')
        if authors[-1] == ' ':
            # so, if the last entry in the authorlist is empty, means
            # it split a ';', which in turn means there are more 
            # authors, need to add that part...
            authors[-1] = td_tags1[1].contents[1].contents[0].encode('UTF-8') + ', COAuth'
        #
        self.authors = [i.split(',') for i in authors]
        self.author = ', '.join(self.authors[0])
        #
        #~ self.
    def __repr__(self):
        return repr([self.author, self.authors, self.title, self.url_abstract_page, self.ads_score, self.links, self.bibcode, self.pubdate])
    def _returnlist_(self):
        return [self.author, self.authors, self.title, self.url_abstract_page, self.ads_score, self.links, self.bibcode, self.pubdate]

#~ # second part of the result entry
#~ title = td_tags1[3].contents[0].string.replace(u'\xa0', u' ').encode()
#~ # still in unicode
#~ # TODO need to convert to normal UTF, not unicode
#~ authors = td_tags1[1].string.replace(u'\xa0', u' ').encode().split(';')
#~ authors = [i.split(',') for i in authors]
#~ author = authors[0]




############################################
######## RETURN SORTABLE OBJECT LIST

############################################
######## HOW TO SORT RESULTS
# needs Python 2.6 at least
#~ from operator import itemgetter, attrgetter
#~ 
#~ # now to sort it, just use one of the keys
#~ # score, high to low
#~ sorted(resultlist, key=attrgetter('author'), reverse=True)
#~ 
#~ # cmp=locale.strcoll new and untested addition
#~ 
#~ # authors alphabetical order first and then by score
#~ # i.e. sort by score if same first author
#~ sorted(resultlist, key=attrgetter('ads_score','authors'), reverse=True)


########################################################################
######## NOTES 

### FIELDS
# bibcode
# title
# authors
# score
# pubdate
# possilbe (quick)links : 
#           A Abstract 
#           C CITATIONS
#           D On-line Data
#           E EJOURNAL
#           F Printable Article
#           G Gif Images
#           H HEP/Spires Information
#           I Author Comments
#           L Library Entries
#           M Multimedia
#           N NED Objects
#           O Associated Articles
#           P PDS datasets
#           R REFERENCES
#           S SIMBAD Objects
#           T TOC
#           U Also read
#           X arXiv e-print
#           Z Abstract Custom




#~ Z39.50 or via URLlib/mecahnise etc

#~ The ADS server supports the following services:
#~ 
#~ Initialization
#~ Search
#~ Present
#~ 
#~ Production Server
#~ Domain name: z3950.adsabs.harvard.edu (131.142.185.23)
#~ Port: 210


#~ Protocol Version
#~ Z39.50-1992 (Version 2)
#~ Options Supported
#~ Search
#~ Present
#~ Preferred Message Size
#~ There is no restriction on message size. However, the ADS will only return a maximum of 500 records at a time.
#~ Maximum Record Size
#~ n/a
#~ ID Authentication
#~ User-id and password are not required by ADS Servers at this time.

#~ 
#~ Result Set Name
#~ No named result sets are supported.
#~ Database Names (case sensitive)
#~ (ADS Server supports searching one database at a time)
#~ Element Set Names
#~ ADS will return either brief, full, or tagged records. Database specific Element Set Names are not supported.
#~ Query
#~ Type-1 only is supported.
#~ Attribute Set ID
#~ Bib-1 only is supported.
#~ Operators Supported:
#~ AND
#~ OR
#~ AND-NOT

#~ """
#~ Simple script to search a Z39.50 target using Python
#~ and PyZ3950. 
#~ """
#~ 
#~ from PyZ3950 import zoom
#~ 
#~ 
#~ ISBNs = ['9781905017799', '9780596513986']
#~ 
#~ conn = zoom.Connection ('z3950.loc.gov', 7090)
#~ conn.databaseName = 'VOYAGER'
#~ conn.preferredRecordSyntax = 'USMARC'
#~ 
#~ for isbn in ISBNs:
    #~ query = zoom.Query ('PQF', '@attr 1=7 %s' % str(isbn))
    #~ res = conn.search (query)
    #~ for r in res:
        #~ print str(r)
#~ 
#~ conn.close ()





