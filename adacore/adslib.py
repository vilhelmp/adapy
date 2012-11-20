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

Need :  o scipy (and numpy)
        o matplotlib
        o mechanize module (standard in Python >2.6/2.7?)
        o urllib2 module (standard Python module)
        o beautiful soup/xml (xml standard in Python >2.6/2.7?)

"""

"""
ADSlib - Python Module to interact with NASA ADS

http://adswww.harvard.edu/
OR mirrors
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

29 Oct 2012 - now only uses mechanize module

02 Oct 2012 - file created!



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

def search(query, advanced=0, **kwargs):
    """
    query       :  Normal string to ADS
                   or dictionary for advanced search
    search_type :  0 - Simple
                   1 - Advanced
    
    """
    
    # advanced search query
    
    # quick search
    #~ simple_q = 'index.html'
    
    #~ types_q = [simple_q, advanced_q]
    
    #~ if search_type='adv':
        #~ type_q = advanced_q
    #~ else:
        #~ type_q = simple_q
    
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


##### For one entry
#~ 
#~ td_tags0 = entry[0].findAll('td')
#~ bibcode = td_tags0[1].findAll('input')[0]['value'].replace(u'\xa0', u' ').encode()
#~ url_abstract_page = td_tags0[1].findAll('a')[0]['href'].replace(u'\xa0', u' ').encode()
#~ ads_score = float(td_tags0[3].contents[0].encode())
#~ pubdate = td_tags0[4].contents[0].string.replace(u'\xa0', u' ').encode()
#~ pubday = pubdate[:2]
#~ pubyear = pubdate[3:]
#~ 
#~ links = dict()
#~ for link in td_tags0[5].findAll('a'):
    #~ links[link.string.encode()] = link['href'].replace(u'\xa0', u' ').encode()
#~ 
#~ td_tags1 = entry[1].findAll('td')
#~ 
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







############################################
######## THE END
############################################

### quick test to see if i works...
#~ w = 0
#~ for i in entries:
    #~ print w
    #~ _Results(i)
    #~ w += 1





### example use
# input the results
#~ res = [
    #~ _Results('Persson1', 'Water1', 100, 'various links', '2012bddaldkjf...00'), _Results('Olof', 'Ammonia', 80, 'other links', '2011nbnflkdajf..00')
    #~ ]






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


# Class that can be sorted when in a list.

#~ class _Results:
    #~ """
    #~ Internal object to represent the results
    #~ """
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
        #~ self.
    #~ def __repr__(self):
        #~ return repr([self.authors, self.title, self.score, self.links, self.bibcode])
    #~ def _returnlist_(self):
        #~ return [self.authors, self.title, self.score, self.links, self.bibcode]




### example use
# input the results
#~ res = [
    #~ _Results('Persson1', 'Water1', 100, 'various links', '2012bddaldkjf...00'), _Results('Olof', 'Ammonia', 80, 'other links', '2011nbnflkdajf..00')
    #~ ]
#~ 
#~ res.append(_Results('Jorgensen','I16293', 100,'linking', '2020320322...'))

# for advanced sorting 
# needs Python 2.6 at least
#~ from operator import itemgetter, attrgetter

# now to sort it, just use one of the keys
# score, high to low
#~ sorted(res, key=attrgetter('score'), cmp=locale.strcoll, reverse=True)

# cmp=locale.strcoll new and untested addition

# authors alphabetical order first and then by score
# i.e. sort by score if same first author
#~ sorted(res, key=attrgetter('score','authors'), cmp=locale.strcoll, reverse=True)

# now to get all the entries one by one in the sorted order and print

#~ print ('Author : {0[0].authors}\n' 
        #~ 'Title : {0[0].title}\n' 
        #~ 'Score : {0[0].score}\n' 
        #~ 'Links : {0[0].links}'.format(li)
        #~ )
#~ 
#~ #OR
#~ 
#~ print ('Author : {0[0]}\n'
        #~ 'Title : {0[1]}\n' 
        #~ 'Score : {0[2]}\n' 
        #~ 'Links : {0[3]}'.format(li[0]._returnlist_())
        #~ )









# Paper object... Noooo dictionary is better?
#~ class Paper:
    #~ def __init__(self, entry, entry_type='bfs'):
        #~ """
        #~ Hack to get all the fields, they do not have any naming of 
        #~ the td tags, or any good way to access DB
        #~ 
        #~ internal object
        #~ 
        #~ entry_type  : what type the 'entry' input is in, default is
                      #~ 'bfs', which is a list of two Beautiful Soup 
                      #~ entries. 
                      #~ available: 'bfs' or 'z39.50'
        #~ """
        #~ if entry_type=='bfs':
            #~ # first part of the result
            #~ td_tags0 = entry[0].findAll('td')
            #~ 
            #~ self.bibcode = str(td_tags0[1].findAll('input')[0]['value'])
            #~ self.url_abstract_page = str(td_tags0[1].findAll('a')[0]['href'])
            #~ self.ads_score = float(td_tags0[3].contents[0])
            #~ self.date = str(td_tags0[4].contents[0])
            #~ 
            #~ for link in td_tags0[5].findAll('a'):
                    #~ par = str(link['href'].split('=')[-1].lower()
                    #~ val = 
                    #~ self.__dict__[par] = )
            #~ 
            #~ # these short links does not all allways exist
            #~ self.url_html_article = 
            #~ self.url_pdf_article = 
            #~ self.url_arxiv_article = 
            #~ self.url_references_in = 
            #~ self.url_simbad_objects = 
            #~ self.url_also_read =
            #~ 
            #~ # second part of the result entry
            #~ self.title = str(entry[1].findAll('td')[3].contents[0])
            #~ # still in unicode
            #~ # TODO need to convert to normal UTF, not unicode
            #~ self.authors = entry[1].findAll('td')[1].contents[0]






# try 1
#~ resultform = HTMLParser.ParseResponse(result, backwards_compat=False)
#~ result.close()
#~ resultform = resultform[0]


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


###########################################
# ERROR CLASSES
class ParError(Exception):
    # input parameter error
    def __init__(self, value):
        """ Parameter Error Class
        Takes the wrong parameter as input.
        """
        self.value = value
    def __str__(self):
        """ Prints a message and the wrong parameter with value """
        s1 = '\nWrong format/number of parameters. You input:\n    '
        s2 = '\nas parameters. Check it/them.'
        return s1+str(self.value)+s2
#
###########################################
# HELP FUNCTIONS
def stylify (s='Test text', f='n', fg='r', bg='d'):
    """

    Sends back the string 'txt' with the correct foreground unicode
    color start and finish (reset color).

        Formatting style of text (f)
        f =
            "n" normal
            "b" bold
            "u" underline
            "l" blinking
            "i" inverse
        Forground color of text (fg)
        fg =
             "k" black
             "r" red
             "g" green
             "y" yellow
             "b" blue
             "m" magenta
             "c" cyan
             "a" gray
             "d" default
             "rand" random
        Background color of text (fg)
        bg =
            "k" black
            "r" red
            "g" green
            "y" yellow
            "b" blue
            "m" magenta
            "c" cyan
            "a" gray
            "d" default


    Changelog :

    *2011/10/24 added fg = "rand" for random foreground color

    """

    # needed them in this order for it to work,
    # styles, fg color, bg color
    format_and_colors = {"n_f": 0, #
                         "b_f": 1, #
                         "u_f": 4,
                         "l_f": 5,
                         "i_f": 7,
                         "k": 30,
                         "r": 31,
                         "g": 32,
                         "y": 33,
                         "b": 34,
                         "m": 35,
                         "c": 36,
                         "a": 37,
                         "d": 39,
                         "k_bg": 40,
                         "r_bg": 41,
                         "g_bg": 42,
                         "y_bg": 43,
                         "b_bg": 44,
                         "m_bg": 45,
                         "c_bg": 46,
                         "a_bg": 47,
                         "d_bg": 49}

    CSI = "\x1B["
    end = CSI+'m'

    if f == 'b' and fg =='a':
        print stylify('\n Warning : This combination of colors/styles does not work\n','b','r','d')
        raise ParError((f,fg,bg))
    bg +='_bg' # append to the list, the "_bg" ending
    f += "_f" # append "_f" to the formatting list

    if fg=="rand":
        from random import randint
        c_tmp = ["k","r","g","y","b","m","c","a","d"]
        fg = c_tmp[randint(0,len(c_tmp)-1)]
    #
    try:
        style = [format_and_colors[f.lower()],
                format_and_colors[fg.lower()],
                format_and_colors[bg.lower()]]
        style = [str(x) for x in style]
        formatted_text = CSI+';'.join(style)+'m'
        formatted_text += s + end
    except KeyError:
        raise ParError((f,fg,bg))

    return formatted_text

####################################################################
####                                                            ####
####                 RUNNING FROM CMD LINE                      ####
####                                                            ####
####################################################################
#~ if __name__ == '__main__':
    #~ # statements that you want to be executed only when the
    #~ # module is executed from the command line
    #~ # (not when importing the code by an import statement)
    #~ # wee hooo...
    #~ try:
        #~ from optparse import OptionParser as op
    #~ except (ImportError):
        #~ print stylify('ImportError',fg='r')+' Make sure you have optparse installed.'
#~ 
    #~ from sys import exit as sysexit
#~ 
    #~ #version
    #~ ver = '1.0beta'
#~ 
    #~ desc="""Script to quickly search the splatalogue compilation
    #~ Magnus Vilhelm Persson
    #~ magnusp@nbi.dk"""
    #~ usage = "Usage: %prog [options]"
    #~ epilog = """----------------------------------------------------------------------------
#~ General information :
#~ 
#~ The script uses the modules 'urllib2' and 'ClientForm' to fetch,
#~ fill in and submit the search form from Splatalogue.net. Then the
#~ results are parsed and displayed. This script can be imported into
#~ existing python code. After import the function has parameters to
#~ send the results to the user as lists or a dictionary for
#~ integration into line identification code, calculations etc.
#~ 
#~ Complete dependency list:
#~ SciPy, urllib2, ClientForm
#~ ----------------------------------------------------------------------------"""
#~ 
    #~ parser = op(usage=usage, description=desc, epilog=epilog, version="%prog " + str(ver))
#~ 
    #~ # the options permitted
    #~ parser.add_option("-f", \
        #~ dest="f",
        #~ help="frequency range given as 'F1 F2', if -w flag given, F2 is the width around F1 to look for line. Mandatory." ,
        #~ metavar="<F1> <F2>",
        #~ nargs=2,
        #~ action="store")
    #~ parser.add_option("-w",
        #~ dest="w",
        #~ help="is the f2 parameter (given in -f) the frequency width?",
        #~ default=False,
        #~ action="store_true")
    #~ parser.add_option("-u",
        #~ dest="u",
        #~ metavar="<UNIT>",
        #~ help="frequency unit, GHz or MHz.",
        #~ action="store")
    #~ parser.add_option("-l",
        #~ dest="l",
        #~ metavar="<LIST1>,<LIST2>,...",
        #~ help="molecular line list database(s) to search. \
        #~ possible values : Lovas, SLAIM, JPL, CDMS, ToyaMA, OSU, Recomb, Lisa, RFI.",
        #~ action="store")
    #~ parser.add_option("-e",
        #~ dest="e",
        #~ metavar="<FROM> <TO> <TYPE>",
        #~ nargs=3,
        #~ help="Energy range, given as 'from to type' where E_type is one of EL_cm1, EU_cm1, EL_K, EU_K.",
        #~ action="store")
    #~ parser.add_option("-t",
        #~ dest="t",
        #~ metavar="<TRANSITION>",
        #~ help="Specify transition e.g. '1-0'.",
        #~ action="store")
    #~ parser.add_option("-i",
        #~ dest="i",
        #~ metavar="<LIMIT> <UNIT>",
        #~ nargs=2,
        #~ help="Line intensity lower limit, given as 'LIMIT UNIT' where UNIT is one of CDMS_JPL, Sijmu2, Aij",
        #~ action="store")
#~ 
    #~ # time to parse
    #~ (opts, args) = parser.parse_args()
#~ 
    #~ # create the search dictionary
    #~ params = {}
#~ 
    #~ # one mandatory argument
    #~ if opts.f == None:
        #~ print stylify('\nError :',fg='r')+' No frequencies input.\n'
        #~ parser.print_help()
        #~ print ''
        #~ sysexit()
    #~ else:
        #~ f1,f2 = opts.f
        #~ if opts.w:
            #~ params['freq'] = float(f1)
            #~ params['dfreq'] = float(f2)
        #~ elif not opts.w:
            #~ params['freq'] = [float(f1),float(f2)]
    #~ if opts.u != None:
        #~ params['funit'] = opts.u
    #~ if opts.l != None:
        #~ l = (opts.l).split(',')
        #~ params['linelist'] = list(l)
    #~ if opts.e != None:
        #~ params['e_from'] = float(opts.e[0])
        #~ params['e_to'] = float(opts.e[1])
        #~ params['e_type'] = opts.e[2]
    #~ if opts.t != None:
        #~ params['transition'] = opts.t
    #~ if opts.i != None:
        #~ params['lill'] = [float(opts.i[0]), opts.i[1]]
    #~ #
    #~ params['display'] = True
    #~ params['send'] = False
    #~ # search!
    #~ splatsearch(**params)



