#! /usr/bin/env python
#
# Simple script to search Splatalogue from cli
#
# Magnus Vilhelm Persson 19.05.2011.



"""
This script queries Splatalogue.net

Returns a data structure with the requested information from Splatalogue.

TODO : page output if it is long

TODO : solve the parsing, the removal of tags...

TODO : get the export-to-file-form directly?

TODO : add so that one can enter the wavelength, just translate to frequency

"""



# import optparser
try:
    from urllib2 import urlopen
except (ImportError):
    print 'You need the module \'urllib2\''
try:
    from ClientForm import ParseResponse
except (ImportError):
    print 'You need the module \'ClientForm\' get it at http://wwwsearch.sourceforge.net/old/ClientForm/'
    print 'If you instead have the \'mechanize\' module (http://wwwsearch.sourceforge.net/mechanize/)\
    contact the programmer for implementation...'

#from BeautifulSoup import BeautifulSoup as bfs
from scipy import array, where, nan, arange

###########################################
# ERRORS

# input parameter error
class ParError(Exception):
     def __init__(self, value):
         self.value = value
     def __str__(self):
         s1 = '\nWrong format/number of parameters. You input:\n    '
         s2 = '\nas parameters. Check it/them.'
         return s1+str(self.value)+s2
#

def splatsearch(**kwargs):
    # get the form from the splatalogue search page
    response = urlopen("http://www.cv.nrao.edu/php/splat/b.php")
    forms = ParseResponse(response, backwards_compat=False)
    response.close()
    form = forms[0]
    
    arg = kwargs
    
    #now change the 'from' and 'to' controls with the frequency range requested
    # if no filename is supplied
    # the options
    #
    #### FREQUENCY
    if arg.has_key('freqs'):
        f1, f2 = arg['freqs']
        form['from'] = str(f1)
        form['to'] = str(f2)
    elif arg.has_key('freq') and type(arg['freq']) == type(1) or type(arg['freq']) == type(1.):
        if type(arg['freq']) == type(1):
            arg['freq'] = float(arg['freq'])
        if not arg.has_key('dfreq'):
            print 'Please either give a frequency interval (freq=[f1,f2])\n\
            OR a center frequency and a bandwidth (freq=f, dfreq=df)'
        elif arg.has_key('dfreq'):
            f1, f2 = arg['freq']+array([-1,1])*arg['dfreq']/2.
        form['from'] = str(f1)
        form['to'] = str(f2)
    else: # if no frequency is given, just run example
        # this is not visible when running from outside python
        # check "if __main__ ..." part
        print 'Example run... setting f1,f2 = 203.406, 203.409'
        form['from'] = '203.406'
        form['to'] = '203.409'
    #### FREQUENCY UNIT
    if arg.has_key('funit') and arg['funit'] in ['GHz', 'MHz']:
        form['frequency_units'] = [arg['funit']]
    elif not arg.has_key('funit'):
        form['frequency_units'] = ['GHz']
    else:
        print 'Only give units as either GHz or MHz, (not ghz or Ghz)'
        raise ParError(arg['funit'])
    
    #get species molecular number, ordered by mass
    # TODO : perhaps be able to search in this one
    #        either by mass or by species, text of chem formula
    # TODO : after getting it, should sort the list of dictionaries
    #        clean it up a bit
    sel_species = [i.attrs for i in form.find_control('sid[]').items]
    #
    # GET EVERYTHING!!!!
    # i.e. select everything in the form
    #Line List Display
    form.find_control("displayLovas").get().selected = True
    form.find_control("displaySLAIM").get().selected = True
    form.find_control("displayJPL").get().selected = True
    form.find_control("displayCDMS").get().selected = True
    form.find_control("displayToyaMA").get().selected = True
    form.find_control("displayOSU").get().selected = True
    form.find_control("displayRecomb").get().selected = True
    form.find_control("displayLisa").get().selected = True
    form.find_control("displayRFI").get().selected = True
    #Line Strength Display
    form.find_control("ls1").get().selected = True
    form.find_control("ls2").get().selected = True
    form.find_control("ls3").get().selected = True
    form.find_control("ls4").get().selected = True
    form.find_control("ls5").get().selected = True
    # Energy Levels
    form.find_control("el1").get().selected = True
    form.find_control("el2").get().selected = True
    form.find_control("el3").get().selected = True
    form.find_control("el4").get().selected = True
    # Miscellaneous
    form.find_control("show_unres_qn").get().selected = True
    form.find_control("show_upper_degeneracy").get().selected = True
    form.find_control("show_molecule_tag").get().selected = True
    form.find_control("show_qn_code").get().selected = True
    
    # 'click' the form
    # need to click the form first
    clicked_form = form.click()
    # then get the results page
    result = urlopen(clicked_form)
    
    #### EXPORTING RESULTS FILE
    # so what I do is that I fetch the first results page,
    # click the form/link to get all hits as a colon separated 
    # ascii table file
    #
    # get the form
    resultform = ParseResponse(result, backwards_compat=False)
    result.close()
    resultform = resultform[0]
    # set colon as dilimeter of the table (could use anything I guess)
    resultform.find_control('export_delimiter').items[2].selected =  True
    resultform_clicked = resultform.click()
    result_table = urlopen(resultform_clicked )
    data = result_table.read()
    result_table.close()
    ### PARSE the table of results ###
    # get each line (i.e. each molecule)
    lines = data.split('\n')
    # get the names of the columns
    column_names = lines[0]
    lines = lines[1:-1]
    column_names = column_names.split(':')
    hits = len(lines)
    print 'Got %d hits!' % hits
    lines = [i.split(':') for i in lines]
    #return column_names
    species, name, cfreq, cfreqerr, mfreq, mfreqerr, res_qns, ures_qns, cdmsjpl_I, \
    Smu2, Sij, log10Aij, lovasAST_I, ELcm, ELK, EUcm, EUK, u_degen, \
    mol_tag, QNr, llist = array(lines).transpose()
    # parse the columns
    # first change empty values to Nan
    cfreq[where(cfreq == '')] = 'nan'
    cfreqerr[where(cfreqerr == '')] = 'nan'
    mfreq[where(mfreq == '')] = 'nan'
    mfreqerr[where(mfreqerr == '')] = 'nan'
    ### Done parsing ###
    
    # create arrays
    cfreqerr = array(cfreqerr, dtype='float')
    cfreq = array(cfreq, dtype='float')
    mfreqerr = array(mfreqerr, dtype='float')
    mfreq = array(mfreq, dtype='float')
    N = arange(hits)+1
    #

    if arg.has_key('display') and arg['display']:
        print '{0:2} {1:10}\t{2:10}\t{3:9}\t{4:6} '.format('N', 'Form', \
        'Freq', 'FErr', 'List')
        for i in arange(len(N)):
            print '{0:2} {1:10}\t{2:10}\t{3:9}\t{4:6} '.format(N[i], \
            species[i], cfreq[i], cfreqerr[i], llist[i])
    if arg.has_key('send') and arg['send']=='dict':
        # TODO : change the output dictionary keys, a bit complicated now
        return {'N': N, 'Chem. Species': species, 'Chem. Name': name, \
        'Comp. Freq': cfreq,'Comp.Freq Err': cfreqerr, \
        'Meas. Freq': mfreq, 'Meas.Freq Err': mfreqerr,  \
        'Res. QNr': res_qns,  'URes. QNr': ures_qns, \
        'CDMS/JPL I': cdmsjpl_I, 'Smu2': Smu2, 'Sij': Sij, \
        'log10Aij': log10Aij, 'Lovas/AST I': lovasAST_I, \
        'EL (cm-1)': ELcm, 'EL (K)': ELK, 'EU (cm-1)': EUcm, \
        'EU (K)': EUK, 'Upper Degeneracy': u_degen, \
        'Molecular Tag': mol_tag, 'Quantum Nr': QNr, \
        'Line List': llist}
    elif arg.has_key('send') and arg['send'] == 'list' or arg['send']:
        print 'Sending:\n\
 Number\n Chemical Species\n Chemical Name\n Computed Frequency\n Computed Frequency Error\n\
 Measured Frequency\n Measured Frequency Error\n Resolved Quantum Numbers\n Uresolved Quantum Numbers\n\
 CDMS/JPL Intensity\n Smu**2\n Sij\n log10(Aij)\n Lovas/AST Intensity'
        return N, species, name, cfreq, cfreqerr, mfreq, mfreqerr, res_qns, ures_qns, cdmsjpl_I, \
        Smu2, Sij, log10Aij, lovasAST_I, ELcm, ELK, EUcm, EUK, u_degen, \
        mol_tag, QNr, llist
    else:
        pass
    
    
if __name__ == '__main__':
    # statements that you want to be executed only when the
    # module is executed from the command line
    # (not when importing the code by an import statement)
    # wee hooo...
    from optparse import OptionParser as op
    
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

    if options.f == None:
        print 'No frequencies input, running example with (-f) f1,f2 = 203.406,203.408'
        parser.print_usage()
        f1 = str(203.406)
        f2 = str(203.408)
    else:
        # the options
        f1,f2 = [x for x in (options.f).split(',')]
    # dont want to send anything back from the CLI
    # only display results
    splatsearch(freqs=[f1,f2], display=1, send=0)


##################################################################
##################################################################
##################################################################
########### OLD CODE with bad workarounds
# workaround, add 'submit=Search' to the request url
#result = urlopen(form.click_request_data()[0]+'&submit=Search').read()
# parse the results



####  PARSING HTML with beautiful soup
# problem because if more than 500 hits, several pages, 
# this only fetechs first page
#
#~ bfs.NESTABLE_TABLE_TAGS['td'].append('sub')
#~ bfs.NESTABLE_TABLE_TAGS['td'].append('sup')
#~ bfs.NESTABLE_TABLE_TAGS['td'].append('font')
#~ 
#~ soup = bfs(result)
#~ result.close()
#~ 
#~ def findreplace(isoup, x,p=0):
    #~ r = isoup.findAll(x)
    #~ if p:
        #~ [i.replaceWith('('+i.renderContents()+')') for i in r]
    #~ else:
        #~ [i.replaceWith(i.renderContents()) for i in r]
#~ 
#~ # these work, but only in this order (!?)
#~ findreplace(soup, 'font')
#~ findreplace(soup, 'sub')
#~ findreplace(soup, 'sup',p=1)
#~ findreplace(soup, 'a')
#~ findreplace(soup, 'br')
#~ 
#~ # find the "Found NNNN lines..." string and crab no hits
#~ # find the "Next>" link
#~ # next link from webpage
#~ # 'http://www.cv.nrao.edu/php/splat/c.php?el1=el1&el2=el2&el3=el3&el4=el4&ls1=ls1&ls2=ls2&ls3=ls3&ls4=ls4&ls5=ls5&displayRecomb=displayRecomb&displayLovas=displayLovas&displaySLAIM=displaySLAIM&displayJPL=displayJPL&displayCDMS=displayCDMS&displayToyaMA=displayToyaMA&displayOSU=displayOSU&displayLisa=displayLisa&displayRFI=displayRFI&data_version=v2.0&no_atmospheric=no_atmospheric&no_potential=no_potential&no_probable=no_probable&from=203.406&to=205&frequency_units=GHz&show_unres_qn=show_unres_qn&show_upper_degeneracy=show_upper_degeneracy&show_molecule_tag=show_molecule_tag&show_qn_code=show_qn_code&submit=1&start=500'
#~ lnk = soup.findAll('a', attrs={'class':'norm'})[-1].attrs[-1][-1]
#~ 
#~ 
#~ 
#~ results_table = soup.find('table',"results")
#~ 
#~ data = [[col.renderContents() for col in row.findAll('td')] for row in results_table.findAll('tr')]
#~ 
#~ # put all the stuff in arrays
#~ data[0][0] = 'N'
#~ col_names = data[0]
#~ data = data[1:]
#~ 
#~ N, species, name, cfreq, mfreq, res_qns, ures_qns, cdmsjpl_I, Smu2, Sij, log10Aij, lovasAST_I, ELcm, ELK, EUcm, EUK, u_degen, mol_tag, QNr, llist = array(data).transpose()
#~ 
#~ 
#~ ## this is a not so good coded part...
#~ # running i.strip('(').strip(')').split()
#~ # does not seemt to work either
#~ cfreq[where(cfreq == '&nbsp;')] = 'nan nan'
#~ mfreq[where(mfreq == '&nbsp;')] = 'nan nan'
#~ 
#~ cfreq = array([i.split(' ') for i in cfreq])
#~ mfreq = array([i.split(' ') for i in mfreq])
#~ 
#~ cfreqerr = cfreq[:,1]
#~ mfreqerr = mfreq[:,1]
#~ cfreq = cfreq[:,0]
#~ mfreq = mfreq[:,0]
#~ 
#~ 
#~ mfreqerr = [i.strip('(').strip(')') for i in  mfreqerr]
#~ cfreqerr = [i.strip('(').strip(')') for i in  cfreqerr]
#~ 
#~ cfreqerr = array(cfreqerr, dtype='float')
#~ cfreq = array(cfreq, dtype='float')
#~ mfreqerr = array(mfreqerr, dtype='float')
#~ mfreq = array(mfreq, dtype='float')
#~ 
#~ print '{0:2} {1:10}\t{2:10}\t{3:9}\t{4:6} '.format('N', 'Form', 'Freq', 'FErr', 'List')
#~ for i in arange(len(N)):
    #~ print '{0:2} {1:10}\t{2:10}\t{3:9}\t{4:6} '.format(N[i], species[i], cfreq[i], cfreqerr[i], llist[i])
#~ 
