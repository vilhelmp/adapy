#! /usr/bin/env python
#
# Simple script to search Splatalogue from cli
#
# Magnus Vilhelm Persson 19.05.2011.



"""
Module with a function to  query Splatalogue.net

TODO : page output if it is long

TODO : solve the parsing, the removal of tags...

TODO : get the export-to-file-form directly?

TODO : add so that one can enter the wavelength, just translate to frequency

TODO : Figure out prettier printing of linelists... colors, web-adresses?

"""
###########################################
# MAIN FUNCTION
def splatsearch(**kwargs):
    """
    This script queries Splatalogue.net
    
    Returns a data structure with the requested information from Splatalogue.
    
    Usage:
        splatsearch(**kwargs)
        
    Options/Arguments:
        Frequency
            freq = [f1, f2] 
                frequency interval to query
            freq = f and dfreq = df
                center frequency (f) and bandwidth (df)
            funit = 'GHz'/'MHz'
                freqyency unit
            Note: low ranked plans to implement searching by
                  wavelength/wavenumber (m,cm,mm/m-1) exists
        Line list
            linelist = ['lovas', 'slaim', 'jpl', 'cdms', 'toyama', 'osu',\
         'recomb', 'lisa', 'rfi']
            or linelists = 'all'
    """
    ### import dependencies
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
    from string import lower
    
    # get the form from the splatalogue search page
    response = urlopen("http://www.cv.nrao.edu/php/splat/b.php")
    forms = ParseResponse(response, backwards_compat=False)
    response.close()
    form = forms[0]
    arg = kwargs
    if arg.has_key('dbg_snd_form'):
        return form
    ####################################################################
    ####                                                            ####
    ####                        PARSE INPUT                         ####
    ####                                                            ####
    ####################################################################
    #### FREQUENCY
    if arg.has_key('freq'):
        if type(arg['freq']) == type([1,2]):
            if len(arg['freq']) == 1:
                raise ParError(arg['freq'])
            f1, f2 = arg['freq']
            form['from'] = str(f1)
            form['to'] = str(f2)
            print 'Got frequency interval from %f to %f' % (f1, f2)
        elif type(arg['freq']) == type(1) or type(arg['freq']) == type(1.):
            if type(arg['freq']) == type(1):
                arg['freq'] = float(arg['freq'])
            if not arg.has_key('dfreq'):
                print 'Please either give a frequency interval (freq=[f1,f2])\n\
                OR a center frequency and a bandwidth (freq=f, dfreq=df)'
                raise ParError('freq='+str(arg['freq'])+' and dfreq=None')
            elif arg.has_key('dfreq'):
                f1, f2 = arg['freq']+array([-1,1])*arg['dfreq']/2.
            else:
                raise ParError(arg['dfreq'])
            form['from'] = str(f1)
            form['to'] = str(f2)
    elif not arg.has_key('freq') and arg.has_key('dfreq'):
        print 'Only delta-frequency (dfreq) given, no frequency to look for'
        raise ParError('freq=None and dfreq='+str(arg['dfreq']))
    else: 
        # if no frequency is given, just run example
        # this is not visible when running from outside python
        # check "if __main__ ..." part
        print colorify('Example run... setting f1,f2 = 203.406, 203.409 GHz',c='m')
        form['from'] = '203.406'
        form['to'] = '203.409'
    #### FREQUENCY UNIT
    if arg.has_key('funit'):
        if lower(arg['funit']) in ['ghz', 'mhz']:
            form['frequency_units'] = [arg['funit']]
        else:
            print 'Onlye '
    elif not arg.has_key('funit'):
        print colorify('No frequency units given, assuming GHz.',c='b')
        form['frequency_units'] = ['GHz']
    #### MOLECULAR SPECIES
    #
    #get species molecular number, ordered by mass
    # TODO : perhaps be able to search in this one
    #        either by mass or by species, text of chem formula
    # TODO : after getting it, should sort the list of dictionaries
    #        clean it up a bit
    sel_species = [i.attrs for i in form.find_control('sid[]').items]
    #### LINE LIST
    # define a reference list of names
    mylinelist = ['lovas', 'slaim', 'jpl', 'cdms', 'toyama', 'osu', \
    'recomb', 'lisa', 'rfi']
    # list of strings with the format that the search form wants
    formcontrol_linelist = ["displayLovas", "displaySLAIM", \
    "displayJPL", "displayCDMS", "displayToyaMA", "displayOSU", \
    "displayRecomb", "displayLisa", "displayRFI"]
    if arg.has_key('linelist'):
        if type(arg['linelist'])==type('string'):
            # if linelist is given as linelist='all'
            if lower(arg['linelist']) == 'all':
                # if we want to set all, just copy mylinelist
                arg['linelist'] = mylinelist
            else:
                print 'Linelist input not understood'
                raise ParError(linelist)
        elif type(arg['linelist'])==type(['list']):
            # get all values to lower case, to accept capitals
            arg['linelist'] = [lower(x) for x in arg['linelist']]
        else:
            print 'Linelist input not understood'
            raise ParError(linelist)
    else:
        # if none given, search with all
        arg['linelist'] = mylinelist

    # now set the linelist search form
    # check for every linelist, if it exists in the input linelist
    for i,j in zip(mylinelist, formcontrol_linelist):
        if i in arg['linelist']:
            form.find_control(j).get().selected = True
        else:
            form.find_control(j).get().selected = False
    # ['Lovas', 'SLAIM', 'JPL', 'CDMS', 'ToyaMA', 'OSU', \
    #'Recomb', 'Lisa', 'RFI']
    # Figure out prettier printing here... colors, web-adresses?
    print 'Using linelist(s): '+', '.join(arg['linelist'])+' '
    ### Energy Range
    #<TextControl(energy_range_from=)>
    #<TextControl(energy_range_to=)>
    #<RadioControl(energy_range_type=[el_cm1, eu_cm1, el_k, eu_k])>
    # range
    
    # type
    
    ### Line Intensity Lower Limits
    #<RadioControl(lill=[on, on, on])>          # only one can be on
    #<TextControl(lill_cdms_jpl=) (disabled)>   # enter energy limit (log)
    #<TextControl(lill_sijmu2=) (disabled)>     # debye**2
    #<TextControl(lill_aij=) (disabled)>        # log

    
    
    #### Line Strength Display
    form.find_control("ls1").get().selected = True
    form.find_control("ls2").get().selected = True
    form.find_control("ls3").get().selected = True
    form.find_control("ls4").get().selected = True
    form.find_control("ls5").get().selected = True
    #### Energy Levels
    form.find_control("el1").get().selected = True
    form.find_control("el2").get().selected = True
    form.find_control("el3").get().selected = True
    form.find_control("el4").get().selected = True
    #### Miscellaneous
    form.find_control("show_unres_qn").get().selected = True
    form.find_control("show_upper_degeneracy").get().selected = True
    form.find_control("show_molecule_tag").get().selected = True
    form.find_control("show_qn_code").get().selected = True
    ####################################################################
    ####                                                            ####
    ####                        GET RESULTS                         ####
    ####                                                            ####
    ####################################################################
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
    result_table = urlopen(resultform_clicked)
    data = result_table.read()
    result_table.close()
    ####################################################################
    ####                                                            ####
    ####                        PARSE RESULT                        ####
    ####                                                            ####
    ####################################################################
    # get each line (i.e. each molecule)
    lines = data.split('\n')
    # get the names of the columns
    column_names = lines[0]
    lines = lines[1:-1]
    column_names = column_names.split(':')
    hits = len(lines)
    if hits == 0:
        print 'No lines found!'
        return None
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
    # create arrays
    cfreqerr = array(cfreqerr, dtype='float')
    cfreq = array(cfreq, dtype='float')
    mfreqerr = array(mfreqerr, dtype='float')
    mfreq = array(mfreq, dtype='float')
    N = arange(hits)+1
    ####################################################################
    ####                                                            ####
    ####                     DISPLAY RESULTS                        ####
    ####                                                            ####
    ####################################################################
    if arg.has_key('display') and arg['display']:
        print '{0:2} {1:15}\t{2:10}\t{3:9}\t{4:6} '.format('N', 'Form', \
        'Freq', 'FErr', 'List')
        for i in arange(len(N)):
            print '{0:2} {1:15}\t{2:10}\t{3:9}\t{4:6} '.format(N[i], \
            species[i], cfreq[i], cfreqerr[i], llist[i])
    if arg.has_key('send'):
        if arg['send']=='dict':
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
        if arg['send'] == 'list' or arg['send']:
            if arg.has_key('silen'):
                if not arg['silent']:
                    print 'Sending:\n\
             Number\n Chemical Species\n Chemical Name\n Computed Frequency\n Computed Frequency Error\n\
             Measured Frequency\n Measured Frequency Error\n Resolved Quantum Numbers\n Uresolved Quantum Numbers\n\
             CDMS/JPL Intensity\n Smu**2\n Sij\n log10(Aij)\n Lovas/AST Intensity'
                elif arg['silent']:
                    pass
            return N, species, name, cfreq, cfreqerr, mfreq, mfreqerr, res_qns, ures_qns, cdmsjpl_I, \
            Smu2, Sij, log10Aij, lovasAST_I, ELcm, ELK, EUcm, EUK, u_degen, \
            mol_tag, QNr, llist
    else:
        pass

def get_mol_species():
    """Function to get the molecular species...
    only started...
    Best to take the list/dictionary as input (more like a matching function)?
    """


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
def colorify (txt, c='r'):
    """ 
    
    Sends back the string 'txt' with the correct foreground unicode 
    color start and finish (reset color).

    Current avaliable colors:
        'r' : red
        'g' : green
        'y' : yellow
        'b' : blue
        'm' : magenta
        'c' : cyan
        'w' : white
    """
    CSI = "\x1B["

    if c=='r':
        # sets, red text color (31)
        start =CSI+'31m'
    elif c=='g':
        # sets, green text color (32)
        start =CSI+'32m'
    elif c=='y':
        # sets, yellow text color (32)
        start =CSI+'33m'
    elif c=='b':
        # sets, blue text color (34)
        start =CSI+'34m'
    elif c=='m':
        # sets, magenta text color (34)
        start =CSI+'35m'
    elif c=='c':
        # sets, cyan text color (34)
        start =CSI+'36m'
    elif c=='w':
        # sets, white text color (34)
        start =CSI+'37m'
    #
    end = CSI+'m'
    return start+txt+end

####################################################################
####                                                            ####
####                 RUNNING FROM CMD LINE                      ####
####                                                            ####
####################################################################
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
    splatsearch(freq=[f1,f2], display=1, send=0)
