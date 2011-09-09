#! /usr/bin/env python
#
#       Simple script to search Splatalogue from cli
#
#       Copyright 2011 Magnus Persson <magnusp@nbi.dk>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       ver. 1.0beta
#
#       Magnus Vilhelm Persson 19.05.2011.



"""
Module with a function to query Splatalogue.net

TODO : page output if it is long, i.e. input a pause and clear screen

TODO : control display-output more

TODO : solve the parsing, the removal of tags...

TODO : get the export-to-file directly from splatalogue?

TODO : add so that one can enter the wavelength/wavenumber in m, cm, m-1, cm-1
        just translate to frequency after input

TODO : figure out prettier printing of linelists... colors, web-adresses?
        so when linelists choosen, print e.g.
        "Line lists choosen:
        JPL         http://jpl.nasa.gov/molecules/etc
        CDMS        http://cologne.de/CDMS/whatever"

TODO : Change dictionary key-names in return-dictionary

TODO : Implement new stylify help function

"""
###########################################
# MAIN FUNCTION
def splatsearch(**arg):
    """
    This script queries Splatalogue.net from the command line
    
    Returns a data structure or displays the requested information 
    from Splatalogue.
    
    Usage:
        splatsearch(**kwargs)
        
    Needs internet connection to work (duh!).
    
    Keyword arguments (kwargs)
        
        Frequency
          o parameter name : freq
                [f1, f2]  or f
                frequency interval to query
          o parameter name : dfreq
                if freq only f, a dfreq must be given
                center frequency (f) and bandwidth (df)
          o parameter name : funit
                possible values : 'GHz' or 'MHz'
                freqyency unit
            Note: low ranked plans to implement searching by
                  wavelength/wavenumber (m,cm,mm/m-1) exists
        
        Line list
          o parameter name : 'linelist'
                type : list
                a list of strings ['item1', 'item2']
                conainting the line list catalogs that you want 
                to search in. Possible entries:
                ['Lovas', 'SLAIM', 'JPL', 'CDMS', 'ToyaMa',\
                            'OSU', 'Recomb', 'Lisa', 'RFI']
                or linelists = 'all' for all of them
            
        Energy range
          o parameter name : e_from
                type : int/float
          o parameter name : e_to
                type : int/float
          o parameter name : e_type
                type :  string
                one of ['el_cm1', 'eu_cm1', 'el_k', 'eu_k']
                unit and type, in cm-1 or K
                if not given, defaults to eu_k
        
        Line Intensity Lower Limit
          o parameter name : lill
                type : list
                list in the form [value, 'type']
                possible types:
                    'cdms_jpl' : CDMS/JPL Intensity given in log10
                     'Sij mu^2' : Sij mu^2 in Debye^2
                         'Aij' : Aij in log10
                example lill = [-5, 'cdms_jpl']
                for CDMS/JPL Intensity of 10^(-5)
        Transition
          o parameter name : transition
                type : string
                example transition = '1-0'
        
        OUTPUT
        Display output
        Frequency
        
        
        
        
        freq measured if exist otherwise computed
        
    """
    #~ arg = kwargs
    #
    if arg.has_key('display'):
        if arg['display'] == 1:
            def print_greeting():
                text_splat_colors = [colorify(i,c='rand') for i in 'SPLATSEARCH']
                text_splat_colors = ''.join(text_splat_colors)
                print "\n    "+"*"*40
                print "    *"+" "*38+"*"
                print "    *\t           "+text_splat_colors+"             *"
                print "    *\t  Splatalogue.net search script    *"
                print "    *\t    Magnus Vilhelm Persson         *"
                print "    *\t        magnusp@nbi.dk             *"
                print "    *"+" "*38+"*"
                print "    "+"*"*40+"\n"
            if not arg.has_key('send'):
                print_greeting()
            elif arg.has_key('send') and not arg['send']:
                print_greeting()
    ### import dependencies
    try:
        from urllib2 import urlopen
    except (ImportError):
        print 'You need the module \'urllib2\''
    try:
        from ClientForm import ParseResponse
    except (ImportError):
        print 'You need the module \'ClientForm\' get it at http://wwwsearch.sourceforge.net/old/ClientForm/'
        print 'If you instead have the newer \'mechanize\' module (http://wwwsearch.sourceforge.net/mechanize/) contact the programmer for implementation...'
    
    #from BeautifulSoup import BeautifulSoup as bfs
    from scipy import array, where, nan, arange
    from string import lower, upper
    
    # get the form from the splatalogue search page
    response = urlopen("http://www.cv.nrao.edu/php/splat/b.php")
    forms = ParseResponse(response, backwards_compat=False)
    response.close()
    form = forms[0]
    
    if arg.has_key('dbg_snd_form'): # for test purposes
        return form
    ####################################################################
    ####                                                            ####
    ####                        PARSE INPUT                         ####
    ####                                                            ####
    ####################################################################
    #
    #### FREQUENCY
    #
    #
    #           No frequency given, what then?
    #           a looooot of hits returned, perhaps pause and ask if user wants to continue??
    #
    #
    #
    if arg.has_key('freq'):
        if type(arg['freq']) == type([1,2]):
            if len(arg['freq']) == 1:
                raise ParError(arg['freq'])
            f1, f2 = arg['freq']
            form['from'] = str(f1)
            form['to'] = str(f2)
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
    elif not arg.has_key('freq') and not arg.has_key('dfreq') and len(arg.keys()) != 0:
        # no frequency given, but other parameters
        tmp = str(raw_input('No frequency limits given, continue? Press Enter to continue, Ctrl+C to abort.'))
        f1 = ''
        f2 = ''
    else: 
        # if no frequency is given, just run example
        # this is not visible when running from outside python
        # check "if __main__ ..." part
        print colorify('Example run... setting f1,f2 = 203.406, 203.409 GHz',c='m')
        form['from'] = '203.406'
        form['to'] = '203.409'
    #
    #### FREQUENCY UNIT
    #
    if arg.has_key('funit'):
        if lower(arg['funit']) in ['ghz', 'mhz']:
            form['frequency_units'] = [arg['funit']]
        else:
            print 'Allowed frequency units : \'GHz\' or \'MHz\''
    elif not arg.has_key('funit'):
        arg['funit'] = 'GHz'
        form['frequency_units'] = ['GHz']
    #
    #### MOLECULAR SPECIES
    #
    #Get species molecular number, ordered by mass
    # TODO : perhaps be able to search in this one
    #        either by mass or by species, text of chem formula
    # TODO : after getting it, should sort the list of dictionaries
    #        clean it up a bit
    # get the avaliable species from the form
    sel_species = [i.attrs for i in form.find_control('sid[]').items]
    #
    #### LINE LIST
    #
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
    # Figure out prettier printing here... 
    #    web-adresses?
    #
    ### Energy Range
    #
    # form['energy_range_from/to'] is a text field in the form
    # while it is called e_from/to in the function
    #
    if arg.has_key('e_from') or arg.has_key('e_to'):
        e_type_ref = ['el_cm1', 'eu_cm1', 'el_k', 'eu_k'] 
        # check that unit is given, and correct
        # or set default (eu_k)
        
        if arg.has_key('e_from'):
            form['energy_range_from'] = str(arg['e_from'])
        if arg.has_key('e_to'):
            form['energy_range_to'] = str(arg['e_to'])
        if arg.has_key('e_from') or arg.has_key('e_to'):
            if arg.has_key('e_type'):
                if lower(arg['e_type']) in e_type_ref:
                    pass
                else:
                    print 'Energy range type keyword \'e_type\' malformed.'
                    raise ParError(arg['e_type'])
                e_type_default = 0
            else:
                e_type_default = 1
                arg['e_type'] = 'eu_k'
            # now set the radio button to the correct value
            form.find_control('energy_range_type').toggle(arg['e_type'])
        if not arg.has_key('e_from') and not arg.has_key('e_to') and arg.has_key('e_type'):
            print 'You gave the Enery range type keyword, but no energy range...'
            raise ParError(arg['e_type'])
    #
    ### Specify Transition
    #    
    if arg.has_key('transition'):
        form['tran'] = str(arg['transition'])
    #
    ### Line Intensity Lower Limits
    if arg.has_key('lill'):
        if lower(arg['lill'][1]) == 'cdms_jpl':
            form.find_control('lill_cdms_jpl').disabled = False
            form['lill_cdms_jpl'] = str(arg['lill'][0])
        elif lower(arg['lill'][1]) == 'sijmu2':
            form.find_control('lill_sijmu2').disabled = False
            form['lill_sijmu2'] = str(arg['lill'][0])
        elif lower(arg['lill'][1]) == 'aij':
            form.find_control('lill_aij').disabled = False
            form['lill_aij'] = str(arg['lill'][0])
    #
    ### FREQUENCY ERROR LIMIT
    #
    
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
    ####               DISPLAY SEARCH PARAMETERS                    ####
    ####                                                            ####
    ####################################################################
    print colorify('** SEARCH PARAMETERS **',c='b')
    if arg.has_key('freq'):
        print colorify('Frequency range :',c='g')+' '+str(f1)+' - '+str(f2)
        print colorify('Frequency unit \t:',c='g')+' '+arg['funit']
    else:
        print 'No frequency range specified'
    print colorify('Line list(s) \t:',c='g')+' '+', '.join(arg['linelist'])
    if arg.has_key('e_from') or arg.has_key('e_to'):
        if arg.has_key('e_from') and not arg.has_key('e_to'):
            print colorify('Energy range \t:',c='g')+'from '+str(arg['e_from'])+'( Type : %s)' % str([arg['e_type'],'def (EU(K))'][e_type_default])
        elif not arg.has_key('e_from') and arg.has_key('e_to'):
            print colorify('Energy range \t:',c='g')+'to '+str(arg['e_to'])+'( Type : %s)' % str([arg['e_type'],'def (EU(K))'][e_type_default])
        else:
            #print colorify('Energy range \t:',c='g')+upper(arg['e_type'][:2])+' from '+str(arg['e_from'])+' to '+str(arg['e_to'])+' 'upper(arg['e_type'][3:])+'( Type : %s)' % str([arg['e_type'],'yes'][e_type_default])
            print colorify('Energy range \t:',c='g')+' %s from %s to %s %s (Type : %s)' % (upper(arg['e_type'][:2]),str(arg['e_from']),str(arg['e_to']),upper(arg['e_type'][3:]), str([arg['e_type'],'yes'][e_type_default]))
    if arg.has_key('lill'):
        if lower(arg['lill'][1]) == 'cdms_jpl':
            print colorify('Line lower lim \t:',c='g')+' 1E('+str(arg['lill'][0])+') - CDMS/JPL Intensity'
        elif lower(arg['lill'][1]) == 'sijmu2':
            print colorify('Line lower lim \t:',c='g')+' '+str(arg['lill'][0])+' Debye^2 - Sijmu^2'
        elif lower(arg['lill'][1]) == 'aij':
            print colorify('Line lower lim \t:',c='g')+' 1E('+str(arg['lill'][0])+') - Aij'
    if arg.has_key('transition'):
        print colorify('Transition \t:',c='g')+' '+arg['transition']
    #~ if arg.has_key(''):
    print ''
    
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
        print '\nNo lines found!'
        return None
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
    # create global frequency array, and a 
    # array telling if it is measured or computed
    # empty arrays
    from scipy import zeros
    freq = zeros(cfreq.shape)
    freqerr  = zeros(cfreqerr.shape)
    freqtype = []
    # use measured frequency if exists
    # otherwise use computed
    for i in arange(hits):
        if str(mfreq[i]) == 'nan':
            freq[i] = cfreq[i]
            freqerr[i] = cfreqerr[i]
            freqtype.append('C')
        else:
            freq[i] = mfreq[i]
            freqerr[i] = mfreqerr[i]
            freqtype.append('M')
    N = arange(hits)+1
    ####################################################################
    ####                                                            ####
    ####                     DISPLAY RESULTS                        ####
    ####                                                            ####
    ####################################################################
    if arg.has_key('display') and arg['display']:
        print colorify('** RESULTS **',c='b')
        print 'Got %s hits!' % colorify(str(hits),c='r')
        print colorify('{0:2} {1:13}\t{2:10}\t{3:9}\t{4:10}\t{5:3}\t{6:6}',c='m').format('N', 'Form', \
        'Freq', 'Smu^2', 'EU(K)','C/M', 'List')
        for i in arange(hits):
            print '{0:2} {1:13}\t{2:10}\t{3:9}\t{4:10}\t{5:3}\t{6:6} '.format(N[i], \
            species[i], freq[i], Smu2[i], EUK[i], freqtype[i], llist[i])
    if arg.has_key('send'):
        if arg['send']=='dict':
            # TODO : change the output dictionary keys, a bit complicated now
            return {'N': N, 'Chem. Species': species, 'Chem. Name': name, \
            'Comp. Freq': cfreq,'Comp.Freq Err': cfreqerr, \
            'Meas. Freq': mfreq, 'Meas.Freq Err': mfreqerr,  \
            'Freq': freq, 'Freq Err': freqerr,  \
            'FreqType': freqtype, \
            'Res. QNr': res_qns,  'URes. QNr': ures_qns, \
            'CDMS/JPL I': cdmsjpl_I, 'Smu2': Smu2, 'Sij': Sij, \
            'log10Aij': log10Aij, 'Lovas/AST I': lovasAST_I, \
            'EL (cm-1)': ELcm, 'EL (K)': ELK, 'EU (cm-1)': EUcm, \
            'EU (K)': EUK, 'Upper Degeneracy': u_degen, \
            'Molecular Tag': mol_tag, 'Quantum Nr': QNr, \
            'Line List': llist}
        if arg['send'] == 'list' or arg['send']:
            if arg.has_key('silent'):
                if not arg['silent']:
                    print 'Sending:\n\
             Number\n Chemical Species\n Chemical Name\n Computed Frequency\n Computed Frequency Error\n\
             Measured Frequency\n Measured Frequency Error\n Resolved Quantum Numbers\n Uresolved Quantum Numbers\n\
             CDMS/JPL Intensity\n Smu**2\n Sij\n log10(Aij)\n Lovas/AST Intensity'
                elif arg['silent']:
                    pass
            return N, species, name, freq, freqerr, freqtype, cfreq, cfreqerr, mfreq, mfreqerr, res_qns, ures_qns, cdmsjpl_I, \
            Smu2, Sij, log10Aij, lovasAST_I, ELcm, ELK, EUcm, EUK, u_degen, \
            mol_tag, QNr, llist
    else:
        pass

def get_mol_species():
    """Function to get the molecular species...
    only started...
    Best to take the list/dictionary as input (more like a matching function)?
    Should return a structure compilant with the form input
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
    elif c=='rand':
        # white not included
        from scipy import rand
        #colors = ['r','g','y','b','m','c']
        codes = ['31m','32m','33m','34m','35m','36m']
        i = int(round((rand()*len(codes)))-1)
        start = CSI+codes[i]
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
    try:
        from optparse import OptionParser as op
    except (ImportError):
        print colorify('ImportError',c='r')+' Make sure you have optparse installed.'
    
    from sys import exit as sysexit
    
    #version
    ver = '1.0beta'
    
    desc="""Script to quickly search the splatalogue compilation                      
    Magnus Vilhelm Persson                                                            
    magnusp@nbi.dk""" 
    usage = "Usage: %prog [options]" 
    epilog = """----------------------------------------------------------------------------
General information :                                                       

The script uses the modules 'urllib2' and 'ClientForm' to fetch,
fill in and submit the search form from Splatalogue.net. Then the
results are parsed and displayed. This script can be imported into
existing python code. After import the function has parameters to
send the results to the user as lists or a dictionary for
integration into line identification code, calculations etc.                        

Complete dependency list:                                                   
SciPy, urllib2, ClientForm                                        
----------------------------------------------------------------------------"""
    
    parser = op(usage=usage, description=desc, epilog=epilog, version="%prog " + str(ver))
    
    # the options permitted
    parser.add_option("-f", \
        dest="f", 
        help="frequency range given as 'F1 F2', if -w flag given, F2 is the width around F1 to look for line. Mandatory." , 
        metavar="<F1> <F2>",
        nargs=2,
        action="store")
    parser.add_option("-w", 
        dest="w", 
        help="is the f2 parameter (given in -f) the frequency width?",
        default=False,
        action="store_true")
    parser.add_option("-u", 
        dest="u",
        metavar="<UNIT>",
        help="frequency unit, GHz or MHz.",
        action="store")
    parser.add_option("-l", 
        dest="l", 
        metavar="<LIST1>,<LIST2>,...",
        help="molecular line list database(s) to search. \
        possible values : Lovas, SLAIM, JPL, CDMS, ToyaMA, OSU, Recomb, Lisa, RFI.",
        action="store")
    parser.add_option("-e", 
        dest="e", 
        metavar="<FROM> <TO> <TYPE>",
        nargs=3,
        help="Energy range, given as 'from to type' where E_type is one of EL_cm1, EU_cm1, EL_K, EU_K.",
        action="store")
    parser.add_option("-t", 
        dest="t",
        metavar="<TRANSITION>",
        help="Specify transition e.g. '1-0'.",
        action="store")
    parser.add_option("-i", 
        dest="i", 
        metavar="<LIMIT> <UNIT>",
        nargs=2,
        help="Line intensity lower limit, given as 'LIMIT UNIT' where UNIT is one of CDMS_JPL, Sijmu2, Aij",
        action="store")
    
    # time to parse
    (opts, args) = parser.parse_args()
    
    # create the search dictionary
    params = {}
    
    # one mandatory argument
    if opts.f == None:
        print colorify('\nError :',c='r')+' No frequencies input.\n'
        parser.print_help()
        print ''
        sysexit()
    else:
        f1,f2 = opts.f
        if opts.w:
            params['freq'] = float(f1)
            params['dfreq'] = float(f2)
        elif not opts.w:
            params['freq'] = [float(f1),float(f2)]
    if opts.u != None:
        params['funit'] = opts.u
    if opts.l != None:
        l = (opts.l).split(',')
        params['linelist'] = list(l)
    if opts.e != None:
        params['e_from'] = float(opts.e[0])
        params['e_to'] = float(opts.e[1])
        params['e_type'] = opts.e[2]
    if opts.t != None:
        params['transition'] = opts.t
    if opts.i != None:
        params['lill'] = [float(opts.i[0]), opts.i[1]]
    #
    params['display'] = True
    params['send'] = False
    # search!
    splatsearch(**params)
    
