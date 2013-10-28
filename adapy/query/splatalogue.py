
SPLAT_FORM_URL = "http://www.cv.nrao.edu/php/splat/c_export.php"
HIT_LIMIT = 2000
SPLATALOGUE_TIMEOUT = 30

"""
NB this script has a export limit of 2000 hits.
Change HIT_LIMIT to accomodate your needs.
"""

__all__ = ['search']

import numpy as np
try:
    from astropy.table import Table
    use_astropy = True
except (ImportError):
    use_astropy = False

import numpy as _np

"""

TODO : improve the parsing of the astropy.table output

TODO : get partition values/function for certain molecule?

TODO : simple escape probability calculation for transitions
        of molecule to estimate strengths?

TODO : be able to search on molecular species
        either just parse the output (grep type)
                          or
        actually search for specific molecule (sid[])

TODO : create pretty printing for on screen results
        in the bin directory


TODO : add help docs
        explain how to use it

        
TODO : clean up the QN strings in the results


"""

#~ The urllib2 module has been split across several modules in Python 3.0
#~ named urllib.request and urllib.error. The 2to3 tool will automatically adapt imports when converting your sources to 3

# used for Python 2: 
# urllib.urlencode
# urllib2.Request
# urllib2.urlopen

try:
    # For Python 3.0 and later
    from urllib.request import urlopen, Request
    from urllib.parse import urlencode
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib import urlencode
    from urllib2 import urlopen, Request
    
###############################################################################

def search( freq = [203.4, 203.42],
            fwidth = None,
            funit = 'GHz',
            linelist = ['lovas', 'slaim', 'jpl', 'cdms', 'toyama', 'osu', 'recomb', 'lisa', 'rfi'],
            efrom = None,
            eto = None,
            eunit = None,    # 'el_cm1', 'eu_cm1', 'el_k', 'eu_k'
            trans = None,
            lill = None,    # line intensity lower limits, 'cdms_jpl', 'sijmu2', 'aij'
             **settings):

        """
            
            settings:
            ----------

            version : 2.0, 1.0, all



            otype   : output type

                        astropy.table - an astropy.table table
            
        """

        # prepare some of the input
        # lowercase insensitive input.
        if type(linelist) == type([]):
            linelist = [i.lower() for i in linelist]
        else:
            lineliest = [linelist.lower()]
        
        #~ _check_input(freq, fwidth, funit, linelist, efrom, eto, eunit,
            #~ transition, lill,**settings)

        # parameters list, to be urlencoded later        
        parameters = []
        
        parameters.extend( _parameters_preamble() )
        parameters.extend( _parse_linelist( linelist ) )
        parameters.extend( _parse_settings( settings ) )
        
        parameters.extend( _parse_frequency(freq, fwidth, funit) )
        parameters.extend( _parse_linelist( linelist ) )
        
        if ((efrom is not None) or (eto is not None)):
            parameters.extend( _parse_erange( efrom, eto, eunit ) )
        if (trans is not None):
            parameters.extend( _parse_transition( trans ) )
        if (lill is not None):
            parameters.extend( _parse_lill( lill ) )
        
        parameters.extend( _parameters_ending() )
        results = _get_results(parameters)

        if settings.has_key('otype'):
            results = _parse_results(results, settings['otype'])
        else:
            results = _parse_results(results, 'astropy.table')
        return results
        
def _parameters_preamble():
    """
        Set the default display parameters, i.e. display
        everything in the result table


        #  Energy level display (triggered)
        # 1 : Elower (cm-1)
        # 2 : Elower (K)
        # 3 : Eupper (cm-1)
        # 4 : Eupper (K)
        el1=el1
        el2=el2
        el3=el3
        el4=el4

        # Line strength display (triggered)
        # 1 : CDMS/JPL Intensity
        # 2 : Sij mu2
        # 3 : Sij
        # 4 : Aij
        # 5 : Lovas/AST
        ls1=ls1
        ls2=ls2
        ls3=ls3
        ls4=ls4
        ls5=ls5

        # Display Unresolved quantum numbers (triggered)
        # always on
        #~ show_unres_qn=show_unres_qn

        # Show upper degeneracy (triggered)
        # always on
        #~ show_upper_degeneracy=show_upper_degeneracy

        # Display Molecule Tag (triggered)
        # always on
        #~ show_molecule_tag=show_molecule_tag

        # No HFS Display (triggered)
        # not included 
        #~ noHFS=noHFS

        # Display HFS Intensity (triggered)
        # always on
        #~ displayHFS=displayHFS

        # Display Quantum Number Code (triggered)
        # always on
        #~ show_qn_code=show_qn_code

        # Display Lab Ref (triggered)
        # always off
        #~ show_lovas_labref=show_lovas_labref

        # Display Obs Ref (triggered)
        # always off
        #~ show_lovas_obsref=show_lovas_obsref

        # Display Ordered Frequency ONLY (triggered)
        #~ show_orderedfreq_only=show_orderedfreq_only

        # Display NRAO Recommended Frequencies (triggered)
        #~ show_nrao_recommended=show_nrao_recommended

        SUMMARY :
        ALWAYS ON-----------------------------------
        Display HFS Intensity  
        Display Unresolved Quantum Numbers  
        Display Upper State Degeneracy 
        Display Molecule Tag  
        Display Quantum Number Code  
        Display NRAO Recommended Frequencies
        E_levels (all)
        Line Strength Display (all)
        Display Ordered Frequency ONLY (only one frequency to parse)
        --------------------------------------------

    """
    returnlist = [
        ('submit', 'Search'),
        ('ls1','ls1'),
        ('ls2','ls2'),
        ('ls3','ls3'),
        ('ls4','ls4'),
        ('ls5','ls5'),
        ('el1', 'el1'),
        ('el2', 'el2'),
        ('el3', 'el3'),
        ('el4', 'el4'),
        ('show_unres_qn', 'show_unres_qn'),
        ('show_upper_degeneracy', 'show_upper_degeneracy'),
        ('show_molecule_tag', 'show_molecule_tag'),
        ('displayHFS', 'displayHFS'),
        ('show_qn_code', 'show_qn_code'),
        #('show_lovas_labref', 'how_lovas_labref'),          # Always OFF
        #('show_lovas_obsref', 'show_lovas_obsref'),         # Always OFF
        #~ ('show_orderedfreq_only', 'show_orderedfreq_only'),
        ('show_nrao_recommended', 'show_nrao_recommended')
        ]
    return returnlist

def _parse_linelist(linelist):
    """
    Only search the requested line lists.
    
    # line list (triggered)
    # def all on
    displayRecomb=displayRecomb 
    displayLovas=displayLovas
    displaySLAIM=displaySLAIM
    displayJPL=displayJPL
    displayCDMS=displayCDMS
    displayToyaMA=displayToyaMA
    displayOSU=displayOSU
    displayLisa=displayLisa
    displayRFI=displayRFI
    """
    returnlist = []
    if 'lovas' in linelist:
        returnlist.append(('displayLovas'  ,'displayLovas'))
    if 'slaim' in linelist:
        returnlist.append(('displaySLAIM'  ,'displaySLAIM'))
    if 'jpl' in linelist:
        returnlist.append(('displayJPL'    ,'displayJPL'))
    if 'cdms' in linelist:
        returnlist.append(('displayCDMS'   ,'displayCDMS'))
    if 'toyama' in linelist:
        returnlist.append(('displayToyaMA' ,'displayToyaMA'))
    if 'osu' in linelist:
        returnlist.append(('displayOSU'    ,'displayOSU'))
    if 'recomb' in linelist:
        returnlist.append(('displayRecomb' ,'displayRecomb'))
    if 'lisa' in linelist:
        returnlist.append(('displayLisa'   ,'displayLisa'))
    if 'rfi' in linelist:
        returnlist.append(('displayRFI'    ,'displayRFI'))

    return returnlist

def _set_bool(settings, key, param, default):
    """
    help function to check the dictionary settings if key exsists,
    and return a on (param, param) or off (empty) tuple depending in
    settings[key] value or to the default (0:off, or 1:on ) value
    """
    if settings.has_key( key ):
        if not settings[key]:   # if its False
            return () # return empty list
        elif settings[key]: # if its True
            return (param, param)
    else: # Else we set it to the default value
        if default: # if default is On (i.e. 1)
            return (param, param)
        elif not default: # if default is Off (i.e. 0)
            return ()

def _parse_settings( settings ):
    """
    set the data release version of the splatalogue compilation

    # data versions (choose)
    # def v2.0
    data_version=v2.0
    or
    data_version=v1.0
    or
    data_version=vall
    """
    returnlist = []
    
    # Data release version
    # first parese the input
    # def v2.0
    if settings.has_key( 'version' ):
        version = settings['version']
    else:
        version = '2.0'
    # now set the parameter
    # def v2.0
    if str(version) in ['2.0', '2', '2.']:
        returnlist.append(('data_version', 'v2.0'))
    elif str(version) in ['1.0', '1', '1.']:
        returnlist.append(('data_version', 'v1.0'))
    elif str(version).lower() in ['all', 'al', 'a']:
        returnlist.append(('data_version', 'vall'))
    else:
        returnlist.append(('data_version', 'vall'))
    # Frequency error limit
    # def off
    # fel=fel
    key = 'felim'
    param = 'fel'
    default = 0
    returnlist.append( _set_bool(settings, key, param, default) )
    # Exclude atmospheric species (triggered)
    # def on
    # no_atmospheric=no_atmospheric
    key = 'no_atm'
    param = 'no_atmospheric'
    default = 1
    returnlist.append( _set_bool(settings, key, param, default) )
    # Show ONLY NRAO Recommended Freq (triggered)
    # def off
    #~ include_only_nrao=include_only_nrao
    key = 'nrao'
    param = 'include_only_nrao'
    default = 0
    returnlist.append( _set_bool(settings, key, param, default) )
    # Exclude potential interstellar species (triggered)
    # def on
    #~ no_potential=no_potential
    key = 'potential'
    param = 'no_potential'
    default = 1
    returnlist.append( _set_bool(settings, key, param, default) )
    # Exclude probable interstellar species (triggered)
    # def on
    #~ no_probable=no_probable
    key = 'probable'
    param = 'no_probable'
    default = 1
    returnlist.append( _set_bool(settings, key, param, default) )
    # Exclude known AST species (triggered)
    # def off
    #~ known=known
    key = 'known'
    param = 'known'
    default = 0
    returnlist.append( _set_bool(settings, key, param, default) )

    # stupid(!) hack to remove empty entries, need to just not add them...
    while 1:
        try:
            returnlist.remove(())
        except (ValueError):
            break
    
    return returnlist

def _parse_frequency(freq, fwidth, funit):
    """
        # frequency
        from=31
        to=31
        frequency_units=GHz
        or 
        frequency_units=MHz
    """
    returnlist = []
    #### FREQUENCY
    # Two casees:
    #   1. A list with length two
    #   2. A integer/float
    if type(freq) == str:
        raise(Exception, 'Wrong format for frequency. Need list or float')
    # First guess : a list of floats with length two
    try:
        returnlist.append( ('from', str(freq[0])) )
        returnlist.append( ('to', str(freq[1])) ) 
    except (IndexError, TypeError):
        # If not a list, should be a float, and fwidth given
        try:
            freq = float(freq)
        except (ValueError):
            raise (Exception, 'Wrong format for frequency. Need list or float')
        if fwidth not in [0, 0.0, None]:
            # with freq and fwidth given, we can calculate start and end
            f1, f2 = freq + _np.array([-1,1]) * fwidth / 2.
            returnlist.append( ('from',  str(f1)) ) 
            returnlist.append( ('to',  str(f2)) )
        else:
            # the fwidth parameter is missing
            raise (Exception, 'The fwidth parameter is missing. '
            'Frequency parameter(s) malformed')
    #### FREQUENCY UNIT
    #
    if funit not in [0, None]:
        if funit.lower() in ['ghz', 'mhz']:
            returnlist.append( ('frequency_units', funit) )
        else:
            print 'Allowed frequency units : \'GHz\' or \'MHz\''
    elif not funit in [0, None]:
        funit = 'GHz'
        returnlist.append( ('frequency_units', 'GHz') )
    return returnlist

def _parse_erange( efrom, eto, eunit ):
    """
        # Energy range (triggered)
        # but if one exists, the energy_range_type must exist
        energy_range_from=10
        energy_range_to=500
        energy_range_type=eu_k
         or
        energy_range_type=el_k
        or
        energy_range_type=el_cm1
        or
        energy_range_type=eu_cm1

    """
    
    returnlist = []
    ### Energy Range
    # form['energy_range_from/to'] is a text field in the form
    # while it is called e_from/to in the function
    
    if efrom == None and eto == None and eunit != None:
        print 'You gave the Enery range type keyword, but no energy range...'
        raise Exception('energy range keywords malformed')
    #~ if (efrom not None) or (eto not None):
    eunit_ref = ['el_cm1', 'eu_cm1', 'el_k', 'eu_k']
        # check that unit is given, and correct
        # or set default (eu_k)
    # set efrom if supplied
    if efrom != None:
        returnlist.append( ('energy_range_from', str(efrom)) )
    # set eto if supplied
    if eto != None:
        returnlist.append( ('energy_range_to',  str(eto)) )
    # check if eunit is given, and tick the corresponding radio
    # button, if none then assume Kelvin
    if eunit != None: #arg.has_key('efrom') or arg.has_key('eto'):
        if eunit.lower() in eunit_ref:
            pass
        else:
            print 'Energy range unit keyword \'eunit\' malformed.'
            raise Exception('eunit keyword malformed')
    else:
        # no value, assuming its in Kelvin (i.e. Eu/kb)
        eunit = 'eu_k'
    # now set the eunit radio button
    returnlist.append( ('energy_range_type', eunit.lower() ) )   
    return returnlist

def _parse_transition( trans ):
    """
        # transition (triggered)
        tran=1-0
    """
    return ('tran', str(trans))
    
def _parse_lill( lill ):
    """
        # line intensity lower limit (triggered)
        #~ lill_cdms_jpl=-5
        #~ or
        #~ lill_sijmu2
        #~ or
        #~ lill_aij
    """
    
    ### Line Intensity Lower Limits
    if lill != None:
        if lill[1].lower() == 'cdms_jpl':
            return ( 'lill_cdms_jpl', str(lill[0]) )
        elif lill[1].lower() == 'sijmu2':
            return ( 'lill_sijmu2', str(lill[0]) )
        elif lill[1].lower() == 'aij':
            return ( 'lill_aij', str(lill[0]) )

def _parameters_ending():

    returnlist = [
    ('export_type','current'),
    ('export_delimiter','colon'),
    ('offset','0'),
    ('limit', str(HIT_LIMIT)),
    ('range','on'),
    ('submit','Export')
    ]
    return returnlist

def _get_results(parameters):
    
    parameters = urlencode(parameters)
    path = SPLAT_FORM_URL  
    req = Request(path, parameters)
    req.add_header("Content-type", "application/x-www-form-urlencoded")
    results = urlopen(req, timeout=SPLATALOGUE_TIMEOUT).read()
    return results

def _parse_results(data, output='astropy.table'):
    """
    Only one output type at the moment, the astropy.table table
    
    """
    #TODO : what if results are empty
    if output == 'astropy.table':
        if not use_astropy:
            #~ print('Astropy not installed, try other output format')
            raise(ImportError('Astropy not installed, try other output format'))
        # get each line (i.e. each molecule)
        rows = data.split('\n')
        # get the names of the columns
        column_names = rows[0]
        column_names = column_names.split(':')
        # clean them up a bit
        for i in _np.arange(len(column_names)):
            column_names[i] = column_names[i].replace('<br>', ' ')
            column_names[i] = column_names[i].replace('<sub>', '_')
            column_names[i] = column_names[i].replace('<sup>', '^')
            column_names[i] = column_names[i].replace('</sup>', '')
            column_names[i] = column_names[i].replace('</sub>', '')
            column_names[i] = column_names[i].replace('&#956;', 'mu')
            column_names[i] = column_names[i].replace('sid[0] is null', '')
            column_names[i] = column_names[i].replace('sid[0] is null', '')
        """
        Column Names should now be:
            ['Species',
             'NRAO Recommended',
             'Chemical Name',
             'Freq-GHz',
             'Freq Err',
             'Meas Freq-GHz',
             'Meas Freq Err',
             'Resolved QNs',
             'Unresolved Quantum Numbers',
             'CDMS/JPL Intensity',
             'S_ijmu^2 (D^2)',
             'S_ij',
             'Log_10 (A_ij)',
             'Lovas/AST Intensity',
             'E_L (cm^-1)',
             'E_L (K)',
             'E_U (cm^-1)',
             'E_U (K)',
             'HFS int',
             'Upper State Degeneracy',
             'Molecule Tag',
             'Quantum Number Code',
             'Linelist']
        """
        rows = rows[1:-1]
        rows = [i.split(':') for i in rows]
        rows = _np.array(rows)
        rows[rows == ''] = -999999
        #~ print column_names
        #~ return rows
        column_dtypes = ['str',        # 'Species',
                        'str',        # 'NRAO Recommended',
                        'str',        # 'Chemical Name',
                        'float',      # 'Freq-GHz',
                        'float',      # 'Freq Err',
                        'float',      # 'Meas Freq-GHz',
                        'float',      # 'Meas Freq Err',
                        'str',        # 'Resolved QNs',
                        'str',        # 'Unresolved Quantum Numbers',
                        'float',      # 'CDMS/JPL Intensity',
                        'float',      # 'S_ijmu^2 (D^2)',
                        'float',      # 'S_ij',
                        'float',      # 'Log_10 (A_ij)',
                        'str',        # 'Lovas/AST Intensity',
                        'float',      # 'E_L (cm^-1)',
                        'float',      # 'E_L (K)',
                        'float',      # 'E_U (cm^-1)',
                        'float',      # 'E_U (K)',
                        'float',      # 'HFS int',
                        'float',      # 'Upper State Degeneracy',
                        'int',        # 'Molecule Tag',
                        'int',        # 'Quantum Number Code',
                        'str']        # 'Linelist']

        funit = str(column_names[3][-3:])
        
        column_units = [None,           # 'Species',
                        None,           # 'NRAO Recommended',
                        None,           # 'Chemical Name',
                        funit,          # 'Freq-GHz',
                        funit,          # 'Freq Err',
                        funit,          # 'Meas Freq-GHz',
                        funit,          # 'Meas Freq Err',
                        None,           # 'Resolved QNs',
                        None,           # 'Unresolved Quantum Numbers',
                        '?',            # 'CDMS/JPL Intensity',
                        'Debye^2',      # 'S_ijmu^2 (D^2)',
                        '?',            # 'S_ij',
                        'log10(s^-1)',  # 'Log_10 (A_ij)',
                        '?',            # 'Lovas/AST Intensity',
                        'cm^-1',        # 'E_L (cm^-1)',
                        'K',            # 'E_L (K)',
                        'cm^-1',        # 'E_U (cm^-1)',
                        'K',            # 'E_U (K)',
                        '?',            # 'HFS int',
                        None,           # 'Upper State Degeneracy',
                        None,           # 'Molecule Tag',
                        None,           # 'Quantum Number Code',
                        None]           # 'Linelist']

        column_names_original = column_names[:]

        #~ column_names = [i.lower() for i in column_names]

        #~ for i in _np.arange(len(column_names)):
            #~ column_names[i] = column_names[i].replace('nrao recommended', 'nrao_rec')
            #~ column_names[i] = column_names[i].replace('chemical name', 'name')
            #~ if 'meas freq err' in column_names[i]:
                #~ column_names[i] = 'mferr'
            #~ elif 'meas freq' in column_names[i]:
                #~ column_names[i] = 'mfreq'
            #~ elif 'freq err' in column_names[i]:
                #~ column_names[i] = 'ferr'
            #~ elif 'freq' in column_names[i]:
                #~ column_names[i] = 'freq'
            #~ column_names[i] = column_names[i].replace('resolved qns', 'resqn')
            #~ column_names[i] = column_names[i].replace('unresolved quantum numbers', 'resqn')

        column_names = ['species',
                 'nrao_rec',
                 'name',
                 'ofreq',
                 'oferr',
                 'mfreq',
                 'mferr',
                 'res_qn',
                 'uresqn',
                 'cdmsjplint',
                 'sijmu2',
                 'Sij',
                 'logaij',
                 'lovasastint',
                 'el_cm',
                 'el_k',
                 'eu_cm',
                 'eu_k',
                 'hfsint',
                 'gu',
                 'tag',
                 'qncode',
                 'list']
        results = Table(data = rows , 
                        names = column_names, 
                        dtypes = column_dtypes)

        
        for i in _np.arange(len(column_units)):
            results.field(i).units = column_units[i]
        return results
    else:
        print('Nothing else than astropy.table output is implemented atm')
        return results


"""

sid[]=

#  Energy level display (triggered)
# 1 : Elower (cm-1)
# 2 : Elower (K)
# 3 : Eupper (cm-1)
# 4 : Eupper (K)
el1=el1
el2=el2
el3=el3
el4=el4

# Line strength display (triggered)
# 1 : CDMS/JPL Intensity
# 2 : Sij mu2
# 3 : Sij
# 4 : Aij
# 5 : Lovas/AST
ls1=ls1
ls2=ls2
ls3=ls3
ls4=ls4
ls5=ls5

# line list (triggered)
# def all on
displayRecomb=displayRecomb 
displayLovas=displayLovas
displaySLAIM=displaySLAIM
displayJPL=displayJPL
displayCDMS=displayCDMS
displayToyaMA=displayToyaMA
displayOSU=displayOSU
displayLisa=displayLisa
displayRFI=displayRFI


# data versions (choose)
# def v2.0
data_version=v2.0
or
data_version=v1.0
or
data_version=vall

# Exclude atmospheric species (triggered)
# def on
no_atmospheric=no_atmospheric

# Exclude potential interstellar species (triggered)
# def on
no_potential=no_potential

# Exclude probable interstellar species (triggered)
# def on
no_probable=no_probable

# Exclude known AST species (triggered)
# def off
known=known

# Show ONLY NRAO Recommended Freq (triggered)
# def off
include_only_nrao=include_only_nrao

# Display Unresolved quantum numbers (triggered)
# def on
show_unres_qn=show_unres_qn

# Show upper degeneracy (triggered)
# def on
show_upper_degeneracy=show_upper_degeneracy

# Display Molecule Tag (triggered)
# def on
show_molecule_tag=show_molecule_tag

# No HFS Display (triggered)
noHFS=noHFS

# Display HFS Intensity (triggered)
displayHFS=displayHFS

# Display Quantum Number Code (triggered)
show_qn_code=show_qn_code

# Display Lab Ref (triggered)
show_lovas_labref=show_lovas_labref

# Display Obs Ref (triggered)
show_lovas_obsref=show_lovas_obsref

# Display Ordered Frequency ONLY (triggered)
show_orderedfreq_only=show_orderedfreq_only

# Display NRAO Recommended Frequencies (triggered)
show_nrao_recommended=show_nrao_recommended


# transition (triggered)
tran=1-0

# frequency
from=31
to=31
frequency_units=GHz
or 
frequency_units=MHz

# line intensity lower limit (triggered)
lill_cdms_jpl=-5
or
lill_sijmu2
or
lill_aij


# Energy range (triggered)
# but if one exists, the energy_range_type must exist
energy_range_from=10
energy_range_to=500
energy_range_type=eu_k
 or
energy_range_type=el_k
or
energy_range_type=el_cm1
or
energy_range_type=eu_cm1




submit=1
"""
