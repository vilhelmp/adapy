#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
#  adacore.py
#
#  Core module/library for Astronomical Data
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
#  version 


from cgsconst import *

#----[ DESCRIPTION ]----
"""
Module with functions and objects to analyze Astronomical Data 
(mainly FITS files)


Need : o scipy (and numpy)
       o mpfit.py (optional, for Gaussian fitting)
       o congridding.py (optional, for regridding)
       (o coords (not yet, but later with ability to plot larger regions))
"""


#----[ CHANGE LOG ]----
"""

* 2012 Oct 16
    Moved data handling and manipulation to adacore.py
    move UI to separate module adavis.py
    and data I/O and handling to adacore.py DONE!


"""



#----[ BLUEPRINTS ]----
"""
Top of the list TODO:

TODO : Clean up functions

TODO : LineID plot

TODO : Able to change restfrequency, or at least the V=0 point

TODO : get rid of all text output in this module, move it to adavis, or 
        just the individual __str__ methods

TODO : method calc_offset to Fits object (calc offset from phase center given
        sexadecimal or sim coord input)

TODO : extend the "box" plotting keyword so that it accepts arbitrary
regions (just like the region keyword) <- perhaps sharing code is a good
idea here(?)

TODO : Cannot plot or deduce levels when some pixels have value 'nan'

TODO : load frequency information from fits file if it is there, and create
velocity array as well
        -> more dynamic fits-file loading

TODO : constants move to a module, update the code to use it.

TODO : implement different lineID plot, make it more object oriented

TODO : update the __str__ method of all the classes to print as much
        info as possible

TODO : Moment 2 maps
        -Check others mom2 function (how to get line center?)
        -Need MPFIT Gaussian fitting?

TODO : Clean up code again, remove font handler, or inactivate it
      All fonts should be Times New Roman, to big of a hassle to change
       to Sans-serif font, for the tick labels mostly.
       (works better in later matplotlib?)

TODO : Check that the set_rc function is used in all functions

TODO : Clean up moment 0/1(/2) function

TODO : Change to handle DataObject in all functions.
       (Use v_sys, dist, name from the Fits data object)
        X Spectra
        X Fit
        X moment0
        X moment1
        O moment2
        O The rest

TODO : How to divide it into a runnable program
        Perhaps use **kwargs i.e. dictionaries more which can
        easily be passed along from e.g. OptParser

 TODO : The Error classes should take two input, which parameter
        is raising the error and the reason

------------------------------------------------------------------------
 Less pressing matters:

TODO : P-V diagram - rotatable
       what happens with the coordinates?
       need to calc an offset, but call them x and y or something

 TODO : RMS and other units, e.g. when using SD data, Kelvin instead of Jy.

 TODO : axes_grid1 backwards compability? (for use in CASA)

 TODO : Tick locators are good for small regions/narrow spectra,
        but not for big/wide
       Perhaps change it with the 'box' keyword
            -> e.g a 10th of the box keyword for major locators
       Implemented a parameter locator=[major,minor] as a first step

 TODO : Function to bin in spatial (2D) regime (congrid map),
         change ra/dec_delt etc

 TODO : Spectra class, use for fitting, lineid, binning

 TODO : Interactive multiple Gaussian fitting, i.e. click at start
        end of each line, calculate moments within it

 TODO : Interactive splatsearch

 TODO : What does the VELREF keyword in the GILDAS fits header mean?
       Check GILDAS manual

 TODO : Check Jes functions how they work and try to merge/replace them.

 TODO : Implement it as command line program
        -> partially implemented, should work with spectra
        Make sure all parameters are reachable through **kwargs
        use Dictionary.get('attribute', ValueIfNoKey)

 TODO : Create a partition function class, to access partition functions
       for different molecules.
         o Example adavis.qrot.sulfurdioxide('18(1,2)-17(1,8)AE', Trot=170)
         o Using the transition to identify it, and print shows the
           available ones
         o Store the Q values for different Tex
         o Able to update from Q and Tex values

"""
########################################################################
# USEFUL STRINGS
KMS = u"km\u00b7s\u207b\u00b9"
########################################################################
# ERROR HANDLING (EXCEPTIONS)
#
# input parameter error
class ParError(Exception):
    def __init__(self, value, reason=None):
        self.value = value
    def __str__(self):
        s1 = 'Parameter(s) \"{0}\" is(are) malformed or missing.'.format(self.value)
        return stylify(s1, fg='r')

class FitsError(Exception):
    def __init__(self, value, reason=None):
        self.value = value
    def __str__(self):
        s1 = 'Error reading fits file: {0}'.format(self.value)
        return stylify(s1, fg='r')

#~ class DataError(Exception):


########################################################################
# GENERAL FUNCTIONS

def calc_frequency(vlsr, freq0):
    """
    vlsr in kms
    freq0 in whatever you want out
    """
    from scipy import constants
    return (1-vlsr/(constants.c*1e-3))*freq0
def calc_vlsr (f,f0):
    """ calc vlsr in km/s from two freq """
    from scipy import constants
    return (1-f/f0)*constants.c*1e-3
def calc_sigma(N,rms,v_cdelt):
    from scipy import sqrt
    return sqrt(N)*rms*abs(v_cdelt)



def get_telescope_diameter(telescope):
    from string import upper
    from scipy import where, array
    name = array(['SMA', 'PDBI', 'JCMT', 'AP-H201-F102', 'APEX', 'IRAM30M', 'ALMA'])
    dia = array([6,     # SMA
                 15,    # PdBI
                 15,    # JCMT
                 12,    # APEX (Johans data, odd code for telescope?)
                 12,    # APEX
                 30,    # IRAM30M
                 12])   # ALMA
    try:
        diameter = dia[where(upper(telescope)==name)][0]
    except IndexError:
        diameter = 1
    return diameter
def get_vals(chvals=None, nvals=None):
    from scipy import array
    from matplotlib.pyplot import close
    #
    # if there is no channels supplied
    if chvals==None:
        chvals = array(
        raw_input('input the limits, comma separated: ').split(','),
                    dtype ='float')
    if nvals==None:
        # ask for noise calculation velocity limits
        try:
            nvals = array(
            raw_input('input the noise limits (velocity). \
            comma separated: ').split(','),
                dtype='float')
        except (ValueError):
            print("Since you did not input any or \
            input was wrong we will guess some values...")
            nvals = None
    close(1)
    return chvals, nvals
def calc_offset(ra,dec,**kwargs):
    """

    ra,dec - string with coordinate
    There are three ways of calling the function:

    1.
    calc_offset(ra,dec,offset=(dX,dY))
    if  "offset" is given, it will
    calculate new coordinates from old ones
    offset - offsets in a tuple/list i.e. (dX,dY), [dX,dY]

    2.
    if no "offset" is given, it will
    assume you have input two values in each ra, and dec as list
    and calculate the offset between those two. Where the first one
    is the reference position

    3.
    if "data" is given as a ADAVIS Fits data object,
    it will calculate the offset from the phase center of that data
    i.e. offset from Fits.ra_crval and Fits.dec_crval
    """
    from numpy import cos, pi, array, sign
    #
    if kwargs.get('offset') != None:
        #ra_inp = ra
        ra = ra.split(':')
        ra = [float(i) for i in ra]
        # convert to decimal number in degrees
        ra_decimal = (ra[0] + ra[1]/60.0 + ra[2]/3600.0)*15.0
        #dec_inp = dec
        dec = dec.split(':')
        dec = [float(i) for i in dec]
        # convert to decimal number
        dec_decimal = sign(dec[0])*(abs(dec[0]) + dec[1]/60.0 + dec[2]/3600.0)
        offset_inp = kwargs['offset']
        offset = array(offset_inp)/3600.0
        cosdec = cos(dec_decimal*pi/180.) # convert to radians
        new_ra = ra_decimal + offset[0]/cosdec
        new_dec = dec_decimal + offset[1]
        print 'Input coordinates:'
        print 'RA:\t{0}\nDEC:\t{1}'.format(parse_ra(ra_decimal,string=1),
                                    parse_dec(dec_decimal,string=1))
        print 'Offset: {}'.format(offset_inp)
        print 'New coordinates:'
        print 'RA:\t{}\nDEC:\t{}'.format(parse_ra(new_ra,string=1),
                                    parse_dec(new_dec,string=1))
    elif kwargs.get('offset') == None:
        if kwargs.get('data') != None:
            # calculate the offset from the phase center
            ralist=[]
            ralist.append(parse_ra(kwargs['data'].ra_crval,string=1))
            ralist.append(ra)
            declist=[]
            declist.append(parse_dec(kwargs['data'].dec_crval,string=1))
            declist.append(dec)
            # calculate the offset from the two lists
        else:
            ralist = ra
            declist = dec
        ra0_inp = ralist[0]
        ra0 = ra0_inp.split(':')
        ra0 = [float(i) for i in ra0]
        # convert to decimal number in degrees
        ra0_decimal = (ra0[0] + ra0[1]/60.0 + ra0[2]/3600.0)*15.0
        ra1_inp = ralist[1]
        ra1 = ra1_inp.split(':')
        ra1 = [float(i) for i in ra1]
        # convert to decimal number in degrees
        ra1_decimal = (ra1[0] + ra1[1]/60.0 + ra1[2]/3600.0)*15.0
        #
        dec0_inp = declist[0]
        dec0 = dec0_inp.split(':')
        dec0 = [float(i) for i in dec0]
        # convert to decimal number
        dec0_decimal = sign(dec0[0])*(abs(dec0[0]) + dec0[1]/60.0 +
                    dec0[2]/3600.0)
        dec1_inp = declist[1]
        dec1 = dec1_inp.split(':')
        dec1 = [float(i) for i in dec1]
        # convert to decimal number
        dec1_decimal = sign(dec1[0])*(abs(dec1[0]) + dec1[1]/60.0 +
                    dec1[2]/3600.0)
        # correction factor
        cosdec = cos(dec0_decimal*pi/180.) # convert to radians
        # calculate offsets
        ra_offset = (ra1_decimal-ra0_decimal)*cosdec
        dec_offset = dec1_decimal-dec0_decimal
        print 'Reference\nRA:\t{0}\nDEC:\t{1}'.format(
                                            parse_ra(ra0_decimal,string=1),
                                            parse_dec(dec0_decimal,string=1))
        print 'Final\nRA:\t{0}\nDEC \t{1}'.format(
                                            parse_ra(ra1_decimal,string=1),
                                            parse_dec(dec1_decimal,string=1))
        print '\nOffset: {0:.4f}, {1:.4f}'.format(ra_offset*3600,
                                                dec_offset*3600)
    else:
        raise(ParError(ra,dec,kwargs))



# TODO : get rid of all text output in this module
def stylify (s='Test text', f='n', fg='r', bg='d'):
    """
    Stylify text - change format, forground color and background color.

    Sends back the string 'txt' with the correct foreground unicode
    color start and finish (reset color).

    Arguments
    -----------------
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
        print stylify('\n Warning : '
        'This combination of colors/styles does not work\n','b','r','d')
        raise ParError((f,fg,bg))
    bg += '_bg' # append to the list, the "_bg" ending
    f += "_f" # append "_f" to the formatting list
    if fg == "rand":
        from random import randint
        c_tmp = ["k", "r", "g", "y", "b", "m", "c", "a", "d"]
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




########################################################################
# DATA/COORDINATE PARSING
# to adacore.py
def parse_ra (ra,string=False):
    """

    Parses a simple float coordinate and returns hours, minutes and seconds
    and vice verse.
    Input
        ra : the right ascention coordinate to be parsed

    ---------------------------------------------------------------------------

                            oOO Changelog OOo

    *2010/09 Funciton created

    *2010/12/18 Added documentation

    *2012/07/03 Added calculation in other direction, if input is a string

    """
    from scipy import array
    if type(ra) != type(''):
        a = ra/15
        hours = int(a)
        b = (a-hours)*60
        minutes = int(b)
        seconds = (b-minutes)*60
        if string:
            return '{0:0>2}:{1:2}:{2:0>4.2f}'.format(hours,minutes,seconds)
        return hours,minutes,seconds
    elif type(ra) == type(''):
        h, m, s = array(ra.split(':')).astype('float')
        h += m/60.+s/3600.
        h *= 15
        if string:
            return '{0:5.2f}'.format(h)
        return h
def parse_dec (dec, string=False):
    """

    Parses a simple float coordinate and returns hours, minutes and seconds.
    Input
        dec: the declination coordinate to be parsed

    ---------------------------------------------------------------------------

                            oOO Changelog OOo

    *2010/09 Funciton created

    *2010/12/18 Added documentation

    *2012/07/03 Added calculation in other direction, if input is a string

    """
    from scipy import array, sign
    if type(dec) != type(''):
        degrees = int(dec)
        b = abs((dec-degrees)*60)
        minutes = int(b)
        seconds = (b-minutes)*60
        if string:
            return '{0:0>2}:{1:2}:{2:0>4.2f}'.format(degrees,minutes,seconds)
        return degrees,minutes,seconds
    elif type(dec) == type(''):
        d, m, s = array(dec.split(':')).astype('float')
        # trying to catch all variants of writing a negative sexadecimal
        if -1 in [sign(i) for i in (d, m, s)] or d == -0:
            d = -1*abs(d) - abs(m/60.) - abs(s/3600.)
        else:
            d += m/60.+s/3600.
        if string:
            return '{0:5.2f}'.format(h)
        return d
def parse_linelist(linelist, out='brief'):
    """

    Parse a linelist/CSV file of the lines, so that it is easier to just
    use a lot of lines from Splatalogue or create a list with the lines.

    Input
        linelist : Can be formatted in three different ways

                1.  The path to a COLON separated table file from
                    Splatalogue,where the first two columns are as

                   ----------------------------
                   |   name   |   frequency   |
                   ----------------------------

                    The lines (rows/entries etc) can be commented out with a
                    hash '#'. Easy to get exported from Splatalogue.net.

                2.  A list in the format

                    linelist = ['name1', frequency1, 'name2', frequency2]

                3. A list of lists in the format

                    linelist = [[list_of_names], [list_of_frequencies]]


    ---------------------------------------------------------------------------

                            oOO Changelog OOo

    *2010/06
        Funciton created

    *2010/12
        Doc written

    TODO : when list is  [['Name'],[203.45654]] it cannot output
           the right format.



    """
    from scipy import size, array # arange, zeros,
    def get_lines(linelist):
        # pars a list of lists scenario
        names = [i for i in linelist[::2]]
        if len(names) == 1:
            return linelist
        freqs = [float(i) for i in linelist[1::2]]
        #arr_len = len(names)+len(freqs)
        #list = zeros((arr_len))
        return [names, array(freqs)]
        #
    if size(linelist[1])==1 and type(linelist)==type([]):
        # if it is just a list of lines   e.g. ['Name',203.45654]
        return get_lines(linelist)
    elif size(linelist[1])>1:
        # if it is a list with two lists in it
        # that is, the second element is (strictly) longer than 1
        # e.g. [['Name1','Name2'],[203.45654,203.45654]]
        return linelist
    elif size(linelist[1])==1 and type(linelist) == type(''):
        # in case it is the path to a file
        # load the table with the load ascii table function loadatbl
        try:
            names, freqs = loadatbl(linelist,
                                    dtype='string', rtype='array')[0:2]
            f = []
            n = []
            for i,j in zip(freqs,names):
                if i != '':
                    f.append(i)
                    n.append(j)
            freqs = array(f)
            names = n
        except ValueError:
            names, freqs = loadatbl(linelist, dtype='string', rtype='array')[0:2]
        #freqs = array([float(x) for x in freqs]).astype('float')
        freqs = array(freqs).astype('float')
        # return list, if array freq will be a string
        # cannot mix types in an numpy.ndarray, works in list
        return [names, freqs]
    else:
        print_warning('Did not parse!')
def get_indices (arr,vals,disp=False):
    """

    Get the indices of all the elements between vals[0] and vals[1].
    Alternatively also between vals[2] and vals[3] if they are given.

    Input:
        arr  : the array in which to look for the elements
        vals : a list with either 2 or 4 values that corresponds
               limits inbetween which the indices of the values

    Optional argument(s):
        disp : Bolean parameter, if True it displays start and end
               index and the number of channels inbetween. Only works
               for value lists of length 2.

    Assumes the values in 'arr' is the mid values and that it is evenly
    spaced for all values.

    ********************** Important! **********************************
    The output indices are Python friendly, i.e. they are 0-based. Take
    when using the indices in other software e.g. GILDAS, MIRIAD, which
    are 1-based.

    --------------------------------------------------------------------

                            oOO Changelog OOo

    *2012/02
        Added more documentation, "important" notice about indexing
    *2011/07
        Removed +1 in the output indices to be compatible with rest of
        module, where Pythons 0-based indexing is used.
    *2010/12
        Doc written
    *2010/06
        Funciton created
    """

    from scipy import concatenate, where, array, diff
    dx = abs(.5*diff(arr)[0])
    if len(vals)==4:
        v1,v2,v3,v4 = vals + array([-1,1,-1,1])*dx
        # if the user wants two velocity areas to calculate noise
        low = where((arr>=v1)*(arr<=v2))[0]
        high = where((arr>=v3)*(arr<=v4))[0]
        channels = concatenate((low,high))
    elif len(vals)==2:
        v1,v2 = vals + array([-1,1])*dx
        #channels = where((arr>=v1)*(arr<v2))[0]+1
        # this is because if +1 it is FITS/Fortran safe
        # changed: removed +1 for consistency in program
        channels = where((arr>=v1)*(arr<=v2))[0]
    #
    if disp and len(vals)==2:
        first, last = channels.min(), channels.max()
        n = last-first+1
        print '\nFirst: %d,\n Last: %d\n Nchan: %d\n' % (first, last, n)
    return channels

#
# Help functions for fitting
# 1D
# to adacore.py
def gauss1d(x, params=None, height=None):
    """

    NB: Use this for plotting, the gauss1dfit gaussian function should be
        faster for fitting.

    returns the function vector of input x with parameters p
    just a gaussian function, perhaps to plot a gaussian after the fitting

    fix is the fixed parameters array

    should have subtracted baseline of data before fitting
    i.e. not height is fitted


    gaussian of the form a[0]*exp(-(X-a[1])**2/(2*a[2]**2)*2*sqrt(2*log(2)))
    is fitted. Here a[2] is the FWHM, we multiply with  2*sqrt(2*log(2)) to
    transform it to sigma.
    For an astronomer the FWHM is much more interesting, right?

    TODO : be able to pass height as optional variable

    Changes:
    *using numpy instead of scipy

    """
    from scipy import exp, array, alen, array, log, sqrt, arange
    #
    # Check the input parameters, since this is for plotting, we check them
    #
    # params should always be a list of parameters
    # if it is a series of tuple(s)
    # ie params=((A1,x1,sig1),(A2,x2,sig2))
    if alen(params[0])==3:
        # just put them in a long array
        params = array([list(i) for i in params]).flatten()
        # or in their respective arrays directly?
        #a = params[arange(0,len(params),3)]
        #x = params[arange(1,len(params),3)]
        #dx = params[arange(2,len(params),3)]
    elif params==None:
        # add code to
        # calculate moments!
        raise ParError(params)
    #
    no_fits = len(params)/3
    #
    # check that all parameters are input
    if len(params)%3 != 0:
        raise ParError(params)
    #
    # 2*sqrt(2*log(2)) to convert to sigma (two is cancled out)
    if height==None:
        gaussian = lambda X, a : (
                    a[0]*exp(-(X-a[1])**2/(a[2]**2)*4*log(2)))
    elif height!=None:
        gaussian = lambda X, a : height + (
                            a[0]*exp(-(X-a[1])**2/(a[2]**2)*4*log(2)))
    #
    def func (X, p):
        S=0
        try:
            for i in arange(0,len(p),3):
                S += gaussian(X, p[i:i+3])
        except IndexError:
            raise ParError(p)
        return S
    #
    #
    return func(x, params)
def fit_gauss1d((X,Y),
                params,
                err = None,
                fixlist = None,
                minbool = None,
                minpar = None,
                maxbool = None,
                maxpar = None,
                tie = None,
                verbose = 1,
                full_output=0):
    """
    X - the coordinates of the x-axis (array)
    Y - the data to be fitted (array)
    params - are the initial guesses of the parameters,
            the parameters are grouped in three, with order
            amplitude, peak position and distribution FWHM
            [[AMPL, POS, FWHM]] or ((AMPL, POS, FWHM))
    err - error of the Y data
    fixlist - a list of which parameters to be fixed during fit
            e.g. fixlist = [[0,0,1],[0,0,1]] to fix the width of two
            gaussians that are being fit
    minbool - if != None, will assume that you want to limit something with a
            minimum e.g. [[False,False,True]] to fix the FWHM
    minpar - so then you must specify what that minimum limit is
              e.g. [[0,0,0.5]] for a minimum FWHM of 0.5
    maxbool - here alseo, if != None you will have to specify a
    maxpar - which shows what that value is, you have to suply a value
            even if you are not limiting it so that the specific limiting
            value can be found
    verbose - 0 to be completely quiet
            1 to output final fit parameters (default)
            2 for full verbosity, output each iteration of the fit process
            AND the final fit parameters
    full_output - output the whole shebang?

    MPFIT status messages interpretation:
        0  Improper input parameters.
        1  Both actual and predicted relative reductions in the sum of squares
           are at most ftol.
        2  Relative error between two consecutive iterates is at most xtol
        3  Conditions for status = 1 and status = 2 both hold.
        4  The cosine of the angle between fvec and any column of the jacobian
           is at most gtol in absolute value.
        5  The maximum number of iterations has been reached.
        6  ftol is too small. No further reduction in the sum of squares is
           possible.
        7  xtol is too small. No further improvement in the approximate
           solution x is possible.
        8  gtol is too small. fvec is orthogonal to the columns of the jacobian
           to machine precision.


    TODO : initial guesses for params
                    o intereactive - point and clic
                    o 1D : add a 1 gaussian guessing alorithm
    TODO : minimum width >= X[1]-X[0]
    """
    #from scipy import optimize
    from numpy import exp, array, log, alen, where, arange
    # hstack, zeros, sqrt, diag
    from mpfit import mpfit
    #from sys import exit as sysexit
    print 'Fitting Gaussians, checking input parameters'
    # flatten the parameters, so it is readable for the gaussian fitting
    params = array(params).flatten()

    pos = params[arange(1,len(params),3)]
    if (pos<X.min()).any() or (pos>X.max()).any():
        print_warning('You are trying to fit a Gaussian outside of the \
    data range\n Not allowed. Exciting without fitting')
        return None
    widths = params[arange(2,len(params),3)]
    if (widths<abs(X[1]-X[0])).any():
        raise ParError("Width of a line can not be less than channel width of data (<{0:3.3} kms).".format(abs(X[1]-X[0])))
    # number of gaussians to fit
    no_fits = (len(params))/3
    #
    # check the total number of params
    if len(params)%3 != 0:
        raise ParError(params.reshape((alen(params)/3,3)))

    # that was the most important parameters,
    # now on to the other
    # FIXLIST
    if fixlist != None:
        fixlist = array(fixlist).flatten()
        # length of fixlist must equal the number of parameters
        if alen(fixlist) != alen(params):
            raise ParError(fixlist)
    elif fixlist == None:
        fixlist = [False, False, False]*no_fits
    #
    # MINIMUM LIMITS
    #
    # first we define a minimum limit for the width
    xwidth = abs(X[1]-X[0]) # it is evenly space
    #
    # first case, it is not None, we are giving limits
    if minbool != None:
        # first the boolean array
        minbool = array(minbool).flatten()
        # then the minumum parameter array
        minpar = array(minpar).flatten()
        if alen(minbool) != alen(params) or alen(minbool) != alen(minpar):
            raise ParError((minbool,minpar))
        #
        # any width that need to be set?
        # we know that we are setting the limits, but which FWHM limits are
        # not set to True?
        index = where(minbool[arange(2,alen(minbool),3)]==False)[0]
        # go back to the full array indices

        index = (index+2)+2*index
        # set the bool array to True at the FWHM
        minbool[index] = True
        # set the FWHM minimum parameter to xwidth
        minpar[index] = xwidth
    #
    elif minbool == None:
        # it is not limited
        # here we still would want to limit the gaussian fwhm to the
        # width of one x axis data unit, can't be thinner than that, can it
        minbool = [False,False,True]*no_fits
        minpar = [0, 0, xwidth]*no_fits
    #
    # MAXIMUM LIMITS
    if maxbool != None:
        maxbool = array(maxbool).flatten()
        # then the maximum parameter array
        maxpar = array(maxpar).flatten()
        if alen(maxbool) != alen(params) or alen(maxbool) != alen(maxpar):
            raise ParError((maxbool,maxpar))
    elif maxbool == None:
        # just setting them to False and zero
        maxbool = [False, False, False]*no_fits
        maxpar = [0,0,0]*no_fits
    #
    if tie != None:
        # parse the tied parameters
        tie = array(tie).flatten()
    elif tie == None:
        tie = ['','','']*no_fits
    #
    print '\nDefining fitting function and fitting'
    #
    ## Defining the fitting function and the error function
    # NB: we are fitting FWHM directly, and its squared, so sqrt disapears
    gaussian = lambda X, a : a[0]*exp(-(X-a[1])**2/(a[2]**2)*4*log(2))
    #
    def fitfunc (X, p):
        S=0
        for i in arange(0,len(p),3):
            S += gaussian(X, p[i:i+3])
        return S
    #
    def errfunc(x,y,err=None):
        if err == None:
            def f(p,fjac=None): return [0,(y-fitfunc(x,p))]
        else:
            def f(p,fjac=None): return [0,(y-fitfunc(x,p))/err]
        return f
    # return([status, (y-model)/err]
    # define the parameter dictionary
    PAR = ['Amplitude','Position','Fwhm']
    fitlist = []
    for i in xrange(alen(params)):
        # create the dictionary to be appended to the list=
        dictline = dict(n  = i,
                        value   = float(params[i]),
                        limits  = [float(minpar[i]),float(maxpar[i])],
                        limited = [minbool[i],maxbool[i]],
                        fixed   = fixlist[i],
                        parname = PAR[i%3],
                        tied = tie[i])
        fitlist.append(dictline)

    #print fitlist

    if verbose == 2:
        # verbose level 2, output all fit iterations and final fit params
        quiet = False
    else:
        quiet = True
    #
    mpfit_out = mpfit(errfunc(X,Y,err=err), parinfo=fitlist, quiet=quiet)
    pfit = mpfit_out.params
    pfit_err = mpfit_out.perror
    chi2 = mpfit_out.fnorm
    #
    # any errors
    if mpfit_out.status == 0:
        raise Exception(mpfit_out.errmsg)

    if verbose in [1,2]:
        print('\n')
        print '*'*40
        print 'Results of fitting:\n'
        j = 1
        for i,val in enumerate(pfit):
            fitlist[i]['value'] = val
            if i in arange(0,alen(pfit),3):
                print 'Fit number: %d' % j
                j+=1
            print u"%s :  %2.3f \u00b1 %2.3f" % (fitlist[i]['parname'],val, pfit_err[i])
        print '*'*40
        print u"\u03C7\u00B2 : %2.3f Red. \u03C7\u00B2 : %2.3f DOF : %d\n" % (chi2,(chi2/len(Y)),(len(Y)-len(pfit)))
    #
    #
    if full_output:
        return pfit, pfit_err, chi2, mpfit_out
    else:
        return pfit, pfit_err
# 2D
# to adacore.py
def gauss2d (a, X, Y):
    """ Gaussian 2D """
    #
    from scipy import cos, sin, exp
    xp = X*cos(a[5]) - Y*sin(a[5])
    xp0 = a[1]*cos(a[5]) - a[2]*sin(a[5])
    yp = X*sin(a[5]) + Y*cos(a[5])
    yp0 = a[1]*sin(a[5]) + a[2]*cos(a[5])
    #
    f = a[0]*exp(-.5*(((xp - xp0)/a[3])**2 + ((yp - yp0)/a[4])**2))
    return f
def gaussfit2d((X,Y,Z), params=None, err=None):
    """

    params= (height, (amplitude, x, y, sig_x, sig_y, rota), (amplitude, x, y, sig_x, sig_y, rota)) etc


        xp = cos(rota) * x - sin(rota) * y
        yp = sin(rota) * x + cos(rota) * y
        (rota should be in degrees)
        g = b + a exp ( - ( ((xp - center_xp)/width_x)**2 +
        ((yp-center_yp)/width_y)**2 ) / 2 )

    TODO : move code over to MPFIT
    TODO : initial guesses for params
            o intereactive - point and clic
            o 1D : add a 1 gaussian guessing alorithm
        if(f1.lt.f2)then
          t = f1
          f1 = f2
          f2 = t
          t = sf1
          sf1 = sf2
          sf2 = t
          p = p + 90
        endif
        p = mod(p,180.)
        if(p.lt.-90) p = p + 180
        if(p.gt. 90) p = p - 180
        end

    TODO : Be able to lock a gaussfit parameter
    """

    from scipy import optimize, exp, hstack, array, cos, sin, indices, ravel
    from scipy import sqrt, arange, pi, ceil, diff
    from sys import exit as sysexit
    #
    #
    if params==None:
        print '--> No parameters given, guestimates\n'
        from scipy.ndimage import center_of_mass as cmass
        from scipy.stats import mode
        y0,x0 = cmass(Z)
        row = Z[int(y0),:]
        xsig = sqrt(abs(((arange(row.size)-x0)**2*row).sum()/row.sum()))*diff(X[0,1:3])
        col = Z[:,int(x0)]
        ysig = sqrt(abs(((arange(col.size)-y0)**2*col).sum()/col.sum()))*diff(Y[1:3,0])
        height = mode(Z.ravel())[0][0]
        amplitude = Z.max() - height
        # guess the PA to be 0
        params = [height, amplitude, X[0,x0], Y[y0,0], xsig, ysig, 0]
        print('done.')
    elif params!=None:
        print '--> Parameters entered, using them\n'
        a = array(params[0])
        b = array([list(i) for i in params[1:]]).flatten()
        # switch to radians!
        index = arange(5,len(b),5)
        for i in index:
            b[i] = b[i]*pi/180
        params = hstack((a,b))
        if len(params[1:])%6 != 0:
            print ' '
            print 'wrong number of input parameters'
            print '(N*6)+1 parameters where N are the number of gaussians fited'
            print '  height,  amplitude, x0, y0, xisg, ysig, angle(deg)'
            print ' '
            return ; sysexit()
    #
    #
    no_fits = len(params[1:])/6

    # List of the indices (a):
    #   0   1   2    3    4      5
    #   B, x0, y0, xsig, ysig, theta
    #gau2d = lambda a, X, Y  : a[0]*exp(\
                    #-((X*cos(a[5])-Y*sin(a[5])) - \
                    #(a[1]*cos(a[5])-a[2]*sin(a[5])) )**2 / (2*a[3]**2)\
                    #-((X*sin(a[5])+Y*cos(a[5])) - \
                    #(a[1]*sin(a[5])+a[2]*cos(a[5])) )**2 / (2*a[4]**2))
    #
    #if no_fits == 1:
        #fitfunc = lambda p, X, Y: p[0] + gau2d(p[1:7],X,Y)
    #if no_fits == 2:
        #fitfunc = lambda p, X, Y: p[0] + gau2d(p[1:7],X,Y)\
                                        #+gau2d(p[7:13],X,Y)
    #if no_fits == 3:
        #fitfunc = lambda p, X, Y: p[0] + gau2d(p[1:7],X,Y)\
                                        #+gau2d(p[7:13],X,Y)\
                                        #+gau2d(p[13:19],X,Y)
    def fitfunc (p, X, Y):
        S = p[0] # first the baseline
        # then the gaussians
        for i in arange(1,len(p[1:]),6):
            S += gauss2d(p[i:i+6], X, Y)
        return S
    # error function
    # List of the indices (p):
    #   0   1   2   3   4      5     6
    #   A,  B, x0, y0, xsig, ysig, theta
    if err == None:
        errfunc = lambda p, X, Y, Z: (fitfunc(p, X, Y) - Z).ravel() # Distance to the target function
    else:
        errfunc = lambda p, X, Y, Z: ((fitfunc(p, X, Y) - Z)/err).ravel() # Distance to the target function
    # the fitting
    p1, success = optimize.leastsq(errfunc, params, args=(X,Y,Z))
    #
    #switch back to degrees and positive width
    # check so degrees betweem 0 and 180
    index = arange(4,len(p1),5) # start at xsig
    for i in index:
        p1[i] = abs(p1[i])
        p1[i+1] = abs(p1[i+1])
        # return degrees
        p1[i+2] = p1[i+2]*180/pi-90
        # removed this to not do it twice (is done in parse_gau2fit)
        #if p1[i+2]<0: p1[i+2] += ceil(abs(p1[i+2]/180))*180
        #if p1[i+2]>180: p1[i+2] -= ceil(abs(p1[i+2]/180))*180
    #
    return p1,success,no_fits
def parse_gau2dfit (fwhm1, fwhm2, sfwhm1, sfwhm2, pa, spa):
    """ Function doc

    fwhm in pixels, pa in degrees
    """
    from scipy import sqrt, log, append
    #from sys import exit as sysexit
    #if abs(DATA.ra_cdelt) != abs(DATA.dec_cdelt):
    #    print 'Not equal scale of x and y pixels'
    #    return ; sysexit()
    ##
    ## convert to asecs
    #delta = abs(DATA.ra_cdelt)
    #fwhm1 = fwhm1*delta
    #fwhm2 = fwhm2*delta
    #sfwhm1 = sfwhm1*delta
    #sfwhm2 = sfwhm2*delta
    #
    # check if major is in fact minor axis
    if fwhm1<fwhm2:
        # switch!
        fwhm1, fwhm2 = (fwhm2, fwhm1)
        sfwhm1, sfwhm2 = (sfwhm2, sfwhm1)
        pa = pa + 90
    #
    # now check the PA
    pa = pa%180
    if pa<-90: pa += 180
    if pa>90: pa -= 180

    return (fwhm1, fwhm2, sfwhm1, sfwhm2, pa, spa)
def gauss2d_decon ((bmaj1, bmin1, theta1, bmaj2, bmin2, theta2), ang='rad'):
    """
    Deconvolves one gaussian  with parameters bmaj1, bmin1, theta1 (major,
    minor, PA) with another (bmaj2,bmin2, theta2)
    all in FWHM and radians (if deg is wanted, set ang='deg')


    uses:
    pi, cos, sin, arctan2, sqrt, min,

    """
    from scipy import pi, cos, sin, arctan2, sqrt
    #
    # check the ang keyword, if deg, go over to radians from deg
    if ang=='deg':
        theta1 *= pi/180
        theta2 *= pi/180
    #
    # define some calculations
    alpha  = (bmaj1*cos(theta1))**2 + (bmin1*sin(theta1))**2 - \
             (bmaj2*cos(theta2))**2 - (bmin2*sin(theta2))**2
    beta   = (bmaj1*sin(theta1))**2 + (bmin1*cos(theta1))**2 - \
             (bmaj2*sin(theta2))**2 - (bmin2*cos(theta2))**2
    gamma  = 2 * ( (bmin1**2-bmaj1**2)*sin(theta1)*cos(theta1) -\
                   (bmin2**2-bmaj2**2)*sin(theta2)*cos(theta2) )
    #
    # calculate the intermediate results
    s = alpha + beta
    t = sqrt((alpha-beta)**2 + gamma**2)
    limit = 0.1*min(bmaj1,bmin1, bmaj2, bmin2)**2
    #
    # now check if result is illigal/close to a point source
    if alpha < 0 or beta < 0 or s < t:
        bmaj, bmin, bpa = [0, 0, 0]
        #
        # now check if result is close to a point source
        tmp_par =.5*(s-t)
        if tmp_par < limit and alpha > -limit and beta > -limit:
            success = 1
        #
        # it was not close to point source, but results are thus illigal
        else:
            success = 2
    #
    # since (if) everything is ok, go ahead and calculate the bmaj, bmin & bpa
    else:
        bmaj = sqrt(.5*(s+t))
        bmin = sqrt(.5*(s-t))
        #
        # bpa
        if (abs(gamma)+abs(alpha-beta)) == 0:
            bpa = 0
        else:
            bpa = 0.5 * arctan2(-gamma,(alpha-beta))
        success = 0
        #
        # go back to degrees if asked for
        if ang=='deg':
            bpa *= 180/pi
    #
    # send back the results
    return (bmaj, bmin, bpa, success)
########################################################################
# DATA HANDLING
# to adacore.py
# classes etc
#
# FITS DATA CLASS
# MAIN class
# Main data object, needs fits file as input
# keywords can be appended to object later
# to complement fits header information
class Fits:
    """
    ------------------------------------------
    Adavis Fits Object (Data object)

    Usage :
    ObjectName = FITS(PathToFitsFile)

    ------------------------------------------

    Should be able to read:
    Map, Cube and SD

    TODO : for DataObject loading, learn it to parse the FITS type INDEX
       that Class can output for -32 Bits, (single dish data)

    TODO : maybe change the velocity attr, so that it follows naming
    similar to ra and dec, e.g. self.v_delt, self.v_arr

    TODO : Create a frequency array as well, much simpler later on then

    """
    def __init__(self, fitsfile, telescope=None, vsys=0, distance=0):
        """

        attributes
        ------------------------
        datatype
            possible values : 'SDSPECT', 'CUBE', 'IMAGE'
        telescope
            supported values : 'SMA', 'PDBI', 'IRAM30M', 'APEX, 'ALMA'
            the name of the telescope
        diameter
            diameter of the telescope
        v_arr
            array with every channels velocity, not corrected for the systemic velocity
        dist
            distance to the source

        TODO : If the rotational matrice is non empty, load data and rotate it (?)
        TODO : What if add a function to grid it to a certain size, say 512x512?
                o So it is possible to combine data with different res.
                    - If extent keyword is given for both datasets, what happens?
                o Better to do this when plotting (se above comment)?
        TODO : More robust loading (detecting the axis etc)?
                o Minimum of different types of data
                    - NAXIS==2 (Freq/Vel & Flux/Intensity) - 1D Spectra
                    - NAXIS==2 (RA & DEC) - 2D map
                    - NAXIS==3 (RA & DEC & Flux/Intensity) 3D Spectral map
                    - Polarization data?
                    (i.e, in SD spectra need to get rid of 3 axes
                        self.d = self.d[0][0][0])
        TODO : make loadcube delete an axis, along with the hdr keywords if all
               the axis keywords/values are empty/null (delet hdr keyword really needed?)
                    o only to loaded data, not save to raw data (fits file)


        OVERALL change:

        make it load the fits info in steps. start with the RA and DEC keywords
        then go over to loading frequency/velocity array, and so on...
        lastly determine what type it is



        """

        #imports
        from pyfits import open as fitsopen
        #from  sys import
        from os.path import getsize
        from scipy import where, array, nan
        from string import upper
        from sys import exit as sysexit
        #

        # create the class, but without any init script.
        # a class (object) the easy way
        print u'Loading fitsfile :  %s ' % stylify(str(fitsfile),fg='g')
        s  = getsize(fitsfile)
        print " Size %0.2f MB" % (s/(1024.*1024.))
        f = fitsopen(fitsfile)
        self.hdr, self.d = f[0].header, f[0].data
        #self.d = self.d[0] # this is if the stokes axis is present,
        # but it should not be there anymore
        f.close()
        # save the fitsfile, perhaps the path too, for updating it
        self.fitsfile = fitsfile
        # the telescope diameter
        # first check if there was keyword sent in
        if telescope!=None:
            self.hdr.update('TELESCOP', telescope)
            self.telescope = str(telescope)
        #
        if self.hdr.has_key('TELESCOP'):
            #~ name = array(['SMA', 'PDBI', 'JCMT', 'AP-H201-F102', 'IRAM30M'])
            #~ dia = array([6, 15, 15, 12, 30])
            #~ try:
                #~ self.diameter = dia[where(upper(self.hdr['TELESCOP'])==name)][0]
            #~ except IndexError, ex:
                #~ self.diameter = 1
            self.diameter = get_telescope_diameter(self.hdr['TELESCOP'])
            self.telescope = self.hdr['TELESCOP']
        else:
            self.diameter= 1
            self.telescope = None
        if self.hdr.has_key('LINE'):
            self.linename = self.hdr['LINE']
        #
        #
        # spectra, 1 : 3 axis and 3rd axis is >1 in size
        # image, 2 : 2 axis (or 3 axis, 3rd is =1 in size)
        # cube,3 : spectral cube
        #
        # simple solution to the common extra empty STOKES axis
        # and sometimes even an extra empty axis for ?
        from numpy import diff, arange
        #~ while self.d.shape[0] == 1:
            #~ self.d = self.d[0]

        try:
            self.restfreq = self.hdr['RESTFREQ'] # in Hertz
        except KeyError:
            try:
                self.restfreq = self.hdr['RESTFRQ'] # in Hertz
            except KeyError:
                print ('No frequency information.')


        if self.hdr['NAXIS']==4 and self.d.shape[0:2] == (1,1):
            self.datatype = ('IMAGE',2)
            self.d = self.d[0][0]
        #naxis = self.hdr['NAXIS']
        #axshape = self.d.shape
        #if axshape[0] array([i>1 for i in a.shape[1:]]).all()
        # an image, only 2 dimensions
        elif self.hdr['NAXIS']==2 and self.hdr['NAXIS1']>1 and self.hdr['NAXIS2']>1:
            #if self.hdr['NAXIS']==3 and self.hdr['NAXIS1']>1 and self.hdr['NAXIS2']>1 and self.hdr['NAXIS3']==1:
            # image, not SD spectra or anything,
            # really 2D and greater extent than 1x1
            self.datatype = ('IMAGE',2)
            pass
        #
        # spectral image cube (extra axis for frequency/velocity)
        elif self.hdr['NAXIS']==3 and self.hdr['NAXIS1']>1 and self.hdr['NAXIS2']>1 and self.hdr['NAXIS3']==1:
            self.datatype = ('IMAGE',2)
            # extra if the continuum image has the freq and width
            self.freq = self.hdr['CRVAL3']
            self.freqwidth = self.hdr['CDELT3']
        # a spectra! the 3rd axis is longer than 1
        elif self.hdr['NAXIS']>=3 and self.hdr['NAXIS3']>1:
            # spectral cube
            # only support for velo-lsr in 3rd axis
            self.datatype = ('CUBE',3)
            # load the third axis
            # need frequency!
            while self.d.shape[0] == 1:
                self.d = self.d[0]
            ##### have to add loading of frequency and calculate velocity
            # UGLY HACK BELOW, BEWARE!
            # need to be changed to a more flexible code...
            hdr_values = [self.hdr[i] for i in self.hdr.keys()]
            if 'VELO' in hdr_values:
                print('VELO')
                velax = str([x for x in self.hdr.keys() if x[:-1]=='CTYPE' and 'VELO' in self.hdr[x]][0][-1:])
                vel_info = True
            elif 'VELO-LSR' in hdr_values:
                print('VELO')
                velax = str([x for x in self.hdr.keys() if x[:-1]=='CTYPE' and 'VELO' in self.hdr[x]][0][-1:])
                vel_info = True
            elif 'VRAD' in hdr_values:
                print('VRAD')
                velax = str([x for x in self.hdr.keys() if x[:-1]=='CTYPE' and 'VRAD' in self.hdr[x]][0][-1:])
                vel_info = True
            else:
                print('No velocity axis defined')
                vel_info = False
            if vel_info:
                self.v_type = velax
                self.v_crpix = self.hdr['CRPIX'+self.v_type]-1
                self.v_crval = self.hdr['CRVAL'+self.v_type]
                self.v_ctype = self.hdr['CTYPE'+self.v_type]
                self.v_cdelt = self.hdr['CDELT'+self.v_type]
                self.v_naxis = self.hdr['NAXIS'+self.v_type]
                self.v_cdeltkms = self.v_cdelt/float(1e3)
            # load frequency and calculate velocity stuff
            if not vel_info:
                from scipy import sign
                if 'FREQ' in hdr_values:
                    freqax = str([x for x in self.hdr.keys() if x[:-1]=='CTYPE' and 'FREQ' in self.hdr[x]][0][-1:])
                    self.f_type = freqax
                    self.f_crpix = self.hdr['CRPIX'+self.f_type]-1
                    self.f_crval = self.hdr['CRVAL'+self.f_type]
                    self.f_ctype = self.hdr['CTYPE'+self.f_type]
                    self.f_cdelt = self.hdr['CDELT'+self.f_type]
                    self.f_naxis = self.hdr['NAXIS'+self.f_type]
                    #self.f_cdeltkms = self.f_cdelt/float(1e3)
                    self.f_arr = ((arange(0,self.f_naxis)-self.f_crpix)*self.f_cdelt+self.f_crval) # in Hz
                else:
                    print('No velocity or frequency axis defined')
                    raise FitsError('Could not load FITS file')
                # velocity
                # UGLY hack warning...
                self.v_crpix = self.f_crpix
                self.v_crval = calc_vlsr (self.f_crval, self.restfreq)*1e3
                self.v_ctype = 'VELO'
                self.v_cdelt = abs(calc_vlsr (self.f_arr[0], self.restfreq)*1e3-calc_vlsr (self.f_arr[1], self.restfreq)*1e3)*sign(self.f_cdelt)*(-1)
                self.v_naxis = self.f_naxis
                self.v_cdeltkms = self.v_cdelt/float(1e3)
                #self.f_arr = ((arange(0,self.f_naxis)-self.f_crpix)*self.f_cdelt+self.f_crval)

            #data = loadvelocity(data, velax)
            # loading velocity information

            # plus one because to start at 0 is wrong, need 1 to v_naxis <- this is WRONG
            # start at 0 because we
            self.v_arr = ((arange(0,self.v_naxis)-self.v_crpix)*self.v_cdelt+self.v_crval)/float(1e3) # so it is i kms
            self.v_rangekms = self.v_arr.max()-self.v_arr.min()
            # not good if vcdletkms has more than 3 significant digits.
            #self.v_arr = self.v_arr.round(3)
            #self.v_cdeltkms = round(self.v_cdeltkms,2)

            #start = self.v_crval-self.v_cdelt*(self.v_crpix-1)
            #stop =  self.v_crval+self.v_cdelt*(self.v_naxis-self.v_crpix)
            #arr = arange(start,stop-1,self.v_cdelt)/float(1e3)
            #print self.v_arr-arr
            # calculate the FOV = 58.4*lambda/D*3600 asec


            self.fov = 58.4*(3.e8/self.restfreq)/float(self.diameter)*3600.
            print 'Field of view: %.2f asecs, for dish size: %.1f m' % (self.fov, self.diameter)
            #print self.veltype, self.v_crpix, self.v_crval, self.v_cdeltkms, self.v_naxis
            print 'Velocity range \t: {0:.2f} km/s'.format(self.v_rangekms)
            print 'Velocity step \t: {0:2.4f} km/s'.format(self.v_cdeltkms)
            #
            # now if we want to have the spectral array as well to use
            if self.hdr.has_key('RESTFREQ'):
                self.restfreq = self.hdr['RESTFREQ']
            #if self.hdr['NAXIS']==4 and  self.hdr['NAXIS4']==1:
            #    self.d = self.d[0]
            #
            # this was just a test
            # SD pointing spectra
        elif self.hdr['NAXIS']>1 and self.hdr['NAXIS2']==1 and self.hdr['NAXIS3']==1:
            self.datatype = ('SDSPECT',1)
            self.v_cdelt = self.hdr['DELTAV']
            self.v_cdeltkms = self.hdr['DELTAV']/float(1e3)
            self.v_crpix = self.hdr['CRPIX1']-1
            self.v_naxis = self.hdr['NAXIS1']
            self.v_crval = self.hdr['VELO-LSR']
            self.v_arr = ((arange(0,self.v_naxis)-self.v_crpix)*self.v_cdelt+self.v_crval)/float(1e3)
            self.restfreq = self.hdr['RESTFREQ'] # in Hertz
            # huh?
            self.fov = 58.4*(3e8/self.restfreq)/(self.diameter)*3600
            if 'BEAMEFF' in self.hdr:
                self.beameff = self.hdr['BEAMEFF']
            if 'FORWEFF' in self.hdr:
                self.forweff = self.hdr['FORWEFF']
            #self.d = self.d[0][0][0] # specific for this data...
        #~ elif 'Miriad fits' in self.hdr['ORIGIN']:
        # below to load CLASS bits -32
        #~ elif 'FITS_rec' in str(type(self.d)) and not self.hdr['NAXIS']:
            #~ self.d = self.d[0][0]
            #~ self.datatype = 'SDSPECT',1
            #~ #
            #~ # self.d.dtype shows SPECTRUM and WAVE for CLASS data
            #~ #
            #~ self.v_cdelt = self.hdr['DELTAV']
            #~ self.v_cdeltkms = self.hdr['DELTAV']/float(1e3)
            #~ self.v_crpix = self.hdr['CRPIX1']-1
            #~ self.v_naxis = self.hdr['NAXIS1']
            #~ self.v_crval = self.hdr['VELO-LSR']
            #~ self.v_arr = ((arange(0,self.v_naxis)-self.v_crpix)*self.v_cdelt+self.v_crval)/float(1e3)
            #~ self.restfreq = self.hdr['RESTFREQ'] # in Hertz
            #~ self.fov = 58.4*(3e8/self.restfreq)/(self.diameter)*3600
        else:
            # if it is not an image or a spectral cube
            print_error('The dimensions of the data is wrong\n at least the header keywords indicate that.\n The data has '+str(self.hdr['NAXIS'])+' axes. \n\n Perhaps use the removeaxis script?\n')
            sysexit()
        print 'Datatype : {0}'.format(self.datatype[0])
        # perhaps check in the header?
        # velref probably at what velocity that middle of spectra is?
        self.v_sys = float(vsys)
        self.dist = float(distance)
        #
        # FREQUENCY ARRAY
        #
        # construct the frequency array!
        # the 3rd axis longer than 1, and 4th axis is the frequency
        # if the data is constructed in gildas
        if self.datatype[0] in ['CUBE', 'SDSPECT']:
            self.v_arr_syscorr = self.v_arr - self.v_sys
        #
        # load the coordinate parameters
        # for the CRPIXNax parameter I take -1 because
        # FITS starts at 1 and Python starts at 0, hence in
        # an array, crpix-1 will show the Python position of the crpix
        # DEC
        decax = str([x for x in self.hdr.keys() if x[:-1]=='CTYPE' and 'DEC' in self.hdr[x]][0][-1:])
        self.dec_cdelt = self.hdr['CDELT'+decax]*3600 # arcs
        self.dec_npix = self.hdr['NAXIS'+decax]
        self.y_npix = self.hdr['NAXIS'+decax]
        self.dec_crpix = self.hdr['CRPIX'+decax]-1
        self.dec_crval = self.hdr['CRVAL'+decax]
        # RA
        raax = str([x for x in self.hdr.keys() if x[:-1]=='CTYPE' and 'RA' in self.hdr[x]][0][-1:])
        self.ra_cdelt = self.hdr['CDELT'+raax]*3600 # arcs
        self.ra_npix = self.hdr['NAXIS'+raax]
        self.x_npix = self.hdr['NAXIS'+raax]
        self.ra_crpix = self.hdr['CRPIX'+raax]-1
        self.ra_crval = self.hdr['CRVAL'+raax]
        if self.datatype[0] in ['CUBE','IMAGE']:
            # create a extent keyword
            #~ ylen, xlen = self.d[0].shape
            #~ ycoords = arange(-ylen/2,ylen/2,1)*self.dec_cdelt
            #~ xcoords = arange(-xlen/2,xlen/2,1)*self.ra_cdelt
            #~ left, right = xcoords[0],xcoords[-1]
            #~ bottom, top = ycoords[0],ycoords[-1]
            #~ extent=(left,right,bottom,top)
            X = array([0,self.ra_npix-1]) # self.*_npix-1 because we're
            Y = array([0,self.dec_npix-1]) # slicing the python-way
            left,right = (X-self.ra_crpix)*self.ra_cdelt
            bottom,top = (Y-self.dec_crpix)*self.dec_cdelt
            self.extent = (left,right,bottom,top)
            #self.extent = (left,right,bottom,top)
            #~ xcoords = arange(-(self.ra_crpix),(self.ra_npix-self.ra_crpix),1)*self.ra_cdelt
            #~ ycoords = arange(-(self.dec_crpix),(self.dec_npix-self.dec_crpix),1)*self.dec_cdelt
            #~ print xcoords[0],xcoords[-1]
            #~ print left,right
            #~ print ycoords[0],ycoords[-1]
            #~ print bottom,top
        try:
            # convert Beam size from degrees to asecs
            self.bmaj = self.hdr['BMAJ']*3600
            self.bmin = self.hdr['BMIN']*3600
            self.bpa = self.hdr['BPA']
        except KeyError, ex:
            msg='Header keywords (bmaj,bmin,bpa) incomplete and not loaded.'
            print_warning(msg)
            #~ self.bmaj = None
            #~ self.bmin = None
            #~ self.bpa = None
            #~ self.gain = None
        if 'BUNIT' in self.hdr:
            self.unit = self.hdr['BUNIT']
            # units
            if 'JY/BEAM' in upper(self.unit):
                self.unitpixel = u"Jy\u00b7beam\u207b\u00b9\u00b7channel\u207b\u00b9"
                self.unitint = u"Jy\u00b7beam\u207b\u00b9\u00b7" + KMS
            elif 'K' in upper(self.unit):
                self.unitpixel = u"K\u00b7channel\u207b\u00b9"
                self.unitint = u"K\u00b7" + KMS
        else:
            print_warning('No beam unit in header.')
            self.unitpixel = "INTENSITY"
            self.unitint = "INTEGRATED-INTENSITY"
        # calculate the GAIN of the observations (interferometric observations?)
        if self.datatype[0] in ['CUBE','SDSPECT'] and hasattr(self,'restfreq') and hasattr(self,'bmin'):
            # gain depends on restfreq being there
            self.gain = 8.168e-25*(self.restfreq)**2*self.bmin*self.bmaj
        #
        # Object name
        self.obj = self.hdr['OBJECT']
    def __str__(self):
        print '\n','='*40
        print ' '*8,'FITS file\n'
        print 'Data type : %s' % str(self.datatype[0])
        if self.datatype[1] in [3]:
            print 'Shape of image cube : {0}'.format(self.d.shape)
        print 'Object : %s' % self.obj
        if hasattr(self,'beameff'):
            print 'Beam Efficiency : {0:3.4f}'.format(self.beameff)
        if hasattr(self,'forweff'):
            print 'Fwd Efficiency : {0:3.4f}'.format(self.forweff)
        #
        print ''
        if self.datatype[0] != 'SDSPECT':
            self.ra_size = abs(self.ra_cdelt)*self.ra_npix
            self.dec_size = abs(self.dec_cdelt)*self.dec_npix
            print 'Spatial size of image\n RA\t: %2.3f asec\n DEC\t: %2.3f asec' % (self.ra_size, self.dec_size)
        print 'Phase center '
        print '  RA : {0}'.format(parse_ra(self.ra_crval,string=1))
        print ' DEC : {0}'.format(parse_dec(self.dec_crval,string=1))
        if hasattr(self,'bmaj') and hasattr(self,'bmin') and hasattr(self,'bpa'):
            print '\nBeam info'
            print ' Beam major axis : {0}'.format(self.bmaj)
            print ' Beam minor axis : {0}'.format(self.bmin)
            print ' Beam position angle : {0}'.format(self.bpa)
        #
        print ''
        if hasattr(self,'restfreq'):
            if (1E-9*self.restfreq)<1:
                freq = 1E-6*self.restfreq
                freq_unit = 'MHz'
            else:
                freq = 1E-9*self.restfreq
                freq_unit = 'GHz'
            print 'Rest frequency of data : {0} {1}'.format(freq,freq_unit)
        return '\n ADAVIS - Fitsfile Object \n'
    def parse_pxlcoord (self, x, y):
        """ Function doc """
        xoffset = (x-self.ra_crpix)*self.ra_cdelt
        yoffset = (y-self.dec_crpix)*self.dec_cdelt
        return xoffset, yoffset
    def parse_region(self, region, f=False):
        """
        Parser for the region parameter, three different possibilities to supply
        the region command:

            o region = [i1, j1, i2, j2]
                The four corners of a square around the object, in offset from
                phase center position.

            o region = [i1, j1, a]
                The center coordinate (i1, j1) and the side (a) of a square around
                the center coordinate (i1, j1).

            o region = [d1, d2]
                Just the square sides length, will be centered on the phase center.

        All the coorinates are given in lenghts and offsets (in asec) from the
        data center as displayed normally in radio data.

        Inspired by the miriad 'region' parameter

        ---------------------------------------------------------------------------

                                oOO Changelog OOo

        *2010/06 Funciton created

        *2010/10(11) Doc written and some errors in the code corrected (+-1 in
        some places)

        *2010/12/09 in len(region)==3, changed the division with an abs()
        array([-region[2],region[2]])/(2*data.ra_cdelt) to abs(2*data.ra_cdelt).
        In len(region)==2 same change, now it is correct, I hope.

        *2010/12/13 the previous "fix" made the len=3 procedure to be erronous.
        corrected it

        *2011/10/03 incorporated into the Fits class

        """
        from scipy import ceil, floor, array
        from sys import exit as sysexit
        if len(region)==4:
            xcheck = region[0]==region[2]
            ycheck = region[1]==region[3]
            #~ if region[0]<region[2]: # if you enter it as in miriad i.e. (-5,-5,5,5)
                #~ reg2 = region[2]
                #~ reg0 = region[0]
                #~ region[0] = reg2
                #~ region[2] = reg0
            #x1, x2 = (data.ra_npix+1)/2 + array([region[0],region[2]])/abs(data.ra_cdelt) + array([0,xcheck])
            #y1, y2 = (data.dec_npix+1)/2+ array([region[1],region[3]])/abs(data.dec_cdelt)+ array([0,ycheck])
            #
            x1, x2 = array([region[0],region[2]])/self.ra_cdelt + self.ra_crpix + array([0,xcheck])
            y1, y2 = array([region[1],region[3]])/self.dec_cdelt + self.dec_crpix + array([0,ycheck])
            #
        elif len(region)==3:
            check = region[2]==0
            #x1, x2 = (data.ra_npix+1)/2 + array([-region[2],region[2]])/(2*abs(data.ra_cdelt)) + region[0]/data.ra_cdelt + array([0,check])
            #y1, y2 = (data.dec_npix+1)/2+ array([-region[2],region[2]])/(2*abs(data.dec_cdelt)) +region[1]/data.dec_cdelt+ array([0,check])
            #
            x1, x2 = self.ra_crpix + region[0]/self.ra_cdelt + array([-region[2],region[2]])/abs(2*self.ra_cdelt) + array([0,check])
            y1, y2 = self.dec_crpix + region[1]/self.dec_cdelt + array([-region[2],region[2]])/abs(2*self.dec_cdelt) + array([0,check])
            #
        elif len(region)==2:
            xcheck = region[0]==0
            ycheck = region[1]==0
            #x1, x2 = (data.ra_npix+1)/2 + array([-1,1])*region[0]/abs(data.ra_cdelt)  + array([0,xcheck])
            #y1, y2 = (data.dec_npix+1)/2+ array([-1,1])*region[1]/abs(data.dec_cdelt) + array([0,ycheck])
            #
            x1, x2 = array([-region[0],region[0]])/(2*abs(self.ra_cdelt)) + self.ra_crpix + array([0,xcheck])
            y1, y2 = array([-region[1],region[1]])/(2*abs(self.dec_cdelt)) + self.dec_crpix + array([0,ycheck])
            #
        elif():
            print ('Error, region keyword malformed')
            sysexit(1)
            #
        # so that we are returning usable pixel coordinates
        if f==False:
            x1,x2,y1,y2 = array([x1,x2,y1,y2]).round().astype('int')
        return x1,x2,y1,y2
    def calc_fov(self):
        # method to calculate FOV after the correct telescope name/diameter
        # has been input and thus correcting the current FOV of
        # the DataObject
        if self.telescope!=None:
            self.diameter = get_telescope_diameter(self.telescope)
        elif self.diameter == 1:
            print 'You have not changed either the diameter of the telescope or the telescope name'
        self.fov = 58.4*(3.e8/self.restfreq)/float(self.diameter)*3600.
    def calc_rms(self, nvals, area):
        from scipy import sqrt,array
        i1,i2,j1,j2 = self.parse_region(area)
        n_channels = get_indices(self.v_arr, nvals)
        # just to find out which channels (start, stop) to print
        if len(nvals)==2:
            n = array([n_channels.min(),n_channels.max()])
            nv = self.v_arr[n]
            print "RMS calculated in intervals {0} ({1}) and region {2}".format(n, nv,nvals,area)
        if len(nvals)==4:
            n_1 = get_indices(self.v_arr,array(nvals)[:2])
            n_1min = min(n_1)
            n_1max = max(n_1)
            n_2 = get_indices(self.v_arr,array(nvals)[2:])
            n_2min = min(n_2)
            n_2max = max(n_2)
            #n = array([n_channels.min(),n_channels.max()])
            #nv = self.v_arr[n]
            print "RMS calculated in intervals {0} and {1} ({2}) and region {3}".format([n_1min,n_1max], [n_2min,n_2max],nvals,area)
        rms_data = self.d[n_channels]
        self.rms = sqrt(((rms_data[:, j1:j2, i1:i2])**2).mean())
        del rms_data
    def add_line(self, name, frequency=None, channels=None, width=None):
        """
        Add identified line(s) to the class

        TODO : update the fits file as well?
        '"""
        try:
            known_lines = self.known_lines
        except AttributeError, ex:
            known_lines = {}
        known_lines[204.38343] = {'name' : 'SO$_2$','frequency' : frequency, 'channels' : channels, 'width' : width}
    #
    # method to change the v_sys
    def change_v_sys (self, v_sys):
        self.v_sys = v_sys
        # now, change the v_arr_syscorr array as well
        if self.datatype[0] in ['CUBE', 'SDSPECT']:
            self.v_arr_syscorr = self.v_arr - self.v_sys
    #
    def change_dist (self, dist):
        self.dist = dist # unit of pc
# UV-FITS DATA CLASS
class Uvfits:
    """
    Reads uv-fits data...



    --------------------------------------------------------------------
    Normal structure of UV-fits data:

    Header : same as for normal fits

    Data : Gropu data

    --------------------------------------------------------------------

    TODO :  Assumes that CRVAL4 is frequency, is that always true?
            Make it more robust, look for the frequency keyword, either as
            "restfreq" or as a "crvalX"

    TODO :  UV fit method
    TODO : __init__ method, print information about:
                - Phase center
                - No correlations
                - No baselines
                - No antennas
                - Telescope
    """
    def __init__(self, uvfitsfile, telescope=None, vsys=0, distance=0):
        """

        Reads the uvfits and calculates useful things, e.g. u,v,w,
        phase and amplitude

        """
        from pyfits import open as pfopen
        from scipy import sqrt, pi, arctan2
        from cgsconst import CC
        f = pfopen(uvfitsfile)

        if f[0].header['NAXIS1'] != 0:
            print "error: this file may not be a UV FITS."
            raise FileError('File format error.')
        f.info()
        try:
            self.hdu = f[0]
        except:
            print "error: cannot open uv data HDU."
        self.hdr = self.hdu.header
        self.data = self.hdu.data
        #f.close() # is this really needed for pyfits file objects?
        """
        The standard unit is to give UU and VV in seconds (??!?)
        So we have to convert to whatever we want.
        """
        # unit nano seconds
        self.u_nsec = self.hdu.data.par(0) * 1.e+9
        self.v_nsec = self.hdu.data.par(1) * 1.e+9
        self.w_nsec = self.hdu.data.par(2) * 1.e+9
        # unit kilo-lambda
        #CC_cm = a.CC*1e2 # light speed in cm/s
        freq = self.hdu.header['CRVAL4'] #TODO
        #lmd = lsp / freq
        # u_klam = uu * CC_cm / (CC_cm/freq)
        self.u_klam = self.hdu.data.par(0) * freq * 1.e-3
        self.v_klam = self.hdu.data.par(1) * freq * 1.e-3
        self.w_klam = self.hdu.data.par(2) * freq * 1.e-3
        # unit meters
        self.u_m = self.hdu.data.par(0) * CC*1e-2
        self.v_m = self.hdu.data.par(1) * CC*1e-2
        self.w_m = self.hdu.data.par(2) * CC*1e-2
        # uv distance
        self.uvdist_nsec= sqrt(self.u_nsec**2 +self.v_nsec**2)
        self.uvdist_klam = sqrt(self.u_klam**2 +self.v_klam**2)
        self.uvdist_m = sqrt(self.u_m**2 +self.v_m**2)
        # visibility data set (COMPLEX)
        visi_index = len(self.hdu.data.parnames)
        if self.hdu.header['NAXIS']  == 7:
            self.visdata = self.hdu.data.par(visi_index)[:,0,0,0,0,0,:]
        #~ self.visdata = self.hdu.data.data[:,0,0,0,0,0,:]
        elif self.hdu.header['NAXIS']  == 6:
            self.visdata = self.hdu.data.par(visi_index)[:,0,0,0,0,:]
        # load the re, im and weight arrays
        self.re = self.visdata[:,0]
        self.im = self.visdata[:,1]
        self.weight = self.visdata[:,2]
        # now calculate the amplitude and phase
        self.amplitude = sqrt(self.re**2 + self.im**2)
        self.phase = arctan2(self.im, self.re) / pi * 180.
        # following 1.0e6 is just for GILDAS, change if needed
        print('NB : Error calculated from weights assuming GILDAS '
        'data(i.e. frequencies in MHz).')
        self.error = 1/sqrt(self.weight*1.0e6)
    def avgamp(self, avg):
        """
        averages amplitude over 'avg' number of uvdist units
        perhaps just 'average', and average everything..?
        """
        return (0,0)
    def __str__():
        return 'Not implemented yet.'

# MOMENTS DATA CLASS
# to adacore.py
# Calculates moment 0 and 1, to use in
# the moment_map function
class Moments:
    """
    Calculate moments
    output : moment object with all info

    calculate the different moments (0/1/2) for a dataset
    append it to the data-class (DataObject) and return it

    idl code
    ; Select which part of the velocity axis to consider
    vx=where(va_l ge 6.0 and va_l le 8.0)
    ; mom0
    mom0 = dmom(i,j,0) = total(line(i,j,vx))
    ; only mom1 if signal high enough (i.e. 3sigma)
    dmom(i,j,1) = total(reform(line(i,j,vx))*va_l(vx))/mom0
    check
    """


    ### not shure if it is allowed to use "FITS" as input here.
    # perhaps need some other input name, for clarity
    def __init__ (self, Fits, chvals, nsig):
        """
        moment class initialiser

        input :


        DONE : Check if I take enough levels, i.e. that the self.maximum
               really is the maximum of the WHOLE 2D array
                -> Yes it is

        """
        from scipy import sqrt, alen, flipud, arange, array, ones, nan
        # -> never use binned array
        # -> never use velocities from/with the v_sys corrected data
        # get the data from the cube
        # copy header for easy acess to stuff
        self.hdr = Fits.hdr
        self.channels = get_indices(Fits.v_arr, chvals)
        imgs = Fits.d[self.channels]
        ## MOMENT 0
        # calculate the moment 0
        self.zero = imgs.sum(axis=0) * abs(Fits.v_cdeltkms)
        #Isum = imgs.sum(axis=0)*abs(Fits.v_cdeltkms) # make a copy for masking <3sigma values in mom1 map̈́
        ## STATISTICS of MOMENT 0 (sigma, min, max, levels)
        # other statistics
        self.sigma = sqrt(alen(imgs)) * Fits.rms * abs(Fits.v_cdeltkms)
        self.minimum = self.zero.min()
        self.maximum = self.zero.max()
        # calculate levels, start at 1 sigma, jump 1 sigma
        # one for positive and one for negative
        # concatenate before displaying if want certain start & jump
        self.levels_neg = -1 * arange(self.sigma, abs(self.minimum) + 2 * self.sigma, self.sigma)
        self.levels_pos = arange(self.sigma, self.maximum + 2 * self.sigma, self.sigma)
        #levels = arange(nsig*moment0_sigma,moment0_max+moment0_sigma,nsjump*moment0_sigma)
        ## MOMENT 1
        # create array that matches imgs array
        # only the velocities that we want, i.e. Fits.v_arr[self.channels]
        velocities_matrix = array([ones(imgs.shape[1:]) * i for i in Fits.v_arr[self.channels]] )
        # removed where(), because a boolean array works fine
        # find out where we calculate the moment 1, i.e. 3 sigma level
        Isum = imgs.sum(axis=0)
        lt_3sigma = (self.zero < (nsig * self.sigma)) * (self.zero > (-1.0 * nsig * self.sigma))
        Isum[ lt_3sigma ] = nan
        Ivsum = (imgs * velocities_matrix).sum(axis=0)
        # calculate the denominator, the sum of all images
        # in our velocity interval

        # MOMENT 1
        #
        # M1 = V(x,y)
        # = sum( v_i * I(x,y,v_i)) / sum(I(x,y,v_i))
        #
        self.one = Ivsum / Isum

        # MOMENT 2
        #
        # M2 = sqrt[ sum( I(x,y,v_i) * (v_i - V(x,y))**2 ) / sum(I(x,y,v_i)) ]
        #
        # M2 = sqrt[ sum( I(x,y,v_i) * (v_i - M1)**2 ) / sum(I(x,y,v_i)) ]
        #
        top = imgs * (velocities_matrix - self.one)**2
        division = abs(top.sum(axis=0) / Isum)
        self.two = sqrt(division)

#
# SPECTRUM DATA CLASS
#
# to adacore.py
class Spectrum:
    """ Class doc """
    ### not shure if it is allowed to use "FITS" as input here.
    # perhaps need some other input name
    def __init__ (self, Fits, **args):
        """
        Class initialiser
        ----
        Arguments (args):
        region = []


        TODO : should the "lines" input accept other than the native type?
        TODO : supply the filename of "writetofile" of identify_lines method

        """
        # copy the header and other useful stuff
        # a bit redundant, don't you think?
        self.hdr = Fits.hdr.copy() # this is a dictionary!
        self.v_arr = Fits.v_arr
        self.v_cdeltkms = Fits.v_cdeltkms
        self.v_cdelt = Fits.v_cdelt
        self.v_crpix = Fits.v_crpix
        self.v_crval = Fits.v_crval
        self.v_sys = Fits.v_sys
        self.restfreq = Fits.restfreq
        self.unitpixel = Fits.unitpixel
        self.unitint = Fits.unitint
        #
        if Fits.datatype[0] == 'SDSPECT':
            print stylify("SD-SPECTRUM - region keyword not doing anything.",fg='y')
            self.d = Fits.d
        elif Fits.datatype[0] in ['CUBE','IMAGE']:
            if 'region' in args:
                pass
            else:
                args['region'] = (0,0,0)
            self.region = args['region']
            x1,x2,y1,y2 = Fits.parse_region(args['region'])
            area_region = ((y2-y1)*(x2-x1))
            self.d = (Fits.d[:,y1:y2,x1:x2].sum(axis=1).sum(axis=1))/float(area_region)
        if hasattr(Fits,'unit'):
            self.unit = Fits.unit
        self.binned = 0 # it is not binned at this point
        if 'lines' in args:
            self.lines = parse_linelist(args['lines'])
            try:
                self.linedata = loadatbl(args['lines'], dtype='string', rtype='array')
            except:
                pass
            #print self.lines[1][0],type(self.lines[1][0])
    def __str__(self):
        if self.binned:
            print 'Spectrum has been binned'
        print 'Extracted from region: {0}'.format(self.region)
        return '<ADAVIS Spectrum Object>'
    def bin_spectrum(self, binning, bintype='mean'):
        from string import lower
        from scipy import alen, arange, array
        from sys import exit as sysexit
        binning = int(binning)
        self.binning = binning
        ##### temporary solution, saving old stuff
        class Original: pass
        Original.d = self.d
        Original.v_arr = self.v_arr
        Original.v_cdelt = self.v_cdelt
        Original.v_cdeltkms = self.v_cdeltkms
        self.Original = Original
        #
        if lower(bintype) == 'resample':
            from congridding import congrid
            # congridding, proper resampling of data
            self.d = congrid(self.d,(alen(self.d)/binning,),centre=True, method='neighbour')
            self.v_arr = congrid(self.v_arr,(alen(self.v_arr)/binning,))
            #
            self.v_cdeltkms = self.v_cdeltkms*binning
            self.v_cdelt = self.v_cdelt*binning
        elif lower(bintype) == 'mean':
            if alen(self.d)%binning!=0:
                print 'Bin has to be evenly devide the number of channels: %d' % alen(self.d)
                sysexit()
            #  Old method - simple binning, just average
            indices = arange(0,alen(self.d),binning)
            self.d = array([self.d[x:x+binning].sum(axis=0)/binning for x in indices])
            self.v_arr = array([self.v_arr[x:x+binning].sum(axis=0)/binning for x in indices])
            #
            self.v_cdeltkms = self.v_cdeltkms*binning
        elif binning == 0 or binning <0:
            print stylify("\nERROR:\n Variable \"bin\" has to be 1 for no binning, or above 1 \n\
            for the number of channels to bin")
        # print out information about the binning
        print '='*40
        print ' '*11,"Binning of data\n"
        print "No channels to bin : %d" % self.binning
        print "Velocity step : %f" % self.v_cdeltkms
        if bintype=='mean':
            print 'Type of binning : Simple mean over selected no. bin channels'
        elif bintype=='resample':
            print 'Type of binning : Resampling - 1D interpolation'
        # set the "binned" flag to True! (i.e. larger than 0)
        # every time we bin, it increases with the number of the binning parameter
        # hence the number of channels that it has been binned is repr
        # by this parameter
        self.binned +=self.binning
    def calc_rms(self,Fits,nvals,region='quarter'):
        from scipy import array, sqrt
        from string import upper
        print '='*40
        print ' '*11,'Noise statistics\n'
        # calculate the rms from the channels in the spectra
        # accounts for it even if it is binned
        # image rms
        # change x1,x2,y1,y2 to quarter region
        # change so that when binning, the rms i calculated
        # x1,x2
        # calc two RMS, one for bin and one without
        # perhaps change this, it just a factor of 1/sqrt(self.binning)?
        # also, change to some other algorithm of deducing the rms
        # i.e. histogram or similar
        #
        """
        TODO : specify arbitrary region (also string "quarter") to
                specify what region to calc rms over

        """
        if Fits.datatype[0] == 'SDSPECT':
            if self.binned>1: # if it is larger than 0
                self.rms = sqrt(((Fits.d[get_indices(Fits.v_arr,nvals)])**2).mean()/self.binned)
                #self.rms = sqrt(((Fits.d[get_indices(Fits.v_arr,nvals),j1:j2,i1:i2])**2).mean()/self.binned)
            else:
                self.rms = sqrt(((Fits.d[get_indices(Fits.v_arr,nvals)])**2).mean())
        else:
            zlen, ylen, xlen = Fits.d.shape
            ydelt = ylen/6
            xdelt = xlen/6
            i1,i2 = xlen/2-xdelt, xlen/2+xdelt
            j1,j2 = ylen/2-ydelt, ylen/2+ydelt
            if self.binned>1: # if it is larger than 0
                self.rms = sqrt(((Fits.d[get_indices(Fits.v_arr,nvals),j1:j2,i1:i2])**2).mean()/self.binned)
            else: # if it is 0
                self.rms = sqrt(((Fits.d[get_indices(Fits.v_arr,nvals),j1:j2,i1:i2])**2).mean())
        self.rms_mjy = self.rms*1e3
        # the sensitivity
        #TODO : if the unit is K km/s how come we divide?
        self.sensitivity = self.rms/sqrt(abs(self.v_cdeltkms))
        # the channels used
        ind = get_indices(self.v_arr,nvals)
        ind1, ind2 = ind.min(), ind.max()
        if self.binned:
            print u' RMS \t\t: {0:4.2f} m{1} (Unbinned: {2:3.1f})'.format(self.rms_mjy, Fits.unitpixel, self.rms_mjy*self.binned**0.5)
        else:
            print u' RMS \t\t: {0:2.3f} m{1}'.format(self.rms_mjy, Fits.unitpixel)
        print u' Sensitivity \t: {0:2.3f} m{1}'.format(self.sensitivity*1e3, Fits.unitint)
        print u' Channels used \t: {0:d}, {1:d} ({2} {3})'.format(ind1, ind2, str(nvals), KMS)
        print u' RMS Region \t: {0} arcsec'.format(region)
        print u' Channel width \t: {0:2.3f} {1}'.format(abs(self.v_cdeltkms), KMS)
        if self.binned:
            print 'Binning parameter : {0}'.format(self.binned)
    def fit_lines(self, **args):
        """
        TODO : calculate errors for the fit etc, look at Kaper et al (1966)\
                and Condon (1996)
        TODO : calculate the sum (integrated) intensity
                and apply Beam_eff? (for SD, in K kms-1)
        TODO : guess parameters if only one gaussian?
        TODO : interactive?
        TODO : streamline what is printed out, and what is saved in
                the Fit-class

        Intensities:
        from
        http://www.iram.es/IRAMES/otherDocuments/manuals/index.html
        and the document
        http://www.iram.es/IRAMES/otherDocuments/manuals/Report/cali_rep_ddo970205.ps


        T*A temperatures are brightness temperatures of an
        equivalent source which fills the entire 2pi sterradians of
        the forward beam pattern of the telescope. To obtain the
        brightness temperature of and equivalent source just filling
        the main beam (def. as main beam brightness temperature
        Tmb), antenna temperatures have to be multiplied bu the
        ratio of the forward and main beam efficiency Beff.

        Tmb = Feff/Beff * T*A


        """
        from scipy import array, sqrt, log, alen, arange, pi, diff
        if 'linefit' not in args:
            # add either interactive fitting or
            # some guessing algorithm
            print 'Need first guesses for the fitting atm.'
            return False
        elif 'linefit' in args: # first guess parameters supplied!
            # first copy the fit-dictionary
            fitting = args['linefit'].copy()
        # start print msgs
        print '=' * 40
        print ' ' * 11, "Line fitting\n"
        ##### check if error is supplied
        if 'error' not in fitting:
            # if no error is supplied, use the calculated
            # main class RMS
            if hasattr(self, 'rms'):
                fitting['error'] = self.rms
            elif  not hasattr(self, 'rms'):
                errmsg = 'No error supplied.\n\
                Either supply with {\'error\': VALUE} or cal\
                calculate with Spect.calc_rms(nvals=[n1,n2,n3,n4])'
                print stylify(errmsg,fg='r')
                raise ParError(fitting)
        elif 'error' in fitting:
            print 'We have error from input, not from Spectrum'
        ##### small useful functions
        fwhmfromsig = 2*sqrt(2*log(2)) # the constant
        fwhm = lambda x: fwhmfromsig*x
        sigma = lambda x: x/fwhmfromsig
        ##### fits 1D gaussian(s) to spectral data
        if 'chvals' in args: # if we supplied limits for the fit
            print args['chvals']
            ch = get_indices(self.v_arr,args['chvals'])
            Fx = self.d[ch]
            X = self.v_arr[ch]
        else: # else use the eveything
            Fx = self.d
            X = self.v_arr
        #
        p = fitting['params'] # parameters for the fit
        #
        if 'limmin' not in fitting:
            fitting['limmin'] = None
            fitting['minpar'] = None
        if 'limmax' not in fitting:
            fitting['limmax'] = None
            fitting['maxpar'] = None
        if 'fixlist' not in fitting:
            fitting['fixlist'] = None
        if 'tie' not in fitting:
            fitting['tie'] = None
        #
        from time import time
        t1 = time()
        #
        fitting_results = fit_gauss1d((X,Fx), params=p, fixlist=fitting['fixlist'],
                minbool=fitting['limmin'], minpar=fitting['minpar'],
                maxbool=fitting['limmax'], maxpar=fitting['maxpar'],
                err=fitting['error'], tie=fitting['tie'], verbose=0, full_output=1)

        if fitting_results==None:
            print(stylify('\n\n No fitting done...',f='b',fg='r'))
        elif fitting_results!=None:
            ###### subclass fit

            class Fit: pass

            Fit.params, Fit.errors, Fit.chi2, Fit.mp = fitting_results
            #~ params, errors, chi2, mp = fitting_results
            print ' Done in %2.3f seconds' % (time()-t1)
            #
            print ' Number of fits : ', alen(Fit.params)/3
            print(' Fit status : ',
                  Fit.mp.status,
                  '(if 0, it should have halted)')
            print(' Chi2 : {0}, reduced : {1}\n'.format(
                                                    Fit.chi2,
                                                    Fit.chi2/float(len(Fx))))
            # now, parse output of fitting and print it out on screen
            Fit.line_widths = []
            Fit.frequencies_fitted = []
            Fit.frequencies_corrected = []
            Fit.vel_shifts = []
            Fit.freq_shifts = []
            Fit.gauss_intensities = []
            Fit.sum_intensities = []
            Fit.sum_intensities2 = []
            Fit.line_names = []
            Fit.peak_intensities = []
            Fit.error_gauss_intensities = []
            #
            Fit.nfits = alen(Fit.params)/3
            j = 0
            from scipy import sqrt
            for i in arange(0,len(Fit.params),3):
                # add 1 because channel 1 is in pos 0
                fwhm = Fit.params[i+2]
                half_fwhm = fwhm/2.
                ampl = Fit.params[i]
                pos = Fit.params[i+1]

                # first figure out the extent of the gaussian (the line)
                # jump half a channel down and up so that it finds the correct channels
                lower_half, upper_half = (pos + array([-1,1])*half_fwhm)
                lower,upper = (pos + array([-1,1])*fwhm)
                lower2,upper2 = (pos + array([-1,1])*fwhm*2)
                #channels = where((velocity>lower)*(velocity<upper))[0]+1
                channels_half = get_indices(self.v_arr,
                                            [lower_half,upper_half])
                channels = get_indices(self.v_arr,[lower,upper])
                channels2 = get_indices(self.v_arr,[lower2,upper2])
                #draw_highlight_box(ax_kms, params[i+1], params[i+2]*3)
                # apply v_sys correction,
                #so that we use v_sys for estimating the correct
                # frequency for the line,
                #especially importat when using line-identification
                frequency_fitted = calc_frequency(pos, self.restfreq/1e9)
                frequency_corrected = calc_frequency(pos -
                                            self.v_sys, self.restfreq/1e9)
                Fit.line_widths.append(fwhm)
                Fit.frequencies_corrected.append(frequency_corrected)
                Fit.frequencies_fitted.append(frequency_fitted)
                Fit.peak_intensities.append(self.d[channels2].max())
                #
                if hasattr(Fits, 'beameff') and hasattr(Fits, 'forweff'):
                    constant = Fits.forweff/Fits.beameff
                else:
                    constant = 1
                #gauss_int = (sqrt(2 * pi) *
                #            sigma(Fit.params[i+2]) *
                #            Fit.params[i])*constant
                gauss_int = ampl*sigma(fwhm)*sqrt(2*pi)*constant
                sum_int = (self.d[channels].sum()*abs(self.v_cdeltkms))*constant
                sum_int2 = (self.d[channels2].sum()*abs(self.v_cdeltkms))*constant
                Fit.gauss_intensities.append(gauss_int)
                Fit.sum_intensities.append(sum_int)
                Fit.sum_intensities2.append(sum_int2)
                Fit.error_gauss_intensities.append(1.064 *
                                sqrt((fwhm*Fit.errors[i])**2 +
                                        (ampl*Fit.errors[i+2])**2))
                if hasattr(self, 'lines'):
                    # frequency shift = rest frequency - measured frequency
                    # velocity shift  =  c (freq_shift/freq._rest)
                    # 'lines' contain the name and frequency of
                    # the identified lines, that is all the
                    # lines we're fitting
                    # TODO : self.lines[1] are strings, change to floats
                    #print self.lines[1][j],type(self.lines[1][j])
                    freq_shift = self.lines[1][j] - frequency_corrected
                    vel_shift = CC*1e-2 * freq_shift/self.lines[1][j] * 1E-3 # in kms
                    Fit.freq_shifts.append(freq_shift)
                    Fit.vel_shifts.append(vel_shift)
                    Fit.line_names.append(self.lines[0][j])
                print  stylify('Fit number : {0}'.format(j+1),fg='g',bg='k')
                intensity_string = u" Intensity: Fit= {0:2.4f}, Data= {1:2.4f} (\u00b1FWHM), {2:2.4f} (\u00b12*FWHM) {3}".format(
                                gauss_int,sum_int,sum_int2,self.unitint)
                print intensity_string
                parameter_string=u" Ampl= {0:2.3f} (\u00b1{1:2.3f}) {2}, Pos= {3:2.3f} (\u00b1{4:2.3f} {5})".format(
                                    Fit.params[i],
                                    Fit.errors[i],
                                    self.unit,
                                    Fit.params[i+1],
                                    Fit.errors[i+1],
                                    KMS)
                print stylify(parameter_string,fg='b')
                print u" Width= {0:2.3f} (\u00b1{1:2.3f}) {2} (FWHM, \u03c3={3:2.3f})".format(Fit.params[i+2],Fit.errors[i+2],KMS,sigma(Fit.params[i+2]))
                # calculate frequency and velocity offset if linelist exists
                if hasattr(self,'lines'):
                    # print offset (freq. & vel.)
                    print "Id molecule : {0}".format(self.lines[0][j])
                    print u"Frequency shift : {0:2.5} GHz Vel shift : {1:5.5} {2}".format(freq_shift,vel_shift,KMS)
                frequency_string =\
                u' Frequency : {0:3.9f} GHz (v_sys corrected)'.format(frequency_corrected)
                print stylify(frequency_string,fg='b')
                print u' FWHM      : {0}, {1} ({2}) (0-based) ([{3:.2f}, {4:.2f}] {5})'.format(channels_half.min(), channels_half.max(),(channels_half.max()-channels_half.min()+1), lower_half, upper_half,KMS)
                print u' \u00b1FWHM   : {0}, {1} ({2}) (0-based) ([{3:.2f}, {4:.2f}] {5})'.format(channels.min(), channels.max(), (channels.max()-channels.min()+1), lower,upper,KMS)
                if self.binned:
                    channels_nobin = get_indices(self.Original.v_arr,[lower,upper])
                    channels_half_nobin = get_indices(self.Original.v_arr,[lower_half,upper_half])
                    print u'Original channels :\n \t FWHM width  : {0}, {1} (\u00b1 1/2FWHM) (0-based)'.format(channels_half_nobin.min(), channels_half_nobin.max())
                    print  u' \t \u00b1FWHM width : {0}, {1} ({2}) (0-based) \n'.format(channels_nobin.min(), channels_nobin.max(),(channels_nobin.max()-channels_nobin.min()+1))
                j+=1
            #
            Fit.line_widths = array(Fit.line_widths)
            print 20*'- '
            print u'Mean FWHM : {0:2.1f} \u00b1{1:2.2f} km\u00b7s\u207b\u00b9\n'.format(Fit.line_widths.mean(),Fit.line_widths.std())
            Fit.xarr = arange(X[0],X[-1],(diff(X)[0]/4))
            # lastly get the Fit class into the Spectrum class (self)
            self.Fit = Fit
    def identify_lines(self, **args):
        # later when the kwargs is implemented...
        # checks that we have done a fit first
        #if not kwargs['linefit']:
        #    print 'If you want lineid you need linefit'
        #    raise ParError(kwargs['lineid'], kwargs['linefit'])
        #
        # only works after a fit has been done, this should be changed
        # so that it works with both a native fit, and with supplied
        # frequencies
        #
        #print args
        if not args.has_key('writetofile'):
            args['writetofile'] = 0
        if not hasattr(self, 'Fit'):
            print_warning('Spectrum not fitted, aborting line id')
            raise ParError('linefit - need to fit lines')
        if hasattr(self, 'lines'):
            print_warning('Line list exists, overwriting...')
        # now import the splatsearch module and other stuff
        import splatsearch as spl
        from scipy import arange, array
        print 'Trying to indentify candidates for the fitted lines.'
        frequency_pairs = []
        if 'nfwhm' in args:
            nfwhm = args['nfwhm']
        else:
            nfwhm = 1.5
        print nfwhm
        for i in arange(0, len(self.Fit.params), 3):
            # correct for v_sys, to get correct frequency for the
            # correct frequency range in the splatalogue search
            # when calculating the "correct" frequency, we subtract v_sys
            vel_lower, vel_upper = (self.Fit.params[i+1] -
                                    self.v_sys +
                                    array([-1,1]) *
                                    self.Fit.params[i+2] *
                                    nfwhm)
            # frequency increases when velocity decreases...
            freq_lower = calc_frequency(vel_upper,self.restfreq/1e9)
            freq_upper = calc_frequency(vel_lower,self.restfreq/1e9)
            frequency_pairs.append([freq_lower,freq_upper])
        list_of_species = []
        list_of_frequencies = []
        number = 1
        lineids = []
        if args['writetofile']:
            with open(args['writetofile'],'w') as f:
                f.write('#Results file for Line ID\n#ADAVIS.py - Magnus Persson\n')
                f.write('#{0:15} {1:10} {2:10} {3:10}   {4:20}\n'.format('Species','Frequency','Smu2','Eu(K)','UResQNr'))
        for i in arange(len(frequency_pairs)):
            tmpspecies, tmpfreq = [],[]
            #~ df=8e-3 # range to find line
            CSI = "\x1b["
            start =CSI+'1m'+CSI+'32m'+CSI+'40m'
            end = CSI+'m'
            print '\n'+start+'Line number : '+str(number)+'\t\t\t\t'+end
            print 'Frequency : {0}  GHz'.format(self.Fit.frequencies_corrected[i])
            # N, species, name, freq, freqerr, freqtype, cfreq, cfreqerr,
            # mfreq, mfreqerr, res_qns (10), ures_qns, cdmsjpl_I, Smu2, Sij,
            # log10Aij, lovasAST_I, ELcm, ELK, EUcm, EUK, u_degen, \
            # mol_tag, QNr, llist
            result = spl.splatsearch(
                        freq=frequency_pairs[i],
                        send=1,
                        display=1,
                        linelist=['jpl','cdms'])
                        #e_to=1000)
            if args['writetofile']:
                with open(args['writetofile'],'a') as f:
                    f.write('\n# Line no. {0}\n'.format(i+1))
            if result!=None:
                species, freq = result[1],result[3]
                smu2, eu = result[13], result[20]
                uresqnr = result[10]
                llist = result[24]
                for j in arange(len(freq)):
                    list_of_species.append(species[j])
                    list_of_frequencies.append(freq[j])
                    tmpspecies.append(species[j])
                    tmpfreq.append(freq[j])
                #~ for i in arange(len(freq)):
                    #~ if i>0 and freq[i]!=freq[i-1]: # remove duplicates
                        #~ list_of_species.append(species[i])
                        #~ list_of_frequencies.append(freq[i])
                    #~ elif i==0:
                        #~ list_of_species.append(species[i])
                        #~ list_of_frequencies.append(freq[i])
                    #~ else:
                        #~ pass
                    if args['writetofile']:
                        with open(args['writetofile'],'a') as f:
                            f.write(
                                '{0:20}  {1: <10}  {2:>10}  {3:>10}'
                                '  {4:<25} {5}\n'.format(
                                species[j],freq[j],smu2[j],eu[j],
                                uresqnr[j],llist[j]))
                lineids.append([tmpspecies,tmpfreq])
            else:
                if args['writetofile']:
                    with open(args['writetofile'],'a') as f:
                        f.write('NO LINES\n')
                print('No lines found...')
            number+=1
            # done now define the linelist
        self.lines = [list_of_species,list_of_frequencies]
        self.lineids = lineids



























