


########################################################################
# UNIT object
# and object that has a unit as well

class Unit(float):
    """
    Constant float type that can have a unit (unit)
    and a description (desc)

    Call : W = Constant(value, 'unit', desc = 'description')
    where 'unit' and desc, are optional.
    Then 'W' just returns the value assigned to it, and can
    be used in calculations as a normal float.
    e.g. W * 2e3 / 5.
    """

    def __new__(self, value, *args, **kwargs):
        # return without the *args and **kwargs
        # to make sure no error is raised in __new__
        return super(Unit, self).__new__(self, value)

    def __init__(self, value, *args, **kwargs):
        # store the different arguments into the class
        self.unit = args[0] if args else None
        self.desc = kwargs.pop('desc', None)




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
# GENERAL FUNCTIONS

def calc_gain(bmin, bmaj, frequency):
    return 8.168E-25 * bmin * bmaj * frequency**2
    
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
def calc_offset(ra, dec, display = True, **kwargs):
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
        if display:
            print 'Input coordinates:'
            print 'RA:\t{0}\nDEC:\t{1}'.format(parse_ra(ra_decimal,string=1),
                                        parse_dec(dec_decimal,string=1))
            print 'Offset: {}'.format(offset_inp)
            print 'New coordinates:'
            print 'RA:\t{}\nDEC:\t{}'.format(parse_ra(new_ra,string=1),
                                        parse_dec(new_dec,string=1))
    elif kwargs.get('offset') == None:
        if kwargs.get('data') != None:
            #~ # calculate the offset from the phase center
            ralist = []
            ralist.append(kwargs['data'].ra_crval)
            #~ ralist.append(parse_ra(kwargs['data'].ra_crval, string=1))
            ralist.append(ra)
            declist = []
            declist.append(kwargs['data'].dec_crval)
            #~ declist.append(parse_dec(kwargs['data'].dec_crval, string=1))
            declist.append(dec)
            #~ # calculate the offset from the two lists
        else:
            ralist = ra
            declist = dec
            #
        #~ ra0_inp = ralist[0]
        #~ ra0 = ra0_inp.split(':')
        #~ ra0 = [float(i) for i in ra0]
        # convert to decimal number in degrees
        #~ ra0_decimal = (ra0[0] + ra0[1]/60.0 + ra0[2]/3600.0)*15.0
        ra0_decimal = ralist[0]
        #~ ra1_inp = ralist[1]
        ra1_decimal = ralist[1]
        #~ ra1 = ra1_inp.split(':')
        #~ ra1 = [float(i) for i in ra1]
        # convert to decimal number in degrees
        #~ ra1_decimal = (ra1[0] + ra1[1]/60.0 + ra1[2]/3600.0)*15.0
        #
        
        #~ dec0_inp = declist[0]
        #~ dec0 = dec0_inp.split(':')
        #~ dec0 = [float(i) for i in dec0]
        # convert to decimal number
        #~ dec0_decimal = sign(dec0[0])*(abs(dec0[0]) + dec0[1]/60.0 +
        dec0_decimal = declist[0]
        
        #~ dec1_inp = declist[1]
        #~ dec1 = dec1_inp.split(':')
        #~ dec1 = [float(i) for i in dec1]
        # convert to decimal number
        #~ dec1_decimal = sign(dec1[0])*(abs(dec1[0]) + dec1[1]/60.0 + dec1[2]/3600.0)
        dec1_decimal = declist[1]
        
        ##### calculate the offset
        # correction factor
        cosdec = cos(dec0_decimal*pi/180.) # convert to radians
        # calculate offsets
        ra_offset = (ra1_decimal - ra0_decimal) * cosdec
        dec_offset = dec1_decimal - dec0_decimal
        
        if display:
            print 'Reference\nRA:\t{0}\nDEC:\t{1}'.format(
                                                parse_ra(ra0_decimal,string=1),
                                                parse_dec(dec0_decimal,string=1))
            print 'Final\nRA:\t{0}\nDEC \t{1}'.format(
                                                parse_ra(ra1_decimal,string=1),
                                                parse_dec(dec1_decimal,string=1))
            print '\nOffset: {0:.4f}, {1:.4f}'.format(ra_offset*3600,
                                                    dec_offset*3600)
        elif not display:
            return ra_offset*3600, dec_offset*3600
    else:
        raise(ParError(ra,dec,kwargs))

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
        return hours, minutes, seconds
    elif type(ra) == type(''):
        h, m, s = array(ra.split(':')).astype('float')
        h += m/60. + s/3600.
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
        print('Did not parse!')
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
    #~ from mpfit import mpfit
    from adapy.fitting.mpfit import mpfit
    #from sys import exit as sysexit
    if verbose:
        print 'Fitting Gaussians, checking input parameters'
    # flatten the parameters, so it is readable for the gaussian fitting
    params = array(params).flatten()

    pos = params[arange(1,len(params),3)]
    if (pos<X.min()).any() or (pos>X.max()).any():
        print('You are trying to fit a Gaussian outside of the \
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
    if verbose:
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
    """ Gaussian 2D 
    a[0] - Amplitude
    a[1] - x grid
    a[2] - y grid
    a[3] - x sigma 
    a[4] - y sigma
    a[5] - PA
    
    
    """
    #
    from scipy import cos, sin, exp
    xp = X*cos(a[5]) - Y*sin(a[5])
    xp0 = a[1]*cos(a[5]) - a[2]*sin(a[5])
    yp = X*sin(a[5]) + Y*cos(a[5])
    yp0 = a[1]*sin(a[5]) + a[2]*cos(a[5])
    #
    f = a[0] * exp(-.5*(((xp - xp0)/a[3])**2 + ((yp - yp0)/a[4])**2))
    return f
def gaussfit2d((X,Y,Z), params=None, err=None, fitheight=0,verbose=0):
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
        if verbose:
            print('--> No parameters given, guestimates\n')
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
        if fitheight:
            params = [height, amplitude, X[0,x0], Y[y0,0], xsig, ysig, 0]
            no_fits = len(params[1:])/6
        elif not fitheight:
            params = [amplitude, X[0,x0], Y[y0,0], xsig, ysig, 0]
            no_fits = len(params)/6
        if verbose:
            print('done.')
    elif params!=None:
        if verbose:
            print '--> Parameters entered, using them\n'
        if fitheight:
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
            no_fits = len(params[1:])/6
        elif not fitheight:
            b = array([list(i) for i in params[1:]]).flatten()
            # switch to radians!
            index = arange(5,len(b),5)
            for i in index:
                b[i] = b[i]*pi/180
            params = b
            if len(params)%6 != 0:
                print ' '
                print 'wrong number of input parameters'
                print '(N*6)+1 parameters where N are the number of gaussians fited'
                print '  height,  amplitude, x0, y0, xisg, ysig, angle(deg)'
                print ' '
                return ; sysexit()
            no_fits = len(params)/6
    #
    #
    

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
        if fitheight:
            S = p[0] # first the baseline
            # then the gaussians
            for i in arange(1,len(p[1:]),6):
                S += gauss2d(p[i:i+6], X, Y)
        elif not fitheight:
            S = 0.0
            # then the gaussians
            for i in arange(0,len(p),6):
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
    if fitheight:
        index = arange(4,len(p1),5) # start at xsig
    else:
        index = arange(3,len(p1),5) # start at xsig
    
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
