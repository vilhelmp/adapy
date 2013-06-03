#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
#  adavis.py
#
#  UI functions to interact with (Radio-) Astronomical Data.
#  Needs the adacore.py module to function.
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
#  version  2.0b
#
#
#
"""
Python module for handling Astronomical data and analysis

Reads fits cubes and spectra.

Plots spectra, mom0/1/2,

Interfaces with

Transphere (Dust continuum model?)
RADEX

RATRAN

Splatalogue

NASA ADS
"""

"""
Python tips!

"A Foolish Consistency is the Hobgoblin of Little Minds"
Python recommendations (http://www.python.org/dev/peps/pep-0008/)

 - UpperCamelCase for class names

 - CAPITALIZED_WITH_UNDERSCORES for constants

 - lowercase_separated_by_underscores for other names (functions, variables etc)


"""


#----[ DESCRIPTION ]----
"""
Script with functions to perform different actions on interferometric/SD
radio fits cubes/arrays.

Need :  o scipy (and numpy)
        o matplotlib
        o mechanize module (standard in Python >2.6/2.7?)
        o urllib2 module (standard Python module)
        o beautiful soup/xml (xml standard in Python >2.6/2.7?)

"""


#----[ CHANGE LOG ]----
"""


* 2012 Oct 16 
    Moved data handling and manipulation to adacore.py
    move UI to separate module adavis.py
    and data I/O and handling to adacore.py DONE!

* 2012
    Moved constants and the radiative transfer scripts to separate
        modules (cgsconst[.py] and pyrt[.py]).

* 2012
    Update constants (AU, LY, etc)
       Created a Constant - object, so that they are more useful
            (can check unit and description of the constant).

* 2012?
    The "extent" keyword is NOW computed exactly. Uses crval/crpix
    FITS keyword to compute it, so it is truly from phase center.
     -> DONE but needs testing.

* 2011/2012
    fixed some ugly (but working) code in write_latex_table

* 2011
    Disabled the general parse_region function, it was wrong.
    use the method method instead or correct the general function?

* 2011?
    Common function for moment maps?
    Started for 0 and 1
    Needs update

* 2011?
    Add ALMA in get_telescope_diameter function


"""


#----[ BLUEPRINTS ]----
#Top of the list 

# TODO : naming scheme: 
#                       function names - whatit_does()
#                       constant - CONSTANT
#                       class name - ClassName

# TODO : move over to import MODULE
#        Then use it as MODULE.function(input)

# TODO method calc_offset to Fits object (calc offset from phase center given
#sexadecimal or sim coord input)

# TODO : extend the "box" plotting keyword so that it accepts arbitrary
#regions (just like the region keyword) <- perhaps sharing code is a good
#idea here(?)

#TODO : Cannot plot or deduce levels when some pixels have value 'nan'

#TODO : load frequency information from fits file if it is there, and create
#velocity array as well
#        -> more dynamic fits-file loading

#TODO : constants move to a module, update the code to use it.

#TODO : implement different lineID plot, make it more object oriented

#TODO : update the __str__ method of all the classes to print as much
#        info as possible

#TODO : Moment 2 maps
#        -Check others mom2 function (how to get line center?)
#        -Need MPFIT Gaussian fitting?

#TODO : Clean up code again, remove font handler, or inactivate it
#      All fonts should be Times New Roman, to big of a hassle to change
#       to Sans-serif font, for the tick labels mostly.
#       (works better in later matplotlib?)

#TODO : Check that the set_rc function is used in all functions

#TODO : Clean up moment 0/1(/2) function

#TODO : Change to handle DataObject in all functions.
#       (Use v_sys, dist, name from the Fits data object)
#        X Spectra
#        X Fit
#        X moment0
#        X moment1
#        O moment2
#        O The rest

#TODO : How to divide it into a runnable program
#        Perhaps use **kwargs i.e. dictionaries more which can
#        easily be passed along from e.g. OptParser

# TODO : The Error classes should take two input, which parameter
#        is raising the error and the reason

#------------------------------------------------------------------------
# Less pressing matters:

#TODO : P-V diagram - rotatable
#       what happens with the coordinates?
#       need to calc an offset, but call them x and y or something

# TODO : RMS and other units, e.g. when using SD data, Kelvin instead of Jy.

# TODO : axes_grid1 backwards compability? (for use in CASA)

# TODO : Tick locators are good for small regions/narrow spectra,
#        but not for big/wide
#       Perhaps change it with the 'box' keyword
#            -> e.g a 10th of the box keyword for major locators
#       Implemented a parameter locator=[major,minor] as a first step

# TODO : Function to bin in spatial (2D) regime (congrid map),
#         change ra/dec_delt etc

# TODO : Spectra class, use for fitting, lineid, binning

# TODO : Interactive multiple Gaussian fitting, i.e. click at start
#        end of each line, calculate moments within it

# TODO : Interactive splatsearch

# TODO : What does the VELREF keyword in the GILDAS fits header mean?
#       Check GILDAS manual

# TODO : Check Jes functions how they work and try to merge/replace them.

# TODO : Implement it as command line program
#        -> partially implemented, should work with spectra
#        Make sure all parameters are reachable through **kwargs
#        use Dictionary.get('attribute', ValueIfNoKey)

# TODO : Create a partition function class, to access partition functions
#       for different molecules.
#         o Example adavis.qrot.sulfurdioxide('18(1,2)-17(1,8)AE', Trot=170)
#         o Using the transition to identify it, and print shows the
#           available ones
#         o Store the Q values for different Tex
#         o Able to update from Q and Tex values

########################################################################
# IMPORTS

#~ from ..libs import cgsconst as _cgs  # if I only need it locally
                                        # import with preceeding "_"


########################################################################
# ASTRONOMY & ASTROPHYSICS definitions for figure sizes
# 255.76535 pt column width equals 88 mm (from aa.doc.pdf section 3.4)
# one inch in mm = 25.4 mm/inch
#one_col_fig_width_pt = 249.448819
#~ class FigSizes():
    #~ def __init__(self):
        #~ ONE_COL_FIG_WIDTH_MM = 88
        #~ TWO_COL_FIG_WIDTH_MM = 180
        #~ SIDE_CAPTION_FIG_WIDTH_MM = 120
        #~ INCHES_PER_PT = 1.0/72.27                               # convert pt to inches
        #~ INCHES_PER_MM = 1/25.4                                  # convert mm to inches
        #~ GOLDEN_MEAN = (5**0.5-1.0)/2.0                          # aesthetic ratio
        #~ ONE_COL_FIG_WIDTH = ONE_COL_FIG_WIDTH_MM*INCHES_PER_MM  # width in inches
        #~ ONE_COL_FIG_HEIGHT = ONE_COL_FIG_WIDTH*GOLDEN_MEAN      # height in inches
        #~ TWO_COL_FIG_WIDTH = TWO_COL_FIG_WIDTH_MM*INCHES_PER_MM
        #~ SIDE_CAPTION_FIG_WIDTH = SIDE_CAPTION_FIG_WIDTH_MM*INCHES_PER_MM
        #~ FIG_SIZE = [ONE_COL_FIG_WIDTH,ONE_COL_FIG_HEIGHT]
from ._figsizes import AandA as _AandA # only to local namespace

from ..adacore import *
########################################################################
# USEFUL STRINGS
_KMS = u"km\u00b7s\u207b\u00b9"

########################################################################
# GENERAL FUNCTIONS


def print_warning(s):
    import sys
    msg1 = stylify('WARNING:',f='b',fg='r')
    msg2 = stylify(' '+s,f='b',fg='k')
    sys.stderr.write(msg1 + msg2 +'\n')
def print_error(s):
    # Dont know if this works as intended
    import sys
    msg1 = stylify('ERROR:',f='b',fg='r')
    msg2 = stylify(' '+s,f='b',fg='k')
    sys.stderr.write(msg1 + msg2)
    sys.exit()
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
def parse_tick_font (font):
    """
    You have to add your formatter yourself with:
    labels = label_X.replace('X',data_fmt)
    where data_fmt  eg = '%g'
    and also:
    formatted_labels = FormatStrFormatter(labels)

    and lastly:
    ax.xaxis.set_major_formatter(formatted_labels)

    """
#    from scipy import array, where
#    if font.has_key('family'):
#        fformatters = array(['$\mathrm{X}$', '$\mathsf{X}$',
#                             '$\mathtt{X}$' , '$\mathit{X}$'])
#        ffamilies = array(['serif', 'sans-serif', 'monospace', 'cursive'])
#        label_X = fformatters[where(font['family']==ffamilies)][0]
#    else:
#        label_X = '$\mathrm{X}$'
    #return label_X
    return 'X'


########################################################################
# PLOTTING HELP FUNCTIONS
# to adavis.py
def draw_fov(ax, data):
    """
    Function to draw the field of view into the
    plot with axes 'ax'.
    Assumes that the axes data is set to arcsec
    """
    from matplotlib.patches import Circle
    cir = Circle((0,0), transform=ax.transData, fill=False, ec='k',
                 lw=1, ls='dashed', radius=data.fov/2)
    ax.add_patch(cir)
def draw_sizebar(ax, data, dist=220, au=200, loc=8):
    """
    distance in pc
    then theta = AU/(distance in pc)
    """
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    # draw horizontal bar with length of 10 in data coordinate
    # 10 arcs * 220pc = 2200
    asb = AnchoredSizeBar(ax.transData,
                            au/float(dist),
                            str(au)+" AU",
                            loc=loc,
                            pad=0.1, borderpad=0.5, sep=5,
                            frameon=False)
    ax.add_artist(asb)
    return asb
def draw_beam(ax, data, loc=3, bpad=0.2, ppad=0.15, box=True):
    """
    function that draws the beam
    the attributes data.bmin, .bmaj and .bpa must exist in the data class
    """
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse
    #from scipy import pi
    # the PA is calculated from N to E, but the plotting is relating the
    # minor axis to E-W line so just add -1*
    ae = AnchoredEllipse(transform=ax.transData,\
                            width=data.bmin,\
                            height=data.bmaj,\
                            angle=-1*data.bpa,\
                            loc=loc,\
                            pad=ppad,\
                            borderpad=bpad,\
                            frameon=box)

    ax.add_artist(ae)
def draw_highlight_box(ax, xpos, xwidth):
    """
    adds a highlight box that spans whole y-axes space and specified data
    coordinates in x
    """
    from matplotlib.patches import Rectangle
    from matplotlib.transforms import blended_transform_factory

    trans = blended_transform_factory(ax.transData, ax.transAxes)
    # We want x to be in data coordinates and y to
    # span from 0..1 in axes coords
    xpos -= xwidth/float(2)
    #print xpos, xwidth
    rect = Rectangle((xpos,0), width=xwidth, height=1,
                             transform=trans, color='yellow',
                             alpha=0.3)

    ax.add_patch(rect)
def put_line_indicator(ax, velocity, spect, xpos, text_string, \
            text_size=5, text_weight='extra bold', text_color='black',\
             lc='b', offset=0):
    """
    Put text in figure, the positon is in axes data coordinates
    additional kwargs are text_size, text_weight, text_color
    oset_bool is for lines at same places
    """
    # calculate position
    line_height = spect.max()*0.08
    try:
        maxval = spect[get_indices(velocity,[xpos-2,xpos+2])].max()
    except ValueError:
        maxval = spect.max()
    except IndexError:
        maxval = spect.max()
    line_ypos = maxval + maxval*0.15
    text_ypos = line_ypos + line_height + line_height*0.5
    text_ypos = text_ypos + len(text_string)*offset*3e-2
    # plot line
    ax.plot([xpos, xpos],[line_ypos,line_ypos+line_height],lc)
    # print text
    ax.text(xpos+0.1,text_ypos,\
    text_string,\
    size=text_size, weight=text_weight, color=text_color,\
    ha='center',va='bottom',rotation='vertical',\
    transform = ax.transData)
def set_rc(font={'family':'serif', 'serif': ['Times New Roman'],
        'size':10},
        quality=[300, 150], latex=True, **kwargs):
    from matplotlib import rc
    ################################
    # setting global rc properties #
    ################################
    rc('text', usetex=latex)
    #~ rc('mathtext',**{'fontset':'custom', 'default':'sf'})
    #~ rc('pdf', fonttype=42)
    rc('savefig', **{'dpi': quality[0]})
    rc('figure',**{'facecolor': '1', 'dpi': quality[1]})
    rc('image', **{'origin': 'lower'})
    rc('text.latex', preamble=[r'\renewcommand{\familydefault}{\sfdefault}'])
    #to set the global font properties
    rc('font', **font)
    # ticksize
    #rc('xtick',**{'minor.size':3, 'major.size':7})
    #rc('ytick',**{'minor.size':3, 'major.size':7})
    # linewidths
    rc('axes', linewidth=0.45)
    rc('patch', linewidth=0.45)
    rc('lines', linewidth=0.45, markeredgewidth=0.45)
def steppify(arr, isX = False, interval = 0):
    """
    Converts an array to double-length for step plotting
    Adapted from a script by Adam Ginsburg 
    (http://casa.colorado.edu/~ginsbura/)
    """
    from scipy import array
    if isX and interval==0:
        interval = (arr[1]-arr[0]) / 2.0
        newarr = array(zip(arr - interval, arr + interval)).ravel()
        return newarr
    else:
        newarr = array([arr,arr]).transpose().flatten()
        return newarr
def calc_levels(sigma, mini, maxi):
    from scipy import arange
    levels_neg = -1 * arange(sigma, abs(mini) + 2 * sigma, sigma)
    levels_pos = arange(sigma, maxi + 2 * sigma, sigma)
    return levels_neg, levels_pos
def plot_contours(ax, data, extent, levels, cp = 'b', cn='r', **kwargs):
    # plot contours with nice control over 
    # the dash size for negative values
    from scipy import where
    i_neg = data<0
    i_pos = data>=0
    d_neg = data.copy()
    d_neg[i_pos] = 0
    d_pos = data.copy()
    d_pos[i_neg] = 0
    
    c_pos = ax.contour(d_pos,
                extent = extent,
                levels = levels,
                colors = cp,
                **kwargs)
    c_neg = ax.contour(d_neg,
                extent = extent,
                levels = levels,
                colors = cn,
                **kwargs)    
    for c in c_neg.collections:
        c.set_dashes([(0, (2.5, 2.5))])
    return c_pos, c_neg
########################################################################
# DATA OUTPUT
def write_latex_table(self, filename):
    """
    What Héctor needs in a Latex table
    -From self.Spect.linedata:
    Molecule
    Transition
    Rest frequency
    Eu (K)
    -From self.Spect.v_cdelt & self.Spect.rms:
    Velocity Resolution
    RMS
    -From self.Spect.Fit:
    Peak intensity (line intensity)
    Line width
    Integrated intensity
    Frequency shift
    Velocity shift

    separated with &-sign
    line ended with \\ (and \n)
    comment, file/table info start with %
    """
    from scipy import arange,array

    # mol_data is a ndarray (4xNlines) that we need to loop through
    # when writing to the file
    (molecules, rest_frequencies, smu2s,
    eus, transitions, catalogs) = self.Spect.linedata
    # format the transitions array entries in latex fashion as:
    # $N1_{n1,n2} --N2_{n1,n2}$
    transitions = array(['$'+
                    i.replace('(','_{').replace(')','}')+
                    '$' for i in transitions])
    gaussian_peaks = self.Spect.Fit.params[::3]
    line_width_error = self.Spect.Fit.errors[2::3]
    line_pos_error = self.Spect.Fit.errors[1::3]
    v_res_kms = self.Spect.v_cdeltkms
    rms = self.Spect.rms
    #peak_gauss_intensities = self.Spect.Fit.params[::3]
    #peak_intensities = self.Spect.Fit.peak_intensities
    line_widths = self.Spect.Fit.line_widths
    int_intensities = self.Spect.Fit.gauss_intensities
    error_int_intensities = self.Spect.Fit.error_gauss_intensities
    #
    freq_shifts = self.Spect.Fit.freq_shifts
    vel_shifts = self.Spect.Fit.vel_shifts

    f = open(filename, 'w')
    f.write('% LaTex table of detected lines\n')
    f.write('% Created from the '
    'file {0}\n'.format(self.fitsfile.split('/')[-1]))
    f.write('% Intensities: I is the peak intensity (Gaussian fit)\n% '
    'while Int I is the integrated intensity (Gaussian fit).\n')
    cols =['Molecule', 'Transition', 'Rest F', 'Eu', 'V res', 'I',
           'FWHM', 'Int I', 'F shift', 'V shift', 'RMS']
    info = ('% Columns:\n%  {0:15}  {1:28}   {2:12} {3:6}   {4:6}   '
    '{5:8}   {6:15} {7:13} {8:13} {9:10} {10} \n'.format(*cols))
    f.write(info)
    unitinfo = [
        '---', '---', 'MHz', 'K', 'km/s', 'mK', 'km/s', 'K km/s', 'MHz',
        'km/s', 'mK']
    units = ('%  {0:17}  {1:28} {2:12} {3:8} {4:8} {5:10} {6:15} '
    '{7:15} {8:13} {9:9} {10} \n').format(*unitinfo)
    f.write(units)
    for i in arange(self.Spect.Fit.nfits):
        #Molecule Transition Rest frequency Eu (K) Velocity Resolution RMS
        #Peak intensity Line width Integrated intensity Frequency shift
        #Velocity shift
        data = (molecules[i], transitions[i], float(rest_frequencies[i])*1e3,
        float(eus[i]), float(v_res_kms),float(gaussian_peaks[i])*1e3,
        float(line_widths[i]),float(line_width_error[i]), float(int_intensities[i]),
        float(error_int_intensities[i]),float(freq_shifts[i])*1e3, float(vel_shifts[i]),float(line_pos_error[i]),
        float(rms)*1e3)
        line = ('{0:15} & {1:28} & {2:8.2f} & {3:8.3f} & {4:4.2f} & '
        '{5:5.1f} & ${6:4.1f} \\pm{7:4.1f}$ & {8:6.3f} \\pm{9:5.3f} & '
        '{10:5.1f}  & {11:5.1f} \\pm{12:3.2f} & {13:2.1f} \\\\ \n'.format(*data))
        f.write(line)
#
########################################################################
# ASCII TABLE HELP FUNCTIONS
# to adatools.py(?)
# table handling help functions. easily reads in
# simple ascii tables with # as comment character
def saveatbl(filename, dataList, names):
    """
    saves a list of data arrays ("dataList") in to a table with
    the coloumn names in "names" as an ASCII file with name
    "filename" (can be a path to a file as well)
    """

    def fileExists(f):
        try:
            file = open(f)
        except IOError:
            return False
        else:
            return True

    if type(dataList) != type([]):
        dataList = list(dataList) # if its not a list, make it one

    if len(dataList) != len(names):
        raise Exception('The number of column names does not match the number of data arrays')

    if fileExists(filename) == False:
        with open(filename,'w') as f: # force creating "filename"
            #first line, the column names
            out = "\t\t".join(['no']+names)
            f.write('#'+out+'\n')
            f.write('\n')

            for i in xrange(len(dataList[0])): # make the index as long as
                # the first array in the list
                out = "\t\t".join([str(i+1)]+[str(x[i]) for x in dataList])
                f.write(out+'\n') # add a line break

    elif fileExists(filename) == True:
        print(' ')
        print('FILE EXISTS ('+str(filename) +') - SKIPPING SAVE')
        print(' ')
def loadatbl(filename, dtype='float', rtype='array',sep=None, c_char=['#', '!', '|', '/']):
    """
    loads a list of data arrays in "filename", returns an array of the
    whole shebang (loads arrays saved with the savetable command). just do:
    (not the number column...)
    a = loadtable(filename)
    and then a[0] for first column, a[1] for second column and so on.
    """
    from scipy import array
    try:
        with open(filename,'r') as f:
            values = []
            for line in f:
                start_test = [line.startswith(x) for x in c_char]
                if True in start_test or not line.strip():
                    continue # skip lines that are comments and empty
                line = line.strip('\n')
                cols = line.split(sep)
                values.append(cols)
    except IOError:
        raise IOError('file ' +str(filename)+' does NOT exist...')
    except ValueError:
        raise ValueError('Trying to convert to '+str(dtype)+' while it is a string\
                        try to change it to \'str\'')
    if rtype=='array':
        try:
            return array(values,dtype=dtype).transpose()
        except ValueError:
            raise ValueError('Rows in text file have different number of columns')
    elif rtype=='native':
        return values
def infoatbl(filename, sep=None, c_char=['#', '!', '|', '/']):
    """
    just returns the lines with comments (ie the column names) from a *.atbl file
    """
    try:
        with open(filename,'r') as f:
            strings = []
            for line in f:
                start_test = [line.startswith(x) for x in c_char]
                if True in start_test:
                    strings.append(line.split(sep))
    except IOError:
        raise IOError('file' + str(filename)+'does not exist...')

    return strings
#
########################################################################
# UV DATA VISUALISATION FUNCTIONS
def plot_uvplane(self, units='klam'):
    import matplotlib.pyplot as pl
    from scipy import where, array
    pl.ion()
    pl.close(1)
    unit = array(['k$\lambda$', 'm', 'nsec'])
    unit_which = where(array(['klam', 'm', 'nsec']) == units)[0][0]
    x = [self.u_klam,
                self.u_m,
                self.u_nsec][unit_which]
    y = [self.v_klam,
                self.v_m,
                self.v_nsec][unit_which]

    fig = pl.figure(1, figsize=((_AandA.ONE_COL_FIG_WIDTH,_AandA.ONE_COL_FIG_HEIGHT*1.5)))
    ax = fig.add_subplot(111)
    ax.plot(x, y, '.b', ms=0.8)
    ax.set_aspect(1)
    ax.set_xlabel(unit[unit_which])
    ax.set_ylabel(unit[unit_which])
    fig.subplots_adjust(bottom=0.14,left=0.16)

def plot_uvdata(self, x='uvdist', xunit='klam', y='amp', overplot=0, avg=0, **kwargs):
    """
    average : average every X unit
              e.g. average=10 and unit klam, means average every
              10 kilolambda togeter.


    TODO : averaging not done correct.
    TODO : labels
    TODO : plot several things in one figure
    TODO : plot UV fits
    TODO : plot several UV sets in one figure
    TODO : function to average
    """
    import matplotlib.pyplot as pl
    pl.ion()
    from scipy import arange, where, zeros, array, sqrt
    set_rc()
    self.avg = avg
    #~ if avg>1:
        #~ if y == 'amp':
            #~ x,y = self.avgamp(avg)
    class Uvplt1: pass

    # first take care of X-axis
    Uvplt1.xunit = xunit
    xunits = array(['k$\lambda$', 'm', 'nsec'])
    xunit_which = where(array(['klam', 'm', 'nsec']) == xunit)[0][0]
    Uvplt1.xunitstr = xunits[xunit_which]
    Uvplt1.x = [self.uvdist_klam,
                self.uvdist_m,
                self.uvdist_nsec][xunit_which]
    # get the Y-axis
    ytype_which = where(array(['amp','pha']) == y)[0][0]
    Uvplt1.y = [self.amplitude, self.phase][ytype_which]
    Uvplt1.yunitstr = array(['Amplitude [Jy]','Phase [Degrees]'])[ytype_which]
    if avg:
        x = arange(Uvplt1.x.min()+avg/2.,\
                              Uvplt1.x.max()-avg/2.,\
                              avg)
        y = zeros(len(x))
        yerr = zeros(len(x)); j = 0
        for i in x:
            ipos = where((Uvplt1.x>=(i-avg/2))*\
                (Uvplt1.x<(i+avg/2)))[0]

            #~ y[j] = Uvplt1.y[ipos].mean()
            y[j] = sqrt(self.re[ipos].mean()**2 + self.im[ipos].mean()**2)
            yerr[j] = Uvplt1.y[ipos].std()
            j += 1
        Uvplt1.x = x
        Uvplt1.y = y
        Uvplt1.yerr = yerr
    if overplot:
        class Uvplt2: pass
        Uvplt2.x = [overplot.uvdist_klam,
                    overplot.uvdist_m,
                    overplot.uvdist_nsec][xunit_which]
        Uvplt2.y = [overplot.amplitude, overplot.phase][ytype_which]
        if avg:
            x = arange(Uvplt2.x.min()+avg/2.,\
                                  Uvplt2.x.max()-avg/2.,\
                                  avg)
            y = zeros(len(x)); j = 0
            yerr = zeros(len(x))
            for i in x:
                ipos = where((Uvplt2.x>=(i-avg/2))*\
                    (Uvplt2.x<(i+avg/2)))[0]
                y[j] = Uvplt2.y[ipos].mean()
                yerr[j] = Uvplt2.y[ipos].std()
                j += 1
            Uvplt2.x = x
            Uvplt2.y = y
            Uvplt2.yerr = yerr
    print Uvplt1.x.shape
    print Uvplt1.y.shape
    class Plot_data: pass
    Plot_data.Uvplt1 = Uvplt1
    pl.close(1)
    fig = pl.figure(1, figsize=((_AandA.ONE_COL_FIG_WIDTH*1.7,_AandA.ONE_COL_FIG_HEIGHT*1.5)))
    ax = fig.add_subplot(111)
    if y == 'amp' and avg == 0:
        err = self.error
    else:
        err = zeros(Uvplt1.y.shape)
    ax.errorbar(Uvplt1.x, Uvplt1.y, err, fmt='.',color='b')
    if overplot:
        ax.plot(Uvplt2.x, Uvplt2.y, '.r', ms=2)
        Plot_data.Uvplt2 = Uvplt2
    ax.set_xlabel('UV distance ['+Uvplt1.xunitstr+']')
    ax.set_ylabel(Uvplt1.yunitstr)
    fig.subplots_adjust(bottom=0.15,right=0.97,top=0.9)
    ax.set_title('UV Plot')
    self.Plot_data = Plot_data
    if 'send' in kwargs:
        return ax


########################################################################
# CUBE DATA VISUALISATION FUNCTIONS
# to adavis.py
def plot_spectrum (self,
    region=[0,0,0,0],
    source = dict(),
    show_freq=False,
    font={'family':'serif', 'serif': ['Times New Roman'],
    'size':8},
    send=False,
    quality=[300, 300],
    plot_adjust= [0.15, 0.17, 0.98, 0.95],
    axspace = [1., 1., 1., 1.],
    ylimits=None,
    telescope=None,
    fsize=(_AandA.FIG_SIZE),
    binning=1,
    bintype='mean',
    **args):
    """
    Plot the spectrum of a DataObject

    List of available parameters
    and their default value
    ----------
    DataObject
    chvals = None
    nvals = None
    region = [0,0,0,0]
    source = dict(v_sys=0)
    show_freq = False
    font = {'family':'serif', 'serif': ['Times New Roman'],'size':8}
    bin = 1
    bintype ='mean'
    linefit = dict(type=None, params=[(0.09, 7.3, 3.6)],
        guess=False, interactive=False, fixlist=None, error=None,
        limmin=None, minpar=None, limmax=None, maxpar=None, tie=None)
    send = False
    quality = [300, 300]
    plot_adjust = [0.15, 0.17, 0.98, 0.95]
    lines = None
    axspace = [1.01, 1.01, 1.01, 1.05]
    ylimits = None
    telescope = None
    fsize = (FIG_SIZE)

    Description of parameters
    ----------
    DataObject :
        Optional

    chvals :
        Default : Default = None

    nvals :
        Default : None

    region :
        Default : [0,0,0,0]

    source :
        Default : dict(v_sys=0)

    show_freq :
        Default : False

    font :
        Default : {'family':'serif', 'serif': ['Times New Roman'],'size':8}

    bin :
        Default : 1

    bintype :
        Default : 'mean'

    linefit :
        Default : dict(type=None, params=[(0.09, 7.3, 3.6)],
        guess=False, interactive=False, fixlist=None, error=None,
        limmin=None, minpar=None, limmax=None, maxpar=None, tie=None,
        lineid=0/dict(nfwhm=1.5))

    send :
        Default : False

    quality : [int. int]
        Default : [300, 300]
        Set the quality of the saved figure (a) and the displayed figure (b), where [a, b]

    plot_adjust :
        Default : [0.15, 0.17, 0.98, 0.95]
        Set the left, bottom, right, top of the axes

    lines :
        Default : None
        draw lines with name 'a' & 'c' and rest frequencies b & d.
        if obj_dict=dict(v_sys=x), or DataObject.v_sys is not null, then it will shift the positions
        of the lines accordingly.

    axspace :
        Default : [1.01, 1.01, 1.01, 1.05]

    ylimits :
        Default : None

    fsize :
        Default : (FIG_SIZE)

    lineid : boolean
        Default : False
        Use splatsearch module to search the splatalogue for lines?


    TODO : What is the RMS when binning and giving a area_region>1?
                    # RMS is not accounting for binning in vel. axis
                    # RMS is not correct when giving a are_region>1
            ie does the end unit change? yes-> fix it
                    # change x1,x2,y1,y2 to quarter region (line 2350)
                    # change so that when binning, the rms i calculated

    TODO : Remove the most weird font-settings

    TODO : RMS units, e.g. when using SD data, Kelvin instead of Jy...

    TODO : The RMS calculation, always for both binned and unbinned data

    #TODO : Fix the axspace implementation, must be a better way

    #TODO : Perhaps be able to supply more for the figure, e.g. figure number



    Remember! : hasattr(data,'vsys')

    """
    # imports
    print 'importing...'
    from scipy import array, where, median, std, sign, arange, alen, vstack, \
                    concatenate, sqrt, log10, exp, log, ceil, floor, diff, \
                    flipud, pi, nan
    #from string import lower
    import matplotlib.pyplot as pl
    from mpl_toolkits.axes_grid1 import AxesGrid
    #from matplotlib.patches import Circle
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from matplotlib import cm, rc, rc_params
    import adapy
    print 'done'
    ####
    #### PARSE INPUT
    ####

    # formatting
    data_fmt = '%g'         # for normal y and x axis
    freq_data_fmt = '%5.2f' # for the frequency array
    label_X = parse_tick_font(font)
    tick_label_formatter = FormatStrFormatter(label_X.replace('X',data_fmt))
    freq_label_formatter = FormatStrFormatter(label_X.replace('X',freq_data_fmt))

    # set the subplot_adjust parameters, if not given a standard set will
    # be used
    pl_left, pl_bottom, pl_right, pl_top = plot_adjust
    # load input parameters
    #~ inputdict = dict(region=region)
    # create the spectrum object
    inputargs = args
    Spect = adapy.Spectrum(self, region=region, **inputargs)
    if binning > 1:
        if bintype == 'resample':
            Spect.bin_spectrum(binning=binning, bintype=bintype)
        else:
            Spect.bin_spectrum(binning=binning)
    # calculate RMS
    #~ print args
    if 'nvals' in args:
        Spect.calc_rms(self, nvals=args['nvals'])
    else:
        pass
    if 'linefit' in args:
        Spect.fit_lines(**args)
    # fit data?
    #~ if 'linefit' in args:
        #~ Spect.fit

    # 1. Extract spectra
    # 2. BIN data
    # 3. Calculate RMS
    # 4. Plot spectrum
    # 5. If fitkey exist : Fit data
    # 6. If lineID in fitkey OR linenameskey exists : Line ID

    # OLD CODE IN COMMENT BELOW
    """
    #self.freqarr = calc_frequency(self.v_arr,self.restfreq)
    #
    # parse the region parameter
    # does it exist?
    # now parse the region keyword
    x1,x2,y1,y2 = parse_region(self, region)
    #
    # now start the plotting
    """

    ###################################
    #      extract the spectra        #
    ###################################
    """
    if self.datatype[0] != 'SDSPECT': # if it is a cube
        area_region = ((y2-y1)*(x2-x1))
        spect = (self.d[:,y1:y2,x1:x2].sum(axis=1).sum(axis=1))/float(area_region)
    elif self.datatype[0] == 'SDSPECT': # if it is a single dish spectra
        area_region = 1
        print stylify("SD-SPECTRUM - region keyword not doing anything.",fg='r')
        spect = self.d
    """
    ####################################
    #          binning of data         #
    ####################################
    """
    # it has to be an integer, for now at least
    if binning != 1:
        binning = int(binning)
        if lower(bintype) == 'resample':
            from congridding import congrid
            # congridding, proper resampling of data
            spect = congrid(spect,(alen(spect)/binning,),centre=True, method='neighbour')
            velocity = congrid(velocity,(alen(velocity)/binning,))
            #
            velocity_delta = self.v_cdeltkms*binning
        elif lower(bintype) == 'mean':
            if alen(spect)%binning!=0:
                print 'Bin has to be evenly devide the number of channels: %d' % alen(spect)
                sysexit()
            #  Old method - simple binning, just average
            indices = arange(0,alen(spect),binning)
            spect = array([spect[x:x+binning].sum(axis=0)/binning for x in indices])
            velocity = array([velocity[x:x+binning].sum(axis=0)/binning for x in indices])
            #
            velocity_delta = self.v_cdeltkms*binning
        elif binning == 0 or binning <0:
            print stylify("\nERROR:\n Variable \"bin\" has to be 1 for no binning, or above 1 \n\
            for the number of channels to bin")
        # print out information about the binning
        print '='*40
        print ' '*11,"Binning of data\n"
        print "No channels to bin : %d" % binning
        print "Velocity step : %f" % velocity_delta
        if bintype=='mean':
            print 'Type of binning : Simple mean over selected no. bin channels'
        elif bintype=='resample':
            print 'Type of binning : Resampling - 1D interpolation'
    """

    ####################################
    #       noise calculation          #
    ####################################
    """
    if 'nvals' in args: #nvals!=None:
    #~ if nvals!=None:
        print '='*40
        print ' '*11,'Noise statistics\n'
        # calculate the rms from the channels in the spectra
        # accounts for it even if it is binned
        #
        # image rms
        # change x1,x2,y1,y2 to quarter region
        # change so that when binning, the rms i calculated
        # x1,x2
        if self.datatype[0] == 'SDSPECT':
            rms = sqrt(((self.d[get_indices(velocity,args['nvals'])])**2).mean()/float(binning))
        else:
            zlen, ylen, xlen = self.d.shape
            ydelt = ylen/6
            xdelt = xlen/6
            i1,i2 = xlen/2-xdelt, xlen/2+xdelt
            j1,j2 = ylen/2-ydelt, ylen/2+ydelt
            rms = sqrt(((self.d[get_indices(velocity,args['nvals']),j1:j2,i1:i2])**2).mean()/float(binning))
        rms_mjy = rms*1e3
        #rms_0 = sqrt(((spect[get_indices(velocity,args['nvals'])])**2).mean())
        #rms_2 = sqrt(((self.d[get_indices(velocity,args['nvals']),:,:])**2).mean())

        #rms_0= rms/sqrt(abs(self.v_cdeltkms))
        #print rms_0
        #print 'rms_0 =', rms_0*1e3
        #print 'rms_2 =', rms_2*1e3
        #
        # the sensitivity
        s = rms/sqrt(abs(velocity_delta))
        s_mjy = s*1e3
        # the channels used
        ind = get_indices(velocity,args['nvals'])
        ind1, ind2 = ind.min(), ind.max()
        print u'RMS \t\t: %2.3f mJy\u00b7beam\u207b\u00b9\u00b7channel\u207b\u00b9' % rms_mjy
        print u' Sensitivity \t: %2.3f mJy\u00b7beam\u207b\u00b9\u00b7km\u207b\u00b9\u00b7s' % s_mjy
        print u' Channels used \t: %d, %d (%s km\u00b7s\u207b\u00b9)' % (ind1, ind2, str(args['nvals']))
        print u' Region \t: %s arcsec' % (str(region))
        print u' Channel width \t: %2.3f km\u00b7s\u207b\u00b9' % abs(velocity_delta)
    """
    ####################################
    #           fitting data           #
    ####################################
    """
    if linefit['type']=='gauss':
        print '='*40
        print ' '*11,"Line fitting\n"
        #
        fit = linefit.copy()
        if not linefit.has_key('error'):
            #~ if nvals!=None:
            if 'nvals' not in args:
                fit['error'] = rms
                print 'No error supplied, but nvals given and rms calculated.\n\
                using rms = %2.3f Jy beam-1 channel-1' % rms
            elif args['nvals']==None:
                errmsg='You have to supply an error to the fit, they\'re kind of important you know.'
                print stylify(errmsg)
                return ; sysexit()

        fwhmfromsig = 2*sqrt(2*log(2)) # the constant
        fwhm = lambda x: fwhmfromsig*x
        sigma = lambda x: x/fwhmfromsig

        # fits 1D gaussian(s) to spectral data
        if 'chvals' in args: # if we supplied limits for the fit
            ch = where((velocity>args['chvals'][0])*(velocity<args['chvals'][1]))[0]
            Fx = spect[ch]
            X = velocity[ch]
        else: # else use the eveything
            Fx = spect
            X = velocity
        #
        #
        p = fit['params'] # parameters for the fit
        #
        if not fit.has_key('limmin'):
            fit['limmin'] = None
            fit['minpar'] = None
        if not fit.has_key('limmax'):
            fit['limmax'] = None
            fit['maxpar'] = None
        if not fit.has_key('fixlist'):
            fit['fixlist'] = None
        if not fit.has_key('tie'):
            fit['tie'] = None
        #
        from time import time
        t1= time()
        #
        fitting_results = fit_gauss1d((X,Fx), params=p, fixlist=fit['fixlist'],
                minbool=fit['limmin'], minpar=fit['minpar'],
                maxbool=fit['limmax'], maxpar=fit['maxpar'],
                err=fit['error'], tie=fit['tie'], verbose=0, full_output=1)

        if fitting_results==None:
            print(stylify('\n\n No fitting done...',f='b',fg='r'))
        elif fitting_results!=None:
            params, errors, chi2, mp = fitting_results
            #
            print ' Done in %2.3f seconds \n' % (time()-t1)
            #
            print ' Number of fits : ', alen(params)/3
            print ' Fit status : ', mp.status, '(if 0, it should have halted)'
            print ' Chi2 : {0}, reduced : {1}\n'.format(chi2,chi2/float(len(Fx)))
            # now, parse output of fitting and print it out on screen
            j = 1
            line_widths = []
            frequencies = []
            #
            # temporary code, move these attributes and methods over to a class (Spectrum?)
            # both the fitting and the initialisation, each spectra class can be fitted individually
            #
            # move fitting up before plotting
            #
            #
            class Tmp: pass
            #Tmp.fitting_results = fitting_results
            Tmp.params = params
            Tmp.errors = errors
            Tmp.chi2 = chi2
            Tmp.nfits = alen(params)/3
            #
            #
            #
            for i in arange(0,len(params),3):
                # add 1 because channel 1 is in pos 0

                half_fwhm = params[i+2]/2.
                fwhm = params[i+2]
                line_widths.append(fwhm)
                # first figure out the extent of the gaussian (the line)
                # jump half a channel down and up so that it finds the correct channels
                lower_half, upper_half = (params[i+1] + array([-1,1])*half_fwhm)
                lower,upper = (params[i+1] + array([-1,1])*fwhm)
                #channels = where((velocity>lower)*(velocity<upper))[0]+1
                channels_half_nobin = get_indices(self.v_arr,[lower_half,upper_half])
                channels_nobin = get_indices(self.v_arr,[lower,upper])
                channels_half = get_indices(velocity,[lower_half,upper_half])
                channels = get_indices(velocity,[lower,upper])
                #draw_highlight_box(ax_kms, params[i+1], params[i+2]*3)
                # apply v_sys correction, so that we use v_sys for estimating the correct
                # frequency for the line, especially importat when using line-identification
                frequency = calc_frequency(params[i+1]-v_sys,self.restfreq/1e9)
                frequencies.append(frequency)
                print  'Fit number : %i' % j
                print  ' Intensity : %2.4f \t=calculate it=' % (sqrt(2*pi)*sigma(params[i+2])*params[i]) # the area under the 1D Gaussian
                print u' Amplitude : %2.3f (\u00b1%2.3f) \t Jy\u00b7Beam\u207b\u00b9' % (params[i],errors[i])
                print u' Position  : %2.3f (\u00b1%2.3f) \t km\u00b7s\u207b\u00b9' % (params[i+1],errors[i+1])
                print u' Width     : %2.3f (\u00b1%2.3f) \t km\u00b7s\u207b\u00b9 (FWHM, \u03c3=%2.3f)' % (params[i+2],errors[i+2], sigma(params[i+2]))
                print  ' Frequency : %3.9f GHz (v_sys corrected)' % frequency
                print u' FWHM     : %d, %d (%d) ([%.2f, %.2f] km/s)'% (channels_half_nobin.min(), channels_half_nobin.max(),(channels_half_nobin.max()-channels_half_nobin.min()+1), lower_half, upper_half)
                print u' \u00b1FWHM   : %d, %d (%d) ([%.2f, %.2f] km/s)' % (channels_nobin.min(), channels_nobin.max(), (channels_nobin.max()-channels_nobin.min()+1), lower,upper)
                if bin!=1:
                    print u' Rebinned channels :\n \t FWHM width  : %d, %d (\u00b1FWHM)' % (channels_half.min(), channels_half.max())
                    print  u' \t \u00b1FWHM width : %d, %d (%d)  \n' % (channels.min(), channels.max(),(channels.max()-channels.min()+1))
                j+=1
            #
            line_widths = array(line_widths)
            print 20*'- '
            print u'Mean FWHM : %2.1f \u00b1%2.2f km\u00b7s\u207b\u00b9' % (line_widths.mean(),line_widths.std())

            #### old code block
            #~ if send: # does this clause do anything? sending this way earlier...
                #~ j = 1
                #~ f = []
                #~ for i in arange(1,len(params),3):
                    #~ nu = calc_frequency(params[i]-v_sys,self.restfreq/1e9)
                    #~ f.append(nu)
                    #~ print '%3.9f' % nu
                    #~ j+=1
            ### end of old code block
            # draw the fit(s) into the figure
            # X is the velocity array
            # xarr has more 3 times more datapoints than velocity array
            # the plotted lines looks smoother and nicer that way
            xarr = arange(X[0],X[-1],(diff(X)[0]/3))
            #~ for i in arange(0,len(params),3):
                #~ lower,upper = (params[i+1] + array([-1,1])*fwhm*4)
                #~ channels = get_indices(xarr,[lower,upper])
                #~ ax_kms.plot(xarr[channels], gauss1d(xarr[channels],params[i:i+3]))
            #~ ax_kms.plot(xarr, gauss1d(xarr,params), color='0.2', lw=1, alpha=0.6)
        #
            if linefit.has_key('lineid'):
                if linefit['lineid']:
                    # later when the args is implemented...
                    # checks that we have done a fit first
                    #if not kwargs['linefit']:
                    #    print 'If you want lineid you need linefit'
                    #    raise ParError(kwargs['lineid'], kwargs['linefit'])
                    import splatsearch as spl
                    print 'Trying to indentify candidates for the fitted lines.'
                    frequency_pairs = []
                    for i in arange(0,len(params),3):
                        # correct for v_sys, to get correct frequency for the
                        # correct frequency range in the splatalogue search
                        # when calculating the "correct" frequency, we subtract v_sys
                        vel_lower, vel_upper = (params[i+1]-v_sys + array([-1,1])*params[i+2]*1.5)
                        # frequency increases when velocity decreases...
                        freq_lower = calc_frequency(vel_upper,self.restfreq/1e9)
                        freq_upper = calc_frequency(vel_lower,self.restfreq/1e9)
                        frequency_pairs.append([freq_lower,freq_upper])
                    list_of_species = []
                    list_of_frequencies = []
                    number = 1
                    for i in arange(len(frequency_pairs)):
                        df=8e-3 # range to find line
                        CSI = "\x1b["
                        start =CSI+'1m'+CSI+'32m'+CSI+'40m'
                        end = CSI+'m'
                        print '\n'+start+'Line number : '+str(number)+'\t\t\t\t'+end
                        print 'Frequency : %f  GHz' % (frequencies[i])
                        result = spl.splatsearch(freq=frequency_pairs[i], send=1, display=1, linelist=['jpl','cdms'], e_to=500)
                        if result!=None:
                            species, freq = result[1],result[3]
                            for i in arange(len(freq)):
                                list_of_species.append(species[i])
                                list_of_frequencies.append(freq[i])
                            #~ for i in arange(len(freq)):
                                #~ if i>0 and freq[i]!=freq[i-1]: # remove duplicates
                                    #~ list_of_species.append(species[i])
                                    #~ list_of_frequencies.append(freq[i])
                                #~ elif i==0:
                                    #~ list_of_species.append(species[i])
                                    #~ list_of_frequencies.append(freq[i])
                                #~ else:
                                    #~ pass

                        number+=1
                    # done now define the linelist
                    lines=[list_of_species,list_of_frequencies]
    """
    ####################################
    #            return data           #
    ####################################
    if send:
        return Spect
    #    return Spect
        #~ if nvals!=None and linefit['type']!=None: # send everything!
        #~ if 'nvals' in args and 'type' in args: # send everything!
            #~ txt = '\n sending you spectra, v_arr, data, noise-spectra'
            #~ print stylify(txt,fg='g')
            #~ return spect, velocity, self, Tmp
        #~ elif linefit['type']=='gauss': # gaussian fitting
            #~ txt = '\n sending you spect, v_arr, data'
            #~ print stylify(txt,fg='g')
            #~ return spect, velocity, self, Tmp
        #~ elif 'nvals' in args and 'type' not in linefit: # no fit but noise
            #~ txt =  '\n sending you spect, v_arr, data, noise-spectra'
            #~ print stylify(txt,fg='g')
            #~ return spect, velocity, self, spect[get_indices(velocity,args['nvals'])]
        #~ else: # well non of them are supplied
            #~ txt =  '\n sending you spectra, v_arr, data'
            #~ print stylify(txt,fg='g')
            #~ return spect, velocity, self
    # use set_rc here!
    #set_rc

    ################################
    # setting global rc properties #
    ################################
    #rc('savefig', **{'dpi': quality[0]})
    #rc('text', usetex=True)
    #rc('figure',**{'facecolor': '1', 'dpi': quality[1]})
    #to set the global font properties
    #rc('font', **font)
    #rc('axes',linewidth=1)
    #rc('lines', linewidth=0.5, markeredgewidth=0.8)
    #rc('patch',linewidth=0.5)
    # ticksize
    #rc('xtick',**{'minor.size':2, 'major.size':4, 'major.pad': 3})
    #rc('ytick',**{'minor.size':2, 'major.size':4, 'major.pad': 1})
    set_rc()
    ### done!
    pl.ion()
    pl.close()
    fig = pl.figure(1,figsize=fsize)
    #fig = pl.figure(1,figsize=(10.5,4.5))
    fig.clf()
    ax_kms = fig.add_subplot(111)
    """
    PLOTTING NOT ACCOUNTING FOR BINNING!!!
    Just overwrite old value and add a flag binned True or False
    DONE?

    """
    # now correct so that negative values are to the left
    if Spect.v_cdelt<0:
        ax_kms.step(flipud(Spect.v_arr), flipud(Spect.d), 'k',  where='mid')
    elif self.v_cdelt>=0:
        ax_kms.step(Spect.v_arr, Spect.d, 'k', where='mid')
    # if we want som space in the figure
    cx1 = axspace[0]
    cx2 = axspace[1]
    cy1 = axspace[2]
    cy2 = axspace[3]
    xmin, xmax = round(Spect.v_arr.min()), round(Spect.v_arr.max())
    ymin, ymax = round(Spect.d.min(),3), round(Spect.d.max(),3)
    sxmin,sxmax = sign(xmin)<0, sign(xmax)<0
    symin,symax = sign(ymin)<0, sign(ymax)<0
    xmin, xmax = xmin*[1/cx1,cx1][sxmin], xmax*[cx2,1/cx2][sxmax]
    ymin, ymax = ymin*[1/cy1,cy1][symin], ymax*[cy2,1/cy2][symax]

    if 'nvals' in args:
        ax_kms.plot([xmin, xmax],[Spect.rms,Spect.rms],'b--', alpha=0.7)
        ax_kms.plot([xmin, xmax],[-Spect.rms,-Spect.rms],'b--', alpha=0.7)
    #
    # plot dotted lines at y=0
    # and x=0 if v_sys is in the range that we plot
    if (Spect.v_arr<self.v_sys).any()*(Spect.v_arr>self.v_sys).any() or (Spect.v_arr==self.v_sys).any():
        ax_kms.plot([xmin-10,xmax+10],[0,0],'g:',[self.v_sys,self.v_sys],[ymin-1,ymax+1],'g:')
    else:
        ax_kms.plot([xmin-10,xmax+10],[0,0],'g:')
    #####################
    #   Ticklocations   #
    #####################
    #xmajorLocator = MultipleLocator(5)
    #xminorLocator = MultipleLocator(2.5)
    #ymajorLocator = MultipleLocator(0.05)
    #yminorLocator = MultipleLocator(0.025)

    #ax_kms.xaxis.set_major_locator(xmajorLocator)
    #ax_kms.xaxis.set_minor_locator(xminorLocator)
    #ax_kms.yaxis.set_major_locator(ymajorLocator)
    #ax_kms.yaxis.set_minor_locator(yminorLocator)
    #
    #~ return Spect
    # just so we se something! :)
    ax_kms.set_ylim(ymin,ymax); ax_kms.set_xlim(xmin,xmax);pl.draw()
    if 'linefit' in args:
        #Spect.fit_lines(linefit=args['linefit'])
        # what if the fitting fails, no "Fit" subclass is created!
        j=1
        for i in arange(0,len(Spect.Fit.params),3):
            lower,upper = (Spect.Fit.params[i+1] + array([-1,1])*Spect.Fit.params[i+2]*4)
            channels = get_indices(Spect.Fit.xarr,[lower,upper])
            ax_kms.plot(Spect.Fit.xarr[channels], gauss1d(Spect.Fit.xarr[channels],Spect.Fit.params[i:i+3]))
            ax_kms.text(Spect.Fit.params[i+1],
                        -5*Spect.rms,
                        str(j),
                        verticalalignment='top',
                        horizontalalignment='center')
            j += 1
            #draw_highlight_box(ax_kms, Spect.Fit.params[i+1], Spect.Fit.params[i+2]*2)
        ax_kms.plot(Spect.Fit.xarr, gauss1d(Spect.Fit.xarr,Spect.Fit.params), color='0.2', lw=1, alpha=0.6)
        pl.draw()
        if args['linefit'].has_key('lineid'):
            if 'lines' in args:
                print_warning('You have supplied a lines keyword, this will be ignored now')
            #if args['linefit']['lineid'].has_key('writetofile'):
            #    writetofile = args['linefit']['lineid']['writetofile']
            #else:
            #    args['linefit']['lineid']['writetofile'] = 0
            Spect.identify_lines(**args['linefit']['lineid'])
    if hasattr(Spect,'lines'):
        print u'Marking the lines, using %2.2f km\u00b7s\u207b\u00b9' % Spect.v_sys
        #~ print Spect.lines
        lines = parse_linelist(Spect.lines)
        #~ print lines
        #~ print Spect.lines, lines
        # add v_sys so that it plots it at the right v_sys
        #~ if len(lines) == 2:
            #~ lines = [lines[0][0], lines[1][0]]
        #~ print Spect.lines
        #~ print lines
        #~ print [float(lines[i+1]) for i in arange(0, len(lines), 2)]
        #~ print lines
        v = array([calc_vlsr(float(lines[1][i])*1e9,self.restfreq)+self.v_sys for i in arange(0, len(lines[1]))])
        colors = ['k']*len(lines[1])
        x = 0
        for i in arange(len(lines[1])):
            # only check for lines behind the new line
            # change this so that it takes the number of channels away
            # check lineid_plot.py script for inspiration
            no_dbl = len(where(array(v)[:i+1].round(0) == round(v[i],0))[0])
            put_line_indicator(ax_kms, Spect.v_arr, Spect.d, v[i], lines[0][i],lc=colors[x], offset=no_dbl, text_color=colors[x])
            #put_line_indicator(ax_kms, velocity, spect, v[j], lines[i],lc='k', offset=no_dbl)
            # verbose output
            #print u'Line %1d : %2.2f\t km\u00b7s\u207b\u00b9 (%2.3f)' % (j+1,v[j],v[j]+obj_dict['vsys'])
            x+=1
        ax_kms.set_ylim(ymin,ymax+ymax*0.4)
    else:
        ax_kms.set_ylim(ymin,ymax)
    ax_kms.set_xlim(xmin,xmax)
    #~ if lines!=[]:
        #~ ax_kms.set_ylim(ymin,ymax+ymax*0.4)
    #~ else:
        #~ ax_kms.set_ylim(ymin,ymax)
    if ylimits!=None:
        ax_kms.set_ylim(ylimits)
    #ax_kms.set_title(self.obj)
    ax_kms.text(0.07,0.87 ,self.obj, transform=ax_kms.transAxes)
    ax_kms.set_xlabel('$v$ [km s$^{-1}$]')
    #~ if self.unit=='Jy/beam':
        #~ ax_kms.set_ylabel('$I$ [Jy beam$^{-1}$]')
    #~ else:
        #~ ax_kms.set_ylabel('$I$ ['+self.unit+']')
    ax_kms.set_ylabel('$I$ ['+self.unit+']')
    # to show the channel numbers
    if show_freq==True:
        # create the frequency axis
        ax_hz = ax_kms.twiny()
        # change the font of the ticks
        # this messes things up!
        #ax_hz.set_xticklabels(ax_hz.get_xticks(),ffont)
        #ax_hz.set_yticklabels(ax_hz.get_yticks(),ffont)

        # to load special fonts
        # rc('text.latex',preamble='\usepackage[bitstream-charter]{mathdesign}')
        ax_hz.xaxis.set_major_formatter(freq_label_formatter)
        # removed the set_minor_locator settings below because there
        #   where to many ticks (wide spectra)
        #ax_hz.xaxis.set_minor_locator(MultipleLocator(0.001))

        x_1, x_2 = ax_kms.get_xlim()
        # i want the frequency in GHz so, divide by 1e9
        ax_hz.set_xlim(calc_frequency(x_1,self.restfreq/1e9), calc_frequency(x_2,self.restfreq/1e9))
        pl_top -= 0.1
        pl.xticks(rotation=40,horizontalalignment ='left')
    #
    #elif lines!=[] and obj_par['vsys']==0:
    #    print('please, input a v_sys!=0')
    #    return ; sysexit()
    #ax_kms.xaxis.set_major_formatter(tick_label_formatter)
    #ax_kms.yaxis.set_major_formatter(tick_label_formatter)

    #ax_kms.set_xticklabels(ax_kms.get_xticks(),ffont)
    #ax_kms.set_yticklabels(ax_kms.get_yticks(),ffont)
    ax_kms.minorticks_on()
    #for tick in ax_kms.xaxis.get_major_ticks():
    #    tick.label.set_family('Arial')
    #pl.show()
    #pl.draw()
    pl.subplots_adjust(left=pl_left, bottom=pl_bottom, right=pl_right, top=pl_top)
    #~ return spect, velocity, self, Tmp
    #
    #
    # LASTLY modify the indata object so that it has the spectrum (and fits)
    self.Spect = Spect

# to adavis.py
def plot_moment_map(self,
                chvals=None,
                nvals=None,
                region=[0,0,0,0],
                filled=True,
                rms_area = [0, 0, 10],
                type=0,
                box=[0,0],
                nx=6,
                nsig=3,
                nsjump=2,
                source = dict(v_sys=0, dist=250),
                font={'family':'serif', 'serif': ['Times New Roman'],
                'size':8},
                fit = dict(params = None, interactive = False),
                send=False,
                quality=[300, 300],
                cbar=True,
                colormap=True,
                plot_adjust= [0.07, 0.06, 0.82, 0.94],
                cpeak=[0,0,'k'],
                ccol='k',
                sbar=dict(au=200),
                locators = [2,1],
                telescope=None,
                negcontours=True,
                fsize=(_AandA.TWO_COL_FIG_WIDTH, _AandA.TWO_COL_FIG_WIDTH*0.7),
                **kwargs):
    # imports
    #import scipy as sp
    from scipy import array, where, median, std, sign, arange, alen, vstack, \
    concatenate, sqrt, log10, ones, sqrt, flipud
    #import pyfits as pf
    from time import sleep
    import matplotlib.pyplot as pl
    #from matplotlib.patches import Circle
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from matplotlib import cm, rc
    #


    """
    1. get data
    2. calc moments
    3. calc stats
    4. calc levels
    5. display
    """
    ################################
    # setting the tick label font  #
    ################################
    # we are using latex formatting here!
    # following sets the ticklabel format string
    data_fmt = '%g'                 # for normal y and x axis
    cbar_data_fmt = '%2.2f'         # for normal y and x axis
    label_X = parse_tick_font(font)
    tick_label_formatter = FormatStrFormatter(label_X.replace('X',data_fmt))
    cbar_tick_label_formatter = FormatStrFormatter(label_X.replace('X',cbar_data_fmt))


    pl_left, pl_bottom, pl_right,pl_top = plot_adjust

    if filled==0:
        pl_right*=1.2
        pl_bottom*=2
    if hasattr(self,'rms'):
        pass
    elif not nvals:
        raise(Exception,'Missing nvals/RMS input for level calculation')
    else: # if nvals was given and self.rms does not exist
        self.calc_rms(nvals, rms_area)
    ### for the plotting box
    ylen, xlen = self.d[0].shape
    ycoords = arange(-ylen/2,ylen/2,1)*self.dec_cdelt
    xcoords = arange(-xlen/2,xlen/2,1)*self.ra_cdelt
    left,right,bottom,top = self.extent
    # set plot boundaries
    if box == [0,0]:
        i1,i2 = left,right
        j1,j2 = bottom,top
    elif box != [0,0]:
        i1,i2 = array([-1,1])*box[0]/2.*sign(self.ra_cdelt)
        j1,j2 = array([-1,1])*box[1]/2.*sign(self.dec_cdelt)
    # calculate the moments (0/1 so far)
    Mom = Moments(self, chvals = chvals, nsig = nsig)
    #
    d1,d2,r1,r2 = self.parse_region([0, 0, box[0]])
    print 'No. sigmas for max in box: ', (Mom.zero[d1:d2,r1:r2].max() / Mom.sigma)

    x1,x2,y1,y2 = self.parse_region(region)

    Mom.region_pxl = (x1,x2,y1,y2)

    # calculate the velocity intervals
    print('remember that the calculation for vmin and vmax uses the region keyword')
    from string import lower
    # MOMENT 1 stats
    arr = array([line for line in Mom.one[y1:y2,x1:x2].flatten() if str(line).lower() != 'nan'])
    Mom.one_sigma = arr.std()
    from scipy import median
    Mom.one_median = median(arr)
    # MOMENT 2 stats
    arr = array([line for line in Mom.two[y1:y2,x1:x2].flatten() if str(line).lower() != 'nan'])
    Mom.two_sigma = arr.std()
    Mom.two_median = median(arr)


    # now create the levels
    if negcontours:
        # the levels_neg/pos arrays, are [1 sigma, 2 sigma, 3 sigma, 4 sigma] and so on
        # therefore, when slicing, we need to first take -1, i.e. slice index 0 is 1 sigma
        # so 3 sigma is slice index 2 -> nsig - 1
        levs = concatenate((
                Mom.levels_neg[flipud(arange(nsig - 1, alen(Mom.levels_neg), nsjump))],
                Mom.levels_pos[arange(nsig - 1, alen(Mom.levels_pos), nsjump)]
                          ))
        levs_contour = concatenate((
                Mom.levels_neg[flipud(arange(nsig - 1, alen(Mom.levels_neg), 2 * nsjump))],
                Mom.levels_pos[arange(nsig - 1, alen(Mom.levels_pos), 2 * nsjump)]
                                  ))
    else:
        #~ levs = arange(nsig*img_sigma,img_max+img_sigma,nsjump*img_sigma)
        levs = Mom.levels_pos[arange(nsig -1, alen(Mom.levels_pos), nsjump)]
        #~ levs_contour = arange(nsig*img_sigma,img_max+img_sigma,2*nsjump*img_sigma)
        levs_contour = Mom.levels_pos[arange(nsig - 1, alen(Mom.levels_pos), 2 * nsjump)]
    #levs = arange(nsig*img_sigma, img_max, nsjump*img_sigma)
    # print some info out
    def print_title(title):
        print '\n{0}'.format('='*70)
        print '{0}{1}'.format(' '*30, title)
        print '{0}'.format('='*70)

    print_title('LINE DATA')
    mjybeamkms = u'mJy\u00b7beam\u207b\u00b9\u00b7km\u00b7s\u207b\u00b9'
    mjybeamchannel = u'mJy\u00b7beam\u207b\u00b9\u00b7channel\u207b\u00b9'
    print '\n Channels \t: {0:3d} to {1:3d}'.format(Mom.channels.min(), Mom.channels.max())
    print u' RMS \t\t: {0:>10.2f}  {1}'.format(1e3*self.rms, mjybeamchannel)
    print u' Sigma \t\t: {0:>10.2f}  {1}'.format(1e3*Mom.sigma, mjybeamkms)
    print u' Map min \t: {0:>10.2f}  {1}'.format(1e3*Mom.minimum, mjybeamkms)
    print u' Map max \t: {0:>10.2f}  {1}'.format(1e3*Mom.maximum, mjybeamkms)
    print u' Start sigma \t: {0:>10.2f}  {1}'.format(1e3*nsig*Mom.sigma, mjybeamkms)
    print u' Sigma step \t: {0:>10.2f}  {1}'.format(1e3*nsjump*Mom.sigma, mjybeamkms)
    #
    # print beam parameters for the data
    print_title('BEAM')
    print 'Line cube:'
    print ' Minor (FWHM)\t: {0:>10.2f}   arcseconds'.format(self.bmin)
    print ' Major (FWHM)\t: {0:>10.2f}   arcseconds'.format(self.bmaj)
    print u' PA \t\t: {0:>10.2f}   degrees (0\u00b0 < \u0398 < 180\u00b0)'.format(self.bpa)
    print u' Gain\t\t: {0:>10.2f}   mJy\u00b7K\u207b\u00b9\n'.format(self.gain*1e3)

    self.Mom = Mom
    if send==True:
        print 'sending you moment class and extents'
        return Mom, array(self.extent).reshape((2,2))

    ####################################################################

    ####                  END CALCULATION BLOCK                     ####

    ####################################################################

    #
    # tick density
    majorLocator = MultipleLocator(locators[0])
    minorLocator = MultipleLocator(locators[1])
    set_rc()
    #
    pl.ion()
    pl.close()
    fig = pl.figure(1, figsize=fsize)
    fig.clf()
    #
    ax1 = fig.add_subplot(221)
    #
    #
    # plot the continuum data
    #~ if cfile != None:
        #~ cont = ax.imshow(contdata.d, cmap=cm.gray_r, extent=(left,right,bottom,top))
        #~ if cbar:
            #~ cb = pl.colorbar(cont)
            #~ cb.ax.set_ylabel(str(contdata.unit))
    #~
        #~ print '-'*40
        #~ print 'Continuum data:'
        #~ print u' Minor (FWHM)\t: %2.3f \tasec' % contdata.bmin
        #~ print u' Major (FWHM)\t: %2.3f  \tasec' % contdata.bmaj
        #~ print ' PA \t\t: %2.3f Degrees (0<theta<180) \n' % contdata.bpa
    #
    #
    if colormap and filled==True:
        colormap = cm.bone_r
        # parse colormap command!
    if colormap == None or filled==False:
        # just contours
        levs = levs.round(3)
        #~ levs_contour = levs_contour.round(3)
        cs = ax1.contour(Mom.zero, levels=levs, colors=ccol, extent=self.extent)
    elif colormap!=None:
        levs = levs.round(3)
        levs_contour = levs_contour.round(3)
        cs1 = ax1.contourf(Mom.zero, levels=levs, cmap=colormap, extent=self.extent)
        #cs2 = ax1.contour(img, cs1.levels[::2], colors=ccol, extent=self.extent)
        #return cs1
        cs2 = ax1.contour(Mom.zero, levels=levs_contour, colors=ccol, extent=self.extent)
        cax = fig.add_axes([0.35, 0.54, 0.02, 0.4]) # left, bottom, width, height
        cbar = pl.colorbar(cs1, cax=cax, ticks=levs_contour, format=cbar_tick_label_formatter)
        cbar.ax.set_ylabel(self.unitint)
        #~ cbar = pl.colorbar(cs1, ticks=levs_contour, format=cbar_tick_label_formatter)#label_X.replace('X','%2.2f'))
        cbar.add_lines(cs2)
        if str(self.unit) == 'Jy/beam':
            cbar.ax.set_ylabel(r'Jy\,beam$^{-1}$')
        else:
            cbar.ax.set_ylabel(str(self.unit))
    else:
        line = ax1.contour(Mom.zero, levels=levs, colors='r', extent=self.extent)

    #ax1.text(0.5,0.5,'test',transform = ax1.transax1es)
    draw_beam(ax1, self)
    #~ draw_fov(ax1, self)

    # check distance key
    if sbar.has_key('dist'):
        dist_mark = sbar['dist']
    elif self.dist != 0:
        dist_mark = self.dist
    else:
        dist_mark = 200
    # check the length of the scale bar
    if sbar.has_key('au'):
        au_mark = sbar['au']
    else:
        au_mark = 200
    print 'Using distance {0} pc to source. Scale bar length {1} AU'.format(dist_mark, au_mark)
    draw_sizebar(ax1, self, dist=dist_mark, au=au_mark)
    # parse the cpeak keyword
    #
    ax1.xaxis.set_major_formatter(tick_label_formatter)
    ax1.yaxis.set_major_formatter(tick_label_formatter)

    ax1.xaxis.set_major_locator(majorLocator)
    ax1.xaxis.set_minor_locator(minorLocator)
    ax1.yaxis.set_major_locator(majorLocator)
    ax1.yaxis.set_minor_locator(minorLocator)

    ax1.set_xlabel('RA Offset ($^{\prime\prime}$)')
    ax1.set_ylabel('Dec Offset ($^{\prime\prime}$)')
    if source.has_key('title'):
        ax1.text(0.05,0.92, source['title'], transform = ax1.transAxes)
    else:
        ax1.text(0.05,0.92, self.obj, transform = ax1.transAxes)
    ax1.set_xlim(i1,i2)
    ax1.set_ylim(j1,j2)
    ax1.set_aspect(1)


    ax2 = fig.add_subplot(222, sharex=ax1, sharey=ax1)

    m = Mom.one_median
    w = Mom.one_sigma
    if type:
        im = ax2.imshow(Mom.one,cmap=cm.jet,vmin=m-2*w,vmax=m+2*w,extent=(left,right,bottom,top),interpolation='nearest')
    elif not type:
        im = ax2.contourf(Mom.one, levels=arange(m-2*w,m+2*w,4*w/20), cmap=cm.jet, extent=(left,right,bottom,top))
    cs2 = ax2.contour(Mom.zero, levels=levs_contour, colors='k', extent=self.extent)
    #~ cbar = pl.colorbar(im,format=cbar_tick_label_formatter) #label_X.replace('X','%2.2f'))
    #~ cbar.ax.set_ylabel(r'km\,s$^{\sf -1}$')

    cax = fig.add_axes([0.837, 0.54, 0.02, 0.4]) # left, bottom, width, height
    cbar = pl.colorbar(im, format='%.2f', cax=cax)
    cbar.ax.set_ylabel(r'km\,s$^{\sf -1}$')


    ax3 = fig.add_subplot(223, sharex=ax1, sharey=ax1)

    m = Mom.two_median
    w = Mom.two_sigma
    if type:
        im = ax3.imshow(Mom.two, cmap=cm.jet,vmin=m-2*w,vmax=m+2*w,extent=(left,right,bottom,top),interpolation='nearest')
    elif not type:
        im = ax3.contourf(Mom.two, levels=arange(m-2*w,m+2*w,4*w/20), cmap=cm.jet, extent=(left,right,bottom,top))
    cs3 = ax3.contour(Mom.zero, levels=levs_contour, colors='k', extent=self.extent)
    #~ cbar = pl.colorbar(im,format=cbar_tick_label_formatter) #label_X.replace('X','%2.2f'))
    #~ cbar.ax.set_ylabel(r'km\,s$^{\sf -1}$')

    cax = fig.add_axes([0.35, 0.06, 0.02, 0.4]) # left, bottom, width, height
    cbar = pl.colorbar(im, format='%.2f', cax=cax)
    cbar.ax.set_ylabel(r'km\,s$^{\sf -1}$')




    #draw_beam(ax, self,box=0)
    #~ draw_fov(ax2, self)

    # check distance key
    if sbar.has_key('dist'):
        dist_mark = sbar['dist']
    elif self.dist != 0:
        dist_mark = self.dist
    else:
        dist_mark = 200
    # check the length of the scale bar
    if sbar.has_key('au'):
        au_mark = sbar['au']
    else:
        au_mark = 200
    print 'Using distance {0} pc to source. Scale bar length {1} AU'.format(dist_mark, au_mark)
    draw_sizebar(ax1, self, dist=dist_mark, au=au_mark)


    if len(cpeak) == 3:
        mark = cpeak[2]
        xmark, ymark = cpeak[0:2]
    elif len(cpeak) >4:
        xmark = cpeak[0:-1:2]
        ymark = cpeak[1:-1:2]
        mark = cpeak[-1]
    else:
        mark = '+k'
    cross1 = ax1.plot(xmark, ymark, mark, ms=6)#, mew=3, alpha=0.9)
    cross2 = ax2.plot(xmark, ymark, mark, ms=6)#  mew=3, alpha=0.9)

    #
    #~ ax2.xaxis.set_major_formatter(tick_label_formatter)
    #~ ax2.yaxis.set_major_formatter(tick_label_formatter)

    #~ ax2.xaxis.set_major_locator(majorLocator)
    #~ ax2.xaxis.set_minor_locator(minorLocator)
    #~ ax2.yaxis.set_major_locator(majorLocator)
    #~ ax2.yaxis.set_minor_locator(minorLocator)

    ax2.set_xlabel('RA Offset ($^{\prime\prime}$)')
    ax2.set_ylabel('Dec Offset ($^{\prime\prime}$)')
    if source.has_key('title'):
        ax2.text(0.05,0.92, source['title'], transform = ax2.transAxes)
    else:
        ax2.text(0.05,0.92, self.obj, transform = ax2.transAxes)
    ax2.set_xlim(i1,i2)
    ax2.set_ylim(j1,j2)
    ax2.set_aspect(1)




    pl.subplots_adjust(left=pl_left, bottom=pl_bottom, right=pl_right, top=pl_top,wspace=0.8)


# temp function to plot continuum
# to adavis.py
def plot_map(self,
                chvals=None,
                nvals=None,
                region=[0,0,0,0],
                filled=True,
                rms_area = [0, 0, 10],
                type=0,
                box=[0,0],
                nx=6,
                nsig=3,
                nsjump=2,
                source = dict(v_sys=0,dist=250),
                font={'family':'serif', 'serif': ['Times New Roman'],
                'size':8},
                fit = dict(params = None, interactive = False),
                send=False,
                quality=[300, 300],
                cbar=True,
                colormap=True,
                plot_adjust= [0.07, 0.06, 0.82, 0.99],
                cpeak=[0,0,'k'],
                ccol='k',
                sbar=dict(au=200),
                locators = [2,1],
                telescope=None,
                negcontours=True,
                fsize=(_AandA.ONE_COL_FIG_WIDTH, _AandA.ONE_COL_FIG_WIDTH),
                **kwargs):
    #imports
    #import scipy as sp
    from scipy import array, where, median, std, sign, arange, alen, vstack, \
    concatenate, sqrt, log10, ones, sqrt, flipud
    #import pyfits as pf
    from time import sleep
    import matplotlib.pyplot as pl
    #from matplotlib.patches import Circle
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from matplotlib import cm, rc
    #
    ################################
    # setting the tick label font  #
    ################################
    # we are using latex formatting here!
    # following sets the ticklabel format string
    data_fmt = '%g'                 # for normal y and x axis
    cbar_data_fmt = '%2.2f'         # for normal y and x axis
    label_X = parse_tick_font(font)
    tick_label_formatter = FormatStrFormatter(label_X.replace('X',data_fmt))
    cbar_tick_label_formatter = FormatStrFormatter(label_X.replace('X',cbar_data_fmt))


    pl_left, pl_bottom, pl_right,pl_top = plot_adjust

    if filled==0:
        pl_right*=1.2
        pl_bottom*=2
    ### for the plotting box
    ylen, xlen = self.d.shape
    ycoords = arange(-ylen/2,ylen/2,1)*self.dec_cdelt
    xcoords = arange(-xlen/2,xlen/2,1)*self.ra_cdelt
    left,right,bottom,top = self.extent
    # set plot boundaries
    if box == [0,0]:
        i1,i2 = left,right
        j1,j2 = bottom,top
    elif box != [0,0]:
        i1,i2 = array([-1,1])*box[0]/2.*sign(self.ra_cdelt)
        j1,j2 = array([-1,1])*box[1]/2.*sign(self.dec_cdelt)
    x1,x2,y1,y2 = self.parse_region(region)

    # print beam parameters for the data
    print '='*40
    print ' '*15,'BEAM(S)'
    print '='*40
    print 'Line cube:'
    print u' Minor (FWHM)\t: %2.3f \tasec' % self.bmin
    print u' Major (FWHM)\t: %2.3f  \tasec' % self.bmaj
    print ' PA \t\t: %2.3f \tDegrees (0<theta<180)' % self.bpa
    #~ print u' Gain\t\t: %2.4f \tJy\u00b7K\u207b\u00b9\n' % self.gain

    # tick density
    majorLocator = MultipleLocator(locators[0])
    minorLocator = MultipleLocator(locators[1])
    set_rc()



    pl.ion()
    pl.close()
    fig = pl.figure(1, figsize=fsize)
    fig.clf()
    #
    ax1 = fig.add_subplot(111)
    #

    #
    # plot the continuum data
    #~ if cfile != None:
        #~ cont = ax.imshow(contdata.d, cmap=cm.gray_r, extent=(left,right,bottom,top))
        #~ if cbar:
            #~ cb = pl.colorbar(cont)
            #~ cb.ax.set_ylabel(str(contdata.unit))
#~
        #~ print '-'*40
        #~ print 'Continuum data:'
        #~ print u' Minor (FWHM)\t: %2.3f \tasec' % contdata.bmin
        #~ print u' Major (FWHM)\t: %2.3f  \tasec' % contdata.bmaj
        #~ print ' PA \t\t: %2.3f Degrees (0<theta<180) \n' % contdata.bpa
    #
    #
    if colormap and filled==True:
        colormap = cm.bone_r
        # parse colormap command!
    if colormap == None or filled==False:
        # just contours
        #~ levs = levs.round(3)
        #~ levs_contour = levs_contour.round(3)
        cs = ax1.contour(self.d, colors=ccol, extent=self.extent)
    elif colormap!=None:
        #~ levs = levs.round(3)
        #~ levs_contour = levs_contour.round(3)
        cs1 = ax1.imshow(self.d, extent=self.extent)
        #cs2 = ax1.contour(img, cs1.levels[::2], colors=ccol, extent=self.extent)
        #return cs1
        #~ cs2 = ax1.contour(Mom.zero, levels=levs_contour, colors=ccol, extent=self.extent)
        #~ cax = fig.add_axes([0.35, 0.185, 0.02, 0.68]) # left, bottom, width, height
        #~ cbar = pl.colorbar(cs1, cax=cax, ticks=levs_contour, format=cbar_tick_label_formatter)
        #~ cbar.ax.set_ylabel(self.unitint)
        cbar = pl.colorbar(cs1) #, ticks=levs_contour, format=cbar_tick_label_formatter)#label_X.replace('X','%2.2f'))
        #~ cbar.add_lines(cs2)
        if str(self.unit).lower() == 'jy/beam':
            cbar.ax.set_ylabel(r'Jy\,beam$^{-1}$')
        else:
            cbar.ax.set_ylabel(str(self.unit))
    else:
        line = ax1.contour(self.d, colors='r', extent=self.extent)

    #ax1.text(0.5,0.5,'test',transform = ax1.transax1es)
    draw_beam(ax1, self)
    #~ draw_fov(ax1, self)

    # check distance key
    if sbar.has_key('dist'):
        dist_mark = sbar['dist']
    elif self.dist != 0:
        dist_mark = self.dist
    else:
        dist_mark = 200
    # check the length of the scale bar
    if sbar.has_key('au'):
        au_mark = sbar['au']
    else:
        au_mark = 200
    print 'Using distance {0} pc to source. Scale bar length {1} AU'.format(dist_mark, au_mark)
    draw_sizebar(ax1, self, dist=dist_mark, au=au_mark)
    # parse the cpeak keyword
    if len(cpeak) == 3:
        mark = cpeak[2]
        xmark, ymark = cpeak[0:2]
    elif len(cpeak) >4:
        xmark = cpeak[0:-1:2]
        ymark = cpeak[1:-1:2]
        mark = cpeak[-1]
    else:
        mark = '+k'
    cross = ax1.plot(xmark, ymark, mark, ms=6)#, mew=3, alpha=0.9)

    #
    ax1.xaxis.set_major_formatter(tick_label_formatter)
    ax1.yaxis.set_major_formatter(tick_label_formatter)

    ax1.xaxis.set_major_locator(majorLocator)
    ax1.xaxis.set_minor_locator(minorLocator)
    ax1.yaxis.set_major_locator(majorLocator)
    ax1.yaxis.set_minor_locator(minorLocator)

    ax1.set_xlabel('RA Offset ($^{\prime\prime}$)')
    ax1.set_ylabel('Dec Offset ($^{\prime\prime}$)')
    if source.has_key('title'):
        ax1.text(0.05,0.92, source['title'], transform = ax1.transAxes)
    else:
        ax1.text(0.05,0.92, self.obj, transform = ax1.transAxes)
    ax1.set_xlim(i1,i2)
    ax1.set_ylim(j1,j2)
    ax1.set_aspect(1)

########################################################################
# OLD (NEEDS CORRECTIONS TO WORK)
# to adavis.py

def plot_moment2 (self,
                cfile=None,
                chvals=None,
                nvals=None,
                region=[0,0,0,0],
                nx=6,
                filled=True,
                box=[0,0],
                nsig=3,
                nsjump=2,
                source = dict(v_sys=0),
                font={'family':'sans-serif', 'sans-serif': ['Arial', 'Helvetica'],
                'size':22, 'weight':'bold'},
                fit = dict(gauss = None, params = None, continuum = False,
                interactive = False),
                send=False,
                quality=[150, 72],
                cbar=True,
                plot_adjust= [0.12, 0.01, 0.74, 0.99],
                cpeak=[0,0,'k'],
                ccol='r',
                sbar=dict(dist=250,au=200),
                locators = [2,1],
                telescope=None,
                type=1):
    """ Function doc """

    # imports
    #import scipy as sp
    from scipy import array, where, median, std, sign, arange, alen, vstack, \
    concatenate, sqrt, log10, ones, nan
    #import pyfits as pf
    import matplotlib.pyplot as pl
    #from matplotlib.patches import Circle
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from matplotlib import cm
    #

    set_rc(font=font, quality=quality)


    ################################
    # setting the tick label font  #
    ################################
    # following sets the ticklabel format string
    data_fmt = '%g'         # for normal y and x axis
    cbar_data_fmt = '%2.2f'         # for normal y and x axis
    label_X = parse_tick_font(font)
    tick_label_formatter = FormatStrFormatter(label_X.replace('X',data_fmt))
    cbar_tick_label_formatter = FormatStrFormatter(label_X.replace('X',cbar_data_fmt))

    pl_left, pl_bottom, pl_right,pl_top = plot_adjust

    velocity = self.v_arr - source['vsys'] #now all the velocities are based on the LSR

    # parse the region, get the area for the spectra
    x1,x2,y1,y2 = parse_region(self, region)
    # if any is limits are None, draw spectra and ask for them/it
    if nvals==None or chvals==None:
        # make it a bit interactive
        # draw spectrum
        plot_spectrum(self,region=region, source=source)
        # ask for chvals and nvals
        # if one is not None, it will just be sent around
        chvals, nvals = get_vals(chvals, nvals)

    # calculate common stuff
    ### for the plotting box
    ylen, xlen = self.d[0].shape
    ycoords = arange(-ylen/2,ylen/2,1)*self.dec_cdelt
    xcoords = arange(-xlen/2,xlen/2,1)*self.ra_cdelt
    #
    # for the extent keyword
    #~ left, right = xcoords[0],xcoords[-1]
    #~ bottom, top = ycoords[0],ycoords[-1]
    left,right,bottom,top = self.extent
    #
    # parse the noise channel values
    if nvals!=None:
        # if noise has been input
        noise_channels = get_indices(velocity, nvals)
    #
    else:
        # if user didnt input noise velocity limits
        # choose channels away from the choosen (red blue) lines
        low = where((velocity>(velocity.min()+10))*(velocity<(chvals[0]-10)))[0]
        high = where((velocity>(chvals[1]+10))*(velocity<(velocity.max()-10)))[0]
        noise_channels = concatenate((low, high))
    #
    # the region to calculate the rms in
    i1,i2 = xlen/2+array([-1,1])*xlen/4
    j1,j2 = ylen/2+array([-1,1])*ylen/4

    # the noise, rms
    noise = self.d[noise_channels]
    rms = sqrt(((noise[:,j1:j2,i1:i2])**2).mean())

    # set boundaries
    if box == [0,0]:
        i1,i2 = left,right
        j1,j2 = bottom,top
    elif box != [0,0]:
        i1,i2 = array([-1,1])*box[0]/2.*sign(self.ra_cdelt)
        j1,j2 = array([-1,1])*box[1]/2.*sign(self.dec_cdelt)

    img_channels = get_indices(velocity,chvals)
    imgs = self.d[img_channels]
    #
    # do the summing
    #
    # M1 = dv * sum(I(a,b,vi))_x1^x2
    # or perhaps without dv to save calc?
    moment0 = imgs.sum(axis=0)*abs(self.v_cdeltkms)
    moment0_sigma = sqrt(alen(imgs))*rms*abs(self.v_cdeltkms)

    """
    Fit gaussian in the chvals interval, using velocity and cut out spectra
    in all pixels in region/box.
    Plot it!



    # levels
    moment0_max = moment0.max()
    moment0_min = moment0.min()
    levels = concatenate((-1*flipud(arange(nsig*moment0_sigma,abs(moment0_min+moment0_sigma),nsjump*moment0_sigma)), arange(nsig*moment0_sigma,moment0_max+moment0_sigma,nsjump*moment0_sigma)))

    # calculate moment 1 map

    velocities = linedata.v_arr[img_channels]
    # calculate the denominator
    Isum = imgs.sum(axis=0)

    minpar =
    limmin =
    maxpar =
    limmax =

    HAVE TO ESTIMATE LINE CENTER FOR ALL PIXELS
fit = dict(type=None, params=[(0.09, 7.3, 3.6)],
                guess=False, interactive=False, fixlist=None, error=None,
                limmin=None, minpar=None, limmax=None, maxpar=None, tie=None)

            if not fit.has_key('error'):
            errmsg='You have to supply an error to the fit, they\'re kind of important you know.'
            print stylify(errmsg)
            return ; sysexit()

        fwhmfromsig = 2*sqrt(2*log(2)) # the constant
        fwhm = lambda x: fwhmfromsig*x
        sigma = lambda x: x/fwhmfromsig

        # fits up to five (5) 1D gaussian(s) to the spectral data
        if chvals!=None: # if we supplied limits for the fit
            ch = where((velocity>v1)*(velocity<v2))[0]
            Fx = spect[ch]
            X = velocity[ch]
        else: # else use the eveything
            Fx = spect
            X = velocity
        #
        #
        p = fit['params'] # parameters for the fit

        #
        if not fit.has_key('limmin'):
            fit['limmin'] = None
            fit['minpar'] = None
        if not fit.has_key('limmax'):
            fit['limmax'] = None
            fit['maxpar'] = None
        if not fit.has_key('fixlist'):
            fit['fixlist'] = None
        if not fit.has_key('tie'):
            fit['tie'] = None
        #
        from time import time
        t1= time()
        #
        params, errors, chi2, mp = fit_gauss1d((X,Fx), params=p, fixlist=fit['fixlist'], \
                minbool=fit['limmin'], minpar=fit['minpar'], \
                maxbool=fit['limmax'], maxpar=fit['maxpar'], \
                err=fit['error'], tie=fit['tie'], verbose=0, full_output=1)
        #
        print ' Done in %2.3f seconds \n' % (time()-t1)

        #
        print ' Number of fits : ', alen(params)/3
        print ' Fit status : ', mp.status, '(if 0, it should have halted)\n'


            'value' - the starting parameter value (but see the START_PARAMS
                             parameter for more information).

            'fixed' - a boolean value, whether the parameter is to be held
                             fixed or not.  Fixed parameters are not varied by
                             MPFIT, but are passed on to MYFUNCT for evaluation.

            'limited' - a two-element boolean array.  If the first/second
                               element is set, then the parameter is bounded on the
                               lower/upper side.  A parameter can be bounded on both
                               sides.  Both LIMITED and LIMITS must be given
                               together.

            'limits' - a two-element float array.  Gives the
                              parameter limits on the lower and upper sides,
                              respectively.  Zero, one or two of these values can be
                              set, depending on the values of LIMITED.  Both LIMITED
                              and LIMITS must be given together.

            'parname' - a string, giving the name of the parameter.  The
                               fitting code of MPFIT does not use this tag in any
                               way.  However, the default iterfunct will print the
                               parameter name if available.

            'step' - the step size to be used in calculating the numerical
                            derivatives.  If set to zero, then the step size is
                            computed automatically.  Ignored when AUTODERIVATIVE=0.

            'mpside' - the sidedness of the finite difference when computing
                              numerical derivatives.  This field can take four
                              values:

                                     0 - one-sided derivative computed automatically
                                     1 - one-sided derivative (f(x+h) - f(x)  )/h
                                    -1 - one-sided derivative (f(x)   - f(x-h))/h
                                     2 - two-sided derivative (f(x+h) - f(x-h))/(2*h)

                             Where H is the STEP parameter described above.  The
                             "automatic" one-sided derivative method will chose a
                             direction for the finite difference which does not
                             violate any constraints.  The other methods do not
                             perform this check.  The two-sided method is in
                             principle more precise, but requires twice as many
                             function evaluations.  Default: 0.

            'mpmaxstep' - the maximum change to be made in the parameter
                                     value.  During the fitting process, the parameter
                                     will never be changed by more than this value in
                                     one iteration.

                                     A value of 0 indicates no maximum.  Default: 0.

            'tied' - a string expression which "ties" the parameter to other
                            free or fixed parameters.  Any expression involving
                            constants and the parameter array P are permitted.
                            Example: if parameter 2 is always to be twice parameter
                            1 then use the following: parinfo(2).tied = '2 * p(1)'.
                            Since they are totally constrained, tied parameters are
                            considered to be fixed; no errors are computed for them.
                            [ NOTE: the PARNAME can't be used in expressions. ]

            'mpprint' - if set to 1, then the default iterfunct will print the
                               parameter value.  If set to 0, the parameter value
                               will not be printed.  This tag can be used to
                               selectively print only a few parameter values out of
                               many.  Default: 1 (all parameters printed)

Future modifications to the PARINFO structure, if any, will involve
     adding dictionary tags beginning with the two letters "MP".
     Therefore programmers are urged to avoid using tags starting with
     the same letters; otherwise they are free to include their own
     fields within the PARINFO structure, and they will be ignored.

     PARINFO Example:
     parinfo = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
                                                                                                    for i in range(5)]
     parinfo[0]['fixed'] = 1
     parinfo[4]['limited'][0] = 1
     parinfo[4]['limits'][0]  = 50.
     values = [5.7, 2.2, 500., 1.5, 2000.]
     for i in range(5): parinfo[i]['value']=values[i]

     A total of 5 parameters, with starting values of 5.7,
     2.2, 500, 1.5, and 2000 are given.  The first parameter
     is fixed at a value of 5.7, and the last parameter is
     constrained to be above 50.




    """
    # create array that matches imgs array, so they can be multiplied
    velocities_matrix = array([ones(imgs.shape[1:])*i for i in velocities])
    # set values which are lower than 3 sigma to 'nan' i.e. mask the array
    # use the moment0 to determine where, but the Isum to set it
    Isum[where(moment0<(nsig*moment0_sigma))] = nan
    # now, calculate the numerator of the division
    Ivsum = (imgs*velocities_matrix).sum(axis=0)
    # division
    moment1 = Ivsum/Isum


    ###
    ###
    ###
    # tick density
    majorLocator = MultipleLocator(locators[0])
    minorLocator = MultipleLocator(locators[1])


    if send==True:
        return {'moment0': moment0,
                'moment1': moment1,
                'levels':levels,
                'moment0_sigma':moment0_sigma,
                'extent':self.extent,
                'lims': (i1,i2,j1,j2),
                'data':linedata}
    #
    #
    # calculating the velocity for vmin and vmax
    #
    #
    # remember it uses the region keyword, so
    print('remember that the calculation for vmin and vmax uses the region keyword')
    from string import lower
    arr = array([line for line in moment1[y1:y2,x1:x2].flatten() if lower(str(line)) != 'nan'])
    w = arr.std()
    from scipy import median
    m = median(arr)

    """

    here you should add a histogram calculation, and fit a gaussian or something
    and take some representative width for vmin to vmax
    """
    pl.ion()
    pl.close()
    if cbar:
        fig = pl.figure(1, (9., 7.))
    elif not cbar:
        fig = pl.figure(1, (9.,7.))
    fig.clf()
    #
    ax = fig.add_subplot(111)
    #
    # print beam parameters for the data
    gain = 8.168e-25*(linedata.restfreq)**2*linedata.bmin*linedata.bmaj
    print '='*40
    print ' '*15,'BEAM(S)'
    print '='*40
    print 'Line cube:'
    print u' Minor (FWHM)\t: %2.3f \tasec' % linedata.bmin
    print u' Major (FWHM)\t: %2.3f  \tasec' % linedata.bmaj
    print ' PA \t\t: %2.3f \tDegrees (0<theta<180)' % linedata.bpa
    print u' Gain\t\t: %2.4f \tJy\u00b7K\u207b\u00b9\n' % gain

    #levs=levs.round(2)
    #cs1 = ax.contourf(img, levels=levs, cmap=cm.bone_r, extent=(left,right,bottom,top))
    #cs2 = ax.contour(cs1, levels=cs1.levels[::2], colors=ccol, extent=(left,right,bottom,top))
    #im = pl.imshow(moment1,vmin=6.79,vmax=7,extent=(left,right,bottom,top))
    if type:
        im = pl.imshow(moment1,cmap=cm.jet,vmin=m-2*w,vmax=m+2*w,extent=self.extent,interpolation='nearest')
    elif not type:
        im = ax.contourf(moment1, levels=arange(m-2*w,m+2*w,4*w/20), cmap=cm.jet, extent=self.extent)
    cs2 = ax.contour(moment0, levels=levels, colors='k', extent=self.extent)
    cbar = pl.colorbar(im,format=cbar_tick_label_formatter) #label_X.replace('X','%2.2f'))

    if str(linedata.unit) == 'Jy/beam':
        cbar.ax.set_ylabel(r'Jy\,beam$^{\sf-1}$')
    else:
        cbar.ax.set_ylabel(str(linedata.unit))

    ax.set_xlim(i1,i2)
    ax.set_ylim(j1,j2)

    #draw_beam(ax, linedata,box=0)
    draw_fov(ax, linedata)

    if 'dist' and 'au' in sbar.keys():
        draw_sizebar(ax,linedata,dist=sbar['dist'],au=sbar['au'])
    elif 'dist' in sbar.keys():
        draw_sizebar(ax,linedata,dist=sbar['dist'])
    else :
        raise ParError(sbar)
    if len(cpeak) == 3:
        mark = cpeak[2]
        xmark, ymark = cpeak[0:2]
    elif len(cpeak) == 5:
        xmark = cpeak[0:-1:2]
        ymark = cpeak[1:-1:2]
        mark = cpeak[-1]
    else:
        mark = '+k'
    cross = ax.plot(xmark, ymark, mark, ms=13, mew=2, alpha=0.8)

    #
    ax.xaxis.set_major_formatter(tick_label_formatter)
    ax.yaxis.set_major_formatter(tick_label_formatter)

    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    #pl.setp(ax.get_xticklabels(),family='sans-serif')
    #pl.setp(ax.get_yticklabels(),family='sans-serif')

    ax.set_xlabel('RA offset ($^{\prime\prime}$)')
    ax.set_ylabel('DEC offset ($^{\prime\prime}$)')

    #fig.suptitle(linedata.obj+' (%s km s$^{-1}$)'% str(linedata.unit))
    #if cfile!=None:
    #    fig.subplots_adjust(bottom=0.08, right=0.77, top=0.90)
    ax.text(0.05,0.93,
    linedata.obj,
    transform = ax.transAxes,
    backgroundcolor='w')
    ax.set_aspect(1)
    #
    #
    pl.subplots_adjust(left=pl_left, bottom=pl_bottom, right=pl_right, top=pl_top)
#
def plot_vgradient(radecfile,
                fitsfile=None,
                rms=None,
                chvals=None,
                nvals=None,
                xy_cont = None,
                region=[0,0,0,0],
                source = dict(v_sys=0),
                freq=False,
                font={'family':'serif', 'serif': ['Times', 'Times New Roman'],
                'size': 20, 'weight':'bold'},
                bin=False,
                fit = dict(type=None, params=[(0.09, 7.3, 3.6)],
                guess=False, interactive=False, fixlist=None, error=None,
                limmin=None, minpar=None, limmax=None, maxpar=None, tie=None),
                send=False,
                quality=[150, 72],
                plot_adjust= [0.12, 0.09, 0.99, 0.99],
                lines = [],\
                axspace = [1.01, 1.01, 1.01, 1.05],
                ylimits=None):
    """



    TODO: add a plot of the fit perpendicular to the gradient!


    """

    # imports
    #import scipy as sp
    from scipy import array, where, median, std, sign, arange, alen, vstack, \
    concatenate, sqrt, log10, ones, nan
    #import pyfits as pf
    import matplotlib.pyplot as pl
    #from matplotlib.patches import Circle
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from matplotlib import cm, rc
    #

    ################################
    # setting global rc properties #
    ################################
    rc('savefig', **{'dpi': 300})
    rc('text', usetex=True)
    rc('savefig', **{'dpi': quality[0]})
    rc('figure',**{'facecolor': '1', 'dpi': quality[1]})

    #to set the global font properties
    rc('font', **font)

    ################################
    # setting the tick label font  #
    ################################
    # we are using latex formatting here!
    # following sets the ticklabel format string
#    data_fmt = '%g'         # for normal y and x axis
#    cbar_data_fmt = '%2.2f'         # for normal y and x axis
#    label_X = parse_tick_font(font)
#    tick_label_formatter = FormatStrFormatter(label_X.replace('X',data_fmt))
#    cbar_tick_label_formatter = FormatStrFormatter(label_X.replace('X',cbar_data_fmt))
    # ticksize
    rc('xtick',**{'minor.size':3, 'major.size':7})
    rc('ytick',**{'minor.size':3, 'major.size':7})

    # linewidths
    rc('axes', linewidth=1)
    rc('lines', linewidth=1.5, markeredgewidth=1)

    pl_left, pl_bottom, pl_right,pl_top = plot_adjust


    #
    # TODO: be able to add RMS directly, instead of loading fits-file all the time
    #
    if fitsfile!=None and nvals!=None:
        # first get the fitsfile
        linedata = Fits(filename,telescope=telescope)
        #
        linedata.v_arr = linedata.v_arr - source['vsys'] #now all the velocities are based on the LSR
         # parse the region, get the area for the spectra
        x1,x2,y1,y2 = parse_region(linedata, region)
    elif fitsfile==None and rms==None:
        print 'If you do not supply a fitsfile or rms value the calculations cannot proceed.'
        sysexit()
    elif nvals==None and rms==None:
        print 'Need eiter a RMS value, or a fitsfile and velocities over\
        which to calculate the RMS.'
        sysexit()
    elif rms!=None:
        pass

    if nvals!=None:
        ## get the nvals and calculate RMS etc
        print '='*40
        print ' '*11,'Noise statistics'
        # calculate the rms from the channels in the spectra
        # accounts for it even if it is binned
        #
        # image rms
        # change x1,x2,y1,y2 to quarter region
        # change so that when binning, the rms i calculated
        # x1,x2
        zlen, ylen, xlen = data.d.shape
        ydelt = ylen/6
        xdelt = xlen/6
        i1,i2 = xlen/2-xdelt, xlen/2+xdelt
        j1,j2 = ylen/2-ydelt, ylen/2+ydelt
        rms = sqrt(((data.d[get_indices(velocity,nvals),j1:j2,i1:i2])**2).mean())
        rms_mjy = rms*1e3

        xlen = data.d.shape
        #rms_0 = sqrt(((spect[get_indices(velocity,nvals)])**2).mean())
        #rms_2 = sqrt(((data.d[get_indices(velocity,nvals),:,:])**2).mean())

        #rms_0= rms/sqrt(abs(data.v_cdeltkms))
        #print rms_0
        #print 'rms_0 =', rms_0*1e3
        #print 'rms_2 =', rms_2*1e3
        #
        # the sensitivity
        s = rms/sqrt(abs(velocity_delta))
        s_mjy = s*1e3
        # the channels used
        ind = get_indices(velocity,nvals)
        ind1, ind2 = ind.min(), ind.max()
        print u'\n RMS \t\t: %2.3f mJy\u00b7beam\u207b\u00b9\u00b7channel\u207b\u00b9' % rms_mjy
        print u' Sensitivity \t: %2.3f mJy\u00b7beam\u207b\u00b9\u00b7km\u207b\u00b9\u00b7s' % s_mjy
        print u' Channels used \t: %d, %d (%s km\u00b7s\u207b\u00b9)' % (ind1, ind2, str(nvals))
        print u' Region \t: %s arcsec offset' % (str(region))
        print u' Channel width \t: %2.3f km\u00b7s\u207b\u00b9' % abs(velocity_delta)
        #ax_kms.plot([xmin, xmax],[3*rms,3*rms],'b--', lw=2, alpha=0.7)
        #ax_kms.plot([xmin, xmax],[-3*rms,-3*rms],'b--', lw=2, alpha=0.7)


    ####################################################################
    # need to load the data into the arrays
    # need ra, d_ra, dec, d_dec, channel numbers, flux from the radecfile

    try:
        f=open(radecfile, mode='r')
    except (IOError):
        print 'File %s does not exist!' % radecfile
        sysexit()
    lines = [line for line in f]
    f.close()
    linebyline = [lines[i].split() for i in arange(len(lines))]
    ch = []
    rapos = []
    d_ra = []
    decpos = []
    d_dec = []
    ch = []
    flux = []
    d_flux = []
    vel = []

    for i in arange(len(linebyline)):
        # have to get the chvals here if they where not input?
        # i.e. if chvals==None
        # if I just use output from GILDAS POINT uv_fit
        if linebyline[i][1:5] == ['No', 'data', 'for', 'channel']:
            ch.append(linebyline[i][5])
            rapos.append(nan)
            d_ra.append(nan)
            decpos.append(nan)
            d_dec.append(nan)
            flux.append(nan)
        elif linebyline[i][2:6] == ['data', 'points', 'for', 'channel']:
            if linebyline[i+4][0] != 'POINT':
                print 'Only supports POINT fits for now'
                sysexit()
            ch.append(linebyline[i][6])
            vel.append(linebyline[i+1][7])
            rapos.append(linebyline[i+4][3])
            d_ra.append(linebyline[i+4][5].strip(')'))
            if linebyline[i+5][0] == 'STOP':
                i+=1
            decpos.append(linebyline[i+5][3])
            d_dec.append(linebyline[i+5][5].strip(')'))
            if linebyline[i+6][0] == 'STOP':
                i+=1
            flux.append(linebyline[i+6][3])
            d_flux.append(linebyline[i+6][5].strip(')'))
    #
    # done getting all the values!

    ch = array(ch).astype('i')
    rapos = array(rapos).astype('float')
    d_ra = array(d_ra).astype('float')
    decpos = array(decpos).astype('float')
    d_dec = array(d_dec).astype('float')
    flux = array(flux).astype('float')
    d_flux = array(d_flux).astype('float')
    vel = array(vel).astype('float')

    # First get all values 3*RMS above
    rms3 = 3*rms
    i_3sigma = where(flux>=rms3)[0]
    #pl.plot(flux)
    #pl.plot([0, len(flux)],[rms3,rms3])
    if chvals==None:
        # take what is given, all of it
        # parse the channel values
        print 'You did not input any velocity values (chvals=[v1, v2])\n\
        So will try to take what is given in the radecpos-file'
        i_vel = arange(0,alen(vel),1).astype('int')
    elif chvals!=None:
        # slice!
        v1, v2 = chvals
        i_vel = where((vel>=v1)*(vel<=v2))[0]
    i = [a for a in i_3sigma if a in i_vel]
    x = rapos[i]
    dx = d_ra[i]
    y = decpos[i]
    dy = d_dec[i]
    z = vel[i]
    channels = ch[i]
    f = flux[i]
    df = d_flux[i]
    # need to slice the arrays according to chvals now and set below
    # do not includ values with flux<=3sigma

    ####################################################################
    """
        find the strongest velocity gradient in data.
        input either ra, and dec lists, or just the STDOUT of uvfit
        in GILDAS MAPPING

        Fit a plane to a cloud of points
        From
        http://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf
    """
    from scipy import pi, arange, tan, array, dot, arctan, alen, meshgrid, arcsin, cos, sin

    from matplotlib import cm
    import matplotlib.pyplot as pl
    from scipy.optimize import leastsq
    pl.ion()

    # TODO: include errors.
    #       i.e. need equation which include varance-covariance matrix for x and y
    #       -> the orthogonal 3D regression
    def fit_plane_to_cloud((x,y,z), err=None):
        """
        Fit a plane Ax + By + C = z, to a cloud of points
        Input:
            (x,y,z) - the coordinates of the points to be fit
            err     - the error, in (dx, dy)

        Output:
            G       - (A, B, C) of the fit, ie the normal is (A, B, -1)
        """
        from scipy.linalg import solve
        #x = rapos
        #y = decpos
        #z = arange(alen(decpos)) # for now, need the exact channel numbers
        # can fit all channels of the cube, and then define velocity to fit over!

        if not x.size == y.size == z.size:
            print 'input not same size'
            import sys
            sys.exit(1)

        a11 = (x**2).sum()
        a12 = (x*y).sum()
        a13 = x.sum()

        a21 = (x*y).sum()
        a22 = (y**2).sum()
        a23 = y.sum()

        a31 = x.sum()
        a32 = y.sum()
        a33 = alen(y)

        A = array([[a11, a21, a31],[a12, a22, a32],[a13, a23, a33]])


        b1 = (x*z).sum()
        b2 = (y*z).sum()
        b3 = z.sum()

        B = array([[b1],[b2],[b3]])


        G = solve(A,B)
        return G
    #
    #
    G = fit_plane_to_cloud((x,y,z))
    grad = G[0:2] # the projected gradient (normal of the plane) onto the xy-plane

    # only in the positive side of the y-axis
    # the direction of the gradient is obvious from figures etc
    # this makes the y-axis the same for all angles
    # we essentailly move the vector up to positive y-values, mirror it
    if grad[1]<0:
        grad[1] = grad[1]*(-1)

    # the angle with which we incline the coordinate system
    theta = arcsin(grad[0]/((grad**2).sum())**.5)[0]

    # the new y axis is called "s" (s,t)-plane
    s_offsets = []
    t_offsets = []
    for x,y in zip(x,y):
        x_prime = x*cos(theta) - y*sin(theta)
        y_prime = x*sin(theta) + y*cos(theta)
        t_offsets.append(x_prime)
        s_offsets.append(y_prime)

    # convert to array
    t_offsets = array(t_offsets,dtype='float')
    s_offsets = array(s_offsets,dtype='float')
    #
    #ab = (G[0]**2+G[1]**2)**.5
    #v = z/ab - G[2]/ab

    ## plot the best fit model
    from scipy import polyfit
    p = polyfit(z,s_offsets,deg=1)
    s = lambda x: p[0]*x + p[1]

    fig1 = pl.figure(1,figsize=(10,5))
    #ax1 = fig1.add_subplot(111)
    ax1 = fig1.add_subplot(121)
    ax1.plot(z,s_offsets,'x')
    ax1.plot(z,s(z))
    if xy_cont!=None:
        #new position in (s,t) coordinate system of the continuum
        s_cont = xy_cont[0]*sin(theta) + xy_cont[1]*cos(theta)
        ax1.plot([z[0],z[-1]],[s_cont, s_cont])
    else:
        print '\nINFO: No continuum position (rel to phase center) given\n\
        ->skipping...\n'
    ax1.set_xlabel('Velocity [km~s$^{-1}$]')
    ax1.set_ylabel('s-offset[$\prime\prime$]')
    #ax.set_xlim(5.8,8.2)
    #ax.set_ylim(-0.205,0.305)
    #fig1.subplots_adjust(left=0.17,bottom=0.14,right=0.98,top=0.97)
    print '\n','*'*40
    print '\t\tResults for the parallell s-axis:\n'
    print u'Velocity gradient: %2.2f km\u00b7s\u207b\u00b9\u00b7arcsec\u207b\u00b9' % p[0]**-1


    ## plot perpendicular to the best fit model
    #from scipy import polyfit
    #p = polyfit(z,t_offsets,deg=1)
    #t = lambda x: p[0]*x + p[1]

    #fig2 = pl.figure(2,figsize=(7,6))
    #ax2 = fig2.add_subplot(111)
    ax2 = fig1.add_subplot(122)
    ax2.plot(z,t_offsets,'x')
    #ax2.plot(z,t(z))
    if xy_cont!=None:
        #new position in (s,t) coordinate system of the continuum
        t_cont = xy_cont[0]*cos(theta) - xy_cont[1]*sin(theta)
        ax2.plot([z[0],z[-1]],[t_cont, t_cont])
    else:
        print '\nINFO: No continuum position (rel to phase center) given\n\
        ->skipping...\n'
    ax2.set_xlabel('Velocity [km~s$^{-1}$]')
    ax2.set_ylabel('t-offset[$\prime\prime$]')
    #ax.set_xlim(5.8,8.2)
    #ax.set_ylim(-0.205,0.305)
    #fig2.subplots_adjust(left=0.17,bottom=0.14,right=0.98,top=0.97)
    fig1.subplots_adjust(left=0.12,bottom=0.14,right=0.97,top=0.95,wspace=0.31)
    print '\n','*'*40
    print '\t\tResults for the perpendicular t-axis:\n'
    print u'Velocity gradient: %2.2f km\u00b7s\u207b\u00b9\u00b7arcsec\u207b\u00b9' % p[0]**-1
#
def plot_pv (self,
                cfile=None,
                chvals=None,
                nvals=None,
                region=[0,0,0,0],
                nx=6,
                filled=True,
                box=[0,0],
                nsig=3,
                nsjump=2,
                object=dict(v_sys=0),
                font= dict(family='serif', serif=['Times', 'Palatino',
                'New Century Schoolbook', 'Bookman', 'Computer Modern Roman'],
                 size=12, weight='bold'),\
                fit = dict(gauss = None, params = None, continuum = False,
                interactive = False),
                send=False,
                quality=[150, 72],
                cbar=True):
    """ Function doc """

    # TODO : change the contour level calculation, use RMS?
    # TODO : being able to specify rotation
    #           perhaps draw a line of where the pv-diag is
    #
    #       1 translate center to the given center (region)
    #       2 rotate by given angle
    #           - new x,y pixel scale
    #       3 cut out box
    #       4 produce pv-diagram as before

    print 'importing modules...'
    from scipy import median,zeros
    from scipy.ndimage import rotate
    from matplotlib import cm, rc
    print 'done'
    #

    ################################
    # setting global rc properties #
    ################################
    rc('savefig', **{'dpi': 300})
    rc('text', usetex=True)
    rc('savefig', **{'dpi': quality[0]})
    rc('figure',**{'facecolor': '1', 'dpi': quality[1]})

    #to set the global font properties
    rc('font', **font)

    rc('axes',linewidth=2)
    rc('lines', linewidth=1)


    ################################
    # setting the tick label font  #
    ################################
    # we are using latex formatting here!
    # following sets the ticklabel format string
    data_fmt = '%g'         # for normal y and x axis
    label_X = parse_tick_font(font)
    tick_label_formatter = FormatStrFormatter(label_X.replace('X',data_fmt))
    # ticksize
    rc('xtick',**{'minor.size':5, 'major.size':7})
    rc('ytick',**{'minor.size':5, 'major.size':7})

    # linewidths
    rc('axes', linewidth=2)
    rc('lines', linewidth=1.2, markeredgewidth=2)

    #
    # first get the fitsfile
    #~ linedata = loadcube(filename)
    linedata = self

    ############## P-V diagram
    v1, v2 = chvals
    pl.ion()
    linedata.v_arr = linedata.v_arr - object['vsys'] #now all the velocities are based on the LSR
    #
    #sum along x-axis ie axis=2
    # parse the region, get the area for the spectra
    x1,x2,y1,y2 = parse_region(linedata, region)
    #x1,x2,y1,y2  = 245,265,245,265
    #axis=2 because I sum along x axis, ie naxis1
    d_arr = linedata.d[:,y1:y2,x1:x2].sum(axis=2)*linedata.v_cdeltkms
    # now we have 2D data that has velocity as Y-axis
    # the increment in arcs along declination-axis

    # statistics
    dstd = d_arr.std()
    dmid = median(d_arr)
    dmin = d_arr.min()
    dmax = d_arr.max()
    up = arange(dmid+dstd,dmax,(dmid+dstd+dmax)/12)
    down = arange(abs(dmid-dstd),abs(dmin),(abs(dmid-dstd)+abs(dmin))/12)*(-1)
    levs = zeros(len(down)+len(up))
    levs[0:len(down)] = flipud(down)
    levs[len(down):len(down)+len(up)] = up
    levs= levs.round(1)
    # get the y-axis (declination) coordinates in arcs
    dec = ((arange(0,linedata.dec_npix,1)-linedata.dec_npix/2.)*linedata.dec_cdelt)[y1:y2]
    # get the x-axis (ra) coordinates in arcs
    ra = ((arange(0,linedata.ra_npix,1)-linedata.ra_npix/2.)*linedata.ra_cdelt)[x1:x2]

    fig1 = pl.figure(1)
    ax1 = pl.axes()
    # transpose it so that the velocity s on the X-axis
    channels = get_indices(linedata.v_arr,[v1,v2])
    im = ax1.contourf(linedata.v_arr[channels], dec, d_arr.transpose()[:, channels], levs)
    cb = pl.colorbar(im)
    pl.xlim(v1,v2) # kms boundaries
    pl.ylim(dec.min(),dec.max())
    pl.xlabel(r'$v_{lsr}$')
    pl.ylabel(r'asec')

    fig2 = pl.figure(2)
    ax2 = pl.axes()
    im = ax2.imshow(linedata.d[channels,y1:y2,x1:x2].sum(axis=0)*linedata.v_cdeltkms,
                    interpolation='nearest',cmap=cm.Greys)
    #im = ax2.pcolor(flipud(ra), dec, linedata.d[channels,y1:y2,x1:x2].sum(axis=0)*linedata.v_cdeltkms)
    cb = pl.colorbar(im)
    ax2.axis('image')
    pl.xlabel(r'asec')
    pl.ylabel(r'asec')
#
def plot_chmap (self,
                chvals=None,
                nvals=None,
                region=[0,0,0,0],
                nx=6,
                filled=True,
                box=[0,0],
                nsig=3,
                nsjump=2,
                source={'vsys':0},
                font={'family':'serif', 'serif': ['Times New Roman', 'Times'],
                'size':8},
                nsum=False,
                plot_adjust= [0.12, 0.09, 0.99, 0.99],
                quality=[300, 150],
                send=False,
                color_map='jet',
                cpeak = [0,0],
                locators = [0.5,0.1],
                levs_global=True,
                fsize=(3, 4)):
    # imports
    #import scipy as sp
    print 'importing...'
    from scipy import array, where, median, std, sign, arange, alen, vstack, \
    concatenate, sqrt, log10, exp, flipud, ceil
    import pyfits as pf
    import matplotlib
    import matplotlib.pyplot as pl
    from matplotlib.pyplot import contour, contourf
    if matplotlib.__version__<'1':
        from mpl_toolkits.axes_grid1 import AxesGrid
    else:
        from mpl_toolkits.axes_grid import AxesGrid
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from matplotlib import cm, rc
    print 'done'
    #


    ################################
    # setting the tick label font  #
    ################################
    # we are using latex formatting here!
    # the following set the ticklabel format string
    data_fmt = '%g'         # for normal y and x axis
    #freq_data_fmt = '%5.5f' # for the frequency array
    label_XY = parse_tick_font(font)
    tick_label_formatter = FormatStrFormatter(label_XY.replace('X',data_fmt))
    #freq_label_formatter = FormatStrFormatter(label_X.replace('X',freq_data_fmt))

    #rc('mathtext',**{'rm':'sans\\-serif'})


    ################################
    # parsing input                #
    ################################
    # set the subplot_adjust parameters, if not given a standard set will
    # be used
    pl_left, pl_bottom, pl_right,pl_top = plot_adjust
    #
    color_map = pl.get_cmap(color_map)

    #sigma levels etc
    nsig = nsig
    nsjump = nsjump
    # now check the type3 plots
    if chvals==None:
        ion=True
    elif chvals!=None:
        ion=False
        v1, v2 = chvals
    else: # input error
        raise ParError(chvals)
    #

    ################################
    #       get data               #
    ################################
    # first get the fitsfile
    #~ linedata = loadcube(filename)
    linedata = self


    # this line below is temporary, fix this in Fits.__init__ method
    #linedata.d = linedata.d[0]
    #linedata.v_arr = linedata.v_arr - source['vsys'] #so that now all the velocities are based on the objects
    #if cfile!=None:
    #    contdata= loadcube(cfile)
    cfile=None
    #
    # parse the region parameter
    x1, x2, y1, y2 = self.parse_region(region)

    if ion==True:
        plot_spectrum(filename, region=region)
        #
        # if it is just a channel map
        v1, v2 = array(raw_input('input the limits, comma separated: ').split(','), dtype ='float')
        if nvals==None:
            # ask for noise calculation velocity limits
            try:
                nvals = array(raw_input('input the noise limits (velocity). comma separated: ').split(','), dtype='float')
                if len(nvals)==4:
                    n1, n2, n3, n4 = nvals
                elif len(nvals)==2:
                    n1, n2 = nvals
            except (ValueError):
                print "Since you did not input any or input was wrong we will guess some values..."
                nvals = None
        pl.close(1)
    #######################################
    #######################################
    #

    # if we are to sum channels
    if nsum!=False:
        """
        blue[x:x+y].sum(axis=0)*abs(data.v_cdeltkms)
        calc_sigma(nchans,rms,v_cdelt)


        # give chvals=[v1,v2] and nsum > 1
        # nsum is number of channels to sum together
        # get the channels
        """
        print '\nParameter \'nsum\' is set to %d.' % nsum
        print stylify('Summing channels for better signal.',fg='g')
        # check if it is a integer
        if type(nsum) != type(1) or nsum<1: raise ParError(nsum)
        #
        channels = get_indices(linedata.v_arr, chvals)
        # can we divide it into the correct number of channels directly?
        print ' *Checking if we have to add some channels for the channel map to even up.'
        rest = alen(channels)%nsum
        print ' *%d channels, with nsum=%d we have a rest of %d' % (alen(channels),nsum,rest)
        to_add = nsum-rest
        if rest == 0:
            print '  Everything is fine, continuing.'
            pass
        elif to_add == 1:
            print '  Adding one channel'
            start = channels.min()
            end = channels.max()+to_add+1 # add one extra for python convention
        elif to_add%2 == 0:
            half = to_add/2
            print '  Going down/up %d in channel indices' % half
            start = channels.min()-half
            end = channels.max()+half+1 # add one extra for python convention
        elif to_add%2 != 0:
            half = (to_add-1)/2
            print '  Going down %d and up %d in channel indices' % (half, half+1)
            start = channels.min()-half
            end = channels.max()+(half+1)+1 # add one extra for python convention
        else:
            raise ParError(nsum)
        # now change the channel indices, so we get more channels
        if rest != 0:
            channels = arange(start,end,1)
        #
        # so now the channels are adapted to evenly divide with nsum
        #
        # getting the data to sum
        indices = arange(0, alen(channels), nsum)
        data = linedata.d[channels]
        velocity = linedata.v_arr[channels]
        # the mean velocity is the middle one
        velocity = array([velocity[x:x+nsum].mean() for x in indices])
        velocity_delta = linedata.v_cdeltkms*nsum
        print u'The velocity interval over which you create the maps is %2.3f to %2.3f km\u207b\u00b9\u00b7s' % (velocity.min()-abs(velocity_delta)/2,velocity.max()+abs(velocity_delta)/2)
        maps = array([data[x:x+nsum].sum(axis=0)*abs(linedata.v_cdeltkms) for x in indices])
        N_channels = alen(channels)/nsum
        print stylify('Done, remember that it can change the velocity interval that you specified in chvals. \n',fg='g')
        #
    else :
        # normal procedure here
        #
        #we are making a channel map
        channels = get_indices(linedata.v_arr, chvals)
        N_channels = alen(channels)
        velocity = linedata.v_arr[channels]
        velocity_delta = linedata.v_cdeltkms
        # units of Jy/beam kms
        data = linedata.d[channels]*abs(velocity_delta)
        maps = data
    #
    if nvals!=None:
        noise_channels = get_indices(linedata.v_arr,nvals)
    elif nvals==None:
        # if user didnt input noise velocity limits
        # choose channels away from the choosen (red blue) lines
        #
        # if you have choosen just a min and max to plot all channels
        low = where((linedata.v_arr>(linedata.v_arr.min()+10))*(linedata.v_arr<(v1-10)))[0]
        high = where((linedata.v_arr>(v2+10))*(linedata.v_arr<(linedata.v_arr.max()-10)))[0]
        noise_channels = concatenate((low, high))

    # calculate common stuff
    ### for the plotting box
    # *_crpix-1 because FITS starts at 1
    # changed to just crpix, because now loading with crpix-1
    # assumes that crpix is in the middle
    zlen, ylen, xlen = data.shape
    ycoords = arange(-(linedata.dec_crpix),(linedata.dec_npix-linedata.dec_crpix),1)*linedata.dec_cdelt
    xcoords = arange(-(linedata.ra_crpix),(linedata.ra_npix-linedata.ra_crpix),1)*linedata.ra_cdelt

    # for the extent keyword
    #~ left, right = xcoords[0],xcoords[-1]
    #~ bottom, top = ycoords[0],ycoords[-1]
    left,right,bottom,top = self.extent
    #
    # calculate the RMS
    noise = linedata.d[noise_channels]
    #
    # the region to calculate the RMS in
    i1,i2 = xlen/2+array([-1,1])*xlen/4
    j1,j2 = ylen/2+array([-1,1])*ylen/4
    #
    # RMS
    # since nsum might not be defined i.e. =False
    # dont really know if it is used later on,
    # so this is just to be sure
    # the rms used above changes with a factor of 1/sqrt(nsunm) for the
    # new channels in the channel map (since they are summed)
    nchansum = [1,nsum][nsum != False]
    rms = sqrt(((noise[:,j1:j2,i1:i2])**2).mean()/nchansum)
    #
    #
    if box == [0,0]:
        x1,x2 = left,right
        y1,y2 = bottom,top
    elif box != [0,0]:
        x1,x2 = array([-1,1])*box[0]/2.*sign(linedata.ra_cdelt)
        y1,y2 = array([-1,1])*box[1]/2.*sign(linedata.dec_cdelt)
    # calculate the correct rms and sigma
    # the rms is per channel, ie with the new
    # channel widht "velocity_delta"
    # so the sigma is just rms*abs(velocity_delta) ie. N=1
    sigma = calc_sigma(1,rms,velocity_delta)
    #
    chans_min = array([i.min() for i in maps])
    chans_max = array([i.max() for i in maps])

    # calculate the levels of the contours
    if levs_global:
        vmin, vmax = chans_min.min(),chans_max.max()
        # levs_stat, i.e. contourf starts and end at global min/max
        levs_stat = concatenate((
                    arange(vmin,-nsig*sigma,nsjump*sigma),
                    arange(nsig*sigma,vmax,nsjump*sigma)))
        levs = [levs_stat for x in arange(N_channels)]
    elif not levs_global:
        # levs, i.e. contour starts and ends at local min/max
        levs = [concatenate((arange(x,-nsig*sigma,nsjump*sigma), arange(nsig*sigma,y,nsjump*sigma))) for x,y in zip(chans_min, chans_max)]

    ##########
    majorLocator = MultipleLocator(locators[0])
    minorLocator = MultipleLocator(locators[1])

    print u'RMS \t: %2.2f mJy\u00b7beam\u207b\u00b9\u00b7channel\u207b\u00b9' % (rms*1e3)
    print u'1 sigma : %2.2f mJy\u00b7beam\u207b\u00b9\u00b7km\u00b7s\u207b\u00b9' % (sigma*1e3)
    print u'\u0394 vel \t : %2.3f km\u00b7s\u207b\u00b9 (channel width in km\u00b7s\u207b\u00b9)' % abs(velocity_delta)
    if send==True:
        if cfile!=None:
            return {'ch_map data': maps, 'levs':levs, 'chans_sig':sigma, 'vels':velocity, 'extent':self.extent, 'data':linedata, 'continuum':contdata}
        if cfile==None:
            return {'ch_map data': maps, 'levs':levs, 'chans_sig':sigma, 'vels':velocity, 'extent':self.extent, 'data':linedata}


    if velocity_delta<0:
        maps = flipud(maps)
        levs.reverse()
        velocity = flipud(velocity)
    #
    def print_vel(ax,x,fc='w'):
        ax.text(0.03,.91,
        str(round(velocity[x],1)),
        fontsize=6,
        bbox=dict(edgecolor='k',facecolor=fc, pad=4, lw=0.5),
        transform = ax.transAxes)


    ################################
    # setting global rc properties #
    ################################
    rc('savefig', **{'dpi': quality[0]})
    rc('text', usetex=1)
    rc('figure',**{'facecolor': '1', 'dpi': quality[1]})
    #to set the global font properties
    rc('font', **font)
    rc('axes',linewidth=0.7)
    rc('patch', linewidth=0.5)
    rc('lines', linewidth=0.3, markeredgewidth=0.6)
    #~ rc('font', family='serif', serif='Times New Roman', size=8)
    #~ rc('text', usetex=True)
    # ticksize
    #~ rc('xtick',**{'minor.size':2, 'major.size':4, 'major.pad': 3})
    rc('ytick',**{'minor.size':2, 'major.size':4, 'major.pad': -20})



    pl.ion()
    pl.close()
    fig = pl.figure(1, fsize)

    fig.clf()

    ny = int(ceil(N_channels/float(nx)))
    if filled==False:
        grid = AxesGrid(fig, 111, # similar to subplot(111)
                nrows_ncols = (ny, nx), # creates nyx6 grid of axes
                axes_pad=0, # pad between axes in inch.
                label_mode = "L",
                share_all=True,
                )
            # plot data contours
        for i in range(N_channels):
            grid[i].set_aspect('equal')
            im = grid[i].contour(maps[i],
                                levs_stat,
                                colors=('k'),
                                extent=self.extent)
    #
    elif filled==True:
        grid = AxesGrid(fig, 111, # similar to subplot(111)
                nrows_ncols = (ny, nx), # creates nx x ny grid of axes
                axes_pad=0, # pad between axes in inch.
                label_mode = "L",
                share_all=True,
                cbar_mode="single",
                cbar_location="top",
                cbar_size="2%",)
        # plot data contours
        for i in range(N_channels):
            grid[i].set_aspect('equal')
            im = grid[i].contourf(maps[i],
                                levs[i],
                                extent=self.extent,
                                cmap=color_map)
        grid.cbar_axes[0].colorbar(im)
        # add a unit to the colorbar
        grid.cbar_axes[0].set_xlabel(str(linedata.unit)+r'\,kms$^{-1}$')
        grid.cbar_axes[0].axis["top"].toggle(label=True, ticks=True, ticklabels=True)
    #
    if len(cpeak) == 3:
        mark = cpeak[2]
        xmark, ymark = cpeak[0:2]
        plot_cpeak = 1
    elif len(cpeak) >4:
        xmark = cpeak[0:-1:2]
        ymark = cpeak[1:-1:2]
        mark = cpeak[-1]
        plot_cpeak = 1
    else:
        # no cpeaks given
        mark = '+k'
        plot_cpeak = 0
        print('Not plotting cpeaks, got error here before')
        
    for i in range(N_channels):
        # plot a cross at xmark, ymark
        if plot_cpeak:
            grid[i].plot(xmark, ymark, mark, ms=6)#, mew=3, alpha=0.9)
        #~ cross2 = ax2ax2.plot(xmark, ymark, mark, ms=6)#  mew=3, alpha=0.9)
        
        #~ plt = grid[i].plot(cpeak[0],cpeak[1],'r+', ms=3, mew=0.7, alpha=0.7)
        # set the locator spacings
        #~ grid[i].xaxis.set_major_locator(majorLocator)
        #~ grid[i].xaxis.set_minor_locator(minorLocator)
        #~ grid[i].yaxis.set_major_locator(majorLocator)
        #~ grid[i].yaxis.set_minor_locator(minorLocator)

        #~ draw_fov(grid[i].axes, linedata)
        if i in [0,1,2]:
            print_vel(grid[i].axes, i,fc='#FFAAAA')
        elif i in [12,13,14]:
            print_vel(grid[i].axes, i,fc='0.8')
        else:
            print_vel(grid[i].axes, i)
    #~ fig.text(0.95,0.14,r'SO$_2$', rotation=90)
    #~ fig.text(0.95,0.86,r'H$_2^{18}$O', rotation=90)

    #~ grid[2].plot([21,-2], [-1.8,-1.3],'-', color='#FFAAAA', lw=3, alpha=0.7 ,clip_on=False)
    #~ grid[2].plot([21,-2], [-1.8-1.7,-1.3-1.7],'-', color='0.5', lw=3, alpha=0.7 ,clip_on=False)
    #~ grid[2].plot([11,-2], [-3.2,-2.9],'-', color='0.5', lw=3, alpha=0.7 ,clip_on=False)
    #~ grid[5].plot([21,-2], [-3.6,-3.2],'-', color='0.5', lw=3, alpha=0.7 ,clip_on=False)
    #~ grid[8].plot([21,-2], [-3.76,-2.53],'-', color='0.5', lw=3, alpha=0.7 ,clip_on=False)
    #~ grid[11].plot([21,-2], [-2.16,-0.93],'-', color='0.5', lw=3, alpha=0.7 ,clip_on=False)
    #~ grid[14].plot([21,-2], [-0.93,0.0],'-', color='0.5', lw=3, alpha=0.7 ,clip_on=False)
    #~ grid[12].axhline(y=-2.64, xmin=0.1, xmax=1.4 ,clip_on=False )

    draw_beam(grid[0].axes, linedata, box=1)

    #
    #grid.axes_llc.set_major_formatter(tick_label_formatter)
    #grid.axes_llc.set_major_formatter(tick_label_formatter)

    pl.subplots_adjust(left=pl_left, bottom=pl_bottom, right=pl_right, top=pl_top)


    grid.axes_llc.set_xlim(x1,x2)
    grid.axes_llc.set_ylim(y1,y2)

    fig.text(0.53,0.02,r'RA offset ($^{\prime\prime}$)', ha='center',va='center')
    fig.text(0.03,0.53,r'Dec offset ($^{\prime\prime}$)', rotation='vertical', va='center', ha='center')
    #fig.suptitle(linedata.obj)
    #pl.show()

    """
    from showmoment the 3range procedure.
    Use this to correct the binning!
        ### blue
        # if divided by 3, any rest?
        brest = len(blue)%3
        inc = (len(blue)-brest)/3
        bincrement = array([inc, inc, inc+brest])
        if data.v_cdelt>0:
            bincrement = flipud(bincrement)
        bindices = arange(0,(len(blue)-brest),bincrement[0])

        # create the first 2 maps
        bmaps = array([blue[x:x+y].sum(axis=0)*abs(data.v_cdeltkms) for x,y in zip(bindices,bincrement)])
        # calculate the center positions and their width
        bwidth = bincrement*abs(data.v_cdeltkms)*.5
        bcenter = [data.v_arr[x:x+y].mean() for x,y in zip(blue_channels[bindices],bincrement)]
        if data.v_cdelt>0:
            bmaps = flipud(bmaps)
            bwidth = flipud(bwidth)
            bcenter = flipud(bcenter)
        # calculate the sigmas for each set of data
        bsigma = calc_sigma(bincrement, rms, data.v_cdeltkms)

        ### red
        rrest = len(red)%3
        inc = (len(red)-rrest)/3

        # now flip the arrays (if v_cdelt<0) so that the highest speeds swallow the
        # most channels, if  multiple of 3 no channels
        rincrement = array([inc, inc, inc+rrest])
        if data.v_cdelt<0: # so that the last channel is still the one with the rest
            rincrement = flipud(rincrement)
        # the two middle indices
        # a bit tricky with the first being the biggest here
        rindices = array([0, rincrement[0], rincrement[0]+rincrement[1]])

        # get the maps, the last one (first if vcelt<0) takes the leftovers
        rmaps = array([red[x:x+y].sum(axis=0)*abs(data.v_cdeltkms) for x,y in zip(rindices,rincrement)])

        #get the widht & center of each channelsum
        # flipud so that the channels are the same in blue and red (ie low med high velocity)
        rwidth = rincrement*abs(data.v_cdeltkms)*.5
        rcenter = array([data.v_arr[x:x+y].mean() for x,y in zip(red_channels[rindices],rincrement)])
        if data.v_cdelt<0:
            rmaps = flipud(rmaps)
            rwidth = flipud(rwidth)
            rcenter = flipud(rcenter)

        rsigma = calc_sigma(rincrement,rms,data.v_cdeltkms)

        ### put them together now
        centers = concatenate((bcenter,rcenter)).round(2)
        widths = concatenate((bwidth,rwidth)).round(2)
        maps_sigma = concatenate((bsigma,rsigma))

        maps = vstack((bmaps,rmaps))

        maps_max = array([i.max() for i in maps])
        maps_min = array([i.min() for i in maps])

        # now create the levels for each image


        levs = [concatenate((arange(x,-nsig*z,nsjump*z),arange(nsig*z,y,nsjump*z))) for x,y,z in zip(maps_min, maps_max, maps_sigma)]
    """
#
def cubetracking (self,box=False, nsum=False):
    import matplotlib.pyplot as pl
    from matplotlib import cm
    from scipy import clip, array, sign, alen, arange
    pl.ion()
    class IndexTracker:
        def __init__(self, ax, data, velocity):
            self.ax = ax
            ax.set_title('Use scroll wheel to navigate images')
            self.X = data
            self.vel = velocity
            #self.data = dobject
            self.slices,cols,rows = data.shape
            self.ind  = self.slices/2
            self.im = ax.imshow(self.X[self.ind,:,:],interpolation='nearest',origin='lower', cmap=cm.jet)
            pl.show()
            self.update()

        def onscroll(self, event):
            #print event.button, event.step
            if event.button=='up':
                self.ind = clip(self.ind+1, 0, self.slices-1)
            elif event.button=='down':
                self.ind = clip(self.ind-1, 0, self.slices-1)
            self.update()

        def update(self):
            self.im.set_data(self.X[self.ind,:,:])
            ax.set_title(r'vel: '+str(round(self.vel[self.ind],2))+' kms$^{-1}$')
            self.im.axes.figure.canvas.draw()
    #
    pl.ioff()
    X = self
    if box!=False:
        i1,i2 = array([-1,1])*box[0]/2.*sign(self.ra_cdelt)
        j1,j2 = array([-1,1])*box[1]/2.*sign(self.dec_cdelt)
        i1, i2, j1, j2 = parse_region(self,[i1,i2,j1,j2])
        data = self.d[:, j1:j2, i1:i2]
    else:
        data = self.d
    if nsum!=False:
        if alen(data)%nsum!=0:
            raise ParError(nsum)
        channels = arange(0,alen(data),1)
        #
        # getting the data to sum
        indices = arange(0, alen(data), nsum)
        #data = linedata.d[channels]
        velocity = self.v_arr
        # the mean velocity is the middle one
        velocity = array([velocity[x:x+nsum].mean() for x in indices])
        velocity_delta = self.v_cdeltkms*nsum
        print u'The velocity interval over which you create the maps is %2.3f to %2.3f km\u207b\u00b9\u00b7s' % (velocity.min()-abs(velocity_delta)/2,velocity.max()+abs(velocity_delta)/2)
        print 'Velocity delta : %2.3f' % velocity_delta
        data = array([data[x:x+nsum].sum(axis=0)*abs(self.v_cdeltkms) for x in indices])
        N_channels = alen(data)/nsum
    else:
        velocity = self.v_arr
    fig = pl.figure(1)
    ax = fig.add_subplot(111)
    #ax = pwg.axes(header=X.hdr)
    tracker = IndexTracker(ax, data, velocity)
    #
    fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
#



########################################################################
# DEPRICATED
# to adavis.py(?, or delete?)
def plot_moment0 (self,
                cfile=None,
                chvals=None,
                nvals=None,
                region=[0,0,0,0],
                nx=6,
                filled=True,
                box=[0,0],
                nsig=3,
                nsjump=2,
                source = dict(dist=250),
                font={'family':'serif', 'serif': ['Times New Roman'],
                'size':8},
                fit = dict(gauss = None, params = None, continuum = False,
                    interactive = False),
                send=False,
                quality=[300, 300],
                cbar=True,
                colormap=True,
                plot_adjust= [0.13, 0.06, 0.75, 0.99],
                cpeak=[0,0,'k'],
                ccol='k',
                sbar=dict(au=200),
                locators = [2,1],
                telescope=None,
                negcontours=True,
                fsize=(_AandA.ONE_COL_FIG_WIDTH, _AandA.ONE_COL_FIG_WIDTH*0.8),
                rms_area=[0,0,10]):
    """

    Function doc

    source
        datatype : dictionary
        possible values : v_sys, dist, title
            v_sys - change the systemic velocity used
            dist - change the distanve used
            title - change the title text of the map drawn
    cpeak, continuum peak

    params = [height, amplitude, x0, y0, width_x, width_y, rota]
    TODO : For drawing size bar, add keywords to control it
    TODO : Change fit-routine to MPFIT
    TODO : Fix the units of the RMS/sensitivity/sigma

    """
    # imports
    #import scipy as sp
    from scipy import array, where, median, std, sign, arange, alen, vstack, \
    concatenate, sqrt, log10, ones, sqrt, flipud
    #import pyfits as pf
    from time import sleep
    import matplotlib.pyplot as pl
    #from matplotlib.patches import Circle
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from matplotlib import cm, rc
    #
    #
    #
    ################################
    # setting the tick label font  #
    ################################
    # we are using latex formatting here!
    # following sets the ticklabel format string
    data_fmt = '%g'                 # for normal y and x axis
    cbar_data_fmt = '%2.2f'         # for normal y and x axis
    label_X = parse_tick_font(font)
    tick_label_formatter = FormatStrFormatter(label_X.replace('X',data_fmt))
    cbar_tick_label_formatter = FormatStrFormatter(label_X.replace('X',cbar_data_fmt))
    #
    pl_left, pl_bottom, pl_right,pl_top = plot_adjust
    if filled==0:
        pl_right*=1.2
        pl_bottom*=2
    # parse the region, get the area for the spectra
    x1,x2,y1,y2 = self.parse_region(region)
    # if any is limits are None, draw spectra and ask for them/it
    # INTERACTIVE
    if hasattr(self,'rms'):
        pass
    elif nvals==None or chvals==None: # if self.rms does not exist and nvals not given
        # draw spectrum
        plot_spectrum(self,region=region, source= source)
        #
        chvals, nvals = get_vals(chvals=chvals, nvals=nvals)
        if nvals==None: # if no nvals where still not given...
            n1 = self.v_arr.min()+5*abs(self.v_cdeltkms)
            n2 = chvals[0]-5*abs(self.v_cdeltkms)
            n3 = chvals[1]+5*abs(self.v_cdeltkms)
            n4 =  self.v_arr.max()-5*abs(self.v_cdeltkms)
            nvals = [n1, n2, n3, n4]
        #
        self.calc_rms(nvals, rms_area)
    else: # if nvals was given and self.rms does not exist
        self.calc_rms(nvals, rms_area)
    #
    ### for the plotting box
    ylen, xlen = self.d[0].shape
    ycoords = arange(-ylen/2,ylen/2,1)*self.dec_cdelt
    xcoords = arange(-xlen/2,xlen/2,1)*self.ra_cdelt
    # for the extent keyword
    #~ left, right = xcoords[0],xcoords[-1]
    #~ bottom, top = ycoords[0],ycoords[-1]
    left,right,bottom,top = self.extent
    # set plot boundaries
    if box == [0,0]:
        i1,i2 = left,right
        j1,j2 = bottom,top
    elif box != [0,0]:
        i1,i2 = array([-1,1])*box[0]/2.*sign(self.ra_cdelt)
        j1,j2 = array([-1,1])*box[1]/2.*sign(self.dec_cdelt)

    #~ moments, [moment0_sigma, moment0_min, moment0_max], img_channels, levels = calc_moments(self, chvals=chvals, nvals=nvals, nsig=nsig, nsjump=nsjump, negcontours=negcontours, rms=rms)
    Mom = Moments(self, chvals=chvals, nsig=nsig)
    #
    d1,d2,r1,r2 = self.parse_region([0,0,box[0]])
    print 'No. sigmas for max in box: ', (Mom.zero[d1:d2,r1:r2].max()/Mom.sigma)

    # now create the levels for the contours
    if negcontours:
        #return flipud(arange(nsig,alen(Mom.levels_neg),nsjump)), arange(nsig,alen(Mom.levels_pos),nsjump), Mom
        #
        # not correct, at least not the negative levels
        #
        levs = concatenate((
                Mom.levels_neg[flipud(arange(nsig-1,alen(Mom.levels_neg),nsjump))],
                Mom.levels_pos[arange(nsig-1,alen(Mom.levels_pos),nsjump)]
                          ))
        levs_contour = concatenate((
                Mom.levels_neg[flipud(arange(nsig-1,alen(Mom.levels_neg),2*nsjump))],
                Mom.levels_pos[arange(nsig-1,alen(Mom.levels_pos),2*nsjump)]
                                  ))
    else:
        #~ levs = arange(nsig*img_sigma,img_max+img_sigma,nsjump*img_sigma)
        levs = Mom.levels_pos[arange(nsig,alen(Mom.levels_pos),nsjump)]
        #~ levs_contour = arange(nsig*img_sigma,img_max+img_sigma,2*nsjump*img_sigma)
        levs_contour = Mom.levels_pos[arange(nsig,alen(Mom.levels_pos),2*nsjump)]
    #levs = arange(nsig*img_sigma, img_max, nsjump*img_sigma)
    # print some info out
    print '\n','='*40
    print ' '*8,'INFORMATION : Line data'
    print '='*40
    print '\n Summing from channel %3d to %3d' % (Mom.channels.min(), Mom.channels.max())
    print ' RMS \t\t: %2.3f \tmJy/beam/channel\n Sigma \t\t: %2.3f \tmJy/beam/km/s' % (1e3*self.rms, 1e3*Mom.sigma)
    print ' Map min/max \t: %2.2f/%2.2f \tmJy/beam/km/s' % (1e3*Mom.minimum, 1e3*Mom.maximum)
    print ' Start sigma \t: %2.3f (%1.1f) \tmJy/beam/km/s\n Sigma step \t: %2.3f (%1.1f) \tmJy/beam/km/s\n' % (1e3*nsig*Mom.sigma, nsig, 1e3*nsjump*Mom.sigma, nsjump)

    #
    # tick density
    majorLocator = MultipleLocator(locators[0])
    minorLocator = MultipleLocator(locators[1])

    if send==True:
        print 'sending you moment class and extents'
        return Mom, (left,right,bottom,top)

    ################################
    # setting global rc properties #
    ################################
    set_rc()



    pl.ion()
    pl.close()
    fig = pl.figure(1, figsize=fsize)
    fig.clf()
    #
    ax = fig.add_subplot(111)
    #
    # print beam parameters for the data

    print '='*40
    print ' '*15,'BEAM(S)'
    print '='*40
    print 'Line cube:'
    print u' Minor (FWHM)\t: %2.3f \tasec' % self.bmin
    print u' Major (FWHM)\t: %2.3f  \tasec' % self.bmaj
    print ' PA \t\t: %2.3f \tDegrees (0<theta<180)' % self.bpa
    print u' Gain\t\t: %2.4f \tJy\u00b7K\u207b\u00b9\n' % self.gain
    #
    # plot the continuum data
    #~ if cfile != None:
        #~ cont = ax.imshow(contdata.d, cmap=cm.gray_r, extent=(left,right,bottom,top))
        #~ if cbar:
            #~ cb = pl.colorbar(cont)
            #~ cb.ax.set_ylabel(str(contdata.unit))
#~
        #~ print '-'*40
        #~ print 'Continuum data:'
        #~ print u' Minor (FWHM)\t: %2.3f \tasec' % contdata.bmin
        #~ print u' Major (FWHM)\t: %2.3f  \tasec' % contdata.bmaj
        #~ print ' PA \t\t: %2.3f Degrees (0<theta<180) \n' % contdata.bpa
    #
    #
    if colormap and filled==True:
        colormap = cm.bone_r
    if colormap == None or filled==False:
        # just contours
        levs = levs.round(3)
        #~ levs_contour = levs_contour.round(3)
        cs = ax.contour(Mom.zero, levs, colors=ccol, extent=self.extent)
    elif cfile == None and colormap!=None:
        levs = levs.round(3)
        levs_contour = levs_contour.round(3)
        cs1 = ax.contourf(Mom.zero, levs, cmap=colormap, extent=self.extent)
        #cs2 = ax.contour(img, cs1.levels[::2], colors=ccol, extent=self.extent)
        #return cs1
        cs2 = ax.contour(Mom.zero, levs_contour, colors=ccol, extent=self.extent)
        cbar = pl.colorbar(cs1, ticks=levs_contour, format=cbar_tick_label_formatter)#label_X.replace('X','%2.2f'))
        cbar.add_lines(cs2)
        if str(self.unit) == 'Jy/beam':
            cbar.ax.set_ylabel(r'Jy\,beam$^{-1}$')
        else:
            cbar.ax.set_ylabel(str(self.unit))
    else:
        line = ax.contour(Mom.zero, levels=levs, colors='r', extent=self.extent)

    #ax.text(0.5,0.5,'test',transform = ax.transAxes)
    draw_beam(ax, self)
    draw_fov(ax, self)

    # check distance key
    if sbar.has_key('dist'):
        dist_mark = sbar['dist']
    elif self.dist != 0:
        dist_mark = self.dist
    else:
        dist_mark = 200
    # check the length of the scale bar
    if sbar.has_key('au'):
        au_mark = sbar['au']
    else:
        au_mark = 200
    print 'Using distance {0} pc to source. Scale bar length {1} AU'.format(dist_mark, au_mark)
    draw_sizebar(ax, self, dist=dist_mark, au=au_mark)
    # parse the cpeak keyword
    if len(cpeak) == 3:
        mark = cpeak[2]
        xmark, ymark = cpeak[0:2]
    elif len(cpeak) >4:
        xmark = cpeak[0:-1:2]
        ymark = cpeak[1:-1:2]
        mark = cpeak[-1]
    else:
        mark = '+k'
    cross = ax.plot(xmark, ymark, mark, ms=6)#, mew=3, alpha=0.9)

    #
    ax.xaxis.set_major_formatter(tick_label_formatter)
    ax.yaxis.set_major_formatter(tick_label_formatter)

    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    #pl.setp(ax.get_xticklabels(),family='sans-serif')
    #pl.setp(ax.get_yticklabels(),family='sans-serif')

    ax.set_xlabel('RA Offset ($^{\prime\prime}$)')
    ax.set_ylabel('Dec Offset ($^{\prime\prime}$)')
    #fig.suptitle(self.obj+' (%s km s$^{-1}$)'% str(self.unit))
    #if cfile!=None:
    #    fig.subplots_adjust(bottom=0.08, right=0.77, top=0.90)
    if source.has_key('title'):
        ax.text(0.05,0.92, source['title'], transform = ax.transAxes)
    else:
        ax.text(0.05,0.92, self.obj, transform = ax.transAxes)
    ax.set_xlim(i1,i2)
    ax.set_ylim(j1,j2)
    ax.set_aspect(1)
    #
    #
    pl.subplots_adjust(left=pl_left, bottom=pl_bottom, right=pl_right, top=pl_top)

    if fit['gauss'] == '2d':
        #from gaussfitter2 import gaussfit
        """

        OLD FITTING ROUTINE, GAUSSFIT2D.PY
        SHOULD BE MOVED TO THE NEW ROUTINE MPFIT

        gaussfit(data,err=None,params=[],autoderiv=1,return_all=0,circle=0,rotate=1,vheight=1)
        params = [height, amplitude, x0, y0, width_x, width_y, rota]

        """
        print '\n-------------------------------\n2D Gaussian(s) fitting...\n'
        from scipy import sqrt, log, pi, tan, rot90, meshgrid

        #a = gaussfit(img) # fits xsig wrong in rel to x0
        if cfile != None and fit['continuum'] == True:
            print' Fitting to continuum data\n'
            Z = contdata.d[y1:y2,x1:x2]
            D = contdata
        else:
            print' Fitting to integrated line data'
            Z = img[y1:y2,x1:x2]
            D = self

        try:
            p = fit['params']
        except KeyError:
            p = None

        #
        # get the data coordinates, to fit over them!
        X, Y = meshgrid(array(xcoords)[x1:x2],array(ycoords)[y1:y2])

        #X = ones(Z.shape)*array(xcoords)[x1:x2]
        #Y = ones(Z.shape)*array(ycoords)[y1:y2]
        #return xcoords, ycoords, Z
        #print X,Y
        #u = indices(Z.shape)[1]
        #v = indices(Z.shape)[0]
        #
        # now do the fitting
        #
        params, success, no_fits = gaussfit2d((X,Y,Z), params=p)
        #
        # what we get in params is now
        # params[0] = the baseline
        # params[1] = peak
        # params[2] = x0
        # params[3] = y0
        # params[4] = xsig
        # params[5] = ysig
        #
        # some help functions to go to/from fwhm and sigma
        fwhmfromsig = 2*sqrt(2*log(2)) # the constant
        fwhm = lambda x: fwhmfromsig*x
        sigma = lambda x: x/fwhmfromsig
        #
        print ' Number of fits\t: %d' % no_fits
        print ' Success \t: %d \t(if over 4, not good)' % success
        print u' Baseline \t: %2.3e\tJy\u00b7beam\u207b\u00b9\u00b7km\u00b7s\u207b\u00b9\n' % params[0]
        #
        j = 1
        for i in arange(1,len(params),6):
            #
            #
            # parse the fit stuff to asec, maj/min axis etc
            # the errors are 0 here, for the future implementation of errors
            fwhm1, fwhm2, sfwhm1, sfwhm2, pa, spa = parse_gau2dfit(fwhm(params[i+3]),fwhm(params[i+4]),0,0,params[i+5],0)
            #
            # peak flux, S_peak
            s_peak = params[i]
            #
            # transform to pixel coordinates, not data coordinates
            # in data coordinates
            xoset, yoset  = (params[i+1], params[i+2])
            # get the pixel coordinates, needed to calc RA & DEC
            # f=True for returning the sub-pixel position
            x0, t1, y0,t2 = parse_region(D, [xoset, yoset, 0],f=True)
            # parse the coordinates into offsets
            #xoset, yoset = parse_pxlcoord(D, x0, y0)

            #xsigpxl = params[i+3]
            #ysigpxl = params[i+4]
            #
            # sigma
            sigmaj = sigma(fwhm1)
            sigmin = sigma(fwhm2)

            #
            # calculate the ra & dec coord
            # convert to degrees
            # ra_crpix is FITS based, i.e. starts at 1
            ra = (x0+1-D.ra_crpix)*D.ra_cdelt/3600 + D.ra_crval
            #ra = (D.ra_npix/2 - D.ra_crpix+1)*D.ra_cdelt/3600 + D.ra_crval + xoset
            a1,a2,a3 = parse_ra(ra)
            dec = (y0+1-D.dec_crpix)*D.dec_cdelt/3600 + D.dec_crval
            b1,b2,b3 = parse_dec(dec)
            #
            # the total flux (keep track of units)
            ftot = (s_peak*fwhm1*fwhm2/(D.bmaj*D.bmaj))

            print 'Fit no. : ', j
            print u' Flux \t\t: %2.5f \tJy\u00b7km\u00b7s\u207b\u00b9' % ftot
            print u' Peak \t\t: %2.5f \tJy\u00b7beam\u207b\u00b9\u00b7km\u00b7s\u207b\u00b9' % s_peak
            print ' RA Offset \t: %2.5f \tasec  RA: %02d:%2d:%2.3f' % (xoset, a1, a2, a3)
            print ' DEC Offset \t: %2.5f \tasec DEC: %02d:%2d:%2.3f' % (yoset, b1, b2, b3)
            print u' Major (FWHM)\t: %2.3f \tasec (\u03c3=%2.3f)' % (fwhm1, sigmaj)
            print u' Minor (FWHM)\t: %2.3f \tasec (\u03c3=%2.3f)' % (fwhm2, sigmin)
            print ' PA \t\t: %2.3f \tdegrees (-90<theta<+90) \n' % pa

            a, b, c, ok = gauss2d_decon ((fwhm1, fwhm2, pa, D.bmaj, D.bmin, D.bpa), ang='deg')

            print 'Deconvolved sizes:'
            print ' Major axes \t: %2.3f (asec)' % a
            print ' Minor axes \t: %2.3f (asec)' % b
            print ' PA \t\t: %2.3f (degrees)\n' % c
            j+=1
#
def plot_moment1 (self,
                cfile=None,
                chvals=None,
                nvals=None,
                region=[0,0,0,0],
                nx=6,
                filled=True,
                box=[0,0],
                nsig=3,
                nsjump=2,
                source = dict(v_sys=0),
                font={'family':'serif', 'serif': ['Times New Roman'],
                'size':17},
                fit = dict(gauss = None, params = None, continuum = False,
                interactive = False),
                send=False,
                quality=[150, 72],
                plot_adjust= [0.12, 0.01, 0.74, 0.99],
                cpeak=[0,0,'k'],
                ccol='r',
                sbar=dict(dist=250,au=200),
                locators = [2,1],
                telescope=None,
                type=1,
                fsize=(9.,7.),
                negcontours=True,):
    """
    Plot moment1 map of an area

    M1 = sum(vi * )



    ddir=substwid('~/')+"work/data/iram/h218o/"

    ; This block simply to read-in the data...
    va_l=[1.0]
    line=readfits(ddir+"line.fits",header)
    map_head,header,xa_l,ya_l,vaxis=va_l

    ; Select which part of the velocity axis to consider
    vx=where(va_l ge 6.0 and va_l le 8.0)

    ; Here we calculate the zeroth and first moments.
    dmom=fltarr(n_elements(xa_l),n_elements(ya_l),2)
    for i=0,n_elements(xa_l)-1 do begin
      for j=0,n_elements(ya_l)-1 do begin
        dmom(i,j,0)=total(line(i,j,vx)) ; Calculate moment-0.
        if dmom(i,j,0) ge 0.1           # if enough signal calculate moment-1
            then
            dmom(i,j,1)=total(reform(line(i,j,vx))*va_l(vx))/dmom(i,j,0)
        else
            dmom(i,j,1)=-99.0
      end
    end

    ; Do the plotting.
    set_plot,'ps'
    device,file=substwid("~/")+"work/ps/iram/iras4b/h218o_mom1.eps",xsize=6.0*0.8/0.6,ysize=6.0,/inch,/color,/cmyk & presentation
    loadct,6
    !P.POSITION=[0.15,0.15,0.75,0.95]
    lvls=0.15*findgen(61)/60.0+6.8

    ; The actual image and contours.
    contour,dmom(*,*,1),xa_l,ya_l,xrange=[1,-1],yrange=[-1,1],levels=lvls,/fill,c_color=250-findgen(61)*4.0,charsize=1.5,xtitle='!5RA offset ["]',ytitle='DEC offset ["]',xstyle=1,ystyle=1
    contour,dmom(*,*,0),xa_l,ya_l,xrange=[1,-1],yrange=[-1,1],levels=(0.05*findgen(20)+0.1),/over

    colorbar,ncolors=250,bottom=1,minrange=min(lvls),maxrange=max(lvls),position=[0.80,0.15,0.85,0.95],/vertical,charsize=1.5,divisions=6,FORMAT='(F5.2)',/right,/reverse
    """


    # imports
    #import scipy as sp
    from scipy import array, where, median, std, sign, arange, alen, vstack, \
    concatenate, sqrt, log10, ones, nan, flipud
    #import pyfits as pf
    import matplotlib.pyplot as pl
    #from matplotlib.patches import Circle
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from matplotlib import cm, rc
    #

    ################################
    # setting the tick label font  #
    ################################
    # we are using latex formatting here!
    # following sets the ticklabel format string
    data_fmt = '%g'         # for normal y and x axis
    cbar_data_fmt = '%2.2f'         # for normal y and x axis
    label_X = parse_tick_font(font)
    tick_label_formatter = FormatStrFormatter(label_X.replace('X',data_fmt))
    cbar_tick_label_formatter = FormatStrFormatter(label_X.replace('X',cbar_data_fmt))

    pl_left, pl_bottom, pl_right,pl_top = plot_adjust

    #
    velocity = self.v_arr

    # parse the channel values
    if chvals!=None:
        v1, v2 = chvals
    # parse the region, get the area for the spectra
    x1,x2,y1,y2 = self.parse_region(region)
    # if any is limits are None, draw spectra and ask for them/it
    # INTERACTIVE
    if nvals==None or chvals==None:
        # draw spectrum
        plot_spectrum(filename,region=region, source= {'vsys': source['vsys']})
        #
        # if there is no channels supplied
        #~ if chvals==None:
            #~ v1, v2 = array(raw_input('input the limits, comma separated: ').split(','), dtype ='float')
        #~ if nvals==None:
            #~ # ask for noise calculation velocity limits
            #~ try:
                #~ nvals = array(raw_input('input the noise limits (velocity). comma separated: ').split(','), dtype='float')
                #~ if len(nvals)==4:
                    #~ n1, n2, n3, n4 = nvals
                #~ elif len(nvals)==2:
                    #~ n1, n2 = nvals
            #~ except (ValueError):
                #~ print "Since you did not input any or input was wrong we will guess some values..."
                #~ nvals = None
        #~ pl.close(1)
        chvals, nvals = get_vals()
    # calculate common stuff
    ### for the plotting box
    ylen, xlen = self.d[0].shape
    ycoords = arange(-ylen/2,ylen/2,1)*self.dec_cdelt
    xcoords = arange(-xlen/2,xlen/2,1)*self.ra_cdelt
    #
    # for the extent keyword
    #~ left, right = xcoords[0],xcoords[-1]
    #~ bottom, top = ycoords[0],ycoords[-1]
    left,right,bottom,top = self.extent
    #
    # parse the noise channel values
    if nvals!=None:
        # if noise has been input
        noise_channels = get_indices(velocity, nvals)
    #
    else:
        # if user didnt input noise velocity limits
        # choose channels away from the choosen (red blue) lines
        low = where((velocity>(velocity.min()+10))*(velocity<(v1-10)))[0]
        high = where((velocity>(v2+10))*(velocity<(velocity.max()-10)))[0]
        noise_channels = concatenate((low, high))
    #
    # the region to calculate the rms in
    i1,i2 = xlen/2+array([-1,1])*xlen/4
    j1,j2 = ylen/2+array([-1,1])*ylen/4

    # the noise, rms
    noise = self.d[noise_channels]
    rms = sqrt(((noise[:,j1:j2,i1:i2])**2).mean())

    # set boundaries
    if box == [0,0]:
        i1,i2 = left,right
        j1,j2 = bottom,top
    elif box != [0,0]:
        i1,i2 = array([-1,1])*box[0]/2.*sign(self.ra_cdelt)
        j1,j2 = array([-1,1])*box[1]/2.*sign(self.dec_cdelt)


    Mom = Moments(self, chvals=chvals, nsig=nsig)

    moment0, moment1, [moment0_sigma, moment0_min, moment0_max], img_channels, levels\
    = Mom.zero, Mom.one, [Mom.sigma, Mom.minimum, Mom.maximum], Mom.channels, Mom.levels_pos
    #~ moments, [moment0_sigma, moment0_min, moment0_max], img_channels, levels


    #~ moments, [moment0_sigma, moment0_min, moment0_max], img_channels, levels
    #~
    #~
    #~ = calc_moments(self, chvals=chvals, nvals=nvals, nsig=nsig, nsjump=nsjump, negcontours=negcontours, rms=rms)

    #~ img_channels = get_indices(velocity,[v1,v2])
    imgs = self.d[img_channels]
    #
    # do the summing
    #
    # M1 = dv * sum(I(a,b,vi))_x1^x2
    # or perhaps without dv to save calc?
    #~ moment0 = imgs.sum(axis=0)*abs(self.v_cdeltkms)

    #~ moment0 = moments[0]

    moment0_sigma = sqrt(alen(imgs))*rms*abs(self.v_cdeltkms)
    moment0_max = moment0.max()
    moment0_min = moment0.min()



    # levels
    #~ if negcontours:
        #~ levels = concatenate((-1*flipud(arange(nsig*moment0_sigma,abs(moment0_min+moment0_sigma),nsjump*moment0_sigma)), arange(nsig*moment0_sigma,moment0_max+moment0_sigma,nsjump*moment0_sigma)))
    #~ else:
        #~ levels = arange(nsig*moment0_sigma,moment0_max+moment0_sigma,nsjump*moment0_sigma)
        #~
    # calculate moment 1 map

    #~ velocities = velocity[img_channels]
    #~ # calculate the denominator
    #~ Isum = imgs.sum(axis=0)
    #~ # create array that matches imgs array, so they can be multiplied
    #~ velocities_matrix = array([ones(imgs.shape[1:])*i for i in velocities])
    #~ # set values which are lower than 3 sigma to 'nan' i.e. mask the array
    #~ # use the moment0 to determine where, but the Isum to set it
    #~ Isum[where(moment0<(nsig*moment0_sigma))] = nan
    #~ # now, calculate the numerator of the division
    #~ Ivsum = (imgs*velocities_matrix).sum(axis=0)
    #~ # division
    #~ moment1 = Ivsum/Isum

    #~ moment1 = moments[1]

    ###
    ###
    ###
    # tick density
    majorLocator = MultipleLocator(locators[0])
    minorLocator = MultipleLocator(locators[1])
    if send==True:
        #~ return {'moment0': moment0,
                #~ 'moment1': moment1,
                #~ 'levels':levels,
                #~ 'moment0_sigma':moment0_sigma,
                #~ 'extent':self.extent,
                #~ 'lims': (i1,i2,j1,j2),
                #~ 'data':self}
        return moment0, moment1, levels, moment0_sigma, self.extent, self
    #
    #
    # calculating the velocity for vmin and vmax
    #
    #
    # remember it uses the region keyword, so
    print('remember that the calculation for vmin and vmax uses the region keyword')
    from string import lower
    arr = array([line for line in moment1[y1:y2,x1:x2].flatten() if str(line).lower() != 'nan'])
    w = arr.std()
    from scipy import median
    m = median(arr)
    """

    here you should add a histogram calculation, and fit a gaussian or something
    and take some representative width for vmin to vmax
    """

    set_rc(font=font, quality=quality)

    # linewidths
    rc('axes', linewidth=1)
    rc('lines', linewidth=1, markeredgewidth=1)

    pl.ion()
    pl.close()
    fig = pl.figure(1, figsize=fsize)
    fig.clf()
    #
    ax = fig.add_subplot(111)
    #
    # print beam parameters for the data
    gain = 8.168e-25*(self.restfreq)**2*self.bmin*self.bmaj
    print '='*40
    print ' '*15,'BEAM(S)'
    print '='*40
    print 'Line cube:'
    print u' Minor (FWHM)\t: %2.3f \tasec' % self.bmin
    print u' Major (FWHM)\t: %2.3f  \tasec' % self.bmaj
    print ' PA \t\t: %2.3f \tDegrees (0<theta<180)' % self.bpa
    print u' Gain\t\t: %2.4f \tJy\u00b7K\u207b\u00b9\n' % gain

    #levs=levs.round(2)
    #cs1 = ax.contourf(img, levels=levs, cmap=cm.bone_r, extent=self.extent)
    #cs2 = ax.contour(cs1, levels=cs1.levels[::2], colors=ccol, extent=self.extent)
    #im = pl.imshow(moment1,vmin=6.79,vmax=7,extent=self.extent)
    if type:
        im = pl.imshow(moment1,cmap=cm.jet,vmin=m-2*w,vmax=m+2*w,extent=(left,right,bottom,top),interpolation='nearest')
    elif not type:
        im = ax.contourf(moment1, levels=arange(m-2*w,m+2*w,4*w/20), cmap=cm.jet, extent=(left,right,bottom,top))
    cs2 = ax.contour(moment0, levels=levels, colors='k', extent=self.extent)
    cbar = pl.colorbar(im,format=cbar_tick_label_formatter) #label_X.replace('X','%2.2f'))

    cbar.ax.set_ylabel(r'km\,s$^{\sf -1}$')



    #draw_beam(ax, self,box=0)
    draw_fov(ax, self)

    if 'dist' and 'au' in sbar.keys():
        draw_sizebar(ax,self,dist=sbar['dist'],au=sbar['au'])
    elif 'dist' in sbar.keys():
        draw_sizebar(ax,self,dist=sbar['dist'])
    else :
        raise ParError(sbar)
    if len(cpeak) == 3:
        mark = cpeak[2]
        xmark, ymark = cpeak[0:2]
    elif len(cpeak) >4:
        xmark = cpeak[0:-1:2]
        ymark = cpeak[1:-1:2]
        mark = cpeak[-1]
    else:
        mark = '+k'
    cross = ax.plot(xmark, ymark, mark, ms=13, mew=3, alpha=0.9)

    #
    ax.xaxis.set_major_formatter(tick_label_formatter)
    ax.yaxis.set_major_formatter(tick_label_formatter)

    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    #pl.setp(ax.get_xticklabels(),family='sans-serif')
    #pl.setp(ax.get_yticklabels(),family='sans-serif')

    ax.set_xlabel('RA Offset ($^{\prime\prime}$)')
    ax.set_ylabel('Dec Offset ($^{\prime\prime}$)')

    #fig.suptitle(self.obj+' (%s km s$^{-1}$)'% str(self.unit))
    #if cfile!=None:
    #    fig.subplots_adjust(bottom=0.08, right=0.77, top=0.90)
    ax.text(0.05,0.93,
    self.obj,
    transform = ax.transAxes,
    backgroundcolor='w')
    ax.set_xlim(i1,i2)
    ax.set_ylim(j1,j2)
    ax.set_aspect(1)
    #
    #
    pl.subplots_adjust(left=pl_left, bottom=pl_bottom, right=pl_right, top=pl_top)
#


###################
#
# CMD line implementation
if __name__ == '__main__':
    # statements that you want to be executed only when the
    # module is executed from the command line
    # (not when importing the code by an import statement)
    # wee hooo...

    print '\nThe ADAVIS program was run from the command line\n'
    #import os
    from optparse import OptionParser as op
    #import pyfits as pf
    from string import lower

    #version
    ver = 1.0

    desc="""Add a keyword and value to fitsheader"""
    usage = "usage: %prog [options] fitsfile"
    epilog = """usage: %prog [options] fitsfile, \'fitsfile\' \n\
    should contain the relative path to the file as well"""

    parser = op(usage=usage,description=desc,epilog=epilog, version="%prog " + str(ver))

    # the options permitted
    parser.add_option("-f", "--function", dest="f", help="function to run on fits file" , metavar="FUNCTION")
    parser.add_option("-v", "--velocities", dest="v", help="between what velocities to plot", metavar="VELOCITIES")
    parser.add_option("-b", "--box", dest="b", help="zoom to a box X,Y size", metavar="BOX")
    parser.add_option("-l", "--list", action="store_true", dest="l", help="list the header", metavar="LIST")
    parser.add_option("-k", "--force", action="store_true", dest="k", help="force overwriting of keyword", metavar="FORCE")

    # time to parse
    (options, args) = parser.parse_args()

    # if no filename is supplied
    if len(args) != 1:
        parser.error("Incorrect number of arguments")
        parser.print_usage()
    # the options
    f = str(options.f)

    def convert_to_list(a):
        return [float(x) for x in (a).split(',')]
    #
    #
    # create a dictionary!
    # have to change the functions, to handle **kwargs first
    # the everything except infile and type of plotting function can be **kwargs
    #

    if options.v != None:
        v = convert_to_list(options.v)
    else:
        v = options.v
    if options.b != None:
        b = convert_to_list(options.b)
    else:
        b = options.b
    listhdr = options.l
    force = options.f
    # and the file...
    filename = str(args[0])
    print v, b
    if lower(f) in ['spectrum', 'sectrum', 'stectrum', 'spectraum', 'spectarum']:
        # tries to pick up miss-spellings
        plot_spectrum(filename)

    if lower(f) in ['moment0','mom0','m0']:
        if b!=None:
            plot_moment0(filename,box=b)
        else:
            plot_moment0(filename)
    if lower(f) in ['moment1','mom1','m1']:
        if b!=None:
            plot_moment1(filename,box=b)
        else:
            plot_moment1(filename)


    #
    a = str(raw_input('Press any key to exit...'))
    # put Jes, optparser here!












