#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
#       magpy.py
#
#       Copyright 2010 Magnus Persson <magnusp@snm.ku.dk>
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
#
#       ver 1.0
#
#
#

"""
Script with functions to perform different actions on
our data.


######################################################
#                                                    #
#  move all suggestions to the respective function   #
#                                                    #
######################################################
    for the main funcitons that is

TODO : Function to bin in spatial (2D) regime (congrid map), change ra/dec_delt etc

TODO : P-V diagram - rotatable

TODO : Moment 2 maps - with MPFIT Gaussian fitting

TODO : Calculate common stuff in a function, use in moment 0/1/2 maps

TODO : What does the VELREF keyword in the header mean? check GILDAS 'fits' manual

TODO : Clean up code again, remove font handler, or inactivate it

TODO : Tick locators are good for small regions/narrow spectra, but not for big/wide
        -Perhaps change it with the 'box' keyword

TODO : Separate the loadcube function, so that one can work with data
       objects instead of a function that allways reads the file

Lastly:
- Check Jes functions how they work and try to merge/replace them.
- How to divide it into a runnable program
- implement it as command line program




"""

#import os
from sys import exit as sysexit
#import scipy as sp
from scipy import arange, array, indices, flipud, concatenate, ceil, alen,\
                clip, pi, sqrt
#import pyfits as pf
from platform import system
if system() == 'Darwin':
    print 'Where on OSX so setting backend to MacOSX'
    # must import and run use before pyplot/pylab
    import matplotlib
    matplotlib.use('MacOSX')

from matplotlib import cm, rc
import matplotlib.pyplot as pl
#import pywcsgrid2 as pwg
#import matplotlib.pyplot as pl
#from mpl_toolkits.axes_grid import AxesGrid

#from matplotlib.patches import Circle
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from string import upper

from matplotlib import rc

rc('image',**{'origin':'lower', 'interpolation':'bilinear'})
rc('savefig', **{'dpi': 300})
rc('text', usetex=True)
rc('figure',**{'facecolor': '1', 'dpi': 72})

# do not know if this is needed, check if eps export is rasterised
# need python-poppler installed
#rc('ps.usedistiller' : 'xpdf')

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
###########################################
# HELP FUNCTIONS


def loadcube(fitsfile, telescope=None):
    """
    TODO : loadcube - if the rotational matrice is non empty, load data and rotate it
    TODO : loadcube -what if add a function to grid it to a certain size, say 512x512?
            o so it is possible to combine data with different res.
            o better to do this when plotting?
    TODO : loadcube - More robust loading (detecting the axis etc)?
            o Minimum of different types of data
                - Naxis 2 (Freq/Vel & Flux/Intensity) - 1D Spectra
                - Naxis 2 (RA & DEC) - 2D map
                - Naxis 3 (RA & DEC & Flux/Intensity) 3D Spectral map
                - Polarization data?
    TODO : make loadcube delete an axis, along with the hdr keywords if all
           the axis keywords/values are empty/null
                o only to loaded data, not save to raw data (fits file)
    """


    #imports
    from pyfits import open as fitsopen
    #from  sys import
    from os.path import getsize
    from scipy import where, array
    #

    # create the class, but without any init script.
    # a class (object) the easy way
    class data: pass
    print u'Loading fitsfile :  %s ' % colorify(str(fitsfile),c='g')
    s  = getsize(fitsfile)
    print " Size %0.2f MB" % (s/(1024*1024.0))
    f = fitsopen(fitsfile)
    data.hdr, data.d = f[0].header, f[0].data
    #data.d = data.d[0] # this is if the stokes axis is present,
    # but it should not be there anymore
    f.close()

    #
    #
    # the telescope diameter
    # first check if there was keyword sent in
    if telescope!=None:
        data.hdr.update('TELESCOP', telescope)
    #
    if data.hdr.has_key('TELESCOP') == 1:
        name = array(['SMA', 'PDB', 'JCMT', 'AP-H201-F102'])
        dia = array([6, 15, 15, 12])
        try:
            data.diameter = dia[where(upper(data.hdr['TELESCOP'])==name)][0]
        except IndexError, ex:
            data.diameter = 1
    else:
        data.diameter= 1
    if data.hdr.has_key('LINE')==1:
        data.linename = data.hdr['LINE']
    #
    #
    # getting the velocity
    def loadvelocity(data, axno):
        data.veltype = axno
        data.vcrpix = data.hdr['CRPIX'+data.veltype]
        data.vcrval = data.hdr['CRVAL'+data.veltype]
        data.vctype = data.hdr['CTYPE'+data.veltype]
        data.vcdelt = data.hdr['CDELT'+data.veltype]
        data.vnaxis = data.hdr['NAXIS'+data.veltype]
        data.vcdeltkms = data.vcdelt/float(1e3)
        data.velarr = ((arange(1,data.vnaxis+1)-data.vcrpix)*data.vcdelt+data.vcrval)/float(1e3) # so it is i kms
        data.velrangekms = data.velarr.max()-data.velarr.min()
        # not good if vcdletkms has more than 3 significant digits.
        #data.velarr = data.velarr.round(3)
        #data.vcdeltkms = round(data.vcdeltkms,2)

        #start = data.vcrval-data.vcdelt*(data.vcrpix-1)
        #stop =  data.vcrval+data.vcdelt*(data.vnaxis-data.vcrpix)
        #arr = arange(start,stop-1,data.vcdelt)/float(1e3)
        #print data.velarr-arr
        # calculate the FOV = 58.4*lambda/D*3600 asec
        data.restfreq = data.hdr['RESTFREQ'] # in Hertz
        data.fov = 58.4*(3.e8/data.restfreq)/float(data.diameter)*3600.
        print 'Field of view: %.2f asecs, for dish size: %.1f m' % (data.fov, data.diameter)
        #print data.veltype, data.vcrpix, data.vcrval, data.vcdeltkms, data.vnaxis
        print 'Velocity range \t: %d' % data.velrangekms
        return data
    # spectra: 3 axis and 3rd axis is >1 in size
    # image: 2 axis (or 3 axis, 3rd is =1 in size)
    #
    from numpy import diff, arange
    #naxis = data.hdr['NAXIS']
    #axshape = data.d.shape
    #if axshape[0] array([i>1 for i in a.shape[1:]]).all()
    # an image, only 2 dimensions
    if data.hdr['NAXIS']==2 and data.hdr['NAXIS1']>1 and data.hdr['NAXIS2']>1:
        #if data.hdr['NAXIS']==3 and data.hdr['NAXIS1']>1 and data.hdr['NAXIS2']>1 and data.hdr['NAXIS3']==1:
        # image, not SD spectra or anything,
        # really 2D and greater extent than 1x1
        data.type = 'IMAGE'
        pass
    #
    # spectral image cube (extra axis for frequency/velocity)
    elif data.hdr['NAXIS']==3 and data.hdr['NAXIS1']>1 and data.hdr['NAXIS2']>1 and data.hdr['NAXIS3']==1:
        data.type = 'IMAGE'
        # extra if the continuum image has the freq and width
        data.freq = data.hdr['CRVAL3']
        data.freqwidth = data.hdr['CDELT3']
        # remove the extra axis in the data
        data.d = data.d[0]
    # a spectra! the 3rd axis is longer than 1
    elif data.hdr['NAXIS']>=3 and data.hdr['NAXIS3']>1:
        # spectral cube
        # only support for velo-lsr in 3rd axis
        data.type = 'CUBE'
        # load the third axis
        velax = str([x for x in data.hdr.keys() if x[:-1]=='CTYPE' and 'VELO' in data.hdr[x]][0][-1:])
        data = loadvelocity(data, velax)
        # now if we want to have the spectral array as well to use
        if data.hdr.has_key('RESTFREQ'):
            data.restfreq = data.hdr['RESTFREQ']
        if data.hdr['NAXIS']==4 and  data.hdr['NAXIS4']==1:
            data.d = data.d[0]
        #
        # this was just a test
        # SD pointing spectra
    elif data.hdr['NAXIS']>1 and data.hdr['NAXIS2']==1 and data.hdr['NAXIS3']==1:
        data.type = 'SDSPECT'
        data.vcdelt = data.hdr['DELTAV']
        data.vcdeltkms = data.hdr['DELTAV']/float(1e3)
        data.vcrpix = data.hdr['CRPIX1']
        data.vnaxis = data.hdr['NAXIS1']
        data.vcrval = data.hdr['VELO-LSR']
        data.velarr = ((arange(1,data.vnaxis+1)-data.vcrpix)*data.vcdelt+data.vcrval)/float(1e3)
        data.restfreq = data.hdr['RESTFREQ'] # in Hertz
        data.fov = 58.4*(3e8/data.restfreq)/(data.diameter)*3600
        data.d = data.d[0][0][0] # specific for this data...


    #
    else:
        # if it is not an image or a spectral cube
        print '\n ERROR\nThe dimensions of the data is wrong\n at least the header keywords indicate that.\n The data has '+str(data.hdr['NAXIS'])+' axes. \n\n Perhaps use the removeaxis script?\n\n'
        sysexit(1)

        #
        # FREQUENCY ARRAY
        #
        # construct the frequency array!
        # the 3rd axis longer than 1, and 4th axis is the frequency
        # if the data is constructed in gildas
        #

    # DEC
    data.dec_cdelt = data.hdr['CDELT2']*3600 # arcs
    data.dec_npix = data.hdr['NAXIS2']
    data.y_npix = data.hdr['NAXIS2']
    data.dec_crpix = data.hdr['CRPIX2']
    data.dec_crval = data.hdr['CRVAL2']
    # RA
    data.ra_cdelt = data.hdr['CDELT1']*3600 # arcs
    data.ra_npix = data.hdr['NAXIS1']
    data.x_npix = data.hdr['NAXIS1']
    data.ra_crpix = data.hdr['CRPIX1']
    data.ra_crval = data.hdr['CRVAL1']
    try:
        # Beam size in asecs
        data.bmaj = data.hdr['BMAJ']*3600
        data.bmin = data.hdr['BMIN']*3600
        data.bpa = data.hdr['BPA']
    except KeyError, ex:
        data.bmaj = None
        data.bmin = None
        data.bpa = None
    try:
        # Data units
        data.unit = data.hdr['BUNIT']
        data.obj = data.hdr['OBJECT']
    except Exception, ex:
        data.unit = None
    #
    # Object name
    data.obj = data.hdr['OBJECT']

    print '\n','='*40
    print ' '*8,'INFORMATION : FITS file'
    print '='*40
    print 'Object : %s' % data.obj
    data.ra_size = abs(data.ra_cdelt)*data.ra_npix
    data.dec_size = abs(data.dec_cdelt)*data.dec_npix
    print 'Spatial size of image\n RA\t: %2.3f asec\n DEC\t: %2.3f asec' % (data.ra_size, data.dec_size)
    ra_h, ra_m, ra_s = parse_ra(data.ra_crval)
    dec_h, dec_m, dec_s = parse_dec(data.dec_crval)
    print 'Phase center '
    print ' RA : %2d:%2d:%2.4f' % (ra_h,ra_m,ra_s)
    print ' DEC : %2d:%2d:%2.4f' % (dec_h,dec_m,dec_s)
    #
    print('done')
    return data
#
def calc_frequency(vlsr, freq0):
    """
    vlsr in kms
    freq0 in whatever you want out
    """
    from scipy import constants
    return (1-vlsr/(constants.c*1e-3))*freq0
#
def calc_vlsr (f,f0):
    """ calc vlsr in km/s from two freq """
    from scipy import constants
    return (1-f/f0)*constants.c*1e-3
#
def draw_beam(ax, data,loc=3, box=True):
    """
    function that draws the beam
    the attributes data.bmin, .bmaj and .bpa must exist in the data class
    """
    from mpl_toolkits.axes_grid.anchored_artists import AnchoredEllipse
    # the PA is calculated from N to E, that is +90 degrees from normal
    ae = AnchoredEllipse(transform=ax.transData,\
                            width=data.bmin,\
                            height=data.bmaj,\
                            angle=data.bpa+90,\
                            loc=loc,\
                            pad=0.15,\
                            borderpad=0.2,\
                            frameon=box)

    ax.add_artist(ae)
#
def colorify (txt, c='r'):
    """
    Sends back the string 'txt' with the correct unicode color
    start and finish.

    Current avaliable colors:
        'r' : red
        'g' : green
        'b' : blue
    """
    CSI = "\x1B["

    if c=='r':
        # sets, red text color (31) and black background (40)
        start =CSI+'31m'
    elif c=='g':
        # sets, green text color (32) and black background (40)
        start =CSI+'32m'
    elif c=='b':
        # sets, blue text color (34) and black background (40)
        start =CSI+'34m'
    #
    end = CSI+'m'
    return start+txt+end
#
def draw_fov(ax, data):
    """
    Function to draw the field of view into the
    plot with axes 'ax'.
    Assumes that the axes data is set to arcsec
    """
    from matplotlib.patches import Circle
    cir = Circle( (0,0), transform=ax.transData, fill=False, ec='k', lw=1, ls='dashed', radius=data.fov/2)
    ax.add_patch(cir)
#
def draw_sizebar(ax, data, dist=220, au=200):
    """
    distance in pc
    then theta = AU/(distance in pc)
    """
    from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar
    # draw horizontal bar with length of 10 in data coordinate
    # 10 arcs * 220pc = 2200
    asb = AnchoredSizeBar(ax.transData,
                            au/float(dist),
                            str(au)+" AU",
                            loc=8,
                            pad=0.1, borderpad=0.5, sep=5,
                            frameon=False)
    ax.add_artist(asb)
#
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
    print xpos, xwidth
    rect = Rectangle((xpos,0), width=xwidth, height=1,
                             transform=trans, color='yellow',
                             alpha=0.3)

    ax.add_patch(rect)
#
def put_line_indicator(ax, velocity, spect, xpos, text_string, \
            text_size='medium', text_weight='extra bold', text_color='black', \
            lw=3, lc='b', offset=0):
    """
    Put text in figure, the positon is in axes data coordinates
    additional kwargs are text_size, text_weight, text_color
    oset_bool is for lines at same places
    """
    # calculate position
    line_height = spect.max()*0.08
    try:
        maxval = spect[get_indices(velocity,[xpos-2,xpos+2])].max()
    except ValueError, ex:
        maxval = spect.max()
    except IndexError, ex:
        maxval = spect.max()
    line_ypos = maxval + maxval*0.15
    text_ypos = line_ypos + line_height + line_height*0.5
    text_ypos = text_ypos + len(text_string)*offset*3e-3
    # plot line
    ax.plot([xpos, xpos],[line_ypos,line_ypos+line_height],lc,\
            lw=lw)
    # print text
    ax.text(xpos+0.1,text_ypos,\
    text_string,\
    size=text_size, weight=text_weight, color=text_color,\
    ha='center',va='bottom',rotation='vertical',\
    transform = ax.transData)
#
def calc_sigma(N,rms,vcdelt):
    return sqrt(N)*rms*abs(vcdelt)
#
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
    from scipy import array, where
    if font.has_key('family'):
        fformatters = array(['$\mathrm{X}$','$\mathsf{X}$','$\mathtt{X}$' , '$\mathit{X}$'])
        ffamilies = array(['serif', 'sans-serif', 'monospace', 'cursive'])
        label_X = fformatters[where(font['family']==ffamilies)][0]
    else:
        label_X = '$\mathrm{X}$'
    #return label_X
    return 'X'
#
def parse_region(data, region, f=False):
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

    """
    from scipy import ceil, floor
    if len(region)==4:
        xcheck = region[0]==region[2]
        ycheck = region[1]==region[3]
        #x1, x2 = (data.ra_npix+1)/2 + array([region[0],region[2]])/abs(data.ra_cdelt) + array([0,xcheck])
        #y1, y2 = (data.dec_npix+1)/2+ array([region[1],region[3]])/abs(data.dec_cdelt)+ array([0,ycheck])
        #
        x1, x2 = array([region[0],region[2]])/data.ra_cdelt+ data.ra_crpix - 1 + array([0,xcheck])
        y1, y2 = array([region[1],region[3]])/data.dec_cdelt+ data.dec_crpix - 1 + array([0,ycheck])
        #
    elif len(region)==3:
        check = region[2]==0
        #x1, x2 = (data.ra_npix+1)/2 + array([-region[2],region[2]])/(2*abs(data.ra_cdelt)) + region[0]/data.ra_cdelt + array([0,check])
        #y1, y2 = (data.dec_npix+1)/2+ array([-region[2],region[2]])/(2*abs(data.dec_cdelt)) +region[1]/data.dec_cdelt+ array([0,check])
        #
        x1, x2 = data.ra_crpix + region[0]/data.ra_cdelt + array([-region[2],region[2]])/abs(2*data.ra_cdelt) - 1 + array([0,check])
        y1, y2 = data.dec_crpix + region[1]/data.dec_cdelt + array([-region[2],region[2]])/abs(2*data.dec_cdelt) - 1 + array([0,check])
        #
    elif len(region)==2:
        xcheck = region[0]==0
        ycheck = region[1]==0
        #x1, x2 = (data.ra_npix+1)/2 + array([-1,1])*region[0]/abs(data.ra_cdelt)  + array([0,xcheck])
        #y1, y2 = (data.dec_npix+1)/2+ array([-1,1])*region[1]/abs(data.dec_cdelt) + array([0,ycheck])
        #
        x1, x2 = array([-region[0],region[0]])/(2*abs(data.ra_cdelt)) + data.ra_crpix - 1 + array([0,xcheck])
        y1, y2 = array([-region[1],region[1]])/(2*abs(data.dec_cdelt)) + data.dec_crpix - 1 + array([0,ycheck])
        #
    elif():
        print ('Error, region keyword malformed')
        sysexit(1)
        #
    # so that we are returning usable pixel coordinates
    if f==False:
        x1,x2,y1,y2 = array([x1,x2,y1,y2]).round().astype('int')
    return x1,x2,y1,y2
#
def parse_pxlcoord (data, x, y):
    """ Function doc """
    xoffset = (x-data.ra_crpix)*data.ra_cdelt
    yoffset = (y-data.dec_crpix)*data.dec_cdelt
    return xoffset, yoffset
#
def parse_ra (ra):
    """

    Parses a simple float coordinate and returns hours, minutes and seconds.
    Input
        ra : the right ascention coordinate to be parsed

    ---------------------------------------------------------------------------

                            oOO Changelog OOo

    *2010/09 Funciton created

    *2010/12/18 Added documentation


    """
    a = ra/15
    hours = int(a)
    b = (a-hours)*60
    minutes = int(b)
    seconds = (b-minutes)*60
    return hours,minutes,seconds
#
def parse_dec (dec):
    """

    Parses a simple float coordinate and returns hours, minutes and seconds.
    Input
        dec: the declination coordinate to be parsed

    ---------------------------------------------------------------------------

                            oOO Changelog OOo

    *2010/09 Funciton created

    *2010/12/18 Added documentation


    """
    degrees = int(dec)
    b = (dec-degrees)*60
    minutes = int(b)
    seconds = (b-minutes)*60
    return degrees,minutes,seconds
#
def parse_linelist(linelist):
    """

    Parses a linelist/CSV file of the lines, so that it is easier to just
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

    """
    from scipy import size, arange, zeros
    def get_lines(linelist):
        # pars a list of lists scenario
        names = linelist[0]
        freqs = linelist[1]
        arr_len = len(names)+len(freqs)
        list = zeros((arr_len))
        return array([[names[i],freqs[i]] for i in arange(len(names))]).ravel()
    #
    if size(linelist[1])==1 and type(linelist)==type([]):
        # if it is just a list of lines
        # it is already in the disired format
        return linelist
    elif size(linelist[1])>1:
        # if it is a list with two lists in it
        # that is, the second element is (strictly) longer than 1
        return get_lines(linelist)
    elif size(linelist[1])==1 and type(linelist)==type(''):
        # in case it is the path to a CSV file
        # load the table with the load ascii table function loadatbl
        names, freqs = loadatbl(linelist, dtype='string', sep=':')[0:3:2]
        #
        freqs = freqs.astype('float64')
        return get_lines([names,freqs])
#
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


    ---------------------------------------------------------------------------

                            oOO Changelog OOo

    *2010/06
        Funciton created

    *2010/12
        Doc written

    """

    from scipy import concatenate, where, array, diff
    dx = abs(.5*diff(arr)[0])
    if len(vals)==4:
        v1,v2,v3,v4 = vals + array([-1,1,-1,1])*dx
        # if the user wants two velocity areas to calculate noise
        low = where((arr>=v1)*(arr<v2))[0] + 1
        high = where((arr>=v3)*(arr<v4))[0] + 1
        channels = concatenate((low,high))
    elif len(vals)==2:
        v1,v2 = vals + array([-1,1])*dx
        #  if input just want one velocity area to calculate noise
        channels = where((arr>=v1)*(arr<v2))[0] + 1
    #
    if disp and len(vals)==2:
        first, last = channels.min(), channels.max()
        n = last-first+1
        print '\nFirst: %d,\n Last: %d\n Nchan: %d\n' % (first, last, n)
    return channels
#
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
    from numpy import exp, array, alen, array, log, sqrt
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
        gaussian = lambda X, a : a[0]*exp(-(X-a[1])**2/(a[2]**2)*sqrt(2*log(2)))
    elif height!=None:
        gaussian = lambda X, a : height + a[0]*exp(-(X-a[1])**2/(a[2]**2)*sqrt(2*log(2)))
    #
    def func (X, p):
        S=0
        try:
            for i in arange(0,len(p),3):
                S += gaussian(X, p[i:i+3])
        except IndexError, ex:
            raise ParError(p)
        return S
    #
    #
    return func(x, params)
#
def fit_gauss1d((X,Y),\
                params,\
                err = None,\
                fixlist = None,\
                minbool = None,\
                minpar = None,\
                maxbool = None,\
                maxpar = None,\
                tie = None,\
                verbose = 1,\
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
            minimum
    minpar - so then you must specify what that minimum limit is
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
        7  xtol is too small. No further improvement in the approximate solution
           x is possible.
        8  gtol is too small. fvec is orthogonal to the columns of the jacobian
           to machine precision.


    TODO : move code over to MPFIT
    TODO : fixing of parameters
    TODO : limits to parameters
    TODO : initial guesses for params
                    o intereactive - point and clic
                    o 1D : add a 1 gaussian guessing alorithm
    TODO : minimum width >= X[1]-X[0]
    """
    #from scipy import optimize
    from numpy import exp, hstack, array, log, sqrt, diag, alen, zeros, \
                        where
    from mpfit import mpfit
    print '='*40
    print '\n Fitting Gaussians'
    #
    ## Checking the input parameters
    #
    print '\nChecking input parameters'
    # flatten the parameters, so it is readable for the gaussian fitting
    params = array(params).flatten()
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
        # it is not limited but,
        # here we still would want to limit the gaussian fwhm to the
        # width of one x axis data unit, can't be thinner than that, can it
        minbool = [False,False,True]*no_fits
        minpar = [0,0,xwidth]*no_fits
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
    # NB: we are fitting FWHM directly
    gaussian = lambda X, a : a[0]*exp(-(X-a[1])**2/(a[2]**2)*sqrt(2*log(2)))
    #
    def fitfunc (X, p):
        S=0
        for i in arange(0,len(p),3):
            S += gaussian(X, p[i:i+3])
        return S
    #
    def errfunc(x,y,err):
        if err == None:
            def f(p,fjac=None): return [0,(y-fitfunc(x,p))]
        else:
            def f(p,fjac=None): return [0,(y-fitfunc(x,p))/err]
        return f
    #
    # define the parameter dictionary
    PAR = ['Amplitude','Position','Fwhm']
    fitlist = []
    for i in xrange(alen(params)):
        # create the dictionary to be appended to the list=
        dictline = dict(n  = i,\
                        value   = params[i],\
                        limits  = [minpar[i],maxpar[i]],\
                        limited = [minbool[i],maxbool[i]],\
                        fixed   = fixlist[i],\
                        parname = PAR[i%3],\
                        tied = tie[i],\
                        err   = 0)
        fitlist.append(dictline)


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
#
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
#
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
#
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
#
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
#
def set_rc(font={'family':'sans-serif', 'sans-serif': ['Arial', 'Helvetica'],
        'size':22, 'weight':'bold'},
        quality=[150, 72]):

    from matplotlib import rc
    ################################
    # setting global rc properties #
    ################################
    rc('text', usetex=True)
    rc('savefig', **{'dpi': quality[0]})
    rc('figure',**{'facecolor': '1', 'dpi': quality[1]})
    #to set the global font properties
    rc('font', **font)
    # ticksize
    rc('xtick',**{'minor.size':3, 'major.size':7})
    rc('ytick',**{'minor.size':3, 'major.size':7})

    # linewidths
    rc('axes', linewidth=2)
    rc('lines', linewidth=1.5, markeredgewidth=1)
#


###########################################
# MAIN FUNCTIONS

def plot_spectrum (filename,
                chvals=None,
                nvals=None,
                region=[0,0,0,0],
                obj_dict = dict(vsys=0),
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
    Function documentation

    Required input:
    filename
        name of fits-file input and if needed the path to it

    Optional parameters:

    chvals = [a, b]

    nvals = [a, b]

    region = [x1]

    obj_dict = {'vsys' : 0}

    freq = True/False (1/0)

    font = {'family':'sans-serif', 'sans-serif': ['Arial', 'Helvetica'],'size': 20, 'weight':'bold'}

    bin = True/False (1/0)

    fit = {'type' : 'gauss', 'params' : ((a1, v1, w1),(a2, v2, w2)),
        'guess' : True/False, interactive : , fixlist : , error : ,
        limmin : , limmax : ,minpar :, limmax : , maxpar :, tie : }

    send = True/False


    quality = [a, b]
        set the quality of the saved figure (a) and the displayed figure (b)
        in dpi normal is quality = [150, 72], for prints, quality = [300, 72]?

    plot_adjust = [a, b, c, d]
        set the left, bottom, right, top of the axes, if none is given the
        standard [0.08, 0.07, 0.98, 0.92] is set.

    lines = ['a', b, 'c', d]
        draw lines with name 'a' & 'c' and rest frequencies b & d.
        if obj_dict=dict(vsys=x) is not null, then it will shift the positions
        of the lines accordingly.

    axspace = []



    TODO : what is the RMS when binning and giving a area_region>1?
            ie does the end unit change? yes-> fix it
                    # change x1,x2,y1,y2 to quarter region (line 2350)
                    # change so that when binning, the rms i calculated
    """
    # imports
    #import scipy as sp
    print 'importing...'
    from scipy import array, where, median, std, sign, arange, alen, vstack, \
                    concatenate, sqrt, log10, exp, log, ceil, floor, diff
    import matplotlib.pyplot as pl
    from mpl_toolkits.axes_grid import AxesGrid
    #from matplotlib.patches import Circle
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from matplotlib import cm, rc, rc_params
    print 'done'
    #
    #font={'family':'serif','serif':['Times', 'Palatino', 'New Century Schoolbook', 'Bookman', 'Computer Modern Roman'],'size':12, 'weight':'bold'}
    #font=['serif','Times',12,'bold'],\
    #font = {'family':'sans-serif','sans-serif':['Arial', 'Helvetica', 'Avant Garde', 'Computer Modern Sans Serif'], 'cursive':['Zapf Chancery']}
    #'monospace':['Courier','Computer Modern Typewriter'], },\
    #from matplotlib.font_manager import FontProperties
    #FontProperties(family=None, style=None, variant=None, weight=None, stretch=None, size=None, fname=None, _init=None)
    ################################
    # setting global rc properties #
    ################################
    rc('savefig', **{'dpi': quality[0]})
    rc('text', usetex=True)
    rc('figure',**{'facecolor': '1', 'dpi': quality[1]})
    #to set the global font properties
    rc('font', **font)
    rc('axes',linewidth=2)
    rc('lines', linewidth=2)
    ################################
    # setting the tick label font  #
    ################################
    # the following set the ticklabel format string
    data_fmt = '%g'         # for normal y and x axis
    freq_data_fmt = '%5.5f' # for the frequency array
    label_X = parse_tick_font(font)
    tick_label_formatter = FormatStrFormatter(label_X.replace('X',data_fmt))
    freq_label_formatter = FormatStrFormatter(label_X.replace('X',freq_data_fmt))
    #rc('mathtext',**{'rm':'sans\\-serif'})
    # ticksize
    rc('xtick',**{'minor.size':4, 'major.size':9, 'major.pad': 6})
    rc('ytick',**{'minor.size':4, 'major.size':9, 'major.pad': 4})
    # linewidths
    rc('axes', linewidth=2)
    rc('lines', linewidth=1, markeredgewidth=2)


    # set the subplot_adjust parameters, if not given a standard set will
    # be used
    pl_left, pl_bottom, pl_right,pl_top = plot_adjust

    if chvals!=None:
        v1, v2 = chvals

    # first get the fitsfile
    data = loadcube(filename)
    data.velarr = data.velarr - obj_dict['vsys'] #now all the velocities are based on the LSR
    #data.freqarr = calc_frequency(data.velarr,data.restfreq)
    #
    # parse the region parameter
    print region
    x1,x2,y1,y2 = parse_region(data, region)
    print x1,x2,y1,y2
    # now start the plotting
    area_region = ((y2-y1)*(x2-x1))
    if data.type != 'SDSPECT':
        spect = (data.d[:,y1:y2,x1:x2].sum(axis=1).sum(axis=1))/(area_region)
    elif data.type == 'SDSPECT':
        spect = data.d
    #
    # binning of data
    if bin!=False:
        # change to congridding?
        from congridding import congrid
        j = int(bin)

        #  Old method - simple binning, just average
        #indices = arange(0,alen(spect),j)
        #spect = array([spect[x:x+j].sum(axis=0)/j for x in indices])
        #velocity = array([data.velarr[x:x+j].sum(axis=0)/j for x in indices])
        #
        # congridding, proper re gridding of data
        #
        spect = congrid(spect,(alen(spect)/bin,),centre=True)
        velocity = congrid(data.velarr,(alen(data.velarr)/j,))
        #
        velocity_delta = data.vcdeltkms*bin
        #
    elif bin==False:
        velocity = data.velarr
        velocity_delta = data.vcdeltkms


    #
    # From here on the velocity array is 'velocity'
    # and the vcdelt is 'velocity_delta'

    #
    pl.ion()
    pl.close()
    fig = pl.figure(1,figsize=(10.5,8))
    #fig = pl.figure(1,figsize=(10.5,4.5))
    fig.clf()
    ax_kms = fig.add_subplot(111)
    # now correct so that negative values are to the left
    cx1 = axspace[0]
    cx2 = axspace[1]
    cy1 = axspace[2]
    cy2 = axspace[3]
    if data.vcdelt<0:
        ax_kms.step(flipud(velocity), flipud(spect), '0.3', lw=2,  where='mid')
        xmin, xmax = round(flipud(velocity)[0])*cx1, round(flipud(velocity)[-1])*cx2
        ymin, ymax = round(min(flipud(spect)),3)*cy1, round(max(flipud(spect)),3)*cy2
    elif data.vcdelt>=0:
        ax_kms.step(velocity, spect, 'k', lw=2, where='mid')
        xmin, xmax = round(velocity[0])*cx1,round(velocity[-1])*cx2
        ymin, ymax = round(min(spect),3)*cy1,round(max(spect),3)*cy2
    if nvals!=None:
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

        #rms_0= rms/sqrt(abs(data.vcdeltkms))
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
        ax_kms.plot([xmin, xmax],[3*rms,3*rms],'b--', lw=2, alpha=0.7)
        ax_kms.plot([xmin, xmax],[-3*rms,-3*rms],'b--', lw=2, alpha=0.7)
    #
    # plot dotted lines at y=0 and x=0
    ax_kms.plot([xmin-10,xmax+10],[0,0],'k:',[0,0],[ymin-1,ymax+1],'k:')
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

    if fit['type']=='gauss':
        """
        TODO : calculate errors for the fit etc, look at Kaper et al (1966)\
                and Condon (1996).
        TODO : guess parameters?
        TODO : calculate mean and std of the FWHM for the fits
        """
        if not fit.has_key('error'):
            errmsg='You have to supply an error to the fit, they\'re kind of important you know.'
            print colorify(errmsg)
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
        # now, parse output of fitting and print it out on screen
        j = 1
        line_widths = []
        for i in arange(0,len(params),3):
            # add 1 because channel 1 is in pos 0
            half_fwhm = params[i+2]/2
            fwhm = params[i+2]
            line_widths.append(fwhm)
            # first figure out the extent of the gaussian (the line)
            # jump half a channel down and up so that it finds the correct channels
            lower_half, upper_half = (params[i+1] + array([-1,1])*half_fwhm)
            lower,upper = (params[i+1] + array([-1,1])*fwhm)
            #
            #channels = where((data.velarr>lower)*(data.velarr<upper))[0]+1
            channels_half = get_indices(data.velarr,[lower_half,upper_half])
            channels = get_indices(data.velarr,[lower,upper])

            print 'Fit number : %i' % j
            print ' Intensity : %2.4f \t=calculate it=' % (sqrt(2*pi)*sigma(params[i+2])*params[i]) # the area under the 1D Gaussian
            print u' Amplitude : %2.3f (\u00b1%2.3f) \t Jy\u00b7Beam\u207b\u00b9' % (params[i],errors[i])
            print u' Position  : %2.3f (\u00b1%2.3f) \t km\u00b7s\u207b\u00b9' % (params[i+1],errors[i+1])
            print u' Width     : %2.3f (\u00b1%2.3f) \t km\u00b7s\u207b\u00b9 (FWHM, \u03c3=%2.3f)' % (params[i+2],errors[i+2], sigma(params[i+2]))
            print ' Frequency : %3.9f GHz' % calc_frequency(params[i+1],data.restfreq/1e9)
            print ' Channels  : %d to %d (2FWHM width : %d to %d (%d) [%.2f to %.2f km/s])\n' % (channels_half.min(), channels_half.max(),channels.min(), channels.max(),(channels.max()-channels.min()+1), lower,upper)
            j+=1
        #
        line_widths = array(line_widths)
        print 20*'- '
        print u'Mean FWHM : %2.1f \u00b1%2.2f km\u00b7s\u207b\u00b9' % (line_widths.mean(),line_widths.std())
        if send:
            j = 1
            f = []
            for i in arange(1,len(params),3):
                nu = calc_frequency(params[i],data.restfreq/1e9)
                f.append(nu)
                print '%3.9f' % nu
                j+=1
        # draw the fit into the figure
        xarr = arange(X[0],X[-1],(diff(X)[0]/3))
        ax_kms.plot(xarr, gauss1d(xarr,params), color='#dd0000', lw=3, alpha=0.8)
    #
    if lines!=[]:
        print u'Marking the lines, using %2.2f km\u00b7s\u207b\u00b9' % obj_dict['vsys']
        if type(lines[-1]) in [type(1),type(1.0)]:
            # to be able to add a frequency shift to the linelist
            # to move the markings
            print ('this is very crude implemented, change/remove if it conflicts.')
            lines_vobs = lines[-1]
            lines = lines[0]
        else:
            # so that we do not have to think about if we created it later
            lines_vobs = 0
        lines = parse_linelist(lines)
        v = array([calc_vlsr(float(lines[i+1])*1e9,data.restfreq) for i in arange(0, len(lines), 2)]) + lines_vobs
        colors8 = ['#95B200','#2A3702','#71AEDB','#6D001D','#4B8618', '#DDB61D', '#DC3A0C', '#003B73']
        colors7 = ['#95B200','#71AEDB','#6D001D','#4B8618', '#DDB61D', '#DC3A0C', '#003B73']
        if len(lines)/2 == 8:
            colors = colors8
        elif len(lines)/2 == 7:
            colors = colors7
        elif len(lines)/2 not in [6, 8]:
            colors = ['k']*len(lines)
        x = 0
        for i,j in zip(arange(0, len(lines), 2), arange(0, len(lines),1)):
            # only check for lines behind the new line
            no_dbl = len(where(array(v)[0:j].round(0) == round(v[j],0))[0])
            put_line_indicator(ax_kms, velocity, spect, v[j], lines[i],lc=colors[x], offset=no_dbl, text_color=colors[x], lw=6)
            #put_line_indicator(ax_kms, velocity, spect, v[j], lines[i],lc='k', offset=no_dbl)
            print u'Line %1d : %2.2f\t km\u00b7s\u207b\u00b9 (%2.3f)' % (j+1,v[j],v[j]+obj_dict['vsys'])
            x+=1
    #
    ax_kms.set_xlim(xmin,xmax)
    if lines!=[]:
        ax_kms.set_ylim(ymin,ymax+ymax*0.4)
    else:
        ax_kms.set_ylim(ymin,ymax)
    if ylimits!=None:
        ax_kms.set_ylim(ylimits)
    #ax_kms.set_title(data.obj)
    ax_kms.set_xlabel('$v_{\mathrm{lsr}}$ [km s$^{-1}$]')
    if data.unit=='Jy/beam':
        ax_kms.set_ylabel('$I$ [Jy beam$^{-1}$]')
    else:
        ax_kms.set_ylabel('$I$ ['+data.unit+']')
    # to show the channel numbers
    if freq==True:
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
        ax_hz.set_xlim(calc_frequency(x_1,data.restfreq/1e9), calc_frequency(x_2,data.restfreq/1e9))
    #
    #elif lines!=[] and obj_par['vsys']==0:
    #    print('please, input a vsys!=0')
    #    return ; sysexit()
    ax_kms.xaxis.set_major_formatter(tick_label_formatter)
    ax_kms.yaxis.set_major_formatter(tick_label_formatter)

    #ax_kms.set_xticklabels(ax_kms.get_xticks(),ffont)
    #ax_kms.set_yticklabels(ax_kms.get_yticks(),ffont)
    ax_kms.minorticks_on()
    #for tick in ax_kms.xaxis.get_major_ticks():
    #    tick.label.set_family('Arial')
    #pl.show()
    #pl.draw()
    pl.subplots_adjust(left=pl_left, bottom=pl_bottom, right=pl_right, top=pl_top)

    #pl.show()

    #pl.ioff()

    if send:
        if nvals!=None and fit['type']!=None: # send everything!
            txt = '\n sending you spectra, data, noise-spectra, axes instance and the fitted frequencies'
            print colorify(txt,c='g')
            return spect, data, spect[get_indices(data.velarr,nvals)], ax_kms, f
        elif nvals==None and fit['type']=='gauss': # no noise channels supplied
            txt = '\n sending you spect, data, axes instance and the fitted frequencies'
            print colorify(txt,c='g')
            return spect, data, ax_kms, f
        elif nvals!=None and fit['type']==None: # no fit but noise
            txt =  '\n sending you spect, data, noise-spectra and axes instance'
            print colorify(txt,c='g')
            return spect, data, spect[get_indices(data.velarr,nvals)], ax_kms
        else: # well non of them are supplied
            txt =  '\n sending you spectra, data, axes instance'
            print colorify(txt,c='g')
            return spect, data, ax_kms
        return
#
def cubetracking (f):
    import matplotlib.pyplot as pl
    pl.ion()
    class IndexTracker:
        def __init__(self, ax, X, data):
            self.ax = ax
            ax.set_title('Use scroll wheel to navigate images')
            self.X = X
            self.data = data
            self.slices,cols,rows = X.shape
            self.ind  = self.slices/2
            self.im = ax.imshow(self.X[self.ind,:,:])
            pl.show()
            self.update()

        def onscroll(self, event):
            #print event.button, event.step
            if event.button=='up':
                self.ind = clip(self.ind+1, 0, self.slices-1)
            else:
                self.ind = clip(self.ind-1, 0, self.slices-1)
            self.update()

        def update(self):
            self.im.set_data(self.X[self.ind,:,:])
            ax.set_title(r'vel: '+str(round(self.data.velarr[self.ind],2))+' kms$^{-1}$')
            self.im.axes.figure.canvas.draw()


    pl.ioff()
    X = loadcube(f)
    fig = pl.figure()
    ax = fig.add_subplot(111)
    #ax = pwg.axes(header=X.hdr)
    tracker = IndexTracker(ax, X.d, X)

    fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
#
def plot_moment0 (filename,
                cfile=None,
                chvals=None,
                nvals=None,
                region=[0,0,0,0],
                nx=6,
                filled=True,
                box=[0,0],
                nsig=3,
                nsjump=2,
                obj_dict = dict(vsys=0,),
                font={'family':'serif', 'serif': ['Times', 'Times New Roman'],
                'size':22, 'weight':'bold'},
                fit = dict(gauss = None, params = None, continuum = False,
                interactive = False),
                send=False,
                quality=[150, 72],
                cbar=True,
                plot_adjust= [0.12, 0.01, 0.74, 0.99],
                cpeak=[0,0,'k'],
                ccol='r',
                sbar=dict(dist=220,au=200),
                locators = [2,1],
                telescope=None):
    """

    Function doc

    params = [height, amplitude, x0, y0, width_x, width_y, rota]
    TODO : For drawing size bar, add keywords to control it
    TODO : Change fit-routine to MPFIT
    TODO : Fix the RMS/sensitivity/sigma calculation units
    cpeak, continuum peak
    """


    # imports
    #import scipy as sp
    from scipy import array, where, median, std, sign, arange, alen, vstack, \
    concatenate, sqrt, log10, ones
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
    data_fmt = '%g'         # for normal y and x axis
    cbar_data_fmt = '%2.2f'         # for normal y and x axis
    label_X = parse_tick_font(font)
    tick_label_formatter = FormatStrFormatter(label_X.replace('X',data_fmt))
    cbar_tick_label_formatter = FormatStrFormatter(label_X.replace('X',cbar_data_fmt))
    # ticksize
    rc('xtick',**{'minor.size':3, 'major.size':7})
    rc('ytick',**{'minor.size':3, 'major.size':7})

    # linewidths
    rc('axes', linewidth=1)
    rc('lines', linewidth=1.5, markeredgewidth=1)

    pl_left, pl_bottom, pl_right,pl_top = plot_adjust

    #
    # first get the fitsfile
    linedata = loadcube(filename,telescope)
    #
    if not obj_dict.has_key('vsys'):
        obj_dict['vsys'] = 0
    linedata.velarr = linedata.velarr - obj_dict['vsys'] #now all the velocities are based on the LSR

    # and the continuum data, if existent
    if cfile != None:
        contdata = loadcube(cfile)
    # parse the channel values
    if chvals!=None:
        v1, v2 = chvals
    # parse the region, get the area for the spectra
    x1,x2,y1,y2 = parse_region(linedata, region)
    # if any is limits are None, draw spectra and ask for them/it
    # INTERACTIVE
    if nvals==None or chvals==None:
        # draw spectrum
        plot_spectrum(filename,region=region, obj_dict= {'vsys': obj_dict['vsys']})
        #
        # if there is no channels supplied
        if chvals==None:
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

    # calculate common stuff
    ### for the plotting box
    ylen, xlen = linedata.d[0].shape
    ycoords = arange(-ylen/2,ylen/2,1)*linedata.dec_cdelt
    xcoords = arange(-xlen/2,xlen/2,1)*linedata.ra_cdelt
    #
    # for the extent keyword
    left, right = xcoords[0],xcoords[-1]
    bottom, top = ycoords[0],ycoords[-1]
    #
    # parse the noise channel values
    if nvals!=None:
        # if noise has been input
        noise_channels = get_indices(linedata.velarr, nvals)
    #
    else:
        # if user didnt input noise velocity limits
        # choose channels away from the choosen (red blue) lines
        low = where((linedata.velarr>(linedata.velarr.min()+10))*(linedata.velarr<(v1-10)))[0]
        high = where((linedata.velarr>(v2+10))*(linedata.velarr<(linedata.velarr.max()-10)))[0]
        noise_channels = concatenate((low, high))
    #
    # the region to calculate the rms in
    i1,i2 = xlen/2+array([-1,1])*xlen/4
    j1,j2 = ylen/2+array([-1,1])*ylen/4

    # the noise, rms
    noise = linedata.d[noise_channels]
    rms = sqrt(((noise[:,j1:j2,i1:i2])**2).mean())

    # set boundaries
    if box == [0,0]:
        i1,i2 = left,right
        j1,j2 = bottom,top
    elif box != [0,0]:
        i1,i2 = array([-1,1])*box[0]/2.*sign(linedata.ra_cdelt)
        j1,j2 = array([-1,1])*box[1]/2.*sign(linedata.dec_cdelt)

    # get the channels
    #img_channels = where((linedata.velarr>v1)*(linedata.velarr<v2))[0]
    img_channels = get_indices(linedata.velarr,[v1,v2])
    imgs = linedata.d[img_channels]
    #
    # do the summing
    #
    # M1 = dv * sum(I(a,b,vi))_x1^x2
    #
    img = imgs.sum(axis=0)*abs(linedata.vcdeltkms)
    img_sigma = calc_sigma(alen(imgs), rms, linedata.vcdeltkms)
    #
    img_max = img.max()
    img_min = img.min()
    #
    # now create the levels for the contours
    levs = concatenate((-1*flipud(arange(nsig*img_sigma,abs(img_min)+img_sigma,nsjump*img_sigma)), arange(nsig*img_sigma,img_max+img_sigma,nsjump*img_sigma)))
    levs_contour = concatenate((-1*flipud(arange(nsig*img_sigma,abs(img_min)+img_sigma,2*nsjump*img_sigma)), arange(nsig*img_sigma,img_max+img_sigma,2*nsjump*img_sigma)))
    print levs
    print levs_contour
    #levs = arange(nsig*img_sigma, img_max, nsjump*img_sigma)
    # print some info out
    print '\n','='*40
    print ' '*8,'INFORMATION : Line data'
    print '='*40
    print '\n Summing from channel %3d to %3d' % (img_channels.min(), img_channels.max())
    print ' RMS \t\t: %2.3f \tmJy/beam/channel\n Sigma \t\t: %2.3f \tmJy/beam/km/s' % (1e3*rms, 1e3*img_sigma)
    print ' Map min/max \t: %2.2f/%2.2f \tmJy/beam/km/s' % (1e3*img_min, 1e3*img_max)
    print ' Start sigma \t: %2.3f (%1.1f) \tmJy/beam/km/s\n Sigma step \t: %2.3f (%1.1f) \tmJy/beam/km/s\n' % (1e3*nsig*img_sigma, nsig, 1e3*nsjump*img_sigma, nsjump)
    #
    # tick density
    majorLocator = MultipleLocator(locators[0])
    minorLocator = MultipleLocator(locators[1])

    if send==True:
        if cfile!=None:
            print 'sending you things'
            return {'img': img, 'levs':levs, 'imgs_sigma':img_sigma, 'extent':(left,right,bottom,top), 'data':linedata, 'continuum':contdata}
        if cfile==None:
			print 'sending you : img, levs, imgs_sigma, extent, data'
			return {'img': img, 'levs':levs, 'imgs_sigma':img_sigma, 'extent':(left,right,bottom,top), 'data':linedata}

    pl.ion()
    pl.close()
    if cbar and cfile!=None:
        fig = pl.figure(1, (9., 7.))
    else:
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
    #
    #
    if cfile != None:
        cont = ax.imshow(contdata.d, cmap=cm.gray_r, extent=(left,right,bottom,top))
        if cbar:
            cb = pl.colorbar(cont)
            cb.ax.set_ylabel(str(contdata.unit))

        print '-'*40
        print 'Continuum data:'
        print u' Minor (FWHM)\t: %2.3f \tasec' % contdata.bmin
        print u' Major (FWHM)\t: %2.3f  \tasec' % contdata.bmaj
        print ' PA \t\t: %2.3f Degrees (0<theta<180) \n' % contdata.bpa
    #
    #
    if cfile == None:
        levs = levs.round(3)
        levs_contour = levs_contour.round(3)
        cs1 = ax.contourf(img, levs, cmap=cm.bone_r, extent=(left,right,bottom,top))
        #cs2 = ax.contour(img, cs1.levels[::2], colors=ccol, extent=(left,right,bottom,top))
        #return cs1 
        cs2 = ax.contour(img, levs_contour, colors=ccol, extent=(left,right,bottom,top))
        cbar = pl.colorbar(cs1, ticks=levs_contour, format=cbar_tick_label_formatter)#label_X.replace('X','%2.2f'))
        cbar.add_lines(cs2)
        if str(linedata.unit) == 'Jy/beam':
            cbar.ax.set_ylabel(r'Jy\,beam$^{-1}$')
        else:
            cbar.ax.set_ylabel(str(linedata.unit))
    else:
        line = ax.contour(img, levels=levs, colors='r', extent=(left,right,bottom,top))
	
	#ax.text(0.5,0.5,'test',transform = ax.transAxes)
    draw_beam(ax, linedata)
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
    ax.set_ylabel('Dec offset ($^{\prime\prime}$)')
    #fig.suptitle(linedata.obj+' (%s km s$^{-1}$)'% str(linedata.unit))
    #if cfile!=None:
    #    fig.subplots_adjust(bottom=0.08, right=0.77, top=0.90)
    if obj_dict.has_key('source'):
        ax.text(0.05,0.93, obj_dict['source'], transform = ax.transAxes)
    else:
        ax.text(0.05,0.93, linedata.obj, transform = ax.transAxes)
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
            D = linedata

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
def plot_moment1 (filename,
                cfile=None,
                chvals=None,
                nvals=None,
                region=[0,0,0,0],
                nx=6,
                filled=True,
                box=[0,0],
                nsig=3,
                nsjump=2,
                obj_dict = dict(vsys=0),
                font={'family':'serif', 'serif': ['Times', 'Times New Roman'],
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
    # we are using latex formatting here!
    # following sets the ticklabel format string
    data_fmt = '%g'         # for normal y and x axis
    cbar_data_fmt = '%2.2f'         # for normal y and x axis
    label_X = parse_tick_font(font)
    tick_label_formatter = FormatStrFormatter(label_X.replace('X',data_fmt))
    cbar_tick_label_formatter = FormatStrFormatter(label_X.replace('X',cbar_data_fmt))


    pl_left, pl_bottom, pl_right,pl_top = plot_adjust


    #
    # first get the fitsfile
    linedata = loadcube(filename,telescope)
    #
    linedata.velarr = linedata.velarr - obj_dict['vsys'] #now all the velocities are based on the LSR

    # parse the channel values
    if chvals!=None:
        v1, v2 = chvals
    # parse the region, get the area for the spectra
    x1,x2,y1,y2 = parse_region(linedata, region)
    # if any is limits are None, draw spectra and ask for them/it
    # INTERACTIVE
    if nvals==None or chvals==None:
        # draw spectrum
        plot_spectrum(filename,region=region, obj_dict= {'vsys': obj_dict['vsys']})
        #
        # if there is no channels supplied
        if chvals==None:
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

    # calculate common stuff
    ### for the plotting box
    ylen, xlen = linedata.d[0].shape
    ycoords = arange(-ylen/2,ylen/2,1)*linedata.dec_cdelt
    xcoords = arange(-xlen/2,xlen/2,1)*linedata.ra_cdelt
    #
    # for the extent keyword
    left, right = xcoords[0],xcoords[-1]
    bottom, top = ycoords[0],ycoords[-1]
    #
    # parse the noise channel values
    if nvals!=None:
        # if noise has been input
        noise_channels = get_indices(linedata.velarr, nvals)
    #
    else:
        # if user didnt input noise velocity limits
        # choose channels away from the choosen (red blue) lines
        low = where((linedata.velarr>(linedata.velarr.min()+10))*(linedata.velarr<(v1-10)))[0]
        high = where((linedata.velarr>(v2+10))*(linedata.velarr<(linedata.velarr.max()-10)))[0]
        noise_channels = concatenate((low, high))
    #
    # the region to calculate the rms in
    i1,i2 = xlen/2+array([-1,1])*xlen/4
    j1,j2 = ylen/2+array([-1,1])*ylen/4

    # the noise, rms
    noise = linedata.d[noise_channels]
    rms = sqrt(((noise[:,j1:j2,i1:i2])**2).mean())

    # set boundaries
    if box == [0,0]:
        i1,i2 = left,right
        j1,j2 = bottom,top
    elif box != [0,0]:
        i1,i2 = array([-1,1])*box[0]/2.*sign(linedata.ra_cdelt)
        j1,j2 = array([-1,1])*box[1]/2.*sign(linedata.dec_cdelt)

    img_channels = get_indices(linedata.velarr,[v1,v2])
    imgs = linedata.d[img_channels]
    #
    # do the summing
    #
    # M1 = dv * sum(I(a,b,vi))_x1^x2
    # or perhaps without dv to save calc?
    moment0 = imgs.sum(axis=0)*abs(linedata.vcdeltkms)
    moment0_sigma = sqrt(alen(imgs))*rms*abs(linedata.vcdeltkms)

    # levels
    moment0_max = moment0.max()
    moment0_min = moment0.min()
    levels = concatenate((-1*flipud(arange(nsig*moment0_sigma,abs(moment0_min+moment0_sigma),nsjump*moment0_sigma)), arange(nsig*moment0_sigma,moment0_max+moment0_sigma,nsjump*moment0_sigma)))

    # calculate moment 1 map

    velocities = linedata.velarr[img_channels]
    # calculate the denominator
    Isum = imgs.sum(axis=0)
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
                'extent':(left,right,bottom,top),
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
        im = pl.imshow(moment1,cmap=cm.jet,vmin=m-2*w,vmax=m+2*w,extent=(left,right,bottom,top),interpolation='nearest')
    elif not type:
        im = ax.contourf(moment1, levels=arange(m-2*w,m+2*w,4*w/20), cmap=cm.jet, extent=(left,right,bottom,top))
    cs2 = ax.contour(moment0, levels=levels, colors='k', extent=(left,right,bottom,top))
    cbar = pl.colorbar(im,format=cbar_tick_label_formatter) #label_X.replace('X','%2.2f'))

    if str(linedata.unit) == 'Jy/beam':
        cbar.ax.set_ylabel(r'Jy\,beam$^{\sf-1}$')
    else:
        cbar.ax.set_ylabel(str(linedata.unit))



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
    ax.set_xlim(i1,i2)
    ax.set_ylim(j1,j2)
    ax.set_aspect(1)
    #
    #
    pl.subplots_adjust(left=pl_left, bottom=pl_bottom, right=pl_right, top=pl_top)




#
def plot_moment2 (filename,
                cfile=None,
                chvals=None,
                nvals=None,
                region=[0,0,0,0],
                nx=6,
                filled=True,
                box=[0,0],
                nsig=3,
                nsjump=2,
                obj_dict = dict(vsys=0),
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
    # we are using latex formatting here!
    # following sets the ticklabel format string
    data_fmt = '%g'         # for normal y and x axis
    cbar_data_fmt = '%2.2f'         # for normal y and x axis
    label_X = parse_tick_font(font)
    tick_label_formatter = FormatStrFormatter(label_X.replace('X',data_fmt))
    cbar_tick_label_formatter = FormatStrFormatter(label_X.replace('X',cbar_data_fmt))


    pl_left, pl_bottom, pl_right,pl_top = plot_adjust


    #
    # first get the fitsfile
    linedata = loadcube(filename,telescope)
    #
    linedata.velarr = linedata.velarr - obj_dict['vsys'] #now all the velocities are based on the LSR

    # parse the channel values
    if chvals!=None:
        v1, v2 = chvals
    # parse the region, get the area for the spectra
    x1,x2,y1,y2 = parse_region(linedata, region)
    # if any is limits are None, draw spectra and ask for them/it
    # INTERACTIVE
    if nvals==None or chvals==None:
        # draw spectrum
        plot_spectrum(filename,region=region, obj_dict= {'vsys': obj_dict['vsys']})
        #
        # if there is no channels supplied
        if chvals==None:
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

    # calculate common stuff
    ### for the plotting box
    ylen, xlen = linedata.d[0].shape
    ycoords = arange(-ylen/2,ylen/2,1)*linedata.dec_cdelt
    xcoords = arange(-xlen/2,xlen/2,1)*linedata.ra_cdelt
    #
    # for the extent keyword
    left, right = xcoords[0],xcoords[-1]
    bottom, top = ycoords[0],ycoords[-1]
    #
    # parse the noise channel values
    if nvals!=None:
        # if noise has been input
        noise_channels = get_indices(linedata.velarr, nvals)
    #
    else:
        # if user didnt input noise velocity limits
        # choose channels away from the choosen (red blue) lines
        low = where((linedata.velarr>(linedata.velarr.min()+10))*(linedata.velarr<(v1-10)))[0]
        high = where((linedata.velarr>(v2+10))*(linedata.velarr<(linedata.velarr.max()-10)))[0]
        noise_channels = concatenate((low, high))
    #
    # the region to calculate the rms in
    i1,i2 = xlen/2+array([-1,1])*xlen/4
    j1,j2 = ylen/2+array([-1,1])*ylen/4

    # the noise, rms
    noise = linedata.d[noise_channels]
    rms = sqrt(((noise[:,j1:j2,i1:i2])**2).mean())

    # set boundaries
    if box == [0,0]:
        i1,i2 = left,right
        j1,j2 = bottom,top
    elif box != [0,0]:
        i1,i2 = array([-1,1])*box[0]/2.*sign(linedata.ra_cdelt)
        j1,j2 = array([-1,1])*box[1]/2.*sign(linedata.dec_cdelt)

    img_channels = get_indices(linedata.velarr,[v1,v2])
    imgs = linedata.d[img_channels]
    #
    # do the summing
    #
    # M1 = dv * sum(I(a,b,vi))_x1^x2
    # or perhaps without dv to save calc?
    moment0 = imgs.sum(axis=0)*abs(linedata.vcdeltkms)
    moment0_sigma = sqrt(alen(imgs))*rms*abs(linedata.vcdeltkms)

    """
    Fit gaussian in the chvals interval, using velocity and cut out spectra
    in all pixels in region/box.
    Plot it!



    # levels
    moment0_max = moment0.max()
    moment0_min = moment0.min()
    levels = concatenate((-1*flipud(arange(nsig*moment0_sigma,abs(moment0_min+moment0_sigma),nsjump*moment0_sigma)), arange(nsig*moment0_sigma,moment0_max+moment0_sigma,nsjump*moment0_sigma)))

    # calculate moment 1 map

    velocities = linedata.velarr[img_channels]
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
            print colorify(errmsg)
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
                'extent':(left,right,bottom,top),
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
        im = pl.imshow(moment1,cmap=cm.jet,vmin=m-2*w,vmax=m+2*w,extent=(left,right,bottom,top),interpolation='nearest')
    elif not type:
        im = ax.contourf(moment1, levels=arange(m-2*w,m+2*w,4*w/20), cmap=cm.jet, extent=(left,right,bottom,top))
    cs2 = ax.contour(moment0, levels=levels, colors='k', extent=(left,right,bottom,top))
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
                obj_dict = dict(vsys=0),
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
    data_fmt = '%g'         # for normal y and x axis
    cbar_data_fmt = '%2.2f'         # for normal y and x axis
    label_X = parse_tick_font(font)
    tick_label_formatter = FormatStrFormatter(label_X.replace('X',data_fmt))
    cbar_tick_label_formatter = FormatStrFormatter(label_X.replace('X',cbar_data_fmt))
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
        linedata = loadcube(filename,telescope)
        #
        linedata.velarr = linedata.velarr - obj_dict['vsys'] #now all the velocities are based on the LSR
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
    
        #rms_0= rms/sqrt(abs(data.vcdeltkms))
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







def plot_pv (filename,
                cfile=None,
                chvals=None,
                nvals=None,
                region=[0,0,0,0],
                nx=6,
                filled=True,
                box=[0,0],
                nsig=3,
                nsjump=2,
                object=dict(vsys=0),
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
    linedata = loadcube(filename)

    ############## P-V diagram
    v1, v2 = chvals
    pl.ion()
    linedata.velarr = linedata.velarr - object['vsys'] #now all the velocities are based on the LSR
    #
    #sum along x-axis ie axis=2
    # parse the region, get the area for the spectra
    x1,x2,y1,y2 = parse_region(linedata, region)
    #x1,x2,y1,y2  = 245,265,245,265
    #axis=2 because I sum along x axis, ie naxis1
    d_arr = linedata.d[:,y1:y2,x1:x2].sum(axis=2)*linedata.vcdeltkms
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
    channels = get_indices(linedata.velarr,[v1,v2])
    im = ax1.contourf(linedata.velarr[channels], dec, d_arr.transpose()[:, channels], levs)
    cb = pl.colorbar(im)
    pl.xlim(v1,v2) # kms boundaries
    pl.ylim(dec.min(),dec.max())
    pl.xlabel(r'$v_{lsr}$')
    pl.ylabel(r'asec')

    fig2 = pl.figure(2)
    ax2 = pl.axes()
    im = ax2.imshow(linedata.d[channels,y1:y2,x1:x2].sum(axis=0)*linedata.vcdeltkms,
                    interpolation='nearest',cmap=cm.Greys)
    #im = ax2.pcolor(flipud(ra), dec, linedata.d[channels,y1:y2,x1:x2].sum(axis=0)*linedata.vcdeltkms)
    cb = pl.colorbar(im)
    ax2.axis('image')
    pl.xlabel(r'asec')
    pl.ylabel(r'asec')
#
def plot_chmap (filename,
                chvals=None,
                nvals=None,
                region=[0,0,0,0],
                nx=6,
                filled=True,
                box=[0,0],
                nsig=3,
                nsjump=2,
                object={'vsys':0},
                font={'family':'serif', 'serif': ['Times', 'Times New Roman'],
                'size':22, 'weight':'bold'},
                nsum=False,
                plot_adjust= [0.12, 0.09, 0.99, 0.99],
                quality=[125, 40],
                send=False,
                color_map='jet',
                cpeak = [0,0],
                locators = [1,0.2],
                levs_global=True):
    # imports
    #import scipy as sp
    print 'importing...'
    from scipy import array, where, median, std, sign, arange, alen, vstack, \
    concatenate, sqrt, log10, exp
    import pyfits as pf
    import matplotlib.pyplot as pl
    from matplotlib.pyplot import contour, contourf
    from mpl_toolkits.axes_grid import AxesGrid
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from matplotlib import cm, rc
    print 'done'
    #

    ################################
    # setting global rc properties #
    ################################
    rc('savefig', **{'dpi': quality[0]})
    rc('text', usetex=1)
    rc('figure',**{'facecolor': '1', 'dpi': quality[1]})
    #to set the global font properties
    rc('font', **font)
    rc('axes',linewidth=2)
    rc('lines', linewidth=1, markeredgewidth=1)
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
    # ticksize
    rc('xtick',**{'minor.size':5, 'major.size':7, 'major.pad': 9})
    rc('ytick',**{'minor.size':5, 'major.size':7, 'major.pad': 5})

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
    linedata = loadcube(filename)
    linedata.velarr = linedata.velarr - object['vsys'] #so that now all the velocities are based on the objects
    #if cfile!=None:
    #    contdata= loadcube(cfile)
    cfile=None
    #
    # parse the region parameter
    x1, x2, y1, y2 = parse_region(linedata, region)

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
        blue[x:x+y].sum(axis=0)*abs(data.vcdeltkms)
        calc_sigma(nchans,rms,vcdelt)


        # give chvals=[v1,v2] and nsum > 1
        # nsum is number of channels to sum together
        # get the channels
        """
        print '\nParameter \'nsum\' is set to %d.' % nsum
        print colorify('Summing channels for better signal.',c='g')
        # check if it is a integer
        if type(nsum) != type(1) or nsum<1: raise ParError(nsum)
        #
        channels = get_indices(linedata.velarr, chvals)
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
        velocity = linedata.velarr[channels]
        # the mean velocity is the middle one
        velocity = array([velocity[x:x+nsum].mean() for x in indices])
        velocity_delta = linedata.vcdeltkms*nsum
        print u'The velocity interval over which you create the maps is %d to %d km\u207b\u00b9\u00b7s' % (velocity.min(),velocity.max())
        maps = array([data[x:x+nsum].sum(axis=0)*abs(linedata.vcdeltkms) for x in indices])
        N_channels = alen(channels)/nsum
        print colorify('Done, remember that it can change the velocity interval that you specified in chvals. \n',c='g')
        #
    else :
        # normal procedure here
        #
        #we are making a channel map
        channels = get_indices(linedata.velarr, chvals)
        N_channels = alen(channels)
        velocity = linedata.velarr[channels]
        velocity_delta = linedata.vcdeltkms
        # units of Jy/beam kms
        data = linedata.d[channels]*abs(velocity_delta)
        maps = data

    #
    if nvals!=None:
        noise_channels = get_indices(linedata.velarr,nvals)
    elif nvals==None:
        # if user didnt input noise velocity limits
        # choose channels away from the choosen (red blue) lines
        #
        # if you have choosen just a min and max to plot all channels
        low = where((linedata.velarr>(linedata.velarr.min()+10))*(linedata.velarr<(v1-10)))[0]
        high = where((linedata.velarr>(v2+10))*(linedata.velarr<(linedata.velarr.max()-10)))[0]
        noise_channels = concatenate((low, high))

    # calculate common stuff
    ### for the plotting box
    zlen, ylen, xlen = data.shape
    ycoords = arange(-(linedata.dec_crpix-1),(linedata.dec_npix-linedata.dec_crpix),1)*linedata.dec_cdelt
    xcoords = arange(-(linedata.ra_crpix-1),(linedata.ra_npix-linedata.ra_crpix),1)*linedata.ra_cdelt

    # for the extent keyword
    left, right = xcoords[0],xcoords[-1]
    bottom, top = ycoords[0],ycoords[-1]
    #
    # calculate the RMS
    noise = linedata.d[noise_channels]
    #
    # the region to calculate the RMS in
    i1,i2 = xlen/2+array([-1,1])*xlen/4
    j1,j2 = ylen/2+array([-1,1])*ylen/4
    #
    # RMS
    rms = sqrt(((noise[:,j1:j2,i1:i2])**2).mean())
    #
    #
    if box == [0,0]:
        x1,x2 = left,right
        y1,y2 = bottom,top
    elif box != [0,0]:
        x1,x2 = array([-1,1])*box[0]/2.*sign(linedata.ra_cdelt)
        y1,y2 = array([-1,1])*box[1]/2.*sign(linedata.dec_cdelt)
    #
    #
    if nsum != False:
        sigma = calc_sigma(nsum, rms, linedata.vcdeltkms)
    else:
        sigma = calc_sigma(1,rms,linedata.vcdeltkms)
    #
    chans_min = array([i.min() for i in maps])
    chans_max = array([i.max() for i in maps])

    # calculate the levels of the contours
    if levs_global:
        f = 0.7
        vmin, vmax = f*chans_min.min(),f*chans_max.max()
        # levs_stat, i.e. contourf starts and end at global min/max
        levs_stat = concatenate((arange(vmin,-nsig*sigma,nsjump*sigma), arange(nsig*sigma,vmax,nsjump*sigma)))
        levs = [levs_stat for x in arange(N_channels)]
    elif not levs_global:
        # levs, i.e. contour starts and ends at local min/max
        levs = [concatenate((arange(x,-nsig*sigma,nsjump*sigma), arange(nsig*sigma,y,nsjump*sigma))) for x,y in zip(chans_min, chans_max)]

    ##########
    majorLocator = MultipleLocator(locators[0])
    minorLocator = MultipleLocator(locators[1])

    print u'RMS \t: %2.3f Jy\u00b7beam\u207b\u00b9\u00b7channel\u207b\u00b9' % rms
    print u'1 sigma : %2.3f Jy\u00b7beam\u207b\u00b9\u00b7km\u00b7s\u207b\u00b9' % sigma
    print u'\u0394 vel \t : %2.3f km\u00b7s\u207b\u00b9' % abs(velocity_delta)
    if send==True:
        if cfile!=None:
            return {'ch_map data': maps, 'levs':levs, 'chans_sig':sigma, 'vels':velocity, 'extent':[left,right,bottom,top], 'data':linedata, 'continuum':contdata}
        if cfile==None:
            return {'ch_map data': maps, 'levs':levs, 'chans_sig':sigma, 'vels':velocity, 'extent':[left,right,bottom,top], 'data':linedata}


    if velocity_delta<0:
        maps = flipud(maps)
        levs.reverse()
        velocity = flipud(velocity)
    #
    def print_vel(ax,x):
        ax.text(0.038,.885,
        str(round(velocity[x],2)),
        size='small',
        bbox=dict(edgecolor='k',facecolor='w',pad=10),
        transform = ax.transAxes)


    pl.ion()
    if nx == 5:
        fig = pl.figure(1, (14., 17.))
    else:
        fig = pl.figure(1, (14., 13.))

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
                                extent=(left,right,bottom,top))
    #
    elif filled==True:
        grid = AxesGrid(fig, 111, # similar to subplot(111)
                nrows_ncols = (ny, nx), # creates 2x3 grid of axes
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
                                extent=(left,right,bottom,top),
                                cmap=color_map)
        grid.cbar_axes[0].colorbar(im)
        # add a unit to the colorbar
        grid.cbar_axes[0].set_xlabel(str(linedata.unit)+r'\,kms$^{-1}$')
        grid.cbar_axes[0].axis["top"].toggle(label=True, ticks=True, ticklabels=True)
    #
    for i in range(N_channels):
        # plot a cross at pointing centre
        plt = grid[i].plot(cpeak[0],cpeak[1],'k+', ms=10)
        # set the locator spacings
        grid[i].xaxis.set_major_locator(majorLocator)
        #grid[i].xaxis.set_minor_locator(minorLocator)
        grid[i].yaxis.set_major_locator(majorLocator)
        #grid[i].yaxis.set_minor_locator(minorLocator)


        draw_fov(grid[i].axes, linedata)
        print_vel(grid[i].axes, i)
    draw_beam(grid[0].axes, linedata, box=1)
    #
    #grid.axes_llc.set_major_formatter(tick_label_formatter)
    #grid.axes_llc.set_major_formatter(tick_label_formatter)

    pl.subplots_adjust(left=pl_left, bottom=pl_bottom, right=pl_right, top=pl_top)


    grid.axes_llc.set_xlim(x1,x2)
    grid.axes_llc.set_ylim(y1,y2)

    fig.text(0.5,0.05,r'RA offset ($^{\prime\prime}$)', ha='center',va='center')
    fig.text(0.05,0.5,r'Dec offset ($^{\prime\prime}$)', rotation='vertical', va='center', ha='center')
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
        if data.vcdelt>0:
            bincrement = flipud(bincrement)
        bindices = arange(0,(len(blue)-brest),bincrement[0])

        # create the first 2 maps
        bmaps = array([blue[x:x+y].sum(axis=0)*abs(data.vcdeltkms) for x,y in zip(bindices,bincrement)])
        # calculate the center positions and their width
        bwidth = bincrement*abs(data.vcdeltkms)*.5
        bcenter = [data.velarr[x:x+y].mean() for x,y in zip(blue_channels[bindices],bincrement)]
        if data.vcdelt>0:
            bmaps = flipud(bmaps)
            bwidth = flipud(bwidth)
            bcenter = flipud(bcenter)
        # calculate the sigmas for each set of data
        bsigma = calc_sigma(bincrement, rms, data.vcdeltkms)

        ### red
        rrest = len(red)%3
        inc = (len(red)-rrest)/3

        # now flip the arrays (if vcdelt<0) so that the highest speeds swallow the
        # most channels, if  multiple of 3 no channels
        rincrement = array([inc, inc, inc+rrest])
        if data.vcdelt<0: # so that the last channel is still the one with the rest
            rincrement = flipud(rincrement)
        # the two middle indices
        # a bit tricky with the first being the biggest here
        rindices = array([0, rincrement[0], rincrement[0]+rincrement[1]])

        # get the maps, the last one (first if vcelt<0) takes the leftovers
        rmaps = array([red[x:x+y].sum(axis=0)*abs(data.vcdeltkms) for x,y in zip(rindices,rincrement)])

        #get the widht & center of each channelsum
        # flipud so that the channels are the same in blue and red (ie low med high velocity)
        rwidth = rincrement*abs(data.vcdeltkms)*.5
        rcenter = array([data.velarr[x:x+y].mean() for x,y in zip(red_channels[rindices],rincrement)])
        if data.vcdelt<0:
            rmaps = flipud(rmaps)
            rwidth = flipud(rwidth)
            rcenter = flipud(rcenter)

        rsigma = calc_sigma(rincrement,rms,data.vcdeltkms)

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



#########################################
# ASCII TABLE HELP FUNCTIONS
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
#
def loadatbl(filename, dtype='float64',sep=None):
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
                if line.startswith('#') or not line.strip():
                    continue # skip lines that are comments (# char) and empty
                cols = line.split(sep)
                values.append(array(cols,dtype=dtype))
    except IOError:
        raise IOError('file ' +str(filename)+' does NOT exist...')
    except ValueError:
        raise ValueError('Trying to convert to '+str(dtype)+' while it is a string\
                        try to change it to \'str\'')

    return array(values,dtype=dtype).transpose()
#
def infoatbl(filename, sep=None):
    """
    just returns the lines with comments (ie the column names) from a *.atbl file
    """
    try:
        with open(filename,'r') as f:
            strings = []
            for line in f:
                if line.startswith('#'):
                    strings.append(line.split(sep))
    except IOError:
        raise IOError('file' + str(filename)+'does not exist...')

    return strings
#

#########################################
# leftovers

def showmoment0 (filename, \
                cfile=None,\
                rvals=None,\
                bvals=None,\
                chvals=None,\
                nvals=None,\
                region=[0,0,0,0],\
                type='mom0',\
                nx=6,\
                filled=True,\
                box=[0,0],\
                nsig=3,\
                nsjump=2,\
                object={'vsys':0},\
                show_channels=False,\
                font=['serif','Times',12,'bold'],\
                bin=False,\
                quality=[125, 50],\
                fit = {'gauss': None, 'params': None, 'guess': False, 'interactive': False},\
                send=False):
    """
    def showmoment0(filename, \ # the line cube to use
                cfile=None, \ # continuum file (same size)
                rvals=None, \ # red velocity values
                bvals=None, \ # blue velocity values
                chvals=None, \ #
                nvals=None, \ # noise velocity limits
                region=[0,0,0,0], \ # miriad region type key
                type='mom0', \ # which plot to plot
                nx=6, \ #
                filled=True, \ #
                box=[0,0], \ # box to use for spectra etc
                nsig=3, \ # lowest contour level in sigmas
                nsjump=2, \ # distance between contour levels in sigmas
                font=['serif','Times',10,'bold'],\
                send=False,\ # send back the computed arrays/classes etc rather than plotting them
                ):

    usage: showmoment0(fitsfile, rvals=[r1,r2], bvals=[b1,b2])

    OR

    showmoment0(filename, cfile='cfilename') for completely interactive mode


    To change the box from where the spectra is plotted, change the
    region keyword, it is given as arcseconds in x and y (ra & dec)
    e.g region=[-2,-2,5,5] takes a box from x,y=-2,-2 to x,y=5,5 asecs
    centerd around the center of the map.

    type='mom0'
    means

    type='3range'
    means that it will create a channel map (3 maps for red
    and blue region respectively) instead of just the blueshifted
    and redshifted in one.

    type='chan'
    means it will create a channel map, i.e. one channel one image

    """

    # imports
    #import scipy as sp
    from scipy import array, where, median, std, sign, arange, alen, vstack, concatenate, sqrt, log10, exp
    import pyfits as pf
    import matplotlib.pyplot as pl
    from mpl_toolkits.axes_grid import AxesGrid
    #from matplotlib.patches import Circle
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from matplotlib import cm, rc
    #

    rc('font',**{'family':font[0],font[0]:[font[1]],'size':font[2], 'weight': font[3]})

    rc('axes',linewidth=2)
    rc('lines', linewidth=1)

    rc('savefig', **{'dpi': quality[0]})
    rc('text', usetex=True)
    rc('figure',**{'facecolor': '1', 'dpi': quality[1]})

    type1 = 'mom0'
    type2 = '3range'
    type3 = 'chan'
    type4 = 'spectra'

    #sigma levels etc
    nsig = nsig
    nsjump = nsjump
    #
    # first check the type1 and type2 plots
    if type in [type1, type2] and rvals==None and bvals==None:
        ion=True
    elif type in [type1, type2] and rvals!=None and bvals!=None:
        ion=False
        r1, r2 = rvals
        b1, b2 = bvals
    elif type in [type4] and chvals==None:
        ion=False
    # now check the type3 plots
    elif type in [type3] and chvals==None:
        ion=True
    elif type in [type3, type4] and chvals!=None:
        ion=False
        v1, v2 = chvals
    else: # input error
        print 'parse error: correct your input'
    #

    # first get the fitsfile
    data = loadcube(filename)
    data.velarr = data.velarr - object['vsys'] #so that now all the velocities are based on the objects
    if cfile!=None:
        continuum = loadcube(cfile)
    #


    if ion==True:
        pl.ion()
        x1,x2,y1,y2 = parse_region(data, region)
        if x1==x2 and y1==y2:
            spect = data.d[:,data.dec_npix/2,data.ra_npix/2]
        else:
            spect = data.d[:,y1:y2,x1:x2].sum(axis=1).sum(axis=1)
        # now correct so that negative values are to the left
        fig = pl.figure(1)
        fig.clf()
        ax = fig.add_subplot(111)
        if data.vcdelt<0:
            ax.step(flipud(data.velarr),flipud(spect),'k')
            xmin, xmax = round(flipud(data.velarr)[0]-1),round(flipud(data.velarr)[-1]+1)
            ymin, ymax = round(min(flipud(spect)*0.95)),round(max(flipud(spect)*1.05))
        elif data.vcdelt>=0:
            ax.step(data.velarr,spect,'k')
            xmin, xmax = round(data.velarr[0]-1),round(data.velarr[-1]+1)
            ymin, ymax = round(min(spect)*0.95),round(max(spect)*1.05)
        ax.plot([xmin,xmax],[0,0],'k:',[0,0],[ymin,ymax],'k:')
        ax.set_xlim(xmin,xmax); pl.ylim(ymin,ymax)
        ax.set_xlabel('km/s'); pl.ylabel('intensity')
        if type==type3:
            # if it is just a channel map
            v1, v2 = array(raw_input('input the limits, comma separated: ').split(','), dtype ='float')
        elif type==type1 or type==type2:
            # if mom0 or 3range
            b1, b2 = array(raw_input('input the blueshifted limits: ').split(','), dtype='float')
            r1, r2 = array(raw_input('input the redshifted limits, comma separated: ').split(','), dtype='float')
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
        pl.ioff()
    #
    if type==type1 or type==type2:
        # we are separating the red and blue channels
        red_channels = where((data.velarr>r1)*(data.velarr<r2))[0]
        blue_channels = where((data.velarr>b1)*(data.velarr<b2))[0]
        red = data.d[red_channels]
        blue = data.d[blue_channels]
    if type==type3:
        #we are making a channel map
        map_channels = where((data.velarr>v1)*(data.velarr<v2))[0]
        channels = data.d[map_channels]
    #
    # calculate common stuff
    ### for the plotting box
    ylen, xlen = data.d[0].shape
    ycoords = arange(-ylen/2,ylen/2,1)*data.dec_cdelt
    xcoords = arange(-xlen/2,xlen/2,1)*data.ra_cdelt

    # for the extent keyword
    left, right = xcoords[0],xcoords[-1]
    bottom, top = ycoords[0],ycoords[-1]

    if nvals!=None:
        # if noise has been input
        if len(nvals)==4:
            n1,n2,n3,n4 = nvals
            # if the user wants two velocity areas to calculate noise
            low = where((data.velarr>n1)*(data.velarr<n2))[0]
            high = where((data.velarr>n3)*(data.velarr<n4))[0]
            noise_channels = concatenate((low,high))
        elif len(nvals)==2:
            n1,n2 = nvals
            #  if input just want one velocity area to calculate noise
            noise_channels = where((data.velarr>n1)*(data.velarr<n2))

    elif nvals==None and type in [type1, type2, type3]:
        # if user didnt input noise velocity limits
        # choose channels away from the choosen (red blue) lines

        if type==type3:
            # if you have choosen just a min and max to plot all channels
            b2 = v2
            r1 = v1

        low = where((data.velarr>(data.velarr.min()+10))*(data.velarr<(r1-10)))[0]
        high = where((data.velarr>(b2+10))*(data.velarr<(data.velarr.max()-10)))[0]
        noise_channels = concatenate((low, high))
    elif nvals==None and type in [type4]:
        noise_channels = array([1,2])
    #

    # the region to calculate the rms in
    i1,i2 = xlen/2+array([-1,1])*xlen/4
    j1,j2 = ylen/2+array([-1,1])*ylen/4

    # the noise, rms
    noise = data.d[noise_channels]
    rms = sqrt(((noise[:,j1:j2,i1:i2])**2).mean())

    if box == [0,0]:
        x1,x2 = left,right
        y1,y2 = bottom,top
    elif box != [0,0]:
        x1,x2 = array([-1,1])*box[0]/2.*sign(data.ra_cdelt)
        y1,y2 = array([-1,1])*box[1]/2.*sign(data.dec_cdelt)
    #
    # MOM0 MAPS
    if type==type1:
        #mom0

        rc('axes',linewidth=2)
        rc('lines', linewidth=1)

        # one map for red and one for blue
        red = red.sum(axis=0)*abs(data.vcdeltkms)
        blue = blue.sum(axis=0)*abs(data.vcdeltkms)

        red_sigma = calc_sigma(len(red_channels),rms,data.vcdeltkms)
        blue_sigma = calc_sigma(len(blue_channels),rms,data.vcdeltkms)

        ############## the moment 0 maps

        # imstat
        #red_3sig = median(red)+2*std(red)
        #blue_3sig = median(blue)+2*std(blue)
        blue_min = blue.min(); blue_max = blue.max()
        #blue_levs = concatenate((arange(blue_min, -nsig*blue_sigma,nsjump*blue_sigma), arange(nsig*blue_sigma, blue_max,nsjump*blue_sigma)))
        blue_levs = arange(nsig*blue_sigma, blue_max,nsjump*blue_sigma)
        red_min = red.min(); red_max = red.max()
        #red_levs = concatenate((arange(red_min, -nsig*red_sigma,nsjump*red_sigma), arange(nsig*red_sigma, red_max,nsjump*red_sigma)))
        red_levs = arange(nsig*red_sigma, red_max,nsjump*red_sigma)
        if send==True:
            if cfile!=None:
                return {'blue': blue, 'blue_levs':blue_levs, 'blue_sigma':blue_sigma, 'red':red, 'red_levs':red_levs, 'red_sigma':red_sigma, 'extent':[left,right,bottom,top], 'data':data, 'continuum':continuum}
            if cfile==None:
                return {'blue': blue, 'blue_levs':blue_levs, 'blue_sigma':blue_sigma, 'red':red, 'red_levs':red_levs, 'red_sigma':red_sigma, 'extent':[left,right,bottom,top], 'data':data}
        pl.ion()
        fig = pl.figure(1)
        fig.clf()


        #ax = pwg.axes(header=data.hdr)d
        ax = fig.add_subplot(111)
        #ax.update_wcsgrid_params(label_density=(4,4))
        #blue_levs = arange(N*blue_sigma, blue_max,2*blue_sigma)
        ax.contour(blue,\
                    blue_levs,\
                    extent=(left,right,bottom,top),\
                    cmap=cm.winter)
        #red_levs = arange(N*red_sigma, red_max,2*red_sigma)
        ax.contour(red,\
                    red_levs,\
                    extent=(left,right,bottom,top), cmap=cm.hot)
        if cfile!=None:
            cont = ax.imshow(continuum.d,cmap=cm.Greys, extent=(left,right,bottom,top)) # remove [0] when fits files are OK
            #cbar = pl.colorbar(cont)
            #cbar.set_label(continuum.unit)
            draw_beam(ax, continuum,loc=4)
        draw_beam(ax,data)
        draw_fov(ax, data)
        draw_sizebar(ax, data,au=1000)
        #ax.plot([0],[0],'k*',ms=10)

        ax.set_xlabel('RA Offset ($^{\prime\prime}$)')
        ax.set_ylabel('DEC Offset ($^{\prime\prime}$)')
        #ax.set_title(data.obj+' (3sig red: ' +str(3*red_sigma.round(2))+'- 3sig blue: '+str(3*blue_sigma.round(2))+data.unit+' kms)')
        #print x1,x2,y1,y2
        ax.axis('image')
        ax.set_xlim(x1,x2)
        ax.set_ylim(y1,y2)
        pl.ioff()
        #pl.show()

    #
    #  3 RANGES
    elif type==type2:
        #
        # 3range
        #

        # 3 maps for red and blue each
        rc('axes',linewidth=1)
        rc('lines', linewidth=1)

        ### blue
        # if divided by 3, any rest?
        brest = len(blue)%3
        inc = (len(blue)-brest)/3
        bincrement = array([inc, inc, inc+brest])
        if data.vcdelt>0:
            bincrement = flipud(bincrement)
        bindices = arange(0,(len(blue)-brest),bincrement[0])

        # create the first 2 maps
        bmaps = array([blue[x:x+y].sum(axis=0)*abs(data.vcdeltkms) for x,y in zip(bindices,bincrement)])
        # calculate the center positions and their width
        bwidth = bincrement*abs(data.vcdeltkms)*.5
        bcenter = [data.velarr[x:x+y].mean() for x,y in zip(blue_channels[bindices],bincrement)]
        if data.vcdelt>0:
            bmaps = flipud(bmaps)
            bwidth = flipud(bwidth)
            bcenter = flipud(bcenter)
        # calculate the sigmas for each set of data
        bsigma = calc_sigma(bincrement, rms, data.vcdeltkms)

        ### red
        rrest = len(red)%3
        inc = (len(red)-rrest)/3

        # now flip the arrays (if vcdelt<0) so that the highest speeds swallow the
        # most channels, if  multiple of 3 no channels
        rincrement = array([inc, inc, inc+rrest])
        if data.vcdelt<0: # so that the last channel is still the one with the rest
            rincrement = flipud(rincrement)
        # the two middle indices
        # a bit tricky with the first being the biggest here
        rindices = array([0, rincrement[0], rincrement[0]+rincrement[1]])

        # get the maps, the last one (first if vcelt<0) takes the leftovers
        rmaps = array([red[x:x+y].sum(axis=0)*abs(data.vcdeltkms) for x,y in zip(rindices,rincrement)])

        #get the widht & center of each channelsum
        # flipud so that the channels are the same in blue and red (ie low med high velocity)
        rwidth = rincrement*abs(data.vcdeltkms)*.5
        rcenter = array([data.velarr[x:x+y].mean() for x,y in zip(red_channels[rindices],rincrement)])
        if data.vcdelt<0:
            rmaps = flipud(rmaps)
            rwidth = flipud(rwidth)
            rcenter = flipud(rcenter)

        rsigma = calc_sigma(rincrement,rms,data.vcdeltkms)

        ### put them together now
        centers = concatenate((bcenter,rcenter)).round(2)
        widths = concatenate((bwidth,rwidth)).round(2)
        maps_sigma = concatenate((bsigma,rsigma))

        maps = vstack((bmaps,rmaps))

        maps_max = array([i.max() for i in maps])
        maps_min = array([i.min() for i in maps])

        # now create the levels for each image


        levs = [concatenate((arange(x,-nsig*z,nsjump*z),arange(nsig*z,y,nsjump*z))) for x,y,z in zip(maps_min, maps_max, maps_sigma)]

        def print_velrange(ax,x):
            """
            pl.text(0.05,0.05,\
            str(round(rcenter[x-1],1))+'$\pm$'+str(round(rwidth[x-1],2))+' kms$^{-1}$',\
            bbox=dict(edgecolor='red',facecolor='red', alpha=0.7),\
            transform = ax.transAxes)

            """
            if x>2:
                fcolor='red'
                #ax.text(0.1,.9,\
                #str(round(centers[x],1))+'$\pm$'+str(round(widths[x],2))+' kms$^{-1}$',\
                #size='medium',\
                #bbox=dict(edgecolor='red',facecolor='red', alpha=0.7),\
                #transform = ax.transAxes)
            else:
                fcolor='blue'
                #ax.text(0.1,.9,\
                #str(round(centers[x],1))+'$\pm$'+str(round(widths[x],2))+' kms$^{-1}$',\
                #size='medium', weight='black',\
                #bbox=dict(edgecolor='blue',facecolor='blue', alpha=0.5),\
                #transform = ax.transAxes)
            ax.text(0.1,.9,\
            str(round(centers[x],1))+'$\pm$'+str(round(widths[x],2))+' kms$^{-1}$',\
            size='medium', weight='extra bold',\
            bbox=dict(edgecolor=fcolor,facecolor=fcolor, alpha=0.7),\
            transform = ax.transAxes)
        #

        # tick density and size
        # these are inserted later in:
        # grid[i].xaxis.set_major_locator(majorLocator)
        # grid[i].xaxis.set_minor_locator(minorLocator)

        majorLocator = MultipleLocator(10)
        minorLocator = MultipleLocator(2)

        from matplotlib import rc
        rc('xtick',**{'minor.size':2.5, 'major.size':5})
        rc('ytick',**{'minor.size':2.5, 'major.size':5})

        if send==True:
            if cfile!=None:
                return {'maps': maps, 'levs':levs, 'maps_sigma':maps_sigma, 'centers':centers, 'widths':widths, 'extent':(left,right,bottom,top), 'data':data, 'continuum':continuum}
            if cfile==None:
                return {'maps': maps, 'levs':levs, 'maps_sigma':maps_sigma, 'centers':centers, 'widths':widths, 'extent':(left,right,bottom,top), 'data':data}

        pl.ion()
        fig = pl.figure(1, (9., 7.))
        fig.clf()


        if cfile!=None:
            # if cfile entered, add cbars
            grid = AxesGrid(fig, 111, # similar to subplot(111)
                            nrows_ncols = (2, 3), # creates 2x3 grid of axes
                            axes_pad=0.2, # pad between axes in inch.
                            label_mode = "L",
                            share_all=True,
                            cbar_mode="single",
                            cbar_location="top",
                            cbar_size="2%",
                            )
        else:
            # if no cfile, no colorbar
            grid = AxesGrid(fig, 111, # similar to subplot(111)
                            nrows_ncols = (2, 3), # creates 2x3 grid of axes
                            axes_pad=0.2, # pad between axes in inch.
                            label_mode = "L",
                            share_all=True,
                            )

        for i in range(6):

            col = ['b','b','b','r','r','r']
            im = grid[i].contour(maps[i], levs[i], colors=col[i], extent=(left,right,bottom,top))
            if cfile!=None:
                # if there is a continuum file
                cont = grid[i].imshow(continuum.d,cmap=cm.Greys, extent=(left,right,bottom,top)) # remove [0] when fits files are OK
                grid.cbar_axes[i].colorbar(cont)
                # add a unit to the colorbar
                grid.cbar_axes[i].set_xlabel(str(continuum.unit))
                grid.cbar_axes[i].axis["top"].toggle(label=True, ticks=True, ticklabels=True)
                draw_beam(grid[i].axes, continuum,loc=4)
                # end
            plt = grid[i].plot([0],[0],'k+', ms=10)
            grid[i].xaxis.set_major_locator(majorLocator)
            grid[i].xaxis.set_minor_locator(minorLocator)
            grid[i].yaxis.set_major_locator(majorLocator)
            grid[i].yaxis.set_minor_locator(minorLocator)
        for i in range(6):
            draw_beam(grid[i].axes, data)
            print_velrange(grid[i],i)
            draw_fov(grid[i].axes, data)
        grid.axes_llc.set_xlim(x1,x2)
        grid.axes_llc.set_ylim(y1,y2)
        #grid.axes_llc.set_aspect(1)
        #grid.axes_llc.set_axes('image')
        grid[4].set_xlabel('RA offset')
        #grid[4].text(-0.5,-0.1,'string')
        #grid.cbar_axes[0].axis["top"].set_title(data.obj)
        #fig.text(0.51,0.03,'RA offset (asec)', ha='center',va='center')
        fig.text(0.05,0.5,'DEC offset (asec)', rotation='vertical', va='center', ha='center')
        #fig.suptitle(data.obj)
        #pl.show()
        pl.ioff()

    #
    # CHANNEL MAP
    elif type==type3:
        # map_channels - the index of the channels used
        # channels - the data
        channels = channels*abs(data.vcdeltkms)
        vels = data.velarr[map_channels]

        # if we are plotting channel maps
        chans_sig = calc_sigma(1,rms,data.vcdeltkms)
        chans_min = array([i.min() for i in channels])
        chans_max = array([i.max() for i in channels])

        # calculate the levels of the contours
        levs = [concatenate((arange(x,-nsig*chans_sig,nsjump*chans_sig), arange(nsig*chans_sig,y,nsjump*chans_sig))) for x,y in zip(chans_min, chans_max)]
        vmin, vmax = chans_min.min(),chans_max.max()
        levs_stat = concatenate((arange(vmin,-nsig*chans_sig,nsjump*chans_sig), arange(nsig*chans_sig,vmax,nsjump*chans_sig)))

        majorLocator = MultipleLocator(1)
        minorLocator = MultipleLocator(0.2)

        from matplotlib import rc
        rc('xtick',**{'minor.size':0, 'major.size':2.5})
        rc('ytick',**{'minor.size':0, 'major.size':2.5})
        rc('axes',linewidth=1)
        rc('lines', linewidth=1)

        if send==True:
            if cfile!=None:
                return {'channels': channels, 'levs':levs, 'chans_sig':chans_sig, 'vels':vels, 'extent':[left,right,bottom,top], 'data':data, 'continuum':continuum}
            if cfile==None:
                return {'channels': channels, 'levs':levs, 'chans_sig':chans_sig, 'vels':vels, 'extent':[left,right,bottom,top], 'data':data}


        if data.vcdelt<0:
            channels = flipud(channels)
            levs.reverse()
            vels = flipud(vels)
        def print_vel(ax,x):
            ax.text(0.08,.9,\
            str(round(vels[x],1))+r'kms$^{-1}$',\
            size='small',\
            bbox=dict(edgecolor='white',facecolor='white', alpha=0.7),\
            transform = ax.transAxes)


        pl.ion()
        fig = pl.figure(1, (14., 9.))

        fig.clf()

        ny = int(ceil(len(map_channels)/float(nx)))

        if filled==False:
            grid = AxesGrid(fig, 111, # similar to subplot(111)
                    nrows_ncols = (ny, nx), # creates nyx6 grid of axes
                    axes_pad=0.2, # pad between axes in inch.
                    label_mode = "L",
                    share_all=True,
                    )
                # plot data contours

            for i in range(len(map_channels)):
                grid[i].set_aspect('equal')
                im = grid[i].contour(channels[i],\
                                    levs[i],\
                                    cmap=cm.Greens,\
                                    extent=(left,right,bottom,top))
        #

        elif filled==True:
            grid = AxesGrid(fig, 111, # similar to subplot(111)
                    nrows_ncols = (ny, nx), # creates 2x3 grid of axes
                    axes_pad=0, # pad between axes in inch.
                    label_mode = "L",
                    share_all=True,
                    cbar_mode="single",
                    cbar_location="top",
                    cbar_size="2%",)
            # plot data contours
            for i in range(len(map_channels)):
                grid[i].set_aspect('equal')
                im = grid[i].contourf(channels[i],\
                                    levs_stat,\
                                    extent=(left,right,bottom,top))
            grid.cbar_axes[0].colorbar(im)
            # add a unit to the colorbar
            grid.cbar_axes[0].set_xlabel(str(data.unit)+r' kms$^{-1}$')
            grid.cbar_axes[0].axis["top"].toggle(label=True, ticks=True, ticklabels=True)
        #

        for i in range(len(map_channels)):
            # plot a cross at pointing centre
            plt = grid[i].plot([0],[0],'k+', ms=5)
            # set the locator spacings
            grid[i].xaxis.set_major_locator(majorLocator)
            grid[i].xaxis.set_minor_locator(minorLocator)
            grid[i].yaxis.set_major_locator(majorLocator)
            grid[i].yaxis.set_minor_locator(minorLocator)
            # beam, fov and velocity info
            draw_beam(grid[i].axes, data, box=False)
            draw_fov(grid[i].axes, data)
            print_vel(grid[i].axes, i)




        grid.axes_llc.set_xlim(x1,x2)
        grid.axes_llc.set_ylim(y1,y2)

        fig.text(0.5,0.05,'RA offset (asec)', ha='center',va='center')
        fig.text(0.05,0.5,'DEC offset (asec)', rotation='vertical', va='center', ha='center')
        fig.suptitle(data.obj)
        #pl.show()
        pl.ioff()
    #
    # SPECTRA
    elif type==type4:
        # spectra
        import matplotlib.transforms as mtransforms
        from mpl_toolkits.axes_grid.parasite_axes import SubplotHost
        #
        x1,x2,y1,y2 = parse_region(data, region)
        pl.ion()
        area_region = ((y2-y1)*(x2-x1))
        spect = data.d[:,y1:y2,x1:x2].sum(axis=1).sum(axis=1)/area_region
        # shift the velocity array so that 0 is at system velocity (if supplied)
        velocity = data.velarr


        # Binning.
        # So
        #
        if send==True:
            return spect, velocity
        if bin!=False:
            j = int(bin)
            indices = arange(0,alen(spect),j)
            spect = array([spect[x:x+j].sum(axis=0)/j for x in indices])
            velocity = array([velocity[x:x+j].sum(axis=0)/j for x in indices])
            # use congridding???
        #
        # set ticklocations
        xmajorLocator = MultipleLocator(5)
        xminorLocator = MultipleLocator(1)
        ymajorLocator = MultipleLocator(0.02)
        yminorLocator = MultipleLocator(0.01)
        #xTopMajorLocator = MultipleLocator(5)
        xTopMajorFormatter = FormatStrFormatter('%5.6f')

        rc('axes',linewidth=2)
        rc('lines', linewidth=1)

        fig = pl.figure(1,figsize=(10.5,8))
        fig.clf()
        ax_kms = fig.add_subplot(111)
        #ax_kms = SubplotHost(fig, 1,1,1)
        # kms_to_hz - kms to hz convertion
        #kms_to_hz =
        # now correct so that negative values are to the left
        if data.vcdelt<0:
            ax_kms.step(flipud(velocity), flipud(spect), 'k', lw=2,  where='mid')
            xmin, xmax = round(flipud(velocity)[0]-data.vcdeltkms), round(flipud(velocity)[-1]+data.vcdeltkms)
            ymin, ymax = round(min(flipud(spect)),2)*1.01, round(max(flipud(spect)),2)*1.01
        elif data.vcdelt>=0:
            ax_kms.step(velocity, spect, 'k', lw=2, where='mid')
            xmin, xmax = round(velocity[0]-1),round(velocity[-1]+1)
            ymin, ymax = round(min(spect),2)*1.01,round(max(spect),2)*1.01
        if nvals!=None:
            #rms = sqrt(((data.d[noise_channels,y1:y2,x1:x2])**2).mean())
            rms = sqrt(((spect[noise_channels])**2).mean())
            #sig = calc_sigma(1/0.12,rms,0.12)
            ax_kms.plot([xmin, xmax],[rms,rms],'b')
        #ax.xaxis.set_major_locator(xmajorLocator)
        #ax.xaxis.set_minor_locator(xminorLocator)
        #ax.yaxis.set_major_locator(ymajorLocator)
        #ax.yaxis.set_minor_locator(yminorLocator)
        ax_kms.plot([xmin,xmax],[0,0],'k:',[0,0],[ymin,ymax],'k:')
        if fit['gauss']=='1d':
            # fits up to five (5) 1D gaussian(s) to the spectral data
            #
            #
            #
            if chvals!=None:
                ch = where((velocity>v1)*(velocity<v2))[0]
                Y = spect[ch]
                X = velocity[ch]
            else:
                Y = spect
                X = velocity
            p = fit['params'] # fitting 2 gaussians
            param, errors, success, no_fits = gaussfit1d((X,Y), params=p)
            ### test
            #fitfunc = lambda p, x: p[0]*cos(2*pi/p[1]*x+p[2]) + p[3]*x # Target function
            #errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
            #p0 = [-15., 0.8, 0., -1.] # Initial guess for the parameters
            #p1, success = optimize.leastsq(errfunc, p0[:], args=(Tx, tX))
            ### test
            ax_kms.plot(X,gaussian(param,X),'r-', lw=2)
        #
        ax_kms.set_xlim(xmin,xmax)
        ax_kms.set_ylim(ymin,ymax)
        #ax.set_title(data.obj)
        ax_kms.set_xlabel('km/s')
        ax_kms.set_ylabel('Intensity ('+data.unit+')')
        # to show the channel numbers
        if show_channels==True:
            ax_hz = ax_kms.twiny()
            ax_hz.xaxis.set_major_formatter(xTopMajorFormatter)
            def frequency (v):
                return (1-v/299792.458)*203.407520000
            #def update_ax_hz(ax_kms):
            #   x_1, x_2 = ax_kms.get_xlim()
            #   ax_hz.set_xlim(frequency(x_1), frequency(x_2))
            #   ax_hz.figure.canvas.draw()

            #ax_hz.xaxis.set_major_locator(xTopMajorLocator = MultipleLocator(5))
            #raise Exception ("\n\n\n This function does not work properly yet\n The channels are not correctly marked\n\n")

            #ax_kms.callbacks.connect("xlim_changed", update_ax_hz)

            from scipy import NaN
            #print data.restfreq

            x_1, x_2 = ax_kms.get_xlim()
            ax_hz.set_xlim(frequency(x_1), frequency(x_2))
            #rax_hz.figure.canvas.draw()

                #ax_hz.plot(flipud(data.freqarr/1e9),arange(alen(data.freqarr))*0)
            #elif data.vcdelt>=0:
            #    ax_hz.set_xlim(data.freqarr[0]/1e9, data.freqarr[-1]/1e9)
            #    ax_hz.figure.canvas.draw()
                #ax_hz.plot(data.freqarr,arange(alen(data.freqarr))*0)
        pl.ioff()
#
def showchannelmap(filename,v_sys):
    ############## moment maps
    """
    make it divide it in to three components
    low, medium and high velocity around
    the systemvelocity

    """
    #vel_maps = [] # empty list to store the maps
    # vcdeltkms is one channel
    # number of channels (vcdelt) to sum for each map
    def maps(v1,v2):


        #increment = Nincr*vcdeltkms # in kms, not m
        # the velocity boundaries have to be an even number
        # pair of velocities determining the boundaries
        #v1=-11; v2=0 # red part
        #v1=9; v2=24 # blue part
        #list_vel = [-11,0,9,24] # pair of velocities determining the  boundaries

        a1 = where((velarr-.49).round()==v1)[0][0]
        a2 = where((velarr+.49).round()==v2)[0][0]
        if a1>a2:
            lower = a2
            upper = a1
        elif a2>a1:
            lower = a1
            upper = a2


        vels = velarr[lower:upper+1]
        if Nincr>1:
            vel_maps = array(d[lower:upper+1]).reshape((upper-lower+1)/Nincr,Nincr,y_npix,x_npix).sum(axis=1)/(vcdeltkms*Nincr)
        elif Nincr==1:
            vel_maps = array(d[lower:upper+1])
        # calculate the center velocity
        vels = velarr[lower:upper+1:Nincr]+Nincr/2*vcdeltkms
        width = Nincr*vcdeltkms
        return vel_maps, vels, width

    #v1=-11; v2=0 # red part
    #v1=9; v2=24 # blue part
    vsys=8.3
    a1,b1,c1 = test(vsys,vsys+14)
    a2,b2,c2 = test(vsys-14,vsys)
    b1 = b1.round(2)-vsys # round of the velocities
    b2 = b2.round(2)-vsys
    levs = arange(.25,2.5,0.1)
    pl.figure(1,figsize=(10,5))
    j=1
    for i in arange(1,len(a1)+1):
        pl.subplot(3,4,j)
        pl.contour(a1[i-1],levs)
        if i in [2,3,4,6,7,8,10,11,12]:
            loc,lab = pl.yticks()
            pl.yticks(loc,[])
        if i in arange(1,9):
            loc,lab = pl.xticks()
            pl.xticks(loc,[])
        pl.text(x1+2,y2-8,str(b1[i-1])+'km/s', bbox=dict(edgecolor='white',facecolor='blue', alpha=0.7))
        pl.axis('image')
        pl.xlim(x1,x2)
        pl.ylim(y1,y2)
        j+=1

    j=len(a1)+1
    for i in arange(0,len(a2)):
        pl.subplot(3,4,j)
        pl.contour(a2[i],levs)
        if j in [2,3,4,6,7,8,10,11,12]:
            loc,lab = pl.yticks()
            pl.yticks(loc,[])
        if j in arange(1,5):
            loc,lab = pl.xticks()
            pl.xticks(loc,[])
        pl.text(x1+2,y2-8,str(b2[i])+'km/s', bbox=dict(edgecolor='white',facecolor='red', alpha=0.7))
        pl.axis('image')
        pl.xlim(x1,x2)
        pl.ylim(y1,y2)
        j+=1

    pl.subplots_adjust(left=0.1, bottom=0.1, right=0.97, top=0.97, wspace=0, hspace=0)

#
# 
if __name__ == '__main__':
    # statements that you want to be executed only when the
    # module is executed from the command line
    # (not when importing the code by an import statement)
    # wee hooo...
    
    print '\nThe ADAVIS program was run from the command line\n'
    import os
    import sys
    from optparse import OptionParser as op
    import pyfits as pf
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
    if options.v != None:
        v = convert_to_list(options.v)
    else:
        v = options.v
    if options.n != None:
        n
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
    #
    a = str(raw_input('Press any key to exit...'))
    # put Jes, optparser here!

