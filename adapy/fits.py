
from .helpers import *



########################################################################
# USEFUL STRINGS
KMS = u"km\u00b7s\u207b\u00b9"


########################################################################
# DATA HANDLING CLASSES
# to adacore.py
# classes etc
#
# FITS DATA CLASS
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
    def __init__(self, fitsfile, telescope=None, vsys=0, distance=0, **kwargs):
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
        f = fitsopen(fitsfile, **kwargs)
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
            self.restfreq = Unit(self.hdr['RESTFREQ'],'Hz' ) # in Hertz
        except KeyError:
            try:
                self.restfreq = Unit(self.hdr['RESTFRQ'],'Hz' ) # in Hertz
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
            # _load_IMAGE_data()
            self.datatype = ('IMAGE',2)
            # extra if the continuum image has the freq and width
            self.freq = self.hdr['CRVAL3']
            self.freqwidth = self.hdr['CDELT3']
        # a spectra! the 3rd axis is longer than 1
        elif self.hdr['NAXIS']>=3 and self.hdr['NAXIS3']>1:
            # spectral cube
            # _load_CUBE_data()
            # only support for velo-lsr in 3rd axis
            self.datatype = ('CUBE',3)
            # load the third axis
            # need frequency!
            while self.d.shape[0] == 1:
                self.d = self.d[0]
            ##### have to add loading of frequency and calculate velocity
            # UGLY HACK BELOW, BEWARE!
            # need to be changed to a more flexible code...
            hdr_values = [self.hdr[i] for i in self.hdr.keys() if not i == 'HISTORY']
            #~ self.hdr_values = hdr_values
            #~ return None
            # need to match VELO, VELO-LSR, VELOCITY and VRAD
            _velname = [i for i in hdr_values if ("VELO" in str(i) or "VRAD" in str(i))]
            if _velname != []:
                _velname = _velname[0]
                velax = str([x for x in self.hdr.keys() if x[:-1]=='CTYPE' and _velname in self.hdr[x]][0][-1:])
                vel_info = True
            elif _velname == []:
                print('No velocity axis defined')
                vel_info = False
            else:
                print('No velocity axis defined')
                vel_info = False
            ##### FUGLY hack END
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
                # END ugly hack
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
            #~ already on line 1480
            #~ if self.hdr.has_key('RESTFREQ'):
                #~ self.restfreq = self.hdr['RESTFREQ']
            #~ elif self.hdr.has_key('RESTFRQ'):
                #~ self.restfreq = self.hdr['RESTFRQ']
            #if self.hdr['NAXIS']==4 and  self.hdr['NAXIS4']==1:
            #    self.d = self.d[0]
            #
            # this was just a test
            # SD pointing spectra
            # single dish
        elif self.hdr['NAXIS']>1 and self.hdr['NAXIS2']==1 and self.hdr['NAXIS3']==1:
            # _load_SD_data(self)
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
        #TODO make npix dynamically determined, so if we cut in the image
        # it updates it, and the crpix
        self.dec_npix = self.hdr['NAXIS'+decax] 
        self.y_npix = self.hdr['NAXIS'+decax]
        #TODO crpix has to be updated when cutting in the image!
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
            #TODO make the extent keyword update dynamically with how
            # many pixels that are there...
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
            self.bmaj = Unit(self.hdr['BMAJ']*3600, 'asecs')
            self.bmin = Unit(self.hdr['BMIN']*3600, 'asecs')
            self.bpa = Unit(self.hdr['BPA'], 'degrees?')
        except KeyError, ex:
            msg='Header keywords (bmaj,bmin,bpa) incomplete and not loaded.'
            print(msg)
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
            print('No beam unit in header.')
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
        TODO : move this to the spectrum class?
        '"""
        try:
            known_lines = self.known_lines
        except AttributeError, ex:
            known_lines = {}
        known_lines[204.38343] = {'name' : 'SO$_2$','frequency' : frequency, 'channels' : channels, 'width' : width}
    #
    def phase_center_string(self):
        ra = parse_ra(self.ra_crval, string = True)
        dec = parse_dec(self.dec_crval, string = True)
        center = [ra, dec]
        return center
    
    def calc_offset(self, InputData):
        #~ ra, dec = self.phase_center_string()
        ra, dec = self.ra_crval, self.dec_crval
        return calc_offset(ra, dec, data = InputData, display = False)
    # method to change the v_sys
    def change_v_sys (self, v_sys):
        self.v_sys = v_sys
        # now, change the v_arr_syscorr array as well
        if self.datatype[0] in ['CUBE', 'SDSPECT']:
            self.v_arr_syscorr = self.v_arr - self.v_sys
    #
    def change_dist (self, dist):
        self.dist = dist # unit of pc

    def box_cut(self,region=[-10,10,-10,10]):
        pass
    

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
        import adapy
        from adapy.libs import cgsconst
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
        self.u_m = self.hdu.data.par(0) * cgsconst.CC*1e-2
        self.v_m = self.hdu.data.par(1) * cgsconst.CC*1e-2
        self.w_m = self.hdu.data.par(2) * cgsconst.CC*1e-2
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



class LoadFits:
    def __init__():
        pass



    def __str__():
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

    # Functions to define things
    def define_vlsr():
        pass

    def define_dpc():
        pass

    # Functions to confine data
    def extract_box():
        pass

    def extract_interval():
        pass

    # Functions to parse coordinates
    def parse_spatial(self,coord, ctype='offset', shape='box'):
        """
            coord       : either [Xi,Xj,Yi,Yj] or
                          [Xi, Yi, R]

            shape       : 'box' or 'circle'
            
            ctype       : 'pixel' - actual pixels
                          'offset' - i.e. arcsec offset
        """
        pass

    def parse_spectral(self, ctype='velocity', otype='where'):
        """
        Parse coordinates along the spectral axis
        in velocity or frequency to indices
            ctype       : coordinate type
                          'velocity', 'vel', 'v'
                          or 'frequency', 'freq', 'f'
            otype       : output type
                          'where' - e.g. array([3,4,5,6,7,8,9])
                          'interval' - e.g. [3,10]
            
        """
        pass
    
    def save_fits():
        pass





