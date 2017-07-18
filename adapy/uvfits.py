



import astropy.units as u
import astropy.constants as co
from astropy.io.fits import open as pfopen
from astropy import wcs
from .uvfits_helpers import *

# check the imports below
from datetime import datetime as _dt
from .helpers import *
from .libs.date import jd2gd
import scipy as _sp
from scipy import sqrt, pi, arctan2
import numpy as _np
import numpy as np
#import adapy




########################################################################
# STRINGS
KMS = u"km\u00b7s\u207b\u00b9"

__all__ = ['Uvfits']

########################################################################
# DATA HANDLING
# UV-FITS DATA CLASS
class Uvfits(object):
    """
    Read uv-fits data
    --------------------------------------------------------------------
    Normal structure of UV-fits data:

    Header : same as for normal fits

    Data : Group data

    --------------------------------------------------------------------

    TODO :  Assumes that CRVAL4 is frequency, is that always true?
            Make it more robust, look for the frequency keyword, either as
            "restfreq" or as a "crvalX"

    TODO :  UV fit method
    TODO : __str__ method, information about:
                - Phase center
                - No correlations
                - No baselines
                - No antennas
                - Telescope (if present)
    """
    def __init__(self, uvfitsfile, telescope=None, vsys=0, distance=0, endian=None, **kwargs):
        """

        Reads the uvfits and calculates useful things, e.g. u,v,w,
        phase and amplitude

        .byteswap().newbyteorder() is applied in various places to
        convert to little endian

        """
        f = pfopen(uvfitsfile, **kwargs)
        self.loadendian = endian
        if f[0].header['NAXIS1'] != 0:
            print "error: this file may not be a UV FITS."
            raise FileError('File format error.')
        #~ f.info()
        try:
            self.hdu = f[0]
        except:
            print "error: cannot open uv data HDU."
        self.hdr = self.hdu.header
        self.data = self.hdu.data
        if self.hdr['NAXIS4'] > 1:
            self.datatype = ('CUBE', 3)
        else:
            self.datatype = ('IMAGE', 2)
        
        # find spectral axis
        axis_types = self.WCS.get_axis_types()
        ax_types = np.array([i['coordinate_type'] for i in axis_types])
        try:
            spec_axis = ('spectral' == ax_types).nonzero()[0][0]
            freq = self.hdu.header['CRVAL{0}'.format(spec_axis+1)]
            # assumes the frequency given in Hz
            self.freq = freq * un.Hz
        except (IndexError):
            print('No spectral axis in header.')
            spec_axis = -1
            self.freq = None
        
        
        if 'RESTFREQ' in self.hdu.header.keys():
            self.restfreq = self.hdu.header['RESTFREQ']
            self.restfreq_unit = self.hdu.header['RESTFREQ'] * u.Hz
        else:
            raise StandardError('No restfrequency found, NEED it!')
        #TODO : Read in velocity and frequency array if present
        """
        The standard unit is to give UU and VV in seconds (??!?)
        So we have to convert to whatever we want.
        """
        # standard storing unit here is kilo-lambdas
        # save a million lines of code!
        u.add_enabled_equivalencies(lambdas_equivalencies(self.restfreq_unit))
        self.u = (self.data.par('UU') * u.s).to(klambdas)
        self.v = (self.data.par('VV') * u.s).to(klambdas)
        self.w = (self.data.par('WW') * u.s).to(klambdas)
        self.uvdist = sqrt(self.u.value**2 + self.v.value**2) * klambdas


        # BASELINE
        self.baseline = self.hdu.data.par('BASELINE').byteswap().newbyteorder()
        # DATES
        self.jdate = self.hdu.data.par('DATE')
        # set date to 1900, Jan, 01, 01:00:00 if date before before this
        self.jdate =self.jdate.clip(2415020.5)
        self.date = _sp.array([jd2gd(i) for i in self.jdate])
        self.date0 = self.date.transpose()
        fields = ['year', 'month', 'day', 'hour', 'minute', 'sec']
        self.date1 = {key:value for key,value in zip(fields, self.date0)}
        # convert to datetime objects
        # LOSES the sub-second resolution
        self.date2 = [_dt(int(i[0]), int(i[1]), int(i[2]), int(i[3]), int(i[4]), int(i[5])) for i in self.date]
        # get number of tracks
        # TODO : rough hack, separate track if diff day is >1
        tmp = _sp.where(_sp.diff(_sp.unique(self.jdate.round(0)))>1)[0]
        self.ntracks = len(tmp)+1
        
        ################################################################
        # NB : need to streamline this.
        # only load the complex visibilities, into a complex array
        # AND then work on that
        
        
        # COMPLEX VISIBILITY
        visi_index = len(self.data.parnames)
        if self.hdu.header['NAXIS']  == 7:
            self.visdata = self.data.par(visi_index)[:,0,0,0,0,0,:].byteswap().newbyteorder()
        #~ self.visdata = self.hdu.data.data[:,0,0,0,0,0,:]
        elif self.hdu.header['NAXIS']  == 6:
            self.visdata = self.data.par(visi_index)[:,0,0,0,0,:].byteswap().newbyteorder()
        # load the re, im and weight arrays
        self.re, self.im, self.wt = self.visdata[:,:].T
        #~ self.re = self.visdata[:,0][:]
        #~ self.im = self.visdata[:,1][:]
        #~ self.wt = self.visdata[:,2][:]
        # complex numbers
        #~ self.comp = self.visdata[:,:2].astype(_np.float64).view(_np.complexfloating)
        #~ self.comp = 1j*self.visdata[:,1][:]
        #~ self.comp += self.visdata[:,0][:]
        #~ self.comp = self.visdata[:,:2].astype(_np.float).view(_np.complex)
        
        # below seems a bit dependent...
        self.cvisi = self.visdata[:,:2].astype(_np.float).view(_np.complex).T[0]
        """
        with complex array, you can do
        amp = np.abs(vis)
        np.angle(vis)   
        vis.real
        vis.imag
        
        """
        # the data is not shifted
        self.isshifted = (False, [0,0])
        # AMPLITUDE 
        self.amp = sqrt(self.re**2 + self.im**2)
        # PHASE
        self.pha = arctan2(self.im, self.re)
        self.pha_deg = self.pha / pi * 180.
        # ERROR / SIGMA
        #TODO : check
        # following 1.0e6 is just for GILDAS, change if needed
        #~ print('NB : Error calculated from weights assuming GILDAS '
        #~ 'data (i.e. frequencies in MHz).')
        self.sigma_alt = 1/sqrt(self.wt*1.0e6)
        # Daniels way of calculating sigma
        # test this first
        self.sigma = _sp.sqrt(0.5 / ( self.wt * float(self.amp.shape[0]) ) )
        #np.sqrt( 0.5/self.wt/float(self.amp.shape[0]) )

    def load_model(self, modelfile, endian = None):
        if endian != None: # overrides the endianness of data loading
            self.Model = Uvfits(modelfile, endian = endian)
        else:
            # make sure it loads the same endian format as data
            self.Model = Uvfits(modelfile, endian = self.loadendian)

    def bin_data_DMC(self,ruv=None, binsize=10):
        """
        Function to bin data
        
        """
        class BinnedDMC:
            pass

        # Bin model
        if self.__dict__.has_key('Model'):
            # If model is loaded, bin that as well
            # to the same bins
            pass
        if ruv is not None:
            uvdist = ruv
        else:
            uvdist = self.uvdist_klam
    
        uvmin = uvdist.min()
        uvmax = uvdist.max()
        # Define the bins
        nbin = int( (uvmax-uvmin)/binsize)+5
        arr_bins = _sp.arange(nbin)
        arr_bins = binsize * arr_bins
        arr_bins1 = 0.5*(arr_bins[1:]+arr_bins[:-1])
        print 'Bin Size: {0}, {1}, {2}'.format(binsize, arr_bins.min(), arr_bins.max())
        # in klamda
        print 'UV Dist: {0:.1f} - {1:.1f} klam'.format(uvmin, uvmax)
        # prepare the data structures to store result in
        uvampdat = _sp.zeros([nbin-1,3,3])  # bins for real, img, amp with n points, mean, dispersion
        expt = _sp.zeros(nbin-1)            # Expected value for no signal
        sn = _sp.zeros(nbin-1)              # signal to noise
        for ibin in range(nbin-1):
            # ibin - index of current working bin
            # to store stuff in the prepared arrays
            minbin = arr_bins[ibin]
            maxbin = arr_bins[ibin]+binsize
            isubs = ( (uvdist < maxbin) == (uvdist >= minbin) ).nonzero()[0]
            npoints = len(isubs)
            if npoints > 0:
                # real
                reals = self.re[isubs]
                wts = self.wt[isubs]
                # Get rid of the negative weights
                wtsubs = (wts >= 0.0).nonzero()[0]
                wts = wts[wtsubs]
                reals = reals[wtsubs]
                npts = int(len(wtsubs))
                # points in each interval(?)
                uvampdat[ibin,0,0] = int(len(wtsubs))
                # mean real value
                uvampdat[ibin,0,1] = (reals).sum()/(npoints)
                uvampdat[ibin,0,2] = _sp.sqrt(( (reals*reals).sum() - (npoints*uvampdat[ibin,0,1]*uvampdat[ibin,0,1]))/(npoints-1))
                # Imaginary
                reals = self.im[isubs]
                uvampdat[ibin,1,0] = int(len(wtsubs))
                uvampdat[ibin,1,1] = (reals).sum()/(npoints)
                uvampdat[ibin,1,2] = _sp.sqrt(( (reals*reals).sum() - (npoints*uvampdat[ibin,1,1]*uvampdat[ibin,1,1]))/(npoints-1))
                # amplitudes
                reals = self.amp[isubs]
                # take real and imaginary part, calculate amplitude
                uvampdat[ibin,2,0] = int(len(wtsubs))
                x = uvampdat[ibin,0,1]
                xerr = uvampdat[ibin,0,2]
                y = uvampdat[ibin,1,1]
                yerr = uvampdat[ibin,1,2]
                temp_amp = _sp.sqrt(x*x + y*y)
                uvampdat[ibin,2,1] = temp_amp
                sigtot = (x*xerr/(temp_amp))**(2.0) + (y*yerr/(temp_amp))**(2.0)
                uvampdat[ibin,2,2] =  _sp.sqrt(sigtot/(npts-2))
                if uvampdat[ibin,2,2] > 0.0:
                    sn[ibin] = temp_amp/uvampdat[ibin,2,2]
                else:
                    sn[ibin] = 0.0
                expt[ibin] = (_sp.sqrt(_sp.pi/2.0))*uvampdat[ibin,2,2]
            else:
                pass
                
        BinnedDMC.nbin = nbin
        BinnedDMC.bins = arr_bins1
        BinnedDMC.uvamp = uvampdat
        BinnedDMC.expt = expt
        self.BinnedDMC = BinnedDMC

    def bin_data(self, ruv=None, binsize=10, nbins=30, ignore_wt=False):
        """
        Function to bin UV data, both vector and scalar average is calculated
        creates Bin.Sca and Bin.Vec objects in self
        
        needs : uvdist_klam
                re, im, wt, amp
        TODO: move core calcs to separate function for reuse
        """
     
        if self.__dict__.has_key('Model'):
            if ruv is not None:
                uvdist = ruv
            else:
                uvdist = self.Model.uvdist_klam
            mre = self.Model.re
            mim = self.Model.im
            if ignore_wt:
                mwt = _sp.ones_like(self.Model.re)
            else:
                mwt = self.Model.wt
            self.Model.BinVec = uv_bin_vector(uvdist, re, im, mwt, binsize=binsize, nbins=nbins)
            self.Model.BinSca = uv_bin_scalar(uvdist, re, im, mwt, binsize=binsize, nbins=nbins)
        if ruv is not None:
            uvdist = ruv
        else:
            uvdist = self.uvdist_klam
        if ignore_wt:
                wt = _sp.ones_like(self.re)
        else:
                wt = self.wt
        self.BinVec = uv_bin_vector(uvdist, self.re, self.im, wt, binsize=binsize, nbins=nbins)
        self.BinSca = uv_bin_scalar(uvdist, self.re, self.im, wt, binsize=binsize, nbins=nbins)
        
    def shift(self, offset):
        if 'isshifted' in self.__dict__.keys():
            print('Data already shifted once, aborting')
            return False
        # create input
        uv = _sp.array([self.u_klam, self.v_klam])
        reim = _sp.array([self.re, self.im])
        # shift/translate the data
        phas, re, im = translate(self.u_klam, self.v_klam, offset)
        # store in self
        self.re = re
        self.im = im
        self.isshifted = (True, offset)
        print('Only shifts the current object')

    def rotate(self, deg, force=False):
        if 'isrotated' in self.__dict__.keys() and not force:
            print('Data already rotated once, aborting')
            return False
        # kilo-lambda
        uv = _sp.array([self.u_klam, self.v_klam])
        self.u_klam, self.v_klam = rotate_field(uv, deg)
        # meters
        uv = _sp.array([self.u_m, self.v_m])
        self.u_m, self.v_m = rotate_field(uv, deg)
        # nano seconds
        #~ uv = _sp.array([self.u_nsec, self.v_nsec])
        #~ self.u_nsec, self.v_nsec = rotate_field(uv, deg)
        # set isrotated to True, so we can check how much it was rotated
        self.isrotated = (True, deg)
        # print('Only rotates the current object.')
        # calculate new uvdistances
        # print('Calculating new uvdistances, better way to implement this,\
        #  need base entities that are used on the fly to calc. other')
        #~ self.uvdist_nsec = _sp.sqrt(self.u_nsec**2 +self.v_nsec**2 + self.w_nsec**2)
        self.uvdist_klam = _sp.sqrt(self.u_klam**2 +self.v_klam**2  + self.w_klam**2)
        self.uvdist_m = _sp.sqrt(self.u_m**2 +self.v_m**2 + self.w_m**2)
        
    def incline(self, deg):
        # print('incline only works on kilo-lamda uv points')
        # klam
        uv = _sp.array([self.u_klam, self.v_klam])
        self.uvdist_klam = incline(uv, deg)
        # set flag isinclined to True and the amount of degrees
        self.isinclined = (True, deg)
        # print('Only inclines the current object, if binned b4, bin again')
        
    def deproject(self, PA, inc):
        uvin = _sp.array([self.u_klam, self.v_klam])
        newuv_klam, ruv_klam = deproject(uvin, PA, inc)
        
        # store original values
        class UV_original: pass
        UV_original.uvdist_klam = self.uvdist_klam[:]
        UV_original.u_klam = self.u_klam[:]
        UV_original.v_klam = self.v_klam[:]
        self.UV_original = UV_original
        
        self.uvdist_klam = ruv_klam
        self.u_klam = newuv_klam[0]
        self.v_klam = newuv_klam[1]
        self.uvdist_lam = ruv_klam*1e3
        self.u_lam = newuv_klam[0]*1e3
        self.v_lam = newuv_klam[1]*1e3
        self.uvdist_m = ruv_klam / (self.freq * 1.0e-3) * co.c.cgs.value * 1.0e-2
        self.u_m = newuv_klam[0] / (self.freq * 1.0e-3) * co.c.cgs.value * 1.0e-2
        self.v_m = newuv_klam[1] / (self.freq * 1.0e-3) * co.c.cgs.value * 1.0e-2
        
        self.isrotated = (True, PA)
        self.isinclined = (True, inc)
        
    def reset_projection(self):
        
        if 'UV_original' in self.__dict__.keys():
            self.uvdist_klam = self.UV_original.uvdist_klam[:]
            self.u_klam = self.UV_original.u_klam[:]
            self.v_klam = self.UV_original.v_klam[:]
            self.uvdist_lam = self.UV_original.uvdist_klam[:]*1e3
            self.u_lam = self.UV_original.u_klam[:]*1e3
            self.v_lam = self.UV_original.v_klam[:]*1e3
            self.uvdist_m = self.UV_original.uvdist_klam[:] / (self.freq * 1.0e-3) * co.c.cgs.value * 1.0e-2
            self.u_m = self.UV_original.u_klam[:] / (self.freq * 1.0e-3) * co.c.cgs.value * 1.0e-2
            self.v_m = self.UV_original.v_klam[:] / (self.freq * 1.0e-3) * co.c.cgs.value * 1.0e-2
        else:
            print self.isrotated
            print self.isinclined
            raise StandardError('No orignal UV points saved.')
            
        
        
    def __str__():
        return 'Not implemented yet...'


#########################
