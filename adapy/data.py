
from .helpers import *


########################################################################
# USEFUL STRINGS
KMS = u"km\u00b7s\u207b\u00b9"



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
        #Isum = imgs.sum(axis=0)*abs(Fits.v_cdeltkms) # make a copy for masking <3sigma values in mom1 map
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
    # perhaps need some other input name -> its ok!
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
        from libs.cgsconst import CC
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
                    print " Id molecule : {0}".format(self.lines[0][j])
                    print u" Frequency shift : {0:2.5} GHz Vel shift : {1:5.5} {2}".format(freq_shift,vel_shift,KMS)
                frequency_string =\
                u' Frequency : {0:3.9f} GHz (v_sys corrected)'.format(frequency_corrected)
                print stylify(frequency_string,fg='b')
                print u' FWHM      : {0}, {1} ({2}) (0-based) ([{3:.2f}, {4:.2f}] {5})'.format(channels_half.min(), channels_half.max(),(channels_half.max()-channels_half.min()+1), lower_half, upper_half,KMS)
                print u' \u00b1FWHM   : {0}, {1} ({2}) (0-based) ([{3:.2f}, {4:.2f}] {5})'.format(channels.min(), channels.max(), (channels.max()-channels.min()+1), lower,upper,KMS)
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
            print('Spectrum not fitted, aborting line id')
            raise ParError('linefit - need to fit lines')
        if hasattr(self, 'lines'):
            print('Line list exists, overwriting...')
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
            result = spl.search(
                        freq=frequency_pairs[i],
                        linelist=['jpl','cdms'])
                        #e_to=1000)
            if args['writetofile']:
                with open(args['writetofile'],'a') as f:
                    f.write('\n# Line no. {0}\n'.format(i+1))
            if result!=None:
                species, freq = results['species'],result[]
                smu2, eu = result[13], result[20]
                uresqnr = result[10]
                llist = result[24]
                #~ species, freq = result[1],result[3]
                #~ smu2, eu = result[13], result[20]
                #~ uresqnr = result[10]
                #~ llist = result[24]
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


