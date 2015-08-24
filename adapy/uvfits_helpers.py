
import astropy.units as u
import astropy.constants as co


########################################################################
# Define new uv-coordinate units

lambdas = u.def_unit('lambdas', format={'format' : r'\lambda'})
klambdas = u.def_unit('kilolambdas', 1e3 * lambdas, format={'format' : r'k\lambda'})

# equivalency (from_unit, to_unit, forward, backward)
def lambdas_equivalencies(restunit):
    try:
        restfreq_hz = restunit.to(u.Hz, equivalencies=u.spectral())
    except (AttributeError):
        raise AttributeError('input \'restfreq\' should be a spectral unit.')
    eq = [
        (lambdas, u.s, lambda x: x/restfreq_hz, lambda x: x*restfreq_hz),
        (lambdas, u.m, lambda x: x/restfreq_hz * co.c.to(u.m/u.s).value, lambda x: x/co.c.to(u.m/u.s).value * restfreq_hz),
        (u.m, u.s, lambda x: x/co.c.to(u.m/u.s).value, lambda x: x*co.c.to(u.m/u.s).value),
        ]
    return eq


########################################################################


# BINNING
def uv_bin_vector(uvdist, re, im, wt, start='zero', binsize=10, nbins=50, weighted=False):
    """

    Vector averaging of amplitudes
    Calculate the binned amplitude and various related things.
    The binned amplitude is calculated from the real and imaginary
    part of the visibilities.
    
    
    zero : where to start the binning, at uvdist = 0, or min(uvdist)

    Description
    ===========
    Vector binning of the amplitude, first bin the real and
    imaginary parts as
    RE = sum(RE_i)/Np (Mean) (not the weighted mean, should be?)
    RE_sig = sqrt( ( sum(RE_i^2) - Np * RE^2 ) / (Np - 1) )
    (-1 because we determined the mean)
    the same goes for the Imaginary part.
    
    the square root of the sum of the squared real and imaginary parts
    A = sqrt( RE^2 + IM^2 )

    and the error propagation of the
    variance of the real and imaginary parts
    A_sig = sqrt( ((RE*RE_sig/A)^2 + (IM*IM_sig/A)^2) / (Np - 2) )
    

    NOTES
    =====
    Some parts divide with possible zeros.
    This operation will create 'nan' entries which messes things
    when running calculations on it, i.e. nan*2 = nan
    therefore I've started using masked arrays.

    
    """
    class Binned_Vector(object):
            pass
    ##### CORE CALC START #####
    if start in ['zero', 0, '0']:
        uvmin = 0.0
    elif start in ['min']:
        uvmin = uvdist.min()
    uvmax = uvmin + binsize * (int(nbins) + 0.5)
    # Define the bins, from uvmin to uvmax
    arr_bins = _sp.arange(_sp.floor(uvmin),
                            _sp.ceil(uvmax),
                            binsize)
    # mid-points of the bins
    arr_bins1 = 0.5*(arr_bins[1:] + arr_bins[:-1])
    minmax = zip(arr_bins[:-1], arr_bins[1:])
    # only choose data with positive weigths
    pos_wt = wt>= 0.0
    def filter_points(i,j):
        # find the indices of data within the limits and
        # with positive weigths
        # i:lower boundary, j:upper boundary
        return ( (uvdist>=i) * (uvdist<j) * pos_wt).nonzero()[0]
    isubs = [filter_points(i,j) for i,j in minmax]
    npoints = _sp.array([len(i) for i in isubs])       # points in each interval
    ###########################################################################
    # Real and Imaginary data binning
    data = _sp.array([re[:], im[:]])
    # mean for Re and Im separately
    if not weighted:
        # changed this
        #~ data_mean = _sp.array([data[:,i].mean(axis=1) for i in isubs])
        # to this instead, got some problems with not explicitly setting
        # value to nan, do not rely on mean or average to take care of it
        # here we use a masked array, which handles the nan much better
        # removed masked arrays again!
        # print('method takes into account nan raw data values')
        data_mean =  _sp.array([_sp.nanmean(data[:,i], axis=1) if j>0 else _sp.array([_sp.nan, _sp.nan]) for i,j in zip(isubs,npoints)])
        #~ check = _sp.array([npoints==0, npoints==0]).T
        #~ data_mean = _sp.ma.masked_array(data_mean, check, fill_value=_sp.nan)
        #~ data_mean = _sp.ma.masked_array(data_mean,_sp.isnan(data_mean), fill_value=_sp.nan)
    elif weighted:
        print ('Calculating weighted average real and imaginary amplitudes')
        print ('method does not take into account nan raw data values')
        # doesn't work if isubs is empty somewhere
        #~ data_mean = _sp.array([_sp.average(data[:,i], axis=1, weights=wt[i]) for i in isubs])
        data_mean = _sp.array([_sp.average(data[:,i], axis=1, weights=wt[i]) if j>0 else _sp.array([_sp.nan, _sp.nan]) for i,j in zip(isubs,npoints)], dtype='float')
        #~ data_mean = _sp.ma.masked_array(data_mean,_sp.isnan(data_mean))
    # Error of real and imaginary data
    #~ data_var = _sp.array([_sp.var(data[:,i], ddof=1, axis=1) if j>0 else _sp.array([_sp.nan, _sp.nan]) for i,j in zip(isubs,npoints)])
    # ddof = j-1, the number of points in bin, minus the parameter determined
    # i.e. the mean.
    data_var = _sp.array([_sp.var(data[:,i], ddof=1, axis=1) if j>0 else _sp.array([_sp.nan, _sp.nan]) for i,j in zip(isubs,npoints)])
    #~ data_var = _sp.ma.masked_array(data_var,_sp.isnan(data_var), fill_value=_sp.nan)

    #~ data_var[data_var==-0] = _sp.nan
    #~ data_var = _sp.ma.masked_array(data_var, data_var==-0, fill_value=_sp.nan)

    #~ data_std = _sp.ma.sqrt(data_var)
    data_std = data_var**0.5 # used sqrt here before, got floating point error
    # Amplitude binning
    amp_mean = ( (data_mean**2).sum(axis=1) )**0.5
    #~ amp_mean = _sp.sqrt( (data_mean**2).sum(axis=1) )
    amp_tmp = amp_mean.reshape((len(amp_mean), 1))
    # calculate the variance of the amplitude
    # error propagation of the variance of the imaginary
    # and real parts
    # if the amp_temp is close to zero, we can end up with
    # inf in variance.
    pars=2
    dof = _sp.array([float(i-pars) if i>0 else 0.0 for i in npoints])
    amp_var = (( ( data_mean * data_var / amp_tmp )**2).sum(axis=1)
        / ( dof ) )**0.5
    amp_std = amp_var**0.5
    # Signal to Noise (SNR), _sp.divide gives 0 when dividing by 0
    amp_snr = _sp.divide( amp_mean, amp_var )
    # expectation value when no signal
    amp_expt = _sp.sqrt( _sp.pi / 2. ) * amp_var
    ###########################################################################
    # get the binned real and imaginary parts
    re_mean, im_mean = data_mean.T
    re_var, im_var = data_var.T
    re_std, im_std = data_std.T
    # Signal to Noise (SNR), _sp.divide gives 0 when dividing by 0
    re_snr, im_snr = _sp.divide(data_mean, data_var).T
    # expectation value when no signal
    re_expt, im_expt = _sp.sqrt( _sp.pi / 2. ) * data_var.T

    ###########################################################################
    ##### CORE CALC END #####

    # store in class        
    #~ Binned_Vector.nbin = nbin
    Binned_Vector.npoints = npoints
    Binned_Vector.bins = arr_bins1
    Binned_Vector.uvdist_klam = arr_bins1
    Binned_Vector.amp = amp_mean
    Binned_Vector.amp_var = amp_var
    Binned_Vector.amp_std = amp_std
    Binned_Vector.data = data_mean
    Binned_Vector.data_var = data_var
    Binned_Vector.data_std = data_std
    Binned_Vector.snr = amp_snr
    Binned_Vector.expt = amp_expt
    # store Re in class        
    Binned_Vector.re = re_mean
    Binned_Vector.re_var = re_var
    Binned_Vector.re_snr = re_snr
    Binned_Vector.re_expt = re_expt
    # store Im in class        
    Binned_Vector.im = im_mean
    Binned_Vector.im_var = im_var
    Binned_Vector.im_snr = im_snr
    Binned_Vector.im_expt = im_expt
    
    return Binned_Vector

def uv_bin_scalar(uvdist, re, im, wt, start='zero', binsize=10, nbins=50, weighted=False):
    """
    Scalar averaging amplitudes
    
    NOTES
    =====
    Some parts divide with possible zeros.
    This operation will create 'nan' entries which messes things
    when running calculations on it, i.e. nan*2 = nan
    therefore I've started using masked arrays.
    
    """
    class Binned_Scalar(object):
            pass
    ##### CORE CALC START #####
    if start in ['zero', 0, '0']:
        uvmin = 0.0
    elif start in ['min']:
        uvmin = uvdist.min()
    uvmax = uvmin + binsize * (int(nbins) + 0.5)
    #~ binsize = int(round(((uvmax-uvmin)/nbins), 0 ))
    # Define the bins, from uvmin to uvmax
    #~ arr_bins = _sp.linspace(uvmin, uvmax, nbins)
    arr_bins = _sp.arange(_sp.floor(uvmin),
                            _sp.ceil(uvmax),
                            binsize)
    #~ print ('{0} bins with {1} binsize'.format(nbins, binsize))
    #~ print len(arr_bins)
    # mid-points of the bins
    arr_bins1 = 0.5*(arr_bins[1:] + arr_bins[:-1])
    minmax = zip(arr_bins[:-1], arr_bins[1:])
    # only choose data with positive weigths
    pos_wt = wt>= 0.0
    def filter_points(i,j):
        # find the indices of data within the limits and
        # with positive weigths
        # i:lower boundary, j:upper boundary
        return ( (uvdist>=i) * (uvdist<j) * pos_wt).nonzero()[0]
    isubs = [filter_points(i,j) for i,j in minmax]
    npoints = _sp.array([len(i) for i in isubs])       # points in each interval
    ###########################################################################
    # AMPLITUDE 
    amp = _sp.sqrt(re**2 + im**2)
    # PHASE
    pha = _sp.arctan2(im, re)
    pha_deg = pha / _sp.pi * 180.
    # put amp and pha in an array
    data = _sp.array([amp, pha])
    # ERROR / SIGMA
    #TODO : check
    # following 1.0e6 is just for GILDAS, change if needed
    #~ print('NB : Error calculated from weights assuming GILDAS '
    #~ 'data (i.e. frequencies in MHz).')
    #~ self.sigma = 1/sqrt(self.wt*1.0e6)
    # Daniels way of calculating sigma
    # test this first
    sigma = _sp.sqrt(0.5 / ( wt * float(amp.shape[0]) ) )
    
    if not weighted:
        data_mean = _sp.array([_sp.nanmean(data[:,i], axis=1) if j>0 else _sp.array([_sp.nan, _sp.nan]) for i,j in zip(isubs,npoints)])
        # above operation create some 'nan' entries which messes things up
        # later on if we're not using a masked array
        #~ data_mean = _sp.ma.masked_array(data_mean,_sp.isnan(data_mean))
    elif weighted:
        print ('Calculating weighted average amplitudes')
        data_mean = _sp.array([_sp.average(data[:,i], axis=1, weights=wt[i]) if j>0 else _sp.array([_sp.nan, _sp.nan])  for i,j in zip(isubs,npoints)])
        #~ data_mean = _sp.ma.masked_array(data_mean,_sp.isnan(data_mean), fill_value=_sp.nan)
    # variance and standard deviation
    data_var = _sp.array([_sp.var(data[:,i], ddof=1, axis=1) if j>0 else _sp.array([_sp.nan, _sp.nan]) for i,j in zip(isubs,npoints)])
    #~ data_var[data_var==-0] = _sp.nan
    #~ data_var = _sp.ma.masked_array(data_var, data_var == -0, fill_value=_sp.nan)
    data_std = data_var**0.5
    
    amp_mean, pha_mean = data_mean.T
    amp_var, pha_var = data_var.T
    amp_std, pha_std = data_std.T
    
    # Signal to Noise (SNR), _sp.divide gives 0 when dividing by 0
    amp_snr, pha_snr = _sp.divide(data_mean, data_var).T
    # expectation value when no signal
    amp_expt, pha_expt = _sp.sqrt( _sp.pi / 2. ) * data_var.T
    
    ###########################################################################
    ##### CORE CALC END #####

    # store in class        
    #~ Binned.nbin = nbin
    Binned_Scalar.npoints = npoints
    Binned_Scalar.bins = arr_bins1
    Binned_Scalar.uvdist_klam = arr_bins1
    Binned_Scalar.data = data_mean
    Binned_Scalar.data_var = data_var
    Binned_Scalar.amp = amp_mean
    Binned_Scalar.amp_var = amp_var
    Binned_Scalar.amp_snr = amp_snr
    Binned_Scalar.amp_expt = amp_expt
    Binned_Scalar.pha = pha_mean
    Binned_Scalar.pha_var = pha_var
    Binned_Scalar.pha_snr = pha_snr
    Binned_Scalar.pha_expt = pha_expt
    
    # send back an object with all the data structures
    return Binned_Scalar

# MANIPULATING, DE-PROJECT, ROTATE, TRANSLATE etc
def translate(uv_klam, reim, offset):
    print ('Has a -1 multiplied here, is this right?')
    phas = -1.0*( ((uv_klam[0])*(offset[0]/pc2au)) +
            ((uv_klam[1])*(offset[1]/pc2au)))*2.0*pi
    re = (reim[0]*_sp.cos(phas)) - (reim[1]*_sp.sin(phas))
    im = (reim[0]*_sp.sin(phas)) + (reim[1]*_sp.cos(phas))
    return phas, re, im

def deproject(uv, PA, inc):
    """
    Rotate and deproject individual visibility coordinates.
    From Hughes et al. (2007) - "AN INNER HOLE IN THE DISK AROUND 
    TW HYDRAE RESOLVED IN 7 mm DUST EMISSION".
    """
    R = ( (uv**2).sum(axis=0) )**0.5
    #~ phi = _sp.arctan(uv[1]/uv[0] - deg2rad(PA))
    phi = _sp.arctan2(uv[1],uv[0]) - deg2rad(PA)
    #~ phi = _sp.arctan2( (uv[1] - deg2rad(PA) * uv[0]) , uv[0])
    newu = R * _sp.cos(phi) * _sp.cos( deg2rad(inc) )
    newv = R * _sp.sin(phi)
    newuv = _sp.array([newu, newv])
    ruv = (newuv**2).sum(axis=0)**.5
    return newuv, ruv

def rotate_field(uv, PA, U_RA_align = True):
    """
    Rotates a coordinate system (UV plane) by PA
    degrees.
    uv : 2-D array with uv[0] U and uv[1] coordinated
    PA : Position Angle, in degrees
    U_RA_align : for ALMA and PdBI the U-axis and RA are aligned
                 and thus one form of the equation must be used
                 While for SMA/CARMA (USA, meh), they are not aligned
                 and thus some sign changes have to impliemented from
                 that presented in Berger & Segransan (2007)
    
    """
    direction =  [-1, 1][int(U_RA_align)]
    u_new = uv[0] * _sp.cos( deg2rad(PA) ) + direction * uv[1] * _sp.sin( deg2rad(PA) )
    v_new = -1 * direction * uv[0] * _sp.sin( deg2rad(PA) ) + uv[1] * _sp.cos( deg2rad(PA) )
    return u_new, v_new
    
def incline(uv, inc):
    #~ ruv = ( uv[0]**2 + (uv[1] * _sp.cos(deg2rad(inc)) )**2  )**.5 
    # the PA goes from North to East in the image plane, and the 
    # Major axis is flipped 90 degrees going from
    # image to UV plane (Major in image in minor in UV)
    ruv = ( uv[0]**2 * _sp.cos(deg2rad(inc))**2 + uv[1]**2  )**.5 
    return ruv

    
