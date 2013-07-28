
### Imports
# python builtins
import os as _os
import sys as _sys
import subprocess as _subprocess

# python extra modules
import scipy as _scipy
from matplotlib import pyplot as _pl

# adapy internal imports
from ...libs import cgsconst as _cgs
from ... import moldata as _moldata
from ...tools import get_colors as _gc

#~ from ...views import set_rc
#~ _pl.rcParams['text.usetex'] = True
#~ _pl.rcParams['savefig.dpi'] = 300
#~ _pl.rcParams['figure.dpi'] = 150
_pl.rcParams['figure.facecolor'] = 'w'
_pl.rcParams['font.size'] = 10
#~ set_rc()

# perhaps put help functions in separate file?
# help.py, extra.py,

### Help functions
def bplanck(nu, T):
    """
        Returns the Spectral Radiance (Planck curve)
    """
    from scipy import exp, constants
    try:
        x = _cgs.HH * nu / (_cgs.KK * T)
        bpl = (2.0 * _cgs.HH * nu**3 / _cgs.CC**2) / (exp(x)-1.0)
        return bpl
    except (ZeroDivisionError):
        return 0

def find_intensity(fitsfile, interval = [], nsig = 3):
    from scipy import array, arange, where, log, pi, meshgrid
    import matplotlib.pyplot as pl
    from pyfits import getdata, getheader
    from adapy.adacore import gaussfit2d
    pl.ion()
    class SkyModel: pass

    data = getdata(fitsfile)
    SkyModel.continuum = data[0]
    #~ data -= data[0] # remove continuum
    SkyModel.data = data - SkyModel.continuum
    header = getheader(fitsfile)
    
    # get some of the stuff from the header
    SkyModel.bunit = header['BUNIT']
    SkyModel.restfreq = header['RESTFREQ']
    
    # create the velocity array, just for fun
    v_cdelt = header['CDELT3']*1e-3     # in km/s
    v_crpix = header['CRPIX3']
    v_crval = header['CRVAL3']
    v_naxis = header['NAXIS3']
    
    v_array = arange(v_naxis) - v_crpix
    v_array *= v_cdelt
    SkyModel.v_array = v_array + v_crval
    SkyModel.v_cdelt = v_cdelt         # in km/s
    
    
    SkyModel.ra_cdelt = header['CDELT1']*3600
    SkyModel.dec_cdelt = header['CDELT2']*3600
    SkyModel.ra_array = ((arange(header['NAXIS1']) - header['CRPIX1']) * header['CDELT1']*3600) + header['CRVAL1']
    SkyModel.dec_array = ((arange(header['NAXIS2']) - header['CRPIX2']) * header['CDELT2']*3600) + header['CRVAL2']
    
    # assume model peak in center
    z, y ,x = SkyModel.data.shape
    SkyModel.spectrum = SkyModel.data[:, y/2, x/2]
    
    if len(SkyModel.data.shape) < 3: # needs to be a cube for analysis to work
        print("Wrong data shape of input fits file")
        return 0
    if interval == []: # if no interval given, need to do it interactively
        from adapy.adacore import fit_gauss1d as gaussfit
        #~ from pylab import ginput
        #~ from matplotlib.widgets import Cursor
        #~ fig = pl.figure()
        #~ ax = fig.add_subplot(111) 
        #~ ax.plot(SkyModel.v_array, SkyModel.spectrum)
        #~ cursor = Cursor(ax, color='red', linewidth=2 )
        #~ print('Click on lower limit')
        #~ x_low  = ginput()[0][0]
        #~ print('Click on upper limit')
        #~ x_high  = ginput()[0][0]
        #~ fig.close()
        #~ SkyModel.mom0 = 
        # simple guesses/assumptions, 
        # perhaps extend to calculate moment0/1 for the position
        # and width of the distribution?
        datamax = SkyModel.spectrum.max()
        # where is the max?
        ii = where(SkyModel.spectrum.max() == SkyModel.spectrum)[0] 
        width_estimate = (SkyModel.v_array.max() - SkyModel.v_array.min()) * 0.2
        # fit a 1D Gaussian
        results_1d = gaussfit((SkyModel.v_array, SkyModel.spectrum), 
                            params=(
                                    datamax,
                                    SkyModel.v_array[ii], 
                                    width_estimate
                                    ),
                            verbose=0
                            )[0]
        SkyModel.results_1d = results_1d
        amplitude_1d = results_1d[2]
        position_1d = results_1d[1]
        fwhm_1d = results_1d[2]
        sigmafromfwhm = 1 / (2 * (2 * log(2))**.5)
        sigma_1d = fwhm_1d * sigmafromfwhm
        
        interval = position_1d + array([-1, 1]) * nsig * sigma_1d
        print("Integration interval : 1D Gaussian fit"
                " (+/- {0} sigma)".format(nsig))
        SkyModel.interval = interval
    else:
        SkyModel.interval = interval
        print("Integration interval : input")
        
    indices = where(
                        (SkyModel.v_array >= interval[0]) * 
                        (SkyModel.v_array <= interval[1])
                        )
    
    SkyModel.zero = SkyModel.data[indices].sum(axis=0) * abs(ModelData.v_cdelt)
    X, Y = meshgrid(arange(header['NAXIS1']),arange(header['NAXIS2']))
    #~ results_2d = gaussfit2d((X, Y, ModelData.zero), params=(0.0, (0.1, 64, 64, 2, 2, 0)))[0]
    results_2d = gaussfit2d((X, Y, SkyModel.zero), fitheight=0)[0]
    # Volume (integral) of the Gaussian
    # V = 2 pi Amp sigma1 sigma2
    SkyModel.amplitude_2d = results_2d[0]
    SkyModel.sigma_x = results_2d[3]
    SkyModel.sigma_y = results_2d[4]
    SkyModel.results_2d = results_2d
    
    SkyModel.intensity = pi * 2 * SkyModel.amplitude_2d * SkyModel.sigma_x * SkyModel.sigma_y
    print('Integrated intensity : {0:.2f} Jy'.format(SkyModel.intensity))
    return SkyModel

def _tex_transition(n1, n2, g1, g2, nu):
    """ 
    Calculate the excitation temperature of transition from 2 to 1
    """
    
    numer = _cgs.HH * nu
    denom = _cgs.KK * _scipy.log(n1 * g2 / (n2 * g1))
    return numer / denom

def _pdfcheck(pdf):
        if not pdf: 
            _plt.ion()
        elif pdf:
            _plt.ioff()

def _pdfsave(pdf, pdfname, **kwargs):
    if pdf:
            _plt.savefig('{0}.pdf'.format(str(pdfname)), kwargs)

def read_molfile(filepath, output='dictionary'):
    return True

########################################################################
# OLD HELP FUNCTIONS

def read_transphereoutput(self, ext = 0):
    from scipy import array
    if ext == 0: ext=''
    filetoread = 'envstruct' + ext + '.dat'
    path_to_file = _os.path.join(self.directory, filetoread)
    with open(path_to_file, 'r') as f:
        lines = f.read().split('\n')
    nr = int(lines[0])
    dat_envstruct = array([i.split() for i in lines[2:nr + 2]], dtype='float')


    filetoread = 'spectrum.dat'
    path_to_file = _os.path.join(self.directory, filetoread)
    with open(path_to_file, 'r') as f:
        lines = f.read().split('\n')
    nr = int(lines[0])
    dat_spectrum = array([i.split() for i in lines[2:nr+2]], dtype='float')

    #~ nr = int(f.readline().strip())
    #~ dat = np.zeros((3,nr),float)
    #~ for ii in range(0,3):
        #~ for jj in range(0,nr):
            #~ dum = f.readline().strip()
            #~ dat[ii,jj] = dum
    #~ f.close()

    class Envstruct:
        r = dat_envstruct[:,0]
        rho_dust = dat_envstruct[:,1]
        temp = dat_envstruct[:,2]

    #~ class Spectrum:
    Envstruct.frequency = dat_spectrum[:,0]
    Envstruct.intensity = dat_spectrum[:,1]
    #~ Envstruct.Spectrum = Spectrum
    #~ self.Envstruct = Envstruct

    #~ import numpy as np
    filetoread = 'convhist.info'
    path_to_file = _os.path.join(self.directory, filetoread)
    f = open(path_to_file, 'r')
    nn = int(f.readline().strip().split()[0])
    f.close()

    # Convergence history
    filetoread = 'convhist.dat'
    path_to_file = _os.path.join(self.directory, filetoread)
    with open(path_to_file, 'r') as f:
        lines = f.read().split('\n')
    nr = int(lines[0].strip())
    if nr == 0: raise Exception('Nothing run, no convergence history.')
    x1 = nr+1

    #These need to depend on value of nr
    dat1 = array([i.split() for i in lines[1:x1]], dtype='float')
    dat2 = array([i.split() for i in lines[x1+1:x1*2]], dtype='float')
    dat3 = array([i.split() for i in lines[x1*2+1:x1*3]], dtype='float')

    dat = array([dat1,dat2,dat3])

    #~ f = open('convhist.dat','r')
    #~ nr = int(f.readline().strip())
    #~ dat = np.zeros((9,nn,nr),float)
    #~ for jj in range(0,nn):
        #~ for kk in range(0,nr):
            #~ dum = f.readline().strip().split()
            #~ if dum == []: dum=f.readline().strip().split()
            #~ dat[0:9,jj,kk]=np.array(dum,dtype=float)
    #~ f.close()

#    if nn gt 1 then idx=[1,2,0] else idx=[1,0]. Note transpose commands not executed...
    class Convhist:
        temp=dat[:,:,0]
        jjme=dat[:,:,1]
        hhme=dat[:,:,2]
        jj=  dat[:,:,3]
        hh=  dat[:,:,4]
        kapt=dat[:,:,5]
        kapj=dat[:,:,6]
        kaph=dat[:,:,7]
        fj=  dat[:,:,8]
    #~ self.Convhist = Convhist

    #~ f = open('envstruct.inp')
    #~ nr = int(f.readline().strip())
    #~ dat = np.zeros((3,nr),float)
    #~ for ii in range(0,nr):
        #~ dum=f.readline().strip().split()
        #~ if dum == []: dum=f.readline().strip().split()
        #~ dat[0:3,ii]=np.array(dum,dtype=float)
    #~ r=dat[0,:]
    #~ f.close()

    #~ convhist={'r': r, 'temp': temp, 'jjme': jjme, 'hhme': hhme, 'jj': jj, 'hh': hh, 'kapt': kapt, 'kapj': kapj, 'kaph': kaph, 'fj': fj}
    #~ self.Envstruct = envstruct
    #~ self.convhist = convhist
    return Envstruct, Convhist

def create_grid(r_in, r_out, nshell, space = 'powerlaw1', end = True):
    # function to create grid
    if space == 'log10':
        from scipy import log10, logspace
        # get the exponent of the start- and
        # stop-radius in input units
        start = [log10(r_in), 0][r_in == 0]
        stop = log10(r_out)
        radii = logspace(start, stop, num=nshell, endpoint=end)
    elif space == "powerlaw1":
        from scipy import arange
        radii = r_in * (r_out/r_in)**(arange(nshell)/(nshell - 1.0))
    elif space == 'linear':
        from scipy import linspace
        # linearly spaced grid
        radii = linspace(r_in, r_out, num=nshell, endpoint=end)
    elif space == 'powerlaw2':
        from scipy import linspace
        # first check if coefficients to the power-law was given
        #~ if 'exp' in kwargs:
            #~ p_exp = kwargs['exp']
        #~ else: # if not, set it to 2, i.e. r^2
            #~ p_exp = 2
        radii = r_in + (r_out - r_in)*(linspace(r_in, r_out, num=nshell, endpoint=end)/(r_out))**2
        #pr_int('Not implemented yet.')
        #raise ParError(spaced)
    else:
        raise Exception(space)
    return radii

# FIXME, does not work tries to import old adavis module
def plot_spectrum(freq, intensity, dpc = 0, jy = 0, pstyle = '', xlog = 1, ylog = 1):
    import sys
    import matplotlib.pyplot as pl
    pl.ion()
    #~ from ..views import set_rc
    #~ set_rc
    xcoord = 1.0e4 * _cgs.CC / freq

    if dpc == 0: sys.exit('Error: distance needs to be set when plotting flux')

    distfact = 1.e0/ (dpc**2)

    if jy != 0:
        lumfact = 1e+23
    else:
        lumfact = freq

    pl.plot(xcoord, distfact * lumfact * intensity, pstyle)
    pl.xlabel(r'$\lambda\, [\mu \mathrm{m}]$')

    if jy != 0:
        pl.ylabel(r'$F_\nu$\, [Jy]')
    else:
        pl.ylabel(r'$\nu F_\nu \, [\mathrm{erg cm}^{-2}\, \mathrm{s}^{-1}]$')

    if xlog == 1: pl.xscale('log')
    if ylog == 1: pl.yscale('log')

def write_ratraninput(self):
    # ugly hack, need to restructure the whole radiative transfer module
    from scipy import log, log10, pi, array, arange
    
    
    input_dir_path = _os.path.join(_os.getcwd(), self.directory)
    # if the directory exists
    if not make_dirs(input_dir_path): # will raise error if things go south (i.e., permissions not correct [I think...])
        print('Directory exists, continuing.')
    
    V = 4 * pi * ((self.ra**3  - self.ra**3 )) / 3     # cm3
    self.M = V * (self.nh + self.ne ) * _cgs.MUH2 * _cgs.MP # g = cm3 * g/cm3
    self.M /= _cgs.MSUN                            # Msun
    
    # to get the column density, integrate over radius r1 to r_10k
    #r_10k * 2 nh2_10k
    
    ################################################################
    #  print info. -> move to __str__ method
    #~ print ('M_10K   : {0:<7.2f} Msun\n'
            #~ 'R_10K   : {1:<7.0f} AU\n'
            #~ 'nH2_10K : {2:<7.1e} cm-3\n'
            #~ 'Y       : {3:<7.0f}\n'
            #~ 'T       : {4:<7.1f} K\n'.format(self.M.sum(),
                                        #~ self.r_10k/_cgs.AU,
                                        #~ self.nh2_10k,
                                        #~ self.Y,
                                        #~ self.temp_10k))
    #~ print 'Constraining the envelope : ', self.r_constraint
    #~ print ('nH2_r1000   : {0:<7.1e} cm-3\n'
            #~ 'T_r1000     : {1:7.1f} K\n'.format(self.nh2_r1000,
                                         #~ self.temp_r1000))
    
    print('printing input files')
    
    """
    id    : shell number
    ra,rb : inner & outer radius (m)
    za,zb : lower & upper height (m) (2D only)
    nh    : density (cm-3) of main collision partner (usually H2)
    nm    : density (cm-3) of molecule
    ne    : density (cm-3) of second collision partner (e.g. electrons)
    tk    : kinetic temperature (K) 
    td    : dust temperature (K)
    te    : electron/second coll. partner temperature (K)
    db    : 1/e half-width of line profile (Doppler b-parameter) (km s-1)
    vr    : radial velocity (km s-1)
    """
    # the model file!
    with open(_os.path.join(self.directory, self.modelfile),'w') as f:
        f.write('# Ratran input file based on Transphere results'+'\n')
        if self.skyonly: 
            f.write('# ... intended for (SKY) continuum calculations only.'+'\n')
        f.write("rmax={0:.5E}\n".format( self.rb[-1] ))       # rmax in METERS (convert from cm i.e. / 100)
        f.write("ncell={0:}\n".format(len(self.rb)))
        f.write("tcmb=2.735\n")
        f.write("columns=id,ra,rb,nh,nm,ne,tk,td,te,db,vr\n")
        f.write("gas:dust={0}\n".format(self.gas2dust))
        if self.skyonly: 
            f.write("kappa={0}\n".format(self.kappa))
        f.write('@\n')
        # r1/r2 in meter (convert from cm)
        for ii in range(0, len(self.rb)):
            test = ("{0:4} "                         #  1 id : shell number
                    "{1:12.5E} "                     #  2 ra : inner radius (m)  
                    "{2:12.5E} "                     #  3 rb : outer radius (m)
                    "{3:12.5E} "                     #  4 nh : density (cm-3) of main coll. partner (usually H2)
                    "{4:12.5E} "                     #  5 nm : density (cm-3) of molecule
                    "{5:12.5E} "                     #  6 ne : density (cm-3) of second coll. partner (e.g. electrons)
                    "{6:12.5E} "                     #  7 tk : kinetic temperature (K) 
                    "{7:12.5E} "                     #  8 td : dust temperature (K)
                    "{8:12.5E} "                     #  9 te : second coll. partner temperature (K)
                    "{9:12.5E} "                     # 10 db : 1/e half-width of line profile (Doppler b-parameter) (km s-1)
                    "{10:12.5E}\n")                  # 11 vr : radial velocity (km s-1)
            # now print the whole shebang
            f.write(test.format(ii + 1,              #  1 id : shell number
                self.ra[ii],                         #  2 ra : inner radius (m)  
                self.rb[ii],                         #  3 rb : outer radius (m)
                self.nh[ii],                         #  4 nh : density (cm-3) of main coll. partner (usually p-H2)
                self.nm[ii],                         #  5 nm : density (cm-3) of molecule
                self.ne[ii],                         #  6 ne : density (cm-3) of second coll. partner (e.g, e^-, o-H2)
                self.tk[ii],                         #  7 tk : kinetic temperature (K) 
                self.td[ii],                         #  8 td : dust temperature (K)
                self.te[ii],                         #  9 te : second coll. partner temperature (K)
                self.db[ii],                         # 10 db : 1/e half-width of line profile (Doppler b-parameter) (km s-1)
                round(self.vr[ii], 15))              # 11 vr : radial velocity (km s-1)
                    )          
                                
                                
                                
            #~ f.write("%4i %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E" % (ii+1, self.r1[ii]/100.0, self.r2[ii]/100.0, self.nh2int[ii], self.nh2int[ii]*self.abund_int[ii], self.tkint[ii], self.tdint[ii], self.db, self.vr[ii])+'\n')
    # the AMC.inp file
    if not self.skyonly:
        if self.molfile == '':
            sys.exit('Error: for AMC calculations the molecular datafile (molfile) needs to be set.')
        with open(_os.path.join(self.directory, "amc.inp"),'w') as f:
            f.write("source={0}\n".format(self.modelfile))
            f.write("outfile=populations.pop\n")
            f.write("molfile={0}\n".format(self.molfile))
            f.write("snr={0}\n".format(self.snr))
            f.write("velo=grid\n") 
            # velo=grid if velocity vector is given in input model?
            f.write("nphot={0}\n".format(self.nphot))
            f.write("kappa={0}\n".format(self.kappa))
            f.write("minpop={0:3.2E}\n".format(self.minpop))
            #~ f.write("seed=1971\n")       # NOT GOOD, avoid!
            f.write("fixset={0:3.2E}\n".format(self.fixset))
            f.write("go\n")
            f.write("q\n")
            f.write("\n")
    
    # the SKY.inp file
    with open(_os.path.join(self.directory, "sky.inp"),'w') as f:
        if self.skyonly:
            f.write("source={0}\n".format(self.modelfile))
        else:
            f.write("source=populations.pop\n")                     # just use the AMC output file (always set to populations.pop above)
        f.write("format={0}\n".format(self.outformat))
        f.write("outfile="+self.outputfile+"\n")
        f.write("trans={0}\n".format(self.trans))
        f.write("pix={0},{1:f},{2},{3}\n".format(self.imsize, self.pixel, self.pxlradius, self.los))
        if self.skyonly:
            f.write("chan=1,1.0\n")
        else:
            f.write("chan={0},{1:f}\n".format(self.chans, self.chwidth))
        f.write("distance={0}\n".format(self.dpc))
        f.write("units={0}\n".format(self.unit))
        f.write("go\n")
        f.write("q\n")
        f.write("\n")

def read_ratraninput(modelfile = "transphere.mdl"):
    """
    Read the input model to Ratran and plots all parameters 
    agains ra/r1 column.
    
        
    #Example model file, can be more lines before @
    
    # Ratran input file based on Transphere results
    rmax=1.757e+15
    ncell=51
    tcmb=2.728
    columns=id,ra,rb,nh,nm,tk,td,db,vr
    gas:dust=75.0
    @
       1 0.000e+00 3.501e+12 0.000e+00 0.000e+00 2.259e+02 2.259e+02 2.400e+00 0.000e+00
     ...
    """
    # imports...
    from scipy import array, where
    
    class Mdl: pass
    
    #with open(modelfile, 'r') as f:
    #    lines = f.read().split('\n')
    #table_start_indicator = where(array(lines) == '@')[0][0]
    ## get the header/preamble with the singel values
    #header = lines[:table_start_indicator]
    ## first get the comments in the header
    #Mdl.comments = [i for i in header if i[0] == "#"]
    ## second get the keyword - value pairs
    #values = [i for i in header if i[0] != "#"]
    #key_val = array([i.split('=') for i in values])
    ## add a checker function for input?
    ## put into class attributes
    #for key, val in zip(key_val[:,0], key_val[:,1]):
    #    if key in ['columns']:
    #        Mdl.__dict__[key] = val.split(',')
    #    elif key in ['ncell']:
    #        Mdl.__dict__[key.replace(':','2')] = int(val)
    #    else:
    #        Mdl.__dict__[key.replace(':','2')] = float(val)
    
    with open(modelfile, 'r') as f:
        line = f.readline()
        Mdl.comments = []
        while not line.startswith('@'):
            if line.startswith('#'):
                Mdl.comments.append(line)
                line = f.readline()
                pass
            else:
                keyval = line.strip('\n').split('=')
                try:
                    setattr(Mdl, keyval[0].replace(':','_'), float(keyval[1]))
                except(ValueError):
                    setattr(Mdl, keyval[0].replace(':','_'), keyval[1])
                line = f.readline()
        Mdl.columns = Mdl.columns.split(',')
        lines = f.readlines()
    
    # now get the whole data table
    #table = lines[table_start_indicator + 1:]
    table = lines
    Mdl.table = array([i.split() for i in table]).astype('float')
    Mdl.table = Mdl.table.transpose()
    # put each column in to its own attribute of the class
    for col,vals in zip(Mdl.columns, Mdl.table):
        Mdl.__dict__[col] = vals
    # convert to AUs
    Mdl.ra_au = Mdl.ra * 100 / _cgs.AU # input in m, so times 100 and devide by cm/AU
    Mdl.rb_au = Mdl.rb * 100 / _cgs.AU
    Mdl.r_au = (Mdl.ra_au + Mdl.rb_au) / 2.
    
    # convert to relative abundances
    Mdl.nh[0] = 1.0
    Mdl.ne[0] = 1.0
    Mdl.nm_rel = Mdl.nm / (Mdl.nh + Mdl.ne)
    Mdl.nh[0] = 0.0
    Mdl.ne[0] = 0.0
    Mdl.nm_rel[0] = 0.0
    return Mdl

# FIXME, does not work tries to import old adavis module
def plot_ratraninput(directory = '', modelfile = "transphere.mdl"):
    import matplotlib.pyplot as pl
    from scipy import arange, array
    pl.ion()
    #~ from adavis import set_rc
    from matplotlib import rc
    rc('axes', linewidth=1)
    rc('patch', linewidth=1.5)
    rc('lines', linewidth=1.5, markeredgewidth=1)
    # read in the model file
    Ratran_mdl = read_ratraninput(_os.path.join(directory, modelfile))
    
    # -2 because we dont want the first two cols, id and ra
    N = len(Ratran_mdl.columns) - 2
    pl.close()
    fig = pl.figure(1)#num = 1, figsize = ())
    plots = dict()
    for (i, dat) in zip(arange(N), Ratran_mdl.table[2:]):
        ax = 'ax{0}'.format(i)
        ylbl = Ratran_mdl.columns[i + 2]
        # now how many subplots do we need?
        if N % 3:
            np = N / 3 + 1 
        else:
            np = N / 3
        pln =  i +1
        plots[ax] = fig.add_subplot(np, 3, pln)
        if ylbl in ['db', 'vr']:
            plots[ax].semilogx(Ratran_mdl.ra, dat, '.')
        else:
            plots[ax].loglog(Ratran_mdl.ra, dat, '.')
        # set the rrights ylabel
        if ylbl in ['nh', 'nm', 'ne']:
            ylbl += ' (cm-3)'
        elif ylbl in ['tk', 'td', 'te']:
            ylbl += ' (K)'
        elif ylbl in ['db', 'vr']:
            ylbl += ' (km s-1)'
        elif ylbl in ['rb']:
            ylbl += ' (m)'
        plots[ax].set_ylabel(ylbl)
        plots[ax].grid()
    [plots['ax{0}'.format(i)].set_xlabel('ra (AU)') for i in arange(N-3, N)]
    #~ [plots['ax{0}'.format(i)].xaxis.set_visible(0) for i in arange(N-3)]
    [plots['ax{0}'.format(i)].xaxis.set_ticklabels([]) for i in arange(N-3)]


    #~ plots['ax{0}'.format(N-1)].set_xlabel('ra (AU)')
    fig.subplots_adjust(left=0.11, right= 0.97, bottom=0.11, top=0.96, wspace=0.43, hspace=0.15)
    return plots, fig
 
def create_molecular_abundance(temperature, 
                                abund_type = 'jump', 
                                Tjump = 100, 
                                Xs = [1E-4, 1E-9],
                                smooth = 0):
    """
    IMPORTANT:
    - assumes there is only one point where T > 100 K
    - and that this point has the highest cell number of all cells with
     T > 100 K
    """
    # calculates the molecular abundance from a predefined type and 
    # options
    # 'jump' abundance creates a jump from Xout to Xin 
    # where (temperature > Tjump)
    #
    # smooth : create jump that spans several 'smooth' numer of cells
    # -> should be an even number...
    #
    from scipy import where, ones, diff, sign, array, arange, exp
    [Xin, Xout] = Xs
    # jump abundance
    if 'jump' in abund_type:
        i100k = where(temperature >= Tjump)[0]
        mol_abundance = ones(len(temperature)) * Xout
        mol_abundance[i100k] = Xin
        if not smooth:
            print('discontinuity')
            # if it should be a discontinuity
            abund_param = False
        if smooth:
            # first make sure it is an even number of cellls to 
            # smooth over
            if smooth % 2 != 0:
                smooth += 1
                print('Argument \'smooth\' needs to be an even number',
                    'adding one.')
            # use sigmoid function to connect the discontinuities        
            # assumes that the abundance is consta chmnt before 
            # and after jump
            ijump = max(i100k)
            cells_x = ijump + array([-1, 1]) *  smooth/2
            Xs_diff = diff(Xs)                  # outer - inner
            d = direction = sign(Xs_diff)[0]    # direction of function
            B = height = abs(Xs_diff)[0]        # height of the sigmoid
            A = constant = max(Xs)              # curve start?-> min value
            x0 = center = ijump+0.5
            a = width = abs(diff(cells_x)[0])/10. # the width needs
                                # to be very small, should investigate
                                # and express more analytically
            sigmoid = lambda x: A + d * B / (1 + exp(-(x - x0) / a))
            #~ y_interp = splineinterpolate1d(x, y, k=2)
            #~ x = arange(cells_x[0],cells_x[1],1);
            #~ y = 1.0 / (1 + exp(-x / 0.1))
            
            #~ for i  in arange(x[0], x[1], 1):
                #~ mol_abundance[i] = y_interp(i)
            #~ mol_abundance[:x[0]] = Xin
            for i in arange(cells_x[0], cells_x[1],1):
                mol_abundance[i] = sigmoid(i)
            abund_param = dict(constant = A, direction = d, height = B, center = x0, width = a)
    # add more abundance types later on
    else:
        raise Exception('No abundance type given.')
    # send back the calculated abundance
    
    return mol_abundance, abund_param

# temporary function
# needs to be more modular
def plot_envstruct(self, mol_abundance = '', mark100k = True, **kawargs):
    if not hasattr(self, 'Envstruct'):
        raise Exception('you havent read in the transphere output')
    import matplotlib.pyplot as pl
    from matplotlib.ticker import ScalarFormatter, LogFormatter
    pl.ion()
    pl.close()
    fig = pl.figure(1, figsize=(8,6))
    ax1 = fig.add_subplot(111)
    pl.grid()
    ax2 = ax1.twinx()
    # Density
    p1 = ax1.loglog(self.Envstruct.r/_cgs.AU, self.n_h2, label='n_H2', **kawargs)
    ax1.set_xlabel('Radius (AU)')
    #~ ax1.set_xscale('log')
    ax1.set_ylabel('Number Density (cm-3)')
    # Temperature
    p2 = ax2.loglog(self.Envstruct.r/_cgs.AU, self.Envstruct.temp, color='r', label='Temp', **NiceLineSettings)

    ax2.yaxis.set_major_formatter(ScalarFormatter())
    
    ax2.set_ylabel('Temp (K)', color='r')
    if mol_abundance != '':
        def make_patch_spines_invisible(ax):
            ax.set_frame_on(True)
            ax.patch.set_visible(False)
            for sp in ax.spines.itervalues():
                sp.set_visible(False)
        ax3 = ax1.twinx()
        #~ ylims = ax3.get_ylim()
        #ax3.set_ylim(-0.05E-7, 1.85E-7)

        #~ p3 = ax3.loglog(self.Envstruct.r/_cgs.AU, mol_abundance, 'g')
        p3 = ax3.semilogx(self.Envstruct.r/_cgs.AU, mol_abundance, color='g', label='Mol Abund', **NiceLineSettings)
        ax3.spines["right"].set_position(("axes", 1.2))
        make_patch_spines_invisible(ax3)
        ax3.spines["right"].set_visible(True)
        ax3.set_ylabel('Rel. Abundance', color='g')
        #ax1.legend([p1, p2, p3], ['Density', 'Temp', 'Rel. Abundance'])
        #~ ax3.yaxis.set_major_formatter(())
        #~ ax3.xticks(['1E-9','1E-8','1E-7','1E-6','1E-5','1E-4'],[1E-9,1E-8,1E-7,1E-6,1E-5,1E-4])
        #~ ax3.set_yticks([1E-9,1E-8,1E-7,1E-6,1E-5,1E-4], minor=True)
        #~ ax3.tick_params(axis='y', direction='in')
        fig.subplots_adjust(right = 0.75)
    
    if mark100k:
        from scipy import where
        # where is the value closest to 100 K?
        i_100k = where(abs(100 - self.Envstruct.temp).round(2) == round(min(abs(100 - self.Envstruct.temp)), 2))[0][0]
        r_100k = self.Envstruct.r[i_100k]/_cgs.AU
        t_100k = self.Envstruct.temp[i_100k]
        ax2.annotate('T = {0:.1f} K\nR = {1:.1f} AU'.format(round(t_100k,2), r_100k),
                xy=(r_100k, t_100k), xycoords='data',
                xytext=(-30, -100), textcoords='offset points', fontsize=12,
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-.2"))
        #~ ax2.plot([r_100k], [t_100k] , 'o',color='r', ms=4, mew=0)
        #~ pl.legend('n_H2', 'Temp', 'Mol Abund')
    #~ else:
    ax1.xaxis.set_major_formatter(ScalarFormatter())
    if mol_abundance == '':
        #Create custom artists
        simArtist = pl.Line2D((0,1),(0,0), color='b')
        anyArtist = pl.Line2D((0,1),(0,0), color='r')
        
        #Create legend from custom artist/label lists
        ax1.legend([simArtist,anyArtist],
                  ['Density', 'Temperature'])
    elif mol_abundance != '':
        #Create custom artists
        simArtist = pl.Line2D((0,1),(0,0), color='b')
        anyArtist = pl.Line2D((0,1),(0,0), color='r')
        molArtist = pl.Line2D((0,1),(0,0), color='g')
        
        #Create legend from custom artist/label lists
        ax1.legend([simArtist, anyArtist, molArtist],
                  ['Density', 'Temperature', 'Mol. abundance'])
    
# test functions
def cleanup_ratran():
    import os
    filelist = ['populations.pop', 'amc.inp', 'sky.inp', 'transphere.mdl']
    os.system('rm -Rf _007')
    os.system('rm -Rf image.conv')
    for f in filelist:
        os.system('rm {0}'.format(f))

def cleanup_transphere():
    import os
    filelist = ['convhist.info', 'external_meanint.inp' , 'spectrum.dat', 'transphere.dat', 'dustopac_1.inp',  'envstruct.dat', 'starinfo.inp', 'transphere.inp', 'convhist.dat',  'dustopac.inp', 'envstruct.inp', 'starspectrum.inp']
    for f in filelist:
        os.system('rm {0}'.format(f))

# ratrarun obsolete. saved to be able to check later
"""
def run_ratran(r = 0.0, rho_dust = 0.0, temp = 0.0, db = 0.0, abund = 0.0, vr = 0.0, tdust = 0.0, dustonly = 0, mdl_file = 'transphere.mdl', dpc = 0.0, imsize = 129, pixel = 0.5, trans = '220.0e9', writeonly = 0, skyonly = 0, molfile='', ncell = 50, outputfile="ratranResult", snr=20, fixset=1e-6, minpop=1e-4, unit='Jypx', opstates=0, gas2dust=100, nphot=1000, temp_limit=10.0, rho_limit=1E-4, pxl_radius = 32, los = 2):

    #~ TODO!!!!: create a class, that saves everything,
    #~ envelope structure, model, etc before running RATRAN
    #~ (type : line or continuum)
    #~ TODO : validate the interpolation
    #~ TODO : check the input files for common errors(?)
    #~ TODO : calculate the column density


    from scipy import zeros, array, logspace, log10, pi, where
    import scipy.interpolate
    import sys
    import os
    import numpy as np
    from time import time

    # rewrite this part when changing to object oriented
    # now it is very hack-ish and non pythonic
        
    # Find the envelope cut off for T and n
    #
    # see if T goes below 10K (def) somewhere in the model
    try:
       ind_T = where(temp<temp_limit)[0].min()
    #~ ind = where(r<(1.2E4*cgs.AU))[0].max()
    except (ValueError): 
        ind_T = False
    # see if n goes below 1E-4 (def) somewhere in the model
    try:
        n_h2 = rho_dust * gas2dust  / _cgs.MUH2 / _cgs.MP
        ind_n = where((n_h2) < rho_limit)[0].min()
    except (ValueError):
        ind_n = False
    # T or n strongest constraints on radius
    # Neither T nor n constrain the radius
    # thus just use the last element
    if ind_n == False and ind_T == False:
        r_constraint = None
        ind = len(r)-1
    # Both constraint, which comes first
    elif ind_n != False and ind_T != False:
        # ind_n comes first
        ind = min((ind_n, int_T))
        # what if both have the same...
        # it will pick T, ok
        r_constraint = ['n', 'T'][ind_n < ind_T]
    elif ind_n != False:
        ind = ind_n
        r_constraint = 'n'
    elif ind_T != False:
        ind = ind_T
        r_constraint = 'T'
    
    r_10k = r[ind]
    print ind
    rho_dust_10k = rho_dust[ind]
    nh2_10k = rho_dust_10k * 100 / _cgs.MUH2 / _cgs.MP
    temp_10k = temp[ind]
    Y = r.max() / r.min()

    ind_r1000 = where(r > 1000 * _cgs.AU)[0].min()
    rho_dust_r1000 = rho_dust[ind_r1000]
    nh2_r1000 = rho_dust_r1000 * 100 / _cgs.MUH2 / _cgs.MP
    temp_r1000 = temp[ind_r1000]

    ###### cut off where T<10 K
    # first we have to remove all cells where T<10 K
    # RATRAN does not work well with them
    r = r[:ind]
    rho_dust = rho_dust[:ind]
    temp = temp[:ind]
    abund = abund[:ind]

    if tdust == 0.0:
        tdust = temp
    # from dust density to number density
    #    g/cm-3 * 100 / g
    # ! rho is dust density !
    nh = rho_dust * 100.0 / _cgs.MUH2 / _cgs.MP

    print 'Writing model in {0}'.format(mdl_file)
    #
    # for the refinement, easiest way is to redefine rx!
    #
    rx = logspace(log10(r[0]), log10(r[-1]), num=ncell+1, endpoint=True)
    rx = np.insert(rx, 0, 0)
    r1 = rx[0:-1]
    r2 = rx[1:]
    #~ r1=np.insert(r[0:-1],0,0)
    #~ r2=np.array(r)

    rr = zeros(ncell+1,float)
    rr[1:] = 10**( (np.log10(r1[1:]) + np.log10(r2[1:])) / 2.0 )
    rr[0] = rr[1]
    # Interpolate the values to 'ncell' cells
    nhf = scipy.interpolate.interp1d(log10(r), log10(nh))
    tkf = scipy.interpolate.interp1d(log10(r), log10(temp))
    tdf = scipy.interpolate.interp1d(log10(r), log10(tdust))
    abund_f = scipy.interpolate.interp1d(log10(r), log10(abund))
    #
    # Convert logarithms to floats
    nhint = 10**nhf(log10(rr))
    tkint = 10**tkf(log10(rr))
    tdint = 10**tdf(log10(rr))
    abund_int = 10**abund_f(np.log10(rr))
    ############################
    #~ nhint_p = nhint*2/4.
    #~ nhint_p[0] = 0.0
    #~ nhint_o = nhint*2/4.
    #~ nhint_o[0] = 0.0
    #~ teint = 10**tkf(log10(rr))
    ############################
    nhint[0] = 0.0

    # mass of it all
    #~ vol=[]
    #~ mass=[]
    # V = 4*pi*r**3/3
    # r in cm (?)
    V = 4 * pi * (r2**3-r1**3) / 3         # cm3
    M = V * nhint * _cgs.MUH2 * _cgs.MP # cm3 * g/cm3 ?
    M /= _cgs.MSUN
    
    # to get the column density, integrate over radius r1 to r_10k
    #r_10k * 2 nh2_10k

    print ('M_10K   : {0:<7.2f} Msun\n'
            'R_10K   : {1:<7.0f} AU\n'
            'nH2_10K : {2:<7.1e} cm-3\n'
            'Y       : {3:<7.0f}\n'
            'T       : {4:<7.1f} K\n'.format(M.sum(),
                                        r_10k/_cgs.AU,
                                        nh2_10k,
                                        Y,
                                        temp_10k))
    print 'Constraining the envelope : ', r_constraint
    print ('nH2_r1000   : {0:<7.1e} cm-3\n'
            'T_r1000     : {1:7.1f} K\n'.format(nh2_r1000,
                                         temp_r1000))
    #~ raise Exception('Test')
    #~ return r2,r1
    #
    # INPUT MODEL
    #
    f = open(mdl_file,'w')
    f.write('# Ratran input file based on Transphere results'+'\n')
    if skyonly == 1:
        f.write('# intended for (SKY) continuum calculations only.'+'\n')
    f.write('rmax={0:-5.3e}\n'.format( max(r2) / 1.0e2 ))
    f.write('ncell=%i'  % (len(r2))+'\n')
    f.write('tcmb=2.728\n')
    #~ f.write('columns=id,ra,rb,nh,nm,tk,td,db,vr\n')
    f.write('columns=id,ra,rb,nh,nm,tk,td,db,vr\n')
    f.write('gas:dust={0}\n'.format(gas2dust))
    if skyonly == 1:
        f.write('kappa=jena,thin,e6\n')
    f.write('@\n')
    for ii in range(0,len(r1)):
        f.write("%4i %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e" % (ii+1, r1[ii]/100.0, r2[ii]/100.0, nhint[ii], nhint[ii]*abund_int[ii], tkint[ii], tdint[ii], db, vr)+'\n')
    #~ for ii in range(0,len(r1)):
        #~ f.write("%4i %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E" % (ii+1, r1[ii]/100.0, r2[ii]/100.0, nhint_p[ii], nhint[ii]*abund_int[ii], nhint_o[ii], tkint[ii], tdint[ii], teint[ii], db, vr)+'\n')
    f.close()
    #
    # AMC input file
    #
    if skyonly == 0:
        if molfile == '': sys.exit('Error: for AMC calculations the molecular datafile (molfile) needs to be set.')
        f = open("amc.inp",'w')
        f.write("source="+mdl_file+'\n')
        f.write("outfile=populations.pop\n")
        f.write("molfile="+molfile+'\n')
        f.write("snr={0}\n".format(snr))
        f.write("nphot={0}\n".format(nphot))
        f.write("kappa=jena,thin,e6\n")
        f.write("minpop={0}\n".format(minpop))
        f.write("seed=1971\n")
        f.write("fixset={0}\n".format(fixset))
        f.write("go\n")
        f.write("q\n")
        f.write(" \n")
        f.close()
    #
    # SKY input file
    #
    f = open("sky.inp",'w')
    if skyonly == 0:
        f.write("source=populations.pop\n")
    else:
        f.write("source="+mdl_file+"\n")
    f.write("format=miriad\n")
    f.write("outfile="+outputfile+"\n")
    f.write("trans="+trans+"\n")
    #~ f.write("pix="+str(imsize)+","+str(pixel)+",32,8\n")
    f.write("pix={0},{1},{2},{3}\n".format(imsize, pixel, pxl_radius, los))
    if skyonly == 0:
        f.write("chan=100,0.2\n")
    else:
        f.write("chan=1,1.0\n")
    f.write("distance="+str(dpc)+"\n")
    f.write("units={0}\n".format(unit))
    f.write("go\n")
    f.write("q\n")
    f.close()
    #
    #
    #
    if writeonly == 0:
        if skyonly == 0:
            #~ print "Starting AMC calculation..."
            t1 = time()
            os.system('amc amc.inp')
            print 'AMC took :  {0:2.2f} hours'.format((time()-t1)/3600.0)
            os.system('alert \"AMC has finished running.\"')
        #~ print "Starting SKY calculation..."
        t1 = time()
        os.system("sky sky.inp")
        print ' SKY took : {0:2.2f} seconds'.format(time()-t1)
        os.system('alert \"SKY has finished running.\"')
"""

def save_ratran(Obj, filename = 'ratranmodel.pickle'):
    # take input object
    # and save it as a dictionary to filename in directory
    import pickle
    # take all the original input parameters and put into a dictionary
    inp = vars(Obj.Input)
    # now a dictionary is a non dynamic structure
    # as opposed to a dynamically created object attribute
    # e.g., with the __dict__ method (-> doesn't work with pickle)
    with ChangeDirectory(Obj.directory):
        with open(filename, 'w') as f:
            pickle.dump(inp, f)

def load_ratran(directory = '', filename = 'ratranmodel.pickle'):
    # load input object from filename in directory
    # and create the ratran object
    import pickle
    with open(_os.path.join(directory, filename), 'r') as f:
        inputdict = pickle.load(f)
    Obj = Ratran(**inputdict)
    # IDEA : add so that it loads the output(?) as well?
    return Obj

def save_transphere(Obj, filename = 'transpheremodel.pickle'):
    # take input object
    # and save it as a dictionary to filename in directory
    import pickle
    # take all the original input parameters and put into a dictionary
    inp = vars(Obj.Input)
    # now, a dictionary is a non dynamic structure
    # as opposed to a dynamically created object attribute
    # e.g., with the __dict__ method (-> doesn't work with pickle)
    with ChangeDirectory(Obj.directory):
        with open(filename, 'w') as f:
            pickle.dump(inp, f)

def load_transphere(directory = '', filename = 'transpheremodel.pickle'):
    # load input object from filename in directory
    # and create the transphere object
    import pickle
    with open(_os.path.join(directory, filename), 'r') as f:
        inputdict = pickle.load(f)
    Obj = Transphere(**inputdict)
    # IDEA : add so that it loads the output(?) as well?
    return Obj

def temp_pop(n1, n2, g1, g2, nu):
    """ 
    Calculate the excitation temperature of transition from 2 to 1
    """
    
    numer = _cgs.HH * nu
    denom = _cgs.KK * _scipy.log(n1 * g2 / (n2 * g1))
    return numer / denom

def calc_nu(N, Q, gu, Eu, T):
    """
    Calculate the column density of the upper energy level
    assuming LTE 
    """
    part1 = N / Q * gu
    part1 = _scipy.exp(-Eu / (_cgs.KK * T))
    return part1 * part2
    
def calc_tau(Nu, Aul, freq, width, Eu, T):
    """
    Calculate the optical depth (tau) for given parameters
    from Goldsmith and Langer (1999)
    """
    part1 = Nu * _cgs.CC**3 * Aul / (8 * _scipy.pi * freq**3 * width) 
    part2 = (_scipy.exp(Eu/T) - 1)
    return part1 * part2

# END OLD HELP FUNCTIONS
########################################################################


class Input(object):pass


class Output(object):
    def __init__(self, popfile = '', directory = '', molfile = ''):
        """
        This class should take a directory and ratran output file
        and just read in whatever it can. Should have reasonable
        assumption.
        
        popfile : just name of popfile
        directory : directory where the Ratran output is located
        molfile : path to the moldata file (if different from that 
                    of popfile(s))
        """
        
        ### Parse directory input
        self._check_directory(directory)
        ##### The popfile is mandatory output of RATRAN so it is 
        ##### kind of needed for this class to be of any use...
        ### Parse popfile input
        self._check_popfile(popfile)
        ### create the full path to the popfile
        self.path_to_popfile = _os.path.join(self.directory, self.popfile)
        ### Read the population output from Ratran
        self._read_population_output()
        
        ##### If there are no history files present, don't load anything.
        ### Check if history files present in directory
        if self._check_history_files(): # if so, read them
            self._read_population_history()
        else:
            pass
        ### Check molfile,
        if self._check_molfile(molfile):
            self.Moldata  = _moldata.Read(self.molfile)
    
    #### checking functions
    def _check_directory(self, directory):
        if directory:   # if directory input
            self.directory = directory
        else:
            self.directory = _os.getcwd()
        
    def _check_popfile(self, popfile):
        if popfile:
            self.popfile = popfile
        else:
            try:
                self.popfile = [i for i in _os.listdir(self.directory) if i[-3:] == 'pop'][0]
            except (IndexError):
                raise Exception ('No popfile given or found, or bad directory input.')
    
    def _check_molfile(self, molfile):
        if molfile: 
            # if moldata file input
            if not _os.path.isfile(molfile):
                print('Molfile does not exists/wrong path.')
                return False
            #~ self.Moldata  = _moldata.Read(molfile)
            self.molfile = molfile
            return True
        elif self.Pop.__dict__.has_key('molfile'): 
            # moldata file in popfile
            print self.Pop.molfile
            if not _os.path.isfile(self.Pop.molfile):
                print('Molfile does not exists/wrong path.')
                return False
            self.molfile = self.Pop.molfile
            return True
        else: # moldata not input, nor in Pop
            # NO moldata read, no Moldata nested class exists
            # raise warning?
            print('Molfile neither input nor in population output.')
            return False
        
    def _check_history_files(self):
        hisfiles = [i for i in _os.listdir(self.directory) if i[-3:] == 'his']
        if not hisfiles:
            return False
        else:
            return True
    
    #### reading functions    
    def _read_population_output(self):
        from scipy import arange, array
        #~ if not directory:   # if no directory was given
            #~ directory = _os.getcwd()
        
        class Pop(object): pass
                
        with open(self.path_to_popfile) as f:
            line = f.readline()
            Pop.comments = []
            while not line.startswith('@'):
                if line.startswith('#'):
                    Pop.comments.append(line)
                    line = f.readline()
                    pass
                else:
                    keyval = line.strip('\n').split('=')
                    try:
                        setattr(Pop, keyval[0].replace(':','_'), float(keyval[1]))
                    except(ValueError):
                        setattr(Pop, keyval[0].replace(':','_'), keyval[1])
                    line = f.readline()
            Pop.columns = Pop.columns.split(',')
            lines = f.readlines()
            
        lines = array([i.strip().split() for i in lines], dtype='float')
        lines = lines.transpose()
        for colname, i in zip(Pop.columns, arange(len(Pop.columns))):
            if colname == 'lp':
                # the lp are all the columns left,
                # the number of columns amount to the numbers of 
                # levels of that particulat molecule (see moldata file)
                setattr(Pop, colname, lines[i:]) 
            else:
                setattr(Pop, colname, lines[i])
        self.r = (Pop.rb + Pop.ra) / 2.
        #[setattr(self, col, i) in zip(self.columns,arange(len(self.columns)-1))]
        self.Pop = Pop
        
    def _read_population_history(self):
        hisfiles = [i for i in _os.listdir(self.directory) if i[-3:] == 'his']
        self.hisfiles = sorted(hisfiles)
        class Pop_his(object): pass
        Pop_his.hisfiles = self.hisfiles
        # copy pasted below, check and refine
        from scipy import loadtxt, zeros, array
        tables = []
        onetime = False
        for f in self.hisfiles:    
            with open(_os.path.join(self.directory, f)) as text:
                lines = text.readlines()
                fname = text.name
            print('Reading file : {0}'.format(f))
            rawhdr = [c for c in lines if len(c)<100 and not c.startswith('@')]
            preamble = [i.strip().split('=') for i in rawhdr if not i.startswith('#')]
            while not onetime:
                for keyval in preamble:
                    try:
                        setattr(Pop_his, keyval[0].replace(':','_'), float(keyval[1]))
                    except(ValueError):
                        setattr(Pop_his, keyval[0].replace(':','_'), keyval[1])
                onetime = True
                Pop_his.preamble = preamble
                #~ self.molfile = rawhdr[]
            dat = [c for c in lines if len(c)>100 and not c.startswith('@')]
            # 'raw' table data
            coldata = loadtxt(dat).transpose()
            popfile = dict()
            popfile.__setitem__('filename', str(fname))
            # get the individual settings written in the header/preamble
            [popfile.__setitem__(i[0], i[1]) for i in preamble]
            # from the settings get the columns present
            # the script assumes that the last column is level populations
            columns = [i for i in popfile['columns'].split(',')]
            [popfile.__setitem__(i, j) for (i,j) in zip(columns[:-1], coldata[:len(columns)-1])]
            # convergence percentage
            convergence = [i for i in rawhdr if i.startswith('#') and 'converged' in i]
            convergence = convergence[0].strip().strip('#AMC:').strip('converged').strip('\% ')
            popfile.__setitem__('convergence', float(convergence))
            # get the level populations table
            popfile.__setitem__('lp', coldata[len(columns)-1:])
            # write the raw data and header to the dictionary as well
            popfile.update(dict(header = rawhdr, rawdata = coldata))   
            # append to list of population file info
            tables.append(popfile)
        Pop_his.pop_tables = tables
        self.Pop_his = Pop_his

    #### plotting functions
    ## for final populations
    def plot_structure(self):
        pass
    
    def plot_populations(self, levels = [], runjump = 10, leveljump = 10):
        pass
    
    def plot_tau(self, trans = [12, 10], width = 1.0E4):
        # ra, rb in m, convert to cm
        # width in cm/s
        itup, itdown = trans[0] - 1, trans[1] - 1
        try:
            transIndex = [i['trans'] for i in self.Moldata.radtrans if i['up'] == 12 and i['down'] == 10][0]
        except (IndexError):
            print('No such transition for molecule : {0}'.format(self.Moldata.molecule))
            print('Check your moldata file at : {0}'.format(self.Moldata.molfile))
            return False
        transIndex -= 1
        
        # First, get the opacity from the model outpupt
        class Tau: pass
        self.Tau = Tau
        # the density of the upper energy level
        # from the model data
        self.Tau.Nu = Nu = (self.Pop.rb - self.Pop.ra) * 100.0 * self.Pop.nm * self.Pop.lp[transIndex]
        self.Tau.Au = Aul = self.Moldata.radtrans[transIndex]['aul']
        self.Tau.freq = freq = self.Moldata.radtrans[transIndex]['freq']
        self.Tau.Eu = Eu = self.Moldata.radtrans[transIndex]['eu']
        self.Tau.T = T = self.Pop.tk
        # now calculate the opacity
        self.Tau.tau = tau = calc_tau(Nu, Aul, freq, width, Eu, T)
        
        radii = (self.Pop.ra + self.Pop.rb)/2. * 100 / _cgs.AU
        
        # plotting
        _pl.ion()
        #~ fig = _pl.figure(num=1, figsize=(3.5,3))
        fig = _pl.figure(num=1)
        ax = fig.add_subplot(111)
        ax.loglog(radii , tau, label=r' - '.join([self.Moldata.get_lvl(i, tex = 1) for i in trans]), color=_gc(1).next(), lw=2, marker='o', ms=3, mew=0)

        
        
        # Second, the calculate the opacity if it is pure LTE conditions
        
        try:
            # first initialise the get_partition method
            # to get the partition function values
            # at all temperatures
            self.Moldata.get_partition()
            plot_lte = True
        except:
            print ('cannot get the partition function values.')
            plot_lte = False
            
        if plot_lte:
            # Nu  = nm * (rb - ra) * gu / Qrot(T) * exp(-Eu/T)
            nm = self.Pop.nm
            ra = self.Pop.ra * 100 # in cm
            rb = self.Pop.rb * 100 # in cm
            gu = self.Moldata.elev[itup]['weight']
            qrot = self.Moldata.qrot(T)
            self.Tau.Nu_lte = Nu_lte = nm * (rb - ra) * gu / qrot * _scipy.exp(-Eu/T)
            self.Tau.tau_lte = tau_lte = calc_tau(Nu_lte, Aul, freq, width, Eu, T)
            _pl.loglog(radii, tau_lte, label=r' LTE', color='#008822', lw=2, marker='o', ms=3, mew=0)
        ax.set_xlabel('Radius [AU]')
        ax.set_ylabel(r'$\tau$')
        ax.legend()
        ax.grid()
    
    def plot_populations(self, levels = [], runjump = 10, leveljump = 10):
        """ 
        Function to plot the populations level for 
        every 'runjump' iteration and every 'leveljump' population 
        level
        """
        from scipy import linspace, log10, logspace, array
        #~ import matplotlib.pyplot as pl; 
        _pl.ion()
        from matplotlib import cm
        #~ from adapy.libs import cgsconst as cgs

        radii = (self.Pop_his.pop_tables[0]['ra'] + self.Pop_his.pop_tables[0]['rb'])/2. * 100 / _cgs.AU
        tables_select = self.Pop_his.pop_tables[::10]
        
        
        # plot with import colorsys, and then colorsys.hsl_to_rgb
        # h : hue = color (0 = red, 120 = green, 240 = blue)
        # s : saturation = hardness of color (0 = grey, 1 = full color)
        # l : light = brightness (0 = black, 1 = white, 0.5 = full color)
        # hsv = hue, saturation, value(brightness)
        
        
        if not levels:
            lp = [i['lp'][::int(leveljump)] for i in self.Pop_his.pop_tables[::int(runjump)]] # every 'leveljump' level of every 'runjump' run
            # plot all iterations
            [[_pl.loglog(radii , i, color=str(c), lw=1, ls=':', marker='o', ms=3, mew=0) for i in j] for (j,c) in zip(lp, linspace(1, 0.4, len(lp)))]
            # plot the last iteration
            #~ [pl.loglog(radii , i, color='k', lw=1.5, marker='o', ms=3, mew=0) for i in self.Pop_his.pop_tables[-1]['lp'][::int(leveljump)]]
        if levels:
            lp = array([[i['lp'][j-1] for j in levels] for i in self.Pop_his.pop_tables[::int(runjump)]])
            # plot all iterations
            [[_pl.loglog(radii , i, color=str(c), lw=1, ls=':', marker='o', ms=3, mew=0) for i in j] for (j,c) in zip(lp, linspace(0.7, 0.2, len(lp)))]
            # plot the last iteration
            #~ [pl.loglog(radii , i, color='k', lw=1.5, marker='o', ms=3, mew=0) for i in self.Pop_his.pop_tables[-1]['lp'][array(levels)-1]]
        # allways plot last history item, in green and the others in dashed grey scale
        #~ pl.loglog(radii , self.pop_tables[-1]['lp'][::int(leveljump)].transpose(),color='g', lw=1, marker='o', ms=2, mew=0)
        # should plot the resulting 
        # a bit to complicated list comprehension... reformat to better
        # programming style...
        [_pl.loglog(radii , self.Pop.lp[i], label=self.Moldata.get_lvl(i+1, tex = 1), color=c, lw=1.5, marker='o', ms=3, mew=0) for (i,c) in zip(array(levels)-1, _gc(len(levels)))]
        _pl.legend()
        _pl.grid()
        
    def plot_tex(self, trans = [12, 10], runjump = 10, history = True):    # specify which transition
        """
        Plots the excitation temperature for every 'runjump' iteration 
        and every 'leveljump' population level
        """
        from scipy import linspace, array
        #~ import matplotlib.pyplot as pl; 
        _pl.ion()
        #~ from matplotlib import cm
        #~ from adapy.libs import cgsconst as cgs
        
        radii = (self.Pop.ra + self.Pop.rb)/2. * 100 / _cgs.AU
        if hasattr(self, 'Pop_his') and history:
            ### get the populations for the levels from the different runs 
            runs = array([i for i in self.Pop_his.pop_tables[::int(runjump)]]) # every 10th run
            lp = [[j['lp'][i-1] for i in trans] for j in runs]
        
        ### 
        gweights = [self.Moldata.elev[trans[0]-1]['weight'], self.Moldata.elev[trans[1]-1]['weight']] 
        print('   Molecular weights : upper {0[0]}  lower {0[1]}'.format(gweights))
        nu = _cgs.CC*abs(self.Moldata.elev[trans[0]-1]['energies'] - self.Moldata.elev[trans[1]-1]['energies'])
        print('   Frequency : {0:3.3f} E09 GHz'.format(nu*1e-9))
        trans_str = [self.Moldata.get_lvl(i) for i in trans]
        #trans_str = [self.elev[trans[0] ,.elev[trans[1]-1]['j']]
        print('   Transition : {0[0]} - {0[1]}'.format(trans_str))
        if hasattr(self, 'Pop_his') and history:
            tex = []
            for levels in lp:
                t = temp_pop(levels[1], levels[0], gweights[1], gweights[0], nu)
                tex.append(t)
            [_pl.loglog(radii, j, color=str(c), lw=1, marker='o', ms=4, mew=0)
                for (j, c) in zip(tex, linspace(0.7, 0, len(tex)))]
        #
        #
        ### get the final tex curve
        tex_final = temp_pop(self.Pop.lp[trans[1]-1], self.Pop.lp[trans[0]-1],
                             gweights[1], gweights[0], nu)
        str_transitions = [self.Moldata.get_lvl(i, tex=1) for i in trans]
        _pl.loglog(radii, tex_final, label=r' - '.join(str_transitions),
                   color=_gc(1).next(), lw=2, marker='o', ms=5, mew=0)
        #
        _pl.loglog([radii[-1], radii[-1]], [0, 1], )
        #
        _pl.legend()
        _pl.grid()
        # tex for all transitions? why do that...
        # perhaps be able to supply several pairs of
        # transitions but not more
        #~ tex = [[temp_pop(ldown, lup, gweights[1], gweights[0], nu) for
            #~ (ldown, lup) in zip(run[1], run[0])] for (run) in lp]
        #~ return tex, lp
        #~ for runlp in tex:
            #~ for

    def plot_radiation(self):
        pass
