#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pyrt.py
#
#  Module to interface with RADEX, RATRAN and Transphere
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
#  version  0.1a
#
#
doc="""
PyRT - Python Radiative Transfer

Module to interface with various radiative transfer tools.

Transphere : Dust continuum radiative transfer

Radex : Point/escape probability

Ratran : 1D radiative transfer code


Dependencies:

cgsconst.py - relevant cgs constants

"""
#----[ CHANGE LOG ]----
"""

* 2012 Nov 26 
    "Transphere.run_transphere" only prints the progress, and not the
    whole convergence stuff, it stores that in a parameter called
    "Transphere.transphere_output"

* 2012 Nov 25 
    Now works against a directory if it is given, and read the results 
    from it.

"""
#------------------------------------------------------------------------
# TODO:
"""
TODO :
        Radex 
            - Fix links to libraries
            - test run with :
                ra_inp = dict(
                f=203.4, 
                tkin=124, 
                h2dens=array([1e-1,1,10])*1e8, 
                oflux=7, 
                lwidth=6, 
                molf='p-h218o.dat')
            - make it run in parallell when arrays are input
        Transphere
            - Prepare input and run Transphere
            - Optimize the writeTransphereInput

"""
import cgsconst as _cgs
import os as _os
import sys as _sys
import subprocess as _subprocess

########################################################################
# GENERAL HELP FUNCTIONS (move to adavis_core)
def check_input(input_dictionary, input_defaults):
    return 0

def make_dirs(path):
    import os
    import errno
    # the makedirs function will raise a EEXIST error 
    # if the directory already exists, if so it returns False
    # if some other error is raise, it will raise that 
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            return False
    return True



class ChangeDirectory:
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = _os.getcwd()
        _os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        _os.chdir(self.savedPath)

# Now you can enter the directory like this:

#~ import subprocess
#~ 
#~ with cd("~/Library"):
   #~ # we are in ~/Library
   #~ subprocess.call("ls")






########################################################################
# MODELING HELP FUNCTIONS
def bplanck(nu,T): # Returns the Spectral Radiance (Planck)
    from scipy import exp, constants
    try:
        x = _cgs.HH*nu/(_cgs.KK*T)
        bpl = (2.0*_cgs.HH*nu**3/_cgs.CC**2)/(exp(x)-1.0)
        return bpl
    except (ZeroDivisionError):
        return 0

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
    from ..views import set_rc
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

def read_ratraninput(model_file = "transphere.mdl"):
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
    
    with open(model_file) as f:
        lines = f.read().split('\n')
    table_start_indicator = where(array(lines) == '@')[0][0]
    # get the header/preamble with the singel values
    header = lines[:table_start_indicator]
    # first get the comments in the header
    Mdl.comments = [i for i in header if i[0] == "#"]
    # second get the keyword - value pairs
    values = [i for i in header if i[0] != "#"]
    key_val = array([i.split('=') for i in values])
    # add a checker function for input?
    # put into class attributes
    for key, val in zip(key_val[:,0], key_val[:,1]):
        if key in ['columns']:
            Mdl.__dict__[key] = val.split(',')
        elif key in ['ncell']:
            Mdl.__dict__[key.replace(':','2')] = int(val)
        else:
            Mdl.__dict__[key.replace(':','2')] = float(val)
    # now get the whole data table
    table = lines[table_start_indicator + 1:]
    Mdl.table = array([i.split() for i in table][:-1]).astype('float')
    Mdl.table = Mdl.table.transpose()
    # put each column in to its own attribute of the class
    for col,vals in zip(Mdl.columns,Mdl.table):
        Mdl.__dict__[col] = vals
    # convert to AUs
    Mdl.ra /= 100*_cgs.AU
    Mdl.rb /= 100*_cgs.AU
    return Mdl

# FIXME, does not work tries to import old adavis module
def plot_ratraninput(model_file = "transphere.mdl"):
    import matplotlib.pyplot as pl
    from scipy import arange, array
    pl.ion()
    #~ from adavis import set_rc
    from matplotlib import rc
    rc('axes', linewidth=1)
    rc('patch', linewidth=1.5)
    rc('lines', linewidth=1.5, markeredgewidth=1)
    # read in the model file
    Ratran_mdl = read_ratraninput(model_file)
    
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
            plots[ax].semilogx(Ratran_mdl.ra, dat)
        else:
            plots[ax].loglog(Ratran_mdl.ra, dat)
        # set the rrights ylabel
        if ylbl in ['nh', 'nm', 'ne']:
            ylbl += ' (cm-3)'
        elif ylbl in ['tk', 'td', 'te']:
            ylbl += ' (K)'
        elif ylbl in ['db', 'vr']:
            ylbl += ' (km s-1)'
        elif ylbl in ['rb']:
            ylbl += ' (AU)'
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
                                Xs = [1E-4, 1E-9]):
    # calculates the molecular abundance from a predefined type and 
    # options
    # 'jump' abundance creates a jump from Xout to Xin 
    # where (temperature > Tjump)
    #
    from scipy import where, ones
    
    [Xin, Xout] = Xs
    # jump abundance
    if 'jump' in abund_type:
        i = where(temperature > Tjump)[0]
        mol_abundance = ones(len(temperature)) * Xout
        mol_abundance[i] = Xin
    # add more abundance types later on
    else:
        raise Exception('No abundance type given.')
    # send back the calculated abundance
    return mol_abundance
# temporary function
# needs to be more modular
def plot_envstruct(self, mol_abundance=0):
    if not hasattr(self, 'Envstruct'):
        raise Exception('you havent read in the transphere output')
    import matplotlib.pyplot as pl
    pl.ion()
    pl.close()
    fig = pl.figure(1, figsize=(8,6))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    # Density
    p1 = ax1.loglog(self.Envstruct.r/_cgs.AU, self.n_h2, label='n_H2')
    ax1.set_xlabel('Radius (AU)')
    #~ ax1.set_xscale('log')
    ax1.set_ylabel('Number Density (cm-3)')
    # Temperature
    p2 = ax2.loglog(self.Envstruct.r/_cgs.AU, self.Envstruct.temp, color='r', label='Temp')
    ax2.set_ylabel('Temp (K)', color='r')
    if mol_abundance != 0:
        def make_patch_spines_invisible(ax):
            ax.set_frame_on(True)
            ax.patch.set_visible(False)
            for sp in ax.spines.itervalues():
                sp.set_visible(False)
        ax3 = ax1.twinx()
        #~ ylims = ax3.get_ylim()
        #ax3.set_ylim(-0.05E-7, 1.85E-7)

        #~ p3 = ax3.loglog(self.Envstruct.r/_cgs.AU, mol_abundance, 'g')
        p3 = ax3.semilogx(self.Envstruct.r/_cgs.AU, mol_abundance, color='g', label='Mol Abund')
        ax3.spines["right"].set_position(("axes", 1.2))
        make_patch_spines_invisible(ax3)
        ax3.spines["right"].set_visible(True)
        ax3.set_ylabel('Rel. Abundance', color='g')
        #ax1.legend([p1, p2, p3], ['Density', 'Temp', 'Rel. Abundance'])
        fig.subplots_adjust(right = 0.75)
    pl.legend()

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

# obsolete! yay, finally
def run_ratran(r = 0.0, rho_dust = 0.0, temp = 0.0, db = 0.0, abund = 0.0, vr = 0.0, tdust = 0.0, dustonly = 0, mdl_file = 'transphere.mdl', dpc = 0.0, imsize = 129, pixel = 0.5, trans = '220.0e9', writeonly = 0, skyonly = 0, molfile='', ncell = 50, outputfile="ratranResult", snr=20, fixset=1e-6, minpop=1e-4, unit='Jypx', opstates=0, gas2dust=100, nphot=1000, temp_limit=10.0, rho_limit=1E-4, pxl_radius = 32, los = 2):
    """

    TODO!!!!: create a class, that saves everything,
    envelope structure, model, etc before running RATRAN
    (type : line or continuum)
    TODO : validate the interpolation
    TODO : check the input files for common errors(?)
    TODO : calculate the column density

    """

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


def save_ratran(Obj, filename = 'ratranmodel.pickle'):
    # take input object
    # and save it as a dictionary to filename in directory
    import pickle
    # take all the original input parameters and put into a dictionary
    inp = vars(Obj.InputParameters)
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
    inp = vars(Obj.InputParameters)
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
######################################################################
## # RADIATIVE TRANSFER / MODELING

#### NOT used
class Radex:
    """
    Input (* = mandatory input)
    *f      : Frequency - frequency of line (in GHz)
    df      : Uncertainty of line position (defaults to 20 MHz)
    *tkin   : Kinetic temperature  - one value or array (Kelvin)
    *h2dens : H2 density one value or array (cm-3)
    tbg     : Background temperature - defaults to 2.73 K
    *oflux  : Observed flux in K km/s, to compare
    *lwidth : Line width (FWHM)
    *molf   : Molecular data file (LAMBDA format)
    ftol    : How close can the flux be for it to be the same? (0.1)
    maxiter : Maximum number of iterations (100)
    eps     : Below what is zero or negat ive line intensity (1E-20)
    cdens   : First guess molecular line density (1e12)
    debug   : Output debug messages? (Boolean)
    silent  : Output any msgs at all? (Boolean)
    """
    def __init__(self, **kwargs):
        """
        The program will take arrays as input.
        So fx if you initialize it with
        Radex(f = 100, tkin = arange(10), h2dens = 1e9, cdens=)

        """
        from scipy import array, arange
        #~ from ..adacore import ParError as _ParError
        #~ from multiprocessing import Pool
        # check input parameters and store in class
        # and raise ParError if wrong/missing
        params = array(['silent', False, 'f', None, 'df', 0.02,
                        'h2dens', None, 'tbg', 2.73, 'oflux', None,
                        'lwidth', None, 'molf', None, 'tol', 0.1,
                        'maxiter', 100, 'eps', 1.0e-20, 'cdinit', 1E12,
                        'debug', False, 'tkin', None])
        for par,stdval in zip(params[0::2],params[1::2]):
            if par in kwargs:
                if par == 'molf':
                    self.molf = kwargs['molf']
                elif par == 'tkin':
                    if type(kwargs[par]) in [type([]), type(arange(1))]:
                        self.tovary = 'tkin'
                        self.tkin = kwargs[par]
                    else:
                        self.tkin = kwargs[par]
                elif par == 'h2dens':
                    if type(kwargs[par]) in [type([]), type(arange(1))]:
                        if hasattr(self, 'tovary'):
                            raise Exception('\'tkin\' and \'h2dens\' ',
                                        'are both arrays, only one can',
                                        ' be varied at a time')
                        self.tovary = 'h2dens'
                        self.h2dens = kwargs[par]
                    else:
                        self.h2dens = kwargs[par]
                else:
                    exec('self.{0} = {1}'.format(par, kwargs[par]))
            else:
                if stdval != None:
                    exec('self.{0} = {1}'.format(par, stdval))
                    if not self.silent:
                        print('Using default value '
                        'of {0} for {1}'.format(stdval, par))
                else:
                    raise Exception(par)
    def __str__(self):
        out_string = (40*'='+'\n'+
                        'Frequency  : {0} GHz\n'.format(self.f)+
                        'Bandwidth  : {0} MHz\n'.format(self.df*1E3)+
                        'Kinetic T  : {0} K\n'.format(self.tkin)+
                        'H2 density : {0} cm-3\n'.format(self.h2dens)+
                        'Tbg        : {0} K\n'.format(self.tbg)+
                        'Obs. Flux  : {0} K km s-1\n'.format(self.oflux)+
                        'Line width : {0} km s-1\n'.format(self.lwidth)+
                        'Moldata    : {0} \n'.format(self.molf)+
                        'Tolerance  : {0} \n'.format(self.tol)+
                        'Max. inter : {0} \n'.format(self.maxiter)+
                        'Eps        : {0} \n'.format(self.eps)+
                        'Init Cdens : {0} \n'.format(self.cdinit)+
                        'Debug      : {0}'.format(self.debug))
        return out_string
    def _create_input(self):
        with open('radex.inp','w') as f:
            f.write('{0}\n'.format(self.molf))
            f.write('radex.out\n')
            f.write('{0} {1}\n'.format((self.f-self.df),
                                        (self.f+self.df)))
            f.write('{0}\n'.format(self.tkin))
            f.write('1\n')
            f.write('H2\n')
            f.write('{0:3.2E}\n'.format(self.h2dens))
            f.write('{0}\n'.format(self.tbg))
            f.write('{0:3.2E}\n'.format(self.cdinit))
            f.write('{0}\n'.format(self.lwidth))
            f.write('0\n')
    def _read_output(self):
        import sys
        with open('radex.out', 'r') as f:
            lines = f.readlines()
        if (lines[-2].split()[-1] != '(erg/cm2/s)'):
            print "Error: Ambiguous line selection. Reduce bandwidth?"
            print "See radex.out for details"
            sys.exit()
        return float(lines[-1].split()[-2])
    def _find_cdens(self):
        import os
        ratio = 0
        iteration = 0
        while (ratio > (1+self.tol)) or (ratio < (1-self.tol)) :
            iteration += 1
            self._create_input()
            os.system('radex < radex.inp > /dev/null')
            mflux  = self._read_output()
            if (mflux < self.eps):
                msg = ("Zero or negative line intensity\n"
                        "See radex.out for details")
                if self.debug:
                    raise Exception(msg)
            if self.debug: print "mflx= ",mflux
            ratio = self.oflux/[mflux,self.oflux][mflux==0]
            self.cdinit *= ratio
            if self.debug: print ratio
            if (iteration > self.maxiter):
                if not self.silent: print "Maximum number of iterations exceeded"
                ratio = 1
        self.mflux = mflux
        return self.cdinit
    def run_model(self):
        from scipy import arange, array
        #
        #
        if self.tovary == 'tkin':
            tkin_inp = self.tkin
            self.cdens = []
            if not self.silent: print 'Running model and varying kinetic temperature.'
            for i in tkin_inp:
                self.tkin = i
                density = self._find_cdens()
                self.cdens.append(density)
                if self.debug: print 'density : ',density
            self.tkin = array(tkin_inp)
        elif self.tovary == 'h2dens':
            h2dens_inp = self.h2dens
            self.cdens = []
            if not self.silent: print 'Running model and varying H2 density.'
            for i in h2dens_inp:
                self.h2dens = i
                density = self._find_cdens()
                if self.debug: print 'density : ',density
                self.cdens.append(density)
            self.h2dens = array(h2dens_inp)
        self.cdens = array(self.cdens)
        if not self.silent: print 'Done!'

class Ratran:
    """
    ----------------------------------------------
    Changelog:
        *20/03/2012 - radial profile grid function

        *14/03/2012 - class created
    """
    def __init__(self, **kwargs):
        """
        r = 0.0, 
        rho_dust = 0.0, 
        temp = 0.0, 
        db = 0.0, 
        abund = 0.0, 
        vr = 0.0, 
        tdust = 0.0, 
        dustonly = 0, 
        modelfile = 'transphere.mdl', 
        dpc = 0.0, 
        imsize = 129, 
        pixel = 0.5, 
        trans = '220.0e9', 
        writeonly = 0, 
        skyonly = 0, 
        molfile='', 
        ncell = 50, 
        outputfile="ratranResult", 
        snr=20, 
        fixset=1e-6, 
        minpop=1e-4, 
        unit='Jypx', 
        opstates=0, 
        gas2dust=100, 
        nphot=1000, 
        temp_limit=10.0, 
        rho_limit=1E-4, 
        pxl_radius = 32, 
        los = 2
        
        
        Only supports one collision partner. If both para- and ortho-H2
        is tabulated in \'molfile\', please check the validity and if
        the gas-to-dust ratio has to be changed.

        Input :

        """

        from scipy import array
        # imports
        import cgsconst as _cgs
        from scipy import zeros, array, logspace, log10, where, pi
        import scipy.interpolate
        import sys
        import os
        import numpy as np
        from time import time

        # Checking input parameters
        params = [
        'r',                            0,          'cm',    'array',   # Radial points
        'rhodust',                      0,       'g/cm3',    'array',   # Dust density
        'molfile',      'ph2-18o-ph2.dat',            '',      'str',   # Name of moldata file
        'modelfile',     'transphere.mdl',            '',      'str',   # Name of model input file
        'outputfile',      'ratranresult',            '',      'str',   # Name of output file
        'db',                         0.0,        'km/s',    'float',   # 1/e half-width of line profile (Doppler b-parameter) (km s-1)
        'abundtype',               'jump',            '',      'str',   # Molecular abundance type ['jump', '?']
        'tjump',                    100.0,           'K',    'float',   # If 'jump' profile, at what T should the jump be
        'xs',                [1E-4, 1E-9],    'relative',     'list',   # If 'jump' profile, what is [inner, outer] relative abundance
        'vr',                         0.0,        'km/s',    'array',   # Radial velocity
        'tdust',                      0.0,           'K',    'array',   # Dust temperature profile
        'dustonly',                 False,            '',      'bool',  # 
        'skyonly',                  False,            '',     'bool',   # 
        'writeonly',                False,            '',     'bool',   # Only write the input files (obsolete)
        'temp',                         0,           'K',    'array',   # 
        'templim',                   10.0,           'K',    'float',   #
        'nh2lim',                    1E-4,      'g/cm-3',    'float',   #
        'trans',                        7,            '',    'mixed',   # Transition number(s) as int, or list
        'dpc',                        0.0,          'pc',    'float',   # Distance to source
        'imsize',                     129,      'pixels',      'int',   # Number of pixels in the output image
        'pixel',                      0.5,    'asec/pxl',    'float',   # Pixel size in arcseconds
        'pxlradius',                   32,            '',      'int',   # Region (in numbers of pixels radius w.r.t. image center) over which to use multiple lines of sight (los)
        'los',                          2,            '',      'int',   # Number of lines of sight
        'unit',                    'Jypx',            '',      'str',   # Output units ['Jypx', 'K', 'Wm2Hzsr']
        'trans',                        7,            '',      'int',   # Transition numbers
        'snr',                       10.0,       'ratio',    'float',   # Requested minimum signal-to-noise
        'fixset',                    1E-6,            '',    'float',   # Convergence requirement for first stage
        'minpop',                    1E-4,            '',    'float',   # Minimum population to include in S/N calculation
        'nphot',                     1000,            '',      'int',   # Number of photons
        'opstates',                 False,            '',     'bool',   # ?
        'ncell',                       20,            '',      'int',   # Number of grid cells
        'gas2dust',                 100.0,       'ratio',    'float',   # Gas to dust ratio to be used in the run
        'directory',       'ratr_model_1',    'dir name',      'str']   # Directory to work in
        # if loadfile is input, drop everythin and just load the file 
        # and, if it exists, the output
        
        print ('Model created with the following parameters:')
        
        class InputParameters: pass
        param_zip = zip(params[0::4],params[1::4], params[2::4], params[3::4])
        for par, stdval, unit, typ in param_zip:
            if par in kwargs: # if input was given, save it
                value = kwargs[par]
                #~ printvalue = kwargs[par]
                #~ if par in ['r', 'rho_dust', 'temp']:
                    #~ printvalue = printvalue[0]
                #~ # print the value that was input for that given parameter
                #~ print '   {0:9} : {1} {2}'.format(par, str(printvalue).strip('\"'), unit)
                
            elif par not in kwargs: # if not input was given, use default
                if stdval != None:
                    # print the default value for the parameter
                    #~ print '   {0:9} : {1:15} {2:8} \t(default)'.format(par, stdval, unit)
                    value = stdval
                else:
                    raise Exception('Wrong parameter input/handling'
                                    'Check input and/or code...')
            # Check what type it should be and store it 
            #string
            if typ == 'str':
                #~ kwargs[par] = '{0}'.format(kwargs[par])
                InputParameters.__dict__[par] = str(value)
            #integer (not sure about 'itypemw' perhaps boolean)
            elif typ == 'int':
                InputParameters.__dict__[par] = int(value)
            #float
            elif typ == 'float':
                InputParameters.__dict__[par] = float(value)
            #bool
            elif typ == 'bool':
                InputParameters.__dict__[par] = bool(value)
            # array
            elif typ == 'array':
                InputParameters.__dict__[par] = array(value)
            elif typ == 'list':
                InputParameters.__dict__[par] = list(value)
            # lastly try to fudge it, see if it works
            elif typ == 'mixed':
                #~ try:
                InputParameters.__dict__[par] = value
                #~ except ValueError:
                    #~ try:
                        #~ InputParameters.__dict__[par] = int(value)
                    #~ except ValueError:
                        #~ InputParameters.__dict__[par] = str(value)
            else:
                InputParameters.__dict__[par] = value
        
        # input parameters contains all the input needed to 
        # create this class again
        self.InputParameters = InputParameters
        
        # copy important parameters to the main class
        self.r = self.InputParameters.r
        self.rhodust = self.InputParameters.rhodust
        self.temp = self.InputParameters.temp
        self.tdust = self.InputParameters.tdust        
        self.directory = self.InputParameters.directory
        self.temp = self.InputParameters.temp
        self.db = self.InputParameters.db
        self.vr = self.InputParameters.vr
        
        
        # If no dust temperature given, assume it is in equilibrium 
        # with the gas
        if self.tdust == 0.0:
            self.tdust = self.temp
        
        # calculate the radial dependence of the molecular
        # abundance depends on what type of abundance type is choosen
        self.abund =  create_molecular_abundance(self.temp, 
                                abund_type = self.InputParameters.abundtype, 
                                Tjump = self.InputParameters.tjump, 
                                Xs = self.InputParameters.xs)
        
        # come up with a nicer solution than this...
        # check that none of the required arrays are not empty, or
        # to short
        if len(self.r)<5:
            raise Exception('neeed array as r, rho_dust, abund and temp')
        if len(self.rhodust)<5:
            raise Exception('neeed array as r, rho_dust, abund and temp')
        if len(self.temp)<5:
            raise Exception('neeed array as r, rho_dust, abund and temp')
        #
        # CHECK if directory exists
        input_dir_path = os.path.join(os.getcwd(), self.directory)
        # if the directory exists
        if not make_dirs(input_dir_path): # will raise error if things go south (i.e., permissions not correct [I think...])
            print('Directory exists, continuing.')
        
        save_ratran(self)
        # rewrite this part when changing to object oriented
        # now it is very hack-ish and non pythonic
        ################################################################
        # MOLECULAR H2 NUMBER DENSITY
        # 
        # TODO : why time gas2dust here? <- changed to times 100
        # what is correct, 100 for the standard assumption
        # rhodust is g/cm3
        # from dust density to number density (cm-3)
        #   cm-3 =  g/cm3 * 100 / (muH2 * g)
        # 100 is gas:dust, but the input gas2dust is only for the 
        # run itself, i.e. only related to the molecules
        self.nh2 = self.rhodust * 100  / (_cgs.MUH2 * _cgs.MP)
        # nh2 is number density of H2
        # rhodust * 100 = rhogas (g/cm3) (all gas H2+He+Metals)
        # MUH2 * MP =  molecular mass (H2+He+Metals) in g
        # so
        # rhogas / molecular mass = number density of H2 in cm-3
        
        ################################################################
        # ENVELOPE CUT OFF LIMIT
        #Find the envelope cut off for T and n
        #
        # if T goes below tepmlim (10K def) somewhere in the model
        try:
           ind_T = where(self.temp < self.InputParameters.templim)[0].min()
        #~ ind = where(r<(1.2E4*cgs.AU))[0].max()
        except (ValueError):
            ind_T = False
        # if n goes below nh2lim (1E-4 def) somewhere in the model
        try:
            ind_n = where((self.nh2) < self.InputParameters.nh2lim)[0].min()
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
        
        # get values at cut off
        self.r_10k = self.r[ind]
        self.rhodust_10k = self.rhodust[ind]
        self.nh2_10k = self.rhodust_10k * 100 / _cgs.MUH2 / _cgs.MP
        self.temp_10k = self.temp[ind]
        #~ self.Y = self.InputParameters.r.max() / self.InputParameters.r.min()
        self.Y = self.r_10k / self.r.min()
        
        ################################################################
        # get values at r = 1000 AU
        ind_r1000 = where(self.r > 1000 * _cgs.AU)[0].min()
        self.rhodust_r1000 = self.rhodust[ind_r1000]
        self.nh2_r1000 =  self.nh2[ind_r1000]
        self.temp_r1000 = self.temp[ind_r1000]
    
        ################################################################
        # cut off where T<10 K OR nh2<1E-4
        # first we have to remove all cells where T<10 K
        # RATRAN does not work well with them
        # after this you use the self.parameter.
        # TODO : perhaps not the best tactics..?
        self.r = self.r[:ind]           
        self.rhodust = self.rhodust[:ind]
        self.temp = self.temp[:ind]
        self.tdust = self.tdust[:ind]
        self.abund = self.abund[:ind]
        self.nh2 = self.nh2[:ind]
        
        ################################################################
        # Refinement
        # for the refinement, easiest way is to redefine rx!
        #
        rx = logspace(log10(self.r[0]), log10(self.r[-1]), num = self.InputParameters.ncell + 1, endpoint = True)
        rx = np.insert(rx, 0, 0)
        self.r1 = rx[0:-1]
        self.r2 = rx[1:]
        #~ r1=np.insert(r[0:-1],0,0)
        #~ r2=np.array(r)
        
        self.rr = zeros(self.InputParameters.ncell + 1, float)
        self.rr[1:] = 10**( (np.log10(self.r1[1:]) + np.log10(self.r2[1:])) / 2.0 )
        self.rr[0] = self.rr[1]
        # Interpolate the values to 'ncell' cells
        self.nh2f = scipy.interpolate.interp1d(log10(self.r), log10(self.nh2))
        self.tkf = scipy.interpolate.interp1d(log10(self.r), log10(self.temp))
        self.tdf = scipy.interpolate.interp1d(log10(self.r), log10(self.tdust))
        self.abund_f = scipy.interpolate.interp1d(log10(self.r), log10(self.abund))
        #
        # Convert logarithms to floats
        self.nh2int = 10**self.nh2f(log10(self.rr))
        self.tkint = 10**self.tkf(log10(self.rr))
        self.tdint = 10**self.tdf(log10(self.rr))
        self.abund_int = 10**self.abund_f(np.log10(self.rr))
        ############################
        #~ nhint_p = nhint*2/4.
        #~ nhint_p[0] = 0.0
        #~ nhint_o = nhint*2/4.
        #~ nhint_o[0] = 0.0
        #~ teint = 10**tkf(log10(rr))
        ############################
        # nh2 needs to start at 0
        self.nh2int[0] = 0.0
    
        # mass of it all
        #~ vol=[]
        #~ mass=[]
        # V = 4*pi*r**3/3
        # r in cm (?)
        V = 4 * pi * (self.r2**3 - self.r1**3) / 3         # cm3
        self.M = V * self.nh2int * _cgs.MUH2 * _cgs.MP # cm3 * g/cm3 ?
        self.M /= _cgs.MSUN
        
        # to get the column density, integrate over radius r1 to r_10k
        #r_10k * 2 nh2_10k
    
        ################################################################
        #  print info. -> move to __str__ method
        print ('M_10K   : {0:<7.2f} Msun\n'
                'R_10K   : {1:<7.0f} AU\n'
                'nH2_10K : {2:<7.1e} cm-3\n'
                'Y       : {3:<7.0f}\n'
                'T       : {4:<7.1f} K\n'.format(self.M.sum(),
                                            self.r_10k/_cgs.AU,
                                            self.nh2_10k,
                                            self.Y,
                                            self.temp_10k))
        print 'Constraining the envelope : ', r_constraint
        print ('nH2_r1000   : {0:<7.1e} cm-3\n'
                'T_r1000     : {1:7.1f} K\n'.format(self.nh2_r1000,
                                             self.temp_r1000))

    def write_input(self):
        """
        id    : shell number
        ra,rb : inner & outer radius (m)
        za,zb : lower & upper height (m) (2D only)
        nh    : density (cm-3) of main collision partner (usually H2)
        nm    : density (cm-3) of molecule
        ne    : density (cm-3) of second collision partner (e.g. electrons)
        tk    : kinetic temperature (K) 
        td    : dust temperature (K)
        te    : electron' temperature (K)
        db    : 1/e half-width of line profile (Doppler b-parameter) (km s-1)
        vr    : radial velocity (km s-1)
        """
        with open(_os.path.join(self.directory, self.InputParameters.modelfile),'w') as f: 
            f.write('# Ratran input file based on Transphere results'+'\n')
            if self.InputParameters.skyonly: 
                f.write('# ... intended for (SKY) continuum calculations only.'+'\n')
            f.write('rmax={0:.5E}\n'.format( self.r2[-1] / 100 ))       # rmax in METERS (convert from cm i.e. / 100)
            f.write('ncell=%i'  % (len(self.r2))+'\n')
            f.write('tcmb=2.728\n')
            f.write('columns=id,ra,rb,nh,nm,tk,td,db,vr\n')
            f.write('gas:dust={0}\n'.format(self.InputParameters.gas2dust))
            if self.InputParameters.skyonly: 
                f.write('kappa=jena,thin,e6\n')
            f.write('@\n')
            # r1/r2 in meter (convert from cm)
            """
            id    : shell number
            ra,rb : inner & outer radius (m)
            za,zb : lower & upper height (m) (2D only)
            nh    : density (cm-3) of main collision partner (usually H2)
            nm    : density (cm-3) of molecule
            ne    : density (cm-3) of second collision partner (e.g. electrons)
            tk    : kinetic temperature (K) 
            td    : dust temperature (K)
            te    : electron' temperature (K)
            db    : 1/e half-width of line profile (Doppler b-parameter) (km s-1)
            vr    : radial velocity (km s-1)
            """
            for ii in range(0, len(self.r1)):
                f.write("%4i %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E" % (ii+1, self.r1[ii]/100.0, self.r2[ii]/100.0, self.nh2int[ii], self.nh2int[ii]*self.abund_int[ii], self.tkint[ii], self.tdint[ii], self.db, self.vr)+'\n')
            
        if not self.InputParameters.skyonly:
            if self.InputParameters.molfile == '':
                sys.exit('Error: for AMC calculations the molecular datafile (molfile) needs to be set.')
            with open(_os.path.join(self.directory, "amc.inp"),'w') as f:
                f.write("source={0}\n".format(self.InputParameters.modelfile))
                f.write("outfile=populations.pop\n")
                f.write("molfile={0}\n".format(self.InputParameters.molfile))
                f.write("snr={0}\n".format(self.InputParameters.snr))
                f.write("nphot={0}\n".format(self.InputParameters.nphot))
                f.write("kappa=jena,thin,e6\n")
                f.write("minpop={0:3.2E}\n".format(self.InputParameters.minpop))
                f.write("seed=1971\n")
                f.write("fixset={0:3.2E}\n".format(self.InputParameters.fixset))
                f.write("go\n")
                f.write("q\n")
                f.write("\n")


        with open(_os.path.join(self.directory, "sky.inp"),'w') as f:
            if not self.InputParameters.skyonly:
                f.write("source=populations.pop\n")
            else:
                f.write("source="+self.InputParameters.modelfile+"\n")
            f.write("format=miriad\n")
            f.write("outfile="+self.InputParameters.outputfile+"\n")
            if type(self.InputParameters.trans) == type([]):
                li = [str(i) for i in self.InputParameters.trans]
                f.write("trans={0}\n".format(','.join(li)))
            else:
                f.write("trans={0}\n".format(self.InputParameters.trans))
            f.write("pix="+str(self.InputParameters.imsize)+","+str(self.InputParameters.pixel)+",32,2\n")
            if not self.InputParameters.skyonly:
                f.write("chan=50,0.2\n")
            else:
                f.write("chan=1,1.0\n")
            f.write("distance="+str(self.InputParameters.dpc)+"\n")
            f.write("units={0}\n".format(self.InputParameters.unit))
            f.write("go\n")
            f.write("q\n")
            f.write("\n")
    def __str__(self):
        return 'Info. about the model'

    def save(self, filename='ratran_input.pick'):
        import pickle
        with open(os.path.join(), 'w') as f:
            pickle.dump(vars(self.InputParameters), f)

    def run(self):
        # run ratran with the setup in the directory
        #
        # catch : ### WARNING
        # AMC: minimum S/N  |  converged  |     photons  |  increase to
        # AMC: -------------|-------------|--------------|-------------
        # and the output after this..
        # e.g.
        # AMC: 1.38711E+00  |      78.43% |       89000  |      116000
        # import modules, subprocess replaces os
        import subprocess
        from time import time, sleep
        import sys
        # 
        if not self.InputParameters.writeonly:
            if not self.InputParameters.skyonly:
                #~ print "Starting AMC calculation..."
                with ChangeDirectory(self.directory):
                    t1 = time()
                    proc = subprocess.call(['amc', 'amc.inp'])
                print('\nAMC took {0:2.1f} seconds'.format((time()-t1)))
            with ChangeDirectory(self.directory):
                t1 = time()
                proc = subprocess.call(['sky', 'sky.inp'])
                print('\nSKY took {0:2.1f} seconds'.format((time()-t1)))
                    
          
        #~ # now get the output as a string, so we can check stuff if
        #~ # we want to
        #~ self.transphere_output = ''.join(trans_out)
        #~ # read in the output-files, for checks, plots etc
        #~ self.Envstruct, self.Convhist = read_transphereoutput(self)

class Transphere:
    """
    Title
    ------------
    **test**
    Class to interface with Transphere - dust continuum
    radiative transfer code.
    
    
    
    ######################## copied info from old "pyCollapse.py"
    # floats
    rstar    = 5.3453e0 * nc.RS  # Stellar radius
    mstar    = 1 * nc.MS         # Stellar mass
    tstar    = 5000.             # Stellar temperature
    rin      = 200. * nc.AU      # Inner radius of shell
    rout     = 1.2e4 * nc.AU     # Outer radius of shell
    r0       = 1e3 * nc.AU       # Reference radius
    rho0     = 3.46e-18          # Density at reference radius ()
    plrho    = -1.5              # Powerlaw for rho
    
    ### Parameters related to control of code
    
    rref     = 500. * nc.AU      # Refinement radius
    
    # integers
    nr       = 200               # Nr of radial shells for grid
    nref     = 100               # Nr of refinement points
    nriter   = 30                # Maximum nr of iterations
    convcrit = 0.00001           # Convergence criterion
    ncst     = 10                # Nr of rays for star
    ncex     = 30                # Nr of rays between star and Rin
    ncnr     = 1                 # Nr of rays per radial grid point
    itypemw  = 1                 # Type of mu weighting
    
    # booleans
    idump    = 1                 # Dump convergence history
    localdust= 0                 # Dust opacity local?
    
    
    output: 
    
    
    """
    def __init__(self, **kwargs):
        # imports
        from scipy import log10, log, arange, linspace, logspace,\
        array
        import os
        
        #
        #~ self.silent = silent
        #
        
        
        #~ class Opacity: pass
        """
        localdust : Dust opacity local?
        silent    : Verbose?
        nriter    : Maximum nr of iterations
        convcrit  : Convergence criterion
        ncst      : Nr of rays for star
        ncex      : Nr of rays between star and Rin
        ncnr      : Nr of rays per radial grid point
        itypemw   : Type of mu weighting
        idump     : Dump convergence history
        """
        # Checking input parameters
        params = array([
        'rin',          20,             'AU',   'float',    # Inner radius of shell
        'rout',         8000,           'AU',   'float',    # Outer radius of shell
        'nshell',       200,            ' ',    'int',      # Number of shells
        'spacing',      'powerlaw1',    ' ',    'str',      # Type of grid spacing
        'nref',         0,              ' ',    'int',      # Refinement shells
        'rref',         0.0,            'AU',   'float',    # Refinement radius
        'rstar',        3.0,            'AU',   'float',    # Stellar radius
        'tstar',        5780,           'K',    'float',    # Stellar temperature
        'mstar',        1,              'MSUN', 'float',    # Stellar mass
        'isrf',         0.0,            ' ',    'mixed',    # Scaling of ISRF
        'tbg',          2.73,           'K',    'float',    # Spectral shape of ISRF (Blackbody equivalent temperature). With -1 the spectrum is read from the file isrf.inp and scaled by ISRF.
        'dpc',          250,            'PC',   'float',    # Distance to source in pc
        'r0',           1.0E3,          'AU',   'float',    # Reference radius
        'plrho',        -1.5,           ' ',    'float',    # Powerlaw for rho
        'rho_type',     'powerlaw1',    ' ',    'str',      # Type of rho dependence with radius ['powerlaw1', 'shu-knee']
        'n0',           2e6,            'cm-3', 'float',    # H2 number density at reference radius
        'gas2dust',     100.0,          ' ',    'float',    # Gas to dust ratio
        'localdust',    False,          ' ',    'bool',     # Dust opacity local?
        'silent',       True,           ' ',    'bool',     # Verbose?
        'nriter',       30,             ' ',    'int',      # Maximum nr of iterations
        'convcrit',     1E-5,           ' ',    'float',    # Convergence criterion
        'ncst',         10,             ' ',    'int',      # Nr of rays for star
        'ncex',         30,             ' ',    'int',      # Nr of rays between star and Rin
        'ncnr',         1,              ' ',    'int',      # Nr of rays per radial grid point
        'itypemw',      1,              ' ',    'int',      # Type of mu weighting
        'idump',        True,           ' ',    'bool',     # Dump convergence history
        'opacfile',     0,              ' ',    'str',      # Filename of opacityfile
        'freqfile',     0,              ' ',    'str',      # Filename of frequency file
        'directory',    '' ,            ' ',    'str'       # Directory to work in
                    ])
        print ('Model created with the following parameters:')
              
########################################################################
        # need to check that the input is correct.
        # strings should be strings, and floats floats
        class InputParameters: pass
        param_zip = zip(params[0::4],params[1::4], params[2::4], params[3::4])
        for par, stdval, unit, typ in param_zip:
            if par in kwargs: # if input was given, save it
                # print the value that was input for that given parameter
                print '   {0:9} : {1} {2}'.format(par, str(kwargs[par]).strip('\"'), unit)
                value = kwargs[par]
            elif par not in kwargs: # if not input was given, use default
                if stdval != None:
                    # print the default value for the parameter
                    print '   {0:9} : {1:7} {2:4} \t(default)'.format(par, stdval, unit)
                    value = stdval
                else:
                    raise Exception('Wrong parameter input/handling'
                                    'Check input and/or code...')
            # Check what type it should be and store it 
            #
            #string
            if typ == 'str':
                #~ kwargs[par] = '{0}'.format(kwargs[par])
                InputParameters.__dict__[par] = str(value)
            #integer (not sure about 'itypemw' perhaps boolean)
            elif typ == 'int':
                InputParameters.__dict__[par] = int(value)
            #float
            elif typ == 'float':
                InputParameters.__dict__[par] = float(value)
            #bool
            elif typ == 'bool':
                InputParameters.__dict__[par] = bool(value)
            # just try to fudge it, see if it works
            elif typ == 'mixed':
                # TODO : change mixed to the same as in the 
                # ratran object
                try:
                    InputParameters.__dict__[par] = int(value)
                except ValueError:
                    InputParameters.__dict__[par] = str(value)
            else:
                InputParameters.__dict__[par] = value
        
        # input parameters contains all the input needed to 
        # create this class again
        self.InputParameters = InputParameters
        
        # copy important parameters to the main class
        self.rin = self.InputParameters.rin
        self.rout = self.InputParameters.rout
        self.rref = self.InputParameters.rref
        self.nref = self.InputParameters.nref
        self.nshell = self.InputParameters.nshell
        self.gas2dust = self.InputParameters.gas2dust
        self.n0 = self.InputParameters.n0        
        self.r0 = self.InputParameters.r0   
        self.rstar = self.InputParameters.rstar
        self.mstar = self.InputParameters.mstar
        self.isrf = self.InputParameters.isrf
        self.tbg = self.InputParameters.tbg
        self.tstar = self.InputParameters.tstar
        self.dpc = self.InputParameters.dpc
        
        self.directory = self.InputParameters.directory
        
        # CHECK if directory exists 
        # create dir if it doesn't exist
        # will raise error if things go south 
        # (i.e., permissions not correct [I think...])
        input_dir_path = os.path.join(os.getcwd(), self.directory)
        if not make_dirs(input_dir_path): 
            print('Directory exists, continuing.')
        
        _subprocess.Popen(['cp', self.InputParameters.freqfile, self.directory])
        
        save_transphere(self)
        
        # Convert distances to CGS units
        self.rin   *= _cgs.AU   # cm
        self.rout  *= _cgs.AU   # cm
        self.rref  *= _cgs.AU   # cm
        self.r0    *= _cgs.AU   # cm
        self.rstar *= _cgs.RSUN # cm
        self.mstar *= _cgs.MSUN # g
        
        if self.rref:
            # if we want refinement, just call the grid function twice
            # using the self.rin, self.rref and self.nref for first grid
            # and self.rref and self.rout and self.nshell for second
            from scipy import concatenate
            print('Creating grid with refinement.')
            inner_grid = create_grid(
                        self.rin, self.rref, self.InputParameters.nref,
                        space=self.InputParameters.spacing, end=True
                                    )
            outer_grid = create_grid(
                        self.InputParameters.rref, self.rout, self.nref,
                        space=self.InputParameters.spacing, end=True
                                    )
            radii = concatenate([inner_grid[:-1], outer_grid]).ravel()
            #~ return inner_grid, outer_grid, radii
        else:
            # if we dont want refinement
			radii = create_grid(self.rin, 
                                self.rout, 
                                self.InputParameters.nshell,
                                space = self.InputParameters.spacing, 
                                end = True)
        #
        # radii is now in cm
        # separate the different radii, the lower/upper bounds
        # and the mid-point of the grid "cells"
        lower = radii[:-1]
        upper = radii[1:]
        self.radii = radii
        self.r_midpt = array((lower+upper)/2)
        self.r_grid = array([array([i,j]) for i,j in zip(lower,upper)])

        #
        # Create the density array
        #
        # calculate the density at all radial points
        if self.InputParameters.rho_type == 'powerlaw1':
            # 1E-2  - gas-to-dust
            # get the DUST density
            #~ self.rho_gas = self.rho0_gas * (self.radii / self.r0)**(self.plrho)
            r_dependence = (self.radii / self.r0)**(self.InputParameters.plrho)
        elif self.InputParameters.rho_type == 'shu_knee':
            #~ r_dependence = (self.radii / self.r0)**(self.plrho)
            #
            # what should the parameter be?
            # where should the knee be?
            #~ if rho_flat == 0:
                #~ rho    = (1/gas2dust)*rho0 * (r/r0)**(plrho) # Gas2dust = 100
            #~ else:
                #~ # two component rho, with a "knee" in the middle
                #~ rho    = (1/gas2dust)*rho0 * (r/r0)**(plrho) # Gas2dust = 100
                #~ for i in range(len(rho)):
                    #~ if rho[i] < (1/gas2dust)*rho_flat:
                        #~ rho[i] = (1/gas2dust)*rho_flat
            pass
        else:
            raise Exception('rho_type parameter not recognised')
        # now also save the dust density, and H2 number density
        #~ self.rho_gas = self.rho0_gas * (self.radii / self.r0)**(self.plrho)
        #~ self.n_h2 = self.rho_gas / (_cgs.MUH2 * _cgs.MP)
        #~ self.rho_dust = self.n_h2 * 1 / self.gas2dust
        #
        # n0 is in nh2 / cm3
        #    g/cm3           nh2 / cm3 * amu/h2 * mass of amu(g)
        self.rho0_gas = self.n0 * _cgs.MUH2 * _cgs.MP   # g * cm-3
        # all gas, MUH2 is including H2+He+Metals
        # apply the radial dependence to the gas density
        self.rho_gas = self.rho0_gas * r_dependence     # g * cm-3
        # apply radial dependence to the H2 number density
        self.n_h2 = self.n0 * r_dependence              # nh2 * cm-3
        # gas to dust is mass relationship
        # g * cm-3 
        self.rho_dust = self.rho0_gas / self.InputParameters.gas2dust * r_dependence 
        
        ################################################################        
        from scipy import array
        from sys import exit as sysexit
        #~ import os
        # if a dust opacity file is given
        # that is it is not tabulated correctly
        # as a dustopac.inp file
        if self.InputParameters.opacfile not in [0, '0']:
            if self.InputParameters.freqfile in [0, '0']:
                sysexit('No frequency file given (freqfile). Cannot proceed. '
                'Need two files: one with frequency, and one with opacity'
                ' at corresponding frequency.')
            with open(self.InputParameters.opacfile, 'r') as f:
                lines = f.read().split('\n')
            # get the "header" info,
            # that is the number of tabulated entries in
            # absorption and scattering each.
            try:
                nf, ns = [int(i) for i in lines[:2][0].split()]
            except:
                errmsg = 'Error parsing Opacity-file. Wrong format.'
                raise Exception(errmsg)
            if nf >= 1:
                # shape of the arrays, for absorption and scattering
                # (nf, ns)
                # so we have to nest the list comprehensions, and
                # call float() on the smallest possible element
                # (i.e. when ns > 1).
                try:
                    cabs = array([[float(i) for i in j.split()] for j in lines[2 : nf+2]])
                    csca = array([[float(i) for i in j.split()] for j in lines[nf+2 : 2*nf+2]])
                except:
                    errmsg = 'Error parsing Opacity-file. Wrong format.'
                    raise Exception(errmsg)
                nrtrange = 1
                trange = [0.e0,0.e0]
            else:
                ########################################################
                # this part not edited, dont have a file like this
                sysexit('This part of the code is not tested'
                        'Please optimize code and check it.')
                nf = ns
                ismooth = 0
                nrtrange = 0
                with open(self.InputParameters.opacfile, 'r') as f:
                    line = f.readline().strip()
                    ns, ismooth, nrtrange = line.split()
                    ns, ismooth, nrtrange = int(ns), int(ismooth), int(nrtrange)
                    if ismooth != 0:
                        import sys
                        sys.exit('Error: Smoothing not yet allowed.')
                    import numpy as np
                    cabs = np.zeros((nf, ns), float)
                    csca = np.zeros((nf, ns, nrtrange), float)
                    dum = 0.e0
                    trange = np.zeros(nrtrange+1, float)
                    
                    for ir in range(0, nrtrange):
                        a, b = f.readline().strip().split()
                        a, b = int(a), int(b)
                        trange[ir] = b
                        for kk in range(0, nf):
                            for iss in range(0, ns):
                                dum = float(f.readline().split())
                                cabs[kk, iss, ir] = dum
                        for kk in range(0,nf):
                            for iss in range(0, ns):
                                dum = float(f.readline().split())
                                csca[kk, iss, ir] = dum
                ########################################################
            with open(self.InputParameters.freqfile, 'r') as f:
                lines = f.read().split('\n')
            nf = int(lines[0])
            freq = array(lines[2:2+nf], 'float')
            wave = 2.9979e14 / freq

            ## TODO : change this, do I even need a dictionary?
            # move the calls to where the variable is defined
            class Opacity: pass
            self.Opacity = Opacity
            
            self.Opacity.ns         = ns
            self.Opacity.nf         = nf
            self.Opacity.freq       = freq
            self.Opacity.wave       = wave
            self.Opacity.cabs       = cabs
            self.Opacity.csca       = csca
            self.Opacity.nrtrange   = nrtrange
            self.Opacity.trange     = trange
            
            self.freq               = freq
            #~ self.opacity = {'ns': ns, 'nf': nf, 'freq': freq, 'wave': wave, 'cabs': cabs, 'csca': csca, 'nrt': nrtrange, 'trange': trange}

        else:
            # always neeed a opacity file and a frequency file
            sysexit('No opacity-file given.')
        ################################################################
        # if local dust opacity
        # the opacity is tabulated with radius
        #
        # have not got a file like this, so i am not changing it

        if self.InputParameters.localdust:
            print ('variable localdust not False/0, this part of the code'
            'is not up to date.')
            sysexit('This part of the code is not tested'
                    'Please optimize code and check it.')
            if 'nr' not in kwargs:
                sysexit('Number of radial shells, \"nr\" not given as input.')
            _os.system('rm -f dustopac_1.inp')
            _os.system('rm -f dustopac.inp') # added this, do I really need
                                            # to remove it too?
            f = open(_os.path.join(self.directory, 'dustopac.inp'), 'w')
            f.write('1               Format number of this file'+'\n')
            f.write('1               Nr of dust species'+'\n')
            f.write('============================================================================'+'\n')
            f.write('-1              Way in which this dust species is read (-1=file)'+'\n')
            f.write('0               Extension of name of dustopac_***.inp file'+'\n')
            f.write('----------------------------------------------------------------------------'+'\n')
            f.close
            f = open(os.path.join(self.directory, 'dustopac_0.inp'), 'w')
            f.write(nr)
            f.write(str(nf)+' 1\n')
            f.write(' ')
            redux=1.e0
            
            for ir in range(0,nr):
                for inu in range(0,nr):
                    f.write(self.Opacity.cabs[inu] * redux)
                for inu in range(0,opacity['nf']):
                    f.write(self.Opacity.csca[inu] * redux)
                f.write(' ')
            f.close
        elif not self.InputParameters.localdust:
            # first remove the standard ratran dust opacity input files
            _os.system('rm -f {0}'.format(_os.path.join(self.directory, 'dustopac_0.inp')))
            _os.system('rm -f {0}'.format(_os.path.join(self.directory, 'dustopac.inp'))) # added this, do I really need
                                            # to remove it too?
            with open(_os.path.join(self.directory, 'dustopac.inp'),'w') as f:
                f.write('1               Format number of this file\n')
                f.write('1               Nr of dust species\n')
                f.write('============================================================================\n')
                f.write('-1              Way in which this dust species is read (-1=file)\n')
                f.write('1               Extension of name of dustopac_***.inp file\n')
                f.write('----------------------------------------------------------------------------\n')
            with open(_os.path.join(self.directory, 'dustopac_1.inp'), 'w') as f:
                f.write(str(self.Opacity.nf)+' 1\n \n')
                for inu in range(0, self.Opacity.nf):
                    f.write(str(self.Opacity.cabs[inu][0])+'\n')
                for inu in range(0, nf):
                    f.write(str(self.Opacity.csca[inu][0])+'\n')
    
    def write_transphereinput(self):
        from scipy import pi, zeros
        import os
        # Transphere input file
        text = ('{0}\n{1}\n{2}\n{3}\n'
                '{4}\n{5}\n{6}\n{7}'.format(2,
                self.InputParameters.nriter,
                self.InputParameters.convcrit,
                self.InputParameters.ncst,
                self.InputParameters.ncex,
                self.InputParameters.ncnr,
                self.InputParameters.itypemw,
                int(self.InputParameters.idump)))
        with open(_os.path.join(self.directory, 'transphere.inp'),'w') as f:
            f.write(text)
        #
        # Make the stellar information file
        # (mstar and tstar are irrelevant; they are there for historical reasons)
        # isn't "tstar" used for planck calc?
        #~ f=open('starinfo.inp','w')
        with open(_os.path.join(self.directory, 'starinfo.inp'),'w') as f:
            f.write('1\n'
                    '{0}\n'
                    '{1}\n'
                    '{2}\n'.format(self.rstar,
                                self.mstar,
                                self.tstar))
        #
        # The stellar spectrum
        #
        #~ f=open('starspectrum.inp','w')
        sspec = (self.rstar / _cgs.PC)**2 * pi * bplanck(self.freq, self.tstar)
        with open(_os.path.join(self.directory, 'starspectrum.inp'), 'w') as f:
            f.write('{0}\n'.format(len(self.freq)))
            for inu in range(0,len(self.freq)):
                f.write('{0:20}\t{1:20}\n'.format(self.freq[inu], sspec[inu]))
        #
        # The exterior spectrum
        #
        if self.tbg == 0.0 or self.isrf == 0:
            bgspec = zeros((len(self.freq)),float)
        elif self.tbg == -1:
            f = open('isrf.inp', 'r')
            nf = int(f.readline().strip())
            bgspec = zeros((len(self.freq)),float)
            for ii in range(0,nf):
                bgspec[ii] = float(f.readline().strip())*self.isrf
            f.close()
        else:
            if self.tbg > 0: bgspec = bplanck(self.freq, self.tbg) * self.isrf

        with open(_os.path.join(self.directory, 'external_meanint.inp'), 'w') as f:
            f.write('{0}\n'.format(len(self.freq)))
            for inu in range(0, len(self.freq)):
                f.write('{0:20}\t{1:20}\n'.format(self.freq[inu], bgspec[inu]))
        #
        # Write the envelope structure
        #
        with open(_os.path.join(self.directory, 'envstruct.inp'),'w') as f:
            f.write(str(len(self.radii))+'\n')
            f.write(' '+'\n')
            # rho_dust has unit g/cm3
            for ir in range(0,len(self.radii)):
                f.write("%13.6E %13.6E %13.6E" % (self.radii[ir], self.rho_dust[ir], 0.e0)+'\n') # ,format='(3(E13.6,1X))'
    
    def run(self):
        import subprocess
        from time import time, sleep
        import sys
        
        #
        # re-check the input files here?
        # at least the format?
        #
        print ('Running Transphere')
        # temporarily change directory, if it has been initiated
        # the directory should exist
        with ChangeDirectory(self.directory): 
            t1 = time()
            proc = subprocess.Popen(['transphere'], 
                                    stdout = subprocess.PIPE, 
                                    stderr = subprocess.STDOUT)
            sys.stdout.write('Iteration no : ')
            sys.stdout.flush() # flush output so we can write again
            trans_out = []
            while True:
                # first : if process is done, break the loop
                if proc.poll() != None: 
                    break
                nextline = proc.stdout.readline()
                trans_out.append(nextline)
                
                if "Iteration" in nextline:
                    # grab the iteration number
                    iter_no = int(nextline[10:])
                    sys.stdout.write('{0:3} \b\b\b\b'.format(iter_no)) # print out a point
                    sys.stdout.flush()
                #~ if "Error" in nextline:
                    # if it is the error of the first iteration
                    # we grab the error of it
                    #~ if iter_no == 1:
                        #~ start_error = float(nextline[13:])
                        #~ first = 0 # now we are not at the first any longer
                    # if it is not the first error, calculate per cent done
                    #~ else:
                        #~ current_error = float(nextline[13:])
                        #~ diff1 = start_error - self.convcrit
                        #~ diff2 = start_error - current_error
                        #~ p_done = (diff2 / diff1 * 100)
                        #~ sys.stdout.write(' {0} : {1}\r'.format(iter_no, p_done)) # print out a point
                        #~ sys.stdout.flush()    # flush output so we can write again
                #~ sleep(0.5)            # wait for 0.5 second
                
                #~ sys.stdout.write('Downloading File FooFile.txt [%d%%]\r'%i)
                #~ ys.stdout.flush()

            print('\nDone in {0:2.1f} seconds'.format((time()-t1)))
        # now get the output as a string, so we can check stuff if
        # we want to
        self.transphere_output = ''.join(trans_out)
        # read in the output-files, for checks, plots etc
        self.Envstruct, self.Convhist = read_transphereoutput(self)
        
########################################################################
# depricated
def radial_profile(self, spaced='log10',
                        nshell=200, rin=200*_cgs.AU, rout=8000*_cgs.AU,
                        **kwargs):
    """
    Method to define the radial shells

    spaced : ['log10', 'powerlaw', 'linear']

    """
    #import scipy
    from scipy import log10, log, arange, linspace, logspace,\
    array
    #from cgsconst import *

    if spaced.lower() == 'log10':
        # get the exponent of the start- and
        # stop-radius in input units
        start = [log10(rin), 0][rin == 0]
        stop = log10(rout)
        print start, stop
        radii = logspace(start, stop,
                            num=nshell, endpoint=True)
    elif spaced.lower() == 'linear':
        # linearly spaced grid
        radii = linspace(rin, rout, num=nshell, endpoint=True)
    elif spaced.lower() == 'powerlaw':
        # first check if coefficients to the power-law was given
        #~ if 'exp' in kwargs:
            #~ p_exp = kwargs['exp']
        #~ else: # if not, set it to 2, i.e. r^2
            #~ p_exp = 2
        radii = rin + (rout-rin)*(linspace(rin, rout,
                        num=nshell, endpoint=True)/(rout))**2
        #print('Not implemented yet.')
        #raise ParError(spaced)
    else:
        raise ParError(spaced)

    lower = radii[:-1]
    upper = radii[1:]
    self.r_midpt = array((lower+upper)/2)
    self.r_grid = array([array([i,j]) for i,j in zip(lower,upper)])
    self.radii = radii

def create_density(self, model=[0,2], n0=1E9):
        """

        Create density model

        ----------------------------------------------------------------
        Input:

            model : list/tuple [type, exponent1 (]), exponent2, breakpoint]
                where 'type':
                0 - power law
                1 - linear
                2 - broken power law, break point
                    changed with rhobreak
                and 'exponent2' and 'breakpoint' is only needed when a two
                different functions of radial dependence are needed.

        """

        from scipy import where
        if model[0] == 0: # Power law
            if len(model) > 2: # Broken power law
                # check the model parameters exponent2 and breakpoint
                # here
                pass
                #break_index = where()
                #self.r_midpt
        elif model[0] == 1: # Linear
            if len(model) > 2: # Broken linear
                # check the model parameters exponent2 and breakpoint
                # here
                pass
        else:
            raise ParError(model)
        # calculate the density at the mid-points in all the shells
        self.rho_h2  = n0*(self.r_midpt[0]/self.r_midpt)**float(model[1])
########################################################################
#
