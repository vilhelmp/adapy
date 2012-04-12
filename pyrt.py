#! /usr/bin/env python
# -*- coding: utf-8 -*-


#
#       pyrt.py
#
#
#       Copyright 2012 Magnus Persson <magnusp@nbi.dk>
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
#       ver 0.1 alpha
#
#
#
doc="""
PyRT - Python Radiative Transfer

Module to interface with various radiative transfer tools.

Transphere : Dust continuum radiative transfer

Radex : Escape probability

Ratran : 1D radiative transfer code


"""



########################################################################
# RADIATIVE TRANSFER / MODELING

class Radex:
    """
    Input (* = mandatory input)
    *f      : Frequency - frequency of line (in GHz)
    df      : Uncertainty of line position (defaults to 20 MHz)
    *tkin   : Kinetic temperature  - one or array (Kelvin)
    *h2dens: H2 density - in cm-3 - one or array (Kelvin)
    tbg     : Background temperature - defaults to 2.73 K
    *oflux  : Observed flux in K km/s, to compare
    *lwidth : Line width (FWHM)
    *molf   : Molecular data file (LAMBDA format)
    ftol    : How close can the flux be for it to be the same? (0.1)
    maxiter : Maximum number of iterations (100)
    eps     : Below what is zero or negative line intensity (1E-20)
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
        import os
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
                    raise ParError(par)
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
    def __init__(self, model=None):
        print 'Not implemented yet'

    def _radial_profile(self, spaced='log10',
                            nshell=16, rin=0, rout=8000,
                            **kwargs):
        """
        Method to define the radial shells

        spaced : ['log10', 'powerlaw', 'linear']

        """
        import scipy
        from scipy import log10, log, arange, linspace, logspace,\
        array

        if spaced.lower() == 'log10':
            # get the exponent of the start- and
            # stop-radius in meters
            start = [log10(rin*AU*1e-2), 0][rin == 0] # AU : cm to m
            radii = logspace(start, log10(rout*AU*1e-2),
                                num=nshell, endpoint=True) # AU : cm to m
        elif spaced.lower() == 'linear':
            # linearly spaced grid
            radii = linspace(rin*AU*1e-2, rout*AU*1e-2, num=nshell,
                                endpoint=True)
        elif spaced.lower() == 'powerlaw':
            # first check if coefficients to the power-law was given
            #~ if 'exp' in kwargs:
                #~ p_exp = kwargs['exp']
            #~ else: # if not, set it to 2, i.e. r^2
                #~ p_exp = 2
            radii = rin + (rout-rin)*AU*1e-2*(linspace(rin*AU*1e-2, rout*AU*1e-2, num=nshell,
                                endpoint=True)/(rout*AU*1e-2))**2
            #print('Not implemented yet.')
            #raise ParError(spaced)
        else:
            raise ParError(spaced)
        lower = radii[:-1]
        upper = radii[1:]
        self.r_midpt = array((lower+upper)/2)
        self.r_grid = array([array([i,j]) for i,j in zip(lower,upper)])


    def _create_density(self, model=[0,2], n0=1E9):
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


class Transphere:
    def __init__(self,**kwargs):
        from scipy import array
        # if a dust opacity file is given
        # that is it is not tabulated correctly
        # as a dustopac.inp file
        if 'opacfile' in kwargs:
            if 'freqfile' not in kwargs:
                sysexit('No frequency file given. Cannot proceed. '
                'Need file with opacity file, at corresponding frequency')
            with open(kwargs['opacfile']) as f:
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
                    cabs = array([[float(i) for i in j.split()] for j in lines[2:nf+2]])
                    csca = array([[float(i) for i in j.split()] for j in lines[nf+2:-1]])
                except:
                    errmsg = 'Error parsing Opacity-file. Wrong format.'
                    raise Exception(errmsg)
                nrtrange=1
                trange=[0.e0,0.e0]
            else:
                ########################################################
                # this part not edited, dont have a file like this
                nf = ns
                ismooth = 0
                nrtrange = 0
                with open(kwargs['opacfile']) as f:
                    line=f.readline().strip()
                    line=f.readline().strip()
                    ns,ismooth,nrtrange=line.split()
                    ns,ismooth,nrtrange=int(ns),int(ismooth),int(nrtrange)
                    if ismooth != 0: import sys; sys.exit('Error: Smoothing not yet allowed.')
                    import numpy as np
                    cabs = np.zeros((nf,ns),float)
                    csca = np.zeros((nf,ns,nrtrange),float)
                    dum  = 0.e0
                    trange=np.zeros(nrtrange+1,float)

                    for ir in range(0,nrtrange):
                        a,b=f.readline().strip().split()
                        a,b=int(a),int(b)
                        trange[ir]=b
                        for kk in range(0,nf):
                            for iss in range(0,ns):
                                dum=float(f.readline().split())
                                cabs[kk,iss,ir] = dum

                        for kk in range(0,nf):
                            for iss in range(0,ns):
                                dum=float(f.readline().split())
                                csca[kk,iss,ir] = dum
                ########################################################
                with open(kwargs['freqfile']) as f:
                    lines = f.read().split('\n')

                nf=int(f.readline().strip())
                #   if nnf != nf: sys.exit("ERROR: frequency file has different nr of points as dustopac file")
                dum=f.readline()
                freq = np.zeros(nf,float)
                wave = np.zeros(nf,float)
                for kk in range(0,nf):
                    dum=float(f.readline().strip())
                    freq[kk] = dum
                    wave[kk] = 2.9979e14 / dum
                f.close

                opacity={'ns': ns, 'nf': nf, 'freq': freq, 'wave': wave, 'cabs': cabs, 'csca': csca, 'nrt': nrtrange, 'trange': trange}

                return opacity




        else:
            print ('No opacity-file given, assuming '
                    'correct format exists.')


    def readopac(nr='1'):
        import numpy as np
        import sys

        if nr == -1: nr=1
        filename = 'dustopac_'+str(nr)+'.inp'
        print "Reading "+filename
        f=open(filename,'r')
        nf,ns = f.readline().strip().split()
        nf,ns=int(nf),int(ns)
        f.readline()
        if nf >= 1:
            cabs = np.zeros((nf,ns),float)
            csca = np.zeros((nf,ns),float)
            dum  = 0.e0
            print ns
            for kk in range(0,nf):
                for iss in range(0,ns):
                    dum=float(f.readline().strip())
                    cabs[kk,iss] = dum

            for kk in range(0,nf):
                for iss in range(0,ns):
                    dum=float(f.readline().strip())
                    csca[kk,iss] = dum

            nrtrange=1
            trange=[0.e0,0.e0]

        else:
            nf=ns
            ismooth=0
            nrtrange=0
            line=f.readline().strip()
            line=f.readline().strip()
            ns,ismooth,nrtrange=line.split()
            ns,ismooth,nrtrange=int(ns),int(ismooth),int(nrtrange)
            if ismooth != 0: sys.exit('Error: Smoothing not yet allowed.')
            cabs = np.zeros((nf,ns),float)
            csca = np.zeros((nf,ns,nrtrange),float)
            dum  = 0.e0
            trange=np.zeros(nrtrange+1,float)

            for ir in range(0,nrtrange):
                a,b=f.readline().strip().split()
                a,b=int(a),int(b)
                trange[ir]=b
                for kk in range(0,nf):
                    for iss in range(0,ns):
                        dum=float(f.readline().split())
                        cabs[kk,iss,ir] = dum

                for kk in range(0,nf):
                    for iss in range(0,ns):
                        dum=float(f.readline().split())
                        csca[kk,iss,ir] = dum

        f.close
        file='frequency.inp'
        f=open(file,'r')
        nf=int(f.readline().strip())
     #   if nnf != nf: sys.exit("ERROR: frequency file has different nr of points as dustopac file")
        dum=f.readline()
        freq = np.zeros(nf,float)
        wave = np.zeros(nf,float)
        for kk in range(0,nf):
            dum=float(f.readline().strip())
            freq[kk] = dum
            wave[kk] = 2.9979e14 / dum
        f.close

        opacity={'ns': ns, 'nf': nf, 'freq': freq, 'wave': wave, 'cabs': cabs, 'csca': csca, 'nrt': nrtrange, 'trange': trange}

        return opacity


    def writeopac(f, localdust, nr):
        return 0
    def findkappa(localdust, f):
        return 0
