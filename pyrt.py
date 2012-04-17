#! /usr/bin/env python
# -*- coding: utf-8 -*-


#
#   pyrt.py
#
#
#   Copyright 2012 Magnus Persson <magnusp@nbi.dk>
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#
#   ver 0.1 alpha
#
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
#------------------------------------------------------------------------
#Top of the list TODO:
"""
TODO :
        RADEX:
            - Ability to run grids
        RATRAN
            - Initialisation routine

        Transphere
            - Prepare input and run Transphere
            - Optimize the writeTransphereInput

"""

########################################################################
# MODELING HELP FUNCTIONS
def bplanck(nu,T): # Returns the Spectral Radiance (Planck)
    from scipy import exp, constants
    import cgsconst as cgs
    try:
        print('Dunno if it is correct yet, check!')
        x = cgs.HH*nu/(cgs.KK*T)
        bpl = (2*cgs.HH*nu**3/cgs.CC**2)/(exp(x-1))
        return bpl
    except (ZeroDivisionError):
        return 0



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
        import adavis
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
                    raise adavis.ParError(par)
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

    def radial_profile(self, spaced='log10',
                            nshell=200, rin=200, rout=8000,
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


class Model:
    def __init__(self, silent=True, **kwargs):
        """
        General model object

        All input parameters have default values

        Input Parameters:
        rin      : Inner radius of shell
        rout     : Outer radius of shell
        nshell   : Number of shells
        spacing  : Type of spacing log10, powerlaw1, linear or powerlaw2
        nref     : Number of refinement shells (not implemented yet)
        rref     : Refinement radius (not implemented yet)
        mstar    : Mass of star
        tstar    : Stellar temperature
        isrf     : Scaling of ISRF
        tbg      : Spectral shape of ISRF (Blackbody equivalent temperature).
                   With -1 the spectrum is read from the file isrf.inp and
                   scaled by ISRF.
        dpc      : Distance in pc
        silent   : Verbose?

        r0       : Reference radius (AU)
        plrhp    : rho exponen in the power
        rho_type : rho dependence on radius
        n0       : H2 number density at reference radius (cm-3)

        ----------------------------------------------------------------

        ----------------------------------------------------------------
        TODO:
        - change the name of rho_type parameter
        """



        # imports
        from scipy import log10, log, arange, linspace, logspace,\
        array
        import cgsconst as cgs
        #
        self.silent = silent

        #

        # Checking input parameters
        params = array(['rin', 20, 'rout', 8000, 'nshell', 200,
                        'spacing', 'powerlaw1' ,'nref', 0 ,'rref', 0,
                        'rstar', 3, 'tstar', 5780, 'mstar', 1,
                        'isrf', 0.0, 'tbg', 2.73, 'dpc', 250,
                        'r0', 1.0E3, 'plrho', -1.5,
                        'rho_type', 'powerlaw1', 'n0', 2e6])
        print ('Model created with the following parameters:')
        for par,stdval in zip(params[0::2],params[1::2]):
            if par in kwargs:
                if par in ['spacing', 'rho_type']:
                    kwargs[par] = '\"{0}\"'.format(kwargs[par])
                print '   {0:7} : {1}'.format(par, kwargs[par].strip('\"'))
                exec('self.{0} = {1}'.format(par, kwargs[par]))
            else:
                if stdval != None:
                    print '   {0:7} : {1:10} (default)'.format(par, stdval)
                    if par in ['spacing', 'rho_type']:
                        exec('self.{0} = {1}'.format(par, '\"'+stdval+'\"'))
                    else:
                        exec('self.{0} = {1}'.format(par, stdval))
                else:
                    raise Exception('Wrong parameter input/handling'
                                    'Check input and/or code...')

        """#~ if not self.nref:
            #~ r = self.rin * (self.rout/self.rin)**(arange(self.nshell)/(self.nshell-1.e0)) # Simple Log grid
        #~ elif self.nref:
                #~ #
                #~ # Log grid with refinement at inner edge
                #~ #
                #~ lgr0 = log10(self.rin)
                #~ lgr1 = log10(self.rout)
                #~ lgrr = log10(self.rref)
                #~
                #~ n1 = self.nref
                #~ n2 = self.nshell - self.nref
                #~ n = arange(self.nshell)-n1
                #~ c = (lgr1-lgrr)/(1.e0*n2)
                #~ b = c/(lgrr-lgr0)
                #~
                #~ r = zeros((self.nshell),float)
                #~
                #~ r[0:n1] = (lgrr-lgr0)*e**(b*n[0:n1])
                #~ r[n1:n1+n2] = (lgrr-lgr0)+c*n[n1:n1+n2]
                #~ r = r-(r[0]+(r[0]-r[1]))
                #~ r = r*(lgr1-lgr0)/(max(r))+lgr0
                #~ r = 10.0**r
        #~ self.r_jes = r

        #~ from scipy import log10, log, arange, linspace, logspace,\
        #~ array
        #from cgsconst import *"""

        if self.spacing.lower() == 'log10':
            # get the exponent of the start- and
            # stop-radius in input units
            start = [log10(self.rin), 0][self.rin == 0]
            stop = log10(self.rout)
            radii = logspace(start, stop, num=self.nshell, endpoint=True)
        elif self.spacing.lower() == 'powerlaw1':
            if not self.nref and not self.rref:
                radii = self.rin * (self.rout/self.rin)**(arange(self.nshell)/(self.nshell-1.e0))
        elif self.spacing.lower() == 'linear':
            # linearly spaced grid
            radii = linspace(self.rin, self.rout, num=self.nshell, endpoint=True)
        elif self.spacing.lower() == 'powerlaw2':
            # first check if coefficients to the power-law was given
            #~ if 'exp' in kwargs:
                #~ p_exp = kwargs['exp']
            #~ else: # if not, set it to 2, i.e. r^2
                #~ p_exp = 2

            radii = self.rin + (self.rout - self.rin)*(linspace(self.rin, self.rout, num=self.nshell, endpoint=True)/(self.rout))**2
            #print('Not implemented yet.')
            #raise ParError(spaced)
        else:
            raise Exception(spaced)
        # convert to cm
        # now to get the
        lower = radii[:-1]
        upper = radii[1:]
        self.radii = radii
        self.r_midpt = array((lower+upper)/2)
        self.r_grid = array([array([i,j]) for i,j in zip(lower,upper)])


        #
        # Create the density array
        #

        from scipy import array
        import cgsconst as cgs
        # Checking input parameters
        params = array([])
        print ('Model created with the following parameters:')
        for par,stdval in zip(params[0::2],params[1::2]):
            if par in kwargs:
                if par == 'rho_type':
                    kwargs['rho_type'] = '\"{0}\"'.format(kwargs['rho_type'])
                print '   {0:8} : {1}'.format(par, kwargs[par].strip('\"'))
                exec('self.{0} = {1}'.format(par, kwargs[par]))
            else:
                if stdval != None:
                    print '   {0:8} : {1:10} (default)'.format(par, stdval)
                    if par == 'rho_type':
                        exec('self.{0} = {1}'.format(par, '\"'+stdval+'\"'))
                    else:
                        exec('self.{0} = {1}'.format(par, stdval))
                else:
                    raise Exception('Wrong parameter input/handling'
                                    'Check input and/or code...')

        # calculate the density at the reference radius
        # number of H2 * 2.3 * 1.67E-27
        self.rho0 = self.n0 * cgs.MUH2 * cgs.MP # in g/cm-3?
        # calculate the density at all radial points
        if self.rho_type == 'powerlaw1':
            self.rho    = 1e-2 * self.rho0 * (self.radii / self.r0)**(self.plrho)
        else:
            raise Exception('rho_type parameter not recognised')



    def __str__(self):
        attributes = ['rin', 'rout', 'nshell',
                        'spacing', 'nref','rref',
                        'rstar', 'tstar',
                        'isrf', 'tbg', 'dpc']
        test = '\n'.join(['{0:7} : {1}'.format(i, getattr(self, i, None)) for i in attributes])
        return test

    #~ def density(self, **kwargs):


class Transphere:
    """
    Class to interface with Transphere - dust continuum
    radiative transfer code.
    """
    def __init__(self, ModelObject, **kwargs):
        # imports
        from scipy import log10, log, arange, linspace, logspace,\
        array
        #
        class Opacity: pass
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
        params = array(['localdust', 0, 'silent', True,
                        'nriter', 30, 'convcrit', 1E-5,
                        'ncst', 10,'ncex', 30,
                        'ncnr', 1,'itypemw', 1,
                        'idump', 1])
        #~ print ('Model created with the following parameters:')
        for par,stdval in zip(params[0::2],params[1::2]):
            if par in kwargs:
                #~ if par == 'spacing':
                    #~ kwargs['spacing'] = '\"{0}\"'.format(kwargs['spacing'])
                print '   {0:10} : {1}'.format(par, kwargs[par].strip('\"'))
                exec('self.{0} = {1}'.format(par, kwargs[par]))
            else:
                if stdval != None:
                    print '   {0:10} : {1:10} (default)'.format(par, stdval)
                    #~ if par == 'spacing':
                        #~ exec('Model.{0} = {1}'.format(par, '\"'+stdval+'\"'))
                    #~ else:
                    exec('self.{0} = {1}'.format(par, stdval))
                else:
                    raise Exception('Wrong parameter input/handling'
                                    'Check input and/or code...')
        self.Model = ModelObject

        """#~ if not self.Model.nref:
            #~ r = self.Model.rin * (self.Model.rout/self.Model.rin)**(arange(self.Model.nshell)/(self.Model.nshell-1.e0)) # Simple Log grid
        #~ elif self.Model.nref:
                #~ #
                #~ # Log grid with refinement at inner edge
                #~ #
                #~ lgr0 = log10(self.Model.rin)
                #~ lgr1 = log10(self.Model.rout)
                #~ lgrr = log10(self.Model.rref)
                #~
                #~ n1 = self.Model.nref
                #~ n2 = self.Model.nshell - self.Model.nref
                #~ n = arange(self.Model.nshell)-n1
                #~ c = (lgr1-lgrr)/(1.e0*n2)
                #~ b = c/(lgrr-lgr0)
                #~
                #~ r = zeros((self.Model.nshell),float)
                #~
                #~ r[0:n1] = (lgrr-lgr0)*e**(b*n[0:n1])
                #~ r[n1:n1+n2] = (lgrr-lgr0)+c*n[n1:n1+n2]
                #~ r = r-(r[0]+(r[0]-r[1]))
                #~ r = r*(lgr1-lgr0)/(max(r))+lgr0
                #~ r = 10.0**r
        #~ self.Model.r_jes = r

        #~ from scipy import log10, log, arange, linspace, logspace,\
        #~ array
        #from cgsconst import *"""

        #~ if Model.spacing.lower() == 'log10':
            #~ # get the exponent of the start- and
            #~ # stop-radius in input units
            #~ start = [log10(Model.rin), 0][Model.rin == 0]
            #~ stop = log10(Model.rout)
            #~ radii = logspace(start, stop, num=Model.nshell, endpoint=True)
        #~ elif Model.spacing.lower() == 'powerlaw1':
            #~ if not Model.nref and not Model.rref:
                #~ radii = Model.rin * (Model.rout/Model.rin)**(arange(Model.nshell)/(Model.nshell-1.e0))
        #~ elif Model.spacing.lower() == 'linear':
            #~ # linearly spaced grid
            #~ radii = linspace(Model.rin, Model.rout, num=Model.nshell, endpoint=True)
        #~ elif Model.spacing.lower() == 'powerlaw2':
            # first check if coefficients to the power-law was given
            #~ if 'exp' in kwargs:
                #~ p_exp = kwargs['exp']
            #~ else: # if not, set it to 2, i.e. r^2
                #~ p_exp = 2

            #~ radii = Model.rin + (Model.rout - Model.rin)*(linspace(Model.rin, Model.rout, num=Model.nshell, endpoint=True)/(Model.rout))**2
            #print('Not implemented yet.')
            #raise ParError(spaced)
        #~ else:
            #~ raise Exception(spaced)

        # now to get the
        #~ lower = radii[:-1]
        #~ upper = radii[1:]
        #~ Model.r_midpt = array((lower+upper)/2)
        #~ Model.r_grid = array([array([i,j]) for i,j in zip(lower,upper)])
        #~ #
        #~ self.Model = Model


    #~ def density(self, **kwargs):
        #~ # Checking input parameters
        #~ params = array([r0, '1.0E3', n0])
        #~ print ('Model created with the following parameters:')
        #~ for par,stdval in zip(params[0::2],params[1::2]):
            #~ if par in kwargs:
                #~ if par == 'spacing':
                    #~ kwargs['spacing'] = '\"{0}\"'.format(kwargs['spacing'])
                #~ print '   {0:7} : {1}'.format(par, kwargs[par].strip('\"'))
                #~ exec('Model.{0} = {1}'.format(par, kwargs[par]))
            #~ else:
                #~ if stdval != None:
                    #~ print '   {0:7} : {1:10} (default)'.format(par, stdval)
                    #~ if par == 'spacing':
                        #~ exec('Model.{0} = {1}'.format(par, '\"'+stdval+'\"'))
                    #~ else:
                        #~ exec('Model.{0} = {1}'.format(par, stdval))
                #~ else:
                    #~ raise Exception('Wrong parameter input/handling'
                                    #~ 'Check input and/or code...')

    def fix_opacity_file(self, **kwargs):
        from scipy import array
        from sys import exit as sysexit
        # if a dust opacity file is given
        # that is it is not tabulated correctly
        # as a dustopac.inp file
        print 'Local dust opacity: ',self.localdust
        if 'opacfile' in kwargs:
            if 'freqfile' not in kwargs:
                sysexit('No frequency file given. Cannot proceed. '
                'Need two files: one with frequency, and one with opacity'
                ' at corresponding frequency.')
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
                    cabs = array([[float(i) for i in j.split()] for j in lines[2 : nf+2]])
                    csca = array([[float(i) for i in j.split()] for j in lines[nf+2 : 2*nf+2]])
                except:
                    errmsg = 'Error parsing Opacity-file. Wrong format.'
                    raise Exception(errmsg)
                nrtrange=1
                trange=[0.e0,0.e0]
            else:
                ########################################################
                # this part not edited, dont have a file like this
                sysexit('This part of the code is not tested'
                        'Please optimize code and check it.')
                nf = ns
                ismooth = 0
                nrtrange = 0
                with open(kwargs['opacfile']) as f:
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
            with open(kwargs['freqfile']) as f:
                lines = f.read().split('\n')
            nf = int(lines[0])
            freq = array(lines[2:2+nf], 'float')
            wave = 2.9979e14 / freq

            ## TODO : change this, do I even need a dictionary?
            # move the calls to where the variable is defined
            self.ns = ns
            self.nf = nf
            self.freq = freq
            self.freq = freq
            self.wave = wave
            self.cabs = cabs
            self.csca = csca
            self.nrtrange = nrtrange
            self.trange = trange

            #~ self.opacity = {'ns': ns, 'nf': nf, 'freq': freq, 'wave': wave, 'cabs': cabs, 'csca': csca, 'nrt': nrtrange, 'trange': trange}

        else:
            # always neeed a opacity file and a frequency file
            sysexit('No opacity-file given.')
        ################################################################
        import os
        # if local dust opacity
        # the opacity is tabulated with radius
        #
        # have not got a file like this, so i am not changing it
        if self.localdust:
            sysexit('This part of the code is not tested'
                    'Please optimize code and check it.')
            if 'nr' not in kwargs:
                sysexit('Number of radial shells, \"nr\" not given as input.')
            os.system('rm -f dustopac_1.inp')
            os.system('rm -f dustopac.inp') # added this, do I really need
                                            # to remove it too?
            f = open('dustopac.inp','w')
            f.write('1               Format number of this file'+'\n')
            f.write('1               Nr of dust species'+'\n')
            f.write('============================================================================'+'\n')
            f.write('-1              Way in which this dust species is read (-1=file)'+'\n')
            f.write('0               Extension of name of dustopac_***.inp file'+'\n')
            f.write('----------------------------------------------------------------------------'+'\n')
            f.close
            f = open('dustopac_0.inp','w')
            f.write(nr)
            f.write(str(opacity['nf'])+' 1\n')
            f.write(' ')
            redux=1.e0
            for ir in range(0,nr):
                for inu in range(0,nr):
                    f.write(opacity['cabs'][inu]*redux)
                for inu in range(0,opacity['nf']):
                    f.write(opacity['csca'][inu]*redux)
                f.write(' ')
            f.close
        elif not self.localdust:
            # first remove the standard ratran dust opacity input files
            os.system('rm -f dustopac_0.inp')
            os.system('rm -f dustopac.inp') # added this, do I really need
                                            # to remove it too?
            with open('dustopac.inp','w') as f:
                f.write('1               Format number of this file\n')
                f.write('1               Nr of dust species\n')
                f.write('============================================================================\n')
                f.write('-1              Way in which this dust species is read (-1=file)\n')
                f.write('1               Extension of name of dustopac_***.inp file\n')
                f.write('----------------------------------------------------------------------------\n')
            with open('dustopac_1.inp','w') as f:
                f.write(str(nf)+' 1\n \n')
                for inu in range(0, nf):
                    f.write(str(cabs[inu][0])+'\n')
                for inu in range(0, nf):
                    f.write(str(csca[inu][0])+'\n')

    def writeTransphereInput(self):
        #import natconst as nc
        #~ import math
        #~ import astroProcs
        #~ import numpy as np
        from scipy import pi, zeros
        import cgsconst as cgs
        # Transphere input file
        f=open('transphere.inp','w')
        f.write(str(2)+'\n')
        f.write(str(self.nriter)+'\n')
        f.write(str(self.convcrit)+'\n')
        f.write(str(self.ncst)+'\n')
        f.write(str(self.ncex)+'\n')
        f.write(str(self.ncnr)+'\n')
        f.write(str(self.itypemw)+'\n')
        f.write(str(self.idump)+'\n')
        f.close()
        #
        # Make the stellar information file
        # (mstar and tstar are irrelevant; they are there for historical reasons)
        #
        f=open('starinfo.inp','w')
        f.write(str(1)+'\n')
        f.write(str(self.Model.rstar)+'\n')
        f.write(str(self.Model.mstar)+'\n')
        f.write(str(self.Model.tstar)+'\n')
        f.close()
        #
        # The stellar spectrum
        #
        f=open('starspectrum.inp','w')
        f.write(str(len(self.freq))+'\n')
        sspec=(self.Model.rstar/cgs.PC)**2*pi*bplanck(self.freq,self.Model.tstar)
        for inu in range(0,len(self.freq)):
            f.write(str(self.freq[inu])+' '+str(sspec[inu])+'\n')
        f.close()
        #
        # The exterior spectrum
        #
        if self.Model.tbg == 0.0 or self.Model.isrf == 0:
            bgspec = zeros((len(self.freq)),float)
        elif self.Model.tbg == -1:
            f = open('isrf.inp','r')
            nf = int(f.readline().strip())
            bgspec = zeros((len(self.freq)),float)
            for ii in range(0,nf):
                bgspec[ii]=float(f.readline().strip())*self.Model.isrf
        else:
            if self.Model.tbg > 0: bgspec = bplanck(self.freq, self.Model.tbg)*self.Model.isrf

        f = open('external_meanint.inp','w')
        f.write(str(len(self.freq))+'\n')
        for inu in range(0,len(self.freq)):
            f.write(str(self.freq[inu])+' '+str(bgspec[inu])+'\n')
        f.close()
        #
        # Write the envelope structure
        #
        f = open('envstruct.inp','w')
        f.write(str(len(self.Model.radii))+'\n')
        f.write(' '+'\n')
        for ir in range(0,len(self.Model.radii)):
            f.write("%13.6E %13.6E %13.6E" % (self.Model.radii[ir], self.Model.rho[ir],0.e0)+'\n') # ,format='(3(E13.6,1X))'
        f.close()

    def runTransphere(self):
        import os
        from time import time
        #
        # check the input files here!
        #
        #~ try:
        print ('Running Transphere...')
        t1 = time()
        os.system('transphere')
        t2 = time()
        print('Done, took : {0:3.3f} seconds'.format((t2-t1)))
        #~ except Exception

        # after run read in the results automatically, and store
        # in object

    ####################################################################
