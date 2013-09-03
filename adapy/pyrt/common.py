
### Imports
# python builtins
import os as _os
import sys as _sys
import subprocess as _subprocess



def test_thick(abund, flux, removelast = 1):
    """
        Test the linearity of some values.
    """
    from scipy import polyfit
    from scipy import arange
    print('removelast={0} points in fit'.format(removelast))
    p = polyfit(abund[:-1*int(removelast)], flux[:-1*int(removelast)], deg=1)
    line = lambda x: p[0]*x + p[1]    
    import matplotlib.pyplot as pl
    pl.ion()


    off = (line(abund[-1]) - flux[-1])/line(abund[-1])*100
    print('last value off by : {0:.2f}% from fit'.format(off))
    
    pl.plot(abund, flux, 'or', label='data')
    X = arange(abund[0]*0.9, abund[-1]*1.1, abund[0]*0.1)
    
    pl.plot(X, line(X),'-g', label='fit')
    pl.legend(loc=4)
