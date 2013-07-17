
def get_qrot():
    q = array([1.26,  3.054, 8.648, 23.361, 64.21, 117.004, 179.639])
    # new values from May 2011 by Brian Druin from JPL
    # for para water only! does not seem reliable?
    #~ q = array([1.0108, 1.2060, 2.2700, 5.8434, 16.0522, 29.2468, 44.8740])
    t = array([9.375, 18.75, 37.5, 75, 150, 225, 300])
    x = arange(t.min(),t.max(),10)
    #plot(t,q,'xr',label='data', lw=2,mew=5)
    #xlabel(r'T$_{ex}$')
    #ylabel(r'Q$_{rot}$')
    # with a Qrot =a + b*T**c from G. Busquet et al 2010.: Deuterated ammonia in IRAS 20293+3952
    from scipy.optimize import leastsq
    fitfunc = lambda p,t: q - (p[0]*t**p[1])
    p0 = [0.05,1.5]
    p = leastsq(fitfunc,p0,args=(t))[0]
    qrot = lambda p,t: p[0]*t**p[1]
    
    return qrot
    
