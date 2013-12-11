
import scipy as _sp

#~ dan = uvfits.uvfitting.uvdata()
#~ dan.read_fits('iras2a_cont.aipsfits')

#~ iras2a = adapy.fits.Uvfits('iras2a_cont.uvfits')
#~ iras2a.load_model('sky_model_sim.uvfits')


# fwhm(image plane) = 182 / fwhm(uv plane)
#       where the image plane fwhm is measured in arcseconds, and the uv
#       plane fwhm is measured in kilowavelengths.
        
#~ from adapy.fitting import mpfit

deg2rad = lambda deg : deg * _sp.pi / 180

rad2deg  = lambda rad : rad * 180 / _sp.pi

asec2lam = 206264.806 # 1 PC in AU, or arcsecond in ???
# 206264.806 is also the number of arcsecs in a circle / 2pi

def phase_pos(uv, dRa, dDec):
    return 2 * _sp.pi * ( uv[0] * dRa / asec2lam + uv[1] * dDec / asec2lam )

def rotate_field(uv, PA, direction = 'CCW'):
    """
    Rotates a coordinate system (UV plane) counter clockwise (CCW) by PA
    degrees or clockwise (CW)
    """
    # The sign between the u and v, determines the direction
    morp = dict([['CCW', -1],[ 'CW', 1]])
    d = morp[direction]
    #~ PA *= d
    u_new = uv[0] * _sp.cos( deg2rad(PA) ) + uv[1] * _sp.sin( deg2rad(PA) )
    v_new = -1 * uv[0] * _sp.sin( deg2rad(PA) ) + uv[1] * _sp.cos( deg2rad(PA) )
    #~ u_new = uv[0] * _sp.cos( deg2rad(PA) ) + d * uv[1] * _sp.sin( deg2rad(PA) )
    #~ v_new = uv[0] * _sp.sin( deg2rad(PA) ) - d * uv[1] * _sp.cos( deg2rad(PA) )
    return u_new, v_new
    
def incline(uv, inc):
    ruv = ( uv[0]**2 + (uv[1] * _sp.cos(deg2rad(inc)) )**2  )**.5 
    return ruv

def uvgauss(amp, size, r):
    """
    Parameters:
    p[0]    p[1]    
    amp     size

    Input:
    r
    coordinate of data point
    
    """
    # Gaussian in the UV plane
    # see Jean-Philippe Berger's "AN INTRODUCTION TO VISIBILITY MODELING"
    # for an explanation of the equation
    # p[3] multiplied by pi/3600 to go over to ??
    # r or p[3] multiplied by pi/180 to go over t? ??
    e = -1 * ( _sp.pi**2 * size * r / (3600. * 180.) )**2 / (4. * _sp.log(2.))
    return amp * _sp.exp(e)

def model(p, uv, env):
    """
    uv = _sp.array([[u],[v]])
    
    env = _sp.array([[env_u],[env_v]])
    
    Parameters:
    p[0]    p[1]    p[2]    p[3]    p[4]    p[5]    p[6]
    dRA     dDec    amp     size    inc     PA      PS
                                                    Point Source
    """
    #p['dRA']   = p[0]
    #p['dDec']  = p[1]
    #p['amp']   = p[2]
    #p['size']  = p[3]
    #p['inc']   = p[4]
    #p['PA']    = p[5]
    #p['PS']    = p[6]
    if p.has_key('inc') and p.has_key('PA'):
        pha = phase_pos( uv, p['dRa'].value, p['dDec'].value)
        uv_new = rotate_field( uv, p['PA'].value )
        ruv = incline( uv_new, p['inc'].value )
    else:
        pha = phase_pos( uv, p['dRa'].value, p['dDec'].value)
        ruv = ( uv[0]**2 + uv[1]**2  )**.5 
    gp = p['PS'].value + uvgauss(p['amp'].value, p['size'].value, ruv)
    # assumes that the point source and the Gaussian have the same
    # position
    gp_reim = gp * _sp.array([ _sp.cos( pha ), _sp.sin( pha )])
    # calculate the total model
    mod = gp_reim + env
    return mod 

def residual_fn(params, uv, vis, env, sigma):
    mod = model(params, uv, env)
    #~ return ((vis[0] - mod[0])**2 + (vis[1] - mod[1])**2) / (3 * sigma)**2
    return ((vis - mod)**2).sum(axis=0) / (3 * sigma)**2

###FITTING
def fit(offset=[0,0], ,
        objects='gaussian,point',
        model=0)
    from lmfit import minimize, Parameters, report_fit
    p_gp = Parameters()

    p_gp.add('dRa',  value=0.02)
    p_gp.add('dDec', value=-0.002)
    p_gp.add('amp',  value=0.3)
    p_gp.add('size',  value=1.8)
    p_gp.add('PS',   value=0.1, min=0, vary=True)

    #~ p_gp.add('inc',  value=14, vary=True)    # Elliptical Gaussian
    #~ p_gp.add('PA',   value=-27, min=-90, max=90)

    #~ p_gp.add('inc',  value=0, vary=False)       # Circular Gaussian
    #~ p_gp.add('PA',   value=0, min=-90, max=90, vary=False) # Circular Gaussian




    uv = _sp.array([iras2a.u_lam, iras2a.v_lam])
    vis = _sp.array([iras2a.re, iras2a.im])
    env = _sp.array([iras2a.Model.re, iras2a.Model.im])
    vis -= env
    env = 0
    sigma = iras2a.sigma

    fitout2 = minimize( residual_fn, p_gp, args=(uv, vis, env, sigma))

    report_fit( fitout2.params )

    return fitout2

import matplotlib.pyplot as plt
plt.ion()

plt.errorbar(iras2a.uvdist_klam, iras2a.amp, yerr=iras2a.sigma, color='b', marker=',', ls='None', label='Data')

plt.plot(iras2a.uvdist_klam, iras2a.Model.amp,'.g', lw=1.5, label='Envelope')

plt.plot([0,iras2a.uvdist_klam.max()], [fitout2.params['PS'].value, fitout2.params['PS'].value],'k', lw=1.5, label='Point Source')

r_array = _sp.arange(0, iras2a.uvdist_klam.max(), 1)
plt.plot(r_array, uvgauss(fitout2.params['amp'].value, fitout2.params['size'].value, r_array*1e3), 'm-', lw=2, label='Gauss')


plt.plot(iras2a.uvdist_klam, uvgauss(fitout2.params['amp'].value, fitout2.params['size'].value, iras2a.uvdist_klam*1e3) + fitout2.params['PS'].value + iras2a.Model.amp, 'y.', lw=1.5, label='All')
#~ plt.plot(iras2a.uvdist_klam, uvgauss(fitout2.params['amp'].value, fitout2.params['size'].value, iras2a.uvdist_klam*1e3) + fitout2.params['PS'].value*2 + iras2a.Model.amp, 'y.', lw=1.5, label='All')

"""
Chisq 110972.457742 113576.742451
dra:  0.0219254467327 0.0038136654398
ddec:  -0.0303820184215 0.00862342369389
flux:  0.288737154888 0.00787698379379
size:  1.47575122688 0.0482558647562
inc:  22.253735034 4.07598166308
pa:  90.6899346751 0.314087006013
point flux [Jy]:  0.0410064387655 0.00349010149164

"""
plt.plot(iras2a.uvdist_klam, uvgauss(0.288737154888, 1.47575122688, iras2a.uvdist_klam*1e3) + 0.0410064387655 + iras2a.Model.amp, '0.35', marker='.', ls='None', lw=1.5, label='DMC')

""" 
 C_GAUSS  R.A.        =     0.02272 (  0.00879)  03:28:55.5718
 C_GAUSS  Dec.        =    -0.00491 (  0.00916)  31:14:37.0951
 C_GAUSS  Flux        =     0.25017 (  0.00306)
 C_GAUSS  F.W.H.P.    =     1.64774 (  0.01730)
 POINT    R.A.        =    -0.02941 (  0.00504)  03:28:55.5677
 POINT    DEC.        =    -0.08751 (  0.00421)  31:14:37.0125
 POINT    FLUX        =     0.04962 (  0.00100)
"""

plt.plot(iras2a.uvdist_klam, uvgauss(0.25017, 1.64774, iras2a.uvdist_klam*1e3) + 0.04962 + iras2a.Model.amp, color='#9999AA', marker='.', ls='None', lw=1, label='Mapping')


plt.legend()

plt.grid()
