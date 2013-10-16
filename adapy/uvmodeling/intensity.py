




## Definition of intensity distributions to fit
def Gaussian(p,x,binsize):
    # Gaussian Distribution
    newx = x[0]*cos(p[3]) + x[1]*sin(p[3])
    newy = x[1]*cos(p[3]) - x[0]*sin(p[3])
    newr = sqrt(newx*newx  + (newy*newy*cos(p[2])*cos(p[2])))
    #print x[0].min()*1e-3, x[0].max()*1e-3
    #print x[1].min()*1e-3, x[1].max()*1e-3
    #dist = sqrt(x[0]*x[0] + x[1]*x[1])*1e-3
    #print dist.min(), dist.max()
    #sys.exit()
    # This gives the correct klam distribution!
    ruv = sqrt(x[0]*x[0] + x[1]*x[1])
    vis = p[0]*exp(-1.0*(pi*p[1]/3600.0*(pi/180.0)*newr)*(pi*p[1]/3600.0*(pi/180.0)*newr)/(4.0*log(2.0)))
    vis1 = bin_mod(ruv*1e-3,vis,binsize=binsize,uvmin=0.0,uvmax=350.0)# bin in klam
    isub = (isnan(vis1)).nonzero()[0]
    if len(isub) > 1:
        vis1[isub] = 0.0
    # Find the zeros in the first few bins
    isub = (vis1 == 0.0).nonzero()[0]
    isub1 = (isub < 40).nonzero()[0]
    vis1[isub[isub1]] = p[0]
    return vis1[:,1]
