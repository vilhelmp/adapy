#! /usr/bin/env python



def plot_hist(data, bins=50, gaufit=True, **kwargs):
    
    import matplotlib.mlab as mlab
    import scipy.stats
    
    (mu, sigma) = scipy.stats.norm.fit(data)
    
    import matplotlib.pyplot as pl
    pl.subplot(211)
    n, bns, patches = pl.hist(data, normed=True, bins=bins, histtype='stepfilled', label='Data', **kwargs)
    y = mlab.normpdf(bns, mu, sigma)
    l = pl.plot(bns, y, 'r--', linewidth=2, label='Gaussian')
    pl.legend()
    pl.xlabel('Data')
    pl.ylabel('Normalized counts')

    pl.subplot(212)
    color = patches[0].get_fc()
    pl.plot(bns[:-1], (y[:-1]-n), '.', color=color)
    
    return n, bns, patches, y, l
