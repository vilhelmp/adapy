
import pyfits, numpy, scipy, math, scipy.ndimage

def copy_header(hdr):

    try:
        del hdr['TELESCOPE']
    except KeyError:
        a = 1
    
    try:
        del hdr['INSTRUMENT']
    except KeyError:
        a = 1
    
    try:
        del hdr['REST FREQUENCY']
    except KeyError:
        a = 1
    
    try:
        del hdr['TRANSITION']
    except KeyError:
        a = 1

    hdr_cardlist = hdr.ascard
    ncards = len(hdr_cardlist)
    
    new_hdr_cardlist = pyfits.CardList()

    for i in range(0, ncards):

        new_hdr_cardlist.append(pyfits.Card( 
                                    copy.deepcopy(hdr_cardlist[i].key), 
                                    copy.deepcopy(hdr_cardlist[i].value), 
                                    copy.deepcopy(hdr_cardlist[i].comment) 
                                            )
                                )
    
    new_hdr = pyfits.Header(new_hdr_cardlist)

    return new_hdr

def smooth_image_xy(hdulist, theta_1, theta_2, mode = 'reflect'):

    fwhmxy = math.sqrt(theta_2**2.0 - theta_1**2.0)

    sigxy = (fwhmxy / (2.0 * math.sqrt(2.0 * math.log(2.0)))) 
            / 
            (numpy.abs(hdulist[0].header['CDELT2']) * 3600.0)

    hdulist_smoothed = pyfits.HDUList(
                                [pyfits.PrimaryHDU(hdulist[0].data)]
                                )
    hdulist_smoothed[0].header = copy_header(hdulist[0].header)

    for i in range(0, hdulist[0].data.shape[0]):

        hdulist_smoothed[0].data[i,:,:] = scipy.ndimage.filters.gaussian_filter(
                                            hdulist[0].data[i,:,:], 
                                            sigxy, 
                                            order = 0, 
                                            output = None, 
                                            mode = mode)
         
    return hdulist_smoothed
