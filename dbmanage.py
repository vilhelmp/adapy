import sys
import pyfits as pf
from scipy import *
import MySQLdb as mdb
from datetime import date
import pywcs
"""
Plan:
- First create function to manually add data that works good. Then automate it.
- Create database model in some modelling software to practice and get a hang of it.


What to do:

TODO :




"""
def db_connect():
    # connect to the server
    try: conn = mdb.connect(host = "webserver",
                            user = "root",
                            passwd = "")
    except mdb.Error, e:
        cursor.close(); conn.close()
        raise Exception("MySQL error %d: %s"% (e.args[0], e.args[1]))
    return conn

def db_query(askme, db='fitsfiles', plot=False):
    # connect to the server
    conn = db_connect()
    cursor = conn.cursor()
    try:
        cursor.execute ("USE "+db)
        cursor.execute (askme)
    except mdb.Error, e:
        cursor.close(); conn.close()
        raise Exception("MySQL error %d: %s"% (e.args[0], e.args[1]))
    row = cursor.fetchall ()
    for x in xrange(0,len(row)):
        line = [str(i) for i in row[x]]
        print("\t".join(line))

    cursor.close(); conn.close()
    if plot==True:
        import matplotlib.pyplot as pl
        filename = row[0][8]
        path = row[0][9]
        f = pf.open(path+filename)
        hdr, d = f[0].header, f[0].data
        d = d[0]  # axis=0 STOKES, 1 Vel/Freq/Chan, 2 Y (DEC), 3 X (RA)
        f.close()

        veltype = str([x for x in hdr.keys() if x[:-1]=='CTYPE' and hdr[x]=='VELO-LSR'][0][-1:])
        vcrpix = hdr['CRPIX'+veltype]
        vcrval = hdr['CRVAL'+veltype]
        vctype = hdr['CTYPE'+veltype]
        vcdelt = hdr['CDELT'+veltype]
        vnaxis = hdr['NAXIS'+veltype]
        vcdeltkms = vcdelt/float(1e3)

        velarr = ((arange(1,vnaxis+1)-vcrpix)*vcdelt+vcrval)/float(1e3)

        dec_cdelt = hdr['CDELT2']*3600
        dec_npix = hdr['NAXIS2']
        y_npix = hdr['NAXIS2']

        ra_cdelt = hdr['CDELT1']*3600
        ra_npix = hdr['NAXIS1']
        x_npix = hdr['NAXIS1']

        x1,x2,y1,y2  = 110,140,70,190
        data = d[:,y1:y2,x1:x2].sum(axis=1).sum(axis=1)
        pl.step(velarr,data,'k')
        xmin, xmax = round(velarr[0]-1),round(velarr[-1]+1)
        ymin, ymax = round(min(data)-10),round(max(data)+10)
        pl.plot([xmin,xmax],[0,0],'k:',[0,0],[ymin,ymax],'k:')
        pl.xlim(xmin,xmax); pl.ylim(ymin,ymax)
        pl.show()




def db_create(what):
    """
    Usage "create_database(name_of_new_database)"
    Creates a database within the standard MySQL server.
    """
    db = 'fitsfiles'
    tbl = 'files'
    """
    Create database
    """
    if what=='db':
        # connect to the server
        print('Creating database '+db)
        conn = db_connect()
        cursor = conn.cursor ()
        try: cursor.execute ("CREATE DATABASE "+db+";")
        except mdb.Error, e:
            cursor.close(); conn.close()
            raise Exception("MySQL error %d: %s"% (e.args[0], e.args[1]))
        cursor.execute ("SHOW DATABASES;")
        row = cursor.fetchall ()
        if (db,) in row: print("Database \'"+db+"\' created.")
    """
    Create table
    """
    if what=='tbl':
        print('Creating table '+tbl)
        conn = db_connect()
        cursor = conn.cursor ()
        try:
            cursor.execute ("USE "+db+";")
            # the fields of the file-table
            # assumes that axis 3 is the frequency/velocity axis

            fields= "object varchar(40),\
            center_ra double(40,15) NOT NULL,\
            center_dec double(40,15) NOT NULL,\
            date_obs date,\
            equinox int(6) NOT NULL,\
            naxis int(1) NOT NULL,\
            vobs double(40,15),\
            observer varchar(100),\
            telescope varchar(50),\
            date_added date NOT NULL,\
            filename varchar(50) NOT NULL,\
            path varchar(100) NOT NULL,\
            extra varchar(200),\
            frequency double(40,15),\
            bandwidth double(40,15),\
            freq_res double(40,15),\
            cdelt1 double(40,15) NOT NULL,\
            crpix1 double(40,15) NOT NULL,\
            crval1 double(40,15) NOT NULL,\
            ctype1 varchar(40) NOT NULL,\
            naxis1 int(5) NOT NULL,\
            cdelt2 double(40,15) NOT NULL,\
            crpix2 double(40,15) NOT NULL,\
            crval2 double(40,15) NOT NULL,\
            ctype2 varchar(40) NOT NULL,\
            naxis2 int(5) NOT NULL,\
            cdelt3 double(40,15) NOT NULL,\
            crpix3 double(40,15) NOT NULL,\
            crval3 double(40,15) NOT NULL,\
            ctype3 varchar(40) NOT NULL,\
            naxis3 int(5) NOT NULL"
            cursor.execute ("CREATE TABLE "+tbl+" ("+fields+");")
        except mdb.Error, e:
            cursor.close(); conn.close()
            raise Exception("MySQL error %d: %s"% (e.args[0], e.args[1]))
        row = cursor.fetchall ()
        for x in [x[0] for x in row]:
            print x
        print 'Created table \"',tbl, '\" in database \"', db,'\"'

    cursor.close(); conn.close()

def db_list():
    """
    Usage "db_list()"
    Lists all databases within the standard MySQL server.
    """
    # connect to the server
    conn = db_connect()
    cursor = conn.cursor()
    try: cursor.execute ("SHOW DATABASES;")
    except mdb.Error, e:
        cursor.close(); conn.close()
        raise Exception("MySQL error %d: %s"% (e.args[0], e.args[1]))
    row = cursor.fetchall ()
    for x in [x[0] for x in row]:
        print x
    cursor.close(); conn.close()

def db_describe(db='fitsfiles', tbl='files'):
    """
    Usage
    """
    # connect to the server
    conn = db_connect()
    cursor = conn.cursor ()
    try:
        cursor.execute('USE '+db+';')
        cursor.execute ("DESCRIBE "+tbl+";")
    except mdb.Error, e:
        cursor.close(); conn.close()
        raise Exception("MySQL error %d: %s"% (e.args[0], e.args[1]))
    row = cursor.fetchall ()
    print("Field\t\t Type\t\t Null\t Key\t Default\t Extra\t")
    print(70*"=")
    for x in xrange(0,len(row)):
        line = [str(i) for i in row[x]]
        if len(line[0])<=7:
            line[0] += "\t"
        if len(line[1])<=6:
            line[1] += "\t"
        print("\t".join(line))
    cursor.close(); conn.close()
    #return row



def db_add_fits(path, filename, db='fitsfiles', tbl='files'):
    """
    test
    """
    #check if it has a trailing slash
    if path[-1] != '/': path += '/'
    try:
        f = pf.open(path+filename)
        hdr, d = f[0].header, f[0].data
        f.close()
    except IOError, e:
        print "PyFits Error:", e
        sys.exit(1)
    axes = [x for x in hdr.keys() if x[:-1]=='CTYPE']

    raax = str([x for x in axes if 'RA' in hdr[x]][0][-1:])
    decax = str([x for x in axes if 'DEC' in hdr[x]][0][-1:])

    velax = str([x for x in axes if 'VELO' in hdr[x]][0][-1:])

    vcrpix    = hdr['CRPIX'+velax]
    vcrval    = hdr['CRVAL'+velax]
    vctype    = hdr['CTYPE'+velax]
    vcdelt    = hdr['CDELT'+velax]
    vnaxis    = hdr['NAXIS'+velax]
    if len(axes) == 4:
        stokax = str([x for x in axes if hdr[x]=='STOKES'][0][-1:])

    restfreq = hdr['RESTFREQ']
    # first and last velocity in cube
    v1 = vcrval-vcrpix*vcdelt
    v2 = (vnaxis-vcrpix)*vcdelt
    vedges = array([v1,v2])
    # first and last frequency in cube
    fedges = (1-vedges/3e8)*restfreq

    #center_ra = hdr['CRVAL'+raax]
    #center_dec = hdr['CRVAL'+decax]
    #hwidth_ra = (hdr['NAXIS'+raax]*abs(hdr['CDELT'+raax]))/2 # what about the edge?
    #hwidth_dec = (hdr['NAXIS'+decax]*abs(hdr['CDELT'+decax]))/2


    object = hdr['OBJECT']
    wcsobj = pywcs.WCS(hdr)
    center_ra, center_dec = wcsobj.wcs_pix2sky(array([[hdr['NAXIS'+raax]/2, hdr['NAXIS'+decax]/2,0,0]]),1)[0,:2]
    date_obs = hdr['DATE-OBS']
    try:
        equinox = hdr['EQUINOX']
    except KeyError, e:
        try:
            equinox = hdr['EPOCH']
        except KeyError, e:
            print "PyFits Error:", e
            print "No equinox keyword, exiting..."
            raise Exception("FITS header keyword missing.")

    naxis = hdr['NAXIS']
    vobs = hdr['VOBS']
    observer = hdr['OBSERVER']
    telescope = hdr['TELESCOP']
    date_added = date.isoformat(date.today())
    filename = filename
    path = path
    extra = 'NULL'
    frequency = mean(fedges)
    bandwidth = abs(diff(fedges)[0])/2
    freq_res = bandwidth/(vnaxis/2)
    cdelt1 = hdr['CDELT1']
    crpix1 = hdr['CRPIX1']
    crval1 = hdr['CRVAL1']
    ctype1 = hdr['CTYPE1']
    naxis1 = hdr['NAXIS1']

    cdelt2 = hdr['CDELT2']
    crpix2 = hdr['CRPIX2']
    crval2 = hdr['CRVAL2']
    ctype2 = hdr['CTYPE2']
    naxis2 = hdr['NAXIS2']

    cdelt3 = hdr['CDELT3']
    crpix3 = hdr['CRPIX3']
    crval3 = hdr['CRVAL3']
    ctype3 = hdr['CTYPE3']
    naxis3 = hdr['NAXIS3']


    fields = [
        object,
        center_ra,
        center_dec,
        date_obs,
        equinox,
        naxis,
        vobs,
        observer,
        telescope,
        date_added,
        filename,
        path,
        extra,
        frequency,
        bandwidth,
        freq_res,
        cdelt1,
        crpix1,
        crval1,
        ctype1,
        naxis1,
        cdelt2,
        crpix2,
        crval2,
        ctype2,
        naxis2,
        cdelt3,
        crpix3,
        crval3,
        ctype3,
        naxis3,
    ]
    out = "\',\'".join([str(x) for x in fields])
    out = "(\'"+out+"\')" # add beginning and trailing (' ')

    # adding file to database
    # connect to database
    conn = db_connect()
    cursor = conn.cursor ()

    try:
        cursor.execute ("USE "+db+";")
        cursor.execute ("INSERT INTO "+tbl+" VALUES "+out+";")
    except mdb.Error, e:
        cursor.close(); conn.close()
        raise Exception("MySQL error %d: %s"% (e.args[0], e.args[1]))
    row = cursor.fetchall ()
    for x in [x[0] for x in row]:
        print x
    cursor.close(); conn.close()
    print 'Added file \"',filename, '\" in table \"', tbl,'\"'



path = '/home/magnusp/work/data/db/n1333-i4a/sma'
filename = 'iras4a.13co.line.fits'
filename = 'iras4a.12co.line.fits'
path = '/home/magnusp/work/data/db/n1333-i2a/sma'
filename = 'iras2a.12co.line.fits'

#db_query(dbname, 'SELECT * FROM files WHERE frequency LIKE \'230%\' AND object LIKE \'iras4a\';')
