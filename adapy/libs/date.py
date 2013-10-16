
def jd2gd(jd):

    """
    
    From http://www.astro.ucla.edu/~ianc/python/_modules/date.html
    
    Task to convert a list of julian dates to gregorian dates
    description at http://mathforum.org/library/drmath/view/51907.html
    Original algorithm in Jean Meeus, "Astronomical Formulae for
    Calculators"
    
    These functions were taken from Enno Middleberg's site of useful
    astronomical python references:
    http://www.astro.rub.de/middelberg/python/python.html

    "Feel free to download, use, modify and pass on these scripts, but
    please do not remove my name from it." --E. Middleberg
    """
    #~ import string
    jd=jd+0.5
    Z=int(jd)
    F=jd-Z
    alpha=int((Z-1867216.25)/36524.25)
    A=Z + 1 + alpha - int(alpha/4)

    B = A + 1524
    C = int( (B-122.1)/365.25)
    D = int( 365.25*C )
    E = int( (B-D)/30.6001 )

    dd = B - D - int(30.6001*E) + F

    if E<13.5:
	mm=E-1

    if E>13.5:
	mm=E-13

    if mm>2.5:
	yyyy=C-4716

    if mm<2.5:
	yyyy=C-4715

    months=["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
    daylist=[31,28,31,30,31,30,31,31,30,31,30,31]
    daylist2=[31,29,31,30,31,30,31,31,30,31,30,31]

    h=int((dd-int(dd))*24)
    mins=int((((dd-int(dd))*24)-h)*60)
    sec=86400*(dd-int(dd))-h*3600-mins*60

    # Now calculate the fractional year. Do we have a leap year?
    if (yyyy%4 != 0):
	days=daylist2
    elif (yyyy%400 == 0):
	days=daylist2
    elif (yyyy%100 == 0):
	days=daylist
    else:
	days=daylist2

    hh = 24.0*(dd % 1.0)
    mins = 60.0*(hh % 1.0)
    sec = 60.0*(mins % 1.0)

    dd =  dd-(dd%1.0)
    hh =  hh-(hh%1.0)
    mins =  mins-(mins%1.0)


    #~ print str(jd)+" = "+str(months[mm-1])+ ',' + str(dd) +',' +str(yyyy)
    #~ print string.zfill(h,2)+":"+string.zfill(mins,2)+":"+string.zfill(sec,2)+" UTC"

    #~ print (yyyy, mm, dd, hh, mins, sec)

    return (yyyy, mm, dd, hh, mins, sec)
