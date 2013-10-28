#
# Magnus Persson < vilhelm.nu >
#
# Script to create a moldata file for (o/p)pH2(18)O-H2, combining
# the JPL catalog file for (o/p)H2(18)O and the existing pH2O-H2 moldata file.
#
# Converted from a script by Suzanne Bisschop for CH3OCH3
#

from scipy import *
import scipy.constants as C
jplpath = '/home/magnusp/work/data/molecules/H2-(18)-O/'
moldatapath = '/home/magnusp/work/data/moldata/'


# Para ('0') or Ortho ('1') water?
ortho = '1'
# Only transitions with v = 0, i.e. only the rotational transitions
# and not the vibrational ones.
v0orv1 = '0'
# it only checks this if you input the combined file
# because it is not included in the specific ortho/para files

# if jpldata is the new files, with 'o' or 'p' appended, then this check is
# not necessary, if so set 'docheck' to False
# if you use the combined file, set to True
docheck = True

# either:
#~ jpldata = 'c020003.cat' # combined o/p H2-18O file from JPL molecular catalog 
# or:
#~ jpldata = ['c020003p.cat', 'c020003o.cat'][int(ortho)] # H2-18O file from JPL molecular catalog 
jpldata = ['c020003.cat', 'c020003.cat'][int(ortho)] # H2-18O file from JPL molecular catalog 
# from:
# http://spec.jpl.nasa.gov/ftp/pub/catalog/catdir.html
# http://spec.jpl.nasa.gov/ftp/pub/catalog/c020003.cat # combined ortho/para
# http://spec.jpl.nasa.gov/ftp/pub/catalog/c020003o.cat # ortho separate
# http://spec.jpl.nasa.gov/ftp/pub/catalog/c020003p.cat # para separate

# LAMBDA moldata file(s)
moldatareference = ['ph2o-h2@daniel.dat' , 'oh2o-h2@daniel.dat'][int(ortho)]
# http://home.strw.leidenuniv.nl/~moldata/H2O.html

# name of the new moldata file
newmoldatafile = ['ph2-18o-h2.dat', 'oh2-18o-h2.dat'][int(ortho)]
# gets written to 'moldatapath'



# Q300 value

# either: 
Q300 = 179.639 # ('old') value for the partition function at 300 K
               # for H2-18O (for ortho-para combo) 
               # i.e. not the new values from May 2011 by Brian Druin @ JPL 
               # [only para water], we did not trust them
# from http://www.cv.nrao.edu/php/splat/species_metadata_displayer.php?species_id=620

# or: (for the ortho/para files?)
#~ pQ300 = 44.8740     # 'new' value for para H2-18O 
#~ oQ300 = 150.8666    # 'new' value for ortho H2-18O

#~ Q300 = [pQ300, oQ300][int(ortho)]

########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
#  1
#
# Read in the JPL catalog and convert quantum numbers into strings(?)
#
# Para ('0') or Ortho ('1') water?
#~ oORp = '0'


#~ f= open(jplpath+'c020003.cat','r')
f = open(jplpath + jpldata, 'r')
data = f.read().split('\n')[:-1]
f.close()
class Catalog:
    pass

def v0or1check(inp, v, doit = True):
    if doit:
        return inp == v
    if not doit:
        return True

Catalog.data = data
Catalog.freq = array([float(i[1:13])*1E6 for i in data if v0or1check(i[73:75].strip(), v0orv1, doit=docheck)])      # 12  (1) in Hz from MHz
Catalog.ufreq = array([float(i[13:21])*1E6 for i in data if v0or1check(i[73:75].strip(), v0orv1, doit=docheck)])    # 8   (2) in Hz from MHz
Catalog.linestrength = array([float(i[21:29]) for i in data if v0or1check(i[73:75].strip(), v0orv1, doit=docheck)]) # 8   (3)
Catalog.df = array([int(i[29:31]) for i in data if v0or1check(i[73:75].strip(), v0orv1, doit=docheck)])             # 2   (4)
# this energy level is the energy relative to the ground state energy level
Catalog.elevel = array([float(i[31:41]) for i in data if v0or1check(i[73:75].strip(), v0orv1, doit=docheck)])       # 10  (5)

Catalog.gup = array([int(i[41:44]) for i in data if v0or1check(i[73:75].strip(), v0orv1, doit=docheck)])            # 3   (6)
Catalog.tag = array([int(i[44:51]) for i in data if v0or1check(i[73:75].strip(), v0orv1, doit=docheck)])            # 7   (7)
Catalog.gid = array([int(i[51:55]) for i in data if v0or1check(i[73:75].strip(), v0orv1, doit=docheck)])            # 4   (8)
# a bit nasty read-in of the quantum level numbers
# could perhaps be done in a more elegant way?
Catalog.qup = array(['_'.join((i[55:57].strip(),i[57:59].strip(),i[59:61].strip())) for i in data if v0or1check(i[73:75].strip(), v0orv1, doit=docheck)])  # 8 (9)
# 5 spaces
#~ qlow = array([i[67:75].strip()[:-2] for i in data if i[73:75].strip()==oORp]) # 8 (10)
Catalog.qlow = array(['_'.join((i[67:69].strip(),i[69:71].strip(),i[71:73].strip())) for i in data if v0or1check(i[73:75].strip(), v0orv1, doit=docheck)])

# 
# create our new moldata file and write the first few lines along
# with the energy levels
#
# open reference moldata file for main isotopologue
#~ f = open(moldatapath+'ph2o-h2@daniel.dat','r')
f = open(moldatapath + moldatareference, 'r')
molref = f.read().split('\n')
f.close()
class Moldata:
    pass

molref = molref[:-1]
Moldata.data = molref
Moldata.n_elevels = int(Moldata.data[5])
Moldata.elevels_d = [i.split() for i in Moldata.data[7:7 + Moldata.n_elevels]]

#
# Read in the data for the main isotopologue for the energy levels
# and convert quantum numbers into strings(?)
#
# data1 is the whole molref table
Moldata.wh2o = array([i[2] for i in Moldata.elevels_d],dtype='float')
Moldata.qh2o = array([i[3] for i in Moldata.elevels_d]) # this is the X_Y_Z values

#
# match the energy levels of p-h2-18o with h2o
#
# NOTE : if it doesnt find the energy level, it skips it.
#
# The same energy levels exists in the Catalog.qup array, only the 
# first ('0_0_0') energy level is missing.
indices = [where(i == Catalog.qlow)[0] for i in Moldata.qh2o if len(where(i == Catalog.qlow)[0]) != 0]


neglected = [where(i == Catalog.qlow)[0] for i in Moldata.qh2o if len(where(i == Catalog.qlow)[0]) == 0]
neglected_ID = [i for i in Moldata.qh2o if len(where(i == Catalog.qlow)[0]) == 0]
if len(neglected)>0:
    print ('IMPORTANT: {0} energy levels were NOT matched, this could be significant.'.format(len(neglected)))
    print ('The ignored levels are :\n {0}'.format('\n '.join(neglected_ID)))

class NewMoldata:
    pass

# Get the energy levels 
# (the matches are the same energy levels, so only use first index)
# NewMoldata.en is what is listen in teh ENERGIES(CM^-1) column in 
# the moldata file
NewMoldata.en = array([Catalog.elevel[i[0]] for i in indices])
NewMoldata.nlevels = len(indices)

# 
#
# Read and convert transitions
#
# main header is 7 lines, header for rad. trans is 3 lines (including title for ! TRANS...)
i_rtrans = Moldata.n_elevels+7+3
# main hdr 7 lines, first line "!NUMBER OF RADIATIVE TRANSITIONS"
Moldata.n_rtrans = int(Moldata.data[Moldata.n_elevels+7+1])
NewMoldata.n_rtrans = Moldata.n_rtrans
# get the radiative transitions table from the H2O moldata file
# i_rtrans is 55 for the ph2o-h2@daniel
Moldata.trans_d = array([array(i.split()) for i in Moldata.data[i_rtrans:i_rtrans + Moldata.n_rtrans]])
# now Moldata.trans_d[:,X] is the whole table
# where X (0=TRANS, 1=UP, 2=LOW, 3= EINSTEIN, 4=FREQ, 5=E_up)

Moldata.up = Moldata.trans_d[:,1].astype('int') # up levels (not python slice-friendly)
Moldata.low = Moldata.trans_d[:,2].astype('int') # low

Moldata.qup = Moldata.qh2o[Moldata.up-1] # since we are slicing in Python
Moldata.qlow = Moldata.qh2o[Moldata.low-1]

#~ NewMoldata.up = empty((Moldata.n_rtrans))
#~ NewMoldata.low = empty((Moldata.n_rtrans))
NewMoldata.up = Moldata.up     # we are checking in this order in the
NewMoldata.low = Moldata.low   # for-loop later on.

NewMoldata.Aij = empty((Moldata.n_rtrans))
NewMoldata.freq = empty((Moldata.n_rtrans))
NewMoldata.el = empty((Moldata.n_rtrans))
NewMoldata.eu = empty((Moldata.n_rtrans))
NewMoldata.A = empty((Moldata.n_rtrans))

#~ Q300 = 179.639 # value for the partition function at 300 K

#Now we need to convert the matching entries in the 'catalog' class (jpl data)
#so that it matches what is expected in the moldata radex (lamda) file.
#
#First the frequency is in Hz

# find where the values are in Catalog.qup/qlow (JPL), in the order of
# Moldata.qup/qlow, i.e. Moldata.up/low
# only one transition, or none
test_result = lambda i: where((Catalog.qup == Moldata.qup[i]) * (Catalog.qlow == Moldata.qlow[i]))[0]

not_matched = []

for i in arange(Moldata.n_rtrans):
    # 2012-12-23 - magnusp: changed here so that 
    # I use the size() function, this is correct when checking 
    # numpy/scipy arrays
    if size(test_result(i)): 
        # if a match IS found, get the 'catalog' value and calculate
        # we need to calculate/get : freq, eu, A
        NewMoldata.freq[i] = Catalog.freq[test_result(i)]              # in Hz already
        # convert El from cm-1 to K
        # the ratio C.h/C.k cancels out the cgs or SI units that would change
        # i.e. if in CGS, it gives the same as if in SI. See below both depends on m2
        NewMoldata.el[i] = Catalog.elevel[test_result(i)] * C.c*1e2 * C.h / C.k 
        # UNIT check:
        # cm-1 * (cm/s) * (m2 kg / s)/(J/K) =
        # = cm-1 * cm * (m2 kg / s2) K / (m2 kg / s2) =
        # = K
        # calculate the Eu from El and freq
        NewMoldata.eu[i] = NewMoldata.el[i] + C.h * NewMoldata.freq[i] / C.k
        # NewMoldata.el[i] because we defined it above
        # UNIT check:
        # [el] K + (m2 kg/s) * s-1 / (J/K) =
        # [el] K + K = K
        #
        # linestrength is in units of mm2 MHz at 300 K
        # is this really correct????
        # Aij = It(300K) * nu**2 * (Qrs / gup) * ( exp(-Elow/T) - exp(-Eup/T) )**(-1) * 2.7964 * 10**(-16)
        # 
        # from http://spec.jpl.nasa.gov/ftp/pub/catalog/doc/catintro.pdf
        # Section 3
        It = 10**(Catalog.linestrength[test_result(i)])
        f = exp(-NewMoldata.el[i]/300.) - exp(-NewMoldata.eu[i]/300.)
        constant = 2.7964E-16 # in s-1
        # freq in Hz, so convert to MHz
        NewMoldata.A[i] = It * (NewMoldata.freq[i]*1E-6)**2  * Q300 / Catalog.gup[test_result(i)] * f**-1 * constant
        # Calculate the Einstein A coefficient with Equation 9 from
        # Pickett et al. (1998) "Submillimeter, millimeter and microwave spectral line catalog"
        # where freq should be in MHz
        # now output this as
        # N NewMoldata.up NewMoldata.low NewMoldata.A NewMoldata.freq NewMoldata.eu
    else:# if NO match is found, use 'moldata' value 
        # we need to calculate/get : freq, eu, A
        NewMoldata.A[i] = Moldata.trans_d[i,3]
        NewMoldata.freq[i] = float(Moldata.trans_d[i,4])*1E9  # in Hz from GHz
        # This frequency is WRONG!, need to take the Catalog.freq that 
        # correponds to the Moldata.up - Moldata.low transition
        #~ NewMoldata.freq[i] = float(Moldata.trans_d[i,4])*1E9  # in Hz from GHz
        NewMoldata.eu[i] = float(Moldata.trans_d[i,5])
        not_matched.append(i)

# 
#
# Convert data into correct units
#
"""
The IDL code for this
;Step 7: Convert data into correct units
;
;
;
;File handles
; 1,'o-h2co.dat', moldata reference file
; 2,'c031503.cat' JPL catalog file
; 3,'o-h2c-13-o.dat' new moldata file
;
el=elevel(outputIndex)              # Catalog.elevel[test_result(i)]
freq=freq(outputIndex)              # Catalog.freq[test_result(i)]
dg=gup(outputIndex)                 # Catalog.qup[test_result(i)]
k=1.38066d-23                       # C.k (J/K)
h=6.62608d-34                       # C.h (m2 kg / s)
c=2.997925d10                       # C.c *1e2(?)
eu=dblarr(n_elements(el))
el=el*c*h/k                         # cm-1 * (cm/s) * (m2 kg / s)/(J/K) =
                                    # = cm-1 * cm * (m2 kg / s2) K / (m2 kg / s2) =
                                    # = K
eu=el+h*freq*1.0e6/k                # [el] + (m2 kg/s) * s-1 / (J K-1) =
                                    # [el] + K
It=dblarr(n_elements(el))
It=10^(linestrength(outputIndex))   # 10**(log10(nm2MHz))
f=dblarr(n_elements(el))
f=exp(-el/300.)-exp(-eu/300.)       # e**(-)
mus=dblarr(n_elements(el))
Q300=2956.0542 ; fill in the value for the partition function at 300 K!
mus=2.40251d4*It*Q300/(freq*f)
A=1.16395d-20*freq^3*mus/dg
freq=freq/1.d3

for i=0,ntr-1 do begin
   printf,3,FORMAT='(I3,3X,I2,3X,I2,6X,E9.3,7X,F11.6,4X,F5.1)',i+1,up(i),low(i),A(i),freq(i),eu(i)
end
"""

# 
#
# Copy collisional data from the main isotopologue
#

# FINALLY : Write stuff to file
# remove old
#~ !rm /home/magnusp/work/data/moldata_radex/ph2-18o-h2.dat
# open the file for writing
f_new = open(moldatapath + newmoldatafile, 'w')
#
# write the header
f_new.writelines(['!MOLECULE\n','para-H2-18O (JPL data)\n'])
f_new.writelines(['!MOLECULAR WEIGHT\n','20.\n'])
f_new.writelines(['!NUMBER OF ENERGY LEVELS\n{0}\n'.format(NewMoldata.nlevels)])
# write the energy level table
f_new.write('!LEVEL + ENERGIES(CM^-1) + WEIGHT + J_Kp_Ko\n')
# what is new here is the "en" array
f_new.writelines(['{0: =3} {1: =13.6f} {2: =6.1f} {3:^14}\n'.format(i,j,k,l) for (i,j,k,l) in zip(arange(1, NewMoldata.nlevels+1, 1), NewMoldata.en, Moldata.wh2o, Moldata.qh2o)])
# write header for radiative transitions
f_new.writelines(['!NUMBER OF RADIATIVE TRANSITIONS\n{0}\n'.format(NewMoldata.n_rtrans)])
f_new.writelines(['!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz)+ E_up(K)\n'])
# write the matched radiative transitions as
# i+1,up(i),low(i),A(i),freq(i),eu(i)
format_array = zip(arange(Moldata.n_rtrans), NewMoldata.up, NewMoldata.low, NewMoldata.A, NewMoldata.freq * 1E-9, NewMoldata.eu)
rtrans_output = ['{0:3}    {1:2}    {2:2}    {3:6.3E}    {4: 13.5f}    {5:5.1f}\n'.format(i+1,j,k,l,m,n) for (i,j,k,l,m,n) in format_array]
f_new.writelines(rtrans_output)

#
# Collision rates for main partner
# Perhaps it should account for the mass difference between 
# 16O and 18O? Or does it already? Gah...
# It just prints what ever is there already
retcode = [f_new.writelines(i+'\n') for i in Moldata.data[55+Moldata.n_rtrans:]]
f_new.close()

print("\n{0} transitions not matched and left out.\n".format(len(not_matched)))
