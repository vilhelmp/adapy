#
# Magnus Persson < vilhelm.nu >
#
# Script to create a moldata file for pH2(18)O-H2, combining
# the JPL catalog file for H2(18)O and the existing H2O moldata file.
#
# Converted from a script by Suzanne Bisschop for CH3OCH3
#

from scipy import *
import scipy.constants as C
jplpath = '/home/magnusp/work/data/molecules/H2-(18)-O/'
moldatapath = '/home/magnusp/work/data/moldata/'

oORp = '0'

jpldata = 'c020003_11122012.cat' #H2-18O file from JPL molecular catalog 
# http://spec.jpl.nasa.gov/ftp/pub/catalog/catdir.html

moldata_reference = ['ph2o-h2@daniel.dat' , 'oh2o-h2@daniel.dat'][int(oORp)]
# http://home.strw.leidenuniv.nl/~moldata/H2O.html

newmoldatafile = ['ph2-18o-h2.dat', 'oh2-18o-h2.dat'][int(oORp)]
# gets written to 'moldatapath'

Q300 = 179.639 # ('old') value for the partition function at 300 K
               # for H2-18O (for ortho-para combo) 
               # i.e. not the new values from May 2011 by Brian Druin @ JPL 
               # [only para water], we did not trust them
# from http://www.cv.nrao.edu/php/splat/species_metadata_displayer.php?species_id=620


#  1
#
# Read in the JPL catalog and convert quantum numbers into strings(?)
#
# Para ('0') or Ortho ('1') water?
#~ oORp = '0'


#~ f= open(jplpath+'c020003.cat','r')
f= open(jplpath + jpldata ,'r')
data = f.read().split('\n')[:-1]
f.close()
class catalog:
    pass

catalog.data = data
catalog.freq = array([float(i[1:13])*1E6 for i in data if i[73:75].strip()==oORp])      # 12  (1) in Hz from MHz
catalog.ufreq = array([float(i[13:21])*1E6 for i in data if i[73:75].strip()==oORp])    # 8   (2) in Hz from MHz
catalog.linestrength = array([float(i[21:29]) for i in data if i[73:75].strip()==oORp]) # 8   (3)
catalog.df = array([int(i[29:31]) for i in data if i[73:75].strip()==oORp])             # 2   (4)
catalog.elevel = array([float(i[31:41]) for i in data if i[73:75].strip()==oORp])       # 10  (5)
catalog.gup = array([int(i[41:44]) for i in data if i[73:75].strip()==oORp])            # 3   (6)
catalog.tag = array([int(i[44:51]) for i in data if i[73:75].strip()==oORp])            # 7   (7)
catalog.gid = array([int(i[51:55]) for i in data if i[73:75].strip()==oORp])            # 4   (8)
# a bit nasty read-in of the quantum level numbers
# could perhaps be done in a more elegant way?
catalog.qup = array(['_'.join((i[55:57].strip(),i[57:59].strip(),i[59:61].strip())) for i in data if i[73:75].strip()==oORp])  # 8 (9)
# 5 spaces
#~ qlow = array([i[67:75].strip()[:-2] for i in data if i[73:75].strip()==oORp]) # 8 (10)
catalog.qlow = array(['_'.join((i[67:69].strip(),i[69:71].strip(),i[71:73].strip())) for i in data if i[73:75].strip()==oORp])

# 
# create our new moldata file and write the first few lines along
# with the energy levels
#
# open reference moldata file for main isotopologue
#~ f = open(moldatapath+'ph2o-h2@daniel.dat','r')
f = open(moldatapath + moldata_reference,'r')
molref = f.read().split('\n')
f.close()
class moldata:
    pass

molref = molref[:-1]
moldata.data = molref
moldata.n_elevels = int(moldata.data[5])
moldata.elevels_d = [i.split() for i in moldata.data[7:7+moldata.n_elevels]]

#
# Read in the data for the main isotopologue for the energy levels
# and convert quantum numbers into strings(?)
#
# data1 is the whole molref table
moldata.wh2o = array([i[2] for i in moldata.elevels_d],dtype='float')
moldata.qh2o = array([i[3] for i in moldata.elevels_d]) # this is the X_Y_Z values

#
# match the energy levels of p-h2-18o with h2o
#
# NOTE : if it doesnt find the energy level, it skips it.
#

indices = [where(i == catalog.qlow)[0] for i in moldata.qh2o if len(where(i == catalog.qlow)[0]) != 0]

neglected = [where(i == catalog.qlow)[0] for i in moldata.qh2o if len(where(i == catalog.qlow)[0]) == 0]
neglected_ID = [i for i in moldata.qh2o if len(where(i == catalog.qlow)[0]) == 0]
if len(neglected)>0:
    print ('IMPORTANT: {0} energy levels were NOT matched, this could be significant.'.format(len(neglected)))
    print ('The ignored levels are :\n {0}'.format('\n '.join(neglected_ID)))

class new_moldata:
    pass

# get the energy levels
new_moldata.en = array([catalog.elevel[i[0]] for i in indices])
new_moldata.nlevels = len(indices)

# 
#
# Read and convert transitions
#
# main header is 7 lines, header for rad. trans is 3 lines (including title for ! TRANS...)
i_rtrans = moldata.n_elevels+7+3
# main hdr 7 lines, first line "!NUMBER OF RADIATIVE TRANSITIONS"
moldata.n_rtrans = int(moldata.data[moldata.n_elevels+7+1])
new_moldata.n_rtrans = moldata.n_rtrans
# get the radiative transitions table from the H2O moldata file
# i_rtrans is 55 for the ph2o-h2@daniel
moldata.trans_d = array([array(i.split()) for i in moldata.data[i_rtrans:i_rtrans + moldata.n_rtrans]])
# now moldata.trans_d[:,X] is the whole table
# where X (0=TRANS, 1=UP, 2=LOW, 3= EINSTEIN, 4=FREQ, 5=E_up)

moldata.up = moldata.trans_d[:,1].astype('int') # up levels (not python slice-friendly)
moldata.low = moldata.trans_d[:,2].astype('int') # low

moldata.qup = moldata.qh2o[moldata.up-1] # since we are slicing in Python
moldata.qlow = moldata.qh2o[moldata.low-1]

# find where the values are in catalog.qup/qlow (JPL), in the order of
# moldata.qup/qlow, i.e. moldata.up/low
test_result = lambda i: where((catalog.qup == moldata.qup[i]) * (catalog.qlow == moldata.qlow[i]))[0]

#~ new_moldata.up = empty((moldata.n_rtrans))
#~ new_moldata.low = empty((moldata.n_rtrans))
new_moldata.up = moldata.up     # we are checking in this order in the
new_moldata.low = moldata.low   # for-loop later on.

new_moldata.Aij = empty((moldata.n_rtrans))
new_moldata.freq = empty((moldata.n_rtrans))
new_moldata.el = empty((moldata.n_rtrans))
new_moldata.eu = empty((moldata.n_rtrans))
new_moldata.A = empty((moldata.n_rtrans))

#~ Q300 = 179.639 # value for the partition function at 300 K

#Now we need to convert the matching entries in the 'catalog' class (jpl data)
#so that it matches what is expected in the moldata radex (lamda) file.
#
#First the frequency is in Hz

not_matched = []

for i in arange(moldata.n_rtrans):
    if not test_result(i): # if NO match is found, use 'moldata' value
        # we need to calculate/get : freq, eu, A
        new_moldata.A[i] = moldata.trans_d[i,3]
        new_moldata.freq[i] = float(moldata.trans_d[i,4])*1E9  # in Hz from GHz
        new_moldata.eu[i] = float(moldata.trans_d[i,5])
        not_matched.append(i)
    else: # if a match IS found, get the 'catalog' value and calculate
        # we need to calculate/get : freq, eu, A
        new_moldata.freq[i] = catalog.freq[test_result(i)]              # in Hz already
        # convert El from cm-1 to K
        new_moldata.el[i] = catalog.elevel[test_result(i)] * C.c*1e2 * C.h / C.k
        # cm-1 * (cm/s) * (m2 kg / s)/(J/K) =
        # = cm-1 * cm * (m2 kg / s2) K / (m2 kg / s2) =
        # = K
        # calculate the Eu from El and freq
        new_moldata.eu[i] = new_moldata.el[i] + C.h * new_moldata.freq[i] / C.k
        # new_moldata.el[i] because we defined it above
        # [el] K + (m2 kg/s) * s-1 / (J/K) =
        # [el] K + K = K
        It = 10**(catalog.linestrength[test_result(i)])
        f = exp(-new_moldata.el[i]/300.) - exp(-new_moldata.eu[i]/300.)
        constant = 2.7964E-16 # in s-1
        new_moldata.A[i] = constant * (new_moldata.freq[i]*1E-6)**2 * It * Q300 /(catalog.gup[i] * f)
        # Calculate the Einstein A coefficient with Equation 9 from
        # Pickett et al. (1998) "Submillimeter, millimeter and microwave spectral line catalog"
        # where freq should be in MHz
        # now output this as
        # N new_moldata.up new_moldata.low new_moldata.A new_moldata.freq new_moldata.eu

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
el=elevel(outputIndex)              # catalog.elevel[test_result(i)]
freq=freq(outputIndex)              # catalog.freq[test_result(i)]
dg=gup(outputIndex)                 # catalog.qup[test_result(i)]
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
f_new.writelines(['!NUMBER OF ENERGY LEVELS\n{0}\n'.format(new_moldata.nlevels)])
# write the energy level table
f_new.write('!LEVEL + ENERGIES(CM^-1) + WEIGHT + J_Kp_Ko\n')
# what is new here is the "en" array
f_new.writelines(['{0: =3} {1: =13.6f} {2: =6.1f} {3:^14}\n'.format(i,j,k,l) for (i,j,k,l) in zip(arange(1,46,1), new_moldata.en, moldata.wh2o, moldata.qh2o)])
# write header for radiative transitions
f_new.writelines(['!NUMBER OF RADIATIVE TRANSITIONS\n{0}\n'.format(new_moldata.n_rtrans)])
f_new.writelines(['!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz)+ E_up(K)\n'])
# write the matched radiative transitions as
# i+1,up(i),low(i),A(i),freq(i),eu(i)
format_array = zip(arange(moldata.n_rtrans), new_moldata.up, new_moldata.low, new_moldata.A, new_moldata.freq * 1E-9, new_moldata.eu)
rtrans_output = ['{0:3}    {1:2}    {2:2}    {3:6.3E}    {4: 13.5f}    {5:5.1f}\n'.format(i+1,j,k,l,m,n) for (i,j,k,l,m,n) in format_array]
f_new.writelines(rtrans_output)

# Collision rates for main partner
# Perhaps it should account for the mass difference between 
# 16O and 18O? Or does it already? Gah...
retcode = [f_new.writelines(i+'\n') for i in moldata.data[55+moldata.n_rtrans:]]
f_new.close()
