
import os as _os


class Read:
    def __init__(self, path_to_file):
        from scipy import arange
        f = open(path_to_file)
        self.molfile = path_to_file
        molref = f.read().split('\n')
        f.close()
        self.molref = list(molref)
        # get the preamble, all the molecular properties
        #~ i = 0
        #~ while molref[i] != '!LEVEL + ENERGIES(cm^-1) + WEIGHT + J, Ka, Kc':
        for i in arange(len(molref)):
            if '!LEVEL + ENERGIES'.lower() in str(molref[i]).lower():
                ilvl = i + 1
            elif '!TRANS + UP + LOW + EINSTEINA'.lower() in str(molref[i]).lower():
                rtlvl = i + 1
            
            # get the variables
            txt = molref[i]
            if txt.strip('!') == 'MOLECULE':
                setattr(self, 'molecule', molref[i+1])
            elif txt.strip('!') == 'MOLECULAR WEIGHT':
                setattr(self, 'molweight', float(molref[i+1]))
            elif txt.strip('!') == 'NUMBER OF ENERGY LEVELS':
                setattr(self, 'n_elevels', int(molref[i+1]))
            elif  txt.strip('!') == 'NUMBER OF RADIATIVE TRANSITIONS':
                setattr(self, 'n_radtrans', int(molref[i+1]))
            #~ elif  txt.strip('!') == '':
                #~ setattr(self, '', )
            #~ elif  txt.strip('!') == '':
            
            # below so it wont loop forever...
            #~ if i > 40:
                #~ print 'reached maximum iterations (40)'
                #~ break
            #~ i += 1
        #~ ilvl = i + 1
        
        # read in the energy levels
        elevels_d = [i.split() for i in molref[ilvl:ilvl + self.n_elevels]]
        self.elev = list([dict([['level', int(i[0])], 
                        ['energies', float(i[1])], 
                        ['weight', float(i[2])], 
                        ['j', str(i[3])]]) for i in elevels_d])
        
        # read in the radiative transitions
        radtrans_d = [i.split() for i in molref[rtlvl:rtlvl+self.n_radtrans]]
        #~ print rtlvl, rtlvl+self.n_radtrans
        self.radtrans = list([dict([['trans', int(i[0])],
                        ['ulev', self.get_lvl(i[1])],
                        ['up', int(i[1])],
                        ['llev', self.get_lvl(i[2])],
                        ['down', int(i[2])],
                        ['aul', float(i[3])],
                        ['freq', float(i[4])*1e9],
                        ['eu', float(i[5])]]) for i in radtrans_d])
        
    def get_lvl(self, level, tex = False):
        if not tex:
            return self.elev[int(level)-1]['j']
        if tex:
            lvl = self.elev[int(level)-1]['j'].split('_')
            return ('${0[0]}_{{{0[1]},{0[2]}}}$'.format(lvl))

    def get_partition(self, 
            part_file_directory = '/home/magnusp/work/data/partition/', 
            part_file = ''):
        #
        from scipy.optimize import leastsq as _ls
        from scipy import loadtxt as _loadtxt
        
        molecule = self.molecule
        
        # if no filename has been input directly
        if not part_file:
            part_file = [i for i in _os.listdir(part_file_directory) if i.strip('.dat') in molecule.lower()][0]
        
        # now load the file
        try:
            t,q = _loadtxt(_os.path.join(part_file_directory, part_file))
        except :
            print('partition file could not be loaded')
            import sys as _sys
            _sys.exit('partition file could not be loaded')       
        
        # define the fitting function
        # in this case [a * T**b]
        fitfunc = lambda p,t: q - (p[0]*t**p[1])
        # initial guess...
        p0 = [0.05,1.5] 
        
        # now least square fitting
        p_fit = _ls(fitfunc, p0, args=(t))[0]
        
        # create the function to get the partition function value for 
        # given temperature
        self.qrot = lambda x: p_fit[0] * x**p_fit[1]
        
    def get_transition(self, transition):
        
        pass

class Transition(object):
    def __init__(self, levels=[12, 10]):
        """
        - input population levels
        
        - calculate tex for transition
        
        - get tau
        
        -> put in moldata calculations?
        
        """
        # store everything necessary for calculations
        
        def _calc_level_pop(self, tex):
            self.level_pop = []
        def _calc_tex(self):
            
            #~ _tex_transition(n1, n2, g1, g2, nu):
            
            self.tex = []
        
        def _calc_tau(self):
            self.tau = []
        
        def _calc_radiation_field(self):
            pass
