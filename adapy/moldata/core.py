
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
            if molref[i] == '!LEVEL + ENERGIES(cm^-1) + WEIGHT + J, Ka, Kc':
                ilvl = i + 1
            elif molref[i] == '!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_u(K)':
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
                        ['Aul', float(i[3])],
                        ['freq', float(i[4])*1e9],
                        ['eu', float(i[5])]]) for i in radtrans_d])
        
        
    def get_lvl(self, level, tex = False):
        if not tex:
            return self.elev[int(level)-1]['j']
        if tex:
            lvl = self.elev[int(level)-1]['j'].split('_')
            return ('${0[0]}_{{{0[1]},{0[2]}}}$'.format(lvl))
