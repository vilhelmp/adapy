
class Read:
    def __init__(self, path_to_file):
        f = open(path_to_file)
        self.molfile = path_to_file
        molref = f.read().split('\n')
        f.close()
        self.n_elevels = int(molref[5])
        elevels_d = [i.split() for i in molref[7:7 + self.n_elevels]]
        self.elev = list([dict([['level', int(i[0])], 
                        ['energies', float(i[1])], 
                        ['weight', float(i[2])], 
                        ['j', str(i[3])]]) for i in elevels_d])
    def get_lvl(self, level, tex = False):
        if not tex:
            return self.elev[int(level)-1]['j']
        if tex:
            lvl = self.elev[int(level)-1]['j'].split('_')
            return ('${0[0]}_{{{0[1]},{0[2]}}}$'.format(lvl))
