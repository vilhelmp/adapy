from scipy import *
import cPickle as pick

#create empty database
db = {}
# add objects
db['iras4b'] = {'restfreq':False, 'line':'HCN', 'observer':'jkj' }
db['iras4a'] = {'restfreq':220, 'line':'13CO', 'observer':'test' }
db['iras2a']= {'restfreq':230.538, 'line':'12CO', 'observer':'test' }


dbfile = open("db.cpick", 'w')
pick.dump(db, dbfile)
dbfile.close()

db = {}
dbfile = []

dbfile= open("db.cpick", 'r')
db = pick.load(dbfile)
dbfile.close()
print db.keys()


def genAscii(database):
	with f as open('ascii.file','w'):
		for i in sorted(database.keys()):
			out = 
	
	
	
	
