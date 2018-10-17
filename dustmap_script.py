import sys
import astropy.units as units
from astropy.coordinates import SkyCoord

from dustmaps.config import config
from dustmaps.bayestar import BayestarQuery
from dustmaps.chen2014 import Chen2014Query
from dustmaps.iphas import IPHASQuery
from dustmaps.marshall import MarshallQuery
from dustmaps.sfd import SFDQuery
from dustmaps.lenz2017 import Lenz2017Query
from dustmaps.bh import BHQuery

config['data_dir'] = '/home/arash/astro_sw/dustmap/'

if len(sys.argv) == 3:
    source_ra = sys.argv[1]
    source_dec = sys.argv[2]
    distance = 8.0
    print('Assuming d = 8 kpc for the 3D maps')
elif len(sys.argv) == 4:
    source_ra = sys.argv[1]
    source_dec = sys.argv[2]
    distance = eval(sys.argv[3])
else:
	print('wrong number of parameters')
	print('Example: $ python dustmap_script.py 15h28m17.53s -58d35m13.9s 2.3')
	sys.exit(1)

#source_ra = '15h28m17.53s'
#source_dec = '-58d35m13.9s'
#distance = 2.3
#print(source_ra, source_dec, distance)

print('Estimating E(B-V) using dust maps:')

coords = SkyCoord(source_ra, source_dec, distance=distance*units.kpc, frame='icrs')

bayestar = BayestarQuery(max_samples=5, version='bayestar2017')
print('## Bayestar (3D):', bayestar(coords, mode='median'))

chen = Chen2014Query()
print('## Chen 2014 (3D):', chen(coords))

iphas = IPHASQuery()
print('## IPHAS (3D):', iphas(coords, mode='median'))

marshal = MarshallQuery()
print('## Marshall 2006 (3D):', marshal(coords))

sfdebv = SFDQuery()
print('## SFD 1998/2011 (2D):', sfdebv(coords))

lenz = Lenz2017Query()
print('## Lenz 2017 (2D):', lenz(coords))


bh = BHQuery()
print('## Burstein 1982 (2D):', bh(coords))

print('\n WARNING: REMEMBER TO USE APPROPRIATE BAND CONVERSION FACTORS.\nSEE http://iopscience.iop.org/article/10.1088/0004-637X/737/2/103/meta#apj398709t6')
