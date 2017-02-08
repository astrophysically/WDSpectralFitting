import os.path
import numpy as np

import astropy
from astropy.table import Table, Column, MaskedColumn

from astroquery.vizier import Vizier
from astroquery.sdss import SDSS

# Get catalog from vizier search
catalog_list = Vizier.find_catalogs('New white dwarf SDSS DR12')

# Make sure we have the correct catalog
print({k:v.description for k,v in catalog_list.items()})

# Get the entire catalog, not just the first 50 rows
Vizier.ROW_LIMIT = -1
catalogs = Vizier.get_catalogs(catalog_list.keys())


# cat = 'J/MNRAS/455/3413/table6'
PMF = catalogs[0]['PMF']


print "Downloading " + str(len(PMF)) + " spectra now..."

for i in np.arange(len(PMF)):

    filename = '../data/spectra-'+str(PMF[i])+'.fits'
    if os.path.exists(filename): continue

    plate = int(str(PMF[i])[0:4])
    mjd = int(str(PMF[i])[5:10])
    fiber = int(str(PMF[i])[11:15])


    spec = SDSS.get_spectra(plate=plate, mjd=mjd, fiberID=fiber)

    spec[0].writeto(filename)

print "...finished downloading spectra."
