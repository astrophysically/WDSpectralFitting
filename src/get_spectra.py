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
Vizier.ROW_LIMIT = 2
catalogs = Vizier.get_catalogs(catalog_list.keys())


# cat = 'J/MNRAS/455/3413/table6'
PMF = catalogs[0]['PMF']
#for i in range(len(PMF)):
#    PMF[i] = PMF[i].decode("utf-8")


print("Downloading " + str(len(PMF)) + " spectra now...") 

for i in np.arange(len(PMF)):

    filename = '../data/spectra-'+str(PMF[i].decode("utf-8"))+'.fits'
    if os.path.exists(filename): continue

    # print("PMF:", PMF[i].decode("utf-8"))
    plate = int(str(PMF[i].decode("utf-8"))[0:4])
    mjd = int(str(PMF[i].decode("utf-8"))[5:10])
    fiber = int(str(PMF[i].decode("utf-8"))[11:15])

    try:
        spec = SDSS.get_spectra(plate=plate, mjd=mjd, fiberID=fiber)
        spec[0].writeto(filename)
    except:
        print("Could not download spectra:", plate, mjd, fiber) 
        pass


print("...finished downloading spectra.") 
