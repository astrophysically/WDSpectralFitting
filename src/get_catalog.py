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

catalogs[0].write("../data/kepler_WDs_DR12.csv")
