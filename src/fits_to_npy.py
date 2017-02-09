import numpy as np
from astropy.io import fits
import glob

# Load files
files = glob.glob("../data/*.fits")

dtype = [("wavelength","f8"), ("flux","f8")]

# Loop through ../data/*.fits to open, grab spectra, then convert to .npy files
for i,file in enumerate(files):
    hdulist = fits.open(file)

    if i%500 == 0: print(i)

    # Load data into useful structure
    wavelength = 10**hdulist[1].data['loglam']
    flux = hdulist[1].data['flux']

    spectrum = np.zeros(len(wavelength), dtype=dtype)
    spectrum['wavelength'] = wavelength
    spectrum['flux'] = flux

    fileout = file[:-4]+"npy"

    np.save(fileout, spectrum)
