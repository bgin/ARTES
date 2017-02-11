import numpy as np
import math, os
from astropy.io import fits

# ------------------------------------------------------------
# User input

atmosphere = 'isotropic'
fitsOutput = 'isotropic.fits'

absorption = 0.0 # Absorption opacity [cm2 g-1]
scattering = 1.0 # Scattering opacity [cm2 g-1]

# Wavelengths from FITS file
fitsWavelength = False
wavelengthFile = ''

# Wavelength range [micron]
manualWavelength = True
wavelengthMin = 1.2
wavelengthMax = 1.2
step = 1.0

# ------------------------------------------------------------

scriptDir =  os.path.dirname(os.path.abspath(__file__))
fitsOutput = scriptDir[:-6]+'input/'+atmosphere+'/opacity/'+fitsOutput
wavelengthFile = scriptDir[:-6]+'input/'+atmosphere+'/opacity/'+wavelengthFile

if fitsWavelength:
    
    hdulist = fits.open(wavelengthFile)
    hdu = hdulist[0].data
    wavelength = hdu[0]
    hdulist.close()

elif manualWavelength:

    wavelength = []
    for i in range(int((wavelengthMax-wavelengthMin)/step)+1):
        wavelength.append(wavelengthMin+float(i)*step)

opacity = np.zeros((4,len(wavelength)))
for i in range(len(wavelength)):

    opacity[0,i] = wavelength[i]
    opacity[1,i] = absorption + scattering
    opacity[2,i] = absorption
    opacity[3,i] = scattering

scatter = np.zeros((180,16,len(wavelength)))

for i in range(len(wavelength)):
    for j in range(180):

        scatter[j,0,i] = 1./(4.*math.pi)

hdulist = fits.HDUList()
hdulist.append(fits.ImageHDU(np.array(opacity), name='opacity'))
hdulist.append(fits.ImageHDU(np.array(scatter), name='scattermatrix'))
hdu = hdulist[0]
hdu.header['COMMENT'] = '1. Wavelength [micron]'
hdu.header['COMMENT'] = '2. Extinction [cm2 g-1]'
hdu.header['COMMENT'] = '3. Absorption [cm2 g-1]'
hdu.header['COMMENT'] = '4. Scattering [cm2 g-1]'
hdulist.writeto(fitsOutput, clobber=True)
hdulist.close()