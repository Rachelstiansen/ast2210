import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# Opening the file in python:

filename = "ADP.2017-03-27T12_08_50.541.fits"
hdu = fits.open(filename)
#hdu.info()

# From hdu.info we get the information below:
"""
Filename: ADP.2017-03-27T12_08_50.541.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU    1345   ()      
  1  DATA          1 ImageHDU        44   (320, 317, 3682)   float32   
  2  STAT          1 ImageHDU        44   (320, 317, 3682)   float32   
"""

data = hdu[1].data
hdr = hdu[1].header
# print(data.shape) 
# has shape (3682, 317, 320)

# Now we produce a map of the data cube by averaging the
# flux values within a range of wavelengths

flux_mean = np.nanmean(data, 0)
# print(flux_mean.shape)
# has shape (317, 320)

"""
plt.figure(figsize=(7, 7)) # Setting the size of the image
plt.title(r"Average flux density in the range $\lambda \in [4750, 9351]$ Å")
im = plt.imshow(np.flip(np.nanmean(data,0),0),cmap="gray", vmin=0,vmax=2137)
plt.colorbar(im, fraction=0.046, pad=0.04, label="Flux density [$10^{-20}$ergs$^{-1}$cm$^{-2}\AA{}^{-1}$]")

plt.savefig("integrated_pix.png")
plt.show()
"""
# collapsing hte data within a specific wavelength range:

lower_indx = 6590
upper_indx = 6610

extracted_data = data[lower_indx:upper_indx] # cuts the data in the wavelength dimension

lambda0 = hdr["CRVAL3"] # extract value of the wavelengths lower boundary (smallest wavelength) [A]
dlambda = hdr["CD3_3"] # wavelength resolution, or delta lambda if you will [A]
len_wave = hdr["NAXIS3"] # Number of datapoints in the wavelength dimension
wavelengths = np.linspace(lambda0, lambda0 + (len_wave-1)*dlambda, len_wave)

lower_boundary=6590;upper_boundary=6610 # upper and lower boundary wavelengths in Angstrom.
lower_indx=np.array(np.where(wavelengths>=lower_boundary))[0,0] # locating the index
upper_indx=np.array(np.where(wavelengths<=upper_boundary))[-1,-1]

"""
hdr["CRVAL3"] = lower_wavelength # Updates the header value for the lower wavelength
hdu = fits.PrimaryHDU(extracted_data,header=hdr) # creates an extension
hdul = fits.HDUList([hdu]) # includes the extension into the extension list
hdul.info() # prints info about the extensions (optional)
hdul.writeto("newfile.fits")
"""

# Determining the average flux
r = 5
center_x = 170
center_y = 150
indx = [center_y-r, center_y+r,center_x-r,center_x+r] # [y0, y1 , x0, x1] # aperture indices

collapsed = np.nansum(data,0) # summing all non-nan values in the spectral dimension
aperture_data = collapsed[indx[0]:indx[1], indx[2]:indx[3]]

mean_flux = np.mean(aperture_data) # taking the mean over both spatial dimensions
sum_flux = np.sum(aperture_data) # summing up all the flux contained in the aperture.
std_dev = np.std(aperture_data)

from astropy.wcs import WCS
plt.rcParams["font.size"] = 14

#extracting the spatial WCS from the header
wcs = WCS(hdr)[0,:,:] # indexing[spectral, vertical, horizontal]

Nx = 166; Ny = 156
Bx = 190; By = 131
Ax = 156; Ay = 211

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(1,1,1, projection=wcs)
ax.set_title(r"Average flux density in the range $\lambda \in [4750, 9351]$ Å")

im = ax.imshow(np.nanmean(data,0),cmap="viridis", vmin=0,vmax=2137) # no longer needs np.flip because of WCS

plt.plot(Nx, Ny, "x", color="red")
plt.text(Nx - 12, Ny - 12, "N", color="red")
plt.plot(Bx, By, "x", color="red")
plt.text(Bx - 12, By - 12, "B", color="red")
plt.plot(Ax, Ay, "x", color="red")
plt.text(Ax - 12, Ay - 12, "A", color="red")

plt.xlabel("RA")
plt.ylabel("Dec")
plt.colorbar(im,fraction=0.046, pad=0.04, label="Flux density [$10^{-20}$ergs$^{-1}$cm$^{-2}\AA{}^{-1}$]")
plt.legend()
plt.savefig("integrated_arcsec.png")
plt.show()
