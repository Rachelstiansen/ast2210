import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.colors as colors
from matplotlib.patches import Rectangle
from astropy import cosmology

plt.rcParams["font.size"] = 11.5


# From hdu.info we get the information below:
"""
Filename: ADP.2017-03-27T12_08_50.541.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU    1345   ()      
  1  DATA          1 ImageHDU        44   (320, 317, 3682)   float32   
  2  STAT          1 ImageHDU        44   (320, 317, 3682)   float32   
"""

Nx = 166; Ny = 156 # center of galaxy
Bx = 191; By = 131
Ax = 135; Ay = 233

# Opening the file in python:
filename = "ADP.2017-03-27T12_08_50.541.fits"
hdu = fits.open(filename)
#hdu.info()
data = hdu[1].data
hdr = hdu[1].header
# extracting the spatial WCS from the header
wcs = WCS(hdr)[0,:,:] # indexing[spectral, vertical, horizontal]

"""
Finding indexes with max flux:
"""
# Define the desired wavelength
desired_wavelength = 6599  # Replace with your desired wavelength in the appropriate units

lambda0 = hdr["CRVAL3"] # extract value of the wavelengths lower boundary (smallest wavelength) [A]
dlambda = hdr["CD3_3"] # wavelength resolution, or delta lambda if you will [A]
len_wave = hdr["NAXIS3"] # Number of datapoints in the wavelength dimension
wavelengths = np.linspace(lambda0, lambda0 + (len_wave-1)*dlambda, len_wave)



# Find the index corresponding to the desired wavelength
index = np.abs(wavelengths - desired_wavelength).argmin()

# Get the 2D spatial slice at the desired wavelength
spatial_slice = data[index, :, :]

# Find the pixel with the highest flux at the desired wavelength
max_flux = np.nanmax(spatial_slice)

max_flux_pixel = np.unravel_index(np.nanargmax(spatial_slice), spatial_slice.shape)

# max_flux_pixel now contains the (x, y) coordinates of the pixel with the highest flux at the desired wavelength
y, x = max_flux_pixel
print(f"Pixel with highest flux at {desired_wavelength} nm: ({x}, {y})")


def integrated_map():
  # Opening the file in python:
  filename = "ADP.2017-03-27T12_08_50.541.fits"
  hdu = fits.open(filename)
  #hdu.info()
  data = hdu[1].data
  hdr = hdu[1].header
  # extracting the spatial WCS from the header
  wcs = WCS(hdr)[0,:,:] # indexing[spectral, vertical, horizontal]
  fig = plt.figure(figsize=(8,8))
  ax = plt.subplot(1,1,1, projection=wcs)
  ax.set_title(r"Total integrated map of NGC 1365 in the range $\lambda \in [4750, 9351]$ Å")

  im = ax.imshow(np.nanmean(data, 0),cmap="jet", norm=colors.LogNorm()) # no longer needs np.flip because of WCS

  plt.plot(Nx, Ny, "x", color="black")
  plt.text(Nx - 12, Ny - 12, "N", color="black")
  plt.plot(Bx, By, "x", color="black")
  plt.text(Bx - 12, By - 12, "B", color="black")
  plt.plot(Ax, Ay, "x", color="black")
  plt.text(Ax - 12, Ay - 12, "A", color="black")

  plt.xlabel("RA (J2000)")
  plt.ylabel("Dec (J2000)")
  plt.colorbar(im, fraction=0.046, pad=0.04, label="Logarithmic flux density [log($10^{-20}$ergs$^{-1}$cm$^{-2}\AA{}^{-1}$)]")
  plt.legend()
  plt.tight_layout()
  plt.savefig("integrated_RA_Dec.png")
  plt.show()

integrated_map()

def spectral_plot(file, title, filename):
    fig = plt.figure(figsize=(8,3))
    data = np.loadtxt(file)
    wav = data[:, 0]
    flux = data[:, 1]
    plt.plot(wav, flux)

    #plt.axvline(5034.2, label="$[OIII]^*N_2$", ymin=0.05, ymax=0.98, alpha=0.85, color="g", linestyle="--")
    #plt.axvline(6598.7, label="$H_\\alpha$", ymin=0.05, ymax=0.98, alpha=0.85, color="blue", linestyle="--")
    #plt.axvline(6584, label="$[NII]^*$", ymin=0.05, ymax=0.98, alpha=0.85, color="c", linestyle="--")
    #plt.axvline(6619.5, ymin=0.05, ymax=0.98, alpha=0.85, color="c", linestyle="--")
    
    #plt.axvline(5907.8, label="$HeI^*$", ymin=0.05, ymax=0.98, alpha=0.85, color="brown", linestyle="--")
    #plt.axvline(6009.7, label="$HeII^*$", ymin=0.05, ymax=0.98, alpha=0.85, color="red", linestyle="--")
    #plt.axvline(6334.8, label="$[OII]$", ymin=0.05, ymax=0.98, alpha=0.85, color="g", linestyle="--")
    #plt.axvline(4986.1, label="$[OIII]^*N_1$", ymin=0.05, ymax=0.98, alpha=0.85, color="blue", linestyle="--")
    #plt.axvline(4887.9, label="$H_\\beta$", ymin=0.05, ymax=0.98, alpha=0.85, color="c", linestyle="--")
    #plt.axvline(6753.8, label="$SII^*$", ymin=0.05, ymax=0.98, alpha=0.85, color="purple", linestyle="--")
    #plt.axvline(6768.2, ymin=0.05, ymax=0.98, alpha=0.85, color="purple", linestyle="--")

    #plt.text(6598.7-50, 1.4e6, "$H\\alpha$")
    #plt.text(6584, 0.9e6, "$[NII]*$")
    #plt.text(6584, 0.9e6, "$[NII]*$")
    #plt.text(6740, 0.3e6, "$SII*$")
    #plt.text(5020, 0.6e6, "$[OIII]*N_2$")
    #plt.text(4970, 0.3e6, "$[OIII]*N_1$")
    #plt.text(4880, 0.3e6, "$H\\beta$")

    # Define the wavelength range
    min_wavelength = 5800 
    max_wavelength = 6000

    # Create a boolean mask to filter wavelengths within the specified range
    mask = (wav >= min_wavelength) & (wav <= max_wavelength)

    # Use the masked flux values to find the index of the maximum flux within the range
    max_flux_index = np.argmin(flux[mask])

    # Get the wavelength corresponding to the maximum flux within the range
    max_flux_wavelength = wav[mask][max_flux_index]

    print("peak wav = ", max_flux_wavelength)

    plt.title(title)
    plt.xlim(4800, 5100)
    #plt.xlim(6500, 6800)
    plt.xlabel("Wavelength [Å]")
    plt.ylabel("Flux density [$10^{-20} erg/s/cm^2/Å$]")
    plt.grid()
    #plt.legend(loc="upper right", ncol=3)
    plt.tight_layout()
    #plt.savefig(filename)
    #plt.show()

#spectral_plot("center_166_156_aperture5.txt", "Spectra extracted from galaxy center", "lines_N.png")
#spectral_plot("center_166_156_aperture5.txt", "Spectra extracted from galaxy center, $\lambda \in [4800, 5100]$", "lines_N_left.png")
#spectral_plot("center_166_156_aperture5.txt", "Spectra extracted from galaxy center, $\lambda \in [6500, 6800]$", "lines_N_right.png")

#spectral_plot("B_191_131_aperture5.txt", "Spectra extraced from point B", "line_B.png")
#spectral_plot("B_191_131_aperture5.txt", "Spectra extraced from point B, $\lambda \in [4800, 5100]$", "line_B_left.png")
#spectral_plot("B_191_131_aperture5.txt", "Spectra extraced from point B, $\lambda \in [6500, 6800]$", "line_B_right.png")

#spectral_plot("A_135_233_aperture5.txt", "Spectra extracted from point A", "line_A.png")
#spectral_plot("A_135_233_aperture5.txt", "Spectra extracted from point A, $\lambda \in [4800, 5100]$", "line_A_left.png")
#spectral_plot("A_135_233_aperture5.txt", "Spectra extracted from point A, $\lambda \in [6500, 6800]$", "line_A_right.png")



def line_map(filename, savename, title):
  hdu = fits.open(filename)
  hdu.info()
  data = hdu[0].data
  #hdr = hdu[0].header
  #wcs = WCS(hdr)[:,:] # indexing[vertical, horizontal]

  fig = plt.figure(figsize=(8,8))
  ax = plt.subplot(1,1,1, projection=wcs)
  ax.set_title(title)
  # cutout = data[5:310, 5:313] # Need to change N, A, B
  
  im = ax.imshow(data, cmap="jet", norm=colors.LogNorm()) # no longer needs np.flip because of WCS
  plt.plot(Nx, Ny, "x", color="black")
  plt.text(Nx - 12, Ny - 12, "N", color="black")
  plt.plot(Bx, By, "x", color="black")
  plt.text(Bx - 12, By - 12, "B", color="black")
  plt.plot(Ax, Ay, "x", color="black")
  plt.text(Ax - 12, Ay - 12, "A", color="black")

  # 320 in x, horizontal and 317 in y vertical
  r = 8
  center_y = 285
  center_x = 25

  indx = [center_y-r, center_y+r, center_x-r,center_x+r]
  ax.add_patch(Rectangle((center_y-r, center_x-r), 16, 16, facecolor="None", edgecolor="black", linewidth=2))
  
  box_data = data[indx[0]:indx[1], indx[2]:indx[3]]
  max_flux = np.max(box_data)
  std_dev = np.std(box_data)
  std = std_dev * np.array([3, 10, 20])
  print("std 1 = ", std_dev)
  print("std 3, 10, 20 = ", std)
  S_N = max_flux/std_dev
  print("S/N = ", S_N)

  labels = ["3$\sigma$", "10$\sigma$", "20$\sigma$"]
  color = ["red", "green", "blue"]

  CS = plt.contour(data[:, :], std, colors=["sienna", "crimson", "lavender"])
  plt.clabel(CS, fontsize=11.5)

  """
  for i in range(len(labels)):
     CS.collections[i].set_label(labels[i])
  """

  plt.xlabel("RA (J2000)")
  plt.ylabel("Dec (J2000)")
  plt.colorbar(im, fraction=0.046, pad=0.04, label="Logarithmic flux density [log($10^{-20}$ergs$^{-1}$cm$^{-2}\AA{}^{-1}$)]")
  plt.tight_layout()
  plt.legend()
  plt.savefig(savename)
  plt.show()

#line_map("OIII.fits", "OIII_map.png", "$[OIII]^*N_2$ emission map, $\lambda \in [5024, 5042]$ Å")
#line_map("H_alpa.fits", "H_alpha_map.png", "$H_\\alpha$ emission map, $\lambda \in [6594, 6608]$ Å")


# Preparing file for BBarolo:
def BBarolo_prep(wav_low, wav_up, newname):
  """
  hdu = fits.open(filename)
  hdu.info()
  data = hdu[0].data
  hdr = hdu[0].header
  """
  lambda0 = hdr["CRVAL3"] # extract value of the wavelengths lower boundary (smallest wavelength) [A]
  dlambda = hdr["CD3_3"] # wavelength resolution, or delta lambda if you will [A]
  len_wave = hdr["NAXIS3"] # Number of datapoints in the wavelength dimension

  wavelengths = np.linspace(lambda0, lambda0 + (len_wave-1)*dlambda, len_wave)

  lower_boundary=wav_low; 
  upper_boundary=wav_up # upper and lower boundary wavelengths in Angstrom.

  lower_indx=np.array(np.where(wavelengths>=lower_boundary))[0,0] # locating the index
  upper_indx=np.array(np.where(wavelengths<=upper_boundary))[-1,-1]

  extracted_data = data[lower_indx:upper_indx] # cuts the data in the wavelength dimension

  hdr["CRVAL3"] = lower_boundary # Updates the header value for the lower wavelength
  hdu = fits.PrimaryHDU(extracted_data, header=hdr) # creates an extension
  hdul = fits.HDUList([hdu]) # includes the extension into the extension list
  hdul.info() # prints info about the extensions (optional)
  hdul.writeto(newname)

#BBarolo_prep(5020, 5048, "bbarolo_OIII.fits")
#BBarolo_prep(6590, 6615, "bbarolo_H_alpha.fits") # 6590-6606, 6598
#BBarolo_prep(6612, 6629, "bbarolo_NII.fits")

#hdu = fits.open("bbarolo_OIII.fits")
#hdu.info()