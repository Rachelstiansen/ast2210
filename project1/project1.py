# Importing packages:
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.interpolate import interp1d
import scipy.constants as const
from tqdm import trange

from scipy.optimize import curve_fit

# Making fancier plots:
import seaborn as sns
sns.color_palette("bright")
sns.set_theme()
sns.set_style("darkgrid")

# Loading data to script
idata = np.load("idata_square.npy") # spectral data (intensity for (x, y, wav))
spect_pos = np.load("spect_pos.npy") # wav value in [Å]

# Points to observe + arrays needed for later
A = np.array([49, 197]); B = np.array([238, 443]); 
C = np.array([397, 213]); D = np.array([466, 52])

point_names = np.array(["A", "B", "C", "D"])
points = np.array([A, B, C, D])

plt.rcParams["image.origin"] = "lower" # Fixes the origin at the bottom left when plotting 
plt.rcParams["image.cmap"] = "hot" # For prettier colormap
plt.rcParams.update({"font.size": 14})

# Finding spectra averaged over whole region
avg = np.average(idata, axis=1)
avg_tot = np.average(avg, axis=0)
# print(avg_tot)

# Useful for sub field of view
x1 = 525; y1 = 325

corner = (x1-1, y1-1)

height = 100
width = 150

# Rectangle for visualizing the sub field of view:
rect = Rectangle(corner, width, height, linewidth=2, edgecolor="white", facecolor="none")


def spectra_line_plot(point, point_name, compare_to_avg=None, gauss_fit=None, show=None, save=None):
    """
    Function for plotting the spectra for points A, B, C and D. It also plots the
    spectra averages over the whole region and ...
    """
    fig, ax = plt.subplots()

    x = point[0]; y = point[1]

    wavelength_spectrum = idata[y-1, x-1, :]

    ax.plot(spect_pos, wavelength_spectrum, ls="--", lw=1, color="red", marker="x", label="Intensity")
    
    ax.set_title(f"Spectra at point {point_name} at ({x}, {y}) px")
    ax.set_xlabel(r"Wavelength $\lambda_i$ [Å]")
    ax.set_ylabel(r"Intensity [$erg^{-1} cm^{-2} sr^{-1} Hz^{-1}$]")

    # Finding average and comparing:
    if compare_to_avg == True:
        ax.plot(spect_pos, avg_tot, color="b", label="Average")

    if gauss_fit == True:
        params, xnew, cov = gaussian_fit(point, average=False)
        a, b, c, d = params
        g = gaussian(xnew, a, b, c, d)
        ax.plot(xnew, g, color="g", label="Gaussian fit")
        #ax.errorbar(spect_pos, wavelength_spectrum)
    
    ax.legend()
    fig.tight_layout()  

    if save == True:
        plt.savefig(f"Spectra_point_{point_name}.png")

    if show == True:
        plt.show()


def intensity_plot(wav_indx, point, add_points=None, add_rect=None, show=None, save=None):
    """
    Function for plotting the intensity image of wavelenght location (wav_indx).
    If add_points=True then the observation points are marked with a circle on the image.
    If add_rect=True then the sub field of view is also added to the image.
    """
    intensity_data = idata[:, :, wav_indx]
    wav_val = spect_pos[wav_indx]

    fig, ax = plt.subplots()

    ax.grid(False)
    im = ax.imshow(intensity_data)
    cbar = fig.colorbar(im)

    cbar.ax.set_ylabel(r"Intensity [$erg^{-1} cm^{-2} sr^{-1} Hz^{-1}$]")
    ax.set_title(f"Intensity image of wavelength $\lambda$ = {wav_val} Å")
    ax.set_xlabel(f"x [idx]")
    ax.set_ylabel(f"y [idx]")

    fig.tight_layout()
    
    if add_points == True:
        for i in range(len(point)):
            x = point[i, 0]; y = point[i, 1]
            circle = plt.Circle((x-1, y-1), radius=15, edgecolor="white", linewidth=2, fill=False)
            ax.add_patch(circle)
            ax.set_aspect("equal")
            ax.text(x+9, y+19, point_names[i], color="white", fontweight="bold")

    if add_rect == True:
        ax.add_patch(rect)

    if save == True:
            plt.savefig(f"Intensity_image_{wav_val}Å.png")
    
    if show == True:
        plt.show()
    

x1, y1 = rect.get_xy()
x2 = x1 + rect.get_width()
y2 = y1 + rect.get_height()

slice_x = slice(x1,x2)
slice_y = slice(y1,y2)

idata_cut = idata[slice_y, slice_x, :]

def sub_field_intensity_plot(wav_indx, points, add_points=None, show=None, save=None):
    """
    Function for plotting the intensity image of the sub field of view.
    """    
    # Sclicing out a small grid
    wav_val = spect_pos[wav_indx]
    intensity_data = idata_cut[:,:, wav_indx]

    fig, ax = plt.subplots()
    
    ax.grid(False)
    im = ax.imshow(intensity_data)
    cbar = fig.colorbar(im)

    cbar.ax.set_ylabel(r"Intensity [$erg^{-1} cm^{-2} sr^{-1} Hz^{-1}$]")
    ax.set_title(f"Intensity image of sub FOV for wavelenght $\lambda$ = {wav_val} Å")
    ax.set_xlabel(f"x [idx]")
    ax.set_ylabel(f"y [idx]")
    """
    if add_points == True:
        for i in range(len(points)):
            x = points[i, 0]; y = points[i, 1]
            circle = plt.Circle((x-1, y-1), radius=15, edgecolor="pink", fill=False)
            ax.add_patch(circle)
            ax.set_aspect("equal")
    """  
    fig.tight_layout()

    if save == True:
            plt.savefig(f"Sub_FOV_intensity.png")

    if show == True:
        plt.show()

def gaussian(lamda, a, b, c, d):
    return a * np.exp(-(lamda-b)**2/(2*c**2)) + d

def gaussian_fit(point=None, average=None):
    """
    Function for estimating a gaussian fit to the spectra for points A, B,
    C and D. First we try to estimate the parameters a, b, c and d from the data,
    then we insert them in the function curve_fit from scipy which gives us even
    more accurate parameters as output. These are then used in the function above (gaussian)
    
    We estimate the gaussian fit both for the average spectra and for the spectra for
    each point.
    """

    if average == True: 
        wavelength_spectrum = avg_tot

    else:
        x = point[0]; y = point[1]
        wavelength_spectrum = idata[y-1, x-1, :]

    # Finding a, b, c, d:
    est_a = np.min(wavelength_spectrum) - np.max(wavelength_spectrum)
    est_b = spect_pos[np.argmin(wavelength_spectrum)] #[np.where(min(wavelength_spectrum) == wavelength_spectrum)[0][0]]
    est_c = 0.05
    est_d = np.max(wavelength_spectrum)
    
    tot_est_params = np.array([est_a, est_b, est_c, est_d])

    params, covariance = curve_fit(gaussian, spect_pos, wavelength_spectrum, tot_est_params)

    x_new = np.linspace(spect_pos[0], spect_pos[-1], 100)

    return params, x_new, covariance


def central_wavelength_est(show=None, save=None):

    params, x_new, cov = gaussian_fit(average=True)
    a, b, c, d = params
    g = gaussian(x_new, a, b, c, d)

    est_central_wav = b
    error = np.sqrt(np.diag(cov))
    error_wav = error[1]
    print(f" error wav = {error_wav}")


    if save == True:
        fig, ax = plt.subplots()
        ax.plot(spect_pos, avg_tot, "x--", color="b", label="Average")
        ax.plot(x_new, g, label="Gaussian fit")
        ax.plot(b, np.min(g), "o")
        ax.text(b, np.min(g)-20, f"{(round(b, 3))} Å " , fontsize=11)
        ax.set_title(f"Estimate of central wavelength of Fe I")
        ax.set_xlabel(f"$\lambda$ [Å]")
        ax.set_ylabel(r"Intensity [$erg^{-1} cm^{-2} sr^{-1} Hz^{-1}$]")
        ax.set_ylim(230, 450)

        ax.legend()
        fig.get_tight_layout()
        plt.savefig("est_central_wav.png")

    if show == True:
        plt.show()

    return est_central_wav


def est_doppler_velocity(point):

    c_vel = const.speed_of_light
    Å = const.angstrom

    lamda_em = 6173 #central_wavelength_est(show=False, save=False)# m/s
    
    params, xnew, cov = gaussian_fit(point)
    a, b, c, d = params

    lamda_obs = b # m/s
    #print("b = ", lamda_obs)
    d_lamda = lamda_obs - lamda_em
    # print("lamda_em = ", lamda_em)
    doppler_vel = (d_lamda / lamda_em) * c_vel
    error = np.sqrt(np.diag(cov))*(c_vel/lamda_em) # sqrt(np.diag())????
    error = error[1]

    return doppler_vel, error

def full_doppler_map(wav_indx):

    X = len(idata[0, :, wav_indx])
    Y = len(idata[:, 0, wav_indx])

    title = "full_FOV"
    
    save_vel = np.zeros((Y, X))

    for i in trange(X):
        for j in trange(Y):
            point = np.array([i, j])
            doppler_vel, error = est_doppler_velocity(point)
            save_vel[j, i] = doppler_vel
    
    np.save(f"dop_map_{title}.npy", save_vel)

def plot_doppler_map(vel_data, sub_FOV=None):
    
    if sub_FOV == True:
        x0 = 525; y0 = 325
        x1 = x0 + 150; y1 = y0 + 100

        vel = vel_data[y0:y1+1, x0:x1+1]
        title = "sub FOV"
    else:
        vel = vel_data
        title = "full FOV"

    fig, ax = plt.subplots()
    ax.grid(False)

    im = ax.imshow(vel, cmap='seismic')
    cbar = fig.colorbar(im)

    cbar.ax.set_ylabel(r'Velocity [m/s]')
    ax.set_title(f"Doppler map of {title}")
    ax.set_xlabel("x [idx]")
    ax.set_ylabel("y [idx]")

    fig.tight_layout()
    plt.savefig(f"dop_map_{title}.png")
    plt.show()

#Only call the following function if this is the first time you run the script.
#full_doppler_map(0)

vel = np.load("dop_map_full_FOV.npy")

def bonus():
    intensity = idata[:, :, 4].flatten()
    doppler_vel = vel.flatten()
    plt.scatter(doppler_vel, intensity, s=1)
    plt.title("Intensity vs Doppler velocity")
    plt.ylabel(r"Intensity [$erg^{-1} cm^{-2} sr^{-1} Hz^{-1}$]")
    plt.xlabel("Doppler velocity [m/s]")
    plt.savefig("scatter.png")
    plt.show()

bonus()
    

"""
FUNCTION CALLING:
"""

"""
for i in range(len(points)):
    spectra_line_plot(points[i], point_names[i], gauss_fit=True, compare_to_avg=True, show=True, save=True)

"""
# intensity_plot(4, points, add_points=True, add_rect=True, show=True, save=True)

# sub_field_intensity_plot(4, points, add_points=True, show=True, save=True)

"""
for i in range(len(points)):
    gaussian_fit(points=True, average=True)
"""

# central_wavelength = central_wavelength_est(show=True, save=True)
# print(central_wavelength)


for i in range(len(points)):
    doppler_vel, error = est_doppler_velocity(points[i])
    print(f"doppler velocity of point {point_names[i]} = {doppler_vel:.3f} +- {error:.3f} m/s")

# plot_doppler_map(vel, sub_FOV=False)
