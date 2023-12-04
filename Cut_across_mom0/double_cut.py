import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy import wcs
import sys

plt.rcParams['font.size'] = 25
# Define the radii you want to plot
radii = [2500, -2500, 3000, -3000, 3500, -3500, 4000, -4000, 4500, -4500, 5000, -5000]

# Open the FITS files
file = sys.argv[1]
file1 = sys.argv[2]
file2 = sys.argv[3]
hdu = fits.open(file)
F = np.loadtxt(file1, dtype=str, usecols=0)

hdu2 = fits.open(file2)


# Function to convert radii to kpc for x-axis labels
def radii_to_kpc(radii):
    return [r / 1000 for r in radii]


xticks = [-20, -10, 0, 10, 20]

# Read the data from the FITS files
data = hdu[0].data

data2 = hdu2[0].data
# Model
# Define the WCS parameters for the data
crpix = [350, 350]  # Original CRPIX values
crval = [0, 0]  # Original CRVAL values
cdelt = [50, -50]  # Original CDELT values
w1 = wcs.WCS(naxis=2)
w1.wcs.crpix = crpix
w1.wcs.crval = crval
w1.wcs.cdelt = cdelt
w1.wcs.ctype = ["LINEAR", "LINEAR"]
w1.wcs.cunit = ["parsec", "parsec"]

# Create subplots for the plots in a 3x4 matrix
fig, ax = plt.subplots(2, 6, sharex=False, sharey=False, figsize=(36, 12))

# Iterate over the radii and create plots
for (i, j) in enumerate(radii):
    row, col = divmod(i, 6)  # Calculate row and column position in the subplot matrix

    # x, y = w1.all_world2pix(j, 13000, 0)  # Convert world coordinates to pixel coordinates
    # x = np.round(x).astype(int)
    dx = j / cdelt[0]
    x = 350 + dx
    x = np.round(x).astype(int)
    print(x, 'model')

    # Create data cuts for the three datasets
    data_cut = data[:, x]
    # Generate radius values for the x-axis
    Y = np.arange(data_cut.shape[0])
    wcs_coords = [w1.wcs_pix2world(x, y, 1) for y in Y]
    radius = [wcs_coords[k][1] for k in Y]
    radius_kpc = radii_to_kpc(radius)
    # Plot the data for the current radius
    j_pos = abs(j)
    J = str(j_pos)
    # print(radius_kpc[120:580])
    ax[row, col].plot(radius_kpc[120:580], data_cut[120:580], label='model')
    ax[row, col].set_yticks([])

    # Save data to text files

    # data_model = np.column_stack((radius, data_cut))
    # np.savetxt("ngc551_s_surface_density_model" + J + ".txt", data_model, delimiter="\t")

# IFU

crpix = [33.2, 32.8]  # Original CRPIX values
crval = [0, 0]  # Original CRVAL values
cdelt = [339.291, -339.291]  # Original CDELT values
w1 = wcs.WCS(naxis=2)
w1.wcs.crpix = crpix
w1.wcs.crval = crval
w1.wcs.cdelt = cdelt
w1.wcs.ctype = ["LINEAR", "LINEAR"]
w1.wcs.cunit = ["parsec", "parsec"]

# Iterate over the radii and create plots
for (i, j) in enumerate(radii):
    row, col = divmod(i, 6)  # Calculate row and column position in the subplot matrix

    hdu1 = fits.open(F[i])
    data1 = hdu1[0].data
    # x, y = w1.all_world2pix(j, 13000, 0)  # Convert world coordinates to pixel coordinates
    dx = j / cdelt[0]
    x = 33.2 + dx
    print(x, 'before')
    x = np.round(x).astype(int)
    print('IFU', x)
    # Create data cuts for the three datasets
    data_cut1 = data1[:, x]

    # Process the data as needed (calculating min and normalization)
    data_cut1 = np.where(data_cut1 == 0, np.nan, data_cut1)
    min_data = np.nanmin(data_cut1)
    norm_data = data_cut1 - min_data

    # Generate radius values for the x-axis
    Y = np.arange(data_cut1.shape[0])
    wcs_coords = [w1.wcs_pix2world(x, y, 1) for y in Y]
    radius = [wcs_coords[k][1] for k in Y]
    radius_kpc = radii_to_kpc(radius)
    # print((radius_kpc))
    # Plot the data for the current radius
    if col % 2 == 0:
        col = col
        ax[row, col].plot(radius_kpc, norm_data, label='right cut')
    else:
        col = col - 1
        ax[row, col].plot(radius_kpc, norm_data, label='left cut')
    # ax[row, col].plot(radius_kpc, norm_data, label='IFU')
    j_pos = abs(j / 1000)
    J = str(j_pos)
    ax[row, col].set_title(J + ' ' + ' kpc ' + ' ' + 'Radial Cut' + ' ' + 'IFU', fontsize=25)

    # ax[row, col].legend()

    # Save data to text files

    # data_obs = np.column_stack((radius, norm_data))
    # np.savetxt("ngc551_s_surface_density_IFU" + J + ".txt", data_obs, delimiter="\t")

crpix = [118., 117.5]  # Original CRPIX values  # 232.1  350
crval = [0, 0]  # Original CRVAL values
cdelt = [209.858, -209.858]  # Original CDELT values  # 75.4 50
w1 = wcs.WCS(naxis=2)
w1.wcs.crpix = crpix
w1.wcs.crval = crval
w1.wcs.cdelt = cdelt
w1.wcs.ctype = ["LINEAR", "LINEAR"]
w1.wcs.cunit = ["parsec", "parsec"]

# Iterate over the radii and create plots
for i, j in enumerate(radii):
    row, col = divmod(i, 6)  # Calculate row and column position in the subplot matrix

    # x, y = w1.all_world2pix(j, 13000, 0)  # Convert world coordinates to pixel coordinates
    dx = j / cdelt[0]
    x = 118. + dx
    print(x, 'before')
    x = np.round(x).astype(int)
    print(x, 'spitzer')
    # Create data cuts for the three datasets
    data_cut2 = data2[:, x]
    data_cut2 = np.where(data_cut2 == 0, np.nan, data_cut2)
    min_data1 = np.nanmin(data_cut2)
    norm_data1 = data_cut2 - min_data1

    # Generate radius values for the x-axis
    Y = np.arange(data_cut2.shape[0])
    wcs_coords = [w1.wcs_pix2world(x, y, 1) for y in Y]
    radius = [wcs_coords[k][1] for k in Y]
    radius_kpc = radii_to_kpc(radius)
    # print(radius_kpc[60:180])
    # Plot the data for the current radius
    # ax[row, col].plot(radius_kpc, norm_data1, label='infrared')

    if col % 2 == 0:
        col = col + 1
        ax[row, col].plot(radius_kpc[60:180], norm_data1[60:180], label='right cut')
    else:
        col = col
        ax[row, col].plot(radius_kpc[60:180], norm_data1[60:180], label='left cut')
    j_pos = abs(j / 1000)
    J = str(j_pos)
    ax[row, col].set_title(J + ' ' + 'kpc' + ' ' + 'Radial Cut' + ' ' + '3.6', fontsize=25)

    # ax[row, col].legend()

    # Save data to text files
    # J = str(j)
    # data_infrared = np.column_stack((radius, norm_data1))
    # np.savetxt("ngc551_s_surface_density_infrared" + J + ".txt", data_infrared, delimiter="\t")
ax[0, 0].set_yticks([0, 200, 400, 600, 800, 1000])
ax[1, 0].set_yticks([0, 200, 400, 600, 800, 1000])
ax[1, 5].legend(frameon=False,fontsize=20)
fig.text(0.5, -0.01, 'Radius (Kpc)', ha='center', fontsize=25)
fig.text(-0.01, 0.5, '$M_\odot\: pc^{-2}$', va='center', rotation='vertical', fontsize=25)
# Adjust subplot layout
plt.tight_layout()

# Show the plots
plt.subplots_adjust(wspace=0., hspace=0.2)
plt.savefig('radial_double_cut_round.pdf', format='pdf', bbox_inches='tight')
plt.show()
