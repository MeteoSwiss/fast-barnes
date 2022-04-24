# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 14:51:40 2021

@author: Bruno Zürcher
"""

###############################################################################
### required code to make Basemap library work - refer to stacktrace ##########

import os
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")

###############################################################################

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.basemap import Basemap
import numpy as np

from ReadArray import readGridded2DArray, readCsvArray

###############################################################################

# one of [ "Naive", "Radius", "Convol", "OptConvol", "NaiveS2", "OptConvolS2" ]
method = "Convol"

# one of [ 0.25, 0.5, 1.0, 2.0, 4.0 ]
sigma = 1.0

# one of [ 54, 218, 872, 3490 ]
numPoints = 3490

# one of [ 4.0, 8.0, 16.0, 32.0, 64.0 ]
resolution = 32.0

# applies only to Convol interpolations: one of [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50 ]
numIter = 4

# applies only to Radius interpolation
radius = 3.717

# show mask
printAlphaChannel = True

###############################################################################

def getInterpolFieldData(dir, filename):
    # read grid data
    (a, y0, x0, stepY, stepX) = readGridded2DArray(dir + '/' + filename)

    # x-axis of grid
    gridX = np.arange(0.0, a.shape[1], 1.0)
    gridX = gridX*stepX + x0
    # y-axis of grid
    gridY = np.arange(0.0, a.shape[0], 1.0)
    gridY = gridY*stepY + y0
    # arrays holding x- and y-coordinates of grid points
    X, Y = np.meshgrid(gridX, gridY)

    return (X, Y, a)


# compose name of data file
dataFileName = method + "_" + str(1.0/resolution) + "_" + str(sigma) + "_" + str(numPoints)
if method == "Radius":
    dataFileName = dataFileName + "_" + str(radius)
if "Convol" in method:
    dataFileName = dataFileName + "_" + str(numIter)
dataFileName = dataFileName + ".gdat"


# resolution applies only for vector graphic; 'c' and 'l' too low
fig = plt.figure(figsize=(16, 12), edgecolor='w')
m = Basemap(projection='cyl', resolution='i',
    llcrnrlat=34.5, urcrnrlat=72, llcrnrlon=-26, urcrnrlon=49)

# default background
m.shadedrelief(scale=1.0)

# other background maps can be downloaded from
# http://www.naturalearthdata.com/downloads/10m-raster-data/10m-cross-blend-hypso/
# m.warpimage(image='C:...PATH.../Maps/HYP_HR/HYP_HR_SR_W_DR.tif')

m.drawcountries(linewidth=0.4)
m.drawparallels([20.0,30.0,40.0,50.0,60.0,70.0,80.0],labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,361.,10.),labels=[0,0,0,1])


# print isolines
(X, Y, Dat) = getInterpolFieldData('./output/grid', dataFileName)
# the levels of the plotted isolines
levels = np.arange(976,1026,2)
cs = m.contour(X, Y, Dat, levels, latlon=True, linewidths=1.0, colors='black')
# plot the line labels
plt.clabel(cs, levels[::2], inline=True, fmt='%d', fontsize=11, colors='black')


if printAlphaChannel:
    (X0, Y0, Dat0) = getInterpolFieldData('./output/grid', dataFileName)
    # define mask
    Dat2 = np.isnan(Dat0)
    Alpha = Dat2*0.35
    cmap = ListedColormap(['black', 'black'])
    m.imshow(Dat2, alpha=Alpha, cmap=cmap)


# read observation data and display it on map
pts = readCsvArray('./input/obs/PressQFF_202007271200_' + str(numPoints) + '.csv')
xcoor = pts[1][:]
ycoor = pts[0][:]
m.scatter(xcoor, ycoor, color='red', s=9, marker='.')

plt.show()

