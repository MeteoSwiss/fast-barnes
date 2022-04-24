# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 16:06:11 2021

@author: Bruno Zürcher
"""

import struct
import numpy as np


def readGridded2DArray(filename):
    """
    Reads an array of (Java-endian) doubles from a binary file, whereas integer height
    and width of array are initially read from file, as well as the grid specific
    parameters y0, x0, stepY, stepX.
    :param filename: The file name.
    :return: The grid data array, y- and x-coordinate of grid, step in y- and x-direction.
    """
    f = open(filename, 'rb')

    # read dimensions
    block = f.read(8)
    (h, w) = struct.unpack('>2l', block)

    # read grid specifics
    block = f.read(32)
    (y0, x0, stepY, stepX) = struct.unpack('>4d', block)

    packStruct = struct.Struct('>' + str(w) + 'd')
    numBytes = 8*w
    container = []
    for j in range(h):
        block = f.read(numBytes)
        data = packStruct.unpack(block)
        container.append(data)

    f.close()
    return (np.array(container), y0, x0, stepY, stepX)


def readCsvArray(filename):
    """
    Reads the csv-file of observation data given linewise by their latitude,
    longitude and observation value.
    The file header line defines the number of csv-lines and the number of values per line.
    """
    f = open(filename, 'r')

    # read height and width
    (h, w) = map(int, f.readline().split(','))

    arr = np.zeros((w, h))
    count = 0
    for line in f:
        res = line.split(',')
        arr[0][count] = float(res[0])
        arr[1][count] = float(res[1])
        arr[2][count] = float(res[2])
        count += 1

    f.close()
    return arr
