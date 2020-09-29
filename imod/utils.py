# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains utils functions for IMOD protocols
"""

import csv
import numpy as np


def formatTransformFile(ts, transformFilePath):
    """This method takes a tilt series and the output transformation file path
    and creates an IMOD-based transform file in the location indicated"""
    tsMatrixTransformList = []
    for ti in ts:
        transform = ti.getTransform().getMatrix().flatten()
        transformIMOD = [transform[0],
                         transform[1],
                         transform[3],
                         transform[4],
                         transform[2],
                         transform[5]]
        tsMatrixTransformList.append(transformIMOD)
    with open(transformFilePath, 'w') as f:
        csvW = csv.writer(f, delimiter='\t')
        csvW.writerows(tsMatrixTransformList)


def formatTransformationMatrix(matrixFile):
    """This method takes the IMOD-based transformation matrix file path
    and returns a 3D matrix containing the transformation matrices for each tilt-image
    belonging to the tilt-series"""
    with open(matrixFile, "r") as matrix:
        lines = matrix.readlines()
    numberLines = len(lines)
    frameMatrix = np.empty([3, 3, numberLines])
    i = 0
    for line in lines:
        values = line.split()
        frameMatrix[0, 0, i] = float(values[0])
        frameMatrix[1, 0, i] = float(values[2])
        frameMatrix[0, 1, i] = float(values[1])
        frameMatrix[1, 1, i] = float(values[3])
        frameMatrix[0, 2, i] = float(values[4])
        frameMatrix[1, 2, i] = float(values[5])
        frameMatrix[2, 0, i] = 0.0
        frameMatrix[2, 1, i] = 0.0
        frameMatrix[2, 2, i] = 1.0
        i += 1
    return frameMatrix


def formatFiducialList(fiducialFilePath):
    """This method takes the IMOD-based fiducial model file path and returns a list
    containing the coordinates of each fiducial for each tilt-image belonging to the
    tilt-series"""
    fiducialList = []
    with open(fiducialFilePath) as f:
        fiducialText = f.read().splitlines()
        for line in fiducialText:
            # Fix IMOD bug: columns merge in coordinates exceed 3 digits
            vector = line.replace('-', ' -').split()
            vector = [round(float(i)) for i in vector]
            fiducialList.append(vector)
    return fiducialList


def formatFiducialResidList(fiducialFilePath):
    """This method takes the IMOD-based fiducial residual model file path and returns a
    list containing the coordinates and residual values of each fiducial for each
    tilt-image belonging to the tilt-series. Since IMOD establishes a float value for each
    coordinate the are parsed to int"""
    fiducialResidList = []
    with open(fiducialFilePath) as f:
        fiducialText = f.read().splitlines()
        for line in fiducialText[1:]:
            # Fix IMOD bug: columns merge in coordinates exceed 3 digits
            vector = line.replace('-', ' -').split()
            fiducialResidList.append([round(float(vector[0])),
                                      round(float(vector[1])),
                                      int(vector[2]),
                                      float(vector[3]),
                                      float(vector[4])])
    return fiducialResidList


def formatAngleFile(inputTs, angleFilePath):
    """This method takes a list containing the angles for each tilt-image belonging to the tilt-series and writes the
    IMOD-based angle file at the given location"""
    angleList = []
    for ti in inputTs:
        angleList.append(ti.getTiltAngle())
    angleList.reverse()
    with open(angleFilePath, 'w') as f:
        f.writelines("%s\n" % angle for angle in angleList)


def formatAngleList(tltFilePath):
    """This method takes the IMOD-based angle file path and returns a list containing
    the angles for each tilt-image belonging to the tilt-series"""
    angleList = []
    with open(tltFilePath) as f:
        tltText = f.read().splitlines()
        for line in tltText:
            angleList.append(float(line))
    angleList.reverse()
    return angleList


def format3DCoordinatesList(coordFilePath, xDim, yDim):
    """This method takes the IMOD-based fiducial coordinates file path and returns a
    list containing each coordinate for each fiducial belonging to the tilt-series"""
    coorList = []
    with open(coordFilePath) as f:
        coorText = f.read().splitlines()
        for line in coorText:
            vector = line.split()
            coorList.append([float(vector[1]) - xDim / 2, float(vector[2]) - yDim / 2, float(vector[3])])
    return coorList
