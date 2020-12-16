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
import pyworkflow.object as pwobj


def formatTransformFile(ts, transformFilePath):
    """ This method takes a tilt series and the output transformation file path and creates an IMOD-based transform
    file in the location indicated. """

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
    """ This method takes an IMOD-based transformation matrix file path and returns a 3D matrix containing the
    transformation matrices for each tilt-image belonging to the tilt-series. """

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
    """ This method takes an IMOD-based fiducial model file path and returns a list containing the coordinates of each
    fiducial for each tilt-image belonging to the tilt-series. """

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
    """ This method takes an IMOD-based fiducial residual model file path and returns a list containing the coordinates
    and residual values of each fiducial for each tilt-image belonging to the tilt-series. Since IMOD establishes a
    float value for each coordinate the are parsed to int. """

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
    """ This method takes a list containing the angles for each tilt-image belonging to the tilt-series and writes the
    IMOD-based angle file at the given location. """

    angleList = []

    for ti in inputTs:
        angleList.append(ti.getTiltAngle())
    angleList.reverse()

    with open(angleFilePath, 'w') as f:
        f.writelines("%s\n" % angle for angle in angleList)


def formatAngleList(tltFilePath):
    """ This method takes an IMOD-based angle file path and returns a list containing the angles for each tilt-image
    belonging to the tilt-series. """

    angleList = []

    with open(tltFilePath) as f:
        tltText = f.read().splitlines()
        for line in tltText:
            angleList.append(float(line))
    angleList.reverse()

    return angleList


def format3DCoordinatesList(coordFilePath, xDim, yDim):
    """ This method takes an IMOD-based fiducial coordinates file path and returns a list containing each coordinate
    for each fiducial belonging to the tilt-series. """

    coorList = []

    with open(coordFilePath) as f:
        coorText = f.read().splitlines()
        for line in coorText:
            vector = line.split()
            coorList.append([float(vector[1]) - xDim / 2, float(vector[2]) - yDim / 2, float(vector[3])])

    return coorList


def formatDefocusFile(defocusFilePath):
    """ This method takes an IMOD-based ctf estimation file path and returns a list containing the defocus information
    from the estimation of the ctf of each tilt-image belonging to the tilt-series. """
    defocusTable = []

    with open(defocusFilePath) as f:
        defocusText = f.readlines()
        for index, line in enumerate(defocusText):
            vector = line.split()
            [float(i) for i in vector]
            if index == 0:
                "Remove last element from the first line (it contains the mode of the estimation run)"
                mode = vector.pop()
            defocusTable.append(vector)

        return defocusTable, mode


def formatDefocusAstigmatismFile(defocusAstigmatismFilePath):
    """ This method takes an IMOD-based ctf estimation file path and returns a list containing the defocus and
    astigmatism information from the estimation of the ctf of each tilt-image belonging to the tilt-series. """
    defocusAstigmatismTable = []

    with open(defocusAstigmatismFilePath) as f:
        defocusText = f.readlines()
        for index, line in enumerate(defocusText):
            vector = line.split()
            [float(i) for i in vector]
            if index == 0:
                "Get mode of estimation from the first line (remove rest of the information form this line)"
                mode = vector.pop()
            else:
                defocusAstigmatismTable.append(vector)

        return defocusAstigmatismTable, mode


def refactorCTFDefocusEstimationInfo(ctfInfoIMODTable):
    """ This method takes a table containing the information of an IMOD-based CTF estimation containing only defocus
    information (5 columns) and produces a new table containing the same information in a format readable for
    Scipion. """

    if len(ctfInfoIMODTable[0]) == 5:
        defocusUDict = {}

        for element in ctfInfoIMODTable:

            " Segregate information from range"
            for index in range(int(element[0]), int(element[1]) + 1):
                if index in defocusUDict.keys():
                    defocusUDict[index].append(pwobj.Float(element[4]))
                else:
                    defocusUDict[index] = [pwobj.Float(element[4])]

    else:
        raise Exception("Misleading file format, CTF estiation with no astigmatism should be 5 columns long")

    return defocusUDict


def refactorCTFDesfocusAstigmatismEstimationInfo(ctfInfoIMODTable):
    """ This method takes a table containing the information of an IMOD-based CTF estimation containing defocus and
    astigmatism information (5 columns) and produces a new table containing the same information in a format readable
    for Scipion. """

    if len(ctfInfoIMODTable[0]) == 7:
            defocusUDict = {}
            defocusVDict = {}
            defocusAngleDict = {}

            for element in ctfInfoIMODTable:

                " Segregate information from range"
                for index in range(int(element[0]), int(element[1]) + 1):

                    # Defocus U info
                    if index in defocusUDict.keys():
                        defocusUDict[index].append(pwobj.Float(element[4]))
                    else:
                        defocusUDict[index] = [pwobj.Float(element[4])]

                    # Defocus V info
                    if index in defocusVDict.keys():
                        defocusVDict[index].append(pwobj.Float(element[5]))
                    else:
                        defocusVDict[index] = [pwobj.Float(element[5])]

                    # Defocus angle info
                    if index in defocusAngleDict.keys():
                        defocusAngleDict[index].append(pwobj.Float(element[6]))
                    else:
                        defocusAngleDict[index] = [pwobj.Float(element[6])]

    else:
        raise Exception("Misleading file format, CTF estiation with astigmatism should be 7 columns long")

    return defocusUDict, defocusVDict, defocusAngleDict

    # # Defocus, astigmatism and phase shift estimation
    # if len(ctfInfoIMODTable[0]) == 8:
    #     defocusUDict = []
    #     defocusVDict = []
    #     defocusAngleDict = []
    #     phaseShiftList = []
    #
    #     return defocusUDict, defocusVDict, defocusAngleDict, phaseShiftList
