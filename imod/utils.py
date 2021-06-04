# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
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
import math
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


def generateIMODFiducialTextFile(landmarkModel, outputFilePath):
    """ This method takes a Scipion LandmarkModel object and generates a text file in the sepecified location in IMOD
    convention that contains the information of the position of each fiducial through the tilt-series.
    :param landmarkModel: landmarkModel Scipion object.
    :param outputFilePath: location where the output file must be saved.
    """

    infoTable = landmarkModel.retrieveInfoTable()
    outputLines = []

    for vector in infoTable:
        outputLines.append("\t%s\t%s\t%s\n" % (vector[0], vector[1], vector[2]))

    with open(outputFilePath, 'w') as f:
        f.writelines(outputLines)


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


def format3DCoordinatesList(coordFilePath):
    """ This method takes an IMOD-based fiducial coordinates file path and returns a list containing each coordinate
    for each fiducial belonging to the tilt-series. """

    coorList = []

    with open(coordFilePath) as f:
        coorText = f.read().splitlines()

        for line in coorText:
            if line != '':
                vector = line.split()
                coorList.append([float(vector[1]), float(vector[2]), float(vector[3])])

    return coorList


def getDefocusFileFlag(defocusFilePath):
    """ This method returns the flag that indicate the information contained in an IMOD defocus file. The flag value "is
    the sum of:
          1 if the file has astigmatism values
          2 if the astigmatism axis angle is in radians, not degrees
          4 if the file has phase shifts
          8 if the phase shifts are in radians, not degrees
         16 if tilt angles need to be inverted to match what the
             program expects (what Ctfplotter would produce)
             with the -invert option
         32 if the file has cut-on frequencies attenuating the phase
             at low frequencies"

             from https://bio3d.colorado.edu/imod/doc/man/ctfphaseflip.html """

    with open(defocusFilePath) as f:
        lines = f.readlines()

    " File contains only defocus information (no astigmatism, no phase shift, no cut-on frequency) "
    if len(lines) == 1:
        return 0
    elif len(lines[1].split()) == 5:
        return 0

    else:
        " File contains more information apart "
        return int(lines[0].split()[0])


def readDefocusFileAsTable(defocusFilePath):
    """ This method takes an IMOD-based ctf estimation file path and returns a table containing the CTF estimation
    information of each tilt-image (per line) belonging to the tilt-series. """
    defocusTable = []

    with open(defocusFilePath) as f:
        lines = f.readlines()

        for index, line in enumerate(lines):
            vector = line.split()
            [float(i) for i in vector]

            if index == 0 and len(lines) == 1:
                "CTF estimation is plain (no astigmatism, no phase shift, no cut-on frequency) and is the first line."
                "Remove last element from the first line (it contains the mode of the estimation run). This case "
                "considers that the estimation file only has one line. "
                vector.pop()
                defocusTable.append(vector)

            elif index == 0 and len(lines[1].split()) == 5:
                "CTF estimation is plain (no astigmatism, no phase shift, no cut-on frequency) and is the first line."
                "Remove last element from the first line (it contains the mode of the estimation run). "
                vector.pop()
                defocusTable.append(vector)

            elif index == 0 and len(lines[1].split()) != 5:
                "CTF estimation is not plain and is the first line."
                "Do not add this line to the table. Only contains flag and format info."
                pass

            else:
                "Any posterior line that is not the first one is added to the table ."
                defocusTable.append(vector)

        return defocusTable


def readCTFEstimationInfoFile(defocusFilePath, flag):
    """ This method takes an IMOD-based file path containing the information associated to a CTF estimation and
    produces a set of dictionaries containing the information of each parameter for each tilt-image belonging to the
    tilt-series. These dictionaries are readable information for Scipion, useful to generate the corresponding output
    CTFTomoSeries object. """

    # Read info as table
    ctfInfoIMODTable = readDefocusFileAsTable(defocusFilePath)

    if flag == 0:
        # Plain estimation
        return refactorCTFDefocusEstimationInfo(ctfInfoIMODTable)

    elif flag == 1:
        # Astigmatism estimation
        return refactorCTFDesfocusAstigmatismEstimationInfo(ctfInfoIMODTable)

    elif flag == 4:
        # Phase-shift estimation
        return refactorCTFDefocusPhaseShiftEstimationInfo(ctfInfoIMODTable)

    elif flag == 5:
        # Astigmatism and phase shift estimation
        return refactorCTFDefocusAstigmatismPhaseShiftEstimationInfo(ctfInfoIMODTable)

    elif flag == 37:
        # Astigmatism, phase shift and cut-on frequency estimation
        return refactorCTFDefocusAstigmatismPhaseShiftCutOnFreqEstimationInfo(ctfInfoIMODTable)

    else:
        raise Exception("Defocus file flag do not supported. Only supported formats corresponding to flags 0, "
                        "1, 4, 5, and 37.")


def refactorCTFDefocusEstimationInfo(ctfInfoIMODTable):
    """ This method takes a table containing the information of an IMOD-based CTF estimation containing only defocus
    information (5 columns) and produces a new dictionary containing the same information in a format readable for
    Scipion. Flag 0. """

    if len(ctfInfoIMODTable[0]) == 5:
        defocusUDict = {}

        for element in ctfInfoIMODTable:

            # Segregate information from range
            for index in range(int(element[0]), int(element[1]) + 1):
                if index in defocusUDict.keys():
                    defocusUDict[index].append(pwobj.Float(element[4]))
                else:
                    defocusUDict[index] = [pwobj.Float(element[4])]

    else:
        raise Exception("Misleading file format, CTF estimation with no astigmatism should be 5 columns long")

    return defocusUDict


def refactorCTFDesfocusAstigmatismEstimationInfo(ctfInfoIMODTable):
    """ This method takes a table containing the information of an IMOD-based CTF estimation containing defocus and
    astigmatism information (7 columns) and produces a set of dictionaries table containing the same information in a
    format readable for Scipion. Flag 1. """

    if len(ctfInfoIMODTable[0]) == 7:
        defocusUDict = {}
        defocusVDict = {}
        defocusAngleDict = {}

        for element in ctfInfoIMODTable:

            # Segregate information from range
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
        raise Exception("Misleading file format, CTF estimation with astigmatism should be 7 columns long")

    return defocusUDict, defocusVDict, defocusAngleDict


def refactorCTFDefocusPhaseShiftEstimationInfo(ctfInfoIMODTable):
    """ This method takes a table containing the information of an IMOD-based CTF estimation containing defocus,
    and phase shift information (6 columns) and produces a new set of dictionaries containing the same
    information in a format readable for Scipion. Flag 4. """

    if len(ctfInfoIMODTable[0]) == 6:
        defocusUDict = {}
        phaseShiftDict = {}

        for element in ctfInfoIMODTable:

            # Segregate information from range
            for index in range(int(element[0]), int(element[1]) + 1):

                # Defocus U info
                if index in defocusUDict.keys():
                    defocusUDict[index].append(pwobj.Float(element[4]))
                else:
                    defocusUDict[index] = [pwobj.Float(element[4])]

                # Phase shift info
                if index in phaseShiftDict.keys():
                    phaseShiftDict[index].append(pwobj.Float(element[5]))
                else:
                    phaseShiftDict[index] = [pwobj.Float(element[5])]

    else:
        raise Exception("Misleading file format, CTF estimation with astigmatism and phase shift should be 6 columns "
                        "long")

    return defocusUDict, phaseShiftDict


def refactorCTFDefocusAstigmatismPhaseShiftEstimationInfo(ctfInfoIMODTable):
    """ This method takes a table containing the information of an IMOD-based CTF estimation containing defocus,
    astigmatism and phase shift information (8 columns) and produces a new set of dictionaries containing the same
    information in a format readable for Scipion. Flag 5. """

    if len(ctfInfoIMODTable[0]) == 8:
        defocusUDict = {}
        defocusVDict = {}
        defocusAngleDict = {}
        phaseShiftDict = {}

        for element in ctfInfoIMODTable:

            # Segregate information from range
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

                # Phase shift info
                if index in phaseShiftDict.keys():
                    phaseShiftDict[index].append(pwobj.Float(element[7]))
                else:
                    phaseShiftDict[index] = [pwobj.Float(element[7])]

    else:
        raise Exception("Misleading file format, CTF estiation with astigmatism and phase shift should be 8 columns "
                        "long")

    return defocusUDict, defocusVDict, defocusAngleDict, phaseShiftDict


def refactorCTFDefocusAstigmatismPhaseShiftCutOnFreqEstimationInfo(ctfInfoIMODTable):
    """ This method takes a table containing the information of an IMOD-based CTF estimation containing defocus,
    astigmatism, phase shift information and cut-on frequency (8 columns) and produces a new set of dictionaries
    containing the same information in a format readable for Scipion. Flag 37. """

    if len(ctfInfoIMODTable[0]) == 9:
        defocusUDict = {}
        defocusVDict = {}
        defocusAngleDict = {}
        phaseShiftDict = {}
        cutOnFreqDict = {}

        for element in ctfInfoIMODTable:

            # Segregate information from range
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

                # Phase shift info
                if index in phaseShiftDict.keys():
                    phaseShiftDict[index].append(pwobj.Float(element[7]))
                else:
                    phaseShiftDict[index] = [pwobj.Float(element[7])]

                # Cut-on frequency info
                if index in cutOnFreqDict.keys():
                    cutOnFreqDict[index].append(pwobj.Float(element[8]))
                else:
                    cutOnFreqDict[index] = [pwobj.Float(element[8])]

    else:
        raise Exception("Misleading file format, CTF estiation with astigmatism, phase shift and cut-on frequency "
                        "should be 8 columns long")

    return defocusUDict, defocusVDict, defocusAngleDict, phaseShiftDict, cutOnFreqDict


def generateDefocusIMODFileFromObject(ctfTomoSeries, defocusFilePath):
    """ This methods takes a ctfTomoSeries object a generate a defocus information file in IMOD formatting containing
    the same information in the specified location. """

    tiltSeries = ctfTomoSeries.getTiltSeries()

    # Check if there is CTF estimation information as list
    if ctfTomoSeries.getFirstItem().hasEstimationInfoAsList():

        flag = ctfTomoSeries.getIMODDefocusFileFlag()

        if flag == 0:
            # Plain estimation
            defocusUDict = generateDefocusUDictionary(ctfTomoSeries)

            # Write IMOD defocus file
            with open(defocusFilePath, 'w') as f:
                lines = []

                for index in defocusUDict.keys():

                    if index + ctfTomoSeries.getNumberOfEstimationsInRange() > len(defocusUDict.keys()):
                        break

                    # Dictionary keys is reversed because IMOD set indexes upside down Scipion (highest index for
                    # the tilt-image with the highest negative angle)
                    newLine = ("%d\t%d\t%.2f\t%.2f\t%d\n" % (
                        index,
                        index + ctfTomoSeries.getNumberOfEstimationsInRange(),
                        round(tiltSeries[index + ctfTomoSeries.getNumberOfEstimationsInRange()].getTiltAngle(), 2),
                        round(tiltSeries[index].getTiltAngle(), 2),
                        int(float(defocusUDict[index][0]))
                    ))

                    lines = [newLine] + lines

                # Finally, add flag to the first line of file
                lines[0] = lines[0][0:-1] + "\t2\n"

                f.writelines(lines)

        elif flag == 1:
            # Astigmatism estimation

            defocusUDict = generateDefocusUDictionary(ctfTomoSeries)
            defocusVDict = generateDefocusVDictionary(ctfTomoSeries)
            defocusAngleDict = generateDefocusAngleDictionary(ctfTomoSeries)

            # Write IMOD defocus file
            with open(defocusFilePath, 'w') as f:
                lines = []

                for index in defocusUDict.keys():

                    if index + ctfTomoSeries.getNumberOfEstimationsInRange() > len(defocusUDict.keys()):
                        break

                    # Dictionary keys is reversed because IMOD set indexes upside down Scipion (highest index for
                    # the tilt-image with the highest negative angle)
                    newLine = ("%d\t%d\t%.2f\t%.2f\t%.1f\t%.1f\t%.2f\n" % (
                        index,
                        index + ctfTomoSeries.getNumberOfEstimationsInRange(),
                        round(tiltSeries[index + ctfTomoSeries.getNumberOfEstimationsInRange()].getTiltAngle(), 2),
                        round(tiltSeries[index].getTiltAngle(), 2),
                        float(defocusUDict[index][0]),
                        float(defocusVDict[index][0]),
                        float(defocusAngleDict[index][0]),
                    ))

                    lines = [newLine] + lines

                # This line is added at the beginning of the file in order to match the IMOD defocus file format.
                lines = ["1\t0\t0.0\t0.0\t0.0\t3\n"] + lines

                f.writelines(lines)

        elif flag == 4:
            # Phase-shift estimation

            defocusUDict = generateDefocusUDictionary(ctfTomoSeries)
            phaseShiftDict = generatePhaseShiftDictionary(ctfTomoSeries)

            # Write IMOD defocus file
            with open(defocusFilePath, 'w') as f:
                lines = []

                for index in defocusUDict.keys():

                    if index + ctfTomoSeries.getNumberOfEstimationsInRange() > len(defocusUDict.keys()):
                        break

                    # Dictionary keys is reversed because IMOD set indexes upside down Scipion (highest index for
                    # the tilt-image with the highest negative angle)
                    newLine = ("%d\t%d\t%.2f\t%.2f\t%.1f\t%.2f\n" % (
                        index,
                        index + ctfTomoSeries.getNumberOfEstimationsInRange(),
                        round(tiltSeries[index + ctfTomoSeries.getNumberOfEstimationsInRange()].getTiltAngle(), 2),
                        round(tiltSeries[index].getTiltAngle(), 2),
                        float(defocusUDict[index][0]),
                        float(phaseShiftDict[index][0]),
                    ))

                    lines = [newLine] + lines

                # This line is added at the beginning of the file in order to match the IMOD defocus file format.
                lines = ["4\t0\t0.0\t0.0\t0.0\t3\n"] + lines

                f.writelines(lines)

        elif flag == 5:
            # Astigmatism and phase shift estimation

            defocusUDict = generateDefocusUDictionary(ctfTomoSeries)
            defocusVDict = generateDefocusVDictionary(ctfTomoSeries)
            defocusAngleDict = generateDefocusAngleDictionary(ctfTomoSeries)
            phaseShiftDict = generatePhaseShiftDictionary(ctfTomoSeries)

            # Write IMOD defocus file
            with open(defocusFilePath, 'w') as f:
                lines = []

                for index in defocusUDict.keys():

                    if index + ctfTomoSeries.getNumberOfEstimationsInRange() > len(defocusUDict.keys()):
                        break

                    # Dictionary keys is reversed because IMOD set indexes upside down Scipion (highest index for
                    # the tilt-image with the highest negative angle)
                    newLine = ("%d\t%d\t%.2f\t%.2f\t%.1f\t%.1f\t%.2f\t%.2f\n" % (
                        index,
                        index + ctfTomoSeries.getNumberOfEstimationsInRange(),
                        round(tiltSeries[index + ctfTomoSeries.getNumberOfEstimationsInRange()].getTiltAngle(), 2),
                        round(tiltSeries[index].getTiltAngle(), 2),
                        float(defocusUDict[index][0]),
                        float(defocusVDict[index][0]),
                        float(defocusAngleDict[index][0]),
                        float(phaseShiftDict[index][0])
                    ))

                    lines = [newLine] + lines

                # This line is added at the beginning of the file in order to match the IMOD defocus file format
                lines = ["5\t0\t0.0\t0.0\t0.0\t3\n"] + lines

                f.writelines(lines)

        elif flag == 37:
            # Astigmatism, phase shift and cut-on frequency estimation

            defocusUDict = generateDefocusUDictionary(ctfTomoSeries)
            defocusVDict = generateDefocusVDictionary(ctfTomoSeries)
            defocusAngleDict = generateDefocusAngleDictionary(ctfTomoSeries)
            phaseShiftDict = generatePhaseShiftDictionary(ctfTomoSeries)
            cutOnFreqDict = generateCutOnFreqDictionary(ctfTomoSeries)

            # Write IMOD defocus file
            with open(defocusFilePath, 'w') as f:
                lines = []

                for index in defocusUDict.keys():

                    if index + ctfTomoSeries.getNumberOfEstimationsInRange() > len(defocusUDict.keys()):
                        break

                    # Dictionary keys is reversed because IMOD set indexes upside down Scipion (highest index for
                    # the tilt-image with the highest negative angle)
                    newLine = ("%d\t%d\t%.2f\t%.2f\t%.1f\t%.1f\t%.2f\t%.2f\t%.4f\n" % (
                        index,
                        index + ctfTomoSeries.getNumberOfEstimationsInRange(),
                        round(tiltSeries[index + ctfTomoSeries.getNumberOfEstimationsInRange()].getTiltAngle(), 2),
                        round(tiltSeries[index].getTiltAngle(), 2),
                        float(defocusUDict[index][0]),
                        float(defocusVDict[index][0]),
                        float(defocusAngleDict[index][0]),
                        float(phaseShiftDict[index][0]),
                        float(cutOnFreqDict[index][0])
                    ))

                    lines = [newLine] + lines

                # This line is added at the beginning of the file in order to match the IMOD defocus file format
                lines = ["37\t0\t0.0\t0.0\t0.0\t3\n"] + lines

                f.writelines(lines)

        else:
            raise Exception("Defocus file flag do not supported. Only supported formats corresponding to flags 0, "
                            "1, 4, 5, and 37.")

    else:
        # There is no information available as list (not an IMOD CTF estimation)

        with open(defocusFilePath, 'w') as f:
            lines = []

            # CtfTomoSeries is iterated inversely because IMOD set indexes upside down Scipion (highest index for
            # the tilt-image with the highest negative angle)
            for ctfTomo in ctfTomoSeries:
                index = ctfTomo.getIndex()

                newLine = ("%d\t%d\t%.2f\t%.2f\t%.1f\t%.1f\t%.2f\n" % (
                    index,
                    index,
                    tiltSeries[index].getTiltAngle(),
                    tiltSeries[index].getTiltAngle(),
                    ctfTomo.getDefocusU(),
                    ctfTomo.getDefocusV(),
                    ctfTomo.getDefocusAngle())
                           )

                lines = [newLine] + lines

            f.writelines(lines)


def generateDefocusUDictionary(ctfTomoSeries):
    """ This method generates a dictionary containing the defocus U estimation information from a ctfTomoSeries
    object. """

    defocusUDict = {}

    for ctfTomo in ctfTomoSeries:
        defocusInfoList = ctfTomo.getDefocusUList() if hasattr(ctfTomo, "_defocusUList") \
            else ctfTomo.getDefocusVList()
        defocusInfoList = defocusInfoList.split(",")

        index = ctfTomo.getIndex().get()

        defocusUDict[index] = defocusInfoList

    return defocusUDict


def generateDefocusVDictionary(ctfTomoSeries):
    """ This method generates a dictionary containing the defocus V estimation information from a ctfTomoSeries
    object. """

    defocusVDict = {}

    for ctfTomo in ctfTomoSeries:
        defocusInfoList = ctfTomo.getDefocusVList()
        defocusInfoList = defocusInfoList.split(",")

        index = ctfTomo.getIndex().get()

        defocusVDict[index] = defocusInfoList

    return defocusVDict


def generateDefocusAngleDictionary(ctfTomoSeries):
    """ This method generates a dictionary containing the defocus angle estimation information from a ctfTomoSeries
    object. """

    defocusAngleDict = {}

    for ctfTomo in ctfTomoSeries:
        defocusAngleList = ctfTomo.getDefocusAngleList()
        defocusAngleList = defocusAngleList.split(",")

        index = ctfTomo.getIndex().get()

        defocusAngleDict[index] = defocusAngleList

    return defocusAngleDict


def generatePhaseShiftDictionary(ctfTomoSeries):
    """ This method generates a dictionary containing the phase shift estimation information from a ctfTomoSeries
    object. """

    phaseShiftDict = {}

    for ctfTomo in ctfTomoSeries:
        phaseShiftList = ctfTomo.getPhaseShiftList()
        phaseShiftList = phaseShiftList.split(",")

        index = ctfTomo.getIndex().get()

        phaseShiftDict[index] = phaseShiftList

    return phaseShiftDict


def generateCutOnFreqDictionary(ctfTomoSeries):
    """ This method generates a dictionary containing the cut-on frequency estimation information from a ctfTomoSeries
    object. """

    cutOnFreqDict = {}

    for ctfTomo in ctfTomoSeries:
        cutOnFreqList = ctfTomo.getCutOnFreqList()
        cutOnFreqList = cutOnFreqList.split(",")

        index = ctfTomo.getIndex().get()

        cutOnFreqDict[index] = cutOnFreqList

    return cutOnFreqDict


def formatGoldBead3DCoordinatesList(coordFilePath):
    """This method takes the IMOD-based gold bead 3D coordinates obtained with find3dbeads program file path and
    returns a list containing each coordinate for each bead belonging to the tilt-series"""

    coorList = []

    with open(coordFilePath) as f:
        coorText = f.read().splitlines()

        for line in coorText:
            vector = line.split()
            coorList.append([float(vector[0]), float(vector[1]), float(vector[2])])

    return coorList


def calculateRotationAngleFromTM(ts):
    """ This method calculates que average tilt image rotation angle from its associated transformation matrix."""
    avgRotationAngle = 0

    for ti in ts:
        tm = ti.getTransform().getMatrix()
        cosRotationAngle = tm[0][0]
        sinRotationAngle = tm[1][0]
        avgRotationAngle += math.degrees(math.atan(sinRotationAngle/cosRotationAngle))

    avgRotationAngle = avgRotationAngle / ts.getSize()

    return avgRotationAngle
