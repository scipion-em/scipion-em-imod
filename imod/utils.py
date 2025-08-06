# *****************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
# *****************************************************************************
"""
This module contains utils functions for IMOD protocols
"""
import logging
import os
import csv
from typing import Union
import numpy as np

import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from imod import Plugin
from tomo.objects import TiltSeries

logger = logging.getLogger(__name__)

def generateIMODFiducialTextFile(landmarkModel, outputFilePath):
    """ This method takes a Scipion LandmarkModel object and
    generates a text file in the sepecified location in IMOD
    convention that contains the information of the position
    of each fiducial through the tilt-series.
    :param landmarkModel: landmarkModel Scipion object.
    :param outputFilePath: location where the output file must be saved.
    """

    infoTable = landmarkModel.retrieveInfoTable()
    outputLines = []

    for vector in infoTable:

        outputLines.append("\t%s\t%s\t%s\t%d\n" % (vector[3], vector[0],
                                                   int(float(vector[1])), int(float(vector[2])) - 1))

    with open(outputFilePath, 'w') as f:
        f.writelines(outputLines)


def generateIMODFidFile(protocol, landmarkModel):
    fiducialTextFile = pwutils.replaceExt(landmarkModel.getFileName(), "txt")
    generateIMODFiducialTextFile(landmarkModel, fiducialTextFile)

    fiducialModelGapPath = pwutils.replaceExt(landmarkModel.getFileName(), "fid")

    if not os.path.exists(fiducialModelGapPath):
        paramsPoint2Model = {
            'inputFile': fiducialTextFile,
            'outputFile': fiducialModelGapPath,
            'image': landmarkModel.getTiltSeries().getFirstItem().getFileName(),
            'size': landmarkModel.getSize() / 10
        }

        # -sp <value> parameter: generate sphere with radius <value>
        argsPoint2Model = "-InputFile %(inputFile)s " \
                          "-OutputFile %(outputFile)s " \
                          "-image %(image)s " \
                          "-zc -ci %(size)d"

        protocol.setStepsExecutor()
        Plugin.runImod(protocol, 'point2model',
                       argsPoint2Model % paramsPoint2Model)

    return fiducialModelGapPath


def formatAngleList(tltFilePath):
    """ This method takes an IMOD-based angle file path and
    returns a list containing the angles for each tilt-image
    belonging to the tilt-series. """

    angleList = []

    with open(tltFilePath) as f:
        tltText = f.read().splitlines()
        for line in tltText:
            angleList.append(float(line))

    return angleList


def generateDefocusIMODFileFromObject(ctfTomoSeries, defocusFilePath,
                                      isRelion=False, inputTiltSeries=None, presentAcqOrders=None):
    """ This method takes a ctfTomoSeries object a generate a
    defocus information file in IMOD formatting containing
    the same information in the specified location. """

    if inputTiltSeries is None:
        tiltSeries = ctfTomoSeries.getTiltSeries()
    else:
        tiltSeries = inputTiltSeries

    logger.info("Trying to generate defocus file at %s" % defocusFilePath)

    # Check if there is CTF estimation information as list
    if ctfTomoSeries.getFirstItem().hasEstimationInfoAsList() and not isRelion:

        logger.debug("Defocus file generated form a list.")
        flag = ctfTomoSeries.getIMODDefocusFileFlag()
        nEstimationsInRange = ctfTomoSeries.getNumberOfEstimationsInRange()

        if flag == 0:
            # Plain estimation

            logger.debug("Flag 0: Plain estimation.")
            defocusUDict = generateDefocusUDictionary(ctfTomoSeries)

            # Write IMOD defocus file
            with open(defocusFilePath, 'w') as f:
                lines = []
                pattern = "%d\t%d\t%.2f\t%.2f\t%d\n"

                for index in sorted(defocusUDict.keys()):
                    if index + nEstimationsInRange > len(defocusUDict.keys()):
                        break

                    # Dictionary keys is reversed because IMOD set indexes
                    # upside down Scipion (highest index for
                    # the tilt-image with the highest negative angle)
                    newLine = (pattern % (
                        index,
                        index + nEstimationsInRange,
                        round(tiltSeries[index].getTiltAngle(), 2),
                        round(tiltSeries[index + nEstimationsInRange].getTiltAngle(), 2),
                        # CONVERT DEFOCUS VALUE TO NANOMETERS (IMOD CONVENTION)
                        int(float(defocusUDict[index + nEstimationsInRange][0]) / 10)
                    ))

                    lines.append(newLine)

                if not isRelion:
                    # Finally, add flag to the first line of file
                    lines[0] = lines[0][0:-1] + "\t2\n"

                f.writelines(lines)

        elif flag == 1:
            # Astigmatism estimation
            logger.debug("Flag 1: Astigmatism estimation.")

            defocusUDict = generateDefocusUDictionary(ctfTomoSeries)
            defocusVDict = generateDefocusVDictionary(ctfTomoSeries)
            defocusAngleDict = generateDefocusAngleDictionary(ctfTomoSeries)

            # Write IMOD defocus file
            with open(defocusFilePath, 'w') as f:
                # This line is added at the beginning of the file in
                # order to match the IMOD defocus file format.
                lines = ["1\t0\t0.0\t0.0\t0.0\t3\n"]

                for index in sorted(defocusUDict.keys()):

                    if index + nEstimationsInRange > len(defocusUDict.keys()):
                        break

                    # Dictionary keys is reversed because IMOD set indexes
                    # upside down Scipion (highest index for
                    # the tilt-image with the highest negative angle)
                    newLine = ("%d\t%d\t%.2f\t%.2f\t%.1f\t%.1f\t%.2f\n" % (
                        index,
                        index + nEstimationsInRange,
                        round(tiltSeries[index].getTiltAngle(), 2),
                        round(tiltSeries[index + nEstimationsInRange].getTiltAngle(), 2),
                        # CONVERT DEFOCUS VALUE TO NANOMETERS (IMOD CONVENTION)
                        float(defocusUDict[index + nEstimationsInRange][0]) / 10,
                        # CONVERT DEFOCUS VALUE TO NANOMETERS (IMOD CONVENTION)
                        float(defocusVDict[index + nEstimationsInRange][0]) / 10,
                        float(defocusAngleDict[index + nEstimationsInRange][0]),
                    ))

                    lines.append(newLine)
                f.writelines(lines)

        elif flag == 4:
            # Phase-shift estimation
            logger.debug("Flag 4: Phase shift estimation.")

            defocusUDict = generateDefocusUDictionary(ctfTomoSeries)
            phaseShiftDict = generatePhaseShiftDictionary(ctfTomoSeries)

            # Write IMOD defocus file
            with open(defocusFilePath, 'w') as f:
                # This line is added at the beginning of the file in
                # order to match the IMOD defocus file format.
                lines = ["4\t0\t0.0\t0.0\t0.0\t3\n"]

                for index in sorted(defocusUDict.keys()):

                    if index + nEstimationsInRange > len(defocusUDict.keys()):
                        break

                    # Dictionary keys is reversed because IMOD set indexes
                    # upside down Scipion (highest index for
                    # the tilt-image with the highest negative angle)
                    newLine = ("%d\t%d\t%.2f\t%.2f\t%.1f\t%.2f\n" % (
                        index,
                        index + nEstimationsInRange,
                        round(tiltSeries[index].getTiltAngle(), 2),
                        round(tiltSeries[index + nEstimationsInRange].getTiltAngle(), 2),
                        # CONVERT DEFOCUS VALUE TO NANOMETERS (IMOD CONVENTION)
                        float(defocusUDict[index + nEstimationsInRange][0]) / 10,
                        float(phaseShiftDict[index + nEstimationsInRange][0]),
                    ))

                    lines.append(newLine)

                f.writelines(lines)

        elif flag == 5:
            # Astigmatism and phase shift estimation
            logger.debug("Flag 5: Astigmatism and phase shift.")

            defocusUDict = generateDefocusUDictionary(ctfTomoSeries)
            defocusVDict = generateDefocusVDictionary(ctfTomoSeries)
            defocusAngleDict = generateDefocusAngleDictionary(ctfTomoSeries)
            phaseShiftDict = generatePhaseShiftDictionary(ctfTomoSeries)

            # Write IMOD defocus file
            with open(defocusFilePath, 'w') as f:
                # This line is added at the beginning of the file in order
                # to match the IMOD defocus file format
                lines = ["5\t0\t0.0\t0.0\t0.0\t3\n"]

                for index in sorted(defocusUDict.keys()):

                    if index + nEstimationsInRange > len(defocusUDict.keys()):
                        break

                    # Dictionary keys is reversed because IMOD set indexes
                    # upside down Scipion (highest index for
                    # the tilt-image with the highest negative angle)
                    newLine = ("%d\t%d\t%.2f\t%.2f\t%.1f\t%.1f\t%.2f\t%.2f\n" % (
                        index,
                        index + nEstimationsInRange,
                        round(tiltSeries[index].getTiltAngle(), 2),
                        round(tiltSeries[index + nEstimationsInRange].getTiltAngle(), 2),
                        # CONVERT DEFOCUS VALUE TO NANOMETERS (IMOD CONVENTION)
                        float(defocusUDict[index + nEstimationsInRange][0]) / 10,
                        # CONVERT DEFOCUS VALUE TO NANOMETERS (IMOD CONVENTION)
                        float(defocusVDict[index + nEstimationsInRange][0]) / 10,
                        float(defocusAngleDict[index + nEstimationsInRange][0]),
                        float(phaseShiftDict[index + nEstimationsInRange][0])
                    ))

                    lines.append(newLine)
                f.writelines(lines)

        elif flag == 37:
            # Astigmatism, phase shift and cut-on frequency estimation
            logger.debug("Flag 37: Astigmatism, phase shift and cut-on frequency estimation.")

            defocusUDict = generateDefocusUDictionary(ctfTomoSeries)
            defocusVDict = generateDefocusVDictionary(ctfTomoSeries)
            defocusAngleDict = generateDefocusAngleDictionary(ctfTomoSeries)
            phaseShiftDict = generatePhaseShiftDictionary(ctfTomoSeries)
            cutOnFreqDict = generateCutOnFreqDictionary(ctfTomoSeries)

            # Write IMOD defocus file
            with open(defocusFilePath, 'w') as f:
                # This line is added at the beginning of the file in order
                # to match the IMOD defocus file format
                lines = ["37\t0\t0.0\t0.0\t0.0\t3\n"]

                for index in sorted(defocusUDict.keys()):

                    if index + nEstimationsInRange > len(defocusUDict.keys()):
                        break

                    # Dictionary keys is reversed because IMOD set indexes
                    # upside down Scipion (highest index for
                    # the tilt-image with the highest negative angle)
                    newLine = ("%d\t%d\t%.2f\t%.2f\t%.1f\t%.1f\t%.2f\t%.2f\t%.4f\n" % (
                        index,
                        index + nEstimationsInRange,
                        round(tiltSeries[index].getTiltAngle(), 2),
                        round(tiltSeries[index + nEstimationsInRange].getTiltAngle(), 2),
                        # CONVERT DEFOCUS VALUE TO NANOMETERS (IMOD CONVENTION)
                        float(defocusUDict[index + nEstimationsInRange][0]) / 10,
                        # CONVERT DEFOCUS VALUE TO NANOMETERS (IMOD CONVENTION)
                        float(defocusVDict[index + nEstimationsInRange][0]) / 10,
                        float(defocusAngleDict[index + nEstimationsInRange][0]),
                        float(phaseShiftDict[index + nEstimationsInRange][0]),
                        float(cutOnFreqDict[index + nEstimationsInRange][0])
                    ))

                    lines.append(newLine)
                f.writelines(lines)

        else:
            raise ValueError("Defocus file flag do not supported. Only "
                             "supported formats corresponding to flags 0, "
                             "1, 4, 5, and 37.")

    else:
        # There is no information available as list (not an IMOD CTF estimation)

        logger.info("Defocus file generated from defocus attributes.")

        with open(defocusFilePath, 'w') as f:
            lines = ["1\t0\t0.0\t0.0\t0.0\t3\n"]
            ind = 1
            for ti in tiltSeries:
                ctfTomo = ctfTomoSeries.getCtfTomoFromTi(ti)
                if ctfTomo:
                    if presentAcqOrders and ctfTomo.getAcquisitionOrder() not in presentAcqOrders:
                        continue
                    tiltAngle = ti.getTiltAngle()
                    newLine = ("%d\t%d\t%.2f\t%.2f\t%.1f\t%.1f\t%.2f\n" % (
                        ind,
                        ind,
                        tiltAngle,
                        tiltAngle,
                        # CONVERT DEFOCUS VALUE TO NANOMETERS (IMOD CONVENTION)
                        ctfTomo.getDefocusU() / 10,
                        # CONVERT DEFOCUS VALUE TO NANOMETERS (IMOD CONVENTION)
                        ctfTomo.getDefocusV() / 10,
                        ctfTomo.getDefocusAngle()))

                    lines.append(newLine)
                    ind += 1
            f.writelines(lines)


def generateDefocusUDictionary(ctfTomoSeries):
    """ This method generates a dictionary containing
     the defocus U estimation information from a ctfTomoSeries
    object. """

    defocusUDict = {}

    for ctfTomo in ctfTomoSeries:
        defocusInfoList = ctfTomo.getDefocusUList() if hasattr(ctfTomo, "_defocusUList") \
            else ctfTomo.getDefocusVList()
        defocusInfoList = defocusInfoList.split(",")

        index = ctfTomo.getIndex()

        defocusUDict[index] = defocusInfoList

    return defocusUDict


def generateDefocusVDictionary(ctfTomoSeries):
    """ This method generates a dictionary containing
    the defocus V estimation information from a ctfTomoSeries
    object. """

    defocusVDict = {}

    for ctfTomo in ctfTomoSeries:
        defocusInfoList = ctfTomo.getDefocusVList()
        defocusInfoList = defocusInfoList.split(",")

        index = ctfTomo.getIndex()

        defocusVDict[index] = defocusInfoList

    return defocusVDict


def generateDefocusAngleDictionary(ctfTomoSeries):
    """ This method generates a dictionary containing the
    defocus angle estimation information from a ctfTomoSeries
    object. """

    defocusAngleDict = {}

    for ctfTomo in ctfTomoSeries:
        defocusAngleList = ctfTomo.getDefocusAngleList()
        defocusAngleList = defocusAngleList.split(",")

        index = ctfTomo.getIndex()

        defocusAngleDict[index] = defocusAngleList

    return defocusAngleDict


def generatePhaseShiftDictionary(ctfTomoSeries):
    """ This method generates a dictionary containing
    the phase shift estimation information from a
    ctfTomoSeries object. """

    phaseShiftDict = {}

    for ctfTomo in ctfTomoSeries:
        phaseShiftList = ctfTomo.getPhaseShiftList()
        phaseShiftList = phaseShiftList.split(",")

        index = ctfTomo.getIndex()

        phaseShiftDict[index] = phaseShiftList

    return phaseShiftDict


def generateCutOnFreqDictionary(ctfTomoSeries):
    """ This method generates a dictionary containing
    the cut-on frequency estimation information from
    a ctfTomoSeries object. """

    cutOnFreqDict = {}

    for ctfTomo in ctfTomoSeries:
        cutOnFreqList = ctfTomo.getCutOnFreqList()
        cutOnFreqList = cutOnFreqList.split(",")

        index = ctfTomo.getIndex()

        cutOnFreqDict[index] = cutOnFreqList

    return cutOnFreqDict





