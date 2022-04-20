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

import os

from pwem.protocols import EMProtocol
from pyworkflow.object import Set, CsvList, Integer, Pointer
from pyworkflow.utils import path

from tomo.protocols import ProtTomoBase
from tomo.protocols.protocol_base import ProtTomoImportFiles
from tomo.objects import SetOfTiltSeries, SetOfTomograms, SetOfCTFTomoSeries, CTFTomoSeries, CTFTomo

from imod import utils
from imod import Plugin


class ProtImodBase(ProtTomoImportFiles, EMProtocol, ProtTomoBase):
    """
    Base class with methods used in the rest of the imod protocols
    """

    def __init__(self, **args):
        ProtTomoImportFiles.__init__(self, **args)

    def defineExecutionPararell(self):
        def __init__(self, **args):
            self.stepsExecutionMode = STEPS_PARALLEL

    def _defineImportParams(self, form):
        """ Method to define import params in protocol form """
        ProtTomoImportFiles._defineImportParams(self, form)

    # --------------------------- CACULUS functions ---------------------------
    def convertInputStep(self, tsObjId, generateAngleFile=True, imodInterpolation=True):
        if isinstance(self.inputSetOfTiltSeries, SetOfTiltSeries):
            ts = self.inputSetOfTiltSeries[tsObjId]
        elif isinstance(self.inputSetOfTiltSeries, Pointer):
            ts = self.inputSetOfTiltSeries.get()[tsObjId]

        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        firstItem = ts.getFirstItem()

        outputTsFileName = os.path.join(tmpPrefix, firstItem.parseFileName())

        """Apply the transformation form the input tilt-series"""
        # Use IMOD newstack interpolation
        if imodInterpolation:
            if firstItem.hasTransform():
                # Generate transformation matrices file
                outputTmFileName = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".xf"))
                utils.formatTransformFile(ts, outputTmFileName)

                # Apply interpolation
                paramsAlignment = {
                    'input': firstItem.getFileName(),
                    'output': outputTsFileName,
                    'xform': os.path.join(tmpPrefix, firstItem.parseFileName(extension=".xf")),
                }

                argsAlignment = "-input %(input)s " \
                                "-output %(output)s " \
                                "-xform %(xform)s " \

                rotationAngleAvg = utils.calculateRotationAngleFromTM(ts)

                # Check if rotation angle is greater than 45º. If so, swap x and y dimensions to adapt output image sizes to
                # the final sample disposition.
                if rotationAngleAvg > 45 or rotationAngleAvg < -45:
                    paramsAlignment.update({
                        'size': "%d,%d" % (firstItem.getYDim(), firstItem.getXDim())
                    })

                    argsAlignment += "-size %(size)s "

                Plugin.runImod(self, 'newstack', argsAlignment % paramsAlignment)

            else:
                path.createLink(firstItem.getLocation()[1], outputTsFileName)

        # Use Xmipp interpolation via Scipion
        else:
            outputTsFileName = os.path.join(tmpPrefix, firstItem.parseFileName())
            ts.applyTransform(outputTsFileName)

        if generateAngleFile:
            """Generate angle file"""
            angleFilePath = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt"))
            ts.generateTltFile(angleFilePath)

        # if generateExtraLink:
        #     """"Generate link to input tilt-series (no longer needed for fiducial viewers)"""
        #     inputTS = os.path.join(extraPrefix, firstItem.parseFileName())
        #     if firstItem.hasTransform():
        #         path.copyFile(outputTsFileName, inputTS)
        #
        #     else:
        #         path.createLink(firstItem.getLocation()[1], inputTS)

    # --------------------------- OUTPUT functions ----------------------------
    def getOutputSetOfTiltSeries(self, inputSet, binning=1):
        """ Method to generate output classes of set of tilt-series"""

        if hasattr(self, "outputSetOfTiltSeries"):
            self.outputSetOfTiltSeries.enableAppend()

        else:
            outputSetOfTiltSeries = self._createSetOfTiltSeries()

            if isinstance(inputSet, SetOfTiltSeries):
                outputSetOfTiltSeries.copyInfo(inputSet)
                outputSetOfTiltSeries.setDim(inputSet.getDim())

            elif isinstance(inputSet, SetOfTomograms):
                outputSetOfTiltSeries.setAcquisition(inputSet.getAcquisition())
                outputSetOfTiltSeries.setSamplingRate(inputSet.getSamplingRate())
                outputSetOfTiltSeries.setDim(inputSet.getDim())

            if binning > 1:
                samplingRate = inputSet.getSamplingRate()
                samplingRate *= self.binning.get()
                outputSetOfTiltSeries.setSamplingRate(samplingRate)

            outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfTiltSeries=outputSetOfTiltSeries)
            self._defineSourceRelation(inputSet, outputSetOfTiltSeries)

        return self.outputSetOfTiltSeries

    def getOutputInterpolatedSetOfTiltSeries(self, inputSet):
        """ Method to generate output interpolated classes of set of tilt-series"""

        if hasattr(self, "outputInterpolatedSetOfTiltSeries"):
            self.outputInterpolatedSetOfTiltSeries.enableAppend()

        else:
            outputInterpolatedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Interpolated')

            if isinstance(inputSet, SetOfTiltSeries):
                outputInterpolatedSetOfTiltSeries.copyInfo(inputSet)
                outputInterpolatedSetOfTiltSeries.setDim(inputSet.getDim())

            elif isinstance(inputSet, SetOfTomograms):
                outputInterpolatedSetOfTiltSeries.setAcquisition(inputSet.getAcquisition())
                outputInterpolatedSetOfTiltSeries.setSamplingRate(inputSet.getSamplingRate())
                outputInterpolatedSetOfTiltSeries.setDim(inputSet.getDim())

            if self.binning > 1:
                samplingRate = inputSet.getSamplingRate()
                samplingRate *= self.binning.get()
                outputInterpolatedSetOfTiltSeries.setSamplingRate(samplingRate)

            outputInterpolatedSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputInterpolatedSetOfTiltSeries=outputInterpolatedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputInterpolatedSetOfTiltSeries)

        return self.outputInterpolatedSetOfTiltSeries

    def getOutputFailedSetOfTiltSeries(self, inputSet):
        if hasattr(self, "outputFailedSetOfTiltSeries"):
            self.outputFailedSetOfTiltSeries.enableAppend()
        else:
            outputFailedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Failed')

            if isinstance(inputSet, SetOfTiltSeries):
                outputFailedSetOfTiltSeries.copyInfo(inputSet)
                outputFailedSetOfTiltSeries.setDim(inputSet.getDim())

            else:
                outputFailedSetOfTiltSeries.setAcquisition(inputSet.getAcquisition())
                outputFailedSetOfTiltSeries.setSamplingRate(inputSet.getSamplingRate())
                outputFailedSetOfTiltSeries.setDim(inputSet.getDim())

            outputFailedSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputFailedSetOfTiltSeries=outputFailedSetOfTiltSeries)
            self._defineSourceRelation(inputSet, outputFailedSetOfTiltSeries)

        return self.outputFailedSetOfTiltSeries

    def getOutputFiducialModelNoGaps(self):
        if hasattr(self, "outputFiducialModelNoGaps"):
            self.outputFiducialModelNoGaps.enableAppend()

        else:
            outputFiducialModelNoGaps = self._createSetOfLandmarkModels(suffix='NoGaps')

            outputFiducialModelNoGaps.copyInfo(self.inputSetOfTiltSeries.get())

            outputFiducialModelNoGaps.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputFiducialModelNoGaps=outputFiducialModelNoGaps)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputFiducialModelNoGaps)

        return self.outputFiducialModelNoGaps

    def getOutputFiducialModelGaps(self):
        if hasattr(self, "outputFiducialModelGaps"):
            self.outputFiducialModelGaps.enableAppend()
        else:
            outputFiducialModelGaps = self._createSetOfLandmarkModels(suffix='Gaps')

            outputFiducialModelGaps.copyInfo(self.inputSetOfTiltSeries.get())

            outputFiducialModelGaps.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputFiducialModelGaps=outputFiducialModelGaps)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputFiducialModelGaps)

        return self.outputFiducialModelGaps

    def getOutputSetOfCoordinates3Ds(self, inputSet=None, outputSet=None):
        if hasattr(self, "outputSetOfCoordinates3D"):
            self.outputSetOfCoordinates3D.enableAppend()

        else:
            outputSetOfCoordinates3D = self._createSetOfCoordinates3D(volSet=outputSet,
                                                                      suffix='Fiducials3D')

            outputSetOfCoordinates3D.setSamplingRate(outputSet.getSamplingRate())
            outputSetOfCoordinates3D.setPrecedents(outputSet)

            outputSetOfCoordinates3D.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfCoordinates3D=outputSetOfCoordinates3D)
            self._defineSourceRelation(inputSet, outputSetOfCoordinates3D)

        return self.outputSetOfCoordinates3D

    def getOutputSetOfTomograms(self, inputSet, binning=1):
        if hasattr(self, "outputSetOfTomograms"):
            self.outputSetOfTomograms.enableAppend()

        else:
            outputSetOfTomograms = self._createSetOfTomograms()

            if isinstance(inputSet, SetOfTomograms):
                outputSetOfTomograms.copyInfo(inputSet)

            elif isinstance(inputSet, SetOfTiltSeries):
                outputSetOfTomograms.setAcquisition(inputSet.getAcquisition())
                outputSetOfTomograms.setSamplingRate(inputSet.getSamplingRate())

            if binning > 1:
                samplingRate = inputSet.getSamplingRate()
                samplingRate *= self.binning.get()
                outputSetOfTomograms.setSamplingRate(samplingRate)

            outputSetOfTomograms.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfTomograms=outputSetOfTomograms)
            self._defineSourceRelation(inputSet, outputSetOfTomograms)

        return self.outputSetOfTomograms

    def getOutputSetOfCTFTomoSeries(self, outputSetName):
        if hasattr(self, outputSetName):
            outputSetOfCTFTomoSeries = getattr(self, outputSetName)
            if outputSetOfCTFTomoSeries is not None:
                outputSetOfCTFTomoSeries.enableAppend()
        else:
            outputSetOfCTFTomoSeries = SetOfCTFTomoSeries.create(self._getPath(),
                                                                 template='CTFmodels%s.sqlite')
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(self._getSetOfTiltSeries(pointer=True))
            outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputSetName: outputSetOfCTFTomoSeries})

    def addCTFTomoSeriesToSetFromDefocusFile(self, inputTs, defocusFilePath):
        """ This method generates a CtfTomoSeries Scipion object from a CTF estimation IMOD .defocus file.
        :param inputTs: tilt series associated to the CTF tomo series to be added.
        :param defocusFilePath: Location of the input .defocus file. """

        defocusFileFlag = utils.getDefocusFileFlag(defocusFilePath)

        tsId = inputTs.getTsId()
        tsObjId = inputTs.getObjId()

        newCTFTomoSeries = CTFTomoSeries()

        newCTFTomoSeries.copyInfo(inputTs)
        newCTFTomoSeries.setTiltSeries(inputTs)
        newCTFTomoSeries.setObjId(tsObjId)
        newCTFTomoSeries.setTsId(tsId)
        newCTFTomoSeries.setIMODDefocusFileFlag(defocusFileFlag)

        # We need to create now all the attributes of this object in order to append it to the set and be
        # able to update it posteriorly.

        newCTFTomoSeries.setNumberOfEstimationsInRange(None)

        output = getattr(self, self.outputSetName)
        output.append(newCTFTomoSeries)

        if defocusFileFlag == 0:
            " Plain estimation "
            defocusUDict = utils.readCTFEstimationInfoFile(defocusFilePath,
                                                           flag=defocusFileFlag)

        elif defocusFileFlag == 1:
            " Astigmatism estimation "
            defocusUDict, defocusVDict, defocusAngleDict = utils.readCTFEstimationInfoFile(defocusFilePath,
                                                                                           flag=defocusFileFlag)

        elif defocusFileFlag == 4:
            " Phase-shift information "
            defocusUDict, phaseShiftDict = utils.readCTFEstimationInfoFile(defocusFilePath,
                                                                           flag=defocusFileFlag)

        elif defocusFileFlag == 5:
            " Astigmatism and phase shift estimation "
            defocusUDict, defocusVDict, defocusAngleDict, phaseShiftDict = \
                utils.readCTFEstimationInfoFile(defocusFilePath,
                                                flag=defocusFileFlag)

        elif defocusFileFlag == 37:
            " Astigmatism, phase shift and cut-on frequency estimation "
            defocusUDict, defocusVDict, defocusAngleDict, phaseShiftDict, cutOnFreqDict = \
                utils.readCTFEstimationInfoFile(defocusFilePath,
                                                flag=defocusFileFlag)

        else:
            raise Exception(
                "Defocus file flag do not supported. Only supported formats corresponding to flags 0, "
                "1, 4, 5, and 37.")

        for index, _ in enumerate(inputTs):
            newCTFTomo = CTFTomo()
            newCTFTomo.setIndex(Integer(index + 1))

            if (index + 1) not in defocusUDict.keys():
                raise Exception("ERROR IN TILT-SERIES %s: NO CTF ESTIMATED FOR VIEW %d, TILT ANGLE %f" % (
                tsId, (index + 1), inputTs[index + 1].getTiltAngle()))

            if defocusFileFlag == 0:
                " Plain estimation "
                newCTFTomo._defocusUList = CsvList(pType=float)
                newCTFTomo.setDefocusUList(defocusUDict[index + 1])

            elif defocusFileFlag == 1:
                " Astigmatism estimation "
                newCTFTomo._defocusUList = CsvList(pType=float)
                newCTFTomo.setDefocusUList(defocusUDict[index + 1])

                newCTFTomo._defocusVList = CsvList(pType=float)
                newCTFTomo.setDefocusVList(defocusVDict[index + 1])

                newCTFTomo._defocusAngleList = CsvList(pType=float)
                newCTFTomo.setDefocusAngleList(defocusAngleDict[index + 1])

            elif defocusFileFlag == 4:
                " Phase-shift information "
                newCTFTomo._defocusUList = CsvList(pType=float)
                newCTFTomo.setDefocusUList(defocusUDict[index + 1])

                newCTFTomo._phaseShiftList = CsvList(pType=float)
                newCTFTomo.setPhaseShiftList(phaseShiftDict[index + 1])

            elif defocusFileFlag == 5:
                " Astigmatism and phase shift estimation "
                newCTFTomo._defocusUList = CsvList(pType=float)
                newCTFTomo.setDefocusUList(defocusUDict[index + 1])

                newCTFTomo._defocusVList = CsvList(pType=float)
                newCTFTomo.setDefocusVList(defocusVDict[index + 1])

                newCTFTomo._defocusAngleList = CsvList(pType=float)
                newCTFTomo.setDefocusAngleList(defocusAngleDict[index + 1])

                newCTFTomo._phaseShiftList = CsvList(pType=float)
                newCTFTomo.setPhaseShiftList(phaseShiftDict[index + 1])

            elif defocusFileFlag == 37:
                " Astigmatism, phase shift and cut-on frequency estimation "
                newCTFTomo._defocusUList = CsvList(pType=float)
                newCTFTomo.setDefocusUList(defocusUDict[index + 1])

                newCTFTomo._defocusVList = CsvList(pType=float)
                newCTFTomo.setDefocusVList(defocusVDict[index + 1])

                newCTFTomo._defocusAngleList = CsvList(pType=float)
                newCTFTomo.setDefocusAngleList(defocusAngleDict[index + 1])

                newCTFTomo._phaseShiftList = CsvList(pType=float)
                newCTFTomo.setPhaseShiftList(phaseShiftDict[index + 1])

                newCTFTomo._cutOnFreqList = CsvList(pType=float)
                newCTFTomo.setCutOnFreqList(cutOnFreqDict[index + 1])

                defocusUDict, defocusVDict, defocusAngleDict, phaseShiftDict, cutOnFreqDict = \
                    utils.readCTFEstimationInfoFile(defocusFilePath,
                                                    flag=defocusFileFlag)

            newCTFTomo.completeInfoFromList()

            newCTFTomoSeries.append(newCTFTomo)

        newCTFTomoSeries.setNumberOfEstimationsInRangeFromDefocusList()

        newCTFTomoSeries.calculateDefocusUDeviation(defocusUTolerance=self.defocusUTolerance)
        newCTFTomoSeries.calculateDefocusVDeviation(defocusVTolerance=self.defocusVTolerance)

        if not (newCTFTomoSeries.getIsDefocusUDeviationInRange() and
                newCTFTomoSeries.getIsDefocusVDeviationInRange()):
            newCTFTomoSeries.setEnabled(False)

        newCTFTomoSeries.write(properties=False)

        output.update(newCTFTomoSeries)
        output.write()

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def iterFiles(self):
        """ Iterate through the files matched with the pattern.
        Provide the fileName and fileId.
        """
        filePaths = self.getMatchFiles()

        filePaths = self._excludeByWords(filePaths)

        for fileName in filePaths:
            if self._idRegex:
                # Try to match the file id from filename
                # this is set by the user by using #### format in the pattern
                match = self._idRegex.match(fileName)
                if match is None:
                    raise Exception("File '%s' doesn't match the pattern '%s'"
                                    % (fileName, self.getPattern()))

                fileId = int(match.group(1))

            else:
                fileId = None

            yield fileName, fileId

    def _excludeByWords(self, files):
        exclusionWords = self.exclusionWords.get()

        if exclusionWords is None:
            return files

        exclusionWordList = exclusionWords.split()

        allowedFiles = []

        for file in files:
            if any(bannedWord in file for bannedWord in exclusionWordList):
                print("%s excluded. Contains any of %s" % (file, exclusionWords))
                continue
            allowedFiles.append(file)

        return allowedFiles
