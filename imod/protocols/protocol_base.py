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

import os

from pyworkflow.object import Set, CsvList, Pointer
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.utils import path
from pwem.protocols import EMProtocol
from tomo.protocols.protocol_base import ProtTomoBase, ProtTomoImportFiles
from tomo.objects import (SetOfTiltSeries, SetOfTomograms, SetOfCTFTomoSeries,
                          CTFTomoSeries, CTFTomo, SetOfTiltSeriesCoordinates)

from .. import Plugin, utils

OUTPUT_TS_COORDINATES_NAME = "TiltSeriesCoordinates"
OUTPUT_FIDUCIAL_NO_GAPS_NAME = "FiducialModelNoGaps"
OUTPUT_FIDUCIAL_GAPS_NAME = "FiducialModelGaps"
OUTPUT_TILTSERIES_NAME = "TiltSeries"
OUTPUT_TS_INTERPOLATED_NAME = "InterpolatedTiltSeries"
OUTPUT_TS_FAILED_NAME = "FailedTiltSeries"
OUTPUT_CTF_SERIE = "CTFTomoSeries"
OUTPUT_TOMOGRAMS_NAME = "Tomograms"
OUTPUT_COORDINATES_3D_NAME = "Coordinates3D"


class ProtImodBase(ProtTomoImportFiles, EMProtocol, ProtTomoBase):
    """
    Base class with methods used in the rest of the imod protocols
    """

    def __init__(self, **args):

        # Possible outputs (synchronize these names with the constants)
        self.TiltSeriesCoordinates = None
        self.FiducialModelNoGaps = None
        self.FiducialModelGaps = None
        self.TiltSeries = None
        self.InterpolatedTiltSeries = None
        self.CTFTomoSeries = None
        self.FailedTiltSeries = None
        self.Tomograms = None
        self.Coordinates3D = None

        ProtTomoImportFiles.__init__(self, **args)

    @classmethod
    def worksInStreaming(cls):
        """ So far none of them work in streaming. Since this inherits from the import they were considered as "streamers". """
        return False

    def defineExecutionPararell(self):

        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineImportParams(self, form):
        """ Method to define import params in protocol form """
        ProtTomoImportFiles._defineImportParams(self, form)

    # --------------------------- CACULUS functions ---------------------------
    def convertInputStep(self, tsObjId, generateAngleFile=True,
                         imodInterpolation=True, doSwap=False):
        """

        :param tsObjId: Tilt series identifier
        :param generateAngleFile:  Boolean(True) to generate IMOD angle file
        :param imodInterpolation: Boolean (True) to interpolate the tilt series with
                                  imod in case there is a TM.
                                  Pass None to cancel interpolation.
        :param doSwap: if applying alignment, consider swapping X/y
        :return:
        """
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

        # .. Interpolation cancelled
        if imodInterpolation is None:
            self.info("Tilt series %s linked." % tsId)
            path.createLink(firstItem.getFileName(), outputTsFileName)

        elif imodInterpolation:
            """Apply the transformation form the input tilt-series"""
            # Use IMOD newstack interpolation
            if firstItem.hasTransform():
                # Generate transformation matrices file
                outputTmFileName = os.path.join(tmpPrefix,
                                                firstItem.parseFileName(extension=".xf"))
                utils.formatTransformFile(ts, outputTmFileName)

                argsAlignment, paramsAlignment = self.getBasicNewstackParams(ts,
                                                                             outputTsFileName,
                                                                             xfFile=outputTmFileName,
                                                                             firstItem=firstItem,
                                                                             doSwap=doSwap)

                self.info("Interpolating tilt series %s with imod" % tsId)
                Plugin.runImod(self, 'newstack', argsAlignment % paramsAlignment)

            else:
                self.info("Linking tilt series %s" % tsId)
                path.createLink(firstItem.getFileName(), outputTsFileName)

        # Use Xmipp interpolation via Scipion
        else:
            self.info("Interpolating tilt series %s with emlib" % tsId)
            ts.applyTransform(outputTsFileName)

        self.info("Tilt series %s available for processing at %s." % (tsId, outputTsFileName))

        if generateAngleFile:
            """Generate angle file"""
            angleFilePath = os.path.join(tmpPrefix,
                                         firstItem.parseFileName(extension=".tlt"))
            ts.generateTltFile(angleFilePath)

    def getBasicNewstackParams(self, ts, outputTsFileName, inputTsFileName=None,
                               xfFile=None, firstItem=None, binning=1, doSwap=False):
        """ Returns basic newstack arguments
        
        :param ts: Title Series object
        :param outputTsFileName: tilt series output file name after newstack
        :param inputTsFileName: Input tilt series file name. Default to firsItem.getFilename()
        :param xfFile: xf file name, if passed, alignment will be generated and used
        :param firstItem: Optional, otherwise it will be taken from ts
        :param binning: Default to 1. to apply to output size
        :param doSwap: Default False.
        
        """
        
        if firstItem is None:
            firstItem = ts.getFirstItem()

        if inputTsFileName is None:
            inputTsFileName = firstItem.getFileName()

        # Apply interpolation
        paramsAlignment = {
            'input': inputTsFileName,
            'output': outputTsFileName,
        }
        argsAlignment = "-input %(input)s " \
                        "-output %(output)s " \
                        "-taper 1,1 "

        if xfFile is not None:
            paramsAlignment['xform'] = xfFile
            argsAlignment += "-xform %(xform)s "

            if doSwap:
                rotationAngle = ts.getAcquisition().getTiltAxisAngle()
                # Check if rotation angle is greater than 45º. If so,
                # swap x and y dimensions to adapt output image sizes to
                # the final sample disposition.
                if 45 < abs(rotationAngle) < 135:
                    paramsAlignment.update({
                        'size': "%d,%d" % (round(firstItem.getYDim()/binning),
                                           round(firstItem.getXDim()/binning))
                    })

                    argsAlignment += "-size %(size)s "

        return argsAlignment, paramsAlignment

    # --------------------------- OUTPUT functions ----------------------------
    def getOutputSetOfTiltSeries(self, inputSet, binning=1):
        """ Method to generate output classes of set of tilt-series"""

        if self.TiltSeries:
            self.TiltSeries.enableAppend()

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

            self._defineOutputs(**{OUTPUT_TILTSERIES_NAME: outputSetOfTiltSeries})
            self._defineSourceRelation(inputSet, outputSetOfTiltSeries)

        return self.TiltSeries

    def getOutputInterpolatedSetOfTiltSeries(self, inputSet):
        """ Method to generate output interpolated classes of set of tilt-series"""

        if self.InterpolatedTiltSeries:
            self.InterpolatedTiltSeries.enableAppend()

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

            self._defineOutputs(**{OUTPUT_TS_INTERPOLATED_NAME: outputInterpolatedSetOfTiltSeries})
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputInterpolatedSetOfTiltSeries)

        return self.InterpolatedTiltSeries

    def getOutputFailedSetOfTiltSeries(self, inputSet):
        if self.FailedTiltSeries:
            self.FailedTiltSeries.enableAppend()
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

            self._defineOutputs(**{OUTPUT_TS_FAILED_NAME: outputFailedSetOfTiltSeries})
            self._defineSourceRelation(inputSet, outputFailedSetOfTiltSeries)

        return self.FailedTiltSeries

    def getOutputFiducialModelNoGaps(self, tiltSeries=None):

        if self.FiducialModelNoGaps:
            self.FiducialModelNoGaps.enableAppend()

        else:

            if tiltSeries is None:
                tiltSeriesPointer = self.inputSetOfTiltSeries
                tiltSeries = tiltSeriesPointer.get()
            else:
                tiltSeriesPointer = tiltSeries

            outputFiducialModelNoGaps = self._createSetOfLandmarkModels(suffix='NoGaps')

            outputFiducialModelNoGaps.copyInfo(tiltSeries)
            outputFiducialModelNoGaps.setSetOfTiltSeries(tiltSeriesPointer)

            outputFiducialModelNoGaps.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_FIDUCIAL_NO_GAPS_NAME: outputFiducialModelNoGaps})
            self._defineSourceRelation(tiltSeriesPointer, outputFiducialModelNoGaps)

        return self.FiducialModelNoGaps

    def getOutputFiducialModelGaps(self):

        if self.FiducialModelGaps:
            self.FiducialModelGaps.enableAppend()
        else:
            outputFiducialModelGaps = self._createSetOfLandmarkModels(suffix='Gaps')

            outputFiducialModelGaps.copyInfo(self.inputSetOfTiltSeries.get())
            outputFiducialModelGaps.setSetOfTiltSeries(self.inputSetOfTiltSeries)
            outputFiducialModelGaps.setStreamState(Set.STREAM_OPEN)


            self._defineOutputs(**{OUTPUT_FIDUCIAL_GAPS_NAME: outputFiducialModelGaps})
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputFiducialModelGaps)

        return self.FiducialModelGaps

    def getOutputSetOfTiltSeriesCoordinates(self, setOfTiltSeries=None):

        if self.TiltSeriesCoordinates:
            self.TiltSeriesCoordinates.enableAppend()

        else:
            outputSetOfCoordinates3D = SetOfTiltSeriesCoordinates.create(self._getPath(),
                                                                         suffix='Fiducials3D')

            outputSetOfCoordinates3D.setSetOfTiltSeries(setOfTiltSeries)
            outputSetOfCoordinates3D.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_TS_COORDINATES_NAME: outputSetOfCoordinates3D})
            self._defineSourceRelation(setOfTiltSeries, outputSetOfCoordinates3D)

        return self.TiltSeriesCoordinates

    def getOutputSetOfCoordinates3Ds(self, inputSet=None, outputSet=None):

        if self.Coordinates3D:
            self.Coordinates3D.enableAppend()

        else:
            outputSetOfCoordinates3D = self._createSetOfCoordinates3D(volSet=outputSet,
                                                                      suffix='Fiducials3D')

            outputSetOfCoordinates3D.setSamplingRate(outputSet.getSamplingRate())
            outputSetOfCoordinates3D.setPrecedents(outputSet)
            outputSetOfCoordinates3D.setBoxSize(32)

            outputSetOfCoordinates3D.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_COORDINATES_3D_NAME: outputSetOfCoordinates3D})
            self._defineSourceRelation(inputSet, outputSetOfCoordinates3D)

        return self.Coordinates3D

    def getOutputSetOfTomograms(self, inputSet, binning=1):

        if self.Tomograms:
            getattr(self, OUTPUT_TOMOGRAMS_NAME).enableAppend()

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

            self._defineOutputs(**{OUTPUT_TOMOGRAMS_NAME: outputSetOfTomograms})
            self._defineSourceRelation(inputSet, outputSetOfTomograms)

        return self.Tomograms

    def getOutputSetOfCTFTomoSeries(self, outputSetName):

        outputSetOfCTFTomoSeries = getattr(self, outputSetName, None)

        if outputSetOfCTFTomoSeries:
            outputSetOfCTFTomoSeries.enableAppend()
        else:
            outputSetOfCTFTomoSeries = SetOfCTFTomoSeries.create(self._getPath(),
                                                                 template='CTFmodels%s.sqlite')
            ts = self._getSetOfTiltSeries(pointer=True)
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(ts)
            outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputSetName: outputSetOfCTFTomoSeries})
            self._defineCtfRelation(outputSetOfCTFTomoSeries, ts.get())

        return outputSetOfCTFTomoSeries

    def addCTFTomoSeriesToSetFromDefocusFile(self, inputTs, defocusFilePath, output):
        """ This method generates a CtfTomoSeries Scipion object
        from a CTF estimation IMOD .defocus file.

        :param inputTs: tilt series associated to the CTF tomo series to be added.
        :param defocusFilePath: Location of the input .defocus file.
        :param output: SetOfCTFTomoSeries to do the append

        """

        defocusFileFlag = utils.getDefocusFileFlag(defocusFilePath)

        tsId = inputTs.getTsId()
        tsObjId = inputTs.getObjId()

        newCTFTomoSeries = CTFTomoSeries()

        newCTFTomoSeries.copyInfo(inputTs)
        newCTFTomoSeries.setTiltSeries(inputTs)
        newCTFTomoSeries.setObjId(tsObjId)
        newCTFTomoSeries.setTsId(tsId)
        newCTFTomoSeries.setIMODDefocusFileFlag(defocusFileFlag)

        # We need to create now all the attributes of this object
        # in order to append it to the set and be
        # able to update it later.

        newCTFTomoSeries.setNumberOfEstimationsInRange(None)

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
                f"Defocus file flag {defocusFileFlag} is not supported. Only supported formats "
                "correspond to flags 0, 1, 4, 5, and 37.")

        excludedViews = inputTs.getExcludedViewsIndex()
        ids = inputTs.getIdSet()
        for index in ids:
            newCTFTomo = CTFTomo()
            newCTFTomo.setIndex(index)

            if index not in defocusUDict.keys() and index not in excludedViews:
                raise Exception("ERROR IN TILT-SERIES %s: NO CTF ESTIMATED FOR VIEW %d, TILT ANGLE %f" % (
                    tsId, index, inputTs[index].getTiltAngle()))

            " Plain estimation (any defocus flag)"
            newCTFTomo._defocusUList = CsvList(pType=float)
            newCTFTomo.setDefocusUList(defocusUDict.get(index, [0.]))

            if defocusFileFlag == 1:
                " Astigmatism estimation "
                newCTFTomo._defocusVList = CsvList(pType=float)
                newCTFTomo.setDefocusVList(defocusVDict.get(index, [0.]))

                newCTFTomo._defocusAngleList = CsvList(pType=float)
                newCTFTomo.setDefocusAngleList(defocusAngleDict.get(index, [0.]))

            elif defocusFileFlag == 4:
                " Phase-shift information "
                newCTFTomo._phaseShiftList = CsvList(pType=float)
                newCTFTomo.setPhaseShiftList(phaseShiftDict.get(index, [0.]))

            elif defocusFileFlag == 5:
                " Astigmatism and phase shift estimation "
                newCTFTomo._defocusVList = CsvList(pType=float)
                newCTFTomo.setDefocusVList(defocusVDict.get(index, [0.]))

                newCTFTomo._defocusAngleList = CsvList(pType=float)
                newCTFTomo.setDefocusAngleList(defocusAngleDict.get(index, [0.]))

                newCTFTomo._phaseShiftList = CsvList(pType=float)
                newCTFTomo.setPhaseShiftList(phaseShiftDict.get(index, [0.]))

            elif defocusFileFlag == 37:
                " Astigmatism, phase shift and cut-on frequency estimation "
                newCTFTomo._defocusVList = CsvList(pType=float)
                newCTFTomo.setDefocusVList(defocusVDict.get(index, [0.]))

                newCTFTomo._defocusAngleList = CsvList(pType=float)
                newCTFTomo.setDefocusAngleList(defocusAngleDict.get(index, [0.]))

                newCTFTomo._phaseShiftList = CsvList(pType=float)
                newCTFTomo.setPhaseShiftList(phaseShiftDict.get(index, [0.]))

                newCTFTomo._cutOnFreqList = CsvList(pType=float)
                newCTFTomo.setCutOnFreqList(cutOnFreqDict.get(index, [0.]))

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

    # --------------------------- UTILS functions -----------------------------
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
