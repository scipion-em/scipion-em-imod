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
import logging

import numpy as np

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from imod.protocols.protocol_base import IN_TS_SET
from pyworkflow.protocol.constants import STEPS_SERIAL
import pwem.objects as data
from tomo.objects import SetOfTiltSeries
from tomo.protocols.protocol_base import ProtTomoImportFiles
from tomo.convert.mdoc import normalizeTSId

from imod import utils
from imod.constants import XF_EXT, OUTPUT_TILTSERIES_NAME
from imod.protocols import ProtImodBase

logger = logging.getLogger(__name__)


class ProtImodImportTransformationMatrix(ProtImodBase, ProtTomoImportFiles):
    """
    Import the transformation matrices assigned to an input set of tilt-series
    """
    _label = 'Import transformation matrix'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}

    def __init__(self, **kwargs):
        ProtImodBase().__init__(**kwargs)
        self.stepsExecutionMode = STEPS_SERIAL
        ProtTomoImportFiles.__init__(self, **kwargs)
        self.matchingTsIds = None
        self.iterFilesDict = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        ProtTomoImportFiles._defineImportParams(self, form)
        ProtTomoImportFiles.addExclusionWordsParam(form)
        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Set of tilt-series on which transformation matrices '
                           'will be assigned.',
                      label='Input set of tilt-series')

        groupMatchBinning = form.addGroup('Match binning')

        groupMatchBinning.addParam('binningTM',
                                   params.IntParam,
                                   default=1,
                                   label='Transformation matrix binning',
                                   help='Binning of the tilt series at which '
                                        'the transformation matrices were '
                                        'calculated.')

        groupMatchBinning.addParam('binningTS',
                                   params.IntParam,
                                   default=1,
                                   label='Tilt-series binning',
                                   help='Binning of the tilt-series.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        matchBinningFactor = self.binningTM.get() / self.binningTS.get()
        for tsId in self.tsDict.keys():
            self._insertFunctionStep(self.generateTransformFileStep, tsId, matchBinningFactor)
            self._insertFunctionStep(self.assignTransformationMatricesStep, tsId)

        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.initializeParsing()
        if self.regEx:
            logger.info("Using regex pattern: '%s'" % self.regExPattern)
            logger.info("Generated glob pattern: '%s'" % self.globPattern)
            self.iterFilesDict = self.getMatchingFilesFromRegEx()
        else:
            dictBaseNames = {}
            for iFile in self.iterFiles():
                # We will Look for basename - tsId or base -name normalized basename - tsId matches. See tomo.convert.mdoc
                # normalizeTSId
                iFname = iFile[0]
                fBaseName = pwutils.removeBaseExt(iFname)
                dictBaseNames[fBaseName] = iFname
                dictBaseNames[normalizeTSId(fBaseName)] = iFname
            self.iterFilesDict = dictBaseNames
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self.getInputSet() if ts.getTsId()
                       in self.iterFilesDict.keys()}  # Use only the ones that are not excluded with the excluded words

    def generateTransformFileStep(self, tsId, matchBinningFactor):
        self.genTsPaths(tsId)
        ts = self.tsDict[tsId]
        tiNum = ts.getSize()
        outputTransformFile = self.getExtraOutFile(tsId, ext=XF_EXT)
        self.debug(f"Matching tsIds: {self.matchingTsIds}")

        tmFilePath = self.iterFilesDict.get(tsId, None)
        if tmFilePath:
            if matchBinningFactor != 1:
                inputTransformMatrixList = utils.formatTransformationMatrix(tmFilePath)
                # Update shifts from the transformation matrix considering
                # the matching binning between the input tilt
                # series and the transformation matrix. We create an empty
                # tilt-series containing only tilt-images
                # with transform information.
                transformMatrixList = []

                for index in range(tiNum):
                    inputTransformMatrix = inputTransformMatrixList[:, :, index]

                    outputTransformMatrix = inputTransformMatrix
                    outputTransformMatrix[0][0] = inputTransformMatrix[0][0]
                    outputTransformMatrix[0][1] = inputTransformMatrix[0][1]
                    outputTransformMatrix[0][2] = inputTransformMatrix[0][2] * matchBinningFactor
                    outputTransformMatrix[1][0] = inputTransformMatrix[1][0]
                    outputTransformMatrix[1][1] = inputTransformMatrix[1][1]
                    outputTransformMatrix[1][2] = inputTransformMatrix[1][2] * matchBinningFactor
                    outputTransformMatrix[2][0] = inputTransformMatrix[2][0]
                    outputTransformMatrix[2][1] = inputTransformMatrix[2][1]
                    outputTransformMatrix[2][2] = inputTransformMatrix[2][2]

                    transformMatrixList.append(outputTransformMatrix)

                utils.formatTransformFileFromTransformList(transformMatrixList, outputTransformFile)

            else:
                pwutils.createLink(tmFilePath, outputTransformFile)

    def assignTransformationMatricesStep(self, tsId):
        ts = self.tsDict[tsId]
        outputTransformFile = self.getExtraOutFile(tsId, ext=XF_EXT)
        output = self.getOutputSetOfTS(self.getInputSet(pointer=True))
        alignmentMatrix = utils.formatTransformationMatrix(outputTransformFile)

        self.copyTsItems(output, ts, tsId,
                         updateTsCallback=self.updateTs,
                         updateTiCallback=self.updateTi,
                         copyId=True,
                         copyTM=False,
                         alignmentMatrix=alignmentMatrix,
                         isSemiStreamified=False)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errorMsg = []
        self.initializeParsing()
        if self.regEx:
            matchingFileDict = self.getMatchingFilesFromRegEx()
            if not matchingFileDict:
                errorMsg.append('No files matching the pattern %s were found.' % self.globPattern)
        else:
            matchingFiles = self.getMatchFiles()
            if matchingFiles:
                tsIdList = self.getInputSet().getTSIds()
                tmFileList = [normalizeTSId(fn) for fn, _ in self.iterFiles()]
                self.matchingTsIds = list(set(tsIdList) & set(tmFileList))
                if not self.matchingTsIds:
                    errorMsg.append("No matching files found.\n\n"
                                    f"\tThe tsIds detected are: {tsIdList}\n"
                                    "\tThe transform files base names detected are: "
                                    f"{tmFileList}")
            else:
                errorMsg.append("Unable to find the files provided:\n\n"
                                f"\t-filePath = {self.filesPath.get()}\n"
                                f"\t-pattern = {self.filesPattern.get()}")

        return errorMsg

    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append(f"Input tilt-series: {self.getInputSet().getSize()}\n"
                           "Transformation matrices assigned: "
                           f"{self.TiltSeries.getSize()}")
        return summary

    # --------------------------- UTILS functions -----------------------------
    def iterFiles(self):
        """ Iterate through the files matched with the pattern.
        Returns the fileName and fileId.
        """
        filePaths = self.getMatchFiles()
        filePaths = self._excludeByWords(filePaths)

        for fileName in filePaths:
            if self._idRegex:
                # Try to match the file id from filename
                # this is set by the user by using #### format in the pattern
                match = self._idRegex.match(fileName)
                if match is None:
                    raise ValueError(f"File {fileName} doesn't match the "
                                     f"pattern '{self.getPattern()}'")
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
                print(f"{file} excluded. Contains any of {exclusionWords}")
                continue
            allowedFiles.append(file)

        return allowedFiles

    def updateTi(self, origIndex, index, tsId, ts, ti, tsOut, tiOut, **kwargs):
        transform = data.Transform()
        alignmentMatrix = kwargs.get("alignmentMatrix")

        if ti.hasTransform():
            previousTransform = ti.getTransform().getMatrix()
            newTransform = alignmentMatrix[:, :, index]
            previousTransformArray = np.array(previousTransform)
            newTransformArray = np.array(newTransform)
            outputTransformMatrix = np.matmul(previousTransformArray, newTransformArray)
            transform.setMatrix(outputTransformMatrix)
        else:
            transform.setMatrix(alignmentMatrix[:, :, index])

        tiOut.setTransform(transform)

    @staticmethod
    def updateTs(tsId, ts, tsOut, **kwargs):
        tsOut.setAlignment2D()

