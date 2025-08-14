# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
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
import csv
import logging
import numpy as np
import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from imod.convert.convert import readXfFile
from imod.protocols.protocol_base import IN_TS_SET
from pwem.objects import Transform
from pyworkflow.protocol.constants import STEPS_SERIAL
from pyworkflow.utils import cyanStr, redStr
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from tomo.protocols.protocol_base import ProtTomoImportFiles
from tomo.convert.mdoc import normalizeTSId
from imod.constants import XF_EXT, OUTPUT_TILTSERIES_NAME
from imod.protocols import ProtImodBase

logger = logging.getLogger(__name__)


class ProtImodImportTransformationMatrix(ProtImodBase, ProtTomoImportFiles):
    """
    Import the transformation matrices assigned to an input set of tilt-series
    """
    _label = 'Import transformation matrix'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}
    stepsExecutionMode = STEPS_SERIAL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.matchingTsIds = None
        self.iterFilesDict = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        ProtTomoImportFiles._defineImportParams(self, form)
        ProtTomoImportFiles.addExclusionWordsParam(form)
        super().addInTsSetFormParam(form)
        form.addParam("override",
                      params.BooleanParam,
                      default=True,
                      help='If True, the imported transformations will override the previous alignments, '
                           'otherwise, the alignments will be combined (alignment matrices multiplied.',
                      label='Override alignments')

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
            self._insertFunctionStep(self.generateTransformFileStep,
                                     tsId,
                                     matchBinningFactor,
                                     needsGPU=False)
            self._insertFunctionStep(self.assignTransformationMatricesStep,
                                     tsId,
                                     needsGPU=False)

        self._insertFunctionStep(self.closeOutputSetsStep,
                                 OUTPUT_TILTSERIES_NAME,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.initializeParsing()
        if self.regEx:
            logger.info(cyanStr("Using regex pattern: '%s'" % self.regExPattern))
            logger.info(cyanStr("Generated glob pattern: '%s'" % self.globPattern))
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
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self.getInputTsSet() if ts.getTsId()
                       in self.iterFilesDict.keys()}  # Use only the ones that are not excluded with the excluded words
        logger.info(cyanStr(f"Matching tsIds: {self.matchingTsIds}"))

    def generateTransformFileStep(self,
                                  tsId: str,
                                  matchBinningFactor: int):
        try:
            self.genTsPaths(tsId)
            ts = self.tsDict[tsId]
            tiNum = ts.getSize()
            outputTransformFile = self.getExtraOutFile(tsId, ext=XF_EXT)
            tmFilePath = self.iterFilesDict.get(tsId, None)
            if tmFilePath:
                if matchBinningFactor != 1:
                    inputTransformMatrixList = readXfFile(tmFilePath)
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

                    self.writeXfFileFromTrList(transformMatrixList, outputTransformFile)

                else:
                    pwutils.createLink(tmFilePath, outputTransformFile)
        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> failed with the exception -> {e}'))

    def assignTransformationMatricesStep(self, tsId: str):
        if tsId in self.failedItems:
            self.addToOutFailedSet(tsId)
        else:
            try:
                ts = self.tsDict[tsId]
                outputTransformFile = self.getExtraOutFile(tsId, ext=XF_EXT)
                outTsSet = self.getOutputSetOfTS(self.getInputTsSet(pointer=True))
                alignmentMatrix = readXfFile(outputTransformFile)
                outTs = TiltSeries()
                outTs.copyInfo(ts)
                outTs.setAlignment2D()
                outTsSet.append(outTs)
                enabledInd = 0
                for ti in ts.iterItems():
                    outTi = TiltImage()
                    outTi.copyInfo(ti)
                    if ti.isEnabled():
                        newTransform = alignmentMatrix[:, :, enabledInd]
                        self.updateTiltImage(ti, outTi, newTransform)
                        enabledInd += 1
                    outTs.append(outTi)
                outTsSet.update(outTs)

            except Exception as e:
                logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))

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
                tsIdList = self.getInputTsSet().getTSIds()
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
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            summary.append(f"Input tilt-series: {self.getInputTsSet().getSize()}\n"
                           "Transformation matrices assigned: "
                           f"{output.getSize()}")
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

    @staticmethod
    def writeXfFileFromTrList(transformMatrixList, transformFilePath):
        """ This method takes a list of Transform matrices and the output
        transformation file path and creates an
        IMOD-based transform file in the location indicated.
        :param transformMatrixList: list of transformation matrices to be written
        :param transformFilePath: output location of the file. """

        tsMatrixTransformList = []

        for tm in transformMatrixList:
            tmFlat = tm.flatten()
            transformIMOD = [tmFlat[0],
                             tmFlat[1],
                             tmFlat[3],
                             tmFlat[4],
                             tmFlat[2],
                             tmFlat[5]]
            tsMatrixTransformList.append(transformIMOD)

        with open(transformFilePath, 'w') as f:
            csvW = csv.writer(f, delimiter='\t')
            csvW.writerows(tsMatrixTransformList)

    def updateTiltImage(self,
                        ti: TiltImage,
                        tiOut: TiltImage,
                        newTransform: np.array) -> None:
        transform = Transform()
        if ti.hasTransform() and not self.override.get():
            previousTransform = ti.getTransform().getMatrix()
            previousTransformArray = np.array(previousTransform)
            newTransformArray = np.array(newTransform)
            outputTransformMatrix = np.matmul(previousTransformArray, newTransformArray)
            transform.setMatrix(outputTransformMatrix)
        else:
            transform.setMatrix(newTransform)
        tiOut.setTransform(transform)

