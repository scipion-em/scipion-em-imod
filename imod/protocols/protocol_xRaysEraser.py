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
import logging
import traceback
from collections import Counter
from os.path import exists
import pyworkflow.protocol.params as params
from pwem.convert.headers import setMRCSamplingRate
from pyworkflow.protocol import STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.utils import Message, cyanStr, redStr
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from imod.protocols import ProtImodBase
from imod.constants import OUTPUT_TILTSERIES_NAME, ODD, EVEN, MOD_EXT, CCDERASER_PROGRAM

logger = logging.getLogger(__name__)


class ProtImodXraysEraser(ProtImodBase, ProtStreamingBase):
    """
    Erase Xrays from aligned tilt-series based on the IMOD procedure.
    More info:
            https://bio3d.colorado.edu/imod/doc/man/ccderaser.html
    """

    _label = 'X-rays eraser'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @classmethod
    def worksInStreaming(cls):
        return True

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        super().addInTsSetFormParam(form)
        form.addParam('peakCriterion',
                      params.FloatParam,
                      default=8.0,
                      label='Peak criterion (in std)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Criterion # of SDs above local mean for erasing '
                           'peak based on intensity (the default is 10 SDs)')
        form.addParam('diffCriterion',
                      params.FloatParam,
                      default=6.0,
                      label='Difference criterion  (in std)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Criterion # of SDs above mean pixel-to-pixel '
                           'difference for erasing a peak based on '
                           'differences (the default is 8 SDs).')
        form.addParam('maximumRadius',
                      params.FloatParam,
                      default=4.2,
                      label='Maximum radius (px)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Maximum radius of peak area to erase (the '
                           'default is 4.2 pixels).')
        form.addParam('bigDiffCriterion',
                      params.IntParam,
                      default=19,
                      label='Big difference criterion  (in std)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='An extra-large peak will be erased only if the '
                           'value for the maximum difference between '
                           'adjacent pixels, averaged over the most extreme '
                           'one-fourth of the pixels in the patch, '
                           'exceeds this criterion, evaluated as the '
                           'number of SDs above the mean absolute difference '
                           'between adjacent pixels in the scan area. The '
                           'default is 19.  This high a value is needed '
                           'to prevent gold erasure on low-noise data sets '
                           'with small gold particles, and a lower value '
                           'may be needed to make extra-large peak removal '
                           'useful.')
        self.addOddEvenParams(form)
        form.addParallelSection(threads=3, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def stepsGeneratorStep(self) -> None:
        self._initialize()
        closeSetStepDeps = []
        inTsSet = self.getInputTsSet()
        outTsSet = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        self.readingOutput(outTsSet)

        while True:
            with self._lock:
                listInTsIds = inTsSet.getTSIds()
                tsToProcessDict = {ts.getTsId(): ts.clone() for ts in inTsSet.iterItems()
                                   if ts.getTsId() not in self.tsIdReadList  # Only not processed tsIds
                                   and ts.getSize() > 0}  # Avoid processing empty TS

            if not inTsSet.isStreamOpen() and Counter(self.tsIdReadList) == Counter(listInTsIds):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetsStep,
                                         OUTPUT_TILTSERIES_NAME,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break

            for tsId, ts in tsToProcessDict.items():
                cInId = self._insertFunctionStep(self.linkTsStep,
                                                 ts,
                                                 prerequisites=[],
                                                 needsGPU=False)
                compId = self._insertFunctionStep(self.eraseXraysStep,
                                                  ts,
                                                  prerequisites=cInId,
                                                  needsGPU=False)
                outId = self._insertFunctionStep(self.createOutputStep,
                                                 ts,
                                                 prerequisites=compId,
                                                 needsGPU=False)
                closeSetStepDeps.append(outId)
                logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                self.tsIdReadList.append(tsId)

            self.refreshStreaming(inTsSet)

    # -------------------------- STEPS functions ------------------------------
    def eraseXraysStep(self, ts: TiltSeries):
        tsId = ts.getTsId()
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'tsId = {tsId} -> Erasing the X-Rays...'))
                with self._lock:
                    firstItem = ts.getFirstEnabledItem()
                inputFile = self.getTmpOutFile(tsId)
                outputFile = self.getExtraOutFile(tsId)
                paramsCcderaser = self.getCcdEraserParamsDict(tsId, inputFile, outputFile)
                self.runProgram(CCDERASER_PROGRAM, paramsCcderaser)
                if self.doOddEven:
                    # Odd
                    logger.info(cyanStr(f'tsId = {tsId} -> Erasing the X-Rays (ODD Tilt-series) ...'))
                    inputFile = firstItem.getOdd()
                    outputFile = self.getExtraOutFile(tsId, suffix=ODD)
                    paramsCcderaser = self.getCcdEraserParamsDict(tsId, inputFile, outputFile)
                    self.runProgram(CCDERASER_PROGRAM, paramsCcderaser)
                    # Even
                    logger.info(cyanStr(f'tsId = {tsId} -> Erasing the X-Rays (EVEN Tilt-series) ...'))
                    inputFile = firstItem.getEven()
                    outputFile = self.getExtraOutFile(tsId, suffix=EVEN)
                    paramsCcderaser = self.getCcdEraserParamsDict(tsId, inputFile, outputFile)
                    self.runProgram(CCDERASER_PROGRAM, paramsCcderaser)
            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> {CCDERASER_PROGRAM} execution failed '
                                    f'with the exception -> {e}'))
                logger.error(traceback.format_exc())

    def createOutputStep(self, ts: TiltSeries):
        tsId = ts.getTsId()
        if tsId in self.failedItems:
            self.addToOutFailedSet(tsId)
            return
        try:
            outTsFile = self.getExtraOutFile(tsId)
            if not exists(outTsFile):
                logger.error(redStr(f'tsId = {tsId} -> Output file {outTsFile} was not generated. Skipping... '))

            setMRCSamplingRate(outTsFile, ts.getSamplingRate())  # Update the apix value in file header
            with self._lock:
                ts = self.getCurrentTs(tsId)
                # Set of tilt-series
                outTsSet = self.getOutputSetOfTS(self.getInputTsSet(pointer=True))
                # Tilt-series
                outTs = TiltSeries()
                outTs.copyInfo(ts)
                outTsSet.append(outTs)
                # Tilt-images
                for ti in ts.iterItems():
                    outTi = TiltImage()
                    outTi.copyInfo(ti)
                    outTi.setFileName(outTsFile)
                    self.setTsOddEven(tsId, outTi, binGenerated=True)
                    outTs.append(outTi)
                # Data persistence
                outTs.write()
                outTsSet.update(outTs)
                outTsSet.write()
                self._store(outTsSet)
                # Close explicitly the outputs (for streaming)
                self.closeOutputsForStreaming()

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    # --------------------------- UTILS functions -----------------------------
    def getCcdEraserParamsDict(self,
                               tsId: str,
                               inputFile: str,
                               outputFile=str) -> dict:
        return {
            "-InputFile": inputFile,
            "-OutputFile": outputFile,
            "-FindPeaks": 1,
            "-PeakCriterion": self.peakCriterion.get(),
            "-DiffCriterion": self.diffCriterion.get(),
            "-GrowCriterion": 4.,
            "-ScanCriterion": 3.,
            "-MaximumRadius": self.maximumRadius.get(),
            "-GiantCriterion": 12.,
            "-ExtraLargeRadius": 8.,
            "-BigDiffCriterion": self.bigDiffCriterion.get(),
            "-AnnulusWidth": 2.0,
            "-XYScanSize": 100,
            "-EdgeExclusionWidth": 4,
            "-PointModel": self.getExtraOutFile(tsId, suffix="fid", ext=MOD_EXT),
            "-BorderSize": 2,
            "-PolynomialOrder": 2,
        }

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            summary.append(f"Input tilt-series: {self.getInputTsSet().getSize()}\n"
                           "X-rays erased output tilt series: "
                           f"{output.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")

        return summary

    def _methods(self):
        methods = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output:
            methods.append(f"The x-rays artifacts have been erased for "
                           f"{output.getSize()} tilt-series using "
                           f"the IMOD *{CCDERASER_PROGRAM}* command.")

        return methods
