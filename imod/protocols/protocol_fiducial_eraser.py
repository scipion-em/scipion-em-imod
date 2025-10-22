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
from os.path import exists
import pyworkflow.protocol.params as params
from imod.protocols.protocol_base import IN_TS_SET
from pwem.convert.headers import setMRCSamplingRate
from pyworkflow.protocol import STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.utils import Message, cyanStr, redStr
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from imod.protocols import ProtImodBase
from imod.constants import OUTPUT_TILTSERIES_NAME, IMODFINDBEADS_PROGRAM, CCDERASER_PROGRAM, ODD, EVEN

logger = logging.getLogger(__name__)


class ProtImodFiducialEraser(ProtImodBase, ProtStreamingBase):
    """
    This protocol will erase the fiducial gold beads present in the tilt
     series images using IMOD procedures.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/imodfindbeads.html
        https://bio3d.colorado.edu/imod/doc/man/ccderaser.html
    """

    _label = 'fiducial eraser'
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
        form.addParam('fidDiameter',
                      params.FloatParam,
                      important=True,
                      default= 10.0,
                      label='Fiducial diameter (nm)')
        form.addParam('erasedDiameter',
                      params.FloatParam,
                      important=True,
                      default=12.0,
                      label='Diameter to erase (nm)')
        form.addParallelSection(threads=3, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def stepsGeneratorStep(self) -> None:
        self._initialize()
        closeSetStepDeps = []
        inTsSet = self.getInputTsSet()
        outTsSet = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        self.readingOutput(outTsSet)

        while True:
            listInTsIds = inTsSet.getTSIds()
            if not inTsSet.isStreamOpen() and self.tsIdReadList == listInTsIds:
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetsStep,
                                         OUTPUT_TILTSERIES_NAME,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break

            for ts in inTsSet.iterItems():
                tsId = ts.getTsId()
                if tsId not in self.tsIdReadList and ts.getSize() > 0:  # Avoid processing empty TS
                    cInId = self._insertFunctionStep(self.linkTsStep,
                                                     tsId,
                                                     prerequisites=[],
                                                     needsGPU=False)
                    beadsId = self._insertFunctionStep(self.imodfindbeadsStep,
                                                       tsId,
                                                       prerequisites=cInId,
                                                       needsGPU=False)
                    eraserId = self._insertFunctionStep(self.ccderaserStep,
                                                        tsId,
                                                        prerequisites=beadsId,
                                                        needsGPU=False)
                    createOutputId = self._insertFunctionStep(self.createOutputStep,
                                                              tsId,
                                                              prerequisites=eraserId,
                                                              needsGPU=False)
                    closeSetStepDeps.append(createOutputId)
                    logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                    self.tsIdReadList.append(tsId)

                self.refreshStreaming(inTsSet)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        super()._initialize()

    def imodfindbeadsStep(self, tsId: str):
        """This step creates a fiducial model"""
        try:
            logger.info(cyanStr(f'tsId = {tsId} -> Finding the fiducials...'))
            with self._lock:
                ts = self.getCurrentTs(tsId)
                firstTi = ts.getFirstItem()

            paramsImodFindBeads = {
                "-inp": firstTi.getFileName(),
                "-o": self._getModelFileName(tsId),
                "-size": self._getFiducialDiameterPx(ts)
            }
            self.runProgram(IMODFINDBEADS_PROGRAM, paramsImodFindBeads)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {IMODFINDBEADS_PROGRAM} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def ccderaserStep(self, tsId: str):
        """This step erase the gold beads from the fiducial model"""
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'tsId = {tsId} -> Erasing the gold beads...'))
                with self._lock:
                    ts = self.getCurrentTs(tsId)
                    firstTi = ts.getFirstItem()

                paramsCCDeraser = {
                    "-ExpandCircleIterations": 3,
                    "-BetterRadius": self._getFiducialDiameterPx(ts) / 2,
                    "-input": firstTi.getFileName(),
                    "-output": self.getExtraOutFile(tsId),
                    "-ModelFile": self._getModelFileName(tsId),
                    "-MergePatches": '',
                    "-ExcludeAdjacent": '',
                    "-CircleObjects": '/',
                    "-SkipTurnedOffPoints": 1,
                    "-PolynomialOrder": -1,
                }

                self.runProgram(CCDERASER_PROGRAM, paramsCCDeraser)

                if self.doOddEven:
                    logger.info(cyanStr(f'tsId = {tsId} -> Erasing the gold beads (ODD Tilt-series) ...'))
                    paramsCCDeraser['-input'] = ts.getOddFileName(),
                    paramsCCDeraser['-output'] = self.getExtraOutFile(tsId, suffix=ODD)
                    self.runProgram(CCDERASER_PROGRAM, paramsCCDeraser)

                    logger.info(cyanStr(f'tsId = {tsId} -> Erasing the gold beads (EVEN Tilt-series) ...'))
                    paramsCCDeraser['-input'] = ts.getEvenFileName(),
                    paramsCCDeraser['-output'] = self.getExtraOutFile(tsId, suffix=EVEN)
                    self.runProgram(CCDERASER_PROGRAM, paramsCCDeraser)

            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> {CCDERASER_PROGRAM} execution failed'
                                    f' with the exception -> {e}'))
                logger.error(traceback.format_exc())

    def createOutputStep(self, tsId: str):
        if tsId in self.failedItems:
            self.addToOutFailedSet(tsId)
            return

        try:
            outTsFile = self.getExtraOutFile(tsId)
            if exists(outTsFile):
                with self._lock:
                    ts = self.getCurrentTs(tsId)
                    # Set of tilt-series
                    outTsSet = self.getOutputSetOfTS(self.getInputTsSet(pointer=True))
                    # Tilt-series
                    outTs = TiltSeries()
                    outTs.copyInfo(ts)
                    outTsSet.append(outTs)
                    setMRCSamplingRate(outTsFile, ts.getSamplingRate())  # Update the apix value in file header
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
            else:
                logger.error(redStr(f'tsId = {tsId} -> Output file {outTsFile} was not generated. Skipping... '))

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    # --------------------------- INFO functions ------------------------------

    # --------------------------- UTILS functions -----------------------------
    def _getFiducialDiameterPx(self, ts: TiltSeries) -> float:
        return 10 * self.fidDiameter.get() / ts.getSamplingRate()

    def _getModelFileName(self, tsId: str) -> str:
        return self._getExtraPath(tsId, f'{tsId}_model')