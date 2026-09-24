# *****************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *              Scipion Team (scipion@cnb.csic.es) [1]
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
import sqlite3
import traceback
from os.path import exists
from imod.protocols.protocol_base import NEWSTACK_PROGRAM
from imod.protocols.protocol_base_preprocess import ProtImodBasePreprocess
from typing_extensions import List
from pwem import getExecStatusDir, appendStreamItem
from pwem.convert.headers import setMRCSamplingRate
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.utils import Message, cyanStr, redStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import SetOfTiltSeries, TiltImage, TiltSeries
from imod.constants import OUTPUT_TILTSERIES_NAME, ODD, EVEN
from tomo.protocols.protocol_base_streaming_tomo import ProtocolBaseStreamingTomo
from tomo.utils import writeTsSidecar

logger = logging.getLogger(__name__)


class ProtImodTsNormalization(ProtImodBasePreprocess, ProtocolBaseStreamingTomo):
    """
    Normalize input tilt-series and change its storing formatting.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/newstack.html

    IMOD tilt series preprocess makes use of the Newstack command.
    In particular, three functionalities are possible:\n

    _1 Binning_: The protocol also allows to bin tilt series. This
    means to reduce the dimensions of the tilt series keeping but
    keeping most of the information. The binning factor or simply
    binning is an integer number and represent the scaling factor
    of the images. Binning 2 means that the original images will
    be twice the binned ones.
    _2 Normalization_: This protocol allows to scale the gray values
    of the images, also called normalization, to a common range or
    mean of density. The most used normalization consists in zero
    mean and standard deviation one.\n

    _3 storage format_: IMOD is able to modify the number of bit of
    the stored data in order to reduce the disc occupancy.

    """

    _label = 'Tilt-series preprocess'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @classmethod
    def worksInStreaming(cls):
        return True

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form, *args):
        form.addSection(Message.LABEL_INPUT)
        super().addInTsSetFormParam(form)
        super()._defineParams(form)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self) -> None:
        inTsSet = self.getInputTsSet()
        if inTsSet.isStreamOpen():
            self._insertFunctionStep(self.stepsGeneratorStep,
                                     prerequisites=[],
                                     needsGPU=False)
        else:
            self._insertNonStreamingSteps()

    def _insertNonStreamingSteps(self):
        closeSetStepDeps = []
        self._initialize()
        inTsSet = self.getInputTsSet()
        tsList = [ts.clone() for ts in inTsSet.iterItems()]
        for ts in tsList:
            self._insertCommonSteps(ts, closeSetStepDeps=closeSetStepDeps)
        self._insertFunctionStep(self.closeOutputSetsStep,
                                 OUTPUT_TILTSERIES_NAME,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    def _insertCommonSteps(self, *stepsInputs, closeSetStepDeps: List[int]) -> None:
        ts = stepsInputs[0]
        convId = self._insertFunctionStep(self.linkTsStep,
                                          ts,
                                          prerequisites=[],
                                          needsGPU=False)
        compId = self._insertFunctionStep(self.generateOutputStackStep,
                                          ts,
                                          prerequisites=convId,
                                          needsGPU=False)
        outId = self._insertFunctionStep(self.createOutputStep,
                                         ts,
                                         prerequisites=compId,
                                         needsGPU=False)
        closeSetStepDeps.append(outId)

    def generateOutputStackStep(self, ts: TiltSeries):
        tsId = ts.getTsId()
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'===> tsId = {tsId}: preprocessing...'))
                norm = self.floatDensities.get()
                paramsDict = self.getBasicNewstackParams(ts,
                                                         self.getTmpOutFile(tsId),
                                                         self.getExtraOutFile(tsId),
                                                         binning=self.binning.get(),
                                                         doNorm=norm != 0)

                paramsDict["-antialias"] = self.antialias.get() + 1
                # Float densities
                if norm > 0:
                    paramsDict["-FloatDensities"] = norm
                    if norm == 2:
                        paramsDict["-MeanAndStandardDeviation"] = f"{self.scaleMean.get()},{self.scaleSd.get()}"
                    elif norm == 4:
                        paramsDict["-ScaleMinAndMax"] = f"{self.scaleMax.get()},{self.scaleMin.get()}"

                if self.getModeToOutput() is not None:
                    paramsDict["-ModeToOutput"] = self.getModeToOutput()

                self.runProgram(NEWSTACK_PROGRAM, paramsDict)

                if self.doOddEven:
                    paramsDict['-input'] = self.getTmpOutFile(tsId, suffix=ODD)
                    paramsDict['-output'] = self.getExtraOutFile(tsId, suffix=ODD)
                    self.runProgram(NEWSTACK_PROGRAM, paramsDict)

                    paramsDict['-input'] = self.getTmpOutFile(tsId, suffix=EVEN)
                    paramsDict['-output'] = self.getExtraOutFile(tsId, suffix=EVEN)
                    self.runProgram(NEWSTACK_PROGRAM, paramsDict)

            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> {NEWSTACK_PROGRAM} execution '
                                    f'failed with the exception -> {e}'))
                logger.error(traceback.format_exc())

    def createOutputStep(self, inTs: TiltSeries):
        tsId = inTs.getTsId()
        if tsId in self.failedItems:
            self.addToOutFailedSet(inTs)
            return

        try:
            outputFn = self.getExtraOutFile(tsId)
            if not exists(outputFn):
                logger.error(redStr(f'tsId = {tsId} -> Output file {outputFn} was not generated. Skipping... '))
                return

            binning = self.binning.get()
            samplingRate = self.getInputTsSet().getSamplingRate()
            if binning > 1:
                samplingRate *= binning


            setMRCSamplingRate(outputFn, samplingRate)  # Update the apix value in file header
            newTs = TiltSeries()
            newTs.copyInfo(inTs)

            inTiltList = inTs.loadTiltImgsInMemory()
            inTiltList.sort(key=lambda item: item.getIndex())
            tiltImages = []
            for inTi in inTiltList:
                newTi = TiltImage()
                newTi.copyInfo(inTi)
                newTi.setFileName(outputFn)
                self.updateTransformMatrix(newTi, binning=binning)
                self.setTsOddEven(tsId, newTi, binGenerated=True)
                tiltImages.append(newTi)
            self._registerOutput(newTs, tiltImages)

            # Streaming sidecar files
            execStatusDir = getExecStatusDir(self)
            if exists(execStatusDir):
                writeTsSidecar(execStatusDir, newTs, tiltImages)
                appendStreamItem(self, tsId)

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    @retry_on_sqlite_lock(log=logger, max_attempts=30, initial_delay=0.5,
                          backoff_factor=1.5, max_delay=15)
    def _registerOutput(self,
                        newTs: TiltSeries,
                        tiltImages: List[TiltImage]):
        with self._lock:
            try:
                binning = self.binning.get()
                # Set of tilt-series
                outTsSet = self.getOutputSetOfTS(self.getInputTsSet(pointer=True), binning)
                # Tilt-series
                outTsSet.append(newTs)
                # Tilt-images
                for newTi in tiltImages:
                    newTs.append(newTi)
                # Data persistence
                newTs.write()
                outTsSet.update(newTs)
                outTsSet.write()
                self._store(outTsSet)

            except sqlite3.OperationalError as e:
                # Release the write lock and reset the in-memory append state so
                # the @retry_on_sqlite_lock retry is a clean, non-hogging redo
                # (covers the later commits -- newTs.write/outTsSet.write -- not
                # just the append phase) and never trips the duplicate-tsId guard.
                self._releaseOutputWriteLock(outTsSet, newTs.getTsId())
                raise e

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            summary.append(f"Input tilt-series: {self.getInputTsSet().getSize()}\n"
                           f"Interpolations applied: {output.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            methods.append(f"{output.getSize()} tilt-series have been "
                           "normalized using the IMOD *newstack* command.")
        return methods

    # --------------------------- UTILS functions -----------------------------
    def getModeToOutput(self):
        parseParamsOutputMode = {
            0: None,
            1: 101,
            2: 0,
            3: 1,
            4: 6,
            5: 2
        }
        return parseParamsOutputMode[self.modeToOutput.get()]

    @staticmethod
    def updateTransformMatrix(ti: TiltImage, binning: int = 1) -> None:
        if ti.hasTransform() and binning != 1:
            transform = ti.getTransform()
            matrix = transform.getMatrix()

            matrix[0][2] /= binning
            matrix[1][2] /= binning

            transform.setMatrix(matrix)
            ti.setTransform(transform)


