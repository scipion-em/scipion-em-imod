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
import sqlite3
import traceback
from os.path import exists
from typing import List

import numpy as np
import pyworkflow.protocol.params as params
from pwem import genExecStatusDir, getExecStatusDir, appendStreamItem
from pwem.convert.headers import setMRCSamplingRate
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.utils import Message, cyanStr, redStr, yellowStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from imod.protocols import ProtImodBase
from imod.constants import (ODD, EVEN, SCIPION_IMPORT, FIXED_DOSE,
                            OUTPUT_TILTSERIES_NAME, MTTFILTER_PROGRAM)
from tomo.protocols.protocol_base_streaming_tomo import ProtocolBaseStreamingTomo
from tomo.utils import sleepRandomly, writeTsSidecar

logger = logging.getLogger(__name__)


class ProtImodDoseFilter(ProtImodBase, ProtocolBaseStreamingTomo):
    """
    Tilt-series dose filtering based on the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/mtffilter.html

    A specialized filter can be applied to perform dose weight-filtering of
    cryoEM images, particularly ones from tilt series.  The filter is as
    described in Grant and Grigorieff, 2015 (DOI: 10.7554/eLife.06980) and
    the implementation follows that in their "unblur" program.  At any
    frequency, the filter follows an exponential decay with dose, where the
    exponential is of the dose divided by 2 times a "critical dose" for
    that frequency.  This critical dose was empirically found to be
    approximated by a * k^b + c, where k is frequency; the values of a, b, c in
    that paper are used by default.
    """

    _label = 'Dose filter'
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
        form.addParam('initialDose',
                      params.FloatParam,
                      default=0.0,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Initial dose (e/Å^2)',
                      help='Dose applied before any of the images in the '
                           'input file were taken; this value will be '
                           'added to all the dose values.')
        form.addParam('inputDoseType',
                      params.EnumParam,
                      choices=['Scipion import', 'Fixed dose'],
                      default=SCIPION_IMPORT,
                      label='Input dose source',
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='Where to find the dose information:\n'
                           '- Scipion import: use the dose provided '
                           'during import of the tilt-series\n'
                           '- Fixed dose: manually input fixed dose '
                           'for each image of the input file, '
                           'in electrons/Å^2.')
        form.addParam('fixedImageDose',
                      params.FloatParam,
                      default=FIXED_DOSE,
                      label='Fixed dose (e/Å^2)',
                      condition='inputDoseType == %i' % FIXED_DOSE,
                      help='Fixed dose for each image of the input file, '
                           'in electrons/square Ångstrom.')
        self.addOddEvenParams(form)
        form.addParallelSection(threads=3, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self) -> None:
        inTsSet = self.getInputTsSet()
        if inTsSet.isStreamOpen():
            self._insertFunctionStep(self.stepsGeneratorStep,
                                     prerequisites=[],
                                     needsGPU=False)
        else:
            self._insertNonStreamingSteps()

    # stepsGeneratorStep is centralized in ProtocolBaseStreamingTomo; the
    # per-protocol hooks it needs (_getStreamingInputTs, _getProcessedTsIds,
    # _getStreamingOutputNames, _streamingInitialize) are provided by ProtImodBase.

    def _insertNonStreamingSteps(self):
        closeSetStepDeps = []
        self._initialize()
        inTsSet = self.getInputTsSet()
        tsList = [ts.clone() for ts in inTsSet.iterItems()]
        for ts in tsList:
            self._insertCommonSteps(ts, closeSetStepDeps)
        self._insertFunctionStep(self._closeOutputSet,
                                 OUTPUT_TILTSERIES_NAME,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    def _insertCommonSteps(self, ts: TiltSeries, closeSetStepDeps: List[int]) -> None:
        cInId = self._insertFunctionStep(self.linkTsStep,
                                         ts,
                                         prerequisites=[],
                                         needsGPU=False)
        compId = self._insertFunctionStep(self.doseFilterStep,
                                          ts,
                                          prerequisites=cInId,
                                          needsGPU=False)
        outId = self._insertFunctionStep(self.createOutputStep,
                                         ts,
                                         prerequisites=[compId],
                                         needsGPU=False)
        closeSetStepDeps.append(outId)

    # --------------------------- STEPS functions -----------------------------
    def doseFilterStep(self, ts: TiltSeries):
        """Apply the dose filter to every tilt series"""
        tsId = ts.getTsId()
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'tsId = {tsId} -> Dose filtering...'))
                firstItem = ts.getFirstEnabledItem()

                progParams = {
                    '-input': self.getTmpOutFile(tsId),
                    '-output': self.getExtraOutFile(tsId),
                    '-PixelSize': ts.getSamplingRate(),
                    '-Voltage': int(ts.getAcquisition().getVoltage()),
                }

                if self.initialDose.get() != 0.0:
                    progParams["-InitialDose"] = self.initialDose.get()

                if self.inputDoseType.get() == SCIPION_IMPORT:
                    outputDoseFilePath = self.getExtraOutFile(tsId, ext="dose")
                    self.generateDoseFile(ts, outputDoseFilePath)
                    progParams["-TypeOfDoseFile"] = 2
                    progParams["-DoseWeightingFile"] = outputDoseFilePath

                elif self.inputDoseType.get() == FIXED_DOSE:
                    progParams["-FixedImageDose"] = self.fixedImageDose.get()

                self.runProgram(MTTFILTER_PROGRAM, progParams)

                if self.doOddEven:
                    # Odd
                    logger.info(cyanStr(f'tsId = {tsId} ODD -> Dose filtering...'))
                    progParams['-input'] = firstItem.getEven()
                    progParams['-output'] = self.getExtraOutFile(tsId, suffix=ODD)
                    self.runProgram(MTTFILTER_PROGRAM, progParams)
                    # Even
                    logger.info(cyanStr(f'tsId = {tsId} EVEN -> Dose filtering...'))
                    progParams['-input'] = firstItem.getOdd()
                    progParams['-output'] = self.getExtraOutFile(tsId, suffix=EVEN)
                    self.runProgram(MTTFILTER_PROGRAM, progParams)

            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> {MTTFILTER_PROGRAM} execution failed '
                                    f'with the exception -> {e}'))
                logger.error(traceback.format_exc())

    def createOutputStep(self, inTs: TiltSeries):
        """Generate output filtered tilt series"""
        tsId = inTs.getTsId()
        if tsId in self.failedItems:
            self.addToOutFailedSet(inTs)
            return

        try:
            outTsFile = self.getExtraOutFile(tsId)
            if not exists(outTsFile):
                logger.error(redStr(f'tsId = {tsId} -> Output file {outTsFile} was not generated. Skipping... '))
                return

            setMRCSamplingRate(outTsFile, inTs.getSamplingRate())  # Update the apix value in file header
            newTs = TiltSeries()
            newTs.copyInfo(inTs)
            self.updateTsAcquisition(newTs)  # Acquisition dose goes to 0 after having been applied

            inTiltList = inTs.loadTiltImgsInMemory()
            inTiltList.sort(key=lambda item: item.getIndex())
            tiltImages = []
            for inTi in inTiltList:
                outTi = TiltImage()
                outTi.copyInfo(inTi)
                outTi.setFileName(outTsFile)
                self.updateTiAcquisition(outTi)
                self.setTsOddEven(tsId, outTi, binGenerated=True)
                tiltImages.append(outTi)
            self._registerOutput(newTs, tiltImages)

            # Streaming only: publish the per-TS metadata sidecar (built from the
            # in-memory ts/tiltImages, no DB read) and the journal id
            execStatusDir = getExecStatusDir(self)
            if exists(execStatusDir):
                writeTsSidecar(execStatusDir, newTs, tiltImages)
                appendStreamItem(self, tsId)

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    # The producer's write competes with concurrent readers of the same
    # tiltseries.sqlite (journal_mode=DELETE => one writer vs many readers, e.g.
    # a chained downstream consumer). Use a more patient retry budget than the
    # default so a transient burst of consumer reads cannot exhaust it.
    @retry_on_sqlite_lock(log=logger, max_attempts=30, initial_delay=0.5,
                          backoff_factor=1.5, max_delay=15)
    def _registerOutput(self,
                        newTs: TiltSeries,
                        tiltImages: List[TiltImage]):
        with self._lock:
            # Set of tilt-series
            outTsSet = self.getOutputSetOfTS(self.getInputTsSet(pointer=True))
            try:
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
                outTsSet.rollbackFailedAppend(newTs.getTsId())
                raise e

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        if self.inputDoseType.get() == SCIPION_IMPORT:
            for ts in self.getInputTsSet():
                if ts.getFirstEnabledItem().getAcquisition().getDosePerFrame() is None:
                    validateMsgs.append(f"{ts.getTsId()} has no dose information stored "
                                        "in Scipion Metadata. To solve this, re-import "
                                        "tilt-series using the mdoc option.")
                    break

        return validateMsgs

    def _summary(self):
        summary = []

        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            summary.append(f"Input tilt-series: {self.getInputTsSet().getSize()}\n"
                           "Dose weighting applied: "
                           f"{output.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")

        return summary

    def _methods(self):
        methods = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            methods.append("The dose-weighting has been applied to "
                           f"{output.getSize()} "
                           "tilt-series using the IMOD *mtffilter* command.")
        return methods

    # --------------------------- UTILS functions -----------------------------
    @staticmethod
    def updateTiAcquisition(tiOut: TiltImage) -> None:
        """Sets the initial and accumulated doses to 0 for a given tilt-image"""
        # Output is dose-weighted
        acq = tiOut.getAcquisition()
        acq.setDoseInitial(0.)
        acq.setAccumDose(0.)
        tiOut.setAcquisition(acq)

    @staticmethod
    def updateTsAcquisition(tsOut: TiltSeries) -> None:
        """Sets the initial and accumulated doses to 0 for a given tilt-series"""
        # Output is dose-weighted
        acq = tsOut.getAcquisition()
        acq.setAccumDose(0.)
        acq.setDoseInitial(0.)
        tsOut.setAcquisition(acq)

    @staticmethod
    def generateDoseFile(ts: TiltSeries, doseFileOutputPath: str) -> None:
        """ This method generates a file containing the dose information
        of a tilt series in the specified location from the accumulated
        dose and dose per tilt. The format is two columns per each tilt image:
         the prior accumulated dose and the image dose
         """
        doseInfoList = []

        for ti in ts.iterItems(iterate=False):
            acq = ti.getAcquisition()
            doseInfoList.append((acq.getAccumDose() - acq.getDosePerFrame(), acq.getDosePerFrame()))

        np.savetxt(doseFileOutputPath, np.asarray(doseInfoList), fmt='%f', delimiter=" ")
