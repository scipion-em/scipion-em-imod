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
from os.path import exists
import numpy as np
import pyworkflow.protocol.params as params
from imod.protocols.protocol_base import IN_TS_SET
from pyworkflow.protocol import STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.utils import Message, cyanStr, redStr
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from imod.protocols import ProtImodBase
from imod.constants import (ODD, EVEN, SCIPION_IMPORT, FIXED_DOSE,
                            OUTPUT_TILTSERIES_NAME, MTTFILTER_PROGRAM)

logger = logging.getLogger(__name__)


class ProtImodDoseFilter(ProtImodBase, ProtStreamingBase):
    """
    Tilt-series dose filtering based on the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/mtffilter.html

    A specialized filter can be applied to perform dose weight-filtering of
    cryoEM images, particularly ones from tilt series.  The filter is as
    described in Grant and Grigorieff, 2015 (DOI: 10.7554/eLife.06980) and
    the implementation follows that in their "unblur" program.  At any fre-
    quency, the filter follows an exponential decay with dose, where the
    exponential is of the dose divided by 2 times a "critical dose" for
    that frequency.  This critical dose was empirically found to be approx-
    imated by a * k^b + c, where k is frequency; the values of a, b, c in
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
        form.addParallelSection(threads=2, mpi=0)

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
                cInId = self._insertFunctionStep(self.linkTsStep,
                                                 tsId,
                                                 prerequisites=[],
                                                 needsGPU=False)
                compId = self._insertFunctionStep(self.doseFilterStep,
                                                  tsId,
                                                  prerequisites=cInId,
                                                  needsGPU=False)
                outId = self._insertFunctionStep(self.createOutputStep,
                                                 tsId,
                                                 prerequisites=[compId],
                                                 needsGPU=False)
                closeSetStepDeps.append(outId)
                logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                self.tsIdReadList.append(tsId)

            self.refreshStreaming(inTsSet)

    # --------------------------- STEPS functions -----------------------------
    def doseFilterStep(self, tsId: str):
        """Apply the dose filter to every tilt series"""
        try:
            logger.info(cyanStr(f'tsId = {tsId} -> Dose filtering...'))
            with self._lock:
                ts = self.getCurrentTs(tsId)

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
                progParams['-input'] = ts.getOddFileName()
                progParams['-output'] = self.getExtraOutFile(tsId, suffix=ODD)
                self.runProgram(MTTFILTER_PROGRAM, progParams)

                progParams['-input'] = ts.getEvenFileName()
                progParams['-output'] = self.getExtraOutFile(tsId, suffix=EVEN)
                self.runProgram(MTTFILTER_PROGRAM, progParams)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {MTTFILTER_PROGRAM} execution failed with the exception -> {e}'))

    def createOutputStep(self, tsId: str):
        """Generate output filtered tilt series"""
        if tsId in self.failedItems:
            self.addToOutFailedSet(tsId)
            return

        try:
            outTsFile = self.getExtraOutFile(tsId)
            if exists(outTsFile):
                with self._lock:
                    ts = self.getCurrentTs(tsId)
                    outTsSet = self.getOutputSetOfTS(self.getInputTsSet(pointer=True))
                    outTs = TiltSeries()
                    outTs.copyInfo(ts)
                    self.updateTsAcquisition(outTs)  # Acquisition dose goes to 0 after having been applied
                    outTsSet.append(outTs)
                    for ti in ts.iterItems():
                        outTi = TiltImage()
                        outTi.copyInfo(ti)
                        outTi.setFileName(self.getExtraOutFile(tsId))
                        self.updateTiAcquisition(outTi)
                        if self.doOddEven:
                            outTi.setOddEven([self.getExtraOutFile(tsId, suffix=ODD),
                                              self.getExtraOutFile(tsId, suffix=EVEN)])
                        else:
                            outTi.setOddEven([])  # the input may have odd/even but the user may have decided not
                            # to consider them in the current execution, so they should be set to empty to avoid
                            # next protocols be confused about having them.
                        outTs.append(outTi)
                    outTs.write()
                    outTsSet.update(outTs)
                    outTsSet.write()
                    self._store(outTsSet)
                    # Close explicitly the outputs (for streaming)
                    for outputName in self._possibleOutputs.keys():
                        output = getattr(self, outputName, None)
                        if output:
                            output.close()
            else:
                logger.error(redStr(f'tsId = {tsId} -> Output file {outTsFile} was not generated. Skipping... '))
        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        if self.inputDoseType.get() == SCIPION_IMPORT:
            for ts in self.getInputTsSet():
                if ts.getFirstItem().getAcquisition().getDosePerFrame() is None:
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
