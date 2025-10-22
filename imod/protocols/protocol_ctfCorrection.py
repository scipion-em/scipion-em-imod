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
import traceback
from os.path import exists
from typing import Set, Tuple
from imod.protocols.protocol_base import IN_CTF_TOMO_SET
from imod.protocols.protocol_base_ts_align import ProtImodBaseTsAlign
from pwem import ALIGN_NONE
import pyworkflow.protocol.params as params
from pwem.convert.headers import setMRCSamplingRate
from pyworkflow.protocol import STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.utils import Message, yellowStr, redStr, cyanStr
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries, CTFTomoSeries
from tomo.utils import getCommonTsAndCtfElements
from imod import utils
from imod.constants import (DEFOCUS_EXT, TLT_EXT, XF_EXT, ODD,
                            EVEN, OUTPUT_TILTSERIES_NAME, CTF_PHASE_FLIP_PROGRAM, OUTPUT_CTF_SERIE)

logger = logging.getLogger(__name__)


class ProtImodCtfCorrection(ProtImodBaseTsAlign, ProtStreamingBase):
    """
    CTF correction of a set of input tilt-series using the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/ctfphaseflip.html

    This program will correct the CTF of an input tilt series by phase
    flipping, with an option to attenuate frequencies near the zeros of the
    CTF.

    Ctfphaseflip corrects each view strip by strip.  A strip is defined as
    an image region whose defocus difference is less than a user specified
    value, the defocus tolerance.  Normally, the strips are vertically ori-
    ented and defocus is assumed to be the same along a vertical line.
    Thus, the tilt series must be aligned so that the tilt axis is vertical
    before applying this correction.  The original thinking was that an
    image region with defocus difference less than the tolerance could be
    considered to have constant defocus and could be corrected as one
    strip.  However, the validity of the correction at the center of the
    strip probably does not depend on whether it contains material beyond
    this focus range, since only vertical lines near or at the center are
    used in the corrected image.  The program may limit the width further
    to reduce computation time, or expand it to retain enough resolution
    between successive zeros in the X direction of frequency space.

    Through most of the image, each strip is corrected based on the defocus
    at the center of the strip.  However, the strips at the left and right
    edges of the image may be corrected repeatedly, at different defocus
    values, in order to extend the correction close enough to the edges of
    the image.


    """

    _label = 'CTF correction'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.ctfTsIdReadList = []

    @classmethod
    def worksInStreaming(cls):
        """ So far none of them work in streaming. """
        return True

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        super().addInTsSetFormParam(form)
        form.addParam(IN_CTF_TOMO_SET,
                      params.PointerParam,
                      label="Input CTF estimation",
                      pointerClass='SetOfCTFTomoSeries',
                      help='Select the CTF estimation for the set '
                           'of tilt-series.')
        form.addParam('defocusTol',
                      params.IntParam,
                      label='Defocus tolerance (nm)',
                      default=200,
                      important=True,
                      help='The value introduced must be the same used for '
                           'CTF estimation with IMOD.\n\n'
                           'Defocus tolerance in nanometers defining the '
                           'center strips. The center strips are taken '
                           'from the central region of a view that has defocus '
                           'difference less than this tolerance. '
                           'These kind of center strips from all views within '
                           'AngleRange are considered to have a '
                           'constant defocus and are used to compute the '
                           'initial CTF after being further tessellated '
                           'into tiles.')
        form.addParam('interpolationWidth',
                      params.IntParam,
                      label='Interpolation width (px)',
                      default=15,
                      important=True,
                      help="The distance in pixels between the center lines "
                           "of two consecutive strips. A pixel inside the "
                           "region between those two center lines resides in "
                           "both strips. As the two strips are corrected "
                           "separately, that pixel will have 2 corrected "
                           "values. The final value for that pixel is a "
                           "linear interpolation of the 2 corrected "
                           "values. If a value of 1 is entered, there is "
                           "no such interpolation. For a value greater "
                           "than one, the entered value will be used "
                           "whenever the strip width is less than 256 "
                           "(i.e., at high tilt), and the value will be "
                           "scaled proportional to the strip width for widths "
                           "above 256.  This scaling keeps the computational "
                           "time down and is reasonable because the defocus "
                           "difference between adjacent wide strips at "
                           "wider intervals is still less than that between "
                           "the narrower strips at high tilt. However, strips "
                           "at constant spacing can still be obtained by "
                           "entering the negative of the desired spacing, "
                           "which disables the scaling of the spacing.")
        form.addHidden(params.USE_GPU,
                       params.BooleanParam,
                       default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation."
                            "Select the one you want to use.")
        form.addHidden(params.GPU_LIST,
                       params.StringParam,
                       default='0',
                       label="Choose GPU IDs",
                       help="GPU ID. To pick the best available one set 0. "
                            "For a specific GPU set its number ID "
                            "(starting from 1).")
        self.addOddEvenParams(form)
        form.addParallelSection(threads=3, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def stepsGeneratorStep(self) -> None:
        """
        This step should be implemented by any streaming protocol.
        It should check its input and when ready conditions are met
        call the self._insertFunctionStep method.
        """
        closeSetStepDeps = []
        self._initialize()
        inTsSet = self.getInputTsSet()
        self.readingOutput(getattr(self, OUTPUT_TILTSERIES_NAME, None))
        inCtfSet = self.getInputCtfSet()
        self.readingOutput(getattr(self, OUTPUT_CTF_SERIE, None), tsIdListName='ctfTsIdList')

        while True:
            listTsIdInput = inTsSet.getTSIds()
            listCtfTsIdInput = inCtfSet.getTSIds()
            if ((not inTsSet.isStreamOpen() and self.tsIdReadList == listTsIdInput) and
                    (not inCtfSet.isStreamOpen() and self.ctfTsIdReadList == listCtfTsIdInput)):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetsStep,
                                         OUTPUT_TILTSERIES_NAME,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break
            closeSetStepDeps = []
            for ts in inTsSet.iterItems():
                tsId = ts.getTsId()
                if tsId not in self.tsIdReadList and ts.getSize() > 0:
                    try:
                        ctf = self.getCurrentCtf(tsId)
                    except Exception as e:
                        logger.info(yellowStr(f'tsId = {tsId} - no corresponding CTF was found...'))
                        logger.error(f'{e}')
                        logger.error(traceback.format_exc())
                        continue
                    if tsId not in self.ctfTsIdReadList and ctf.getSize() > 0:  # Avoid processing empty TS (before the Tis are added)
                        pidConvert = self._insertFunctionStep(self.convertInStep,
                                                              tsId,
                                                              prerequisites=[],
                                                              needsGPU=False)
                        pidProcess = self._insertFunctionStep(self.ctfCorrection,
                                                              tsId,
                                                              prerequisites=pidConvert,
                                                              needsGPU=True)
                        pidCreateOutput = self._insertFunctionStep(self.createOutputStep,
                                                                   tsId,
                                                                   prerequisites=pidProcess,
                                                                   needsGPU=False)
                        closeSetStepDeps.append(pidCreateOutput)
                        logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                        self.tsIdReadList.append(tsId)
                        self.ctfTsIdReadList.append(tsId)

            self.refreshStreaming(inTsSet)
            self.refreshStreaming(inCtfSet)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        super()._initialize()
        tsSet = self.getInputTsSet()
        self.sRate = tsSet.getSamplingRate()
        self.acq = tsSet.getAcquisition()

    def convertInStep(self, tsId: str):
        try:
            self.genTsPaths(tsId)
            with self._lock:
                ts = self.getCurrentTs(tsId)
                ctf = self.getCurrentCtf(tsId)
            presentAcqOrders = getCommonTsAndCtfElements(ts, ctf)
            # Generate the defocus file
            self._generateDefocusFile(ts, ctf, presentAcqOrders=presentAcqOrders)
            # Generate the alignment files
            super().convertInputStep(tsId, presentAcqOrders=presentAcqOrders)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'{e}'))
            logger.error(traceback.format_exc())

    def ctfCorrection(self, tsId: str):
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'tsId = {tsId}: correcting the CTF...'))
                with self._lock:
                    ts = self.getCurrentTs(tsId)

                paramsCtfPhaseFlip = {
                    "-InputStack": self.getTmpOutFile(tsId),
                    "-AngleFile": self.getExtraOutFile(tsId, ext=TLT_EXT),
                    "-OutputFileName": self.getExtraOutFile(tsId),
                    "-DefocusFile": self.getExtraOutFile(tsId, ext=DEFOCUS_EXT),
                    "-Voltage": int(self.acq.getVoltage()),
                    "-SphericalAberration": self.acq.getSphericalAberration(),
                    "-DefocusTol": self.defocusTol.get(),
                    "-PixelSize": self.sRate / 10,  # nm
                    "-AmplitudeContrast": self.acq.getAmplitudeContrast(),
                    "-InterpolationWidth": self.interpolationWidth.get()
                }

                if self.usesGpu():
                    paramsCtfPhaseFlip["-UseGPU"] = self.getGpuList()[0]
                    paramsCtfPhaseFlip["-ActionIfGPUFails"] = "2,2"

                if ts.hasAlignment():
                    paramsCtfPhaseFlip["-TransformFile"] = self.getExtraOutFile(tsId, ext=XF_EXT)

                self.runProgram(CTF_PHASE_FLIP_PROGRAM, paramsCtfPhaseFlip)

                if self.doOddEven:
                    # ODD
                    paramsCtfPhaseFlip["-InputStack"] = self.getTmpOutFile(tsId, suffix=ODD)
                    paramsCtfPhaseFlip["-OutputFileName"] = self.getExtraOutFile(tsId, suffix=ODD)
                    self.runProgram(CTF_PHASE_FLIP_PROGRAM, paramsCtfPhaseFlip)

                    # EVEN
                    paramsCtfPhaseFlip["-InputStack"] = self.getTmpOutFile(tsId, suffix=EVEN)
                    paramsCtfPhaseFlip["-OutputFileName"] = self.getExtraOutFile(tsId, suffix=EVEN)
                    self.runProgram(CTF_PHASE_FLIP_PROGRAM, paramsCtfPhaseFlip)

            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> {CTF_PHASE_FLIP_PROGRAM} execution failed '
                                    f'with the exception -> {e}'))
                logger.error(traceback.format_exc())

    def createOutputStep(self, tsId: str):
        if tsId in self.failedItems:
            self.addToOutFailedSet(tsId)
            return
        try:
            outputFn = self.getExtraOutFile(tsId)
            if exists(outputFn):
                with self._lock:
                    ts = self.getCurrentTs(tsId)
                    ctf = self.getCurrentCtf(tsId)
                    presentAcqOrders = getCommonTsAndCtfElements(ts, ctf)
                    # Set of tilt-series
                    inTsSetPointer = self.getInputTsSet(pointer=True)
                    outTsSet = self.getOutputSetOfTS(inTsSetPointer)
                    # Tilt-series
                    outTs = self._createOutputTiltSeries(ts, presentAcqOrders)
                    outTsSet.append(outTs)
                    # Tilt-images
                    tiList, angleMin, angleMax = self._processTiltImages(ts, presentAcqOrders, outputFn)
                    self._updateAcquisition(ts, ctf, presentAcqOrders, outTs, tiList, angleMin, angleMax)
                    setMRCSamplingRate(outputFn, ts.getSamplingRate())  # Update the apix value in file header
                    # Data persistence
                    outTs.write()
                    outTsSet.update(outTs)
                    outTsSet.write()
                    self._store(outTsSet)
                    # Close explicitly the outputs (for streaming)
                    self.closeOutputsForStreaming()
            else:
                logger.error(f'tsId = {tsId} -> Output file {outputFn} was not generated. Skipping... ')
        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    # --------------------------- UTILS functions -----------------------------
    def _generateDefocusFile(self, ts: TiltSeries,
                             ctf: CTFTomoSeries,
                             presentAcqOrders: Set[int]) -> None:
        tsId = ts.getTsId()
        self.debug(f"tsId = {tsId} -> Generating defocus file...")
        defocusFilePath = self.getExtraOutFile(tsId, ext=DEFOCUS_EXT)
        utils.genDefocusFileFromScipion(ctf,
                                                defocusFilePath,
                                                inputTiltSeries=ts,
                                                presentAcqOrders=presentAcqOrders)

    @staticmethod
    def _createOutputTiltSeries(ts: TiltSeries, presentAcqOrders: Set[int]) ->TiltSeries:
        outTs = TiltSeries()
        outTs.copyInfo(ts)
        outTs.setAlignment(ALIGN_NONE)
        outTs.setAnglesCount(len(presentAcqOrders))
        outTs.setCtfCorrected(True)
        outTs.setInterpolated(True)
        outTs.getAcquisition().setTiltAxisAngle(0.)  # 0 because TS is aligned
        return outTs

    def _processTiltImages(self,
                           ts: TiltSeries,
                           presentAcqOrders: Set[int],
                           outputFn: str) -> Tuple[list, float, float]:
        tsId = ts.getTsId()
        angleMin, angleMax = 999, -999
        tiList = []
        for index, inTi in enumerate(ts.iterItems()):
            if inTi.getAcquisitionOrder() in presentAcqOrders:
                outTi = TiltImage()
                outTi.copyInfo(inTi, copyTM=False)
                acq = inTi.getAcquisition()
                acq.setTiltAxisAngle(0.)  # Is interpolated
                outTi.setAcquisition(acq)
                outTi.setFileName(outputFn)
                self.setTsOddEven(tsId, outTi, binGenerated=True)
                # Update the acquisition of the TS. The accumDose, angle min and angle max for the re-stacked TS, as
                # these values may change if the removed tilt-images are the first or the last, for example.
                tiAngle = outTi.getTiltAngle()
                angleMin = min(tiAngle, angleMin)
                angleMax = max(tiAngle, angleMax)
                tiList.append(outTi)

        return tiList, angleMin, angleMax

    @staticmethod
    def _updateAcquisition(ts: TiltSeries,
                           ctf: CTFTomoSeries,
                           presentAcqOrders: Set[int],
                           outTs: TiltSeries,
                           tiList: list,
                           angleMin: float,
                           angleMax: float) -> None:
        if len(presentAcqOrders) != max(len(ts), len(ctf)):
            # Update the acquisition minAngle and maxAngle values of the tilt-series
            acq = outTs.getAcquisition()
            acq.setAngleMin(angleMin)
            acq.setAngleMax(angleMax)
            acq.setAccumDose(0)
            acq.setDoseInitial(0)
            outTs.setAcquisition(acq)
            # Update the acquisition minAngle and maxAngle values of each tilt-image acq while preserving their
            # specific accum and initial dose values
            for tiOut in tiList:
                tiAcq = tiOut.getAcquisition()
                tiAcq.setAngleMin(angleMin)
                tiAcq.setAngleMax(angleMax)
                tiAcq.setAccumDose(0)
                tiAcq.setDoseInitial(0)
                outTs.append(tiOut)
            outTs.setAnglesCount(len(outTs))
        else:
            for tiOut in tiList:
                outTs.append(tiOut)

    # --------------------------- INFO functions ------------------------------
    def _warnings(self):
        warnings = []
        for ts in self.getInputTsSet():
            if not ts.hasAlignment():
                warnings.append(f"Input tilt-series {ts.getTsId()} does not have "
                                "alignment information! The recommended workflow is to "
                                "estimate CTF on raw tilt-series and then here "
                                "provide tilt-series with alignment "
                                "(non-interpolated).")
                break

        return warnings

    def _summary(self):
        summary = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            summary.append(f"Input tilt-series: {self.getInputTsSet().getSize()}\n"
                           f"CTF corrections applied: {output.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            methods.append(f"{output.getSize()} tilt-series have been "
                           "CTF corrected using the IMOD *ctfphaseflip* program.")
        return methods
