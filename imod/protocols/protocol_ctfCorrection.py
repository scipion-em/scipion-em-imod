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
from typing import Set, Tuple, List
from imod.protocols.protocol_base import IN_CTF_TOMO_SET, ProtImodBase
from pwem import ALIGN_NONE, getExecStatusDir, appendStreamItem
import pyworkflow.protocol.params as params
from pwem.convert.headers import setMRCSamplingRate
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.utils import Message, yellowStr, redStr, cyanStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries, CTFTomoSeries
from tomo.protocols.protocol_base_streaming_tomo import ProtocolBaseStreamingTomo
from tomo.utils import getCommonTsAndCtfElements, writeTsSidecar, getTsIdsIntersection, getTsIdsDicts
from imod import utils
from imod.constants import (DEFOCUS_EXT, TLT_EXT, XF_EXT, ODD,
                            EVEN, OUTPUT_TILTSERIES_NAME, CTF_PHASE_FLIP_PROGRAM)

logger = logging.getLogger(__name__)


class ProtImodCtfCorrection(ProtImodBase, ProtocolBaseStreamingTomo):
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
        self.sRate = None
        self.acq = None
        self.tsDict = None
        self.ctfDict = None

    @classmethod
    def worksInStreaming(cls):
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
    def _insertAllSteps(self) -> None:
        tsSet = self.getInputTsSet()
        ctfSet = self.getInputCtfSet()
        if tsSet.isStreamOpen() or ctfSet.isStreamOpen():
            self._insertFunctionStep(self.stepsGeneratorStep,
                                     prerequisites=[],
                                     needsGPU=False)
        else:
            self._insertNonStreamingSteps()

    # Streaming Hooks ############################
    def _streamingInitialize(self) -> None:
        super()._initialize()
        tsSet = self.getInputTsSet()
        self.sRate = tsSet.getSamplingRate()
        self.acq = tsSet.getAcquisition()

    def _getStreamingInputSets(self):
        return [self.getInputTsSet(), self.getInputCtfSet()]

    def _discoverReadyWork(self, tsIds, inputSets):
        # Rebuild the ready TS and CTF series from their OWN producers' sidecars
        # (no live-DB read) and join by tsId. A tsId whose CTF is not yet
        # materialisable is skipped and retried next cycle.
        tsDict = self.getInputTsSet().fetchNewItems(tsIds)
        ctfDict = self.getInputCtfSet().fetchNewItems(tsIds)
        work = {}
        for tsId, ts in tsDict.items():
            ctf = ctfDict.get(tsId)
            if ctf is None:
                logger.info(yellowStr(f'tsId = {tsId} - no corresponding CTF found yet, retrying...'))
                continue
            # Payload is (ts, ctf) to match _insertCommonSteps(self, ts, ctf,
            # closeSetStepDeps); the common acq. orders are (re)computed there, so
            # they are NOT part of the payload (a 3rd element would be unpacked onto
            # the keyword-only closeSetStepDeps by the base loop).
            work[tsId] = (ts, ctf)
        return work

    # End of streaming hooks #####################

    def _insertNonStreamingSteps(self):
        closeSetStepDeps = []
        self._initialize()
        for tsId in self.tsDict.keys():
            ts = self.tsDict[tsId]
            ctf = self.ctfDict[tsId]
            self._insertCommonSteps(ts, ctf, closeSetStepDeps)
        self._insertFunctionStep(self._closeOutputSet,
                                 OUTPUT_TILTSERIES_NAME,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    def _insertCommonSteps(self, ts: TiltSeries, ctf: CTFTomoSeries, closeSetStepDeps: List[int]) -> None:
        presentAcqOrders = getCommonTsAndCtfElements(ts, ctf)
        pidConvert = self._insertFunctionStep(self.convertInStep,
                                              ts, ctf, presentAcqOrders,
                                              prerequisites=[],
                                              needsGPU=False)
        pidProcess = self._insertFunctionStep(self.ctfCorrection,
                                              ts,
                                              prerequisites=pidConvert,
                                              needsGPU=True)
        pidCreateOutput = self._insertFunctionStep(self.createOutputStep,
                                                   ts, ctf, presentAcqOrders,
                                                   prerequisites=pidProcess,
                                                   needsGPU=False)
        closeSetStepDeps.append(pidCreateOutput)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        super()._initialize()
        tsSet = self.getInputTsSet()
        ctfSet = self.getInputCtfSet()
        self.sRate = tsSet.getSamplingRate()
        self.acq = tsSet.getAcquisition()
        commonTsIds = getTsIdsIntersection(tsSet, ctfSet)
        self.tsDict, self.ctfDict = getTsIdsDicts(tsSet, ctfSet, present_ts_ids=commonTsIds)

    def convertInStep(self,
                      ts: TiltSeries,
                      ctf: CTFTomoSeries,
                      presentAcqOrders: Set[int]):
        tsId = ts.getTsId()
        try:
            self.genTsPaths(tsId)
            # Generate the defocus file
            self._generateDefocusFile(ts, ctf, presentAcqOrders=presentAcqOrders)
            # Generate the alignment files
            super().convertInputStep(ts, presentAcqOrders=presentAcqOrders)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'{e}'))
            logger.error(traceback.format_exc())

    def ctfCorrection(self, ts: TiltSeries):
        tsId = ts.getTsId()
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'tsId = {tsId}: correcting the CTF...'))
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
                    gpuId = self._stepsExecutor.getGpuList()
                    paramsCtfPhaseFlip["-UseGPU"] = gpuId[0]
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

    def createOutputStep(self,
                         ts: TiltSeries,
                         ctf: CTFTomoSeries,
                         presentAcqOrders: Set[int]):
        tsId = ts.getTsId()
        if tsId in self.failedItems:
            self.addToOutFailedSet(ts)
            return

        try:
            self.createOutTs(ts, ctf, presentAcqOrders)

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output '
                                f'with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    # --------------------------- UTILS functions -----------------------------
    def createOutTs(self,
                    inTs: TiltSeries,
                    ctf: CTFTomoSeries,
                    presentAcqOrders: Set[int]):
        tsId = inTs.getTsId()
        outputFn = self.getExtraOutFile(tsId)
        if not exists(outputFn):
            logger.error(f'tsId = {tsId} -> Output file {outputFn} was not generated. Skipping... ')
            return
        setMRCSamplingRate(outputFn, self.sRate)  # Update the apix value in file header
        # Tilt-series
        outTs = self._createOutputTiltSeries(inTs, presentAcqOrders)
        # Tilt-images
        tiList, angleMin, angleMax = self._processTiltImages(inTs, presentAcqOrders, outputFn)

        self._registerOutput(inTs, outTs, ctf, tiList, presentAcqOrders, angleMin, angleMax)

        # Streaming only: publish the per-TS metadata sidecar (built from the
        # in-memory ts/tiltImages, no DB read) and the journal id
        execStatusDir = getExecStatusDir(self)
        if exists(execStatusDir):
            writeTsSidecar(execStatusDir, outTs, tiList)
            appendStreamItem(self, tsId)

    @retry_on_sqlite_lock(log=logger)
    def _registerOutput(self,
                        inTs: TiltSeries,
                        newTs: TiltSeries,
                        ctf: CTFTomoSeries,
                        tiltImages: List[TiltImage],
                        presentAcqOrders: Set[int],
                        angleMin: float,
                        angleMax: float) -> None:
        with self._lock:
            # Set of tilt-series
            inTsSetPointer = self.getInputTsSet(pointer=True)
            outTsSet = self.getOutputSetOfTS(inTsSetPointer)
            try:
                # Tilt-series
                outTsSet.append(newTs)
                # Tilt-images
                for newTi in tiltImages:
                    newTs.append(newTi)
                # Data persistence
                self._updateAcquisition(inTs, ctf, presentAcqOrders, newTs, tiltImages, angleMin, angleMax)
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

    def _generateDefocusFile(self,
                             ts: TiltSeries,
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
    def _createOutputTiltSeries(ts: TiltSeries, presentAcqOrders: Set[int]) -> TiltSeries:
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
