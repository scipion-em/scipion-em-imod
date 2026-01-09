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
from typing import Tuple
import pyworkflow.protocol.params as params
from imod.convert.convert import genXfFile
from imod.protocols.protocol_base import BINNING_FACTOR, NEWSTACK_PROGRAM
from pwem import ALIGN_NONE
from pwem.convert.headers import setMRCSamplingRate
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.utils import Message, redStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from imod.protocols import ProtImodBase
from imod.constants import (XF_EXT, ODD, EVEN, OUTPUT_TS_INTERPOLATED_NAME)

logger = logging.getLogger(__name__)


class ProtImodApplyTransformationMatrix(ProtImodBase):
    """
    Compute the interpolated tilt-series from its transform matrix.
    The protocol makes use of the IMod command newstack
    More info:
        https://bio3d.colorado.edu/imod/doc/man/newstack.html

    Generally, the tilt series has an associated transformation matrix
    which contains the alignment information. The transformation matrix
    is usually associated but not applied to avoid to accumulate interpolation
    errors during the image processing. This protocol allows to apply
    the transformation matrix to the tilt series
    """
    _label = 'Apply transformation'
    _possibleOutputs = {OUTPUT_TS_INTERPOLATED_NAME: SetOfTiltSeries}
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsDict = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        super().addInTsSetFormParam(form)
        form.addParam(BINNING_FACTOR, params.IntParam,
                      default=1,
                      label='Binning for the interpolated',
                      help='Binning to be applied to the interpolated tilt-series '
                           'in IMOD convention.\nBinning is a scaling factor '
                           'given by an integer greater than 1. IMOD uses ordinary '
                           'binning (with antialiasing filter) to reduce images in '
                           'size by the given factor. The value of a binned pixel '
                           'is the average of pixel values in each block of pixels '
                           'being binned. Binning is applied before all other image '
                           'transformations.')
        form.addParam('taperInside',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=True,
                      label='Taper inwards from the edge?',
                      help='When the image is transformed areas with no information '
                           'are filled in (e.g. because of rotation). '
                           'Decide whether tapering is done inwards or outwards '
                           'from the edge.')
        form.addParam('linear',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=False,
                      label='Linear interpolation?',
                      help='From newstack man page: Use linear instead of cubic '
                           'interpolation to transform images. Linear interpolation '
                           'is more suitable when images are very noisy, but cubic '
                           'interpolation will preserve fine detail better when '
                           'noise is not an issue.')
        self.addOddEvenParams(form)
        form.addParallelSection(threads=2, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []
        for tsId in self.tsDict.keys():
            compId = self._insertFunctionStep(self.computeAlignmentStep,
                                              tsId,
                                              prerequisites=[],
                                              needsGPU=False)
            outId = self._insertFunctionStep(self.createOutputStep,
                                             tsId,
                                             prerequisites=[compId],
                                             needsGPU=False)
            closeSetStepDeps.append(outId)
        self._insertFunctionStep(self.closeOutputSetsStep,
                                 OUTPUT_TS_INTERPOLATED_NAME,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    # --------------------------- STEPS functions ------------------------------
    def _initialize(self):
        super()._initialize()
        self.tsDict = {ts.getTsId(): ts.clone() for ts in self.getInputTsSet().iterItems()}

    def computeAlignmentStep(self, tsId: str):
        try:
            logger.info(f"tsId = {tsId}: computing the alignment...")
            ts = self.tsDict[tsId]
            firstItem = ts.getFirstEnabledItem()
            self.genTsPaths(tsId)
            inputFile = firstItem.getFileName()
            outputFile = self.getExtraOutFile(tsId)
            # Gen the xf alignment file
            xfFile = self.getExtraOutFile(tsId, ext=XF_EXT)
            genXfFile(ts, xfFile)
            # Get the excluded views indices
            tsExcludedIndices = ts.getTsExcludedViewsIndices(ts.getTsPresentAcqOrders())
            # Get the doSwap argument value
            doSwap = self.getNewstackDoSwap(firstItem, xfFile)
            self._runNewstack(ts, inputFile, outputFile, xfFile, doSwap, tsExcludedIndices)
            if self.doOddEven:
                # Odd
                logger.info(f"tsId = {tsId} ODD: computing the alignment...")
                inputFile = ts.getOddFileName()
                outputFile = self.getExtraOutFile(tsId, suffix=ODD)
                self._runNewstack(ts, inputFile, outputFile, xfFile, doSwap, tsExcludedIndices)
                # Even
                logger.info(f"tsId = {tsId} EVEN: computing the alignment...")
                inputFile = ts.getEvenFileName()
                outputFile = self.getExtraOutFile(tsId, suffix=EVEN)
                self._runNewstack(ts, inputFile, outputFile, xfFile, doSwap, tsExcludedIndices)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {NEWSTACK_PROGRAM} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def createOutputStep(self, tsId: str):
        ts = self.tsDict[tsId]
        if tsId in self.failedItems:
            self.addToOutFailedSet(ts)
            return

        try:
            outputFn = self.getExtraOutFile(tsId)
            if exists(outputFn):
                self._registerOutput(ts)
            else:
                logger.error(redStr(f'tsId = {tsId} -> Output file {outputFn} was not generated. Skipping... '))

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    @retry_on_sqlite_lock(log=logger)
    def _registerOutput(self, ts: TiltSeries):
        with self._lock:
            # Set of tilt-series
            outTsSet = self._prepareOutputSet()
            # Tilt-series
            outTs = self._createOutputTiltSeries(ts)
            outTsSet.append(outTs)
            # Tilt-images
            tiList, angleMin, angleMax, accumDose, initialDose = self._processTiltImages(ts)
            self._updateAcquisition(ts, outTs, tiList, angleMin, angleMax, accumDose, initialDose)
            # Data persistence
            outTs.write()
            outTsSet.update(outTs)
            outTsSet.write()
            self._store(outTsSet)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []
        for ts in self.getInputTsSet():
            if not ts.hasAlignment():
                validateMsgs.append("Some tilt-series from the input set "
                                    "are missing a transformation matrix.")
                break
        return validateMsgs

    def _summary(self):
        summary = []
        output = getattr(self, OUTPUT_TS_INTERPOLATED_NAME, None)
        if output is not None:
            summary.append(f"Input tilt-series: {self.getInputTsSet().getSize()}\n"
                           f"Interpolations applied: {output.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        output = getattr(self, OUTPUT_TS_INTERPOLATED_NAME, None)
        if output is not None:
            methods.append("The interpolation has been computed for "
                           f"{output.getSize()} "
                           "tilt-series using the IMOD *newstack* command.")
        return methods

    # --------------------------- UTILS functions -----------------------------
    def _runNewstack(self,
                    ts: TiltSeries,
                    inputFile: str,
                    outputFile: str,
                    xfFile: str,
                    doSwap: bool,
                    tsExcludedIndices: set) -> None:

        # Generate the command file for newstack
        paramsDict = self.getBasicNewstackParams(ts,
                                                 inputFile,
                                                 outputFile,
                                                 xfFile=xfFile,
                                                 tsExcludedIndices=tsExcludedIndices,
                                                 binning=self.binning.get(),
                                                 doSwap=doSwap,
                                                 doTaper=True)
        paramsDict["-taper"] = "1,1" if self.taperInside else "1,0"
        if self.linear:
            paramsDict["-linear"] = ""
        self.runProgram(NEWSTACK_PROGRAM, paramsDict)

    @staticmethod
    def updateTiltSeries(tsOut: TiltSeries) -> None:
        tsOut.setInterpolated(True)
        tsOut.setAlignment(ALIGN_NONE)
        tsOut.getAcquisition().setTiltAxisAngle(0.)  # 0 because TS is aligned

    def _prepareOutputSet(self) -> SetOfTiltSeries:
            return self.getOutputSetOfTS(
                self.getInputTsSet(pointer=True),
                binning=self.binning.get(),
                attrName=OUTPUT_TS_INTERPOLATED_NAME,
                suffix="Interpolated")

    def _createOutputTiltSeries(self, ts: TiltSeries) -> TiltSeries:
        outTs = TiltSeries()
        outTs.copyInfo(ts)
        self.updateTiltSeries(outTs)
        return outTs

    def _processTiltImages(self, ts: TiltSeries) -> Tuple[list, float, float, float, float]:
        tsId = ts.getTsId()
        angleMin, angleMax = 999, -999
        accumDose, initialDose = 0, 999
        tiFn = self.getExtraOutFile(tsId)
        # Update the sampling rate in the file header
        setMRCSamplingRate(tiFn, ts.getSamplingRate() * self.binning.get())
        tiList = []

        for ti in ts.iterItems(orderBy=TiltImage.INDEX_FIELD):
            if not ti.isEnabled():
                continue

            tiAngle = ti.getTiltAngle()
            angleMin = min(tiAngle, angleMin)
            angleMax = max(tiAngle, angleMax)
            accumDose = max(ti.getAcquisition().getAccumDose(), accumDose)
            initialDose = min(ti.getAcquisition().getDoseInitial(), initialDose)

            outTi = TiltImage()
            outTi.copyInfo(ti)
            outTi.setFileName(tiFn)
            outTi.getAcquisition().setTiltAxisAngle(0.)
            outTi.setTransform(None)
            self.setTsOddEven(tsId, outTi, binGenerated=True)

            tiList.append(outTi)

        return tiList, angleMin, angleMax, accumDose, initialDose


    @staticmethod
    def _updateAcquisition(ts: TiltSeries,
                           outTs: TiltSeries,
                           tiList: list,
                           angleMin: float,
                           angleMax: float,
                           accumDose: float,
                           initialDose: float) -> None:
        # Update the acquisition minAngle and maxAngle values of the tilt-series
        if ts.hasExcludedViews():
            acq = outTs.getAcquisition()
            acq.setAngleMin(angleMin)
            acq.setAngleMax(angleMax)
            acq.setAccumDose(accumDose)
            acq.setDoseInitial(initialDose)
            outTs.setAcquisition(acq)
            # Update the acquisition minAngle and maxAngle values of each tilt-image acq while preserving their
            # specific accum and initial dose values
            for tiOut in tiList:
                tiAcq = tiOut.getAcquisition()
                tiAcq.setAngleMin(angleMin)
                tiAcq.setAngleMax(angleMax)
                tiOut.setAcquisition(tiAcq)
                outTs.append(tiOut)

            outTs.setAnglesCount(len(outTs))
        else:
            for tiOut in tiList:
                outTs.append(tiOut)

