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
import numpy as np
from typing_extensions import List
import pyworkflow.protocol.params as params
from imod.convert.convert import readXfFile
from imod.protocols.protocol_base import ProtImodBase
from imod.protocols.protocol_base_xcorr_fidmodel import ProtImodBaseXcorrFidModel
from pwem import genExecStatusDir, getExecStatusDir, appendStreamItem
from pwem.objects import Transform
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.utils import Message, cyanStr, redStr, yellowStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from imod.constants import (TLT_EXT, PREXF_EXT, PREXG_EXT,
                            OUTPUT_TILTSERIES_NAME,
                            OUTPUT_TS_INTERPOLATED_NAME, XFTOXG_PROGRM, TILT_XCORR_PROGRAM)
from tomo.protocols.protocol_base_streaming_tomo import ProtocolBaseStreamingTomo
from tomo.utils import sleepRandomly, writeTsSidecar

logger = logging.getLogger(__name__)


class ProtImodXcorrPrealignment(ProtImodBase, ProtImodBaseXcorrFidModel, ProtocolBaseStreamingTomo):
    """
    Tilt-series cross correlation alignment based on the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/tiltxcorr.html

    Tiltxcorr uses cross-correlation to find an initial translational
    alignment between successive images of a tilt series.  For a given pair
    of images, it stretches the image with the larger tilt angle perpendic-
    ular to the tilt axis, by an amount equal to the ratio of the cosines
    of the two tilt angles (cosine stretch).  The stretched image is corre-
    lated with the other image, and the position of the peak of the corre-
    lation indicates the relative shift between the images.

    xftoxg program is used after tiltxcorr:
    More info:
        https://bio3d.colorado.edu/imod/doc/man/xftoxg.html

   Xftoxg takes a list of transformations (f) from each section to the
   previous one, and computes a list of xforms (g) to apply to each sec-
   tion to obtain a single consistent set of alignments.  Transforms can
   be simple 6-component linear transforms or warping transformations.
    """

    _label = 'Coarse prealignment'
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
        form.addParam('cumulativeCorr',
                      params.BooleanParam,
                      default=False,
                      label='Use cumulative correlation?',
                      help='The program will take the image at zero tilt as the first'
                           'reference, and correlate it with the image at the next '
                           'most negative tilt. It will then add the aligned image '
                           'to the first reference to make the reference for the next '
                           'tilt. At each tilt, the reference will be the sum of '
                           'images that have already been aligned. When the most '
                           'negative tilt angle is reached, the procedure is repeated '
                           'from the zero-tilt view to more positive tilt angles.')
        form.addParam('Trimming parameters', params.LabelParam,
                      label='Tilt axis angle detected from import. In case another '
                            'value is desired please adjust the number below.')
        form.addParam('tiltAxisAngle',
                      params.FloatParam,
                      allowsNull=True,
                      label='Tilt axis angle (degrees)',
                      help='The tilt axis angle is the tilt axis rotation relative '
                           'to the Y axis of the image. If it was not properly set '
                           'in the import of the tilt series, or the imported'
                           'information is not correct you have the chance to '
                           'correct at in this point. Usually, it will be 90 '
                           'degrees less than the RotationAngle in a system with '
                           'no axis inversions.')
        trimming = form.addGroup('Trimming parameters',
                                 expertLevel=params.LEVEL_ADVANCED)
        self.addTrimingParams(trimming,
                              pxTrimCondition=False,
                              correlationCondition=True,
                              levelType=params.LEVEL_ADVANCED)
        self.filteringParametersForm(form,
                                     condition=True,
                                     levelType=params.LEVEL_ADVANCED)
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

    # Streaming Hooks ############################
    def _streamingInitialize(self):
        self._initialize()

    def _getStreamingInputSets(self):
        return [self.getInputTsSet()]

    def _getProcessedTsIds(self):
        return self.tsIdReadList

    def _getStreamingOutputNames(self):
        return OUTPUT_TILTSERIES_NAME

    # End of streaming hooks #####################

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
        convId = self._insertFunctionStep(self.convertInStep,
                                          ts,
                                          prerequisites=[],
                                          needsGPU=False)
        compId = self._insertFunctionStep(self.computeXcorrStep,
                                          ts,
                                          prerequisites=[convId],
                                          needsGPU=False)
        outId = self._insertFunctionStep(self.createOutputStep,
                                         ts,
                                         prerequisites=[compId],
                                         needsGPU=False)
        closeSetStepDeps.append(outId)

    # --------------------------- STEPS functions -----------------------------
    def convertInStep(self, ts: TiltSeries):
        self.convertInputStep(ts, ts.getTsPresentAcqOrders())

    def computeXcorrStep(self, ts: TiltSeries):
        """Compute transformation matrix for each tilt series. """
        tsId = ts.getTsId()
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'tsId = {tsId} -> Correcting the translations with {TILT_XCORR_PROGRAM}...'))
                tiltAxisAngle = self.getTiltAxisOrientation(ts)

                paramsXcorr = {
                    "-input": self.getTmpOutFile(tsId),
                    "-output": self.getExtraOutFile(tsId, ext=PREXF_EXT),
                    "-tiltfile": self.getExtraOutFile(tsId, ext=TLT_EXT),
                    "-RotationAngle": tiltAxisAngle,
                    "-FilterSigma1": self.filterSigma1.get(),
                    "-FilterSigma2": self.filterSigma2.get(),
                    "-FilterRadius1": self.filterRadius1.get(),
                    "-FilterRadius2": self.filterRadius2.get()
                }

                if self.cumulativeCorr:
                    paramsXcorr["-CumulativeCorrelation"] = ""

                doTrim = any([getattr(self, attr).hasValue() for
                              attr in ["xmin", "xmax", "ymin", "ymax"]])
                if doTrim:
                    xdim, ydim, _ = ts.getDim()
                    xmin, xmax = self.xmin.get() or 0, self.xmax.get() or xdim - 1
                    ymin, ymax = self.ymin.get() or 0, self.ymax.get() or ydim - 1

                    paramsXcorr["-xminmax"] = f"{xmin},{xmax}"
                    paramsXcorr["-yminmax"] = f"{ymin},{ymax}"

                self.runProgram(TILT_XCORR_PROGRAM, paramsXcorr)

                paramsXftoxg = {
                    "-input": self.getExtraOutFile(tsId, ext=PREXF_EXT),
                    "-goutput": self.getExtraOutFile(tsId, ext=PREXG_EXT),
                    "-NumberToFit": 0
                }
                self.runProgram(XFTOXG_PROGRM, paramsXftoxg)

            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> {TILT_XCORR_PROGRAM} or {XFTOXG_PROGRM} execution '
                                    f'failed with the exception -> {e}'))
                logger.error(traceback.format_exc())

    def createOutputStep(self, inTs: TiltSeries):
        tsId = inTs.getTsId()
        if tsId in self.failedItems:
            self.addToOutFailedSet(inTs)
            return

        try:
            outTsFile = self.getExtraOutFile(tsId, ext=PREXG_EXT)
            if not exists(outTsFile):
                logger.error(redStr(f'tsId = {tsId} -> Output file {outTsFile} was not generated. Skipping... '))
                return

            outTs = TiltSeries()
            outTs.copyInfo(inTs)
            tAx = self.getTiltAxisOrientation(inTs)
            outTs.getAcquisition().setTiltAxisAngle(tAx)
            outTs.setAlignment2D()

            aliMatrixStack = readXfFile(outTsFile)
            inTiltList = inTs.loadTiltImgsInMemory()
            inTiltList.sort(key=lambda item: item.getIndex())
            tiltImages = []
            stackIndex = 0
            for inTi in inTiltList:
                outTi = TiltImage()
                outTi.copyInfo(inTi)
                if inTi.isEnabled():
                    self.updateTiltImage(outTi, stackIndex, aliMatrixStack, tAx)
                    stackIndex += 1
                else:
                    self.updateDisabledTi(outTi)
                tiltImages.append(outTi)

            self._registerOutput(outTs, tiltImages)

            # Streaming only: publish the per-TS metadata sidecar (built from the
            # in-memory ts/tiltImages, no DB read) and the journal id
            execStatusDir = getExecStatusDir(self)
            if exists(execStatusDir):
                writeTsSidecar(execStatusDir, outTs, tiltImages)
                appendStreamItem(self, tsId)

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    @retry_on_sqlite_lock(log=logger)
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
    def _summary(self):
        summary = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            summary.append(f"Input tilt-series: {self.getInputTsSet().getSize()}\n"
                           "Transformation matrices calculated: "
                           f"{output.getSize()}")

            interpTS = getattr(self, OUTPUT_TS_INTERPOLATED_NAME, None)
            if interpTS is not None:
                summary.append("Interpolated tilt-series: "
                               f"{interpTS.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            methods.append(f"The transformation matrix has been calculated for "
                           f"{output.getSize()} tilt-series using "
                           f"the IMOD *{TILT_XCORR_PROGRAM}* command.")
        return methods

    # --------------------------- UTILS functions ------------------------------
    def getTiltAxisOrientation(self, ts: TiltSeries) -> float:
        if self.tiltAxisAngle.hasValue():
            return self.tiltAxisAngle.get()
        else:
            return ts.getAcquisition().getTiltAxisAngle()

    @staticmethod
    def updateTiltImage(outTi: TiltImage,
                        stackIndex: int,
                        aliMatrixStack: np.array,
                        tiltAxisAngle: float) -> None:
        transform = Transform()
        newTransform = aliMatrixStack[:, :, stackIndex]
        newTransformArray = np.array(newTransform)

        if outTi.hasTransform():
            previousTransform = outTi.getTransform().getMatrix()
            previousTransformArray = np.array(previousTransform)
            outputTransformMatrix = np.matmul(newTransformArray, previousTransformArray)
            transform.setMatrix(outputTransformMatrix)
        else:
            transform.setMatrix(newTransformArray)

        outTi.setTransform(transform)
        if tiltAxisAngle:
            outTi.getAcquisition().setTiltAxisAngle(tiltAxisAngle)

    @staticmethod
    def updateDisabledTi(outTi: TiltImage) -> None:
        transform = Transform()
        trMatrix = np.eye(3)
        transform.setMatrix(trMatrix)
        outTi.setTransform(transform)
