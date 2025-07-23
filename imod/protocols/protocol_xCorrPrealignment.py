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
import time
from os.path import exists
import numpy as np
import pyworkflow.protocol.params as params
from imod.convert import genXfFile, genTltFile
from imod.protocols.protocol_base import IN_TS_SET
from imod.protocols.protocol_base_xcorr_fidmodel import ProtImodBaseXcorrFidModel
from pwem.objects import Transform
from pyworkflow.utils import Message, cyanStr
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from imod import utils
from imod.constants import (TLT_EXT, PREXF_EXT, PREXG_EXT,
                            OUTPUT_TILTSERIES_NAME,
                            OUTPUT_TS_INTERPOLATED_NAME, XF_EXT)

logger = logging.getLogger(__name__)

TILT_XCORR_PROGRAM = 'tiltxcorr'
XFTOXG_PROGRM = 'xftoxg'


class ProtImodXcorrPrealignment(ProtImodBaseXcorrFidModel):
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

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)

        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-series to be prealigned')

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
        form.addParallelSection(threads=2, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []
        inTsSet = self.getInputSet()
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
                if tsId not in self.tsIdReadList and ts.getSize() > 0:
                    convId = self._insertFunctionStep(self.convertInputStep,
                                                      tsId,
                                                      prerequisites=[],
                                                      needsGPU=False)
                    compId = self._insertFunctionStep(self.computeXcorrStep,
                                                      tsId,
                                                      prerequisites=[convId],
                                                      needsGPU=False)
                    outId = self._insertFunctionStep(self.createAliTsStep,
                                                     tsId,
                                                     prerequisites=[compId],
                                                     needsGPU=False)
                    closeSetStepDeps.append(outId)
                    logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                    self.tsIdReadList.append(tsId)

            time.sleep(10)
            if inTsSet.isStreamOpen():
                with self._lock:
                    inTsSet.loadAllProperties()  # refresh status for the streaming

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self, tsId, **kwargs):
        self.genTsPaths(tsId)
        with self._lock:
            ts = self.getCurrentItem(tsId)
            firstTi = ts.getFirstItem()

        # Generate the tlt file
        tltFile = self.getExtraOutFile(ts.getTsId(), ext=TLT_EXT)
        genTltFile(ts,
                   tltFile,
                   ignoreExcludedViews=True)  # The xcorr program can exclude them itself

        if firstTi.hasTransform():  # Apply a previous alignment if it exists
            logger.info(cyanStr(f'tsId = {tsId} -> Previous alignment detected. Applying...'))
            xfFile = self.getExtraOutFile(ts.getTsId(), ext=XF_EXT)
            genXfFile(ts,
                      xfFile,
                      ignoreExcludedViews=True)  # The xcorr program can exclude them itself
            doSwap = True  # A transformation will be applied
            self.runNewStackBasic(ts,
                                  xfFile=xfFile,
                                  doSwap=doSwap,
                                  ignoreExcludedViews=True)  # The xcorr program can exclude them itself

        else:  # Link it, so the input file expected by xcorr is in the same place in both sides of the "if"
            outTsFn, _, _ = self.getTmpFileNames(ts)
            self.linkTs(firstTi.getFileName(), outTsFn)

    def computeXcorrStep(self, tsId):
        """Compute transformation matrix for each tilt series. """
        try:
            logger.info(cyanStr(f'tsId = {tsId} -> Correcting the translations with {TILT_XCORR_PROGRAM}...'))
            with self._lock:
                ts = self.getCurrentItem(tsId)
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

            # Excluded views
            excludedViews = ts.getTsExcludedViewsIndices(ts.getTsPresentAcqOrders())
            if excludedViews:
                paramsXcorr["-SkipViews"] = ",".join(map(str, excludedViews))

            self.runProgram(TILT_XCORR_PROGRAM, paramsXcorr)

            paramsXftoxg = {
                "-input": self.getExtraOutFile(tsId, ext=PREXF_EXT),
                "-goutput": self.getExtraOutFile(tsId, ext=PREXG_EXT),
                "-NumberToFit": 0
            }
            self.runProgram(XFTOXG_PROGRM, paramsXftoxg)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(f'tsId = {tsId} -> {TILT_XCORR_PROGRAM} or {XFTOXG_PROGRM} execution '
                         f'failed with the exception -> {e}')

    def createAliTsStep(self, tsId):
        """ Generate tilt-series with the associated transform matrix """
        if tsId in self.failedItems:
            self.addToOutFailedSet(tsId)
        else:
            try:
                outputFn = self.getExtraOutFile(tsId, ext=PREXG_EXT)
                if exists(outputFn):
                    with self._lock:
                        ts = self.getCurrentItem(tsId)
                        tAx = self.tiltAxisAngle.get()
                        aliMatrixStack = utils.formatTransformationMatrix(outputFn)
                        outTsSet = self.getOutputSetOfTS(self.getInputSet(pointer=True),
                                                         tiltAxisAngle=tAx)
                        outTs = TiltSeries()
                        outTs.copyInfo(ts)
                        outTs.getAcquisition().setTiltAxisAngle(self.getTiltAxisOrientation(ts))
                        outTs.setAlignment2D()
                        outTsSet.append(outTs)
                        for ti in ts.iterItems(orderBy=TiltImage.INDEX_FIELD):
                            outTi = TiltImage()
                            outTi.copyInfo(ti)
                            if ti.isEnabled():
                                stackIndex = ti.getIndex() - 1
                                self.updateTiltImage(outTi, stackIndex, aliMatrixStack, tAx)
                            else:
                                self.updateDisabledTi(outTi)
                            outTs.append(outTi)
                        outTs.write()
                        outTsSet.update(outTs)
                        outTsSet.write()
                        self._store(outTsSet)
                else:
                    logger.error(f'tsId = {tsId} -> Output file {outputFn} was not generated. Skipping... ')
            except Exception as e:
                logger.error(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... ')

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append(f"Input tilt-series: {self.getInputSet().getSize()}\n"
                           "Transformation matrices calculated: "
                           f"{self.TiltSeries.getSize()}")

            interpTS = getattr(self, OUTPUT_TS_INTERPOLATED_NAME, None)
            if interpTS is not None:
                summary.append("Interpolated tilt-series: "
                               f"{interpTS.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.TiltSeries:
            methods.append(f"The transformation matrix has been calculated for "
                           f"{self.TiltSeries.getSize()} tilt-series using "
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
