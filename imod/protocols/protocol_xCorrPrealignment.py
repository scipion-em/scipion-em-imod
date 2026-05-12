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
import numpy as np
import pyworkflow.protocol.params as params
from imod.convert.convert import readXfFile
from imod.protocols.protocol_base import ProtImodBase
from imod.protocols.protocol_base_xcorr_fidmodel import ProtImodBaseXcorrFidModel
from pwem.objects import Transform
from pyworkflow.protocol import ProtStreamingBase, STEPS_PARALLEL
from pyworkflow.utils import Message, cyanStr, redStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from imod.constants import (TLT_EXT, PREXF_EXT, PREXG_EXT,
                            OUTPUT_TILTSERIES_NAME,
                            OUTPUT_TS_INTERPOLATED_NAME, XFTOXG_PROGRM, TILT_XCORR_PROGRAM)

logger = logging.getLogger(__name__)


class ProtImodXcorrPrealignment(ProtImodBase, ProtImodBaseXcorrFidModel, ProtStreamingBase):
    """
    Performs coarse prealignment of electron tomography tilt-series using
    cross-correlation methods from the IMOD package. The protocol estimates
    translational relationships between successive tilt images in order to
    establish an initial geometric alignment before more refined alignment
    or tomographic reconstruction procedures are applied. More info:
        https://bio3d.colorado.edu/imod/doc/man/tiltxcorr.html
        https://bio3d.colorado.edu/imod/doc/man/xftoxg.html

    AI Generated:

    Coarse Prealignment (ProtImodXcorrPrealignment) - User Manual
        Overview

        The Coarse Prealignment protocol provides an initial alignment of
        tilt-series images using IMOD cross-correlation procedures. Its
        primary objective is to estimate the relative translational shifts
        between neighboring tilt views so that the complete tilt-series can
        be brought into a consistent coordinate system before refinement and
        reconstruction. This stage is one of the fundamental early steps in
        electron tomography workflows because inaccurate initial alignment
        can negatively affect all downstream processing steps.

        In practical cryo-electron tomography workflows, images acquired at
        different tilt angles are affected by shifts introduced during data
        collection. These shifts arise from mechanical stage motion,
        specimen drift, beam-induced movement, or imperfections in tracking
        during acquisition. The protocol compensates for these effects by
        determining how consecutive tilt images relate spatially to one
        another.

        The alignment strategy is based on cross-correlation between
        neighboring tilt images. Because images acquired at different tilt
        angles are geometrically distorted relative to one another, the
        protocol accounts for angular stretching effects before estimating
        image displacements. This improves the robustness of the alignment,
        especially for high tilt angles where geometric distortions become
        more severe. More information:
        https://bio3d.colorado.edu/imod/doc/man/tiltxcorr.html

        Biological Context and Importance

        In electron tomography, the quality of the final tomogram strongly
        depends on the accuracy of the tilt-series alignment. Coarse
        prealignment serves as the first global correction stage and is
        typically followed by fiducial-based refinement or patch tracking
        refinement. Without a reliable initial alignment, later refinement
        stages may converge poorly or produce distorted reconstructions.

        For biological samples containing large cellular structures,
        membranes, organelles, or macromolecular assemblies, this protocol
        provides a robust starting point even when fiducial markers are not
        yet fully optimized. It is especially useful in workflows where
        large datasets must be processed efficiently before detailed local
        refinement.

        The protocol is also important for datasets acquired under low-dose
        conditions. In these situations, images often contain substantial
        noise, and establishing reliable image correspondence becomes more
        difficult. Cross-correlation prealignment improves image coherence
        and stabilizes subsequent processing steps.

        General Workflow

        The protocol receives one or more tilt-series as input and computes
        a consistent set of translational alignments across all tilt views.
        The resulting aligned geometry is stored as transformation matrices
        associated with each image in the tilt-series.

        The workflow is designed to support streaming acquisition
        environments. As new tilt-series become available during data
        collection, they can be processed progressively without waiting for
        the complete dataset to finish acquisition. This behavior is
        particularly useful in automated microscopy facilities and
        high-throughput cryo-electron tomography pipelines.

        The protocol preserves the original identity and metadata of the
        tilt-series while adding alignment information required for later
        tomographic reconstruction and refinement stages.

        Cumulative Correlation Strategy

        An optional cumulative correlation strategy can be used to improve
        alignment stability in challenging datasets. Instead of correlating
        each image only with its direct neighbor, the protocol progressively
        builds a reference from already aligned images. This approach can
        improve robustness when neighboring tilt images differ strongly due
        to large angular changes or low signal-to-noise ratios.

        From a biological perspective, cumulative correlation is especially
        useful for thick specimens, crowded cellular environments, or noisy
        cryo-focused ion beam lamellae where direct neighboring image
        correlation may become unstable. However, users should be aware that
        cumulative strategies may also propagate early alignment errors if
        the initial correlations are inaccurate.

        Tilt Axis Orientation

        Accurate tilt axis orientation is one of the most important physical
        parameters in electron tomography alignment. The protocol allows the
        tilt axis angle to be inherited from acquisition metadata or
        manually adjusted when necessary.

        Incorrect tilt axis orientation can introduce systematic alignment
        artifacts that become increasingly visible at high tilt angles.
        Biologically, such errors may distort membrane geometries, elongate
        macromolecular structures, or reduce interpretability of cellular
        organization in reconstructed tomograms.

        When acquisition metadata is uncertain or inconsistent, users should
        carefully verify the tilt axis orientation before proceeding to
        reconstruction. Even small angular inaccuracies may reduce the
        quality of the final tomogram.

        Filtering and Correlation Optimization

        The protocol provides filtering parameters that help optimize the
        cross-correlation procedure. These controls influence how image
        frequencies contribute to alignment estimation and can significantly
        affect robustness under different experimental conditions.

        For noisy datasets, stronger filtering may improve stability by
        suppressing high-frequency noise. For high-quality datasets,
        preserving more image detail may produce more precise alignment.
        Biological specimens with weak contrast, such as vitrified cellular
        material, often benefit from careful tuning of these filtering
        parameters.

        Excessive filtering, however, may remove biologically meaningful
        structural information and reduce alignment precision. In practice,
        moderate filtering combined with visual inspection usually provides
        the best balance.

        Trimming and Region Selection

        The protocol allows alignment to focus on selected image regions.
        This is particularly useful when portions of the image contain
        artifacts, contamination, empty ice, grid bars, or strongly varying
        background intensity.

        Restricting the correlation to biologically relevant regions can
        substantially improve alignment quality. For example, in cellular
        tomography it is often advantageous to exclude carbon edges or empty
        regions while focusing on the specimen-containing area.

        Proper region selection becomes especially important for large field
        of view acquisitions where only a subset of the image contains
        meaningful biological signal.

        Outputs and Interpretation

        The protocol produces tilt-series associated with transformation
        matrices representing the estimated translational alignment for each
        tilt image. These transformations establish a coherent geometric
        framework that can be used directly in subsequent reconstruction and
        refinement workflows.

        The aligned outputs remain linked to the original acquisition
        information while incorporating the newly estimated geometric
        corrections. This allows downstream protocols to interpret the data
        consistently throughout the tomography processing pipeline.

        Disabled or excluded images are preserved within the dataset to
        maintain acquisition consistency, although they are not used for the
        active alignment estimation.

        Streaming and High-Throughput Processing

        The protocol is designed for modern automated tomography pipelines
        where data may arrive continuously during acquisition. Streaming
        execution allows processing to begin immediately as tilt-series are
        collected, reducing delays between acquisition and reconstruction.

        This capability is especially valuable in facility environments,
        overnight automated collection sessions, or large screening
        experiments where rapid feedback on alignment quality is important.

        Practical Recommendations

        For most biological datasets, the default alignment settings provide
        a reliable starting point. Users should first evaluate the resulting
        image coherence visually before attempting more advanced refinement
        stages.

        If alignment appears unstable at high tilt angles, cumulative
        correlation and moderate filtering often improve robustness.
        Likewise, carefully restricting the correlation region to specimen
        areas may significantly improve results for heterogeneous or noisy
        datasets.

        When processing fiducial-based tomography experiments, this protocol
        is commonly used as a preparatory stage before bead tracking and
        fine alignment refinement. For fiducial-free workflows, the quality
        of the coarse alignment becomes even more important because later
        refinement relies heavily on the consistency established at this
        stage.

        Final Perspective

        Coarse prealignment is a foundational stage in electron tomography
        processing because it establishes the first coherent geometric model
        of the tilt-series. Although it is considered an initial alignment
        step, its quality strongly influences the success of later
        refinement, reconstruction, segmentation, and biological
        interpretation.

        Careful attention to tilt axis orientation, filtering behavior,
        image quality, and correlation regions is essential for obtaining
        reliable alignments that support accurate biological conclusions.
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
    def stepsGeneratorStep(self) -> None:
        self._initialize()
        closeSetStepDeps = []
        inTsSet = self.getInputTsSet()
        outTsSet = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        self.readingOutput(outTsSet)

        while True:
            with self._lock:
                inTsIds = set(inTsSet.getTSIds())

            if not inTsSet.isStreamOpen() and Counter(self.tsIdReadList) == Counter(inTsIds):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetsStep,
                                         OUTPUT_TILTSERIES_NAME,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break

            nonProcessedTsIds = inTsIds - set(self.tsIdReadList)
            tsToProcessDict = {tsId: ts.clone() for ts in inTsSet.iterItems()
                               if (tsId := ts.getTsId()) in nonProcessedTsIds  # Only not processed tsIds
                               and ts.getSize() > 0}  # Avoid processing empty TS
            for tsId, ts in tsToProcessDict.items():
                convId = self._insertFunctionStep(self.convertInStep,
                                                  ts,
                                                  prerequisites=[],
                                                  needsGPU=False)
                compId = self._insertFunctionStep(self.computeXcorrStep,
                                                  ts,
                                                  prerequisites=[convId],
                                                  needsGPU=False)
                outId = self._insertFunctionStep(self.createAliTsStep,
                                                 ts,
                                                 prerequisites=[compId],
                                                 needsGPU=False)
                closeSetStepDeps.append(outId)
                logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                self.tsIdReadList.append(tsId)

            self.refreshStreaming(inTsSet)

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

    def createAliTsStep(self, ts: TiltSeries):
        """ Generate tilt-series with the associated transform matrix """
        tsId = ts.getTsId()
        if tsId in self.failedItems:
            self.addToOutFailedSet(ts)
            return

        try:
            outputFn = self.getExtraOutFile(tsId, ext=PREXG_EXT)
            if exists(outputFn):
                self._registerOutput(ts, outputFn)
            else:
                logger.error(redStr(f'tsId = {tsId} -> Output file {outputFn} was not generated. Skipping... '))

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    @retry_on_sqlite_lock(log=logger)
    def _registerOutput(self, ts: TiltSeries, outputFn: str):
        tsId = ts.getTsId()
        with self._lock:
            tAx = self.getTiltAxisOrientation(ts)
            aliMatrixStack = readXfFile(outputFn)
            inTsSetPointer = self.getInputTsSet(pointer=True)
            # Set of tilt-series
            outTsSet = self.getOutputSetOfTS(inTsSetPointer,
                                             tiltAxisAngle=tAx)
            # Tilt-series
            outTs = TiltSeries()
            outTs.copyInfo(ts)
            outTs.getAcquisition().setTiltAxisAngle(self.getTiltAxisOrientation(ts))
            outTs.setAlignment2D()
            outTsSet.append(outTs)
            # Tilt-images
            stackIndex = 0
            for ti in ts.iterItems(orderBy=TiltImage.INDEX_FIELD):
                outTi = TiltImage()
                outTi.copyInfo(ti)
                if ti.isEnabled():
                    self.updateTiltImage(outTi, stackIndex, aliMatrixStack, tAx)
                    stackIndex += 1
                else:
                    self.updateDisabledTi(outTi)
                self.setTsOddEven(tsId, outTi, binGenerated=False)
                outTs.append(outTi)
            # Data persistence
            outTs.write()
            outTsSet.update(outTs)
            outTsSet.write()
            self._store(outTsSet)
            # Close explicitly the outputs (for streaming)
            self.closeOutputsForStreaming()

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
