# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * National Center of Biotechnology, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
# **************************************************************************
import logging
import subprocess
import traceback
import typing
from collections import Counter
from imod import Plugin
from imod.constants import OUTPUT_TILTSERIES_NAME, TLT_EXT, PATCH_TRACKING, FIDUCIAL_MODEL, \
    BRT_ENV_NAME
from imod.protocols.protocol_base_ts_align import ProtImodBaseTsAlign
from imod.protocols.protocol_fiducialAlignment import TILT_ALIGN_PROGRAM
from pyworkflow.constants import BETA
from pyworkflow.protocol import EnumParam, IntParam, GT, STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.utils import Message, cyanStr, redStr
from tomo.objects import SetOfTiltSeries, TiltSeries

logger = logging.getLogger(__name__)

class ProtImodBRT(ProtImodBaseTsAlign, ProtStreamingBase):
    """
    Performs automatic alignment of tilt-series datasets using the IMOD
    batchruntomo (https://bio3d.colorado.edu/imod/doc/man/batchruntomo.html) wrapper made
    by Team Tomo (yet-another-imod-wrapper https://teamtomo.org/teamtomo-site-archive/)
    workflow and TeamTomo integration utilities.

    AI Generated:

    BatchRunTomo Alignment (ProtImodBRT) — User Manual
        Overview

        The BatchRunTomo Alignment protocol performs automatic alignment of
        cryo-electron tomography tilt-series using the IMOD batchruntomo
        framework integrated through TeamTomo utilities. Its main objective is
        to generate geometrically aligned tilt-series suitable for tomographic
        reconstruction and downstream structural analysis.

        In cryo-electron tomography workflows, tilt-series alignment is one of
        the most critical preprocessing stages because the quality of the final
        tomogram depends strongly on the accuracy of the geometric corrections.
        Misalignment between projection images introduces blurring, artifacts,
        and loss of structural detail that can compromise biological
        interpretation. This protocol automates the alignment process while
        remaining compatible with high-throughput and streaming acquisition
        environments.

        The protocol is designed for both routine facility processing and more
        demanding biological studies involving cellular tomography, in situ
        structural biology, membrane systems, or macromolecular complexes.

        Inputs and General Workflow

        The protocol requires one or more tilt-series datasets together with
        their associated tilt geometry information. Each tilt-series is
        processed independently and produces an aligned version that can later
        be reconstructed into a tomogram.

        The alignment procedure is intended to compensate for shifts,
        distortions, and geometric inconsistencies introduced during data
        acquisition. Correct alignment ensures that corresponding structural
        features remain spatially consistent across projection views, allowing
        reliable three-dimensional reconstruction.

        Because the protocol supports streaming operation, it is especially
        suitable for modern automated microscopy pipelines where datasets are
        processed progressively as they become available. This capability is
        highly valuable in large cryo-electron tomography facilities and
        high-throughput screening workflows.

        Alignment Strategies

        The protocol provides two biologically relevant alignment strategies:
        fiducial-based alignment and patch-tracking alignment. The appropriate
        choice depends primarily on sample preparation quality and the presence
        of fiducial markers.

        Fiducial alignment is generally considered the preferred approach when
        fiducial gold beads are present and clearly visible throughout the
        tilt-series. In this strategy, the alignment relies on the stable
        tracking of fiducial markers across projection angles. This often
        yields highly accurate geometric correction and is commonly used in
        high-resolution tomography workflows.

        From a biological perspective, fiducial alignment is especially useful
        for purified complexes, lamellae, and cellular preparations where
        fiducials are evenly distributed and remain visible across the full
        angular range. Choosing an appropriate fiducial diameter is important
        because it should approximately reflect the physical size of the gold
        markers used during specimen preparation.

        Patch tracking provides an alternative strategy for datasets lacking
        fiducials or containing poorly visible markers. Instead of relying on
        external reference particles, the alignment is derived from intrinsic
        image features distributed throughout the specimen.

        This method is often advantageous for dense cellular samples, crowded
        intracellular environments, or specimens where fiducials interfere with
        biological interpretation. Patch tracking can also simplify workflows
        by eliminating the need for fiducial deposition during sample
        preparation.

        Patch Size and Overlap Considerations

        When patch tracking is selected, the protocol allows control over the
        patch size and overlap between neighboring regions. These parameters
        influence both alignment robustness and computational cost.

        Larger patches generally provide more stable tracking in noisy datasets
        because they contain more structural information. However, excessively
        large patches may average together regions with different local motion
        or deformation. Smaller patches can better capture local image
        behavior, although they may become unstable in low-contrast regions.

        The overlap percentage determines how strongly neighboring tracking
        regions interact spatially. Higher overlap values usually improve
        continuity and robustness at the expense of increased computational
        effort. In challenging biological datasets containing heterogeneous
        cellular structures, moderate to high overlap often produces more
        reliable alignment behavior.

        Streaming and High-Throughput Processing

        One of the major strengths of this protocol is its compatibility with
        streaming data acquisition. New tilt-series can be aligned
        automatically as they are produced, enabling near real-time processing
        during microscope operation.

        This capability is particularly valuable in modern cryo-electron
        tomography facilities, where immediate feedback can help identify
        acquisition problems early in an experiment. Researchers can therefore
        evaluate alignment quality, specimen behavior, or imaging conditions
        without waiting for the completion of an entire data collection
        session.

        Outputs and Their Interpretation

        The protocol produces aligned tilt-series datasets that preserve the
        original acquisition identity while incorporating the computed
        geometric corrections. These aligned datasets are intended for direct
        use in tomographic reconstruction workflows and downstream structural
        analysis.

        Accurate alignment improves tomogram quality by reducing reconstruction
        artifacts and preserving structural consistency across projection
        angles. This directly affects the interpretability of biological
        features such as membranes, cytoskeletal assemblies, ribosomes, viral
        particles, or macromolecular complexes embedded within cells.

        Practical Recommendations

        In routine biological workflows, fiducial alignment is generally the
        preferred starting point whenever high-quality fiducial markers are
        available. Proper fiducial distribution across the field of view often
        provides the most robust and accurate geometric correction.

        Patch tracking is highly valuable for fiducial-free workflows or
        difficult cellular specimens where fiducials are sparse, obscured, or
        absent. In these cases, adjusting patch size and overlap can
        significantly improve alignment stability.

        Researchers should visually inspect alignment quality before proceeding
        to tomographic reconstruction, particularly in datasets with strong
        specimen deformation, low contrast, or limited angular coverage.

        Final Perspective

        For most cryo-electron tomography studies, tilt-series alignment is a
        foundational step that strongly influences the quality and biological
        reliability of the final reconstruction. Careful selection between
        fiducial-based and patch-based alignment strategies allows researchers
        to adapt the workflow to the characteristics of the specimen, imaging
        conditions, and experimental objectives while maintaining compatibility
        with modern automated tomography pipelines.
    """
    """Automatic tilt-series alignment using IMOD's batchruntomo.
    """

    _label = "teamtomo/batchruntomo"
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}
    stepsExecutionMode = STEPS_PARALLEL
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @classmethod
    def worksInStreaming(cls):
        return True

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        super().addInTsSetFormParam(form)
        form.addParam('alignMode', EnumParam,
                      important=True,
                      choices=['Fiducial alignment', 'Patch tracking'],
                      default=FIDUCIAL_MODEL,
                      label='Alignment mode')
        form.addParam('fidSize', IntParam,
                      condition=f'alignMode == {FIDUCIAL_MODEL}',
                      default=10,
                      validators=[GT(0)],
                      label='Fiducial diameter (nm)')
        patchTrackingCond = f'alignMode == {PATCH_TRACKING}'
        form.addParam('patchSize', IntParam,
                      condition=patchTrackingCond,
                      default=500,
                      validators=[GT(0)],
                      label='Patch side-length (A)')
        form.addParam('patchOverlapPercent', IntParam,
                      condition=patchTrackingCond,
                      default=80,
                      validators=[GT(0)],
                      label='Patch overlap percent',
                      help='Percentage of tile-length to overlap on each side.')
        form.addParallelSection(threads=3, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def stepsGeneratorStep(self) -> None:
        """
        This step should be implemented by any streaming protocol.
        It should check its input and when ready conditions are met
        call the self._insertFunctionStep method.
        """
        closeSetStepDeps = []
        inTsSet = self.getInputTsSet()
        self.readingOutput(getattr(self, OUTPUT_TILTSERIES_NAME, None))

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
                cInputId = self._insertFunctionStep(self.convertInStep, ts,
                                                    prerequisites=[],
                                                    needsGPU=False)
                predFidId = self._insertFunctionStep(self.runBRT, ts,
                                                     prerequisites=cInputId,
                                                     needsGPU=False)
                cOutId = self._insertFunctionStep(self.createOutputStep, ts,
                                                  prerequisites=predFidId,
                                                  needsGPU=False)
                closeSetStepDeps.append(cOutId)
                logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                self.tsIdReadList.append(tsId)

            self.refreshStreaming(inTsSet)

    # --------------------------- STEPS functions -----------------------------
    def convertInStep(self, ts: TiltSeries):
        presentAcqOrders = ts.getTsPresentAcqOrders()
        super().convertInputStep(ts, presentAcqOrders=presentAcqOrders)

    def runBRT(self, ts: TiltSeries):
        tsId = ts.getTsId()
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'tsId = {tsId}: aligning...'))
                self.genTsPaths(tsId)
                args = self._getFiducialAliCmd(ts) if self.alignMode.get() == FIDUCIAL_MODEL \
                    else self._getPatchTrackingCmd(ts)
                Plugin.runBRT(self, args)
            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> {TILT_ALIGN_PROGRAM} execution failed with the exception -> {e}'))
                logger.error(traceback.format_exc())

    def createOutputStep(self, ts: TiltSeries):
        tsId = ts.getTsId()
        if tsId in self.failedItems:
            self.addToOutFailedSet(ts)
            return
        try:
            self.createOutTs(ts, self.getInputTsSet(pointer=True))

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    # --------------------------- INFO functions ------------------------------

    # --------------------------- UTILS functions -----------------------------
    def _getCommonCmd(self, ts: TiltSeries):
        tsId = ts.getTsId()
        cmd = [
            f'--tilt-series {self.getTmpOutFile(tsId)}',
            f'--tilt-angles {self.getExtraOutFile(tsId, ext=TLT_EXT)}',
            f'--output-directory {self._getExtraPath(tsId)}',
            f'--pixel-size {ts.getSamplingRate():.3f}',
            f'--nominal-rotation-angle {ts.getAcquisition().getTiltAxisAngle():.2f}',
            f'--basename {tsId}'
        ]
        return ' '.join(cmd)

    def _getFiducialAliCmd(self, ts: TiltSeries):
        cmd = [
            'fiducials',
            self._getCommonCmd(ts),
            f'--fiducial-size {self.fidSize.get()}',
        ]
        return ' '.join(cmd)

    def _getPatchTrackingCmd(self, ts: TiltSeries):
        cmd = [
            'patch-tracking',
            self._getCommonCmd(ts),
            f'--patch-size {self.patchSize.get()}',
            f'--patch-overlap-percentage {self.patchOverlapPercent.get()}'
        ]
        return ' '.join(cmd)

    def getTltFilePath(self, tsId):
        return self.getExtraOutFile(tsId, suffix="fid", ext=TLT_EXT)
