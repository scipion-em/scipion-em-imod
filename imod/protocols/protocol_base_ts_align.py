# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
# *
# * [1] National Center of Biotechnology, CSIC, Spain
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
from os import stat
from os.path import exists
from typing import Tuple, List
import numpy as np
from imod.constants import XF_EXT
from imod.convert.convert import readXfFile
from imod.protocols import ProtImodBase
from imod.utils import formatAngleList
from pwem.objects import Transform
from pyworkflow.object import Pointer
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import TiltSeries, TiltImage

logger = logging.getLogger(__name__)


IDENTITY_MATRIX = np.eye(3)  # Store in memory instead of multiple creation


class ProtImodBaseTsAlign(ProtImodBase):
    """
    Provides a common framework for generating aligned tilt-series data
    from tomography alignment workflows. The protocol is designed to
    standardize how alignment results are interpreted, transferred, and
    preserved across different IMOD-based alignment strategies so that
    downstream tomographic reconstruction and analysis can proceed within
    a consistent geometric reference system.

    AI Generated:

    Base Tilt-Series Alignment Framework (ProtImodBaseTsAlign) — User Manual
        Overview

        The ProtImodBaseTsAlign protocol acts as a foundational alignment
        layer for tomography workflows that produce aligned tilt-series
        datasets. Its primary purpose is to ensure that image transformations,
        tilt geometry, and acquisition metadata are consistently propagated
        into a new aligned dataset that can be safely used for reconstruction,
        visualization, subtomogram averaging, or further processing.

        In cryo-electron tomography, tilt-series alignment is one of the
        most biologically important preprocessing stages because it defines
        the spatial consistency between all projection images acquired during
        the tilt experiment. Accurate alignment directly affects the quality
        of the final tomogram, the interpretability of structural features,
        and the reliability of any downstream quantitative analysis.

        Biological Context

        During tilt-series acquisition, unavoidable specimen drift, beam-
        induced motion, and mechanical inaccuracies introduce positional
        inconsistencies between projection images. Alignment protocols aim
        to compensate for these effects by estimating geometric corrections
        that place every tilt image into a coherent coordinate system.

        This framework is intended for workflows where alignment parameters
        have already been estimated by specialized algorithms such as fiducial
        tracking, patch tracking, or cross-correlation methods. The role of
        the protocol is to consolidate those alignment results into a properly
        formatted and internally consistent aligned tilt-series object.

        Preservation of Experimental Geometry

        One of the key objectives of the protocol is to preserve the biological
        and experimental meaning of the original acquisition. Each tilt image
        retains its associated tilt angle and acquisition identity while being
        updated with corrected geometric transformations. This allows the
        aligned dataset to remain fully traceable to the original experiment.

        From a practical cryo-ET perspective, maintaining accurate tilt-angle
        assignments is essential because tomographic reconstruction algorithms
        depend on the exact angular relationship between projections. Even
        small inconsistencies can produce reconstruction artifacts, elongation
        effects, or reduced interpretability in the final tomogram.

        Handling of Partial or Irregular Datasets

        The protocol is suitable for datasets where some projections may be
        excluded, disabled, or unavailable. In realistic experimental workflows,
        individual images are sometimes removed because of excessive motion,
        contamination, charging, or poor image quality. The framework preserves
        dataset continuity while maintaining valid transformations for the
        remaining projections.

        This capability is biologically important because excluding damaged
        projections often improves reconstruction quality more than attempting
        to retain every acquired image. The resulting aligned tilt-series
        therefore reflects a balance between geometric completeness and data
        quality.

        Transformation Management

        The protocol supports workflows where previous transformations may
        already exist. This is particularly relevant in iterative processing
        pipelines where multiple alignment refinements are performed over time.
        Instead of discarding earlier corrections, the framework preserves
        and combines transformation information so that geometric refinements
        accumulate consistently throughout processing.

        In biological applications, iterative refinement commonly occurs when
        users progressively improve alignment accuracy after fiducial cleaning,
        exposure filtering, or motion correction. Maintaining transformation
        continuity helps prevent inconsistencies between processing stages.

        Streaming and Large-Scale Processing

        The framework is compatible with streaming-oriented workflows in which
        aligned datasets are generated progressively while acquisition or
        preprocessing is still ongoing. This behavior is particularly valuable
        in high-throughput cryo-ET facilities where rapid feedback during data
        collection can improve microscope usage and experimental decision-making.

        For biological users working with large datasets, streaming support
        enables early inspection of alignment quality before the full experiment
        has completed. This can help identify problems such as fiducial loss,
        incorrect tilt-axis estimation, or stage instability at an early stage.

        Outputs and Their Interpretation

        The protocol produces aligned tilt-series datasets in which every tilt
        image is associated with updated geometric transformations and validated
        acquisition metadata. The resulting outputs are intended to serve as
        reliable inputs for tomographic reconstruction and all subsequent stages
        of structural interpretation.

        Biologically, the quality of the aligned tilt-series strongly influences
        the interpretability of reconstructed cellular environments, macromolecular
        assemblies, membrane organization, and subtomogram averaging results.
        Poor alignment propagates directly into reconstruction artifacts, whereas
        stable alignment improves resolution, contrast, and structural fidelity.

        Practical Recommendations

        In routine cryo-ET practice, users should carefully inspect aligned
        projections before reconstruction, paying particular attention to
        fiducial consistency, projection continuity, and residual drift.
        Datasets containing severe alignment instabilities should ideally be
        corrected before reconstruction rather than compensated for afterward.

        When iterative alignment refinement is performed, maintaining consistent
        geometry and metadata across all processing stages is essential for
        reliable biological interpretation. Users should also verify that tilt
        angles remain physically meaningful and correspond correctly to the
        acquisition geometry.

        Final Perspective

        For most cryo-electron tomography workflows, alignment is not simply
        a technical correction step but the process that establishes the spatial
        foundation of the entire experiment. Reliable management of geometric
        transformations, tilt angles, and acquisition consistency is therefore
        essential for producing biologically interpretable tomograms and robust
        downstream structural analyses.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- STEPS functions -----------------------------

    # --------------------------- UTILS functions -----------------------------
    def getTltFilePath(self, tsId: str):
        # To be defined by the child classes
        pass

    @retry_on_sqlite_lock(log=logger)
    def createOutTs(self,
                    ts: TiltSeries,
                    inTsSetPointer: Pointer) -> None:
        tsId = ts.getTsId()
        xfFile = self.getExtraOutFile(tsId, suffix="fid", ext=XF_EXT)
        if exists(xfFile) and stat(xfFile).st_size != 0:
            tltFile = self.getTltFilePath(tsId)
            aliMatrix = readXfFile(xfFile)
            tiltAngles = formatAngleList(tltFile)
            # Set of tilt-series
            with self._lock:
                outTsSet = self.getOutputSetOfTS(inTsSetPointer)
                # Tilt-series
                outTs = TiltSeries()
                outTs.copyInfo(ts)
                outTs.setAlignment2D()
                outTsSet.append(outTs)
                # Tilt-images
                stackIndex = 0
                for ti in ts.iterItems(orderBy=TiltImage.INDEX_FIELD):
                    outTi = TiltImage()
                    outTi.copyInfo(ti)
                    if ti.isEnabled():
                        tiltAngle, newTransformArray = self._getTrDataEnabled(stackIndex,
                                                                              aliMatrix,
                                                                              tiltAngles)
                        stackIndex += 1
                    else:
                        tiltAngle, newTransformArray = self._getTrDataDisabled(ti)
                    self._updateTiltImage(ti, outTi, newTransformArray, tiltAngle)
                    self.setTsOddEven(tsId, outTi, binGenerated=False)
                    outTs.append(outTi)
                # Data persistence
                outTs.write()
                outTsSet.update(outTs)
                outTsSet.write()
                self._store(outTsSet)
                # Close explicitly the outputs (for streaming)
                self.closeOutputsForStreaming()
        else:
            logger.error(f'tsId = {tsId} -> Output file {xfFile} was not generated or is empty. Skipping... ')

    @staticmethod
    def _getTrDataEnabled(stackIndex: int,
                          alignmentMatrix: np.ndarray,
                          tiltAngleList: List[float]) -> Tuple[float, np.ndarray]:
        newTransform = alignmentMatrix[:, :, stackIndex]
        newTransformArray = np.array(newTransform)
        tiltAngle = float(tiltAngleList[stackIndex])
        return tiltAngle, newTransformArray

    @staticmethod
    def _getTrDataDisabled(ti: TiltImage) -> Tuple[float, np.ndarray]:
        tiltAngle = ti.getTiltAngle()
        if ti.hasTransform():
            newTransformArray = ti.getTransform().getMatrix()
        else:
            newTransformArray = IDENTITY_MATRIX
        return tiltAngle, newTransformArray

    @staticmethod
    def _updateTiltImage(ti: TiltImage,
                         outTi: TiltImage,
                         newTransformArray: np.ndarray,
                         tiltAngle: float) -> None:
        transform = Transform()
        if ti.hasTransform():
            previousTransform = ti.getTransform().getMatrix()
            previousTransformArray = np.array(previousTransform)
            outputTransformMatrix = np.matmul(newTransformArray, previousTransformArray)
            transform.setMatrix(outputTransformMatrix)
        else:
            transform.setMatrix(newTransformArray)

        outTi.setTransform(transform)
        outTi.setTiltAngle(tiltAngle)
