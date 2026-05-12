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
import os
import sqlite3
import traceback
from collections import Counter
from os.path import exists
from typing import Union
import pyworkflow.protocol.params as params
from imod.convert.convert import fiducialModel2List, fidResidualModel2List
from imod.protocols.protocol_base_ts_align import ProtImodBaseTsAlign
from pyworkflow.object import Pointer
from pyworkflow.protocol import STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.utils import Message, cyanStr, redStr, yellowStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import (LandmarkModel, SetOfLandmarkModels, SetOfTiltSeries,
                          TiltSeries)
from imod.constants import (TLT_EXT, XF_EXT, FID_EXT, TXT_EXT, XYZ_EXT,
                            MOD_EXT, SFID_EXT, OUTPUT_TILTSERIES_NAME,
                            OUTPUT_FIDUCIAL_NO_GAPS_NAME,
                            OUTPUT_TS_INTERPOLATED_NAME,
                            OUTPUT_TS_COORDINATES_NAME, ALIGNLOG_PROGRAM, TILT_ALIGN_PROGRAM, MODEL2POINT_PROGRAM,
                            POINT2MODEL_PROGRAM)

logger = logging.getLogger(__name__)

# Rotation solution types
NO_ROTATION = 0
ONE_ROTATION = 1
GROUP_ROTATIONS = 2
ALL_ROTATIONS = 3
ROT_SOLUTION_CHOICES = ['No rotation',
                        'One rotation',
                        'Group rotations',
                        'Solve for all rotations']

# Magnification solution types
FIXED_MAG = 0
GROUP_MAGS = 1
ALL_MAGS = 2
MAG_SOLUTION_CHOICES = ['Fixed magnification at 1.0',
                        'Group magnifications',
                        'Solve for all magnifications']

# Tilt angle solution types
FIXED_TILT = 0
GROUP_TILTS = 1
ALL_EXCEPT_MIN = 2
TILT_SOLUTION_CHOICES = ['Fixed tilt angles',
                         'Group tilt angles',
                         'Solve for all except minimum tilt']

# Distortion solution types
DIST_DISABLED = 0
DIST_FULL_SOLUTION = 1
DIST_SKEW_ONLY = 2
DISTORTION_SOLUTION_CHOICES = ['Disabled',
                               'Full solution',
                               'Skew only']


class ProtImodFiducialAlignment(ProtImodBaseTsAlign, ProtStreamingBase):
    """
    Constructs fiducial models and aligns electron tomography tilt-series
    using the IMOD fiducial alignment workflow. The protocol estimates
    geometric transformations between tilted projections so that all views
    are brought into a consistent tomographic coordinate system suitable
    for high-quality 3D reconstruction and downstream structural analysis.
    More info: https://bio3d.colorado.edu/imod/doc/man/tiltalign.html

    AI Generated:

    Fiducial Alignment (ProtImodFiducialAlignment) - User Manual
        Overview

        The Fiducial Alignment protocol performs alignment of tilt-series
        using fiducial markers detected across multiple projection images.
        In electron tomography workflows, this step is essential because it
        determines the geometric relationship between all tilted views before
        tomographic reconstruction. Accurate alignment directly influences the
        quality, interpretability, and resolution of the final 3D volume.

        The protocol is based on the IMOD tiltalign procedure, a widely used
        approach in cryo-electron tomography and cellular electron microscopy.
        It estimates rotations, translations, tilt geometry, magnification
        variations, and optional distortions affecting the acquisition process.
        The resulting alignment model provides a mathematically consistent
        description of how each projection image relates to the specimen and
        microscope geometry.

        Biological and Experimental Context

        Fiducial markers are typically gold beads or similar high-contrast
        particles distributed across the sample before data acquisition.
        Because these markers remain visible through multiple tilt angles,
        they provide reliable reference points for estimating specimen motion
        and imaging distortions during tomography acquisition.

        In practical biological applications, fiducial alignment is required
        before reconstructing cellular environments, organelles, membrane
        systems, macromolecular complexes, or in situ structures. The protocol
        is particularly valuable in cryo-electron tomography workflows where
        accurate geometric correction is necessary to preserve structural
        detail and avoid reconstruction artifacts.

        The protocol supports datasets where fiducial markers are not visible
        in every projection. This flexibility is biologically important because
        fiducials may disappear at high tilts, become occluded by dense sample
        regions, or lose contrast under challenging imaging conditions.

        Inputs and Workflow

        The protocol requires a set of fiducial landmark models associated
        with one or more tilt-series. Each fiducial corresponds to a tracked
        marker observed across several projection images. The alignment process
        uses these trajectories to estimate the global imaging geometry and
        correct positional inconsistencies between views.

        In streaming workflows, tilt-series can be processed progressively as
        new data become available. This is especially useful in automated
        microscopy facilities or high-throughput tomography pipelines where
        acquisition and processing occur simultaneously.

        The protocol automatically prepares fiducial information for IMOD
        processing, computes alignment solutions, and generates updated
        landmark models together with transformed tilt-series geometry. The
        outputs are suitable for subsequent tomographic reconstruction and
        visualization workflows.

        Fiducial Geometry and Surface Modeling

        The protocol supports both single-surface and dual-surface fiducial
        configurations. In many tomography experiments, fiducial beads may be
        deposited on one side of the specimen support film, while in other
        preparations they may exist on both surfaces. Correctly specifying
        this configuration improves geometric stability and alignment accuracy.

        For thicker biological samples or lamella preparations, dual-surface
        configurations can substantially improve robustness because fiducials
        distributed at different depths provide stronger geometric constraints.
        However, the selected option should remain consistent with the original
        fiducial detection strategy to avoid inconsistencies in the alignment
        model.

        Rotation, Tilt, and Magnification Modeling

        The protocol provides several strategies for solving rotational,
        tilt-angle, and magnification parameters. Users may completely fix
        these variables, solve them independently for every projection, or
        group neighboring views together to stabilize the estimation process.

        Grouped solutions are often biologically advantageous because nearby
        tilt images usually share similar acquisition geometry. Reducing the
        number of free variables can improve robustness, particularly in noisy
        datasets or when fiducial coverage is incomplete.

        Solving all parameters independently provides maximum flexibility and
        may be useful for difficult datasets with strong geometric variability.
        However, unconstrained solutions may also become unstable if the number
        of fiducials is limited or unevenly distributed.

        Distortion and Stretch Correction

        In addition to rigid alignment, the protocol can estimate image
        distortions such as stretching and skew. These corrections compensate
        for imperfections in the imaging system, specimen deformation, or
        projection geometry.

        Distortion correction becomes especially important in large fields of
        view, thick cellular specimens, or datasets collected under demanding
        acquisition conditions. Correcting these effects can significantly
        improve reconstruction fidelity and reduce systematic geometric errors.

        Depending on the dataset quality and microscope stability, users may
        choose to disable distortion correction entirely or estimate grouped
        distortion parameters to balance flexibility and robustness.

        Robust Fitting and Error Handling

        The protocol incorporates robust fitting strategies that reduce the
        influence of unreliable fiducial measurements. Fiducials with unusually
        large residual errors may be down-weighted or excluded from the final
        alignment solution.

        This behavior is particularly valuable in biological datasets where
        fiducials can occasionally be misidentified, partially obscured, or
        affected by contamination and radiation damage. Robust fitting improves
        overall stability and helps prevent isolated tracking errors from
        degrading the complete alignment.

        Residual information generated during processing can also help users
        identify problematic regions of the dataset, evaluate alignment quality,
        and decide whether additional fiducial curation is necessary.

        Outputs and Interpretation

        The protocol produces aligned tilt-series geometry together with updated
        fiducial models containing refined landmark positions and residual
        information. These outputs define the corrected spatial relationship
        between all projections and are directly usable for tomographic
        reconstruction procedures.

        Gap-filled fiducial models may also be generated. These models provide
        continuous fiducial trajectories even when some markers are absent in
        particular views, improving compatibility with downstream tomography
        software and reconstruction algorithms.

        From a biological perspective, high-quality alignment is essential for
        preserving fine structural details in reconstructed tomograms. Poor
        alignment can lead to blurred membranes, distorted macromolecular
        assemblies, and loss of interpretability in cellular environments.

        Practical Recommendations

        For most cryo-electron tomography experiments, reliable fiducial
        detection and careful marker distribution are the most important
        prerequisites for successful alignment. Fiducials should ideally be
        well distributed across the field of view and visible through a broad
        tilt range.

        Grouped parameter estimation is often a good starting point because it
        balances robustness and flexibility. More aggressive parameter solving
        strategies should generally be reserved for datasets with strong
        geometric variability or known acquisition inconsistencies.

        When processing thick specimens or dual-surface fiducial preparations,
        enabling the two-surface option can substantially improve alignment
        quality. Users should also inspect residual statistics carefully, since
        unusually high residuals may indicate tracking problems or distorted
        acquisition geometry.

        Final Perspective

        Fiducial alignment is one of the central stages of electron tomography
        processing because it establishes the geometric foundation for all
        downstream reconstruction and interpretation. Accurate alignment not
        only improves visual reconstruction quality but also enhances the
        reliability of quantitative biological conclusions derived from
        tomographic data.

        Careful fiducial tracking, appropriate geometric modeling, and
        thoughtful parameter selection are therefore essential for obtaining
        biologically meaningful tomograms suitable for structural and cellular
        analysis.
    """

    _label = 'Fiducial alignment'
    _possibleOutputs = {
        OUTPUT_TILTSERIES_NAME: SetOfTiltSeries,
        OUTPUT_FIDUCIAL_NO_GAPS_NAME: SetOfLandmarkModels
    }
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @classmethod
    def worksInStreaming(cls):
        return True

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)

        form.addParam('inputSetOfLandmarkModels',
                      params.PointerParam,
                      pointerClass='SetOfLandmarkModels',
                      important=True,
                      label='Fiducial model')

        form.addParam('twoSurfaces',
                      params.BooleanParam,
                      default=False,
                      label='Assume beads on two surfaces?',
                      help="Track fiducials differentiating in which side of the sample are located.\n"
                           "IMPORTANT: It is highly recommended to match the option selected in the "
                           "generation of the fiducial models. In case they  do not match, it is not "
                           "intended to fail but could be missing the whole potential of the algorithm. "
                           "In case the algorithm used for the calculation of the fiducial models does "
                           "not consider this option it is algo recommended to set this "
                           "option to 'No'.")

        form.addSection('Global variables')

        form.addParam('rotationSolutionType',
                      params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=ROT_SOLUTION_CHOICES,
                      default=GROUP_ROTATIONS,
                      label='Rotation solution type',
                      help='Type of rotation solution: See rotOption in tiltalign IMOD command \n'
                           '* No rotation: The in-plane rotation will not be estimated\n'
                           '* One rotation: To solve for a single rotation variable \n'
                           '* Group rotations: Group views to solve for fewer rotations variables. Automapping of '
                           'rotation variables linearly changing values\n'
                           '* Solve for all rotations: for each view having an independent rotation\n')

        form.addParam('groupRotationSize',
                      params.IntParam,
                      default=5,
                      condition='rotationSolutionType == %i' % GROUP_ROTATIONS,
                      label='Group size',
                      help='Default group size when automapping rotation variables')

        form.addParam('magnificationSolutionType',
                      params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=MAG_SOLUTION_CHOICES,
                      default=FIXED_MAG,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Magnification solution type',
                      help='Type of magnification solution: See MagOption in tiltaling IMOD command\n'
                           '* Fixed magnification: Do not solve magnification. This fixes all magnifications at 1.0.\n'
                           '* Group magnifications: Group views to solve for fewer magnifications variables. '
                           'Automapping of variables (linearly changing values)  \n'
                           '* Solve for all magnifications: to vary all magnifications  of each view independently\n')

        form.addParam('groupMagnificationSize',
                      params.IntParam,
                      default=4,
                      condition='magnificationSolutionType == %i' % GROUP_MAGS,
                      label='Group size',
                      help='Group size when automapping magnification variables')

        form.addParam('tiltAngleSolutionType',
                      params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=TILT_SOLUTION_CHOICES,
                      default=GROUP_TILTS,
                      label='Tilt angle solution type',
                      help='Type of tilt angle solution: See TiltOption in tiltalign IMOD command\n'
                           ' * Fixed tilt angles: To fix all tilt angles at their initial (input) values \n'
                           ' * Group tilt angles: To automap groups of tilt angles (linearly changing values) \n'
                           ' * Solve for all except minimum tilt:to solve for all tilt angles except for the view '
                           'at minimum tilt \n')

        form.addParam('groupTiltAngleSize',
                      params.IntParam,
                      default=5,
                      condition='tiltAngleSolutionType == %i' % GROUP_TILTS,
                      label='Group size',
                      help='Average default group size when automapping tilt variables')

        form.addParam('distortionSolutionType',
                      params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=DISTORTION_SOLUTION_CHOICES,
                      default=DIST_DISABLED,
                      label='Distortion solution type',
                      help='Type of skew solution:'
                           '* 0 to fix all skew angles at 0.0 \n'
                           '* 1 to vary all skew angles independently\n '
                           '* 2 to specify a mapping of skew variables, or\n '
                           '* 3 or 4 for automapping of variables (3 for linearly changing values or 4 for values all '
                           'the same within a group)..')

        form.addParam('xStretchGroupSize',
                      params.IntParam,
                      default=7,
                      condition='distortionSolutionType == %i' % DIST_FULL_SOLUTION,
                      label='X stretch group size',
                      help='Basic grouping size for X stretch')

        form.addParam('skewGroupSize',
                      params.IntParam,
                      default=11,
                      condition='distortionSolutionType in [%i, %i]' % (DIST_FULL_SOLUTION, DIST_SKEW_ONLY),
                      label='Skew group size',
                      help='Size of the skew group')
        form.addParallelSection(threads=3, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def stepsGeneratorStep(self) -> None:
        closeSetStepDeps = []
        inSetOfLandmarks = self.getInputSetOfLandmarks()
        inTsSet = self._getInTsSet()
        outTsSet = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        self.readingOutput(outTsSet)

        while True:
            with self._lock:
                inTsIds = set(inTsSet.getTSIds())

            if not inSetOfLandmarks.isStreamOpen() and Counter(self.tsIdReadList) == Counter(inTsIds):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetsStep,
                                         [OUTPUT_TILTSERIES_NAME, OUTPUT_FIDUCIAL_NO_GAPS_NAME],
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break

            nonProcessedTsIds = inTsIds - set(self.tsIdReadList)
            landmarkModelToProcessDict = {tsId: lMk.clone() for lMk in inSetOfLandmarks.iterItems()
                                          if (tsId := lMk.getTsId()) in nonProcessedTsIds  # Only not processed tsIds
                                          and lMk.getSize() > 0}  # Avoid processing empty landmark models
            tsToProcessDict = {tsId: ts.clone() for ts in inTsSet.iterItems()
                               if (tsId := ts.getTsId()) in nonProcessedTsIds  # Only not processed tsIds
                               and ts.getSize() > 0}  # Avoid processing empty TS
            for tsId, lMk in landmarkModelToProcessDict.items():
                ts = tsToProcessDict.get(tsId, None)
                if not ts:
                    logger.info(yellowStr(f'tsId = {tsId} - no corresponding TS to the landmark model was found...'))
                    continue
                cInId = self._insertFunctionStep(self.convertInStep,
                                                 lMk,
                                                 ts,
                                                 prerequisites=[],
                                                 needsGPU=False)
                fidAliId = self._insertFunctionStep(self.computeFiducialAlignmentStep,
                                                    lMk,
                                                    ts,
                                                    prerequisites=cInId,
                                                    needsGPU=False)
                p2mId = self._insertFunctionStep(self.translateFiducialPointModelStep,
                                                 tsId,
                                                 prerequisites=fidAliId,
                                                 needsGPU=False)
                cOutId = self._insertFunctionStep(self.createOutputStep,
                                                  lMk,
                                                  ts,
                                                  prerequisites=p2mId,
                                                  needsGPU=False)
                closeSetStepDeps.append(cOutId)
                logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                self.tsIdReadList.append(tsId)

            self.refreshStreaming(inSetOfLandmarks)

    # --------------------------- STEPS functions -----------------------------
    def convertInStep(self, lMk: LandmarkModel, ts: TiltSeries):
        tsId = lMk.getTsId()
        self.convertInputStep(ts, ts.getTsPresentAcqOrders())
        if tsId in self.failedItems:
            return
        try:
            # TODO: Check the non-imod fiducials combined with excluded views
            # If the fiducial model was not computed with imod, it is converted into the expected format here
            currentModelFn = lMk.getModelName()
            if not currentModelFn or not exists(currentModelFn):
                logger.info(cyanStr(f'tsId = {tsId}: not IMOD fiducial model detected. Generating...'))
                fidFn = self._getExtraPath(tsId, tsId + '_imod.fid')
                tsFn = self.getTmpOutFile(tsId)
                fnTxt = self._getExtraPath(tsId, tsId + '_points.txt')
                self._convertTxt2Fid(currentModelFn, fnTxt)
                paramsPoint2Model = {
                    "-InputFile": fnTxt,
                    "-OutputFile": fidFn,
                    "-image": tsFn
                }
                self.runProgram(POINT2MODEL_PROGRAM, paramsPoint2Model)
        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {POINT2MODEL_PROGRAM} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def computeFiducialAlignmentStep(self, lMk: LandmarkModel, ts: TiltSeries):
        tsId = lMk.getTsId()
        if tsId in self.failedItems:
            return
        try:
            logger.info(cyanStr(f'tsId = {tsId}: aligning...'))
            currentModelFn = lMk.getModelName()
            tsFn = self.getTmpOutFile(tsId)
            paramsTiltAlign = {"-ModelFile": currentModelFn, "-ImageFile": tsFn, "-ImagesAreBinned": 1,
                               "-UnbinnedPixelSize": ts.getSamplingRate() / 10,
                               "-OutputModelFile": self.getExtraOutFile(tsId, suffix="fidxyz", ext=MOD_EXT),
                               "-OutputResidualFile": self.getExtraOutFile(tsId, suffix="resid", ext=TXT_EXT),
                               "-OutputFidXYZFile": self.getExtraOutFile(tsId, suffix="fid", ext=XYZ_EXT),
                               "-OutputTiltFile": self.getExtraOutFile(tsId, suffix="interpolated", ext=TLT_EXT),
                               "-OutputXAxisTiltFile": self.getExtraOutFile(tsId, ext="xtilt"),
                               "-OutputTransformFile": self.getExtraOutFile(tsId, suffix="fid", ext=XF_EXT),
                               "-OutputFilledInModel": self.getExtraOutFile(tsId, suffix="noGaps", ext=FID_EXT),
                               "-RotationAngle": ts.getAcquisition().getTiltAxisAngle(),
                               "-TiltFile": self.getExtraOutFile(tsId, ext=TLT_EXT), "-AngleOffset": 0.0,
                               "-RotOption": self.getRotationType(),
                               "-RotDefaultGrouping": self.groupRotationSize.get(),
                               "-TiltOption": self.getTiltAngleType(),
                               "-TiltDefaultGrouping": self.groupTiltAngleSize.get(), "-MagReferenceView": 1,
                               "-MagOption": self.getMagnificationType(),
                               "-MagDefaultGrouping": self.groupMagnificationSize.get(),
                               "-XStretchOption": self.getStretchType(), "-SkewOption": self.getSkewType(),
                               "-XStretchDefaultGrouping": self.xStretchGroupSize.get(),
                               "-SkewDefaultGrouping": self.skewGroupSize.get(), "-BeamTiltOption": 0,
                               "-XTiltOption": 0, "-XTiltDefaultGrouping": 2000, "-ResidualReportCriterion": 3.0,
                               "-SurfacesToAnalyze": self.getSurfaceToAnalyze(), "-MetroFactor": 0.25,
                               "-MaximumCycles": 1000, "-KFactorScaling": 1.0, "-NoSeparateTiltGroups": 1,
                               "-AxisZShift": 0.0, "-ShiftZFromOriginal": 1, "-TargetPatchSizeXandY": '700,700',
                               "-MinSizeOrOverlapXandY": '0.5,0.5', "-MinFidsTotalAndEachSurface": '8,3',
                               "-FixXYZCoordinates": 0, "-RobustFitting": "",
                               "2>&1 | tee ": self._getExtraPath(tsId, "align.log")}

            self.runProgram(TILT_ALIGN_PROGRAM, paramsTiltAlign)
            self.runProgram(ALIGNLOG_PROGRAM, {'-s': "> taSolution.log"},
                            cwd=self._getExtraPath(tsId))

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {TILT_ALIGN_PROGRAM} or {ALIGNLOG_PROGRAM} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def translateFiducialPointModelStep(self, tsId: str):
        if tsId in self.failedItems:
            return
        try:
            # Check that previous steps have been completed satisfactorily
            noGapsFid = self.getExtraOutFile(tsId, suffix="noGaps", ext=FID_EXT)
            if exists(noGapsFid):
                paramsNoGapModel2Point = {
                    "-InputFile": noGapsFid,
                    "-OutputFile": self.getExtraOutFile(tsId, suffix="noGaps_fid", ext=TXT_EXT)
                }
                self.runProgram(MODEL2POINT_PROGRAM, paramsNoGapModel2Point)
        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {MODEL2POINT_PROGRAM} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def createOutputStep(self, lMk: LandmarkModel, ts: TiltSeries):
        tsId = lMk.getTsId()
        if tsId in self.failedItems:
            return
        try:
            with self._lock:
                # Fiducial model
                self.createOutModel(lMk)
                # Tilt-series
                self.createOutTs(ts, self._getInTsSet(pointer=True))

        except sqlite3.OperationalError:
            # Let the decorator retry
            raise

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output '
                                f'with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        fidModelNoGaps = getattr(self, OUTPUT_FIDUCIAL_NO_GAPS_NAME, None)
        tsCoords = getattr(self, OUTPUT_TS_COORDINATES_NAME, None)

        if fidModelNoGaps is not None:
            summary.append("Fiducial models generated with no gaps: "
                           f"{fidModelNoGaps.getSize()}")

        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            summary.append("Transformation matrices updated from the "
                           f"input tilt-series: {output.getSize()}")

        if tsCoords is not None:
            summary.append("Fiducial 3D coordinates calculated: "
                           f"{tsCoords.getSize()}")

        if not summary:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        fidModelNoGaps = getattr(self, OUTPUT_FIDUCIAL_NO_GAPS_NAME, None)
        if fidModelNoGaps is not None:
            methods.append("Solved fiducials alignment for "
                           f"{fidModelNoGaps.getSize()} "
                           "tilt-series using IMOD *tiltalign* command.")

        return methods

    def _warnings(self):
        pass

    # --------------------------- UTILS functions -----------------------------
    def getRotationType(self):
        return {
            0: 0,
            1: -1,
            2: 3,
            3: 1
        }[self.rotationSolutionType.get()]

    def getMagnificationType(self):
        return {
            0: 0,
            1: 3,
            2: 1
        }[self.magnificationSolutionType.get()]

    def getTiltAngleType(self):
        return {
            0: 0,
            1: 5,
            2: 2
        }[self.tiltAngleSolutionType.get()]

    def getSkewType(self):
        # NOTE: the code used in this method is equivalent to the following (in Python 3.10+)
        # valor = self.distortionSolutionType.get()
        # if valor == 0:
        #     return 0
        # elif valor in (1, 2):
        #     return 3
        # else:
        #     raise KeyError(f"Valor inesperado: {valor}")
        return {
            0: 0,
            1: 3,
            2: 3
        }[self.distortionSolutionType.get()]

    def getStretchType(self):
        return {
            0: 0,
            2: 0,
            1: 3
        }[self.distortionSolutionType.get()]

    def getSurfaceToAnalyze(self) -> int:
        return 2 if self.twoSurfaces else 1

    def getInputSetOfLandmarks(self, pointer: bool = False) -> Union[Pointer, SetOfLandmarkModels]:
        return self.inputSetOfLandmarkModels.get() if not pointer else self.inputSetOfLandmarkModels

    def _getInTsSet(self, pointer: bool = False) -> Union[Pointer, SetOfTiltSeries]:
        return self.getInputSetOfLandmarks().getSetOfTiltSeries(pointer=pointer)

    def getTltFilePath(self, tsId) -> str:
        return self.getExtraOutFile(tsId, suffix="interpolated", ext=TLT_EXT)

    @retry_on_sqlite_lock(log=logger)
    def createOutModel(self, lMk: LandmarkModel) -> None:
        tsId = lMk.getTsId()
        outputFn = self.getExtraOutFile(tsId, suffix="noGaps_fid", ext=TXT_EXT)
        fiducialNoGapFilePath = self.getExtraOutFile(tsId, suffix="noGaps_fid", ext=TXT_EXT)
        if exists(outputFn) and exists(fiducialNoGapFilePath):
            with self._lock:
                output = self.getOutputFiducialModel(self._getInTsSet(pointer=True),
                                                     attrName=OUTPUT_FIDUCIAL_NO_GAPS_NAME,
                                                     suffix="NoGaps")
                fiducialModelNoGapPath = self.getExtraOutFile(tsId,
                                                              suffix="noGaps",
                                                              ext=FID_EXT)
                landmarkModelNoGapsFilePath = self.getExtraOutFile(tsId,
                                                                   suffix="noGaps",
                                                                   ext=SFID_EXT)
                landmarkModelNoGapsResidPath = self.getExtraOutFile(tsId,
                                                                    suffix="resid",
                                                                    ext=TXT_EXT)

                landmarkModelNoGaps = LandmarkModel(tsId=tsId,
                                                    fileName=landmarkModelNoGapsFilePath,
                                                    modelName=fiducialModelNoGapPath,
                                                    size=lMk.getSize(),
                                                    hasResidualInfo=True)

                fiducialNoGapList = fiducialModel2List(fiducialNoGapFilePath)
                fiducialNoGapsResidList = fidResidualModel2List(landmarkModelNoGapsResidPath)
                prevTiltIm = 0
                chainId = 0
                indexFake = 0
                firstExec = True

                for fiducial in fiducialNoGapList:
                    if (int(float(fiducial[2])) <= prevTiltIm) or firstExec:
                        chainId += 1
                        firstExec = False
                    prevTiltIm = int(float(fiducial[2]))

                    if indexFake < len(fiducialNoGapsResidList) and fiducial[2] == \
                            fiducialNoGapsResidList[indexFake][2]:
                        landmarkModelNoGaps.addLandmark(xCoor=fiducial[0],
                                                        yCoor=fiducial[1],
                                                        tiltIm=fiducial[2] + 1,
                                                        chainId=chainId,
                                                        xResid=fiducialNoGapsResidList[indexFake][3],
                                                        yResid=fiducialNoGapsResidList[indexFake][4])
                        indexFake += 1

                    else:
                        landmarkModelNoGaps.addLandmark(xCoor=fiducial[0],
                                                        yCoor=fiducial[1],
                                                        tiltIm=fiducial[2] + 1,
                                                        chainId=chainId,
                                                        xResid=float('nan'),
                                                        yResid=float('nan'))

                output.append(landmarkModelNoGaps)
                output.update(landmarkModelNoGaps)
                output.write(output)
                self._store(output)
                self.closeOutputsForStreaming()
        else:
            logger.error(redStr(f'tsId = {tsId} -> Output file {outputFn} was not generated. Skipping... '))

    @staticmethod
    def _convertTxt2Fid(input_txt_file, output_fid_file):
        """
        input_txt_file: x y view chain
        - x, y are the coordinates of the fiducials
        - vView is the number of the image (0-indexed)
        - chain is ethe chain number or object (0-indexed)
        """
        fiducials_by_chain = {}
        with open(input_txt_file, 'r') as f_in:
            for line in f_in:
                parts = line.strip().split()
                try:
                    x = float(parts[0])
                    y = float(parts[1])
                    view = int(parts[2])
                    chain = int(parts[3])

                    if chain not in fiducials_by_chain:
                        fiducials_by_chain[chain] = []
                    fiducials_by_chain[chain].append((view, x, y))
                except ValueError:
                    print(f"Warning: Skipping line, it is not in the proper format: {line.strip()}")

        with open(output_fid_file, 'w') as f_out:
            # IMOD expects a number of chains at the beginning of the file
            for chain_id in sorted(fiducials_by_chain.keys()):
                fiducials = fiducials_by_chain[chain_id]
                for view, x, y in fiducials:
                    # IMOD coordinates ar integers * 100
                    # and the views are  0-indexed.
                    f_out.write(f"{1}\t{chain_id}\t{x}\t{y}\t{float(view - 1)}\n")
