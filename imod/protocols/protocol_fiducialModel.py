# *****************************************************************************
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
from os.path import exists, abspath
import pyworkflow.protocol.params as params
from imod.convert.convert import fiducialModel2List
from imod.protocols.protocol_base_ts_align import ProtImodBaseTsAlign
from imod.protocols.protocol_base_xcorr_fidmodel import ProtImodBaseXcorrFidModel
import pyworkflow.utils.path as path
from imod.constants import (TLT_EXT, XF_EXT, FID_EXT, TXT_EXT, SEED_EXT,
                            SFID_EXT, OUTPUT_FIDUCIAL_GAPS_NAME,
                            FIDUCIAL_MODEL, PATCH_TRACKING, PT_FRACTIONAL_OVERLAP, PT_NUM_PATCHES, AUTOFIDSEED_PROGRAM,
                            BEADTRACK_PROGRAM, IMODCHOPCONTS_PROGRAM, MODEL2POINT_PROGRAM)
from imod.protocols.protocol_xCorrPrealignment import TILT_XCORR_PROGRAM
from pyworkflow.object import Set
from pyworkflow.protocol import STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.utils import Message, cyanStr, redStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import TiltSeries, SetOfLandmarkModels, LandmarkModel

logger = logging.getLogger(__name__)


class ProtImodFiducialModel(ProtImodBaseTsAlign, ProtImodBaseXcorrFidModel, ProtStreamingBase):
    """
    Constructs fiducial tracking models and performs tilt-series alignment
    using IMOD-based workflows for electron tomography applications.
    The protocol is designed to identify, track, and organize fiducial
    markers across tilt-series images so they can later support accurate
    tomographic reconstruction, geometric correction, and alignment
    refinement. It supports both traditional fiducial bead tracking and
    correlation-based patch tracking approaches, allowing the protocol to
    adapt to datasets with visible fiducials as well as datasets where
    alignment must rely on intrinsic image features. More info:
    https://bio3d.colorado.edu/imod/doc/man/autofidseed.html
    https://bio3d.colorado.edu/imod/doc/man/beadtrack.html
    https://bio3d.colorado.edu/imod/doc/man/model2point.html

    AI Generated:

    Generate Fiducial Model (ProtImodFiducialModel) - User Manual
        Overview

        The Generate Fiducial Model protocol creates fiducial tracking
        models for electron tomography tilt-series using IMOD alignment
        procedures. Its primary purpose is to establish reliable
        correspondences between projections acquired at different tilt
        angles so that the geometry of the acquisition can be accurately
        reconstructed during tomographic processing.

        In electron tomography, fiducial markers such as gold beads are
        frequently deposited on the specimen before data collection.
        These markers act as stable reference points that remain visible
        across the entire tilt-series. By tracking their positions
        through the projections, the protocol determines how images must
        be shifted, rotated, and aligned relative to one another. This
        alignment stage is one of the most critical steps in tomographic
        reconstruction because even small tracking errors can propagate
        into severe reconstruction artifacts.

        The protocol supports two complementary alignment strategies.
        The first strategy is fiducial-based alignment, where discrete
        fiducial markers are automatically detected and tracked through
        the tilt-series. The second strategy is patch tracking, where
        local image regions are correlated between views in order to
        estimate motion and alignment without relying on visible beads.
        This flexibility makes the protocol suitable for a broad range
        of biological specimens and acquisition conditions.

        Fiducial-Based Tracking

        The fiducial-based workflow is intended for datasets containing
        visible fiducial markers, most commonly colloidal gold beads.
        The protocol first generates an initial seed model and then
        iteratively tracks the fiducials through the tilt-series to
        build a consistent landmark model.

        From a biological perspective, fiducial quality strongly
        influences the final tomographic reconstruction. Fiducials
        should ideally be well distributed throughout the field of
        view, clearly separated from one another, and visible across a
        large tilt range. Poorly distributed fiducials often produce
        local distortions or unstable alignments.

        The fiducial diameter parameter should approximately match the
        physical size of the beads used during specimen preparation.
        Accurate bead size estimation improves particle localization and
        tracking robustness, particularly in noisy cryo-electron
        tomography datasets.

        The protocol also allows tracking fiducials on two specimen
        surfaces independently. This option is especially relevant for
        thicker samples where fiducials may exist on both sides of the
        support film. Distinguishing between surfaces can improve
        geometric consistency and reduce alignment ambiguity.

        Automatic Seed Generation

        Automatic seed generation simplifies the creation of initial
        fiducial models by identifying candidate fiducials before the
        tracking stage begins. This approach is particularly useful for
        large datasets where manual fiducial picking would be
        impractical.

        The number of fiducials selected should balance coverage and
        stability. Too few fiducials may provide insufficient geometric
        information, while too many weak or poorly defined fiducials may
        introduce unstable trajectories. In practice, users generally
        benefit from selecting a moderate number of high-quality beads
        distributed across the specimen.

        Once the initial seed model is generated, the protocol can
        optionally reuse the tracked model as a refined seed and repeat
        the tracking stage. This iterative refinement often improves
        alignment quality because the updated model contains more
        accurate fiducial positions than the original seed estimate.

        Sobel Refinement and Fiducial Centering

        The protocol includes optional Sobel-based refinement for
        improving fiducial centering accuracy. This strategy is
        particularly beneficial for noisy cryo-electron tomography data,
        where fiducials may appear blurred or exhibit low contrast.

        In practical biological workflows, enhanced fiducial centering
        improves trajectory smoothness and reduces local alignment
        residuals. This can significantly improve downstream
        reconstruction quality, especially in datasets collected under
        low-dose conditions.

        The sigma parameter associated with Sobel filtering controls the
        scale of the refinement relative to fiducial size. Larger values
        emphasize broader structures, whereas smaller values preserve
        fine localization accuracy. Cryo datasets with substantial noise
        often require stronger smoothing to stabilize detection.

        Patch Tracking Workflow

        Patch tracking provides an alternative alignment strategy when
        fiducials are absent, sparse, or unreliable. Instead of tracking
        beads, the protocol tracks local image patches using
        cross-correlation methods.

        This approach is particularly important for biological samples
        where fiducials cannot be introduced efficiently or where the
        specimen itself contains sufficient structural signal for
        alignment. Examples include large cellular tomograms or thick
        cryo-focused ion beam lamellae.

        The size of the tracking patches determines the balance between
        local sensitivity and robustness. Larger patches improve signal
        stability but may fail to capture local deformations or
        heterogeneous motions. Smaller patches are more sensitive to
        local features but can become unstable in noisy regions.

        Users may define patch overlap or specify the total number of
        patches across the image. Overlapping patches generally improve
        continuity and robustness, especially in datasets with uneven
        signal distribution. However, increasing overlap also increases
        computational cost.

        Subpixel Correlation Refinement

        The protocol supports iterative refinement to achieve subpixel
        alignment accuracy. This refinement is especially important in
        high-resolution tomography workflows where small geometric
        inaccuracies can significantly degrade the reconstruction.

        In biological practice, subpixel refinement is most valuable
        when reconstructing macromolecular complexes or cellular
        structures requiring precise spatial interpretation. Accurate
        alignment contributes directly to improved contrast, sharper
        structural features, and more reliable downstream subtomogram
        analysis.

        Filtering and Correlation Parameters

        Several advanced filtering parameters allow users to optimize
        alignment stability for challenging datasets. Frequency filtering
        can suppress high-frequency noise while preserving the broader
        structural information needed for reliable correlation.

        In cryo-electron tomography, filtering becomes particularly
        important because low-dose imaging conditions produce inherently
        noisy projections. Proper filtering often improves tracking
        consistency and reduces false correlations.

        Trimming and border exclusion parameters help avoid unstable
        image regions near the detector edges or acquisition artifacts.
        Restricting the analysis to well-behaved image areas can
        significantly improve alignment robustness.

        Streaming and High-Throughput Processing

        The protocol is designed to operate in streaming mode, allowing
        incoming tilt-series to be processed continuously as they become
        available. This capability is especially valuable in automated
        acquisition pipelines and facility-scale cryo-electron
        tomography environments.

        Streaming execution allows users to monitor alignment quality
        during data collection instead of waiting until the end of the
        acquisition session. Early detection of problematic datasets can
        save substantial microscope time and improve experimental
        efficiency.

        Outputs and Their Interpretation

        The primary output of the protocol is a set of fiducial landmark
        models containing the tracked fiducial trajectories and
        associated gap information. These landmark models are intended
        for downstream tilt-series alignment and tomographic
        reconstruction workflows.

        Gap-containing trajectories indicate fiducials that were not
        successfully tracked in all projections. Small gaps are common
        and often acceptable, especially at high tilt angles where
        fiducials may become obscured or distorted. However, excessive
        gaps may indicate poor fiducial visibility, low signal quality,
        or unstable tracking conditions.

        Biologically meaningful alignment requires trajectories that are
        spatially consistent across the tilt-series and well distributed
        throughout the specimen area. Users are encouraged to visually
        inspect the tracking quality before reconstruction.

        Practical Recommendations

        For most fiducial-based workflows, users should begin with bead
        diameters closely matching the experimental fiducials and a
        moderate number of tracking points. Well-distributed fiducials
        generally produce more reliable alignments than densely clustered
        fiducials concentrated in a small region.

        Sobel refinement is particularly beneficial for cryogenic
        datasets with weak contrast, whereas conventional resin-embedded
        tomography data may often perform well with default settings.

        Patch tracking should be considered when fiducials are absent or
        unreliable, but users should verify that the specimen contains
        sufficient structural texture to support stable correlations.

        When alignment quality appears unstable, adjusting patch sizes,
        overlap parameters, or filtering settings frequently provides
        substantial improvements.

        Final Perspective

        Fiducial modeling and tilt-series alignment are foundational
        steps in electron tomography because they determine the geometric
        consistency of the final reconstruction. Accurate tracking of
        fiducials or image features directly influences reconstruction
        quality, structural interpretability, and downstream biological
        conclusions.

        Careful parameter selection, proper fiducial distribution, and
        thoughtful validation of tracking results are essential for
        obtaining reliable tomographic reconstructions suitable for
        biological analysis.
    """

    _label = 'Generate fiducial model'
    _possibleOutputs = {OUTPUT_FIDUCIAL_GAPS_NAME: SetOfLandmarkModels}
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @classmethod
    def worksInStreaming(cls):
        return True

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        form.addParam('typeOfModel',
                      params.EnumParam,
                      choices=["Make seed and Track", "Patch Tracking"],
                      default=FIDUCIAL_MODEL,
                      important=True,
                      label='Model generation')
        super().addInTsSetFormParam(form)
        self._patchTrackingForm(form, 'typeOfModel == %i' % PATCH_TRACKING)
        self._fiducialSeedForm(form, 'typeOfModel == %i' % FIDUCIAL_MODEL)
        form.addParallelSection(threads=3, mpi=0)

    def _patchTrackingForm(self, form, condition, levelType=params.LEVEL_NORMAL):
        patchtrack = form.addGroup('Patch Tracking',
                                   expertLevel=levelType,
                                   condition=condition)

        patchtrack.addParam('sizeOfPatches', params.NumericListParam,
                            label='Size of the patches (X,Y)',
                            default='680 680',
                            expertLevel=levelType,
                            help="Size of the  patches to track by correlation. "
                                 "In imod documentation "
                                 "(tiltxcorr: SizeOfPatchesXandY)")

        patchtrack.addParam('patchLayout',
                            params.EnumParam,
                            choices=['Fractional overlap of patches',
                                     'Number of patches'],
                            default=PT_FRACTIONAL_OVERLAP,
                            label='Patch layout',
                            display=params.EnumParam.DISPLAY_HLIST)

        patchtrack.addParam('overlapPatches',
                            params.NumericListParam,
                            default='0.33 0.33',
                            condition='patchLayout == %i' % PT_FRACTIONAL_OVERLAP,
                            label='Fractional overlap of the patches (X,Y)',
                            help="Fractional overlap in X and Y to track by correlation. "
                                 "In imod documentation"
                                 "(tiltxcorr: OverlapOfPatchesXandY)")

        patchtrack.addParam('numberOfPatches',
                            params.NumericListParam,
                            condition='patchLayout == %i' % PT_NUM_PATCHES,
                            label='Number of patches (X,Y)',
                            help="Number of patches in X and Y of the patches. "
                                 "In imod documentation"
                                 "(tiltxcorr: NumberOfPatchesXandY)")

        patchtrack.addParam('iterationsSubpixel',
                            params.IntParam,
                            default=1,
                            label='Iterations to increase subpixel accuracy',
                            help="Number of iteration of each correlation to reduce "
                                 "interpolation of the peak position"
                                 "In imod documentation: (tiltxcorr: IterateCorrelations)")

        self.addTrimingParams(patchtrack, pxTrimCondition=True,
                              correlationCondition=True,
                              levelType=params.LEVEL_ADVANCED)
        self.filteringParametersForm(form, condition=condition,
                                     levelType=params.LEVEL_ADVANCED)
        form.getParam("filterRadius2").setDefault(0.125)
        form.getParam("filterSigma2").setDefault(0.03)

    @staticmethod
    def _fiducialSeedForm(form, condition, levelType=params.LEVEL_NORMAL):
        seedModel = form.addGroup('"Make seed and Track', expertLevel=levelType, condition=condition)
        seedModel.addParam('fiducialDiameter',
                           params.FloatParam,
                           label='Fiducial diameter (nm)',
                           default=10.,
                           important=True,
                           help="Fiducials diameter to be tracked for alignment.")

        seedModel.addParam('twoSurfaces',
                           params.BooleanParam,
                           default=False,
                           label='Find beads on two surfaces?',
                           help="Track fiducials differentiating in which side "
                                "of the sample are located.")

        seedModel.addParam('numberFiducial',
                           params.IntParam,
                           label='Number of fiducials',
                           default=10,
                           help="Number of fiducials to be tracked for alignment.")

        seedModel.addParam('doTrackWithModel', params.BooleanParam,
                           default=True,
                           label="Track with fiducial model as seed",
                           help="Turn the tracked model into new seed and "
                                "repeat tracking.")

        groupGlobalVariables = form.addGroup('Filter variables',
                                             expertLevel=params.LEVEL_ADVANCED, condition=condition)

        groupGlobalVariables.addParam('refineSobelFilter',
                                      params.BooleanParam,
                                      default=True,
                                      label='Refine center with Sobel filter?',
                                      expertLevel=params.LEVEL_ADVANCED,
                                      help='Use edge-detecting Sobel filter '
                                           'to refine the bead positions.')

        groupGlobalVariables.addParam('scalableSigmaForSobelFilter',
                                      params.FloatParam,
                                      default=0.12,
                                      condition='refineSobelFilter',
                                      label='Sobel sigma relative to bead size',
                                      expertLevel=params.LEVEL_ADVANCED,
                                      help='Sigma for gaussian kernel filtering '
                                           'of single beads before Sobel '
                                           'filtering, as fraction of bead '
                                           'diameter. The default sigma is 0.5 '
                                           'pixels regardless of bead size. '
                                           'A value of around 0.12 diameters is '
                                           'needed for higher noise (eg. cryo) '
                                           'data.')

    # -------------------------- INSERT steps functions -----------------------
    def stepsGeneratorStep(self) -> None:
        self._initialize()
        closeSetStepDeps = []
        inTsSet = self.getInputTsSet()
        outTsSet = getattr(self, OUTPUT_FIDUCIAL_GAPS_NAME, None)
        self.readingOutput(outTsSet)

        while True:
            with self._lock:
                inTsIds = set(inTsSet.getTSIds())

            if not inTsSet.isStreamOpen() and Counter(self.tsIdReadList) == Counter(inTsIds):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetsStep,
                                         OUTPUT_FIDUCIAL_GAPS_NAME,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break

            nonProcessedTsIds = inTsIds - set(self.tsIdReadList)
            tsToProcessDict = {tsId: ts.clone() for ts in inTsSet.iterItems()
                               if (tsId := ts.getTsId()) in nonProcessedTsIds  # Only not processed tsIds
                               and ts.getSize() > 0}  # Avoid processing empty TS
            for tsId, ts in tsToProcessDict.items():
                pId = self._insertFunctionStep(self.convertInStep,
                                               ts,
                                               prerequisites=[],
                                               needsGPU=False)

                if self.typeOfModel.get() == FIDUCIAL_MODEL:
                    pId = self._insertFunctionStep(self.generateFiducialSeedStep,
                                                   ts,
                                                   prerequisites=pId,
                                                   needsGPU=False)
                    pId = self._insertFunctionStep(self.generateFiducialModelStep,
                                                   ts,
                                                   prerequisites=pId,
                                                   needsGPU=False)
                else:
                    pId = self._insertFunctionStep(self.xcorrStep,
                                                   tsId,
                                                   prerequisites=pId,
                                                   needsGPU=False)
                    pId = self._insertFunctionStep(self.chopcontsStep,
                                                   tsId,
                                                   prerequisites=pId,
                                                   needsGPU=False)

                pId = self._insertFunctionStep(self.translateFiducialPointModelStep,
                                               tsId,
                                               prerequisites=pId,
                                               needsGPU=False)
                pId = self._insertFunctionStep(self.computeOutputModelsStep,
                                               ts,
                                               prerequisites=pId,
                                               needsGPU=False)
                closeSetStepDeps.append(pId)
                logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                self.tsIdReadList.append(tsId)

            self.refreshStreaming(inTsSet)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        tsSet = self.getInputTsSet()
        self.sRate = tsSet.getSamplingRate()
        self.acq = tsSet.getAcquisition()

    def convertInStep(self, ts: TiltSeries):
        self.convertInputStep(ts, ts.getTsPresentAcqOrders())

    def generateFiducialSeedStep(self, ts: TiltSeries):
        tsId = ts.getTsId()
        if tsId in self.failedItems:
            return
        try:
            logger.info(cyanStr(f'tsId = {tsId}: generating the fiducial seeds...'))
            self.generateTrackCom(ts)
            trackFile = f'{tsId}_track.com'
            if exists(self.getExtraOutFile(tsId, suffix="track", ext="com")):
                paramsAutofidseed = {
                    "-TrackCommandFile": trackFile,
                    "-MinSpacing": 0.85,
                    "-AdjustSizes": "",
                    "-PeakStorageFraction": 1.0,
                    "-TargetNumberOfBeads": self.numberFiducial.get(),
                }
                if self.twoSurfaces:
                    paramsAutofidseed["-TwoSurfaces"] = ""

                self.runProgram(AUTOFIDSEED_PROGRAM, paramsAutofidseed, cwd=self._getExtraPath(tsId))

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {AUTOFIDSEED_PROGRAM} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def generateFiducialModelStep(self, ts: TiltSeries):
        tsId = ts.getTsId()
        if tsId in self.failedItems:
            return
        try:
            logger.info(cyanStr(f'tsId = {tsId}: generating the fiducial model...'))
            paramsBeadtrack = self.genBeadTrackParams(ts)
            self.runProgram(BEADTRACK_PROGRAM, paramsBeadtrack)

            if self.doTrackWithModel:
                # repeat tracking with the current model as seed
                path.copyFile(paramsBeadtrack['-InputSeedModel'],
                              self.getExtraOutFile(tsId, suffix="orig", ext=SEED_EXT))
                path.moveFile(paramsBeadtrack['-OutputModel'],
                              paramsBeadtrack['-InputSeedModel'])

                self.runProgram(BEADTRACK_PROGRAM, paramsBeadtrack)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {BEADTRACK_PROGRAM} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def xcorrStep(self, tsId):
        if tsId in self.failedItems:
            return
        try:
            logger.info(cyanStr(f'tsId = {tsId}: executing the {TILT_XCORR_PROGRAM}...'))
            angleFilePath = self.getExtraOutFile(tsId, ext=TLT_EXT)
            xfFile = self.getExtraOutFile(tsId, ext=XF_EXT)
            borders = self.pxTrim.getListFromValues(caster=str)
            sizePatches = self.sizeOfPatches.getListFromValues(caster=str)

            paramsTiltXCorr = {
                "-InputFile": self.getTmpOutFile(tsId),
                "-OutputFile": self.getExtraOutFile(tsId, suffix="pt", ext=FID_EXT),
                "-RotationAngle": self.acq.getTiltAxisAngle(),
                "-TiltFile": angleFilePath,
                "-FilterRadius2": self.filterRadius2.get(),
                "-FilterSigma1": self.filterSigma1.get(),
                "-FilterSigma2": self.filterSigma2.get(),
                "-BordersInXandY": ",".join(borders),
                "-IterateCorrelations": self.iterationsSubpixel.get(),
                "-SizeOfPatchesXandY": ",".join(sizePatches),
                "-PrealignmentTransformFile": xfFile,
                "-ImagesAreBinned": 1,
            }

            if self.patchLayout.get() == PT_FRACTIONAL_OVERLAP:
                patchesXY = self.overlapPatches.getListFromValues(caster=str)
                paramsTiltXCorr["-OverlapOfPatchesXandY"] = ",".join(patchesXY)
            else:
                numberPatchesXY = self.numberOfPatches.getListFromValues(caster=str)
                paramsTiltXCorr["-NumberOfPatchesXandY"] = ",".join(numberPatchesXY)

            self.runProgram(TILT_XCORR_PROGRAM, paramsTiltXCorr)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {TILT_XCORR_PROGRAM} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def chopcontsStep(self, tsId: str):
        if tsId in self.failedItems:
            return
        try:
            logger.info(cyanStr(f'tsId = {tsId}: executing {IMODCHOPCONTS_PROGRAM}...'))
            paramschopconts = {
                "-InputModel": self.getExtraOutFile(tsId, suffix="pt", ext=FID_EXT),
                "-OutputModel": self.getExtraOutFile(tsId, suffix="gaps", ext=FID_EXT),
                "-MinimumOverlap": 4,
                "-AssignSurfaces": 1,
                "-LengthOfPieces": -1
            }
            self.runProgram(IMODCHOPCONTS_PROGRAM, paramschopconts)
        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {IMODCHOPCONTS_PROGRAM} execution '
                                 f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def translateFiducialPointModelStep(self, tsId: str):
        if tsId in self.failedItems:
            return
        try:
            logger.info(cyanStr(f'tsId = {tsId}: executing {MODEL2POINT_PROGRAM}...'))
            gapsFidFile = self.getExtraOutFile(tsId, suffix='gaps', ext=FID_EXT)
            paramsGapModel2Point = {
                "-InputFile": gapsFidFile,
                "-OutputFile": self.getExtraOutFile(tsId, suffix="gaps_fid", ext=TXT_EXT),
                "-ObjectAndContour": ''
            }
            self.runProgram(MODEL2POINT_PROGRAM, paramsGapModel2Point)
        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {MODEL2POINT_PROGRAM} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def computeOutputModelsStep(self, ts: TiltSeries):
        """ Create the output set of landmark models with gaps. """
        tsId = ts.getTsId()
        if tsId in self.failedItems:
            self.addToOutFailedSet(ts)
            return

        try:
            outputFn = self.getExtraOutFile(tsId, suffix='gaps', ext=FID_EXT)
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
            output = self.getOutputFiducialModel(self.getInputTsSet(pointer=True),
                                                 attrName=OUTPUT_FIDUCIAL_GAPS_NAME,
                                                 suffix="Gaps")
            landmarkModelGapsFilePath = self.getExtraOutFile(tsId,
                                                             suffix='gaps',
                                                             ext=SFID_EXT)
            fiducialModelGapTxtPath = self.getExtraOutFile(tsId,
                                                           suffix="gaps_fid",
                                                           ext=TXT_EXT)

            fiducialGapList = fiducialModel2List(fiducialModelGapTxtPath)
            fiducialDiameter = self.fiducialDiameter.get() * 10  # From nm to angstroms

            landmarkModelGaps = LandmarkModel(tsId=tsId,
                                              tiltSeriesPointer=ts,
                                              fileName=landmarkModelGapsFilePath,
                                              modelName=outputFn,
                                              size=fiducialDiameter,
                                              hasResidualInfo=False)
            landmarkModelGaps.setTiltSeries(ts)

            for index, fiducial in enumerate(fiducialGapList):
                landmarkModelGaps.addLandmark(xCoor=fiducial[2],
                                              yCoor=fiducial[3],
                                              tiltIm=fiducial[4] + 1,
                                              chainId=fiducial[1],
                                              xResid=0,
                                              yResid=0)

            output.append(landmarkModelGaps)
            output.update(landmarkModelGaps)
            output.write(output)
            self._store(output)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        fidModelGaps = getattr(self, OUTPUT_FIDUCIAL_GAPS_NAME, None)
        if fidModelGaps is not None:
            summary.append(f"Input tilt-series: {self.getInputTsSet().getSize()}\n"
                           "Fiducial models (with gaps) generated: "
                           f"{fidModelGaps.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        fidModelGaps = getattr(self, OUTPUT_FIDUCIAL_GAPS_NAME, None)
        if fidModelGaps is not None:
            methods.append(f"The fiducial model (with gaps) has been computed for "
                           f"{fidModelGaps.getSize()} tilt-series using "
                           f"the IMOD *{BEADTRACK_PROGRAM}* command.")
        return methods

    # --------------------------- UTILS functions -----------------------------
    def generateTrackCom(self, ts: TiltSeries):
        tsId = ts.getTsId()
        logger.info(cyanStr(f'tsId = {tsId}: generating the tracking command file...'))
        fiducialDiameterPixel, boxSizeXandY, scaling = self.getFiducialParams()
        # Absolute paths are because the cwd will be self._getExtraPath(tsId) to avoid the generic dir
        # autofidseed.dir (one for each tilt-series processed, but in the same place and with the same name
        # unless the cwd is changed) to be generated, deleted or used when processing in parallel
        paramsDict = {
            'imageFile': abspath(self.getTmpOutFile(tsId)),
            'inputSeedModel': abspath(self.getExtraOutFile(tsId, ext=SEED_EXT)),
            'outputModel': abspath(self.getExtraOutFile(tsId, suffix="gaps", ext=FID_EXT)),
            'tiltFile': abspath(self.getExtraOutFile(tsId, ext=TLT_EXT)),
            'rotationAngle': self.acq.getTiltAxisAngle(),
            'fiducialDiameter': fiducialDiameterPixel,
            'samplingRate': self.sRate / 10,
            'scalableSigmaForSobelFilter': self.scalableSigmaForSobelFilter.get(),
            'boxSizeXandY': boxSizeXandY,
            'distanceRescueCriterion': 10 * scaling,
            'postFitRescueResidual': 2.5 * scaling,
            'maxRescueDistance': 2.5 * scaling,
            'minDiamForParamScaling': 12.5,
            'deletionCriterionMinAndSD': f"{0.04 * scaling:0.3f},2.0",
        }

        hasAlign = ts.hasAlignment()
        self.createBeadTrackTemplate(ts, paramsDict, hasAlign)

    def createBeadTrackTemplate(self,
                                ts: TiltSeries,
                                paramsDict: dict,
                                hasAlignment: bool = False):

        tsId = ts.getTsId()
        trackFilePath = self.getExtraOutFile(tsId, suffix="track", ext="com")
        template = """# Command file for running BEADTRACK
#
####CreatedVersion####4.9.12
#
# For beads lighter than background, add a line with "LightBeads"
#
# To restrict tilt alignment to a subset of views, add a line with:
# "MaxViewsInAlign #_of_views"
#
# To exclude views, add a line "SkipViews view_list"; with the list of views
#
# To specify sets of views to be grouped separately in automapping, add a line
# "SeparateGroup view_list" with the list of views, one line per group
#
$beadtrack -StandardInput
LowPassCutoffInverseNm  0.3
ImageFile	%(imageFile)s
ImagesAreBinned	1
InputSeedModel	%(inputSeedModel)s
OutputModel	%(outputModel)s
RotationAngle	%(rotationAngle).2f 
TiltFile	%(tiltFile)s
TiltDefaultGrouping	7
MagDefaultGrouping	5
RotDefaultGrouping	1
PixelSize   %(samplingRate)f
BeadDiameter	%(fiducialDiameter).2f
FillGaps
MaxGapSize	5
RoundsOfTracking	2
#
# Set this to 1 to track in local areas
LocalAreaTracking	1
LocalAreaTargetSize	1000
MinBeadsInArea	8
MinOverlapBeads	5
#
# CONTROL PARAMETERS FOR EXPERTS, EXPERIMENTATION, OR SPECIAL CASES
#
# minimum range of tilt angles for finding axis and for finding tilts
MinViewsForTiltalign	4
MinTiltRangeToFindAxis	10.0
MinTiltRangeToFindAngles	20.0
BoxSizeXandY	%(boxSizeXandY)d,%(boxSizeXandY)d
MaxBeadsToAverage	4
# points and minimum for extrapolation
PointsToFitMaxAndMin	7,3
# fraction of mean, and # of SD below mean: density criterion for rescue
DensityRescueFractionAndSD	0.6,1.0
# distance criterion for rescue
DistanceRescueCriterion	%(distanceRescueCriterion).2f
# relaxation of criterion for density and distance rescues
RescueRelaxationDensityAndDistance	0.7,0.9
# distance for rescue after fit
PostFitRescueResidual	%(postFitRescueResidual).2f
# relaxation of density criterion, maximum radius to search
DensityRelaxationPostFit	0.9
MaxRescueDistance	%(maxRescueDistance).2f
# Max and min residual changes to use to get mean and SD change
ResidualsToAnalyzeMaxAndMin	9,5
# minimum residual difference, criterion # of sd's
DeletionCriterionMinAndSD	%(deletionCriterionMinAndSD)s
MinDiamForParamScaling %(minDiamForParamScaling).1f
"""

        if self.refineSobelFilter:
            template += "\nSobelFilterCentering"
            template += "\nScalableSigmaForSobel   %(scalableSigmaForSobelFilter)f"

        if hasAlignment:
            XfFileName = abspath(self.getExtraOutFile(tsId, ext=XF_EXT))
            template += f"\nPrealignTransformFile {XfFileName}"

        with open(trackFilePath, 'w') as f:
            f.write(template % paramsDict)

    def getFiducialParams(self):
        """ Precalculate a few params. """
        fiducialDiameterPixel = self.fiducialDiameter.get() / (self.sRate / 10)
        # formulas from bin/copytomos
        boxSizeXandY = max(3.3 * fiducialDiameterPixel + 2, 2 * fiducialDiameterPixel + 20, 32)
        boxSizeXandY = min(512, 2 * int(boxSizeXandY / 2))
        scaling = fiducialDiameterPixel / 12.5 if fiducialDiameterPixel > 12.5 else 1

        return fiducialDiameterPixel, boxSizeXandY, scaling

    def getOutputFiducialModelGaps(self, inputSet):
        if self.FiducialModelGaps:
            self.FiducialModelGaps.enableAppend()
        else:
            fidModelGaps = self._createSetOfLandmarkModels(suffix='Gaps')

            fidModelGaps.copyInfo(inputSet)
            fidModelGaps.setSetOfTiltSeries(inputSet)
            fidModelGaps.setHasResidualInfo(False)
            fidModelGaps.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_FIDUCIAL_GAPS_NAME: fidModelGaps})
            self._defineSourceRelation(inputSet, fidModelGaps)

        return self.FiducialModelGaps

    def genBeadTrackParams(self, ts: TiltSeries) -> dict:
        tsId = ts.getTsId()
        fiducialDiameterPixel, boxSizeXandY, scaling = self.getFiducialParams()

        paramsBeadtrack = {
            "-ImageFile": self.getTmpOutFile(tsId),
            "-InputSeedModel": self.getExtraOutFile(tsId, ext=SEED_EXT),
            "-OutputModel": self.getExtraOutFile(tsId, suffix="gaps", ext=FID_EXT),
            "-TiltFile": self.getExtraOutFile(tsId, ext=TLT_EXT),
            "-RotationAngle": ts.getAcquisition().getTiltAxisAngle(),
            "-PixelSize": self.sRate / 10,
            "-ImagesAreBinned": 1,
            "-TiltDefaultGrouping": 7,
            "-MagDefaultGrouping": 5,
            "-RotDefaultGrouping": 1,
            "-MinViewsForTiltalign": 4,
            "-BeadDiameter": fiducialDiameterPixel,
            "-FillGaps": 1,
            "-MaxGapSize": 5,
            "-MinTiltRangeToFindAxis": 10.0,
            "-MinTiltRangeToFindAngles": 20.0,
            "-BoxSizeXandY": f"{boxSizeXandY},{boxSizeXandY}",
            "-RoundsOfTracking": 2,
            "-LocalAreaTracking": 1,
            "-LocalAreaTargetSize": 1000,
            "-MinBeadsInArea": 8,
            "-MinOverlapBeads": 5,
            "-MaxBeadsToAverage": 4,
            "-PointsToFitMaxAndMin": '7,3',
            "-DensityRescueFractionAndSD": '0.6,1.0',
            "-DistanceRescueCriterion": 10 * scaling,
            "-RescueRelaxationDensityAndDistance": '0.7,0.9',
            "-PostFitRescueResidual": 2.5 * scaling,
            "-DensityRelaxationPostFit": 0.9,
            "-MaxRescueDistance": 2.5 * scaling,
            "-ResidualsToAnalyzeMaxAndMin": '9,5',
            "-DeletionCriterionMinAndSD": f"{0.04 * scaling:0.3f},2.0",
            "-MinDiamForParamScaling": 12.5,
            "-LowPassCutoffInverseNm": 0.3
        }

        if self.refineSobelFilter:
            paramsBeadtrack["-SobelFilterCentering"] = ""
            paramsBeadtrack["-ScalableSigmaForSobel"] = self.scalableSigmaForSobelFilter.get()

        # # Excluded views
        # if ts.hasExcludedViews():
        #     excludedViews = ts.getTsExcludedViewsIndices(ts.getTsPresentAcqOrders())
        #     paramsBeadtrack["-SkipViews"] = ','.join(map(str, excludedViews))

        if ts.hasAlignment():
            paramsBeadtrack["-prexf"] = self.getExtraOutFile(tsId, ext=XF_EXT)

        return paramsBeadtrack
