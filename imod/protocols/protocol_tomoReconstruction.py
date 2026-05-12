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
import pyworkflow.protocol.params as params
from pyworkflow.protocol import ProtStreamingBase
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.utils import Message, cyanStr, redStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import Tomogram, SetOfTomograms, TiltSeries
from imod import Plugin
from imod.protocols import ProtImodBase
from imod.constants import (TLT_EXT, ODD, EVEN, MRC_EXT,
                            OUTPUT_TOMOGRAMS_NAME, TRIMVOL_PROGRAM, TILT_PROGRAM)

logger = logging.getLogger(__name__)


class ProtImodTomoReconstruction(ProtImodBase, ProtStreamingBase):
    """
    Reconstructs three-dimensional tomograms from aligned tilt-series
    using the IMOD tomography workflow. The protocol transforms a set
    of two-dimensional projection images acquired at multiple tilt
    angles into volumetric reconstructions suitable for structural
    interpretation, segmentation, visualization, and subtomogram
    analysis in cryo-electron tomography experiments. More info:
        https://bio3d.colorado.edu/imod/doc/man/tilt.html

    AI Generated:

    Tomogram Reconstruction (ProtImodTomoReconstruction) - User Manual
        Overview

        The Tomogram Reconstruction protocol generates three-dimensional
        tomograms from aligned tilt-series using the IMOD reconstruction
        framework. Its primary objective is to recover the volumetric
        organization of biological specimens from a collection of
        two-dimensional projections acquired across a range of tilt
        angles. This process is one of the central stages in cryo-electron
        tomography because it transforms aligned experimental images into
        interpretable three-dimensional cellular or macromolecular volumes.

        In practical biological workflows, tomographic reconstruction is
        typically performed after motion correction, preprocessing,
        alignment, and tilt-series refinement. The resulting tomograms
        provide the structural context needed for visualization of organelles,
        membranes, viral particles, protein complexes, and cellular
        architecture in situ. The protocol is designed for both routine
        reconstructions and more advanced tomography pipelines where
        reconstruction quality strongly affects downstream interpretation.

        More information:
        https://bio3d.colorado.edu/imod/doc/man/tilt.html

        Inputs and Reconstruction Strategy

        The protocol requires aligned tilt-series as input. These tilt-series
        should already contain reliable geometric alignment information,
        including tilt angles and image transformations. Accurate alignment
        is essential because reconstruction quality directly depends on the
        consistency of the projection geometry across the complete angular
        range.

        The reconstruction process assumes that the specimen was imaged by
        rotation around a tilt axis with limited deviations from an ideal
        acquisition geometry. Under these conditions, the protocol combines
        the information from all projections into a volumetric representation
        of the specimen. Proper preprocessing and alignment before this stage
        are particularly important because reconstruction artifacts are often
        amplified when inconsistencies remain unresolved.

        Biological users commonly reconstruct tomograms for visualization of
        intracellular organization, extraction of subtomograms, membrane
        tracing, particle picking, or quantitative spatial analysis. The
        protocol therefore acts as the bridge between image alignment and
        structural interpretation.

        Tomogram Dimensions and Volume Geometry

        One of the most important parameters is the tomogram thickness,
        which defines the reconstructed size along the beam direction.
        Choosing an appropriate thickness is biologically important because
        values that are too small may truncate relevant structures, whereas
        unnecessarily large volumes increase computation time and storage
        requirements while adding mostly empty regions.

        The protocol also allows reconstruction of a reduced tomogram width.
        This option is especially useful when the biological region of
        interest occupies only a central portion of the field of view.
        Cropping the reconstructed width can substantially reduce storage
        needs and accelerate downstream analysis while preserving the most
        relevant structural information.

        In practice, users should ensure that the selected dimensions fully
        encompass the specimen across all tilt angles. Failure to do so may
        result in incomplete reconstructions or clipping of biologically
        meaningful regions.

        Reconstruction Offsets and Spatial Adjustments

        Advanced users may apply shifts or offsets during reconstruction.
        These adjustments modify the spatial positioning of reconstructed
        slices and can help optimize the placement of the specimen inside
        the final tomogram volume.

        Shifts along the reconstruction axes are particularly useful when
        the region of interest is not centered in the acquisition images.
        Correcting these offsets during reconstruction can improve the
        usability of the final volume and simplify downstream subtomogram
        extraction or segmentation workflows.

        The protocol also supports tilt angle and tilt axis offsets. These
        parameters provide additional flexibility when acquisition geometry
        requires refinement. Small corrections may improve reconstruction
        quality when the imported acquisition metadata does not perfectly
        represent the experimental setup.

        From a biological perspective, careful adjustment of these values
        becomes especially important for high-resolution tomography studies
        or for datasets containing elongated or asymmetric structures.

        Radial Filtering and Noise Control

        Tomographic reconstruction inevitably amplifies noise because the
        available angular information is incomplete. To mitigate this issue,
        the protocol provides radial filtering controls that regulate the
        balance between signal preservation and noise suppression.

        The radial filtering parameters define how high-frequency information
        is attenuated during reconstruction. Conservative filtering generally
        produces smoother tomograms with reduced noise, whereas weaker
        filtering preserves finer structural details but may increase grainy
        artifacts and reconstruction streaks.

        For many biological applications, moderate filtering provides the
        best compromise between interpretability and structural preservation.
        However, datasets intended for subtomogram averaging or quantitative
        analysis may benefit from more careful optimization of the filtering
        behavior.

        Users should remember that excessive filtering may erase biologically
        meaningful high-resolution information, while insufficient filtering
        may complicate segmentation or interpretation because of elevated
        noise levels.

        Super-Sampling and Artifact Reduction

        The protocol supports super-sampling strategies to reduce
        reconstruction artifacts generated by discrete backprojection.
        These artifacts often appear as directional rays or streaks that
        follow the acquisition geometry.

        Super-sampling improves reconstruction smoothness by performing the
        reconstruction at a finer intermediate sampling before reducing the
        data back to the desired resolution. In practice, this option can
        produce cleaner tomograms with improved visual continuity and reduced
        directional artifacts.

        The improvement is especially relevant in datasets intended for
        subtomogram averaging or detailed structural interpretation. Although
        higher super-sampling values increase computational cost, they may
        enhance the visual quality of the reconstructed volume and reduce
        reconstruction-induced artifacts.

        SIRT-Like Filtering

        The protocol also provides an option that approximates the visual
        characteristics of iterative SIRT reconstructions through filtering
        approaches. This strategy can produce tomograms with smoother density
        transitions and reduced high-frequency noise compared to conventional
        weighted backprojection.

        Biologically, SIRT-like reconstructions are often preferred for
        visualization, segmentation, and manual annotation because cellular
        structures can appear more continuous and easier to interpret.
        However, smoother reconstructions may also reduce apparent sharpness,
        so users should evaluate the balance between interpretability and
        preservation of fine structural details.

        GPU and Computational Considerations

        The reconstruction workflow supports both CPU and GPU execution.
        GPU acceleration is particularly valuable for large tomographic
        datasets because reconstruction is computationally intensive and
        often represents one of the slowest stages in cryo-electron
        tomography processing pipelines.

        Depending on the reconstruction strategy and hardware availability,
        GPU execution can dramatically reduce processing time while preserving
        reconstruction quality. This capability is especially important in
        facility-scale workflows or streaming environments where many
        tilt-series must be processed continuously.

        Streaming Reconstruction Workflows

        The protocol is compatible with streaming operation, allowing
        tomograms to be reconstructed progressively as tilt-series become
        available. This capability is highly useful during automated data
        acquisition sessions because users can evaluate reconstruction
        quality in near real time.

        Streaming workflows provide early feedback regarding alignment
        quality, specimen thickness, contamination, and acquisition
        performance. This enables rapid troubleshooting during microscope
        sessions and improves the efficiency of large-scale tomography
        projects.

        Odd-Even Reconstructions

        The protocol can optionally reconstruct odd and even subsets of
        projections independently. These paired reconstructions are useful
        for quality assessment, validation strategies, and downstream
        resolution estimation approaches commonly employed in tomography
        workflows.

        Separating odd and even projections may also support advanced
        denoising or validation procedures where independent reconstructions
        are required to estimate reproducibility or structural consistency.

        Outputs and Biological Interpretation

        After execution, the protocol produces reconstructed tomograms
        associated with the input tilt-series. These tomograms preserve
        the acquisition metadata and spatial information necessary for
        downstream analysis and visualization.

        The reconstructed volumes can subsequently be used for subtomogram
        averaging, segmentation, membrane tracing, particle picking,
        structural annotation, or correlation with fluorescence microscopy
        data. Their interpretability depends strongly on alignment quality,
        reconstruction parameters, angular coverage, and specimen thickness.

        Biological users should interpret reconstructed densities carefully,
        particularly in regions affected by the missing wedge, low contrast,
        or strong noise amplification. Structural continuity and apparent
        resolution may vary substantially across different regions of the
        tomogram.

        Practical Recommendations

        In routine biological practice, it is generally advisable to begin
        with moderate reconstruction thickness and default filtering values,
        then evaluate the resulting tomograms visually before performing
        aggressive optimization. Overly sharp filtering parameters may
        emphasize noise, whereas excessive smoothing may obscure relevant
        structural details.

        Super-sampling by a factor of two is often a reasonable compromise
        between artifact reduction and computational efficiency. Larger
        values are typically reserved for demanding applications involving
        subtomogram averaging or detailed structural interpretation.

        When acquisition metadata are uncertain, users should carefully
        inspect reconstructed tomograms for signs of tilt-axis errors,
        elongation artifacts, or geometric distortions. Small corrections
        to offsets or reconstruction geometry can significantly improve
        the biological quality of the final volume.

        Final Perspective

        Tomographic reconstruction is one of the defining steps in
        cryo-electron tomography because it transforms aligned projection
        images into interpretable three-dimensional biological volumes.
        The quality of the final tomogram depends not only on reconstruction
        settings but also on the accuracy of the preceding alignment and
        preprocessing stages. Careful optimization of reconstruction
        geometry, filtering behavior, and artifact reduction strategies is
        therefore essential for obtaining biologically meaningful tomograms
        suitable for structural and cellular interpretation.
    """

    _label = 'Tomo reconstruction'
    _possibleOutputs = {OUTPUT_TOMOGRAMS_NAME: SetOfTomograms}
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsReadList = []

    @classmethod
    def worksInStreaming(cls):
        return True

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        super().addInTsSetFormParam(form)
        form.addParam('tomoThickness',
                      params.IntParam,
                      default=300,
                      label='Tomogram thickness (voxels)',
                      important=True,
                      help='Size in voxels of the tomogram along the z '
                           'axis (beam direction).')
        form.addParam('tomoWidth',
                      params.IntParam,
                      default=0,
                      label='Tomogram width (voxels)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Number of pixels to cut out in X, centered on the middle in X. '
                           'Leave 0 for default X.')
        lineShift = form.addLine('Tomogram shift (px)',
                                 expertLevel=params.LEVEL_ADVANCED,
                                 help="This entry allows one to shift the reconstructed"
                                      " slice in X or Z before it is output.  If the "
                                      " X shift is positive, the slice will be shifted to "
                                      " the right, and the output will contain the left "
                                      " part of the whole potentially reconstructable area. "
                                      " If the Z shift is positive, the slice is shifted "
                                      " upward. The Z entry is optional and defaults to 0 when "
                                      " omitted.")
        lineShift.addParam('tomoShiftX',
                           params.FloatParam,
                           default=0,
                           label=' in X ')
        lineShift.addParam('tomoShiftZ',
                           params.FloatParam,
                           default=0,
                           label=' in Z ')
        lineoffSet = form.addLine('Offset (deg) of the ',
                                  expertLevel=params.LEVEL_ADVANCED,
                                  help="* Tilt angle offset: pply an angle offset in "
                                       "degrees to all tilt angles. This offset "
                                       "positively rotates the reconstructed sections "
                                       "anticlockwise.\n* Tilt axis offset: Apply an "
                                       "offset to the tilt axis in a stack of full-sized "
                                       "projection images, cutting the X-axis at NX/2 + "
                                       "offset instead of NX/2.")
        lineoffSet.addParam('angleOffset',
                            params.FloatParam,
                            default=0,
                            label='Tilt angles ',
                            help='Apply an angle offset in degrees to all tilt '
                                 'angles. This offset positively rotates the '
                                 'reconstructed sections anticlockwise.')
        lineoffSet.addParam('tiltAxisOffset',
                            params.FloatParam,
                            default=0,
                            label='Tilt axis',
                            help='Apply an offset to the tilt axis in a stack of '
                                 'full-sized projection images, cutting the '
                                 'X-axis at NX/2 + offset instead of NX/2. The '
                                 'DELXX entry is optional and defaults to 0 '
                                 'when omitted.')
        form.addParam('superSampleFactor',
                      params.IntParam,
                      default=2,
                      label='Super-sampling factor',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Compute slices in pixels smaller by this factor to reduce artifacts.'
                           'Super-sampling refers to computing the back projection in a slice '
                           'larger by an integer factor in each dimension, which is done here in '
                           'one of two ways: by interpolating the projection data at smaller '
                           'intervals during backprojection, or by expanding the input lines by '
                           'the factor using sync interpolation (padding in Fourier space).  '
                           'Super-sampling will reduce the rays along the projection angles that '
                           'appear to reflect from the edges of a Fourier transform of an X/Z slice.'
                           'These rays result from back-projecting into discrete pixels and represent'
                           'inappropriate information generated by the transitions between successive '
                           'pixels along a backprojection ray.  Super-sampling by 2 will remove most '
                           'of these rays, especially oblique ones.  The additional benefit (amount '
                           'of change in the image) of going from 2 to 3 is about 10% as large as '
                           'that of super-sampling by 2; it is about 3% going from 3 to 4 and 1.5% '
                           'going from 4 to 8 (the highest allowed value). The super-sampled slice '
                           'is reduced to the original resolution by cropping its Fourier transform. '
                           'Super-sampling alone with the first method does not reduce the rays in '
                           'the corners of the Fourier transform past 0.5 cycle/pixel. These rays '
                           'are  also inappropriate since they originate from  lower-frequency '
                           'information in the projection images, so they are removed from the '
                           'cropped Fourier transform before inverting it. This removal has an '
                           'added benefit about 1/3 as large as the benefit from supersampling by 2. '
                           'The effect of these removals is a subtle change in the noise, and a '
                           'benefit may show up only with subvolume averaging.\n'
                           'The additional effect of expanding the input lines is to avoid attenu-'
                           'ating frequencies past half Nyquist (0.25/pixel) any more than is '
                           'achieved with the falloff of the radial filter.  This will be noticeable '
                           'in a tomogram and would be particularly helpful if setting the radial '
                           'filter cutoff higher than the default for subvolume averaging.\n')
        form.addParam('fakeInteractionsSIRT',
                      params.IntParam,
                      default=0,
                      label='Iterations of a SIRT-like equivalent filter',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Modify the radial filter to produce a '
                           'reconstruction equivalent to the one produced by '
                           'the given number of iterations of SIRT. The '
                           'Gaussian filter is applied at the high-frequency '
                           'end of the filter. The functioning of this filter '
                           'is described in: \n\t'
                           'https://bio3d.colorado.edu/imod/doc/man/tilt.html'
                           'This entry corresponds to the imod parameter – FakeSIRTiterations')
        groupRadialFrequencies = form.addGroup('Radial filtering',
                                               help='This entry controls low-pass filtering with the radial weighting '
                                                    'function. The radial weighting function is linear away from the '
                                                    'origin out to the distance in reciprocal space specified by the '
                                                    'first value, followed by a Gaussian fall-off determined by the '
                                                    'second value.',
                                               expertLevel=params.LEVEL_ADVANCED)
        groupRadialFrequencies.addParam('radialFirstParameter',
                                        params.FloatParam,
                                        default=0.35,
                                        label='Cutoff linear region',
                                        help='Linear region of the filter. This is the first'
                                             'parameter of the -RADIAL flag in IMOD tilt command.'
                                             'If the cutoff is greater than 1 the distances are interpreted '
                                             'as pixels in Fourier space; otherwise they are treated as '
                                             'frequencies in cycles/pixel, which range from 0 to 0.5. '
                                             'Use a cutoff of 0.5 for no low-pass filtering.')
        groupRadialFrequencies.addParam('radialSecondParameter',
                                        params.FloatParam,
                                        default=0.035,
                                        label='Radial fall-off',
                                        help='Gaussian fall-off parameter. The sigma (or standard deviation) '
                                             'of the Gaussian is the second value times 0.707.')
        form.addHidden(params.USE_GPU,
                       params.BooleanParam,
                       default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation. "
                            "Select the one you want to use.")
        form.addHidden(params.GPU_LIST,
                       params.StringParam,
                       default='0',
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="GPU ID. To pick the best available one set 0. "
                            "For a specific GPU set its number ID "
                            "(starting from 1).")
        self.addOddEvenParams(form, isTomogram=True)
        form.addParallelSection(threads=3, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def stepsGeneratorStep(self) -> None:
        """
        This step should be implemented by any streaming protocol.
        It should check its input and when ready conditions are met
        call the self._insertFunctionStep method.
        """
        self._initialize()
        tomoWidth = self.tomoWidth.get()
        inTsSet = self.getInputTsSet()
        outTomoSet = getattr(self, OUTPUT_TOMOGRAMS_NAME, None)
        self.readingOutput(outTomoSet)
        widthWarnTsIds = []
        closeSetStepDeps = []

        while True:
            with self._lock:
                inTsIds = set(inTsSet.getTSIds())

            if not inTsSet.isStreamOpen() and Counter(self.tsReadList) == Counter(inTsIds):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetsStep,
                                         OUTPUT_TOMOGRAMS_NAME,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break

            nonProcessedTsIds = inTsIds - set(self.tsReadList)
            tsToProcessDict = {tsId: ts.clone() for ts in inTsSet.iterItems()
                               if (tsId := ts.getTsId()) in nonProcessedTsIds  # Only not processed tsIds
                               and ts.getSize() > 0}  # Avoid processing empty TS
            for tsId, ts in tsToProcessDict.items():
                xDim = ts.getXDim()
                if tomoWidth > xDim:
                    tomoWidth = 0
                    widthWarnTsIds.append(tsId)
                cId = self._insertFunctionStep(self.convertInStep,
                                               ts,
                                               prerequisites=[],
                                               needsGPU=False)
                recId = self._insertFunctionStep(self.computeReconstructionStep,
                                                 tsId,
                                                 tomoWidth,
                                                 prerequisites=cId,
                                                 needsGPU=True)
                cOutId = self._insertFunctionStep(self.createOutputStep,
                                                  ts,
                                                  prerequisites=recId,
                                                  needsGPU=False)
                closeSetStepDeps.append(cOutId)
                logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                self.tsReadList.append(tsId)

            self.refreshStreaming(inTsSet)

    # --------------------------- STEPS functions -----------------------------
    def convertInStep(self, ts: TiltSeries):
        super().convertInputStep(ts, presentAcqOrders=ts.getTsPresentAcqOrders())

    def computeReconstructionStep(self, tsId: str, tomoWidth: int):
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'===> tsId = {tsId}: reconstructing the tomogram...'))
                outTsFn, outTsOddFn, outTsEvenFn = self.getTmpFileNames(tsId)
                # run tilt
                paramsTilt = {
                    "-InputProjections": outTsFn,
                    "-OutputFile": self.getTmpOutFile(tsId, ext=MRC_EXT),
                    "-TILTFILE": self.getExtraOutFile(tsId, ext=TLT_EXT),
                    "-THICKNESS": self.tomoThickness.get(),
                    "-FalloffIsTrueSigma": 1,
                    "-RADIAL": f"{self.radialFirstParameter.get()},{self.radialSecondParameter.get()}",
                    "-SHIFT": f"{self.tomoShiftX.get()},{self.tomoShiftZ.get()}",
                    "-OFFSET": f"{self.angleOffset.get()},{self.tiltAxisOffset.get()}",
                    "-MODE": 1,
                    "-PERPENDICULAR": "",
                    "-AdjustOrigin": "",
                    "SuperSampleFactor": self.superSampleFactor.get()
                }

                if self.fakeInteractionsSIRT.get() != 0:
                    paramsTilt["-FakeSIRTiterations"] = self.fakeInteractionsSIRT.get()

                if self.usesGpu():
                    paramsTilt["-UseGPU"] = self.getGpuList()[0]
                    paramsTilt["-ActionIfGPUFails"] = "2,2"

                self.runProgram(TILT_PROGRAM, paramsTilt)

                # run trimvol
                trimVolOpts = "-rx "
                if tomoWidth > 0:
                    trimVolOpts += f" -nx {tomoWidth}"

                paramsTrimVol = {
                    'input': self.getTmpOutFile(tsId, ext=MRC_EXT),
                    'output': self.getExtraOutFile(tsId, ext=MRC_EXT),
                    'options': trimVolOpts
                }

                argsTrimvol = "%(options)s %(input)s %(output)s"
                Plugin.runImod(self, TRIMVOL_PROGRAM, argsTrimvol % paramsTrimVol)

                oddEvenTmp = [[], []]
                if self.doOddEven:
                    # Odd
                    paramsTilt['-InputProjections'] = outTsOddFn
                    oddEvenTmp[0] = self.getTmpOutFile(tsId, suffix=ODD, ext=MRC_EXT)
                    paramsTilt['-OutputFile'] = oddEvenTmp[0]
                    self.runProgram(TILT_PROGRAM, paramsTilt)

                    paramsTrimVol['input'] = oddEvenTmp[0]
                    paramsTrimVol['output'] = self.getExtraOutFile(tsId, suffix=ODD, ext=MRC_EXT)
                    Plugin.runImod(self, TRIMVOL_PROGRAM, argsTrimvol % paramsTrimVol)

                    # Even
                    paramsTilt['-InputProjections'] = outTsEvenFn
                    oddEvenTmp[1] = self.getTmpOutFile(tsId, suffix=EVEN, ext=MRC_EXT)
                    paramsTilt['-OutputFile'] = oddEvenTmp[1]
                    self.runProgram(TILT_PROGRAM, paramsTilt)

                    paramsTrimVol['input'] = oddEvenTmp[1]
                    paramsTrimVol['output'] = self.getExtraOutFile(tsId, suffix=EVEN, ext=MRC_EXT)
                    Plugin.runImod(self, TRIMVOL_PROGRAM, argsTrimvol % paramsTrimVol)

            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> {TILT_PROGRAM} or {TRIMVOL_PROGRAM} execution '
                                    f'failed with the exception -> {e}'))
                logger.error(traceback.format_exc())

    def createOutputStep(self, ts: TiltSeries):
        tsId = ts.getTsId()
        if tsId in self.failedItems:
            self.addToOutFailedSet(ts)
            return

        try:
            outputFn = self.getExtraOutFile(tsId, ext=MRC_EXT)
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
            # Set of tomograms
            outTomoSet = self.getOutputSetOfTomograms(self.getInputTsSet(pointer=True))
            # Tomogram
            outTomo = Tomogram(tsId=tsId)
            outTomo.copyInfo(ts)
            outTomo.setFileName(outputFn)
            self.setTomoOddEven(tsId, outTomo)
            # Set default tomogram origin
            outTomo.setOrigin(newOrigin=None)
            shiftX = self.tomoShiftX.get()
            shiftZ = self.tomoShiftZ.get()
            if shiftX or shiftZ:
                sRate = outTomo.getSamplingRate()
                x, y, z = outTomo.getShiftsFromOrigin()
                shiftXang = shiftX * sRate
                shiftZang = shiftZ * sRate
                outTomo.setShiftsInOrigin(x=x - shiftXang, y=y, z=z - shiftZang)
            # Data persistence
            outTomoSet.append(outTomo)
            outTomoSet.update(outTomo)
            outTomoSet.write()
            self._store(outTomoSet)
            # Close explicitly the outputs (for streaming)
            self.closeOutputsForStreaming()

    # --------------------------- UTILS functions ---------------------------

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        output = getattr(self, OUTPUT_TOMOGRAMS_NAME, None)
        if output is not None:
            summary.append(f"Input tilt-series: {self.getInputTsSet().getSize()}\n"
                           f"Tomograms reconstructed: {output.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        output = getattr(self, OUTPUT_TOMOGRAMS_NAME, None)
        if output is not None:
            methods.append("The reconstruction has been computed for "
                           f"{output.getSize()} tilt-series using "
                           "the IMOD *tilt* command.")
        return methods
