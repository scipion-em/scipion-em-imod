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
from pwem.convert.headers import setMRCSamplingRate
from pyworkflow.protocol import STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.utils import Message, cyanStr, redStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from imod.protocols import ProtImodBase
from imod.constants import OUTPUT_TILTSERIES_NAME, IMODFINDBEADS_PROGRAM, CCDERASER_PROGRAM, ODD, EVEN

logger = logging.getLogger(__name__)


class ProtImodFiducialEraser(ProtImodBase, ProtStreamingBase):
    """
    Removes fiducial gold beads from tilt-series images using IMOD-based
    image restoration procedures. The protocol generates corrected tilt
    images in which fiducial markers are erased to reduce reconstruction
    artifacts and improve downstream tomographic analysis. More info:
    https://bio3d.colorado.edu/imod/doc/man/imodfindbeads.html
    https://bio3d.colorado.edu/imod/doc/man/ccderaser.html

    AI Generated:

    Fiducial Eraser (ProtImodFiducialEraser) - User Manual
        Overview

        The Fiducial Eraser protocol removes fiducial gold beads from
        electron tomography tilt-series after alignment procedures have
        been completed. In cryo-electron tomography workflows, fiducial
        markers are essential for accurate alignment because they provide
        stable reference points across tilted views. However, once the
        alignment stage is finished, these beads may become undesirable
        because they introduce strong artificial densities into the final
        reconstruction.

        This protocol is particularly useful when preparing datasets for
        segmentation, visualization, subtomogram analysis, or publication
        figures where fiducial markers could obscure biologically relevant
        structures. By replacing the fiducial regions with locally modeled
        image information, the protocol produces cleaner tilt-series that
        preserve the biological specimen while minimizing contamination
        from gold particles.

        Biological Context and Motivation

        Fiducial gold beads are commonly deposited onto cryo-EM grids to
        facilitate tilt-series alignment. Although extremely valuable during
        motion correction and geometric alignment, they often become a source
        of reconstruction artifacts due to their high contrast and strong
        scattering signal. These artifacts may propagate into tomograms as
        streaks or bright densities that interfere with interpretation.

        For biological studies involving membranes, cytoskeletal assemblies,
        ribosomes, or crowded intracellular environments, the presence of
        fiducials can complicate visualization and automated analysis. Their
        removal is therefore recommended before performing segmentation,
        particle picking, denoising, or machine-learning-based analyses.

        Inputs and General Workflow

        The protocol requires an aligned tilt-series as input. It first
        identifies fiducial gold beads and then removes them from the
        images using interpolation and local image reconstruction methods.
        The corrected tilt-series preserves the original acquisition geometry
        while replacing bead-containing regions with smoothly reconstructed
        surrounding signal.

        The workflow is designed to operate automatically and can process
        multiple tilt-series in streaming environments. This makes the
        protocol particularly useful in facility-scale cryo-ET pipelines
        where data are continuously acquired and processed.

        Fiducial Detection

        The first stage identifies fiducial markers based on their expected
        diameter. Accurate specification of fiducial size is biologically
        important because it determines which image features are interpreted
        as gold beads. If the diameter is underestimated, portions of the
        fiducials may remain visible after correction. If it is overestimated,
        surrounding biological signal may be unnecessarily modified.

        In most practical situations, the fiducial diameter should correspond
        to the nominal bead size used during sample preparation. Typical
        cryo-ET experiments employ fiducials ranging from approximately
        5 to 20 nanometers depending on the specimen and microscope settings.

        Erasing the Fiducials

        After fiducials are identified, the protocol removes them from the
        images using interpolation procedures designed to preserve local
        image continuity. The erased region is generally chosen to be slightly
        larger than the bead itself to ensure that residual high-density
        signal is fully eliminated.

        From a biological perspective, careful selection of the erasing
        diameter is important. Excessively small values may leave visible
        fiducial remnants, while excessively large values may remove nearby
        structural information from the specimen. In crowded cellular
        environments, conservative values are often preferable to minimize
        alteration of neighboring densities.

        Odd and Even Tilt-Series Processing

        In workflows involving odd-even tilt-series separation, the protocol
        can process both subsets independently. This is particularly relevant
        for advanced tomographic reconstruction strategies that estimate
        reproducibility, resolution, or noise properties by comparing
        independently processed halves of the dataset.

        Maintaining consistent fiducial removal across odd and even datasets
        helps preserve compatibility with downstream validation procedures
        and subtomogram refinement workflows.

        Outputs and Their Interpretation

        The protocol produces corrected tilt-series in which fiducial markers
        have been removed while preserving the acquisition geometry and image
        organization of the original data. The resulting datasets can be used
        directly for tomographic reconstruction and downstream biological
        interpretation.

        The corrected images should be visually inspected to ensure that
        fiducials have been adequately removed without introducing obvious
        interpolation artifacts. Small residual traces may remain in difficult
        regions, particularly when fiducials overlap dense biological
        structures or lie near image boundaries.

        Practical Recommendations

        In routine cryo-ET practice, it is generally advisable to remove
        fiducials only after alignment quality has been verified. Fiducials
        should remain available during troubleshooting and alignment
        optimization because they provide essential geometric references.

        When selecting the erasing diameter, users should begin with values
        slightly larger than the fiducial diameter and inspect the results
        visually. Dense cellular samples may require more conservative
        settings to avoid altering nearby biological signal.

        For datasets intended for segmentation, visualization, or machine
        learning applications, fiducial erasure often provides a substantial
        improvement in image interpretability and reduces the risk of
        artificial feature detection.

        Final Perspective

        Fiducial removal is an important finishing step in many electron
        tomography workflows. Although fiducials are indispensable during
        alignment, their persistence in reconstructed tomograms may interfere
        with biological interpretation and computational analysis. Careful
        removal of these markers helps produce cleaner and more biologically
        meaningful datasets suitable for visualization, reconstruction, and
        quantitative structural analysis.
    """

    _label = 'fiducial eraser'
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
        form.addParam('fidDiameter',
                      params.FloatParam,
                      important=True,
                      default= 10.0,
                      label='Fiducial diameter (nm)')
        form.addParam('erasedDiameter',
                      params.FloatParam,
                      important=True,
                      default=12.0,
                      label='Diameter to erase (nm)')
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
                cInId = self._insertFunctionStep(self.linkTsStep,
                                                 ts,
                                                 prerequisites=[],
                                                 needsGPU=False)
                beadsId = self._insertFunctionStep(self.imodfindbeadsStep,
                                                   ts,
                                                   prerequisites=cInId,
                                                   needsGPU=False)
                eraserId = self._insertFunctionStep(self.ccderaserStep,
                                                    ts,
                                                    prerequisites=beadsId,
                                                    needsGPU=False)
                createOutputId = self._insertFunctionStep(self.createOutputStep,
                                                          ts,
                                                          prerequisites=eraserId,
                                                          needsGPU=False)
                closeSetStepDeps.append(createOutputId)
                logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                self.tsIdReadList.append(tsId)

                self.refreshStreaming(inTsSet)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        super()._initialize()

    def imodfindbeadsStep(self, ts: TiltSeries):
        """This step creates a fiducial model"""
        tsId = ts.getTsId()
        try:
            logger.info(cyanStr(f'tsId = {tsId} -> Finding the fiducials...'))
            firstTi = ts.getFirstItem()

            paramsImodFindBeads = {
                "-inp": firstTi.getFileName(),
                "-o": self._getModelFileName(tsId),
                "-size": self._getFiducialDiameterPx(ts)
            }
            self.runProgram(IMODFINDBEADS_PROGRAM, paramsImodFindBeads)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {IMODFINDBEADS_PROGRAM} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def ccderaserStep(self, ts: TiltSeries):
        """This step erase the gold beads from the fiducial model"""
        tsId = ts.getTsId()
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'tsId = {tsId} -> Erasing the gold beads...'))
                firstTi = ts.getFirstItem()

                paramsCCDeraser = {
                    "-ExpandCircleIterations": 3,
                    "-BetterRadius": self._getFiducialDiameterPx(ts) / 2,
                    "-input": firstTi.getFileName(),
                    "-output": self.getExtraOutFile(tsId),
                    "-ModelFile": self._getModelFileName(tsId),
                    "-MergePatches": '',
                    "-ExcludeAdjacent": '',
                    "-CircleObjects": '/',
                    "-SkipTurnedOffPoints": 1,
                    "-PolynomialOrder": -1,
                }

                self.runProgram(CCDERASER_PROGRAM, paramsCCDeraser)

                if self.doOddEven:
                    logger.info(cyanStr(f'tsId = {tsId} -> Erasing the gold beads (ODD Tilt-series) ...'))
                    paramsCCDeraser['-input'] = ts.getOddFileName(),
                    paramsCCDeraser['-output'] = self.getExtraOutFile(tsId, suffix=ODD)
                    self.runProgram(CCDERASER_PROGRAM, paramsCCDeraser)

                    logger.info(cyanStr(f'tsId = {tsId} -> Erasing the gold beads (EVEN Tilt-series) ...'))
                    paramsCCDeraser['-input'] = ts.getEvenFileName(),
                    paramsCCDeraser['-output'] = self.getExtraOutFile(tsId, suffix=EVEN)
                    self.runProgram(CCDERASER_PROGRAM, paramsCCDeraser)

            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> {CCDERASER_PROGRAM} execution failed'
                                    f' with the exception -> {e}'))
                logger.error(traceback.format_exc())

    def createOutputStep(self, ts: TiltSeries):
        tsId = ts.getTsId()
        if tsId in self.failedItems:
            self.addToOutFailedSet(ts)
            return

        try:
            outTsFile = self.getExtraOutFile(tsId)
            if exists(outTsFile):
                self._registerOutput(ts, outTsFile)
            else:
                logger.error(redStr(f'tsId = {tsId} -> Output file {outTsFile} was not generated. Skipping... '))

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    @retry_on_sqlite_lock(log=logger)
    def _registerOutput(self, ts: TiltSeries, outTsFile: str):
        tsId = ts.getTsId()
        with self._lock:
            # Set of tilt-series
            outTsSet = self.getOutputSetOfTS(self.getInputTsSet(pointer=True))
            # Tilt-series
            outTs = TiltSeries()
            outTs.copyInfo(ts)
            outTsSet.append(outTs)
            setMRCSamplingRate(outTsFile, ts.getSamplingRate())  # Update the apix value in file header
            # Tilt-images
            for ti in ts.iterItems():
                outTi = TiltImage()
                outTi.copyInfo(ti)
                outTi.setFileName(outTsFile)
                self.setTsOddEven(tsId, outTi, binGenerated=True)
                outTs.append(outTi)
            # Data persistence
            outTs.write()
            outTsSet.update(outTs)
            outTsSet.write()
            self._store(outTsSet)
            # Close explicitly the outputs (for streaming)
            self.closeOutputsForStreaming()

    # --------------------------- INFO functions ------------------------------

    # --------------------------- UTILS functions -----------------------------
    def _getFiducialDiameterPx(self, ts: TiltSeries) -> float:
        return 10 * self.fidDiameter.get() / ts.getSamplingRate()

    def _getModelFileName(self, tsId: str) -> str:
        return self._getExtraPath(tsId, f'{tsId}_model')