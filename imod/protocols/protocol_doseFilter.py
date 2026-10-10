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
from pwem.convert.headers import setMRCSamplingRate
from pyworkflow.protocol import STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.utils import Message, cyanStr, redStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from imod.protocols import ProtImodBase
from imod.constants import (ODD, EVEN, SCIPION_IMPORT, FIXED_DOSE,
                            OUTPUT_TILTSERIES_NAME, MTTFILTER_PROGRAM)

logger = logging.getLogger(__name__)


class ProtImodDoseFilter(ProtImodBase, ProtStreamingBase):
    """
    Applies dose-dependent filtering to cryo-electron tomography tilt-series
    using the IMOD dose-weighting strategy. The protocol compensates for
    radiation damage accumulated during acquisition in order to preserve
    high-resolution information and improve the interpretability of tomographic
    reconstructions.
    More info:
        https: // bio3d.colorado.edu / imod / doc / man / mtffilter.html

    AI Generated:

    Dose Filter (ProtImodDoseFilter) — User Manual
        Overview

        The Dose Filter protocol performs radiation damage compensation on
        cryo-electron tomography tilt-series using the IMOD implementation of
        dose weighting. In cryo-EM experiments, biological samples are
        progressively damaged as electrons accumulate during imaging. High
        spatial frequencies are particularly sensitive to this effect, meaning
        that later images in a tilt-series often contain less reliable
        structural information than earlier ones.

        The purpose of this protocol is to reduce the negative impact of
        accumulated dose by attenuating damaged frequency components according
        to experimentally validated dose-response models. The resulting
        tilt-series are therefore more suitable for downstream tomographic
        reconstruction, subtomogram averaging, particle detection, and
        structural interpretation.

        In practical biological workflows, dose filtering is considered a
        standard preprocessing step for modern cryo-electron tomography data,
        especially when working with high-resolution subtomogram averaging or
        structurally heterogeneous samples. Proper dose compensation improves
        contrast preservation and helps maintain meaningful high-frequency
        signal throughout the reconstruction pipeline.

        Biological Context and Importance

        Radiation damage is one of the central limitations in cryo-electron
        microscopy. As the specimen receives increasing electron exposure,
        delicate structural details become progressively degraded. This effect
        is cumulative across a tilt-series because each projection contributes
        additional dose to the same biological specimen.

        Dose filtering addresses this problem by weighting the contribution of
        each image according to the amount of damage already accumulated at the
        time of acquisition. Earlier tilts generally preserve more structural
        detail, while later tilts contribute proportionally less at high
        resolution.

        For biological users, this correction is especially important when
        analyzing macromolecular complexes, membrane proteins, viral particles,
        or cellular environments where subtle structural features are required
        for interpretation. Without dose weighting, reconstructions may appear
        blurred, noisy, or depleted in high-resolution detail.

        Inputs and Dose Information

        The protocol requires an input set of tilt-series together with their
        acquisition metadata. The most important information is the electron
        dose associated with each tilt image. Accurate dose information allows
        the protocol to estimate the extent of radiation damage accumulated
        throughout the acquisition process.

        Two main workflows are supported. In the first, dose values are
        imported directly from acquisition metadata already stored in the
        project. This is the preferred approach because it reflects the real
        experimental conditions used during data collection. In the second,
        users may provide a fixed dose value that is assumed to be identical
        for all tilt images. This option is useful when metadata are incomplete
        or unavailable.

        An optional initial dose parameter can also be included to account for
        exposure received before the recorded tilt-series began. This becomes
        relevant in workflows where preview exposures, focusing procedures, or
        other acquisition steps contributed additional irradiation to the
        specimen.

        Dose Weighting Strategy

        The protocol follows a biologically motivated filtering model in which
        the contribution of high-frequency information decreases progressively
        with accumulated exposure. Frequencies that are more sensitive to
        radiation damage are attenuated more strongly as dose increases.

        From a practical perspective, this means that later projections still
        contribute useful low-frequency structural information, while their
        less reliable high-frequency signal is suppressed. The resulting
        balance improves reconstruction quality without discarding valuable
        low-resolution information.

        This strategy is particularly important in tomography because large
        tilt-series may involve substantial cumulative exposure. Biological
        samples with flexible domains, membrane regions, or weakly scattering
        densities often benefit significantly from proper dose compensation.

        Streaming and High-Throughput Processing

        The protocol supports streaming operation, making it suitable for
        facility-scale cryo-ET pipelines and automated acquisition workflows.
        Newly imported tilt-series can be processed continuously as they become
        available, allowing rapid integration into preprocessing pipelines.

        This capability is particularly valuable in modern high-throughput
        cryo-electron tomography facilities, where large numbers of tilt-series
        are collected during extended microscope sessions. Streaming operation
        enables dose correction to occur in parallel with ongoing acquisition,
        reducing turnaround time for reconstruction and analysis.

        Odd and Even Tilt-Series Processing

        The protocol can also process odd and even image subsets separately
        when these data are available. This functionality is especially useful
        for workflows involving resolution estimation, independent half-set
        validation, or subtomogram averaging strategies that rely on split-data
        processing.

        Maintaining separate odd and even dose-weighted outputs allows users to
        preserve statistical independence during downstream refinement and
        validation procedures. This is important for preventing overfitting and
        ensuring reliable structural interpretation.

        Outputs and Their Interpretation

        The protocol produces dose-weighted tilt-series that retain the
        original acquisition geometry while incorporating radiation damage
        compensation. These filtered tilt-series are intended to replace the
        original unfiltered data in subsequent reconstruction workflows.

        The resulting datasets generally exhibit improved preservation of
        biologically meaningful signal, especially at intermediate and high
        spatial frequencies. Reconstructions derived from dose-weighted data
        often display clearer macromolecular boundaries, improved contrast, and
        more reliable structural features.

        Because the dose compensation has already been applied, the output
        metadata are updated accordingly so downstream processing stages treat
        the data as fully corrected inputs.

        Practical Recommendations

        In most biological workflows, using experimentally imported dose
        metadata is strongly recommended because it provides the most accurate
        representation of the acquisition conditions. Fixed-dose approaches are
        acceptable when metadata are unavailable, but they may not capture
        variations introduced during acquisition.

        Users should ensure that microscope voltage and pixel size metadata are
        correct before processing, since these parameters influence the
        filtering behavior and the interpretation of radiation damage.

        Dose filtering is typically performed early in the preprocessing
        workflow, before tomographic reconstruction. Applying it consistently
        across all datasets helps improve comparability between reconstructions
        and enhances the robustness of downstream subtomogram analysis.

        Final Perspective

        For cryo-electron tomography users, dose weighting is more than a
        technical preprocessing correction. It directly influences the quality
        and interpretability of reconstructed biological structures by
        compensating for one of the most fundamental limitations of electron
        microscopy: radiation damage. Accurate dose modeling, reliable
        acquisition metadata, and consistent preprocessing practices are key
        factors for obtaining high-quality tomographic reconstructions suitable
        for meaningful biological interpretation.
    """

    _label = 'Dose filter'
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
        form.addParam('initialDose',
                      params.FloatParam,
                      default=0.0,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Initial dose (e/Å^2)',
                      help='Dose applied before any of the images in the '
                           'input file were taken; this value will be '
                           'added to all the dose values.')
        form.addParam('inputDoseType',
                      params.EnumParam,
                      choices=['Scipion import', 'Fixed dose'],
                      default=SCIPION_IMPORT,
                      label='Input dose source',
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='Where to find the dose information:\n'
                           '- Scipion import: use the dose provided '
                           'during import of the tilt-series\n'
                           '- Fixed dose: manually input fixed dose '
                           'for each image of the input file, '
                           'in electrons/Å^2.')
        form.addParam('fixedImageDose',
                      params.FloatParam,
                      default=FIXED_DOSE,
                      label='Fixed dose (e/Å^2)',
                      condition='inputDoseType == %i' % FIXED_DOSE,
                      help='Fixed dose for each image of the input file, '
                           'in electrons/square Ångstrom.')
        self.addOddEvenParams(form)
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
                compId = self._insertFunctionStep(self.doseFilterStep,
                                                  ts,
                                                  prerequisites=cInId,
                                                  needsGPU=False)
                outId = self._insertFunctionStep(self.createOutputStep,
                                                 ts,
                                                 prerequisites=[compId],
                                                 needsGPU=False)
                closeSetStepDeps.append(outId)
                logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                self.tsIdReadList.append(tsId)

            self.refreshStreaming(inTsSet)

    # --------------------------- STEPS functions -----------------------------
    def doseFilterStep(self, ts: TiltSeries):
        """Apply the dose filter to every tilt series"""
        tsId = ts.getTsId()
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'tsId = {tsId} -> Dose filtering...'))
                firstItem = ts.getFirstEnabledItem()

                progParams = {
                    '-input': self.getTmpOutFile(tsId),
                    '-output': self.getExtraOutFile(tsId),
                    '-PixelSize': ts.getSamplingRate(),
                    '-Voltage': int(ts.getAcquisition().getVoltage()),
                }

                if self.initialDose.get() != 0.0:
                    progParams["-InitialDose"] = self.initialDose.get()

                if self.inputDoseType.get() == SCIPION_IMPORT:
                    outputDoseFilePath = self.getExtraOutFile(tsId, ext="dose")
                    self.generateDoseFile(ts, outputDoseFilePath)
                    progParams["-TypeOfDoseFile"] = 2
                    progParams["-DoseWeightingFile"] = outputDoseFilePath

                elif self.inputDoseType.get() == FIXED_DOSE:
                    progParams["-FixedImageDose"] = self.fixedImageDose.get()

                self.runProgram(MTTFILTER_PROGRAM, progParams)

                if self.doOddEven:
                    # Odd
                    logger.info(cyanStr(f'tsId = {tsId} ODD -> Dose filtering...'))
                    progParams['-input'] = firstItem.getEven()
                    progParams['-output'] = self.getExtraOutFile(tsId, suffix=ODD)
                    self.runProgram(MTTFILTER_PROGRAM, progParams)
                    # Even
                    logger.info(cyanStr(f'tsId = {tsId} EVEN -> Dose filtering...'))
                    progParams['-input'] = firstItem.getOdd()
                    progParams['-output'] = self.getExtraOutFile(tsId, suffix=EVEN)
                    self.runProgram(MTTFILTER_PROGRAM, progParams)

            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> {MTTFILTER_PROGRAM} execution failed '
                                    f'with the exception -> {e}'))
                logger.error(traceback.format_exc())

    def createOutputStep(self, ts: TiltSeries):
        """Generate output filtered tilt series"""
        tsId = ts.getTsId()
        if tsId in self.failedItems:
            self.addToOutFailedSet(ts)
            return

        try:
            outTsFile = self.getExtraOutFile(tsId)
            if not exists(outTsFile):
                logger.error(redStr(f'tsId = {tsId} -> Output file {outTsFile} was not generated. Skipping... '))
                return

            setMRCSamplingRate(outTsFile, ts.getSamplingRate())  # Update the apix value in file header
            self._registerOutput(ts, outTsFile)

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    @retry_on_sqlite_lock(log=logger)
    def _registerOutput(self, ts: TiltSeries, outTsFile: str):
        tsId = ts.getTsId()
        with self._lock:
            # Set of tilt-series
            outTsSet = self.getOutputSetOfTS(self.getInputTsSet(pointer=True))
            outTs = TiltSeries()
            outTs.copyInfo(ts)
            self.updateTsAcquisition(outTs)  # Acquisition dose goes to 0 after having been applied
            outTsSet.append(outTs)
            # Tilt-images
            for ti in ts.iterItems():
                outTi = TiltImage()
                outTi.copyInfo(ti)
                outTi.setFileName(outTsFile)
                self.updateTiAcquisition(outTi)
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
    def _validate(self):
        validateMsgs = []

        if self.inputDoseType.get() == SCIPION_IMPORT:
            for ts in self.getInputTsSet():
                if ts.getFirstEnabledItem().getAcquisition().getDosePerFrame() is None:
                    validateMsgs.append(f"{ts.getTsId()} has no dose information stored "
                                        "in Scipion Metadata. To solve this, re-import "
                                        "tilt-series using the mdoc option.")
                    break

        return validateMsgs

    def _summary(self):
        summary = []

        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            summary.append(f"Input tilt-series: {self.getInputTsSet().getSize()}\n"
                           "Dose weighting applied: "
                           f"{output.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")

        return summary

    def _methods(self):
        methods = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            methods.append("The dose-weighting has been applied to "
                           f"{output.getSize()} "
                           "tilt-series using the IMOD *mtffilter* command.")
        return methods

    # --------------------------- UTILS functions -----------------------------
    @staticmethod
    def updateTiAcquisition(tiOut: TiltImage) -> None:
        """Sets the initial and accumulated doses to 0 for a given tilt-image"""
        # Output is dose-weighted
        acq = tiOut.getAcquisition()
        acq.setDoseInitial(0.)
        acq.setAccumDose(0.)
        tiOut.setAcquisition(acq)

    @staticmethod
    def updateTsAcquisition(tsOut: TiltSeries) -> None:
        """Sets the initial and accumulated doses to 0 for a given tilt-series"""
        # Output is dose-weighted
        acq = tsOut.getAcquisition()
        acq.setAccumDose(0.)
        acq.setDoseInitial(0.)
        tsOut.setAcquisition(acq)

    @staticmethod
    def generateDoseFile(ts: TiltSeries, doseFileOutputPath: str) -> None:
        """ This method generates a file containing the dose information
        of a tilt series in the specified location from the accumulated
        dose and dose per tilt. The format is two columns per each tilt image:
         the prior accumulated dose and the image dose
         """
        doseInfoList = []

        for ti in ts.iterItems(iterate=False):
            acq = ti.getAcquisition()
            doseInfoList.append((acq.getAccumDose() - acq.getDosePerFrame(), acq.getDosePerFrame()))

        np.savetxt(doseFileOutputPath, np.asarray(doseInfoList), fmt='%f', delimiter=" ")
