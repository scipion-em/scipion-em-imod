# *****************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
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
from imod.protocols.protocol_base import NEWSTACK_PROGRAM
from imod.protocols.protocol_base_preprocess import ProtImodBasePreprocess
from pwem.convert.headers import setMRCSamplingRate
from pyworkflow.protocol import STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.utils import Message, cyanStr, redStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import SetOfTiltSeries, TiltImage, TiltSeries
from imod.constants import OUTPUT_TILTSERIES_NAME, ODD, EVEN

logger = logging.getLogger(__name__)


class ProtImodTsNormalization(ProtImodBasePreprocess, ProtStreamingBase):
    """
    Preprocesses cryo-electron tomography tilt-series by applying
    normalization, binning, and storage format conversion using IMOD
    utilities. The protocol is intended to standardize tilt-series
    before downstream reconstruction, alignment, denoising, or
    subtomogram analysis workflows. More info:
        https://bio3d.colorado.edu/imod/doc/man/newstack.html

    AI Generated:

    Tilt-Series Preprocess (ProtImodTsNormalization) - User Manual
        Overview

        The Tilt-Series Preprocess protocol provides a unified
        environment for preparing aligned or partially processed
        tilt-series before subsequent tomographic analysis. Its main
        objective is to improve dataset consistency, reduce storage
        requirements, and optimize image quality for later stages of
        cryo-electron tomography workflows. In practice, this protocol
        is commonly used immediately after motion correction or tilt
        alignment and before reconstruction, CTF correction, particle
        picking, or subtomogram averaging.

        The protocol combines three major preprocessing operations:
        binning, normalization of image intensities, and modification
        of storage precision. These operations are frequently required
        in facility pipelines and large-scale tomography projects
        where datasets may originate from different acquisition
        sessions, microscopes, or processing conditions.

        Inputs and General Workflow

        The protocol operates on a set of tilt-series and produces a
        new processed tilt-series dataset suitable for downstream
        tomographic procedures. The resulting data preserve the
        geometric and acquisition relationships of the original
        dataset while adapting the sampling, intensity scaling, and
        storage representation to the desired processing conditions.

        In biological practice, preprocessing is often one of the most
        influential stages for determining reconstruction quality.
        Properly normalized and appropriately binned tilt-series
        generally improve stability in later alignment and
        reconstruction steps, particularly when handling noisy
        cryo-electron tomography data acquired under low-dose
        conditions.

        Binning and Sampling Considerations

        Binning reduces the dimensions of the tilt-series images while
        preserving the dominant structural information. This operation
        is especially important in cryo-electron tomography because
        raw datasets are frequently very large and computationally
        demanding. Reducing image size decreases memory usage,
        accelerates reconstruction, and simplifies exploratory
        analysis.

        From a biological perspective, binning represents a compromise
        between resolution and computational efficiency. Low binning
        factors preserve more structural detail and are generally
        preferred for high-resolution subtomogram averaging or
        structural interpretation. Higher binning factors reduce noise
        and improve computational speed, making them useful during
        early exploratory processing, rapid alignment procedures, or
        large screening projects.

        In many workflows, users first reconstruct heavily binned
        tomograms for quick inspection and later return to minimally
        binned data for final refinement and interpretation.

        Image Normalization

        The normalization functionality adjusts image intensity
        distributions to achieve more homogeneous density statistics
        across the tilt-series. This step is biologically important
        because cryo-electron tomography datasets often contain
        substantial variations in illumination, detector response,
        exposure conditions, or ice thickness.

        Standard normalization strategies commonly aim to produce
        images with zero-centered densities and controlled variance.
        Such normalization improves numerical stability during
        alignment, reconstruction, and machine learning applications.
        It also facilitates visual interpretation by reducing strong
        intensity disparities between projections.

        Different normalization modes may be appropriate depending on
        the biological objective. Conservative normalization is often
        sufficient for standard reconstruction workflows, whereas more
        carefully controlled scaling may be beneficial when datasets
        are intended for quantitative analysis, denoising pipelines,
        or neural-network-based segmentation approaches.

        Biological users should remember that normalization changes
        image intensity statistics but does not recover missing
        structural information. Excessive scaling or inappropriate
        normalization ranges may amplify artifacts or suppress weak
        biological signal.

        Storage Format Optimization

        The protocol also allows modification of the numerical storage
        format used for the processed tilt-series. This capability is
        particularly useful for reducing disk usage and improving data
        management efficiency in large tomography facilities or
        high-throughput projects.

        Lower-precision storage formats may significantly decrease
        storage requirements while maintaining sufficient information
        for many intermediate processing tasks. However, aggressive
        reduction of numerical precision can limit downstream
        quantitative analyses or reduce the ability to recover weak
        high-resolution features.

        In practice, reduced storage precision is commonly acceptable
        for exploratory visualization, rapid alignment, or temporary
        intermediate datasets, whereas high-precision formats are
        generally preferred for final reconstruction and
        publication-quality analyses.

        Interpolation and Geometric Consistency

        During preprocessing, geometric consistency between tilt
        images is preserved so that the processed tilt-series remain
        compatible with downstream alignment and reconstruction
        procedures. This is particularly important when tilt-series
        already contain alignment information or transformations
        generated during previous processing stages.

        In cryo-electron tomography, maintaining correct geometric
        relationships is essential because even small inconsistencies
        may propagate into reconstruction artifacts, reduced
        resolution, or inaccurate subtomogram localization.

        Odd and Even Tilt-Series Processing

        The protocol can also operate on odd and even subsets of the
        tilt-series when these datasets are available. This feature is
        especially valuable in workflows involving resolution
        estimation, denoising validation, or independent half-set
        processing strategies.

        Processing odd and even datasets consistently ensures that
        downstream validation procedures remain statistically reliable
        and biologically meaningful.

        Streaming and High-Throughput Workflows

        The protocol is designed to support streaming execution,
        allowing preprocessing to begin while data acquisition or
        upstream processing is still ongoing. This capability is
        particularly useful in automated cryo-electron tomography
        facilities where rapid feedback and continuous data handling
        are essential.

        Streaming workflows enable preprocessing to proceed in
        parallel with acquisition, thereby reducing idle computation
        time and accelerating the transition toward reconstruction and
        interpretation stages.

        Outputs and Their Interpretation

        The protocol produces a new set of processed tilt-series with
        updated sampling properties, normalized image statistics, and
        optionally modified storage precision. These outputs preserve
        the biological content of the original data while adapting the
        datasets for improved computational handling and downstream
        analysis.

        Biologically, preprocessing does not create new structural
        information. Instead, it improves consistency, interpretability,
        and computational accessibility of the acquired signal.

        Practical Recommendations

        For most routine cryo-electron tomography workflows, moderate
        binning combined with standard normalization provides a good
        balance between computational efficiency and preservation of
        structural detail. Exploratory reconstruction and alignment
        procedures are often substantially faster on preprocessed
        datasets.

        Users pursuing high-resolution subtomogram averaging should
        avoid excessive binning and aggressive precision reduction,
        particularly during final refinement stages. Conversely,
        strongly binned datasets are highly useful for rapid visual
        inspection, segmentation prototyping, and quality control.

        When handling heterogeneous datasets collected under varying
        imaging conditions, normalization is particularly important
        for ensuring stable downstream behavior and reducing unwanted
        variability unrelated to biological structure.

        Final Perspective

        In modern cryo-electron tomography workflows, preprocessing is
        not simply a technical convenience but an essential stage that
        strongly influences computational efficiency, reconstruction
        quality, and interpretability. Appropriate choices for
        binning, normalization, and storage precision should always
        reflect the biological question, the expected resolution, and
        the computational resources available for the project.
    """

    _label = 'Tilt-series preprocess'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @classmethod
    def worksInStreaming(cls):
        return True

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form, *args):
        form.addSection(Message.LABEL_INPUT)
        super().addInTsSetFormParam(form)
        super()._defineParams(form)

    # -------------------------- INSERT steps functions -----------------------
    def stepsGeneratorStep(self) -> None:
        """
        This step should be implemented by any streaming protocol.
        It should check its input and when ready conditions are met
        call the self._insertFunctionStep method.
        """
        self._initialize()
        closeSetStepDeps = []
        binning = self.binning.get()
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
                convId = self._insertFunctionStep(self.linkTsStep,
                                                  ts,
                                                  prerequisites=[],
                                                  needsGPU=False)
                compId = self._insertFunctionStep(self.generateOutputStackStep,
                                                  ts,
                                                  binning,
                                                  prerequisites=convId,
                                                  needsGPU=False)
                outId = self._insertFunctionStep(self.createOutputStep,
                                                 ts,
                                                 binning,
                                                 prerequisites=compId,
                                                 needsGPU=False)
                closeSetStepDeps.append(outId)
                logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                self.tsIdReadList.append(tsId)

            self.refreshStreaming(inTsSet)


    # --------------------------- STEPS functions -----------------------------
    def generateOutputStackStep(self, ts: TiltSeries, binning: int):
        tsId = ts.getTsId()
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'===> tsId = {tsId}: preprocessing...'))
                norm = self.floatDensities.get()
                paramsDict = self.getBasicNewstackParams(ts,
                                                         self.getTmpOutFile(tsId),
                                                         self.getExtraOutFile(tsId),
                                                         binning=binning,
                                                         doNorm=norm != 0)

                paramsDict["-antialias"] = self.antialias.get() + 1
                # Float densities
                if norm > 0:
                    paramsDict["-FloatDensities"] = norm
                    if norm == 2:
                        paramsDict["-MeanAndStandardDeviation"] = f"{self.scaleMean.get()},{self.scaleSd.get()}"
                    elif norm == 4:
                        paramsDict["-ScaleMinAndMax"] = f"{self.scaleMax.get()},{self.scaleMin.get()}"

                if self.getModeToOutput() is not None:
                    paramsDict["-ModeToOutput"] = self.getModeToOutput()

                self.runProgram(NEWSTACK_PROGRAM, paramsDict)

                if self.doOddEven:
                    paramsDict['-input'] = self.getTmpOutFile(tsId, suffix=ODD)
                    paramsDict['-output'] = self.getExtraOutFile(tsId, suffix=ODD)
                    self.runProgram(NEWSTACK_PROGRAM, paramsDict)

                    paramsDict['-input'] = self.getTmpOutFile(tsId, suffix=EVEN)
                    paramsDict['-output'] = self.getExtraOutFile(tsId, suffix=EVEN)
                    self.runProgram(NEWSTACK_PROGRAM, paramsDict)

            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> {NEWSTACK_PROGRAM} execution '
                                    f'failed with the exception -> {e}'))
                logger.error(traceback.format_exc())

    def createOutputStep(self, ts: TiltSeries, binning: int):
        tsId = ts.getTsId()
        if tsId in self.failedItems:
            self.addToOutFailedSet(ts)
            return

        try:
            outputFn = self.getExtraOutFile(tsId)
            if not exists(outputFn):
                logger.error(redStr(f'tsId = {tsId} -> Output file {outputFn} was not generated. Skipping... '))
                return

            samplingRate = self.getInputTsSet().getSamplingRate()
            if binning > 1:
                samplingRate *= binning
            setMRCSamplingRate(outputFn, samplingRate)  # Update the apix value in file header
            self._registerOutput(ts, binning, outputFn)

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    @retry_on_sqlite_lock(log=logger)
    def _registerOutput(self, ts: TiltSeries, binning: int, outputFn: str):
        tsId = ts.getTsId()
        with self._lock:
            outTsSet = self.getOutputSetOfTS(self.getInputTsSet(pointer=True), binning)
            outTs = TiltSeries()
            outTs.copyInfo(ts)
            outTsSet.append(outTs)
            for ti in ts.iterItems(orderBy=TiltImage.INDEX_FIELD):
                outTi = TiltImage()
                outTi.copyInfo(ti)
                outTi.setFileName(outputFn)
                self.updateTransformMatrix(outTi, binning=binning)
                self.setTsOddEven(tsId, outTi, binGenerated=True)
                outTs.append(outTi)
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
                           f"Interpolations applied: {output.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            methods.append(f"{output.getSize()} tilt-series have been "
                           "normalized using the IMOD *newstack* command.")
        return methods

    # --------------------------- UTILS functions -----------------------------
    def getModeToOutput(self):
        parseParamsOutputMode = {
            0: None,
            1: 101,
            2: 0,
            3: 1,
            4: 6,
            5: 2
        }
        return parseParamsOutputMode[self.modeToOutput.get()]

    @staticmethod
    def updateTransformMatrix(ti: TiltImage, binning: int = 1) -> None:
        if ti.hasTransform() and binning != 1:
            transform = ti.getTransform()
            matrix = transform.getMatrix()

            matrix[0][2] /= binning
            matrix[1][2] /= binning

            transform.setMatrix(matrix)
            ti.setTransform(transform)


