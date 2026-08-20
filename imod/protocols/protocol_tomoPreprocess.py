# *****************************************************************************
# *
# * Authors:     Scipion Team [1]
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
from os.path import exists
from imod.protocols.protocol_base import NEWSTACK_PROGRAM
from imod.protocols.protocol_base_preprocess import ProtImodBasePreprocess
from pwem.convert.headers import setMRCSamplingRate
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.utils import Message, moveFile, cyanStr, redStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import Tomogram, SetOfTomograms
from imod.constants import OUTPUT_TOMOGRAMS_NAME, MRC_EXT, ODD, EVEN, BINVOL_PROGRAM

logger = logging.getLogger(__name__)


class ProtImodTomoNormalization(ProtImodBasePreprocess):
    """
    Preprocesses tomograms by applying normalization, isotropic binning,
    and storage format conversion in order to optimize tomographic datasets
    for visualization, analysis, reconstruction refinement, and downstream
    cryo-electron tomography workflows. The protocol is designed to improve
    data consistency, reduce storage requirements, and adapt tomograms to
    computational or biological requirements while preserving the essential
    structural information contained in the original reconstructions. More info:
    https://bio3D.colorado.edu/imod/doc/newstack.html
    https://bio3d.colorado.edu/imod/doc/man/binvol.html

    AI Generated:

    Tomogram Normalization and Preprocessing (ProtImodTomoNormalization) —
    User Manual

        Overview

        The Tomogram Normalization and Preprocessing protocol prepares one or
        more tomograms for subsequent cryo-electron tomography analysis by
        applying common preprocessing operations such as intensity
        normalization, isotropic binning, and output format conversion. These
        operations are essential in many biological workflows because raw or
        reconstructed tomograms often differ in voxel size, density scale,
        storage precision, or computational requirements.

        In practical cryo-ET projects, preprocessing is commonly performed
        before segmentation, subtomogram averaging, visualization, denoising,
        machine learning analysis, or quantitative comparison between datasets.
        Proper preprocessing improves consistency across experiments and helps
        ensure that downstream procedures operate under stable and biologically
        meaningful conditions.

        Inputs and General Workflow

        The protocol requires one or more tomograms as input. Each tomogram is
        processed independently, allowing large collections of reconstructed
        volumes to be standardized in a reproducible manner. Depending on the
        selected configuration, the protocol can normalize voxel intensities,
        reduce tomogram size through isotropic binning, and modify the storage
        precision of the output files.

        These preprocessing operations can be combined in a single workflow,
        making the protocol suitable both for rapid exploratory processing and
        for more demanding production-level cryo-electron tomography pipelines.

        Binning and Resolution Considerations

        One of the central preprocessing operations is isotropic binning.
        Binning reduces the dimensions of the tomogram by averaging neighboring
        voxels into larger sampling units. This operation decreases file size,
        reduces memory consumption, and accelerates downstream computations.

        From a biological perspective, binning is often useful during early
        exploratory analysis, template matching, manual inspection, or
        segmentation of large cellular tomograms. Moderate binning can improve
        signal-to-noise characteristics and facilitate interpretation of broad
        structural features. However, aggressive binning inevitably reduces
        high-resolution information and may obscure fine molecular details.

        A common practical strategy is to begin processing with moderately
        binned tomograms for rapid analysis and later return to higher
        resolution datasets for detailed structural interpretation or
        publication-quality processing.

        Density Normalization

        Density normalization adjusts voxel intensity distributions so that
        tomograms share a more consistent statistical scale. This is
        particularly important when combining datasets originating from
        different microscopes, acquisition conditions, reconstruction
        procedures, or preprocessing pipelines.

        The most common normalization strategy in cryo-electron tomography is
        scaling the data toward zero mean and unit variance. This standard
        normalization improves compatibility with many downstream computational
        procedures, including machine learning approaches, template matching,
        subtomogram classification, and denoising algorithms.

        Alternative scaling approaches may also be useful when preparing data
        for visualization or for software environments expecting specific
        intensity ranges. Biological users should nevertheless remain aware
        that normalization changes the numerical interpretation of voxel values,
        even though it does not alter the underlying structural content.

        Storage Format Optimization

        Tomographic datasets can occupy very large amounts of disk space,
        particularly in large-scale cellular tomography studies. The protocol
        therefore allows modification of the storage precision used to encode
        voxel densities. Reducing the storage precision can substantially lower
        storage and transfer requirements, especially when handling large
        collections of tomograms.

        In many biological workflows, reduced precision formats are adequate
        for visualization, exploratory analysis, or intermediate processing.
        However, when quantitative analysis or high-resolution interpretation
        is required, users should carefully evaluate whether reduced precision
        could affect subtle density variations.

        Odd and Even Dataset Handling

        The preprocessing workflow can preserve independent odd and even
        tomographic datasets when such information is available. Maintaining
        these separated datasets is important in workflows involving resolution
        validation, noise estimation, subtomogram averaging, or independent
        refinement procedures.

        Preserving odd and even volumes throughout preprocessing helps ensure
        that downstream validation strategies remain statistically meaningful.

        Outputs and Their Interpretation

        The protocol produces a new set of tomograms containing the selected
        preprocessing modifications. Depending on the configuration, the
        outputs may differ from the original datasets in voxel size, intensity
        scale, storage precision, or overall dimensions.

        Despite these transformations, the biological identity and acquisition
        metadata of each tomogram are preserved, allowing the outputs to remain
        fully compatible with downstream cryo-electron tomography procedures.

        Practical Recommendations

        For large tomograms or cellular datasets, moderate isotropic binning is
        often an effective first step because it substantially reduces
        computational cost while preserving overall structural organization.
        Users performing exploratory analysis or segmentation commonly benefit
        from working initially with binned datasets.

        Standardized normalization is generally recommended when combining
        tomograms from multiple experiments or when using machine learning and
        quantitative analysis methods. Consistent intensity scaling improves
        comparability and algorithmic stability.

        Storage precision reduction should be applied carefully. While it may
        dramatically reduce disk usage, excessive reduction in precision can
        potentially limit downstream quantitative interpretation.

        Final Perspective

        Tomogram preprocessing is not merely a technical preparation stage but
        an important component of reliable cryo-electron tomography analysis.
        Appropriate normalization, binning, and storage optimization can
        greatly improve computational efficiency and workflow consistency while
        preserving the biological interpretability of the reconstructed
        structures. Thoughtful preprocessing choices help balance data quality,
        computational cost, and downstream analytical goals.
    """

    _label = 'tomo preprocess'
    _possibleOutputs = {OUTPUT_TOMOGRAMS_NAME: SetOfTomograms}
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tomoDict = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form, *args):
        form.addSection(Message.LABEL_INPUT)
        super().addInTomoSetFormParam(form)
        super()._defineParams(form, isTomogram=True)
        form.addParallelSection(threads=1, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        binning = self.binning.get()
        norm = self.floatDensities.get()
        runNewStack = norm != 0 or (self.getModeToOutput() is not None)
        closeSetStepDeps = []

        for tsId in self.tomoDict.keys():
            compId = self._insertFunctionStep(self.preprocessStep,
                                              tsId,
                                              runNewStack,
                                              binning,
                                              prerequisites=[],
                                              needsGPU=False)
            outId = self._insertFunctionStep(self.generateOutputStep,
                                             tsId,
                                             binning,
                                             prerequisites=compId,
                                             needsGPU=False)
            closeSetStepDeps.append(outId)

        self._insertFunctionStep(self.closeOutputSetsStep,
                                 OUTPUT_TOMOGRAMS_NAME,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        inTomosSet = self.getInputTomoSet()
        self.doOddEven = self.applyToOddEven(inTomosSet)
        self.tomoDict = {tomo.getTsId(): tomo.clone() for tomo in inTomosSet.iterItems()}

    def preprocessStep(self, tsId: str,
                       runNewstack: bool,
                       binning: int):
        try:
            self.genTsPaths(tsId)
            tomo = self.tomoDict[tsId]
            outputFile = self.getExtraOutFile(tsId, ext=MRC_EXT)

            # Operations that don't concern the 3D binning
            if runNewstack:
               self._runNewstack(tomo, outputFile)

            # Run binvol for 3D binning (more precise for volumes than newstack)
            if binning != 1:
                self._runBinvol(tomo, outputFile, binning, runNewstack)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {NEWSTACK_PROGRAM} or {BINVOL_PROGRAM} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def generateOutputStep(self,
                           tsId: str,
                           binning: int):
        tomo = self.tomoDict[tsId]
        if tsId in self.failedItems:
            self.addToOutFailedSet(tomo)
            return

        try:
            outputFn = self.getExtraOutFile(tsId, ext=MRC_EXT)
            if exists(outputFn):
                self._registerOutput(tomo, binning, outputFn)
            else:
                logger.error(redStr(f'tsId = {tsId} -> Output file {outputFn} was not generated. Skipping... '))

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    @retry_on_sqlite_lock(log=logger)
    def _registerOutput(self, tomo: Tomogram, binning: int, outputFn: str):
        tsId = tomo.getTsId()
        with self._lock:
            # Set of tomograms
            outTomoSet = self.getOutputSetOfTomograms(self.getInputTomoSet(pointer=True),
                                                      binning=binning)
            setMRCSamplingRate(outputFn, outTomoSet.getSamplingRate())  # Update the apix value in file header
            # Tomogram
            outTomo = Tomogram(tsId=tsId)
            outTomo.copyInfo(tomo)
            outTomo.setFileName(outputFn)
            if binning > 1:
                outTomo.setSamplingRate(tomo.getSamplingRate() * binning)
                # Fix the mrc tomogram
                outTomo.fixMRCVolume(setSamplingRate=True)
                # Set default tomogram origin
                outTomo.setOrigin(newOrigin=None)
            self.setTomoOddEven(tsId, outTomo)
            # Data persistence
            outTomoSet.append(outTomo)
            outTomoSet.updateDim()
            outTomoSet.update(outTomo)
            outTomoSet.write()
            self._store(outTomoSet)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        output = getattr(self, OUTPUT_TOMOGRAMS_NAME, None)
        if output is not None:
            summary.append(f"Input tilt-series: {self.getInputTomoSet().getSize()}\n"
                           f"Interpolations applied: {output.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        output = getattr(self, OUTPUT_TOMOGRAMS_NAME, None)
        if output is not None:
            methods.append(f"{output.getSize()} tomograms have been "
                           "preprocessed using the IMOD *binvol* command.")
        return methods

    def _validate(self):
        errorMsg = []
        binning = self.binning.get()
        norm = self.floatDensities.get()
        runNewStack = norm != 0 or (self.getModeToOutput() is not None)
        if binning == 1 and not runNewStack:
            errorMsg.append('No operation will be carried out with the parameters introduced.')
        return errorMsg

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

    def _runNewstack(self,
                     tomo: Tomogram,
                     outputFile: str) -> None:
        tsId = tomo.getTsId()
        oddEvenOutput = [[], []]
        logger.info(cyanStr(f'tsId = {tsId}: Preprocessing...'))

        norm = self.floatDensities.get()
        # Important not to bin at this step
        paramsNewstack = self.getBasicNewstackParams(tomo,
                                                     tomo.getFileName(),
                                                     outputFile,
                                                     doNorm=norm != 0)

        paramsNewstack["-antialias"] = self.antialias.get() + 1
        # Float densities
        if norm > 0:
            paramsNewstack["-FloatDensities"] = norm
            if norm == 2:
                paramsNewstack["-MeanAndStandardDeviation"] = f"{self.scaleMean.get()},{self.scaleSd.get()}"
            elif norm == 4:
                paramsNewstack["-ScaleMinAndMax"] = f"{self.scaleMax.get()},{self.scaleMin.get()}"

        if self.getModeToOutput() is not None:
            paramsNewstack["-ModeToOutput"] = self.getModeToOutput()
        self.runProgram(NEWSTACK_PROGRAM, paramsNewstack)

        if self.doOddEven:
            evenFn, oddFn = sorted(tomo.getHalfMaps(asList=True))
            # Odd
            logger.info(cyanStr(f'tsId = {tsId} ODD: Preprocessing...'))
            paramsNewstack['-input'] = oddFn
            oddEvenOutput[0] = self.getExtraOutFile(tsId, suffix=ODD, ext=MRC_EXT)
            paramsNewstack['-output'] = oddEvenOutput[0]
            self.runProgram(NEWSTACK_PROGRAM, paramsNewstack)
            # Even
            logger.info(cyanStr(f'tsId = {tsId} EVEN: Preprocessing...'))
            paramsNewstack['-input'] = evenFn
            oddEvenOutput[1] = self.getExtraOutFile(tsId, suffix=EVEN, ext=MRC_EXT)
            paramsNewstack['-output'] = oddEvenOutput[1]
            self.runProgram(NEWSTACK_PROGRAM, paramsNewstack)

    def _runBinvol(self,
                   tomo: Tomogram,
                   outputFile: str,
                   binning: int,
                   runNewstack: bool) -> None:
        tsId = tomo.getTsId()
        logger.info(cyanStr(f'tsId = {tsId}: Executing {BINVOL_PROGRAM}:.'))
        oddEvenOutput = [[], []]
        inputEven, inputOdd = sorted(tomo.getHalfMaps(asList=True)) if self.doOddEven else None, None
        inputTomoPath = tomo.getFileName()

        if runNewstack:
            # Move previous outputs to tmp as intermediate results and use them as inputs for binvol
            tmpPath = self.getTmpOutFile(tsId, ext=MRC_EXT)
            moveFile(outputFile, tmpPath)
            inputTomoPath = tmpPath
            if self.doOddEven:
                inputOdd, inputEven = (self.getTmpOutFile(tsId, suffix=ODD, ext=MRC_EXT),
                                       self.getTmpOutFile(tsId, suffix=EVEN, ext=MRC_EXT))
                moveFile(oddEvenOutput[0], inputOdd)
                moveFile(oddEvenOutput[1], inputEven)

        paramsBinvol = {
            '-input': inputTomoPath,
            '-output': outputFile,
            '-binning': binning,
            '-antialias': self.antialias.get() + 1
        }

        self.runProgram(BINVOL_PROGRAM, paramsBinvol)

        if self.doOddEven:
            # Odd
            logger.info(cyanStr(f'tsId = {tsId} ODD: Executing {BINVOL_PROGRAM}:.'))
            paramsBinvol['-input'] = inputOdd
            paramsBinvol['-output'] = self.getExtraOutFile(tsId, suffix=ODD, ext=MRC_EXT)
            self.runProgram(BINVOL_PROGRAM, paramsBinvol)
            # Even
            logger.info(cyanStr(f'tsId = {tsId} EVEN: Executing {BINVOL_PROGRAM}:.'))
            paramsBinvol['-input'] = inputEven
            paramsBinvol['-output'] = self.getExtraOutFile(tsId, suffix=EVEN, ext=MRC_EXT)
            self.runProgram(BINVOL_PROGRAM, paramsBinvol)
