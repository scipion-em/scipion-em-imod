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
from os.path import exists
import pyworkflow.protocol.params as params
from imod.protocols.protocol_base import IN_TOMO_SET, NEWSTACK_PROGRAM
from imod.protocols.protocol_base_preprocess import ProtImodBasePreprocess
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.utils import Message, moveFile, cyanStr, redStr
from tomo.objects import Tomogram, SetOfTomograms
from imod.constants import OUTPUT_TOMOGRAMS_NAME, MRC_EXT, ODD, EVEN, BINVOL_PROGRAM

logger = logging.getLogger(__name__)


class ProtImodTomoNormalization(ProtImodBasePreprocess):
    """
    Normalize input tomogram and change its storing formatting.

    More info:
        https://bio3D.colorado.edu/imod/doc/newstack.html
        https://bio3d.colorado.edu/imod/doc/man/binvol.html

    IMOD tilt series preprocess makes use of the Newstack and
    binvol commands. In particular, three functionalities are possible:\n

    _1 Binning_: Binvol will bin down a volume in all three dimensions,
    with the binning done isotropically. Binning means summing (actually
    averaging) all the values in a block of voxels (e.g., 2x2x2
    or 1x1x3) in the input volume to create one voxel in the output volume.
    The output file will have appropriately larger pixel spacings
    in its header.\n
    _2 Normalization_: This protocol allows to scale the gray values
    of the tomograms, also called normalization, to a common range or
    mean of density. The most used normalization consists in zero
    mean and standard deviation one.\n

    _3 storage format_: IMOD is able to modify the number of bit of
    the stored data in order to reduce the disc occupancy.
    """

    _label = 'Tomo preprocess'
    _possibleOutputs = {OUTPUT_TOMOGRAMS_NAME: SetOfTomograms}
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form, *args):
        form.addSection(Message.LABEL_INPUT)
        form.addParam(IN_TOMO_SET,
                      params.PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Input set of tomograms')
        super()._defineParams(form, isTomogram=True)
        form.addParallelSection(threads=1, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        binning = self.binning.get()
        norm = self.floatDensities.get()
        runNewStack = norm != 0 or (self.getModeToOutput() is not None)
        closeSetStepDeps = []

        for tomo in self.getInputTomoSet():
            tsId = tomo.getTsId()
            compId = self._insertFunctionStep(self.preprocessStep,
                                              tsId,
                                              runNewStack,
                                              binning,
                                              prerequisites=[],
                                              needsGPU=False)
            outId = self._insertFunctionStep(self.generateOutputStep,
                                             tsId,
                                             runNewStack,
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
        self.doOddEven = self.applyToOddEven(self.getInputTomoSet())

    def preprocessStep(self, tsId: str,
                       runNewstack: bool,
                       binning: int):
        try:
            logger.info(cyanStr(f'tsId = {tsId} -> Preprocessing...'))
            self.genTsPaths(tsId)
            with self._lock:
                tomo = self.getCurrentTomo(tsId)
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

    def generateOutputStep(self, tsId: str,
                           runNewstack: bool,
                           binning: int):
        if tsId in self.failedItems:
            self.addToOutFailedSet(tsId, inputsAreTs=False)
        else:
            try:
                outputFn = self.getExtraOutFile(tsId, ext=MRC_EXT)
                if exists(outputFn):
                    with self._lock:
                        tomo = self.getCurrentTomo(tsId)
                        # Set of tomograms
                        outTomoSet = self.getOutputSetOfTomograms(self.getInputTomoSet(pointer=True),
                                                                  binning=binning)
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
                        if self.doOddEven:
                            halfMapsList = [self.getExtraOutFile(tsId, suffix=ODD, ext=MRC_EXT),
                                            self.getExtraOutFile(tsId, suffix=EVEN, ext=MRC_EXT)]
                            outTomo.setHalfMaps(halfMapsList)
                        # Data persistence
                        outTomoSet.append(outTomo)
                        outTomoSet.updateDim()
                        outTomoSet.update(outTomo)
                        outTomoSet.write()
                        self._store(outTomoSet)
                else:
                    logger.error(redStr(f'tsId = {tsId} -> Output file {outputFn} was not generated. Skipping... '))
            except Exception as e:
                logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))

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
            oddFn, evenFn = tomo.getHalfMaps().split(',')
            # Odd
            paramsNewstack['-input'] = oddFn
            oddEvenOutput[0] = self.getExtraOutFile(tsId, suffix=ODD, ext=MRC_EXT)
            paramsNewstack['-output'] = oddEvenOutput[0]
            self.runProgram(NEWSTACK_PROGRAM, paramsNewstack)
            # Even
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
        oddEvenOutput = [[], []]
        logger.info(cyanStr(f'tsId = {tsId}. Executing {BINVOL_PROGRAM}:.'))

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

        else:
            tomoFn = tomo.getFileName()
            inputTomoPath = tomoFn
            if self.doOddEven:
                inputOdd, inputEven = tomo.getHalfMaps().split(',')

        paramsBinvol = {
            '-input': inputTomoPath,
            '-output': outputFile,
            '-binning': binning,
            '-antialias': self.antialias.get() + 1
        }

        self.runProgram(BINVOL_PROGRAM, paramsBinvol)

        if self.doOddEven:
            # Odd
            paramsBinvol['-input'] = inputOdd
            paramsBinvol['-output'] = self.getExtraOutFile(tsId, suffix=ODD, ext=MRC_EXT)
            self.runProgram(BINVOL_PROGRAM, paramsBinvol)
            # Even
            paramsBinvol['-input'] = inputEven
            paramsBinvol['-output'] = self.getExtraOutFile(tsId, suffix=EVEN, ext=MRC_EXT)
            self.runProgram(BINVOL_PROGRAM, paramsBinvol)
