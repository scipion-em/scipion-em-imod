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

import pyworkflow.protocol.params as params
import pyworkflow.utils.path as pwpath
from imod.protocols.protocol_base import IN_TOMO_SET
from imod.protocols.protocol_base_preprocess import ProtImodBasePreprocess
from pyworkflow.utils import Message
from tomo.objects import Tomogram, SetOfTomograms
from imod.constants import OUTPUT_TOMOGRAMS_NAME, MRC_EXT, ODD, EVEN


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
    averaging) all of the values in a block of voxels (e.g., 2x2x2
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

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        form.addParam(IN_TOMO_SET,
                      params.PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Input set of tomograms')
        super()._defineParams(form, isTomogram=True)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        binning = self.binning.get()
        norm = self.floatDensities.get()
        runNewStack = norm != 0 or (self.getModeToOutput() is not None)
        closeSetStepDeps = []

        for tsId in self.tomoDict.keys():
            compId = self._insertFunctionStep(self.preprocessStep, tsId, runNewStack, binning,
                                              prerequisites=[])
            outId = self._insertFunctionStep(self.generateOutputStep, tsId, runNewStack, binning,
                                             prerequisites=[compId])
            closeSetStepDeps.append(outId)

        self._insertFunctionStep(self.closeOutputSetsStep, prerequisites=closeSetStepDeps)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.tomoDict = {tomo.getTsId(): tomo.clone() for tomo in self.inputSetOfTomograms.get()}
        self.oddEvenFlag = self.applyToOddEven(self.getInputSet())

    def preprocessStep(self, tsId, runNewstack, binning):
        try:
            self.genTsPaths(tsId)
            tomo = self.tomoDict[tsId]
            tomoFn = tomo.getFileName()
            outputFile = self.getExtraOutFile(tsId, ext=MRC_EXT)

            # Run newstack
            norm = self.floatDensities.get()
            # Important not to bin at this step
            paramsNewstack = self.getBasicNewstackParams(tomo, outputFile,
                                                         firstItem=tomo,
                                                         binning=1,
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

            oddEvenOutput = [[], []]

            if runNewstack:
                self.runProgram("newstack", paramsNewstack)

                if self.oddEvenFlag:
                    oddFn, evenFn = tomo.getHalfMaps().split(',')
                    paramsNewstack['-input'] = oddFn
                    oddEvenOutput[0] = self.getExtraOutFile(tsId, suffix=ODD, ext=MRC_EXT)
                    paramsNewstack['-output'] = oddEvenOutput[0]
                    self.runProgram("newstack", paramsNewstack)

                    paramsNewstack['-input'] = evenFn
                    oddEvenOutput[1] = self.getExtraOutFile(tsId, suffix=EVEN, ext=MRC_EXT)
                    paramsNewstack['-output'] = oddEvenOutput[1]
                    self.runProgram("newstack", paramsNewstack)

            # Run binvol
            if binning != 1:
                if runNewstack:
                    # Move previous outputs to tmp
                    tmpPath = self.getTmpOutFile(tsId, ext=MRC_EXT)
                    pwpath.moveFile(outputFile, tmpPath)
                    inputTomoPath = tmpPath

                    if self.oddEvenFlag:
                        inputOdd, inputEven = (self.getTmpOutFile(tsId, suffix=ODD, ext=MRC_EXT),
                                               self.getTmpOutFile(tsId, suffix=EVEN, ext=MRC_EXT))
                        pwpath.moveFile(oddEvenOutput[0], inputOdd)
                        pwpath.moveFile(oddEvenOutput[1], inputEven)

                else:
                    inputTomoPath = tomoFn
                    if self.oddEvenFlag:
                        inputOdd, inputEven = tomo.getHalfMaps().split(',')

                paramsBinvol = {
                    '-input': inputTomoPath,
                    '-output': outputFile,
                    '-binning': binning,
                    '-antialias': self.antialias.get() + 1
                }

                self.runProgram('binvol', paramsBinvol)

                if self.oddEvenFlag:
                    paramsBinvol['-input'] = inputOdd
                    paramsBinvol['-output'] = self.getExtraOutFile(tsId, suffix=ODD, ext=MRC_EXT)
                    self.runProgram('binvol', paramsBinvol)

                    paramsBinvol['-input'] = inputEven
                    paramsBinvol['-output'] = self.getExtraOutFile(tsId, suffix=EVEN, ext=MRC_EXT)
                    self.runProgram('binvol', paramsBinvol)

        except Exception as e:
            self._failedTomos.append(tsId)
            self.error(f'Preprocessing failed for tsId {tsId} -> {e}')

    def generateOutputStep(self, tsId, runNewstack, binning):
        tomo = self.tomoDict[tsId]
        with self._lock:
            if tsId in self._failedItems:
                self.createOutputFailedSet(tomo)
            else:
                output = self.getOutputSetOfTomograms(self.getInputSet(pointer=True), binning)
                newTomogram = Tomogram()
                newTomogram.copyInfo(tomo)

                if not runNewstack and binning == 1:
                    newTomogram.setLocation(tomo.getFileName())
                else:
                    newTomogram.setLocation(self.getExtraOutFile(tsId, ext=MRC_EXT))

                if binning > 1:
                    newTomogram.setSamplingRate(tomo.getSamplingRate() * binning)
                    # Fix the mrc tomogram
                    newTomogram.fixMRCVolume(setSamplingRate=True)
                    # Set default tomogram origin
                    newTomogram.setOrigin(newOrigin=None)

                if self.oddEvenFlag:
                    halfMapsList = [self.getExtraOutFile(tsId, suffix=ODD, ext=MRC_EXT),
                                    self.getExtraOutFile(tsId, suffix=EVEN, ext=MRC_EXT)]
                    newTomogram.setHalfMaps(halfMapsList)

                output.append(newTomogram)
                output.updateDim()
                output.update(newTomogram)
                output.write()
                self._store(output)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.Tomograms:
            summary.append(f"Input tilt-series: {self.getInputSet().getSize()}\n"
                           f"Interpolations applied: {self.Tomograms.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.Tomograms:
            methods.append(f"{self.Tomograms.getSize()} tomograms have been "
                           "preprocessed using the IMOD *binvol* command.")
        return methods

    # --------------------------- UTILS functions -----------------------------
    def getInputSet(self, pointer=False):
        return (self.inputSetOfTomograms.get() if
                not pointer else self.inputSetOfTomograms)

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
