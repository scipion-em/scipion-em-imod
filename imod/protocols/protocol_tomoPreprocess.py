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
from tomo.objects import Tomogram, SetOfTomograms

from imod.protocols.protocol_base import ProtImodBase
from imod.constants import OUTPUT_TOMOGRAMS_NAME, MRC_EXT, ODD, EVEN


class ProtImodTomoNormalization(ProtImodBase):
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

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTomograms',
                      params.PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Input set of tomograms')

        form.addParam('binning',
                      params.IntParam,
                      default=1,
                      label='Binning',
                      important=True,
                      help='Binning is an scaling factor for the output tomograms. '
                           'Must be an integer greater than 1. IMOD uses ordinary'
                           'binning to reduce tomograms in size by the given factor. '
                           'The value of a binned pixel is the average of pixel '
                           'values in each block of pixels being binned. Binning '
                           'is applied before all')

        form.addParam('floatDensities',
                      params.EnumParam,
                      choices=['No adjust',
                               'range between min and max',
                               'scaled to common mean and standard deviation',
                               'shifted to a common mean without scaling',
                               'shifted to mean and rescaled to a min and max'],
                      default=2,
                      label='Adjust densities mode',
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='Adjust densities of sections individually:\n'
                           '-Default: no adjustment performed\n'
                           
                           '-Range between min and max: This option will scale the gray values'
                           'to be in a range given by a minimum and a maximum values.'
                           'This is the mode 1 in newstack flag -floatDensities.\n'
                           
                           '-Scaled to common mean and standard deviation: This is the most '
                           'common normalization procedure. The new tomograms will have'
                           'a mean and a standard deviation introduced by the user. Generaly,'
                           'a zero mean and a standard deviation one is a good choice.'
                           'This is the mode 2 in newstack flag -floatDensities.\n'
                           
                           '-Shifted to a common mean without scaling: This option only'
                           'add an offset to the gray values of the tomograms. The offset will'
                           'be calculated such as the new tomograms will present a mean gray value'
                           'introduced by the user. This is the mode 3 in newstack flag '
                           'floatDensities.\n'
                           
                           '-shifted to mean and rescaled to a min and max: In this case, an '
                           'offset is added to the tomograms in order to achieve a mean gray value'
                           ' then they are rescale the resulting minimum and maximum densities '
                           'to the Min and Max values specified. This is the mode 4 in newstack'
                           ' flag -floatDensities.\n')

        groupMeanSd = form.addGroup('Mean and SD',
                                    condition='floatDensities==2',
                                    help='Scale all images to the given mean '
                                         'and standard deviation. This option '
                                         'implies -float 2 and is incompatible '
                                         'with all other scaling options. If no '
                                         'values are set, mean=0 and SD=1 by default')

        groupMeanSd.addParam('meanSdToggle',
                             params.EnumParam,
                             choices=['Yes', 'No'],
                             default=1,
                             label='Set mean and SD?',
                             display=params.EnumParam.DISPLAY_HLIST,
                             help='Set mean and SD values')

        groupMeanSd.addParam('scaleMean',
                             params.FloatParam,
                             default=0,
                             label='Mean',
                             help='Mean value for the rescaling')

        groupMeanSd.addParam('scaleSd',
                             params.FloatParam,
                             default=1,
                             label='SD',
                             help='Standard deviation value for the rescaling')

        groupScale = form.addGroup('Scaling values',
                                   condition='floatDensities==4')

        groupScale.addParam('scaleMax',
                            params.FloatParam,
                            default=255,
                            label='Max.',
                            help='Maximum value for the rescaling')

        groupScale.addParam('scaleMin',
                            params.FloatParam,
                            default=0,
                            label='Min.',
                            help='Minimum value for the rescaling')

        form.addParam('modeToOutput',
                      params.EnumParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      choices=['default', '4-bit', '8-bit', 'signed 16-bit',
                               'unsigned 16-bit', '32-bit float'],
                      default=0,
                      label='Storage data type',
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='The storage mode of the output file. The '
                           'default is the mode of the first input file, '
                           'except for a 4-bit input file, where the default '
                           'is to output as bytes')

        form.addParam('scaleRangeToggle',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      condition="floatDensities==1 or floatDensities==3",
                      default=0,
                      label='Set scaling range values?',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='This option will rescale the densities of all '
                           'sections by the same factors so that the original '
                           'minimum and maximum density will be mapped '
                           'to the Min and Max values that are entered')

        form.addParam('scaleRangeMax',
                      params.FloatParam,
                      condition="(floatDensities==1 or floatDensities==3) and scaleRangeToggle==0",
                      default=255,
                      label='Max.',
                      help='Maximum value for the rescaling')

        form.addParam('scaleRangeMin',
                      params.FloatParam,
                      condition="(floatDensities==1 or floatDensities==3) and scaleRangeToggle==0",
                      default=0,
                      label='Min.',
                      help='Minimum value for the rescaling')

        form.addParam('antialias',
                      params.EnumParam,
                      choices=['None', 'Blackman', 'Triangle', 'Mitchell',
                               'Lanczos 2 lobes', 'Lanczos 3 lobes'],
                      default=5,
                      label='Antialias method:',
                      expertLevel=params.LEVEL_ADVANCED,
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='Type of antialiasing filter to use when reducing images.\n'
                           'The available types of filters are:\n\n'
                           'None - Antialias will not be applied\n'
                           'Blackman - fast but not as good at antialiasing as slower filters\n'
                           'Triangle - fast but smooths more than Blackman\n'
                           'Mitchell - good at antialiasing, smooths a bit\n'
                           'Lanczos 2 lobes - good at antialiasing, less smoothing than Mitchell\n'
                           'Lanczos 3 lobes - slower, even less smoothing but more risk of ringing\n'
                           'The default is Lanczos 3 as of IMOD 4.7. Although '
                           'many people consider Lanczos 2 the best compromise '
                           'among the various factors, that sentiment may be '
                           'based on images of natural scenes where there are '
                           'sharp edges.')

        form.addParam('processOddEven',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=True,
                      label='Process odd/even?',
                      help='If True, the full tilt series and the associated odd/even '
                           'tilt series will be processed. The transformations applied '
                           'to the odd/even tilt series will be exactly the same.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        binning = self.binning.get()
        norm = self.floatDensities.get()
        runNewStack = norm != 0 or (self.getModeToOutput() is not None)

        for tsId in self.tomoDict.keys():
            self._insertFunctionStep(self.preprocessStep,
                                     tsId, runNewStack, binning)
            self._insertFunctionStep(self.generateOutputStep,
                                     tsId, runNewStack, binning)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.tomoDict = {tomo.getTsId(): tomo.clone() for tomo in self.inputSetOfTomograms.get()}

    def preprocessStep(self, tsId, runNewstack, binning):
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

        if norm != 0:
            paramsNewstack["-FloatDensities"] = norm

            if norm == 2:
                if self.meanSdToggle.get() == 0:
                    paramsNewstack["-MeanAndStandardDeviation"] = f"{self.scaleMean.get()}, {self.scaleSd.get()}"

            elif norm == 4:
                paramsNewstack["-ScaleMinAndMax"] = f"{self.scaleMax.get()}, {self.scaleMin.get()}"

            else:
                if self.scaleRangeToggle.get() == 0:
                    paramsNewstack["-ScaleMinAndMax"] = f"{self.scaleRangeMax.get()}, {self.scaleRangeMin.get()}"

        if self.getModeToOutput() is not None:
            paramsNewstack["-ModeToOutput"] = self.getModeToOutput()

        oddEvenOutput = [[], []]

        if runNewstack:
            self.runProgram("newstack", paramsNewstack)

            if self.applyToOddEven(tomo):
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

                if self.applyToOddEven(tomo):
                    inputOdd, inputEven = (self.getTmpOutFile(tsId, suffix=ODD, ext=MRC_EXT),
                                           self.getTmpOutFile(tsId, suffix=EVEN, ext=MRC_EXT))
                    pwpath.moveFile(oddEvenOutput[0], inputOdd)
                    pwpath.moveFile(oddEvenOutput[1], inputEven)

            else:
                inputTomoPath = tomoFn
                if self.applyToOddEven(tomo):
                    inputOdd, inputEven = tomo.getHalfMaps().split(',')

            paramsBinvol = {
                '-input': inputTomoPath,
                '-output': outputFile,
                '-binning': binning,
                '-antialias': self.antialias.get() + 1
            }

            self.runProgram('binvol', paramsBinvol)

            if self.applyToOddEven(tomo):
                paramsBinvol['-input'] = inputOdd
                paramsBinvol['-output'] = self.getExtraOutFile(tsId, suffix=ODD, ext=MRC_EXT)
                self.runProgram('binvol', paramsBinvol)

                paramsBinvol['-input'] = inputEven
                paramsBinvol['-output'] = self.getExtraOutFile(tsId, suffix=EVEN, ext=MRC_EXT)
                self.runProgram('binvol', paramsBinvol)

    def generateOutputStep(self, tsId, runNewstack, binning):
        tomo = self.tomoDict[tsId]
        output = self.getOutputSetOfTomograms(self.inputSetOfTomograms.get(),
                                              binning)
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
        else:
            newTomogram.copyAttributes(tomo, '_origin')

        if self.applyToOddEven(tomo):
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
            summary.append(f"Input tilt-series: {self.inputSetOfTomograms.get().getSize()}\n"
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

    def applyToOddEven(self, tomo):
        """ Reimplemented from base class for the tomogram case. """
        return self.processOddEven and tomo.hasHalfMaps()