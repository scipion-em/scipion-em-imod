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

import os

from pyworkflow import BETA
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from tomo.objects import Tomogram, SetOfTomograms

from .. import Plugin
from .protocol_base import ProtImodBase, EXT_MRC_ODD_NAME, EXT_MRC_EVEN_NAME, OUTPUT_TOMOGRAMS_NAME


class ProtImodTomoNormalization(ProtImodBase):
    """
    Normalize input tomogram and change its storing formatting.
    More info:
        https://bio3D.colorado.edu/imod/doc/newstack.html
        https://bio3D.colorado.edu/imod/doc/binvol.html
    """

    _label = 'Tomo preprocess'
    _devStatus = BETA
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
                      help='Binning to be applied to the normalized tomograms '
                           'in IMOD convention. Volumes will be binned by the '
                           'given factor. Must be an integer bigger than 1')

        form.addParam('floatDensities',
                      params.EnumParam,
                      choices=['default', '1', '2', '3', '4'],
                      default=0,
                      label='Adjust densities mode',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Adjust densities of sections individually:\n'
                           '-Default: no adjustment performed\n'
                           '-Mode 1: sections fill the data range\n'
                           '-Mode 2: sections scaled to common mean and standard deviation.\n'
                           '-Mode 3: sections shifted to a common mean without scaling\n'
                           '-Mode 4: sections shifted to a common mean and then '
                           'rescale the resulting minimum and maximum densities '
                           'to the Min and Max values specified')
        form.addParam('modeToOutput',
                      params.EnumParam,
                      choices=['default', '4-bit', 'byte', 'signed 16-bit',
                               'unsigned 16-bit', '32-bit float'],
                      default=0,
                      label='Storage data type',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Apply one density scaling to all sections to map '
                           'current min and max to the given Min and Max. The '
                           'storage mode of the output file. The default is '
                           'the mode of the first input file, except for a '
                           '4-bit input file, where the default is to output '
                           'as bytes')

        form.addParam('scaleRangeToggle',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      condition="floatDensities==0 or floatDensities==1 or floatDensities==3",
                      default=1,
                      label='Set scaling range values?',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='This option will rescale the densities of all '
                           'sections by the same factors so that the original '
                           'minimum and maximum density will be mapped '
                           'to the Min and Max values that are entered')

        form.addParam('scaleRangeMax',
                      params.FloatParam,
                      condition="(floatDensities==0 or floatDensities==1 or floatDensities==3) and scaleRangeToggle==0",
                      default=255,
                      label='Max.',
                      help='Maximum value for the rescaling')

        form.addParam('scaleRangeMin',
                      params.FloatParam,
                      condition="(floatDensities==0 or floatDensities==1 or floatDensities==3) and scaleRangeToggle==0",
                      default=0,
                      label='Min.',
                      help='Minimum value for the rescaling')

        form.addParam('antialias',
                      params.EnumParam,
                      choices=['None', 'Blackman', 'Triangle', 'Mitchell',
                               'Lanczos 2', 'Lanczos 3'],
                      default=5,
                      label='Antialias method:',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Type of antialiasing filter to use when reducing images.\n'
                           'The available types of filters are:\n\n'
                           'None\n'
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

        form.addParam('processOddEven',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=True,
                      label='Process odd/even?',
                      help='If True, the full tilt series and the associated odd/even tilt series will be processed. '
                           'The transformations applied to the odd/even tilt series will be exactly the same.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        for tomo in self.inputSetOfTomograms.get():
            self._insertFunctionStep(self.generateOutputStackStep,
                                     tomo.getObjId())
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def generateOutputStackStep(self, tsObjId):
        tomo = self.inputSetOfTomograms.get()[tsObjId]
        location = tomo.getFileName()
        fileName = os.path.splitext(location)[0]

        extraPrefix = self._getExtraPath(os.path.basename(fileName))
        tmpPrefix = self._getTmpPath(os.path.basename(fileName))
        path.makePath(extraPrefix)
        path.makePath(tmpPrefix)

        runNewstack = False

        paramsNewstack = {
            'input': location,
            'output': os.path.join(extraPrefix, os.path.basename(location)),
            'imagebinned': 1.0,
        }

        argsNewstack = "-input %(input)s " \
                       "-output %(output)s " \
                       "-imagebinned %(imagebinned)s "

        if self.floatDensities.get() != 0:
            runNewstack = True
            argsNewstack += " -FloatDensities " + str(self.floatDensities.get())

            if self.floatDensities.get() == 2:
                if self.meanSdToggle.get() == 0:
                    argsNewstack += " -MeanAndStandardDeviation " + str(self.scaleMean.get()) + "," + \
                                    str(self.scaleSd.get())

            elif self.floatDensities.get() == 4:
                argsNewstack += " -ScaleMinAndMax " + str(self.scaleMax.get()) + "," + str(self.scaleMin.get())

            else:
                if self.scaleRangeToggle.get() == 0:
                    argsNewstack += " -ScaleMinAndMax " + str(self.scaleRangeMax.get()) + "," + \
                                    str(self.scaleRangeMin.get())

        if self.getModeToOutput() is not None:
            runNewstack = True
            argsNewstack += " -ModeToOutput " + str(self.getModeToOutput())

        oddEvenOutput = [[], []]

        if runNewstack:
            Plugin.runImod(self, 'newstack', argsNewstack % paramsNewstack)

            if self.applyToOddEven(tomo):
                oddFn, evenFn = tomo.getHalfMaps().split(',')
                paramsNewstack['input'] = oddFn
                oddEvenOutput[0] = os.path.join(extraPrefix, tomo.getTsId() + EXT_MRC_ODD_NAME)
                paramsNewstack['output'] = oddEvenOutput[0]
                Plugin.runImod(self, 'newstack', argsNewstack % paramsNewstack)
                paramsNewstack['input'] = evenFn
                oddEvenOutput[1] = os.path.join(extraPrefix, tomo.getTsId() + EXT_MRC_EVEN_NAME)
                paramsNewstack['output'] = oddEvenOutput[1]
                Plugin.runImod(self, 'newstack', argsNewstack % paramsNewstack)

        binning = self.binning.get()

        if binning != 1:
            if runNewstack:
                baseLoc = os.path.basename(location)
                tmpPath = os.path.join(tmpPrefix, baseLoc)
                path.moveFile(os.path.join(extraPrefix, baseLoc), tmpPath)
                inputTomoPath = tmpPath

                if self.applyToOddEven(tomo):
                    path.moveFile(oddEvenOutput[0], tmpPath)
                    path.moveFile(oddEvenOutput[1], tmpPath)
                    inputTomoPath = tmpPath
                    inputOdd, inputEven = (os.path.join(tmpPrefix, tomo.getTsId() + EXT_MRC_ODD_NAME),
                                           os.path.join(tmpPrefix, tomo.getTsId() + EXT_MRC_EVEN_NAME))
            else:
                inputTomoPath = location
                if self.applyToOddEven(tomo):
                    inputOdd, inputEven = tomo.getHalfMaps().split(',')


            paramsBinvol = {
                'input': inputTomoPath,
                'output': os.path.join(extraPrefix, os.path.basename(location)),
                'binning': binning,
                'antialias': self.antialias.get() + 1
            }

            argsBinvol = "-input %(input)s " \
                         "-output %(output)s " \
                         "-binning %(binning)d "\
                         "-antialias %(antialias)d "

            Plugin.runImod(self, 'binvol', argsBinvol % paramsBinvol)

            if self.applyToOddEven(tomo):
                paramsBinvol['input'] = inputOdd
                paramsBinvol['output'] = os.path.join(extraPrefix, tomo.getTsId() + EXT_MRC_ODD_NAME)
                Plugin.runImod(self, 'binvol', argsBinvol % paramsBinvol)
                paramsBinvol['input'] = inputEven
                paramsBinvol['output'] = os.path.join(extraPrefix, tomo.getTsId() + EXT_MRC_EVEN_NAME)
                Plugin.runImod(self, 'binvol', argsBinvol % paramsBinvol)

        output = self.getOutputSetOfTomograms(self.inputSetOfTomograms.get(),
                                              binning)

        newTomogram = Tomogram()
        newTomogram.copyInfo(tomo)
        newTomogram.setTsId(tomo.getTsId())

        if not runNewstack and binning == 1:
            newTomogram.setLocation(location)
        else:
            location = os.path.join(extraPrefix, os.path.basename(location))
            newTomogram.setLocation(location)

        if binning > 1:
            sr = tomo.getSamplingRate() * binning

            newTomogram.setSamplingRate(sr)

            # Set default tomogram origin
            newTomogram.setOrigin(newOrigin=None)

        else:
            newTomogram.copyAttributes(tomo, '_origin')

        if self.applyToOddEven(tomo):
            halfMapsList = [os.path.join(extraPrefix, tomo.getTsId() + EXT_MRC_ODD_NAME),
                            os.path.join(extraPrefix, tomo.getTsId() + EXT_MRC_EVEN_NAME)]
            newTomogram.setHalfMaps(halfMapsList)

        output.append(newTomogram)
        output.updateDim()
        output.update(newTomogram)

        output.write()
        self._store()

    def closeOutputSetsStep(self):
        self.Tomograms.setStreamState(Set.STREAM_CLOSED)
        self.Tomograms.write()
        self._store()

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

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.Tomograms:
            summary.append("Input tilt-series: %d\nInterpolations applied: %d"
                           % (self.inputSetOfTomograms.get().getSize(),
                              self.Tomograms.getSize()))
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.Tomograms:
            methods.append("%d tomograms have been normalized using the """
                           "IMOD *binvol* command.\n"
                           % (self.Tomograms.getSize()))
        return methods

    def applyToOddEven(self, tomo):
        """ Reimplemented from base class for the tomogram case. """
        return self.processOddEven and tomo.hasHalfMaps()
