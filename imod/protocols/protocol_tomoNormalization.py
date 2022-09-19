# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
# **************************************************************************

import os
from pyworkflow import BETA
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.objects import Transform
from pwem.emlib.image import ImageHandler
from tomo.objects import Tomogram
from imod import Plugin
from imod.protocols.protocol_base import ProtImodBase


class ProtImodTomoNormalization(ProtImodBase):
    """
    Normalize input tomogram and change its storing formatting.
    More info:
        https://bio3D.colorado.edu/imod/doc/newstack.html
        https://bio3D.colorado.edu/imod/doc/binvol.html
    """

    _label = 'Tomo normalization'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTomograms',
                      params.PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Input set of tomograms')

        form.addParam('binning',
                      params.FloatParam,
                      default=1.0,
                      label='Binning',
                      important=True,
                      help='Binning to be applied to the normalized tomograms in IMOD convention. Volumes will be '
                           'binned by the given factor. Must be an integer bigger than 1')

        form.addParam('floatDensities',
                      params.EnumParam,
                      choices=['default', '1', '2', '3', '4'],
                      default=0,
                      label='Adjust densities mode',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Adjust densities of sections individually:\n'
                           '-Default: no adjustment performed\n'
                           '-Mode 1: sections fill the data range\n'
                           '-Mode 2: sections scaled to common mean and standard deviation.\n'
                           '-Mode 3: sections shifted to a common mean without scaling\n'
                           '-Mode 4: sections shifted to a common mean and then rescale the resulting minimum and '
                           'maximum densities to the Min and Max values specified')

        form.addParam('modeToOutput',
                      params.EnumParam,
                      choices=['default', '4-bit', 'byte', 'signed 16-bit', 'unsigned 16-bit', '32-bit float'],
                      default=0,
                      label='Storage data type',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Apply one density scaling to all sections to map current min and max to the given Min and '
                           'Max. The storage mode of the output file. The default is the mode of the first input file, '
                           'except for a 4-bit input file, where the default is to output as bytes')

        form.addParam('scaleRangeToggle',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      condition="floatDensities==0 or floatDensities==1 or floatDensities==3",
                      default=1,
                      label='Set scaling range values',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='This option will rescale the densities of all sections by the same '
                           'factors so that the original minimum and maximum density will be mapped '
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

        groupMeanSd = form.addGroup('Mean and SD',
                                    condition='floatDensities==2',
                                    help='Scale all images to the given mean and standard deviation. This option '
                                         'implies -float 2 and is incompatible with all other scaling options. If no '
                                         'values are set, mean=0 and SD=1 by default')

        groupMeanSd.addParam('meanSdToggle',
                             params.EnumParam,
                             choices=['Yes', 'No'],
                             default=1,
                             label='Set mean and SD',
                             important=True,
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

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for tomo in self.inputSetOfTomograms.get():
            self._insertFunctionStep(self.generateOutputStackStep, tomo.getObjId())
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions ----------------------------
    def generateOutputStackStep(self, tsObjId):
        tomo = self.inputSetOfTomograms.get()[tsObjId]
        location = tomo.getLocation()[1]
        fileName, fileExtension = os.path.splitext(location)

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

        if runNewstack:
            Plugin.runImod(self, 'newstack', argsNewstack % paramsNewstack)

        if self.binning.get() != 1:
            if runNewstack:
                path.moveFile(os.path.join(extraPrefix, os.path.basename(location)),
                              os.path.join(tmpPrefix, os.path.basename(location)))
                inputTomoPath = os.path.join(tmpPrefix, os.path.basename(location))
            else:
                inputTomoPath = location

            paramsBinvol = {
                'input': inputTomoPath,
                'output': os.path.join(extraPrefix, os.path.basename(location)),
                'binning': self.binning.get(),
            }

            argsBinvol = "-input %(input)s " \
                         "-output %(output)s " \
                         "-binning %(binning)d "

            Plugin.runImod(self, 'binvol', argsBinvol % paramsBinvol)

        output = self.getOutputSetOfTomograms(self.inputSetOfTomograms.get(), self.binning.get())

        newTomogram = Tomogram()
        newTomogram.copyInfo(tomo)
        newTomogram.setTsId(tomo.getTsId())

        if not runNewstack and self.binning.get() == 1:
            newTomogram.setLocation(location)
        else:
            location = os.path.join(extraPrefix, os.path.basename(location))
            newTomogram.setLocation(location)

        if self.binning.get() > 1:
            sr = tomo.getSamplingRate() * int(self.binning.get())

            newTomogram.setSamplingRate(sr)

            # Set default tomogram origin
            newTomogram.setOrigin(tomo.getOrigin(force=True))

        else:
            newTomogram.copyAttributes(tomo, '_origin')

        output.append(newTomogram)
        output.update(newTomogram)
        output.write()
        self._store()

    def closeOutputSetsStep(self):
        self.Tomograms.setStreamState(Set.STREAM_CLOSED)
        self.Tomograms.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------
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

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if self.Tomograms:
            summary.append("Input Tilt-Series: %d.\nInterpolations applied: %d.\n"
                           % (self.inputSetOfTomograms.get().getSize(),
                              self.Tomograms.getSize()))
        else:
            summary.append("Output not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.Tomograms:
            methods.append("%d tomograms have been normalized using the IMOD newstack program.\n"
                           % (self.Tomograms.getSize()))
        else:
            methods.append("Output not ready yet.")
        return methods
