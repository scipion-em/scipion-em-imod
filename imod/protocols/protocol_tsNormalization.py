# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
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
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
from imod import Plugin
from pwem.emlib.image import ImageHandler


class ProtImodTSNormalization(EMProtocol, ProtTomoBase):
    """
    Normalize input tilt-series and change its storing formatting.
    More info:
        https://bio3D.colorado.edu/imod/doc/etomoTutorial.html
    """

    _label = 'tilt-series normalization'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        form.addParam('binning',
                      params.FloatParam,
                      default=1.0,
                      label='Binning',
                      important=True,
                      help='Binning to be applied to the normalized tilt-series in IMOD convention. Images will be '
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
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('generateOutputStackStep', ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        """Apply the transformation form the input tilt-series"""
        outputTsFileName = os.path.join(tmpPrefix, ts.getFirstItem().parseFileName())
        ts.applyTransform(outputTsFileName)

    def generateOutputStackStep(self, tsObjId):
        outputNormalizedSetOfTiltSeries = self.getOutputNormalizedSetOfTiltSeries()

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsNewstack = {
            'input': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName()),
            'output': os.path.join(extraPrefix, ts.getFirstItem().parseFileName(suffix="_norm")),
            'bin': int(self.binning.get()),
            'imagebinned': 1.0,
        }

        argsNewstack = "-input %(input)s " \
                       "-output %(output)s " \
                       "-bin %(bin)d " \
                       "-imagebinned %(imagebinned)s "

        if self.floatDensities.get() != 0:
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
            argsNewstack += " -ModeToOutput " + str(self.getModeToOutput())

        Plugin.runImod(self, 'newstack', argsNewstack % paramsNewstack)

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputNormalizedSetOfTiltSeries.append(newTs)

        if self.binning > 1:
            newTs.setSamplingRate(ts.getSamplingRate() * int(self.binning.get()))

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(index + 1, (os.path.join(extraPrefix, tiltImage.parseFileName(suffix="_norm"))))
            if self.binning > 1:
                newTi.setSamplingRate(tiltImage.getSamplingRate() * int(self.binning.get()))
            newTs.append(newTi)

        ih = ImageHandler()
        x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
        newTs.setDim((x, y, z))

        newTs.write(properties=False)
        outputNormalizedSetOfTiltSeries.update(newTs)
        outputNormalizedSetOfTiltSeries.updateDim()
        outputNormalizedSetOfTiltSeries.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputNormalizedSetOfTiltSeries(self):
        if hasattr(self, "outputNormalizedSetOfTiltSeries"):
            self.outputNormalizedSetOfTiltSeries.enableAppend()
        else:
            outputNormalizedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Normalized')
            outputNormalizedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputNormalizedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            if self.binning > 1:
                samplingRate = self.inputSetOfTiltSeries.get().getSamplingRate()
                samplingRate *= self.binning.get()
                outputNormalizedSetOfTiltSeries.setSamplingRate(samplingRate)
            outputNormalizedSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputNormalizedSetOfTiltSeries=outputNormalizedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputNormalizedSetOfTiltSeries)
        return self.outputNormalizedSetOfTiltSeries

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
        if hasattr(self, 'outputNormalizedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nInterpolations applied: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputNormalizedSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputNormalizedSetOfTiltSeries'):
            methods.append("%d tilt-series have been normalized using the IMOD newstack program.\n"
                           % (self.outputNormalizedSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
