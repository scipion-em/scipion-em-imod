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
import numpy as np
import imod.utils as utils
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase


class ProtTSNormalization(EMProtocol, ProtTomoBase):
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
                      label='Input set of tilt-Series')

        form.addParam('binning',
                      params.FloatParam,
                      default=1.0,
                      label='Binning',
                      help='Binning to be applied to the interpolated tilt-series. '
                           'Must be a integer bigger than 1')

        form.addParam('floatDensities',
                      params.EnumParam,
                      choices=['1', '2', '3', '4'],
                      default=0,
                      label='Adjust densities mode',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Adjust densities of sections individually:\n'
                           '-Mode 1: sections fill the data range. Default option\n'
                           '-Mode 2: sections scaled to common mean and standard deviation\n'
                           '-Mode 3: sections shifted to a common mean without scaling\n'
                           '-Mode 4: sections shifted to a common mean and then rescale\n'
                           'the resulting minimum and maximum densities to the Min and Max values specified')

        form.addParam('modeToOutput',
                      params.EnumParam,
                      choices=['default', '4-bit', 'byte', 'signed 16-bit', 'unsigned 16-bit', '32-bit float'],
                      default=0,
                      label='Storage data type',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='The storage mode of the output file. The default is the mode of the first input file, '
                           'except for a 4-bit input file, where the default is to output as bytes')

        form.addParam('scalingToggle',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Rescale densities',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='The storage mode of the output file. The default is the mode of the first input file, '
                           'except for a 4-bit input file, where the default is to output as bytes')

        group = form.addGroup('Scaling values',
                              condition='(scalingToggle==0) | (floatDensities==3)')

        group.addParam('scaleMax', params.FloatParam,
                       default=255,
                       label='Max.',
                       help='Maximum value for the rescaling')

        group.addParam('scaleMin', params.FloatParam,
                       default=0,
                       label='Min.',
                       help='Minimum value for the rescaling')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('generateTransformFileStep', ts.getObjId())
            self._insertFunctionStep('generateOutputStackStep', ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)
        outputTsFileName = os.path.join(tmpPrefix, "%s.st" % tsId)

        """Apply the transformation form the input tilt-series"""
        ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)
        ts.generateTltFile(angleFilePath)

    def generateOutputStackStep(self, tsObjId):
        outputInterpolatedSetOfTiltSeries = self.getOutputInterpolatedSetOfTiltSeries()

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputInterpolatedSetOfTiltSeries.append(newTs)
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        if ts.getFirstItem().hasTransform():
            paramsAlignment = {
                'input': os.path.join(tmpPrefix, '%s.st' % tsId),
                'output': os.path.join(extraPrefix, '%s.st' % tsId),
                'bin': int(self.binning.get()),
                'imagebinned': 1.0}

            argsAlignment = "-input %(input)s " \
                            "-output %(output)s " \
                            "-bin %(bin)d " \
                            "-imagebinned %(imagebinned)s"
            self.runJob('newstack', argsAlignment % paramsAlignment)

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(index + 1, (os.path.join(extraPrefix, '%s.st' % tsId)))
            if self.binning > 1:
                newTi.setSamplingRate(tiltImage.getSamplingRate() * int(self.binning.get()))
            newTs.append(newTi)

        if self.binning > 1:
            newTs.setSamplingRate(ts.getSamplingRate() * int(self.binning.get()))

        newTs.write()
        outputInterpolatedSetOfTiltSeries.update(newTs)
        outputInterpolatedSetOfTiltSeries.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputNormalizedSetOfTiltSeries(self):
        if not hasattr(self, "outputNormalizedSetOfTiltSeries"):
            outputNormalizedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Normalized')
            outputNormalizedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputNormalizedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            if self.binning > 1:
                samplingRate = self.inputSetOfTiltSeries.get().getSamplingRate()
                samplingRate *= self.binning.get()
                outputNormalizedSetOfTiltSeries.setSamplingRate(samplingRate)
            self._defineOutputs(outputNormalizedSetOfTiltSeries=outputNormalizedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputNormalizedSetOfTiltSeries)
        return self.outputNormalizedSetOfTiltSeries

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
            methods.append("The interpolation has been computed for %d "
                           "Tilt-series using the IMOD newstack program.\n"
                           % (self.outputNormalizedSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
