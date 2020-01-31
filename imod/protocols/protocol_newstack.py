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


class ProtNewstack(EMProtocol, ProtTomoBase):
    """
    Compute the interpolated tilt-series from its transform matrix.
    More info:
        https://bio3D.colorado.edu/imod/doc/etomoTutorial.html
    """

    _label = 'newstack'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-Series')

        form.addParam('binning', params.FloatParam,
                      default=1.0,
                      label='Binning',
                      help='Binning to be applied to the interpolated tilt-series. '
                           'Must be a integer bigger than 1')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('generateTransformFileStep')
        self._insertFunctionStep('generateOutputStackStep')

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            extraPrefix = self._getExtraPath(tsId)
            tmpPrefix = self._getTmpPath(tsId)
            path.makePath(tmpPrefix)
            path.makePath(extraPrefix)
            inputTsFileName = ts.getFirstItem().getLocation()[1]
            outputTsFileName = os.path.join(tmpPrefix, "%s.st" % tsId)

            """Create link to input stack"""
            path.createLink(inputTsFileName, outputTsFileName)

    def generateTransformFileStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            extraPrefix = self._getExtraPath(tsId)
            if not ts.getFirstItem().hasTransform():
                inputTsFilePath = ts.getFirstItem().getLocation()[1]
                outputTsFilePath = os.path.join(extraPrefix, '%s.st' % tsId)
                path.createLink(inputTsFilePath, outputTsFilePath)
            else:
                utils.formatTransformFile(ts, os.path.join(extraPrefix, "%s.prexg" % tsId))


    def generateOutputStackStep(self):
        outputInterpolatedSetOfTiltSeries = self.getOutputInterpolatedSetOfTiltSeries()

        for ts in self.inputSetOfTiltSeries.get():
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
                    'xform': os.path.join(extraPrefix, "%s.prexg" % tsId),
                    'bin': int(self.binning.get()),
                    'mode': 0,
                    'float': 2,
                    'imagebinned': 1.0}

                argsAlignment = "-input %(input)s " \
                                "-output %(output)s " \
                                "-xform %(xform)s " \
                                "-bin %(bin)d " \
                                "-mode %(mode)s " \
                                "-float %(float)s " \
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
    def getOutputInterpolatedSetOfTiltSeries(self):
        if not hasattr(self, "outputInterpolatedSetOfTiltSeries"):
            outputInterpolatedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Interpolated')
            outputInterpolatedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputInterpolatedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            if self.binning > 1:
                samplingRate = self.inputSetOfTiltSeries.get().getSamplingRate()
                samplingRate *= self.binning.get()
                outputInterpolatedSetOfTiltSeries.setSamplingRate(samplingRate)
            self._defineOutputs(outputInterpolatedSetOfTiltSeries=outputInterpolatedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputInterpolatedSetOfTiltSeries)
        return self.outputInterpolatedSetOfTiltSeries

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nInterpolations applied: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            methods.append("The interpolation has been computed for %d "
                           "Tilt-series using the IMOD newstack program.\n"
                           % (self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods