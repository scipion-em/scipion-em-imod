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
from pwem.emlib.image import ImageHandler
from imod import Plugin


class ProtImodCtfEstimationInteractive(EMProtocol, ProtTomoBase):
    """
    CTF estimation of a set of input tilt-series using the IMOD procedure and ctfplotter GUI.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/ctfplotter.html
    """

    _label = 'CTF estimation interactive'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series',
                      help='This should be a raw stack, not an aligned stack, because the interpolation used to make '
                           'an aligned stack attenuates high frequencies and the noise power spectra would no longer '
                           'match.')

        form.addParam('expectedDefocus',
                      params.FloatParam,
                      default=8000,
                      label='Expected defocus',
                      important=True,
                      help='Expected defocus at the tilt axis in nanometers, with a positive value for underfocus. '
                           'The frequency of the first zero of the CTF curve is first computed based on this expected '
                           'defocus.  The segments of the CTF curve of the input stack around that frequency are '
                           'selected to be fitted.')

        form.addParam('axisAngle',
                      params.FloatParam,
                      default=0.0,
                      label='Axis angle',
                      important=True,
                      help='Specifies how much the tilt axis deviates from vertical (Y axis). This angle is in degrees.'
                           ' It follows the right hand  rule and counter-clockwise is positive.')

        form.addParam('skipAstigmaticViews',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Skip astigmatic phase views',
                      expertLevel=params.LEVEL_ADVANCED,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Skip or break views only when finding astigmatism or phase shift')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('ctfEstimation', ts.getObjId())
            self._insertFunctionStep('createOutputStep', ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)
        outputTsFileName = os.path.join(extraPrefix, '%s_ctfEstimated.st' % tsId)
        outputTltFileName = os.path.join(tmpPrefix, '%s.rawtlt' % tsId)

        """Apply the transformation form the input tilt-series"""
        ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        ts.generateTltFile(outputTltFileName)

    def ctfEstimation(self, tsObjId):
        """Run ctfplotter IMOD program"""
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsCtfPlotter = {
            'inputStack': os.path.join(extraPrefix, '%s_ctfEstimated.st' % tsId),
            'angleFile': os.path.join(tmpPrefix, '%s.rawtlt' % tsId),
            'defocusFile': os.path.join(extraPrefix, '%s.defocus' % tsId),
            'axisAngle': self.axisAngle.get(),
            'pixelSize': self.inputSetOfTiltSeries.get().getSamplingRate() / 10,
            'expectedDefocus': self.expectedDefocus.get(),
            'voltage': self.inputSetOfTiltSeries.get().getAcquisition().getVoltage(),
            'sphericalAberration': self.inputSetOfTiltSeries.get().getAcquisition().getSphericalAberration(),
            'amplitudeContrast': self.inputSetOfTiltSeries.get().getAcquisition().getAmplitudeContrast(),
            'angleRange': '0,0',
        }

        argsCtfPlotter = "-InputStack %(inputStack)s " \
                         "-AngleFile %(angleFile)s " \
                         "-DefocusFile %(defocusFile)s " \
                         "-AxisAngle %(axisAngle)f " \
                         "-PixelSize %(pixelSize)f " \
                         "-ExpectedDefocus %(expectedDefocus)f " \
                         "-Voltage %(voltage)d " \
                         "-SphericalAberration %(sphericalAberration)f " \
                         "-AngleRange %(angleRange)s "

        Plugin.runImod(self, 'ctfplotter', argsCtfPlotter % paramsCtfPlotter)

    def createOutputStep(self, tsObjId):
        outputCtfEstimatedSetOfTiltSeries = self.getOutputCtfEstimatedSetOfTiltSeries()

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputCtfEstimatedSetOfTiltSeries.append(newTs)

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(index + 1, (os.path.join(extraPrefix, '%s_ctfEstimated.st' % tsId)))
            newTs.append(newTi)

        newTs.write()
        outputCtfEstimatedSetOfTiltSeries.update(newTs)
        outputCtfEstimatedSetOfTiltSeries.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputCtfEstimatedSetOfTiltSeries(self):
        if not hasattr(self, "outputCtfEstimatedSetOfTiltSeries"):
            outputCtfEstimatedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='CtfEstimated')
            outputCtfEstimatedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputCtfEstimatedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            self._defineOutputs(outputCtfEstimatedSetOfTiltSeries=outputCtfEstimatedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputCtfEstimatedSetOfTiltSeries)
        return self.outputCtfEstimatedSetOfTiltSeries

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputCtfEstimatedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nnumber of CTF estimated: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputCtfEstimatedSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputCtfEstimatedSetOfTiltSeries'):
            methods.append("%d Tilt-series CTF have been estimated using the IMOD ctfplotter software.\n"
                           % (self.outputCtfEstimatedSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
