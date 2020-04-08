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


class ProtCtfEstimation(EMProtocol, ProtTomoBase):
    """
    CTF estimation of a set of input tilt-series using the IMOD procedure.
    More info:
        https://bio3D.colorado.edu/imod/doc/etomoTutorial.html
    """

    _label = 'CTF correction'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-Series')

        form.addParam('defocusTol',
                      params.IntParam,
                      label='Defocus tolerance',
                      default='200',
                      important=True,
                      help="Defocus tolerance in nanometers defining the center strips. The center strips are taken "
                           "from the central region of a view that has defocus difference less than this tolerance. "
                           "These kind of center strips from all views within AngleRange are considered to have a "
                           "constant defocus and are used to compute the initial CTF after being further tessellated "
                           "into tiles.")

        form.addParam('interpolationWidth',
                      params.IntParam,
                      label='Interpolation Width',
                      default='15',
                      important=True,
                      help="The distance in pixels between the center lines of two consecutive strips. A pixel inside "
                           "the region between those two center lines resides in both strips. As the two strips are "
                           "corrected separately, that pixel will have 2 corrected values. The final value for that "
                           "pixel is a linear interpolation of the 2 corrected values. If a value of 1 is entered, "
                           "there is no such interpolation.  For a value greater than one, the entered value will be "
                           "used whenever the strip width is less than 256 (i.e., at high tilt), and the value will be "
                           "scaled proportional to the strip width for widths above 256.  This scaling keeps the "
                           "computational time down and is reasonable because the defocus difference between adjacent "
                           "wide strips at wider intervals is still less than that between the narrower strips at high "
                           "tilt.  However, strips at constant spacing can still be obtained by entering the negative "
                           "of the desired spacing, which disables the scaling of the spacing.")

        form.addParam('expectedDefocus',
                      params.FloatParam,
                      default=8000.0,
                      label='Expected defocus',
                      important=True,
                      help='Expected defocus at the tilt axis in nanometers, with a positive value for underfocus.  '
                           'The frequency of the first zero of the CTF curve is first computed based on this expected '
                           'defocus.  The segments of the CTF curve of the input stack around that frequency are '
                           'selected to be fitted.')

        form.addParam('tileSize',
                      params.IntParam,
                      default=256,
                      label='Tile size',
                      important=True,
                      help='The tile size each strip will be tessellated into. The size is in pixels and the tiles are '
                           'square. Each view is first divided into strips that are considered to have constant '
                           'defocus.')

        groupAngleRange = form.addGroup('Angle range',
                                        help='This entry sets the starting and ending tilt angles for the initial '
                                             'analysis and is a required entry.  Views with a tilt angle within this '
                                             'range are used to compute the CTF curve.')

        groupAngleRange.addParam('minAngleRange',
                                 params.FloatParam,
                                 default=-20.0,
                                 label='Minimum angle range',
                                 help='Lower limit of the angle range.')

        groupAngleRange.addParam('maxAngleRange',
                                 params.FloatParam,
                                 default=20.0,
                                 label='Maximum angle range',
                                 help='Upper limit of the angle range.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('ctfEstimation', ts.getObjId())
            self._insertFunctionStep('ctfCorrection', ts.getObjId())
            self._insertFunctionStep('createOutputStep', ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)
        outputTsFileName = os.path.join(tmpPrefix, "%s.st" % tsId)
        outputTltFileName = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)

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
            'inputStack': os.path.join(tmpPrefix, '%s.st' % tsId),
            'angleFile': os.path.join(tmpPrefix, '%s.rawtlt' % tsId),
            'defocusFile': os.path.join(extraPrefix, '%s.defocus' % tsId),
            'axisAngle': 0.0,
            'pixelSize': self.inputSetOfTiltSeries.get().getSamplingRate(),
            'expectedDefocus': self.expectedDefocus.get(),
            'angleRange': str(self.minAngleRange.get()) + "," + str(self.maxAngleRange.get()),
            'voltage': self.inputSetOfTiltSeries.get().getAcquisition().getVoltage(),
            'sphericalAberration': self.inputSetOfTiltSeries.get().getAcquisition().getSphericalAberration(),
            'amplitudeContrast': self.inputSetOfTiltSeries.get().getAcquisition().getAmplitudeContrast(),
            'defocusTol': self.defocusTol.get(),
            'psResolution': 101,
            'tileSize': self.tileSize.get(),
            'leftDefTol': 2000.0,
            'rightDefTol': 2000.0,
        }

        argsCtfPlotter = "-InputStack %(inputStack)s " \
                         "-AngleFile %(angleFile)s " \
                         "-DefocusFile %(defocusFile)s " \
                         "-AxisAngle %(axisAngle)f " \
                         "-PixelSize %(pixelSize)f " \
                         "-ExpectedDefocus %(expectedDefocus)f " \
                         "-AutoFitRangeAndStep %(angleRange)s " \
                         "-Voltage %(voltage)d " \
                         "-SphericalAberration %(sphericalAberration)f " \
                         "-AmplitudeContrast %(amplitudeContrast)f " \
                         "-DefocusTol %(defocusTol)d " \
                         "-PSResolution %(psResolution)d " \
                         "-TileSize %(tileSize)d " \
                         "-LeftDefTol %(leftDefTol)f " \
                         "-RightDefTol %(rightDefTol)f "

        self.runJob('ctfplotter', argsCtfPlotter % paramsCtfPlotter + " -config /home/fede/Downloads/ctf-sirt/F20.cfg ")

    def ctfCorrection(self, tsObjId):
        """Run ctfphaseflip IMOD program"""

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsCtfPhaseFlip = {
            'inputStack': os.path.join(tmpPrefix, '%s.st' % tsId),
            'angleFile': os.path.join(tmpPrefix, '%s.rawtlt' % tsId),
            'outputFileName': os.path.join(extraPrefix, '%s_ctfCorrected.st' % tsId),
            'defocusFile': os.path.join(extraPrefix, '%s.defocus' % tsId),
            'voltage': self.inputSetOfTiltSeries.get().getAcquisition().getVoltage(),
            'sphericalAberration': self.inputSetOfTiltSeries.get().getAcquisition().getSphericalAberration(),
            'defocusTol': self.defocusTol.get(),
            'pixelSize': self.inputSetOfTiltSeries.get().getSamplingRate(),
            'amplitudeContrast': self.inputSetOfTiltSeries.get().getAcquisition().getAmplitudeContrast(),
            'interpolationWidth': self.interpolationWidth.get(),
        }

        argsCtfPhaseFlip = "-InputStack %(inputStack)s " \
                           "-AngleFile %(angleFile)s " \
                           "-OutputFileName %(outputFileName)s " \
                           "-DefocusFile %(defocusFile)s " \
                           "-Voltage %(voltage)d " \
                           "-SphericalAberration %(sphericalAberration)f " \
                           "-DefocusTol %(defocusTol)d " \
                           "-PixelSize %(pixelSize)f " \
                           "-AmplitudeContrast %(amplitudeContrast)f " \
                           "-InterpolationWidth %(interpolationWidth)d "

        self.runJob('ctfphaseflip', argsCtfPhaseFlip % paramsCtfPhaseFlip)

    def createOutputStep(self, tsObjId):
        outputCtfCorrectedSetOfTiltSeries = self.getOutputCtfCorrectedSetOfTiltSeries()

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputCtfCorrectedSetOfTiltSeries.append(newTs)

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(index + 1, (os.path.join(extraPrefix, '%s_ctfCorrected.st' % tsId)))
            newTs.append(newTi)

        newTs.write()
        outputCtfCorrectedSetOfTiltSeries.update(newTs)
        outputCtfCorrectedSetOfTiltSeries.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputCtfCorrectedSetOfTiltSeries(self):
        if not hasattr(self, "outputCtfCorrectedSetOfTiltSeries"):
            outputCtfCorrectedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='CtfCorrected')
            outputCtfCorrectedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputCtfCorrectedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            self._defineOutputs(outputCtfCorrectedSetOfTiltSeries=outputCtfCorrectedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputCtfCorrectedSetOfTiltSeries)
        return self.outputCtfCorrectedSetOfTiltSeries

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputCtfCorrectedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nCTF corrections applpied applied: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputCtfCorrectedSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputCtfCorrectedSetOfTiltSeries'):
            methods.append("%d Tilt-series have been CTF corrected using the IMOD newstack program.\n"
                           % (self.outputCtfCorrectedSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
