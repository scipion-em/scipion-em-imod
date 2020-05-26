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


class ProtCtfEstimation(EMProtocol, ProtTomoBase):
    """
    CTF estimation of a set of input tilt-series using the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/ctfplotter.html
    """

    _label = 'CTF estimation'

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

        form.addParam('defocusTol',
                      params.FloatParam,
                      label='Defocus tolerance',
                      default=200,
                      important=True,
                      help='Defocus tolerance in nanometers defining the center strips. The center strips are taken '
                           'from the central region of a view that has defocus difference less than this tolerance. '
                           'These kind of center strips from all views within AngleRange are considered to have a '
                           'constant defocus and are used to compute the initial CTF after being further tessellated '
                           'into tiles.')

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

        form.addParam('leftDefTol',
                      params.FloatParam,
                      label='Left defocus tolerance',
                      default=2000,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Defocus tolerance in nanometers for strips to the left of the center strip.")

        form.addParam('rightDefTol',
                      params.FloatParam,
                      label='Right defocus tolerance',
                      default=2000,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Defocus tolerance in nanometers for strips to the right of the center strip.")

        form.addParam('tileSize',
                      params.IntParam,
                      label='Tile size',
                      default=256,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="The tile size each strip will be tessellated into. The size is in pixels and the tiles are "
                           "square.  Each view is first divided into strips that are considered to have constant "
                           "defocus.")

        groupAngleRange = form.addGroup('Autorefinement angle settings',
                                        help='This entry sets the range of angles in each fit and the step size between'
                                             'angles for the initial autofitting of the whole tilt-series.')

        groupAngleRange.addParam('angleStep',
                                 params.FloatParam,
                                 default=2,
                                 label='Angle step',
                                 help='Step size between ranges. A value of zero for the step will make it fit to each '
                                      'single image separately, regardless of the value for the range.')

        groupAngleRange.addParam('angleRange',
                                 params.FloatParam,
                                 condition="angleStep != 0",
                                 default=120,
                                 label='Angle range',
                                 help='Size of the angle range in which the CTF is estimated.')

        groupFrequencyRange = form.addGroup('Autorefinement frequency range',
                                            expertLevel=params.LEVEL_ADVANCED,
                                            help='Starting and ending frequencies of range to fit in power spectrum. '
                                                 'The two values will be used to set the "X1 starts" and "X2 ends" '
                                                 'fields in the fitting dialog.')

        groupFrequencyRange.addParam('startFreq',
                                     params.FloatParam,
                                     default=0.0,
                                     label='Start',
                                     help='Starting frequency. "X1 starts". "X2 ends".')

        groupFrequencyRange.addParam('endFreq',
                                     params.FloatParam,
                                     default=0.0,
                                     label='End',
                                     help='Ending frequency. ')

        form.addParam('extraZerosToFit',
                      params.FloatParam,
                      label='Extra zeros to fit',
                      default=0.0,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="By default, the ending frequency of the fitting range is set to the expected location of "
                           "the second zero.  With this entry, the range will be extended by the given multiple of the "
                           "interval between first and seconds zeros.  For example, entries of 1 and 2 will fit "
                           "approximately to the third and fourth zeros, respectively.  An entry of more than 0.5 will "
                           "trigger fitting to two exponentials, which is important for fitting multiple peaks between "
                           "zeros.")

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
            'autoFitRangeAndStep': str(self.angleRange.get()) + "," + str(self.angleStep.get()),
            'voltage': self.inputSetOfTiltSeries.get().getAcquisition().getVoltage(),
            'sphericalAberration': self.inputSetOfTiltSeries.get().getAcquisition().getSphericalAberration(),
            'amplitudeContrast': self.inputSetOfTiltSeries.get().getAcquisition().getAmplitudeContrast(),
            'defocusTol': self.defocusTol.get(),
            'psResolution': 101,
            'leftDefTol': self.leftDefTol.get(),
            'rightDefTol': self.rightDefTol.get(),
            'tileSize': self.tileSize.get(),
        }

        argsCtfPlotter = "-InputStack %(inputStack)s " \
                         "-AngleFile %(angleFile)s " \
                         "-DefocusFile %(defocusFile)s " \
                         "-AxisAngle %(axisAngle)f " \
                         "-PixelSize %(pixelSize)f " \
                         "-ExpectedDefocus %(expectedDefocus)f " \
                         "-AutoFitRangeAndStep %(autoFitRangeAndStep)s " \
                         "-Voltage %(voltage)d " \
                         "-SphericalAberration %(sphericalAberration)f " \
                         "-AmplitudeContrast %(amplitudeContrast)f " \
                         "-DefocusTol %(defocusTol)d " \
                         "-PSResolution %(psResolution)d " \
                         "-LeftDefTol %(leftDefTol)f " \
                         "-RightDefTol %(rightDefTol)f " \
                         "-tileSize %(tileSize)d " \
                         "-SaveAndExit "

        if self.startFreq.get() != 0 or self.endFreq.get() != 0:
            paramsCtfPlotter.update({
                'startFreq': self.startFreq.get(),
                'endFreq': self.endFreq.get()
            })

            argsCtfPlotter += "-FrequencyRangeToFit %(startFreq)f,%(endFreq)f "

        if self.extraZerosToFit.get() != 0:
            paramsCtfPlotter.update({
                'extraZerosToFit': self.extraZerosToFit.get(),
            })

            argsCtfPlotter += "-ExtraZerosToFit %(extraZerosToFit)f "

        Plugin.runImod(self, 'ctfplotter', argsCtfPlotter % paramsCtfPlotter)

    def createOutputStep(self, tsObjId):
        outputCtfEstimatedSetOfTiltSeries = self.getOutputCtfEstimatedSetOfTiltSeries()

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        ctfInfoList = self.readCtfOutputFile(os.path.join(extraPrefix, "%s.defocus" % tsId))

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
    def readCtfOutputFile(self, ctfOutputFile):
        cttInfoList = []

        return cttInfoList

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
