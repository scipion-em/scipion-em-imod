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


class ProtImodCtfCorrection(EMProtocol, ProtTomoBase):
    """
    CTF correction of a set of input tilt-series using the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/ctfphaseflip.html
    """

    _label = 'CTF correction'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('protCtfEstimation',
                      params.PointerParam,
                      label="IMOD CTF estimation run",
                      pointerClass='ProtImodCtfEstimation',
                      help='Select the previous IMOD CTF estimation run.')

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
                           "tilt. However, strips at constant spacing can still be obtained by entering the negative "
                           "of the desired spacing, which disables the scaling of the spacing.")

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self.inputSetOfTiltSeries = self.protCtfEstimation.get().outputCtfEstimatedSetOfTiltSeries
        for ts in self.inputSetOfTiltSeries:
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('ctfCorrection', ts.getObjId())
            self._insertFunctionStep('createOutputStep', ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        """Apply the transformation form the input tilt-series"""
        outputTsFileName = os.path.join(tmpPrefix, ts.getFirstItem().parseFileName())
        self.protCtfEstimation.get().inputSetOfTiltSeries.get()[tsObjId].applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, ts.getFirstItem().parseFileName(extension=".tlt"))
        self.protCtfEstimation.get().inputSetOfTiltSeries.get()[tsObjId].generateTltFile(angleFilePath)

    def ctfCorrection(self, tsObjId):
        """Run ctfphaseflip IMOD program"""

        ts = self.inputSetOfTiltSeries[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsCtfPhaseFlip = {
            'inputStack': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName()),
            'angleFile': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName(extension=".tlt")),
            'outputFileName': os.path.join(extraPrefix, ts.getFirstItem().parseFileName()),
            'defocusFile': os.path.join(
                self.protCtfEstimation.get()._getExtraPath(tsId),
                self.protCtfEstimation.get().inputSetOfTiltSeries.get()[tsObjId].getFirstItem().parseFileName(extension=".defocus")),
            'voltage': self.inputSetOfTiltSeries.getAcquisition().getVoltage(),
            'sphericalAberration': self.inputSetOfTiltSeries.getAcquisition().getSphericalAberration(),
            'defocusTol': self.protCtfEstimation.get().defocusTol.get(),
            'pixelSize': self.inputSetOfTiltSeries.getSamplingRate()/10,
            'amplitudeContrast': self.inputSetOfTiltSeries.getAcquisition().getAmplitudeContrast(),
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

        Plugin.runImod(self, 'ctfphaseflip', argsCtfPhaseFlip % paramsCtfPhaseFlip)

    def createOutputStep(self, tsObjId):
        outputCtfCorrectedSetOfTiltSeries = self.getOutputCtfCorrectedSetOfTiltSeries()

        ts = self.inputSetOfTiltSeries[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputCtfCorrectedSetOfTiltSeries.append(newTs)

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(index + 1,
                              (os.path.join(extraPrefix, tiltImage.parseFileName())))
            newTs.append(newTi)

        newTs.write(properties=False)
        outputCtfCorrectedSetOfTiltSeries.update(newTs)
        outputCtfCorrectedSetOfTiltSeries.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputCtfCorrectedSetOfTiltSeries(self):
        if hasattr(self, "outputCtfEstimatedSetOfTiltSeries"):
            self.outputCtfEstimatedSetOfTiltSeries.enableAppend()
        else:
            outputCtfCorrectedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='CtfCorrected')
            outputCtfCorrectedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries)
            outputCtfCorrectedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.getDim())
            outputCtfCorrectedSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputCtfCorrectedSetOfTiltSeries=outputCtfCorrectedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputCtfCorrectedSetOfTiltSeries)
        return self.outputCtfCorrectedSetOfTiltSeries

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        if not hasattr(self.protCtfEstimation.get(), 'outputCtfEstimatedSetOfTiltSeries'):
            validateMsgs = "You need to generate an estimation of the CTF to calculate its correction"

        return validateMsgs

    def _summary(self):
        summary = []
        if hasattr(self, 'outputCtfCorrectedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nCTF corrections applied: %d.\n"
                           % (self.protCtfEstimation.get().outputCtfEstimatedSetOfTiltSeries.getSize(),
                              self.outputCtfCorrectedSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputCtfCorrectedSetOfTiltSeries'):
            methods.append("%d Tilt-series have been CTF corrected using the IMOD ctfphaseflip software.\n"
                           % (self.outputCtfCorrectedSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
