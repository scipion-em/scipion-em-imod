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
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
from imod import Plugin
from imod import utils


class ProtImodCtfCorrection(EMProtocol, ProtTomoBase):
    """
    CTF correction of a set of input tilt-series using the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/ctfphaseflip.html
    """

    _label = 'CTF correction'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      label="Input tilt-series",
                      pointerClass='SetOfTiltSeries',
                      help='Select the set of tilt-series to be CTF corrected.')

        form.addParam('inputSetOfCtfTomoSeries',
                      params.PointerParam,
                      label="input tilt-series CTF estimation",
                      pointerClass='SetOfCTFTomoSeries',
                      help='Select the CTF estimation from the set of tilt-series.')

        form.addParam('defocusTol',
                      params.FloatParam,
                      label='Defocus tolerance (nm)',
                      default=200,
                      important=True,
                      help='The value introduced must be the same used for CTF estimation in case ti has been '
                           'performed with IMOD. \n\n'
                           'Defocus tolerance in nanometers defining the center strips. The center strips are taken '
                           'from the central region of a view that has defocus difference less than this tolerance. '
                           'These kind of center strips from all views within AngleRange are considered to have a '
                           'constant defocus and are used to compute the initial CTF after being further tessellated '
                           'into tiles.')

        form.addParam('interpolationWidth',
                      params.IntParam,
                      label='Interpolation Width (pixels)',
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

        form.addHidden(params.USE_GPU,
                       params.BooleanParam,
                       default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation.\
                       Select the one you want to use.")

        form.addHidden(params.GPU_LIST,
                       params.StringParam,
                       default='0',
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="GPU ID. To pick the best available one set 0. For a specific GPU set its number ID.")

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('ctfCorrection', ts.getObjId())
            self._insertFunctionStep('createOutputStep', ts.getObjId())
        self._insertFunctionStep('closeOutputSetsStep')

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        ctfTomoSeries = self.inputSetOfCtfTomoSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        """Apply the transformation form the input tilt-series"""
        outputTsFileName = os.path.join(tmpPrefix, ts.getFirstItem().parseFileName())
        ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, ts.getFirstItem().parseFileName(extension=".tlt"))
        ts.generateTltFile(angleFilePath)

        """Generate defocus file"""
        defocusFilePath = os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".defocus"))
        utils.generateDefocusIMODFileFromObject(ctfTomoSeries, defocusFilePath)

    def ctfCorrection(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        """Run ctfphaseflip IMOD program"""
        paramsCtfPhaseFlip = {
            'inputStack': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName()),
            'angleFile': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName(extension=".tlt")),
            'outputFileName': os.path.join(extraPrefix, ts.getFirstItem().parseFileName()),
            'defocusFile': os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".defocus")),
            'voltage': self.inputSetOfTiltSeries.get().getAcquisition().getVoltage(),
            'sphericalAberration': self.inputSetOfTiltSeries.get().getAcquisition().getSphericalAberration(),
            'defocusTol': self.defocusTol.get(),
            'pixelSize': self.inputSetOfTiltSeries.get().getSamplingRate()/10,
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

        if self.usesGpu():
            paramsCtfPhaseFlip.update({
                "useGPU": self.getGpuList()[0]
            })

            argsCtfPhaseFlip += "-UseGPU %(useGPU)d "

        Plugin.runImod(self, 'ctfphaseflip', argsCtfPhaseFlip % paramsCtfPhaseFlip)

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
            newTi.setLocation(index + 1,
                              (os.path.join(extraPrefix, tiltImage.parseFileName())))
            newTs.append(newTi)

        newTs.write(properties=False)
        outputCtfCorrectedSetOfTiltSeries.update(newTs)
        outputCtfCorrectedSetOfTiltSeries.write()
        self._store()

    def closeOutputSetsStep(self):
        self.getOutputCtfCorrectedSetOfTiltSeries().setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputCtfCorrectedSetOfTiltSeries(self):
        if hasattr(self, "outputCtfCorrectedSetOfTiltSeries"):
            self.outputCtfCorrectedSetOfTiltSeries.enableAppend()
        else:
            outputCtfCorrectedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='CtfCorrected')
            outputCtfCorrectedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputCtfCorrectedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            outputCtfCorrectedSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputCtfCorrectedSetOfTiltSeries=outputCtfCorrectedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries.get(), outputCtfCorrectedSetOfTiltSeries)
        return self.outputCtfCorrectedSetOfTiltSeries

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        if self.inputSetOfTiltSeries.get().getSize() != self.inputSetOfCtfTomoSeries.get().getSize():
            validateMsgs.append("Input set of tilt-series and input set of CTF tomo estimations must contain the "
                                "same number of elements.")

        return validateMsgs

    def _summary(self):
        summary = []
        if hasattr(self, 'outputCtfCorrectedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nCTF corrections applied: %d.\n"
                           % (self.inputSetOfCtfTomoSeries.get().getSize(),
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
