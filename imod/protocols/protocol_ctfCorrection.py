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
import tomo.objects as tomoObj

from .. import Plugin, utils
from .protocol_base import ProtImodBase


class ProtImodCtfCorrection(ProtImodBase):
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
                      help='Select the set of tilt-series to be '
                           'CTF corrected. Usually this will be the '
                           'tilt-series with alignment information.')

        form.addParam('inputSetOfCtfTomoSeries',
                      params.PointerParam,
                      label="Input CTF estimation",
                      pointerClass='SetOfCTFTomoSeries',
                      help='Select the CTF estimation for the set '
                           'of tilt-series.')

        form.addParam('defocusTol',
                      params.FloatParam,
                      label='Defocus tolerance (nm)',
                      default=200,
                      important=True,
                      help='The value introduced must be the same used for '
                           'CTF estimation with IMOD.\n\n'
                           'Defocus tolerance in nanometers defining the '
                           'center strips. The center strips are taken '
                           'from the central region of a view that has defocus '
                           'difference less than this tolerance. '
                           'These kind of center strips from all views within '
                           'AngleRange are considered to have a '
                           'constant defocus and are used to compute the '
                           'initial CTF after being further tessellated '
                           'into tiles.')

        form.addParam('interpolationWidth',
                      params.IntParam,
                      label='Interpolation width (px)',
                      default='15',
                      important=True,
                      help="The distance in pixels between the center lines "
                           "of two consecutive strips. A pixel inside the "
                           "region between those two center lines resides in "
                           "both strips. As the two strips are corrected "
                           "separately, that pixel will have 2 corrected "
                           "values. The final value for that pixel is a "
                           "linear interpolation of the 2 corrected "
                           "values. If a value of 1 is entered, there is "
                           "no such interpolation. For a value greater "
                           "than one, the entered value will be used "
                           "whenever the strip width is less than 256 "
                           "(i.e., at high tilt), and the value will be "
                           "scaled proportional to the strip width for widths "
                           "above 256.  This scaling keeps the computational "
                           "time down and is reasonable because the defocus "
                           "difference between adjacent wide strips at "
                           "wider intervals is still less than that between "
                           "the narrower strips at high tilt. However, strips "
                           "at constant spacing can still be obtained by "
                           "entering the negative of the desired spacing, "
                           "which disables the scaling of the spacing.")

        form.addHidden(params.USE_GPU,
                       params.BooleanParam,
                       default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation."
                            "Select the one you want to use.")

        form.addHidden(params.GPU_LIST,
                       params.StringParam,
                       default='0',
                       label="Choose GPU IDs",
                       help="GPU ID. To pick the best available one set 0. "
                            "For a specific GPU set its number ID "
                            "(starting from 1).")

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep(self.convertInputStep, ts.getObjId())
            self._insertFunctionStep(self.generateDefocusFile, ts.getObjId())
            self._insertFunctionStep(self.ctfCorrection, ts.getObjId())
            self._insertFunctionStep(self.createOutputStep, ts.getObjId())
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self, tsObjId):
        # Considering swapXY is required to make tilt axis vertical
        super().convertInputStep(tsObjId, doSwap=True)

    def generateDefocusFile(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        self.debug(f"Generating defocus file for {tsObjId} (ObjId), {tsId} (TsId)")

        # Compose the defocus file path
        defocusFilePath = self.getDefocusFileName(ts)

        """Generate defocus file"""
        ctfTomoSeries = self.getCtfTomoSeriesFromTsId(self.inputSetOfCtfTomoSeries.get(), tsId)
        utils.generateDefocusIMODFileFromObject(ctfTomoSeries, defocusFilePath)

    def getDefocusFileName(self, ts):
        """ Returns the path of the defocus filename based on
         the tilt series and creates the folder/s"""

        tmpPrefix = self._getTmpPath(ts.getTsId())
        path.makePath(tmpPrefix)
        defocusFn = ts.getFirstItem().parseFileName(extension=".defocus")
        defocusFilePath = os.path.join(tmpPrefix, defocusFn)
        return defocusFilePath

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
            'defocusFile': self.getDefocusFileName(ts),
            'voltage': self.inputSetOfTiltSeries.get().getAcquisition().getVoltage(),
            'sphericalAberration': self.inputSetOfTiltSeries.get().getAcquisition().getSphericalAberration(),
            'defocusTol': self.defocusTol.get(),
            'pixelSize': self.inputSetOfTiltSeries.get().getSamplingRate() / 10,
            'amplitudeContrast': self.inputSetOfTiltSeries.get().getAcquisition().getAmplitudeContrast(),
            'interpolationWidth': self.interpolationWidth.get()
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
            argsCtfPhaseFlip += f"-UseGPU {self.getGpuList()[0]} " \
                                "-ActionIfGPUFails 2,2 "

        if ts.getFirstItem().hasTransform():
            paramsCtfPhaseFlip['xformFile'] = os.path.join(tmpPrefix, ts.getFirstItem().parseFileName(extension=".xf"))
            argsCtfPhaseFlip += "-TransformFile %(xformFile)s "

        Plugin.runImod(self, 'ctfphaseflip', argsCtfPhaseFlip % paramsCtfPhaseFlip)

    def createOutputStep(self, tsObjId):
        inputTs = self.inputSetOfTiltSeries.get()
        output = self.getOutputSetOfTiltSeries(inputTs)
        hasAlign = inputTs.getFirstItem().getFirstItem().hasTransform()

        ts = inputTs[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        newTs.setCtfCorrected(True)
        output.append(newTs)

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True, copyTM=False)
            acq = tiltImage.getAcquisition()
            if hasAlign:
                acq.setTiltAxisAngle(0.)
            newTi.setAcquisition(acq)
            newTi.setLocation(index + 1,
                              (os.path.join(extraPrefix,
                                            tiltImage.parseFileName())))
            newTs.append(newTi)

        if hasAlign:
            acq = newTs.getAcquisition()
            acq.setTiltAxisAngle(0.)  # 0 because TS is aligned
            newTs.setAcquisition(acq)

        newTs.write(properties=False)
        output.update(newTs)
        output.write()
        self._store()

    def closeOutputSetsStep(self):
        self.TiltSeries.setStreamState(Set.STREAM_CLOSED)
        self.TiltSeries.write()
        self._store()

    # --------------------------- UTILS functions -----------------------------
    def getCtfTomoSeriesFromTsId(self, setOfCtfTomoSeries, tsId):
        for ctfTomoSeries in self.inputSetOfCtfTomoSeries.get():
            if tsId == ctfTomoSeries.getTsId():
                return ctfTomoSeries

    # --------------------------- INFO functions ------------------------------
    def _warnings(self):
        warnings = []
        ts = self.inputSetOfTiltSeries.get()
        if not ts.getFirstItem().getFirstItem().hasTransform():
            warnings.append("Input tilt-series do not have alignment "
                            "information! The recommended workflow is to "
                            "estimate CTF on raw tilt-series and then here "
                            "provide tilt-series with alignment "
                            "(non-interpolated).")

        return warnings

    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append("Input tilt-series: %d\nCTF corrections applied: %d"
                           % (self.inputSetOfCtfTomoSeries.get().getSize(),
                              self.TiltSeries.getSize()))
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.TiltSeries:
            methods.append("%d tilt-series have been CTF corrected "
                           "using the IMOD *ctfphaseflip* program."
                           % (self.TiltSeries.getSize()))
        return methods
