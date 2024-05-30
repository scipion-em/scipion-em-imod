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

from pwem.objects import Transform
from pyworkflow import BETA
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
from pwem.emlib.image import ImageHandler
import tomo.objects as tomoObj

from .. import Plugin
from .protocol_base import ProtImodBase, MRCS_EXT


class ProtImodTomoProjection(ProtImodBase):
    """
    Re-project a tomogram given a geometric description (axis and angles).
    More info:
        https://bio3d.colorado.edu/imod/doc/man/xyzproj.html

    This program will compute projections of a tomogram at a series of
    tilts around either the X, the Y or the Z axis.\n

    A projection along a ray line is simply the average of the pixels in
    the block along that line.  However, rather than taking the values of
    the pixels that lie near the ray, interpolation is used to sample den-
    sity at points evenly spaced at one pixel intervals along the ray.
    """

    _label = 'Tomo projection'
    _devStatus = BETA

    AXIS_X = 0
    AXIS_Y = 1
    AXIS_Z = 2

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTomograms',
                      params.PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Input set of tomograms to be projected')

        line = form.addLine('Tilt angles  (deg)',
                            help='Starting, ending, and increment tilt angle.  Enter the same value for '
                                 'starting and ending angle to get only one image')

        line.addParam('minAngle',
                      params.FloatParam,
                      default=-60.0,
                      label='Minimum rotation')

        line.addParam('maxAngle',
                      params.FloatParam,
                      default=60.0,
                      label='Maximum rotation')

        line.addParam('stepAngle',
                      params.FloatParam,
                      default=2.0,
                      label='Step angle')

        form.addParam('rotationAxis',
                      params.EnumParam,
                      choices=['X', 'Y', 'Z'],
                      default=self.AXIS_Y,
                      label='Rotation axis for projection',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Axis to tilt around (X, Y, or Z). Y axis usually '
                           'corresponds to the typical rotation axis '
                           'acquisition.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        for tsId in self.tomoDict.keys():
            self._insertFunctionStep(self.projectTomogram, tsId)
            self._insertFunctionStep(self.generateOutputStackStep, tsId)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.tomoDict = {tomo.getTsId(): tomo.clone() for tomo in self.inputSetOfTomograms.get()}

    def projectTomogram(self, tsId):
        self.genTsPaths(tsId)
        tomo = self.tomoDict[tsId]

        paramsXYZproj = {
            'input': tomo.getFileName(),
            'output': self.getExtraOutFile(tsId, ext=MRCS_EXT),
            'axis': self.getRotationAxis(),
            'angles': str(self.minAngle.get()) + ',' +
                      str(self.maxAngle.get()) + ',' +
                      str(self.stepAngle.get()),
        }

        argsXYZproj = "-input %(input)s " \
                      "-output %(output)s " \
                      "-axis %(axis)s " \
                      "-angles %(angles)s "

        Plugin.runImod(self, 'xyzproj', argsXYZproj % paramsXYZproj)

    def generateOutputStackStep(self, tsId):
        tomo = self.tomoDict[tsId]
        output = self.getOutputSetOfTiltSeries(self.inputSetOfTomograms.get())
        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.setTsId(tomo.getTsId())
        newTs.setAcquisition(tomo.getAcquisition())

        # Add origin to output tilt-series
        output.append(newTs)
        tiltAngleList = self.getTiltAngleList()
        sRate = self.inputSetOfTomograms.get().getSamplingRate()
        for index in range(self.getProjectionRange()):
            newTi = tomoObj.TiltImage()
            newTi.setTiltAngle(tiltAngleList[index])
            newTi.setTsId(tsId)
            newTi.setAcquisitionOrder(index + 1)
            newTi.setLocation(index + 1,
                              self.getExtraOutFile(tsId, ext=MRCS_EXT))
            newTi.setSamplingRate(sRate)
            newTs.append(newTi)

        ih = ImageHandler()
        x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
        newTs.setDim((x, y, z))

        # Set origin to output tilt-series
        origin = Transform()
        origin.setShifts(x / -2. * self.inputSetOfTomograms.get().getSamplingRate(),
                         y / -2. * self.inputSetOfTomograms.get().getSamplingRate(),
                         0)

        newTs.setOrigin(origin)
        newTs.write(properties=False)

        output.update(newTs)
        output.write()
        self._store()

    def closeOutputSetsStep(self):
        self.TiltSeries.setStreamState(Set.STREAM_CLOSED)
        self.TiltSeries.write()
        self._store()

    # --------------------------- UTILS functions -----------------------------
    def getRotationAxis(self):
        parseParamsRotationAxis = {
            self.AXIS_X: 'X',
            self.AXIS_Y: 'Y',
            self.AXIS_Z: 'Z',
        }
        return parseParamsRotationAxis[self.rotationAxis.get()]

    def getProjectionRange(self):
        return int((self.maxAngle.get() - self.minAngle.get()) / self.stepAngle.get()) + 1

    def getTiltAngleList(self):
        tiltAngleList = []

        angle = self.minAngle.get()
        while angle <= self.maxAngle.get():
            tiltAngleList.append(angle)
            angle += self.stepAngle.get()

        return tiltAngleList

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        if self.minAngle.get() > self.maxAngle.get():
            validateMsgs.append("ERROR: Maximum angle of rotation bigger "
                                "than minimum.")

        if (self.maxAngle.get() - self.minAngle.get()) < self.stepAngle.get():
            validateMsgs.append("ERROR: Angle step of rotation bigger than range.")

        if self.stepAngle.get() < 0:
            validateMsgs.append("ERROR: Angle step of rotation must be "
                                "bigger than zero.")

        return validateMsgs

    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append("Input tomograms: %d\nTilt-series generated: %d"
                           % (self.inputSetOfTomograms.get().getSize(),
                              self.TiltSeries.getSize()))
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.TiltSeries:
            methods.append("%d tilt-series have been generated by projecting "
                           "the input tomogram using IMOD *xyzproj* command.\n"
                           % (self.TiltSeries.getSize()))
        return methods
