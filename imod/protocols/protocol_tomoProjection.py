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

from pwem.objects import Transform
from pyworkflow import BETA
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
from imod import Plugin
from pwem.emlib.image import ImageHandler


class ProtImodTomoProjection(EMProtocol, ProtTomoBase):
    """
    Re-project a tomogram given a geometric description (axis and angles).
    More info:
        https://bio3d.colorado.edu/imod/doc/man/xyzproj.html
    """

    _label = 'tomo projection'
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
                      label='Input set of tomograms')

        form.addParam('minAngle',
                      params.FloatParam,
                      default=-60.0,
                      label='Minimum angle of rotation',
                      important=True,
                      help='Minimum angle of the projection range')

        form.addParam('maxAngle',
                      params.FloatParam,
                      default=60.0,
                      label='Maximum angle of rotation',
                      important=True,
                      help='Maximum angle of the projection range')

        form.addParam('stepAngle',
                      params.FloatParam,
                      default=2.0,
                      label='Step angle of rotation',
                      important=True,
                      help='Step angle of the projection range')

        form.addParam('rotationAxis',
                      params.EnumParam,
                      choices=['X', 'Y', 'Z'],
                      default=self.AXIS_Y,
                      label='Rotation axis for projection',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Axis to tilt around (X, Y, or Z). Y axis usually corresponds to the typical rotation axis '
                           'acquisition.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for tomo in self.inputSetOfTomograms.get():
            self._insertFunctionStep('projectTomogram', tomo.getObjId())
            self._insertFunctionStep('generateOutputStackStep', tomo.getObjId())
        self._insertFunctionStep('closeOutputSetsStep')

    # --------------------------- STEPS functions ----------------------------
    def projectTomogram(self, tomoObjId):
        tomo = self.inputSetOfTomograms.get()[tomoObjId]

        tomoId = os.path.splitext(os.path.basename(tomo.getFileName()))[0]

        extraPrefix = self._getExtraPath(tomoId)
        tmpPrefix = self._getTmpPath(tomoId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        paramsXYZproj = {
            'input': tomo.getFileName(),
            'output': os.path.join(extraPrefix, os.path.basename(tomo.getFileName())),
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

    def generateOutputStackStep(self, tomoObjId):
        tomo = self.inputSetOfTomograms.get()[tomoObjId]

        tomoId = os.path.splitext(os.path.basename(tomo.getFileName()))[0]

        extraPrefix = self._getExtraPath(tomoId)

        outputProjectedSetOfTiltSeries = self.getOutputProjectedSetOfTiltSeries()
        
        newTs = tomoObj.TiltSeries(tsId=tomoId)
        newTs.copyInfo(tomo)
        newTs.setTsId(tomoId)

        # Add origin to output tilt-series
        origin = Transform()

        outputProjectedSetOfTiltSeries.append(newTs)

        tiltAngleList = self.getTiltAngleList()

        for index in range(self.getProjectionRange()):
            newTi = tomoObj.TiltImage()
            newTi.setTiltAngle(tiltAngleList[index])
            newTi.setTsId(tomoId)
            newTi.setLocation(index + 1, os.path.join(extraPrefix, os.path.basename(tomo.getFileName())))
            newTi.setSamplingRate(self.inputSetOfTomograms.get().getSamplingRate())
            newTs.append(newTi)

        ih = ImageHandler()
        x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
        newTs.setDim((x, y, z))

        # Set origin to output tilt-series
        origin.setShifts(x / -2. * self.inputSetOfTomograms.get().getSamplingRate(),
                         y / -2. * self.inputSetOfTomograms.get().getSamplingRate(),
                         0)
        newTs.setOrigin(origin)

        newTs.write(properties=False)

        outputProjectedSetOfTiltSeries.update(newTs)
        outputProjectedSetOfTiltSeries.updateDim()
        outputProjectedSetOfTiltSeries.write()
        self._store()

    def closeOutputSetsStep(self):
        self.getOutputProjectedSetOfTiltSeries().setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputProjectedSetOfTiltSeries(self):
        if hasattr(self, "outputProjectedSetOfTiltSeries"):
            self.outputProjectedSetOfTiltSeries.enableAppend()
        else:
            outputProjectedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Projected')
            outputProjectedSetOfTiltSeries.copyInfo(self.inputSetOfTomograms.get())
            outputProjectedSetOfTiltSeries.setDim(self.inputSetOfTomograms.get().getDim())
            outputProjectedSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputProjectedSetOfTiltSeries=outputProjectedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTomograms, outputProjectedSetOfTiltSeries)
        return self.outputProjectedSetOfTiltSeries

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

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        if self.minAngle.get() > self.maxAngle.get():
            validateMsgs.append("ERROR: Maximum angle of rotation bigger than minimum.")

        if (self.maxAngle.get() - self.minAngle.get()) < self.stepAngle.get():
            validateMsgs.append("ERROR: Angle step of rotation bigger than range.")

        if self.stepAngle.get() < 0:
            validateMsgs.append("ERROR: Angle step of rotation mus be bigger than zero.")

        return validateMsgs

    def _summary(self):
        summary = []
        if hasattr(self, 'outputProjectedSetOfTiltSeries'):
            summary.append("Input Tomograms: %d.\nTilt-series generated: %d.\n"
                           % (self.inputSetOfTomograms.get().getSize(),
                              self.outputProjectedSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputProjectedSetOfTiltSeries'):
            methods.append("%d tilt-series have been generated by projecting the input tomogram.\n"
                           % (self.outputProjectedSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
