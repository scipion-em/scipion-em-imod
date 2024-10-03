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

from imod.protocols.protocol_base import IN_TOMO_SET
from pwem.objects import Transform
import pyworkflow.protocol.params as params
from pwem.emlib.image import ImageHandler as ih
from pyworkflow.utils import Message
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries

from imod.protocols import ProtImodBase
from imod.constants import OUTPUT_TILTSERIES_NAME


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
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}

    AXIS_X = 0
    AXIS_Y = 1
    AXIS_Z = 2

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        form.addParam(IN_TOMO_SET,
                      params.PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Input set of tomograms to be projected')

        line = form.addLine('Tilt angles  (deg)',
                            help='Starting, ending, and increment tilt angle. '
                                 'Enter the same value for starting and '
                                 'ending angle to get only one image.')

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

        form.addParallelSection(threads=4, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []
        for tsId in self.tomoDict.keys():
            compId = self._insertFunctionStep(self.projectTomogram, tsId, prerequisites=[])
            outId = self._insertFunctionStep(self.generateOutputStackStep, tsId, prerequisites=[compId])
            closeSetStepDeps.append(outId)

        self._insertFunctionStep(self.closeOutputSetsStep, prerequisites=closeSetStepDeps)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.tomoDict = {tomo.getTsId(): tomo.clone() for tomo in self.inputSetOfTomograms.get()}

    def projectTomogram(self, tsId):
        try:
            self.genTsPaths(tsId)
            tomo = self.tomoDict[tsId]

            paramsXYZproj = {
                '-input': tomo.getFileName(),
                '-output': self.getExtraOutFile(tsId),
                '-axis': self.getRotationAxis(),
                '-angles': ",".join(map(str, [self.minAngle.get(),
                                              self.maxAngle.get(),
                                              self.stepAngle.get()]))
            }

            self.runProgram('xyzproj', paramsXYZproj)

        except Exception as e:
            self._failedItems.append(tsId)
            self.error(f'xyzproj execution failed for tsId {tsId} -> {e}')

    def generateOutputStackStep(self, tsId):
        tomo = self.tomoDict[tsId]
        with self._lock:
            if tsId in self._failedItems:
                self.createOutputFailedSet(tomo)
            else:
                outputFn = self.getExtraOutFile(tsId)
                if os.path.exists(outputFn):
                    output = self.getOutputSetOfTS(self.getInputSet(pointer=True))
                    newTs = TiltSeries(tsId=tsId)
                    acq = tomo.getAcquisition()
                    acq.setAngleMin(self.minAngle.get())
                    acq.setAngleMax(self.maxAngle.get())
                    acq.setStep(self.stepAngle.get())
                    acq.setTiltAxisAngle(0)
                    newTs.setAcquisition(acq)
                    output.append(newTs)

                    tiltAngleList = self.getTiltAngleList()
                    sRate = tomo.getSamplingRate()
                    for index in range(self.getProjectionRange()):
                        newTi = TiltImage(tsId=tsId,
                                          tiltAngle=tiltAngleList[index],
                                          acquisitionOrder=index + 1)
                        newTi.setLocation(index + 1, outputFn)
                        newTi.setSamplingRate(sRate)
                        newTi.setAcquisition(acq)
                        newTs.append(newTi)

                    x, y, z, _ = ih.getDimensions(outputFn)
                    newTs.setDim((x, y, z))
                    newTs.setAnglesCount(len(newTs))

                    # Set origin to output tilt-series
                    origin = Transform()
                    origin.setShifts(x / -2. * sRate, y / -2. * sRate, 0)
                    newTs.setOrigin(origin)

                    output.update(newTs)
                    self._store(output)
                else:
                    self.createOutputFailedSet(tomo)

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
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            summary.append(f"Input tomograms: {self.getInputSet().getSize()}\n"
                           f"Tilt-series generated: {output.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            methods.append(f"{output.getSize()} tilt-series have been "
                           "generated by projecting the input tomogram using "
                           "IMOD *xyzproj* command.")
        return methods

    # --------------------------- UTILS functions -----------------------------
    def getInputSet(self, pointer=False):
        return self.inputSetOfTomograms.get() if not pointer else self.inputSetOfTomograms

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
