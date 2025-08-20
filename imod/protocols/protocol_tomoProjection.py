# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
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
import logging
from os.path import exists
from typing import List
from pwem.objects import Transform
import pyworkflow.protocol.params as params
from pwem.emlib.image import ImageHandler as ih
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.utils import Message, redStr, cyanStr
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries
from imod.protocols import ProtImodBase
from imod.constants import OUTPUT_TILTSERIES_NAME, XYZPROJ_PROGRAM

logger = logging.getLogger(__name__)


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
    stepsExecutionMode = STEPS_PARALLEL

    AXIS_X = 0
    AXIS_Y = 1
    AXIS_Z = 2

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        super().addInTomoSetFormParam(form)
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
        form.addParallelSection(threads=1, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        closeSetStepDeps = []
        for tomo in self.getInputTomoSet():
            tsId = tomo.getTsId()
            compId = self._insertFunctionStep(self.projectTomogram,
                                              tsId,
                                              prerequisites=[],
                                              needsGPU=False)
            outId = self._insertFunctionStep(self.generateOutputStackStep,
                                             tsId,
                                             prerequisites=compId,
                                             needsGPU=False)
            closeSetStepDeps.append(outId)

        self._insertFunctionStep(self.closeOutputSetsStep,
                                 OUTPUT_TILTSERIES_NAME,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def projectTomogram(self, tsId: str):
        try:
            logger.info(cyanStr(f'tsId = {tsId} -> Projecting the tomogram...'))
            self.genTsPaths(tsId)
            with self._lock:
                tomo = self.getCurrentTomo(tsId)

            paramsXYZproj = {
                '-input': tomo.getFileName(),
                '-output': self.getExtraOutFile(tsId),
                '-axis': self.getRotationAxis(),
                '-angles': ",".join(map(str, [self.minAngle.get(),
                                              self.maxAngle.get(),
                                              self.stepAngle.get()]))
            }

            self.runProgram(XYZPROJ_PROGRAM, paramsXYZproj)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {XYZPROJ_PROGRAM} execution '
                                f'failed with the exception -> {e}'))

    def generateOutputStackStep(self, tsId: str):
        if tsId in self.failedItems:
            self.addToOutFailedSet(tsId, inputsAreTs=False)
        else:
            try:
                outputFn = self.getExtraOutFile(tsId)
                if exists(outputFn):
                    tomo = self.getCurrentTomo(tsId)
                    with self._lock:
                        # Set of tilt-series
                        outTsSet = self.getOutputSetOfTS(self.getInputTomoSet(pointer=True))
                        # Tilt-series
                        outTs = TiltSeries(tsId=tsId)
                        acq = tomo.getAcquisition()
                        acq.setAngleMin(self.minAngle.get())
                        acq.setAngleMax(self.maxAngle.get())
                        acq.setStep(self.stepAngle.get())
                        acq.setTiltAxisAngle(0)
                        outTs.setAcquisition(acq)
                        outTsSet.append(outTs)
                        # Tilt-images
                        tiltAngleList = self.getTiltAngleList()
                        sRate = tomo.getSamplingRate()
                        for index in range(self.getProjectionRange()):
                            newTi = TiltImage(tsId=tsId,
                                              tiltAngle=tiltAngleList[index],
                                              acquisitionOrder=index + 1)
                            newTi.setLocation(index + 1, outputFn)
                            newTi.setSamplingRate(sRate)
                            newTi.setAcquisition(acq)
                            outTs.append(newTi)

                        x, y, z, _ = ih.getDimensions(outputFn)
                        outTs.setDim((x, y, z))
                        outTs.setAnglesCount(len(outTs))

                        # Set origin to output tilt-series
                        origin = Transform()
                        origin.setShifts(x / -2. * sRate, y / -2. * sRate, 0)
                        outTs.setOrigin(origin)

                        # Data persistence
                        outTs.write()
                        outTsSet.update(outTs)
                        outTsSet.write()
                        self._store(outTsSet)
                else:
                    logger.error(redStr(f'tsId = {tsId} -> Output file {outputFn} was not generated. Skipping... '))
            except Exception as e:
                logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))

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
            summary.append(f"Input tomograms: {self.getInputTsSet().getSize()}\n"
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
    def getRotationAxis(self) -> str:
        parseParamsRotationAxis = {
            self.AXIS_X: 'X',
            self.AXIS_Y: 'Y',
            self.AXIS_Z: 'Z',
        }
        return parseParamsRotationAxis[self.rotationAxis.get()]

    def getProjectionRange(self):
        return int((self.maxAngle.get() - self.minAngle.get()) / self.stepAngle.get()) + 1

    def getTiltAngleList(self) -> List[float]:
        tiltAngleList = []

        angle = self.minAngle.get()
        while angle <= self.maxAngle.get():
            tiltAngleList.append(angle)
            angle += self.stepAngle.get()

        return tiltAngleList
