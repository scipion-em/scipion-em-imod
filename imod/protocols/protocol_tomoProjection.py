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
import traceback
from os.path import exists
from typing import List

from pwem.convert.headers import setMRCSamplingRate
from pwem.objects import Transform, Pointer
import pyworkflow.protocol.params as params
from pwem.emlib.image import ImageHandler as ih
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.utils import Message, redStr, cyanStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries, CTFTomoSeries, CTFTomo, Tomogram
from imod.protocols import ProtImodBase
from imod.constants import OUTPUT_TILTSERIES_NAME, OUTPUT_CTF_SERIE, XYZPROJ_PROGRAM

logger = logging.getLogger(__name__)


class ProtImodTomoProjection(ProtImodBase):
    """
    Generates synthetic tilt-series projections from reconstructed tomograms
    using the IMOD projection framework. The protocol creates a new series of
    two-dimensional projection images by virtually rotating a tomogram around
    a selected axis and sampling the volume at defined angular intervals. This
    procedure is commonly used in electron tomography workflows to simulate
    projection data, validate reconstruction strategies, benchmark alignment
    procedures, generate synthetic datasets for method development, or study
    the geometric behavior of tomographic reconstructions under controlled
    acquisition conditions. More info:
    https://bio3d.colorado.edu/imod/doc/man/xyzproj.html

    AI Generated:

    Tomogram Projection (ProtImodTomoProjection) — User Manual
        Overview

        The Tomogram Projection protocol generates artificial tilt-series data
        from one or more tomograms by computing projections at user-defined
        rotation angles. In cryo-electron tomography, this type of operation is
        particularly useful when evaluating reconstruction quality, simulating
        acquisition geometries, testing processing pipelines, or preparing
        benchmark datasets for algorithm development.

        Rather than reconstructing a tomogram from projections, this protocol
        performs the inverse conceptual operation: it starts from a 3D volume
        and produces a set of 2D projections corresponding to different tilt
        orientations. The resulting synthetic tilt-series can then be processed
        as if they were experimentally acquired data.

        Inputs and General Workflow

        The protocol requires one or more tomograms as input. These tomograms
        define the three-dimensional density maps that will be projected into
        two-dimensional images. Each input tomogram is processed independently,
        producing a corresponding tilt-series with a defined angular geometry.

        During projection, the user specifies the minimum angle, maximum angle,
        and angular increment. Together, these parameters determine the angular
        coverage and the total number of generated projections. Wide angular
        ranges provide more complete tomographic sampling, while narrower ranges
        may better reproduce realistic experimental limitations.

        The generated projections preserve the sampling characteristics of the
        original tomogram, allowing downstream procedures to operate under
        physically meaningful geometric conditions.

        Choice of Rotation Axis

        One of the most important biological and geometric parameters is the
        selection of the rotation axis. The protocol allows projection around
        the X, Y, or Z axis of the tomogram. In most electron tomography
        experiments, the Y axis typically corresponds to the conventional tilt
        axis used during acquisition, making it the most common choice for
        realistic simulation workflows.

        Rotating around alternative axes can nevertheless be useful for
        visualization studies, geometric validation, educational purposes, or
        algorithm testing under non-standard acquisition geometries. The choice
        of axis directly influences how structural features appear across the
        generated tilt-series.

        Angular Sampling Considerations

        The angular range and angular increment strongly affect the realism and
        interpretability of the generated projections. Small angular increments
        produce denser sampling and smoother reconstruction conditions but also
        increase computational cost and storage requirements. Larger increments
        reduce the number of generated images and may simulate sparse acquisition
        conditions commonly encountered in dose-limited experiments.

        Similarly, the total angular coverage determines the extent of missing
        information in reciprocal space. Limited tilt ranges reproduce the
        missing wedge effect characteristic of experimental tomography, while
        broader ranges provide more isotropic sampling.

        From a biological perspective, realistic angular sampling is important
        when synthetic datasets are intended to mimic actual cryo-electron
        tomography experiments.

        Outputs and Their Interpretation

        The protocol produces one tilt-series for each input tomogram. Each
        generated tilt-series contains projection images associated with their
        corresponding tilt angles and acquisition geometry. The outputs are
        suitable for subsequent alignment, reconstruction, denoising, or
        algorithm benchmarking workflows.

        In addition to the projection images themselves, the protocol also
        generates associated acquisition metadata so that downstream tomography
        procedures can interpret the synthetic datasets consistently. This makes
        the outputs especially useful for validating complete cryo-ET pipelines
        under controlled conditions.

        Practical Recommendations

        For realistic simulation studies, it is generally advisable to select
        angular ranges and increments similar to those used during experimental
        acquisition. Typical cryo-electron tomography datasets often span
        approximately -60 to +60 degrees with increments between 1 and 3
        degrees, although the optimal choice depends on the biological sample
        and imaging strategy.

        When testing reconstruction robustness or alignment algorithms, users
        may intentionally reduce angular coverage or increase angular spacing in
        order to reproduce challenging experimental conditions. Conversely,
        dense angular sampling may be useful for methodological benchmarking or
        visualization purposes.

        The Y-axis rotation mode is usually the preferred configuration for
        generating datasets that resemble standard tomography experiments.
        Alternative axes are more appropriate for specialized geometric studies
        or software validation tasks.

        Final Perspective

        Synthetic projection generation is an important methodological tool in
        cryo-electron tomography because it enables controlled experimentation
        without the variability of experimental acquisition. By generating
        tilt-series directly from reconstructed tomograms, this protocol allows
        researchers to evaluate algorithms, reproduce acquisition geometries,
        and better understand the relationship between three-dimensional
        biological structures and their two-dimensional projections.
    """

    _label = 'Tomo projection'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}
    stepsExecutionMode = STEPS_PARALLEL

    AXIS_X = 0
    AXIS_Y = 1
    AXIS_Z = 2

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tomoDict = None

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
        self._initialize()
        for tsId in self.tomoDict.keys():
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
    def _initialize(self):
        self.tomoDict = {tomo.getTsId(): tomo.clone() for tomo in self.getInputTomoSet().iterItems()}

    def projectTomogram(self, tsId: str):
        try:
            logger.info(cyanStr(f'tsId = {tsId} -> Projecting the tomogram...'))
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

            self.runProgram(XYZPROJ_PROGRAM, paramsXYZproj)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {XYZPROJ_PROGRAM} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def generateOutputStackStep(self, tsId: str):
        tomo = self.tomoDict[tsId]
        if tsId in self.failedItems:
            self.addToOutFailedSet(tomo)
            return

        try:
            outputFn = self.getExtraOutFile(tsId)
            if exists(outputFn):
                setMRCSamplingRate(outputFn, tomo.getSamplingRate())  # Update the apix value in file header
                self._registerOutput(tomo, outputFn)
            else:
                logger.error(redStr(f'tsId = {tsId} -> Output file {outputFn} was not generated. Skipping... '))

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    @retry_on_sqlite_lock(log=logger)
    def _registerOutput(self, tomo: Tomogram, outputFn = str):
        tsId = tomo.getTsId()
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
            # Generate fake CTFs
        ctfSet = self.getOutputSetOfCTFTomoSeries(Pointer(self, extended=OUTPUT_TILTSERIES_NAME),
                                                  OUTPUT_CTF_SERIE)
        ctfSerie = CTFTomoSeries()
        ctfSerie.copyInfo(outTs)
        ctfSerie.setTiltSeries(outTs)
        ctfSerie.setTsId(tsId)
        ctfSet.append(ctfSerie)

        tiltAngleList = self.getTiltAngleList()
        sRate = tomo.getSamplingRate()
        for slice in range(1, self.getProjectionRange() + 1):
            newTi = TiltImage(tsId=tsId,
                              tiltAngle=tiltAngleList[slice - 1],
                              acquisitionOrder=slice)
            newTi.setLocation(slice, outputFn)
            newTi.setSamplingRate(sRate)
            newTi.setAcquisition(acq)
            outTs.append(newTi)

            # Fake CTF
            tiltCTF = CTFTomo(index=slice, acqOrder=slice)
            tiltCTF.setDefocusU(0)
            tiltCTF.setDefocusV(0)
            tiltCTF.setDefocusAngle(0)
            tiltCTF.setResolution(1)
            ctfSerie.append(tiltCTF)

        x, y, z, _ = ih.getDimensions(outputFn)
        outTs.setDim((x, y, z))
        outTs.setAnglesCount(len(outTs))

        # Set origin to output tilt-series
        origin = Transform()
        origin.setShifts(x / -2. * sRate, y / -2. * sRate, 0)
        outTs.setOrigin(origin)

        # Data persistence
        outTsSet.update(outTs)
        self._store(outTsSet)
        ctfSet.update(ctfSerie)
        self._store(outTsSet, ctfSet)

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
