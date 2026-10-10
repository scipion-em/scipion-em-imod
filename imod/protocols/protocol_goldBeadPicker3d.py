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

import pyworkflow.protocol.params as params
from imod.protocols.protocol_fiducialModel import MODEL2POINT_PROGRAM
from pyworkflow.object import Pointer, Set
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.utils import Message, cyanStr, redStr
from tomo.objects import SetOfCoordinates3D, Coordinate3D
import tomo.constants as constants
from imod.protocols import ProtImodBase
from imod.constants import XYZ_EXT, MOD_EXT, OUTPUT_COORDINATES_3D_NAME

logger = logging.getLogger(__name__)

FINDBEADS3D_PROGRAM = 'findbeads3d'

# Beads color options
DARK_BEADS = 0
LIGHT_BEADS = 1


class ProtImodGoldBeadPicker3d(ProtImodBase):
    """
    Detects and extracts three-dimensional fiducial gold bead coordinates
    from tomograms using the IMOD fiducial detection workflow. More info:
    https://bio3d.colorado.edu/imod/doc/man/findbeads3d.html

    AI Generated:

    Gold Bead Picker 3D (ProtImodGoldBeadPicker3d) — User Manual

        Overview

        The Gold Bead Picker 3D protocol identifies fiducial gold beads
        within reconstructed tomograms and generates a set of
        three-dimensional coordinates representing their spatial
        positions. Fiducial markers are commonly used in cryo-electron
        tomography experiments to provide stable reference points for
        tilt-series alignment and geometric calibration. Accurate bead
        detection is therefore essential for reliable tomographic
        reconstruction and downstream structural interpretation.

        In practical cryo-ET workflows, gold beads act as artificial
        landmarks embedded within the specimen support. Because these
        particles produce strong and easily recognizable contrast, they
        can be tracked consistently across tilt images and reconstructed
        volumes. The protocol is designed to automate this detection
        process in reconstructed tomograms, reducing manual annotation
        effort and improving reproducibility in large datasets.

        Biological Context and Fiducial Markers

        Fiducial gold particles are not part of the biological sample
        itself but serve as geometric references during image alignment.
        Their proper identification directly influences the quality of
        tilt-series registration, tomogram reconstruction, and spatial
        consistency across the dataset. Inaccurate fiducial detection
        may propagate into reconstruction artifacts, distorted cellular
        structures, or reduced resolution in subtomogram averaging
        workflows.

        In many biological experiments, fiducials are distributed on
        both sides of the specimen support film to improve geometric
        stability during alignment. Ideally, beads should be well
        separated, clearly visible, and evenly distributed across the
        field of view. Aggregated or overlapping fiducials can reduce
        alignment accuracy and complicate automatic detection.

        Inputs and Detection Strategy

        The protocol requires reconstructed tomograms as input together
        with an approximate fiducial diameter expressed in pixels. This
        parameter is biologically important because it determines the
        expected scale of the fiducial particles inside the tomographic
        volume. Selecting a value close to the true bead diameter
        improves robustness and reduces the likelihood of detecting
        contaminants or structural noise.

        The protocol also allows users to define whether fiducials
        appear as dark particles on a bright background or as bright
        particles on a darker environment. Correct contrast selection
        is essential because tomographic datasets may differ depending
        on acquisition conditions, reconstruction procedures, or image
        inversion conventions.

        Detection Sensitivity and Peak Selection

        The fiducial detection process relies on identifying local
        density features consistent with spherical gold particles. Users
        can control the sensitivity of this process through parameters
        that regulate the minimum accepted signal strength and the
        minimum spacing between detected particles.

        Lower detection thresholds may recover weak or partially
        obscured fiducials but can also increase false positives caused
        by contamination, ice artifacts, or dense biological material.
        Higher thresholds generally improve specificity but may miss
        faint beads in low-contrast regions of the tomogram.

        The spacing parameter is especially relevant in crowded
        fiducial environments. Closely spaced detections may correspond
        either to neighboring beads or to duplicate detections of the
        same particle. Appropriate spacing constraints help maintain
        biologically meaningful fiducial distributions while avoiding
        redundant coordinate assignments.

        Outputs and Workflow Integration

        After execution, the protocol generates a set of three-
        dimensional fiducial coordinates associated with the input
        tomograms. These coordinates can be used for visualization,
        alignment validation, geometric inspection, or integration into
        additional tomography processing workflows.

        The resulting coordinate sets preserve the spatial relationship
        between fiducials and their corresponding tomograms, allowing
        direct interpretation within reconstruction and visualization
        environments. In iterative alignment workflows, the detected
        coordinates may also support refinement procedures or quality
        control analyses.

        Practical Recommendations

        In routine cryo-electron tomography processing, users should
        begin with fiducial diameter values that closely reflect the
        nominal bead size used during sample preparation. Visual
        inspection of the detected coordinates is strongly recommended
        to confirm that fiducials are correctly identified and that
        contaminants or biological structures are not mistakenly
        selected.

        Datasets containing highly crowded fiducials, strong ice
        contamination, or low contrast may require adjustment of
        sensitivity parameters to achieve stable results. In some
        challenging biological specimens, balancing detection
        sensitivity and specificity becomes essential for maintaining
        reconstruction quality.

        It is also advisable to verify that fiducials are spatially
        distributed throughout the tomogram volume rather than
        concentrated in a limited region. Well-distributed fiducials
        generally provide more reliable geometric support for alignment
        and reconstruction procedures.

        Final Perspective

        Automatic fiducial detection is a foundational step in many
        cryo-electron tomography workflows because it establishes the
        geometric reference framework used throughout tomographic
        reconstruction. Reliable gold bead identification improves
        alignment precision, enhances reconstruction consistency, and
        supports accurate biological interpretation of three-
        dimensional cellular and molecular structures.
    """

    _label = 'Gold bead picker 3D'
    _possibleOutputs = {OUTPUT_COORDINATES_3D_NAME: SetOfCoordinates3D}
    stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        super().addInTomoSetFormParam(form)
        form.addParam('beadDiameter',
                      params.FloatParam,
                      label='Fiducials diameter (px)',
                      default=18,
                      help="Diameter of beads in pixels.")
        form.addParam('beadsColor',
                      params.EnumParam,
                      choices=['Dark', 'Light'],
                      label='Bead contrast',
                      default='0',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Contrast of the gold beads:\n'
                           '-Dark: beads are dark on light background.\n'
                           '-Light: beads are light on dark background.')
        form.addParam('minRelativeStrength',
                      params.FloatParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Minimum relative strength',
                      default=0.05,
                      help='Minimum relative peak strength for keeping a '
                           'peak in the analysis.  The square root of the '
                           'specified value is used for comparing with the '
                           'square root of peak strength, for compatibility '
                           'with existing command files. The default is 0.05, '
                           'which corresponds to a relative square root peak '
                           'strength of 0.22. Too many weak peaks can prevent '
                           'a dip from showing up in the smoothed histogram '
                           'of strengths.  If the program fails to find a '
                           'histogram dip, one strategy is to try raising '
                           'this value.')
        form.addParam('minSpacing',
                      params.FloatParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Minimum spacing',
                      default=0.9,
                      help='Minimum spacing between peaks as a fraction of '
                           'the bead size. When two peaks are closer '
                           'than this distance apart, the weaker one is '
                           'eliminated unless the -both option is entered. '
                           'The default is 0.9. A value less than 1 is '
                           'helpful for picking both beads in a pair.')
        form.addParallelSection(threads=1, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        allOutputId = []
        self._initialize()
        for tsId in self.tomoDict.keys():
            pickId = self._insertFunctionStep(self.pickGoldBeadsStep,
                                              tsId,
                                              prerequisites=[])
            outputID = self._insertFunctionStep(self.createOutputStep,
                                                tsId,
                                                prerequisites=pickId)
            allOutputId.append(outputID)

        self._insertFunctionStep(self.closeOutputSetsStep,
                                 OUTPUT_COORDINATES_3D_NAME,
                                 prerequisites=allOutputId)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.tomoDict = {tomo.getTsId(): tomo.clone() for tomo in self.getInputTomoSet()}

    def pickGoldBeadsStep(self, tsId):
        try:
            logger.info(cyanStr(f'tsId = {tsId} -> Finding the gold beads...'))
            self.genTsPaths(tsId)
            tomo = self.tomoDict[tsId]

            # Run findbeads3d
            paramsFindbeads3d = {
                "-InputFile": tomo.getFileName(),
                "-OutputFile": self.getExtraOutFile(tsId, ext=MOD_EXT),
                "-BeadSize": self.beadDiameter.get(),
                "-MinRelativeStrength": self.minRelativeStrength.get(),
                "-MinSpacing": self.minSpacing.get(),
                "-StorageThreshold": 0.0
            }

            if self.beadsColor.get() == 1:
                paramsFindbeads3d["-LightBeads"] = ""

            self.runProgram(FINDBEADS3D_PROGRAM, paramsFindbeads3d)

            # Run model2point
            paramsModel2Point = {
                "-InputFile": self.getExtraOutFile(tsId, ext=MOD_EXT),
                "-OutputFile": self.getExtraOutFile(tsId, ext=XYZ_EXT),
            }
            self.runProgram(MODEL2POINT_PROGRAM, paramsModel2Point)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {FINDBEADS3D_PROGRAM} or {MODEL2POINT_PROGRAM} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def createOutputStep(self, tsId):
        tomo = self.tomoDict[tsId]
        if tsId in self.failedItems:
            self.addToOutFailedSet(tomo)
        else:
            try:
                coordFilePath = self.getExtraOutFile(tsId, ext=XYZ_EXT)
                if exists(coordFilePath):
                    beadDiam = self.beadDiameter.get()
                    coordList = self.formatGoldBead3DCoordinatesList(coordFilePath)
                    output = self.getOutputSetOfCoordinates3Ds(self.getInputTomoSet(pointer=True))
                    output.setBoxSize(beadDiam)

                    for element in coordList:
                        newCoord3D = Coordinate3D()
                        newCoord3D.setVolume(tomo)
                        newCoord3D.setX(element[0], constants.BOTTOM_LEFT_CORNER)
                        newCoord3D.setY(element[1], constants.BOTTOM_LEFT_CORNER)
                        newCoord3D.setZ(element[2], constants.BOTTOM_LEFT_CORNER)

                        output.append(newCoord3D)
                        output.update(newCoord3D)

                    # Data persistence
                    output.write()
                    self._store(output)
                else:
                    logger.error(redStr(f'tsId = {tsId} -> Output file {coordFilePath} was not generated. Skipping... '))
            except Exception as e:
                logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
                logger.error(traceback.format_exc())

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        coords3D = getattr(self, OUTPUT_COORDINATES_3D_NAME, None)
        if coords3D is not None:
            summary.append(f"Input tomograms: {self.getInputTsSet().getSize()}\n"
                           "Output coordinates 3D: "
                           f"{coords3D.getSize()}")
        return summary

    # --------------------------- UTILS functions -----------------------------
    def getOutputSetOfCoordinates3Ds(self,
                                     inTomosSetPointer: Pointer) -> SetOfCoordinates3D:
        coords3D = getattr(self, OUTPUT_COORDINATES_3D_NAME, None)
        if coords3D is not None:
            coords3D.enableAppend()
        else:
            inTomoSet = inTomosSetPointer.get()
            coords3D = self._createSetOfCoordinates3D(volSet=inTomoSet,
                                                      suffix='Fiducials3D')
            coords3D.setSamplingRate(inTomoSet.getSamplingRate())
            coords3D.setPrecedents(inTomoSet)
            coords3D.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_COORDINATES_3D_NAME: coords3D})
            self._defineSourceRelation(inTomosSetPointer, coords3D)

        return coords3D

    @staticmethod
    def formatGoldBead3DCoordinatesList(coordFilePath: str) -> List[List[float]]:
        """This method takes the IMOD-based gold bead 3D
        coordinates obtained with find3dbeads program file
        path and returns a list containing each coordinate
        for each bead belonging to the tilt-series"""

        coorList = []

        with open(coordFilePath) as f:
            coorText = f.read().splitlines()

            for line in coorText:
                vector = line.split()
                coorList.append([float(vector[0]), float(vector[1]), float(vector[2])])

        return coorList



