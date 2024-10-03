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
import pyworkflow.protocol.params as params
from imod.protocols.protocol_base import IN_TOMO_SET
from pyworkflow.utils import Message
from tomo.objects import SetOfCoordinates3D, Coordinate3D
import tomo.constants as constants
from imod import utils
from imod.protocols import ProtImodBase
from imod.constants import XYZ_EXT, MOD_EXT, OUTPUT_COORDINATES_3D_NAME

# Beads color options
DARK_BEADS = 0
LIGHT_BEADS = 1


class ProtImodGoldBeadPicker3d(ProtImodBase):
    """
    3-dimensional gold bead picker using the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/findbeads3d.html
    """

    _label = 'Gold bead picker 3D'
    _possibleOutputs = {OUTPUT_COORDINATES_3D_NAME: SetOfCoordinates3D}

    def __init__(self, **args):
        super().__init__(**args)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        form.addParam(IN_TOMO_SET,
                      params.PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Input set of tomograms',
                      help='Input set of tomograms from which gold beads '
                           'will be picked. A tomogram needs to be thicker '
                           'than normal because the program cannot find '
                           'beads too close to the surfaces of a tomogram.')

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

        form.addParallelSection(threads=4, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        allOutputId = []
        self._initialize()
        for tsId in self.tomoDict.keys():
            pickId = self._insertFunctionStep(self.pickGoldBeadsStep, tsId, prerequisites=[])
            outputID = self._insertFunctionStep(self.createOutputStep, tsId, prerequisites=[pickId])
            allOutputId.append(outputID)

        self._insertFunctionStep(self.closeOutputSetsStep, prerequisites=allOutputId)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.tomoDict = {tomo.getTsId(): tomo.clone() for tomo in self.getInputSet()}

    def pickGoldBeadsStep(self, tsId):
        try:
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

            self.runProgram('findbeads3d', paramsFindbeads3d)

            # Run model2point
            paramsModel2Point = {
                "-InputFile": self.getExtraOutFile(tsId, ext=MOD_EXT),
                "-OutputFile": self.getExtraOutFile(tsId, ext=XYZ_EXT),
            }
            self.runProgram('model2point', paramsModel2Point)

        except Exception as e:
            self._failedItems.append(tsId)
            self.error(f"findbeads3d or model2point execution failed for tsId {tsId} -> {e}")

    def createOutputStep(self, tsId):
        tomo = self.tomoDict[tsId]
        with self._lock:
            if tsId in self._failedItems:
                self.createOutputFailedSet(tomo)
            else:
                coordFilePath = self.getExtraOutFile(tsId, ext=XYZ_EXT)
                if os.path.exists(coordFilePath):
                    beadDiam = self.beadDiameter.get()
                    coordList = utils.formatGoldBead3DCoordinatesList(coordFilePath)
                    output = self.getOutputSetOfCoordinates3Ds(self.getInputSet(pointer=True),
                                                               self.getInputSet())
                    output.setBoxSize(beadDiam)

                    for element in coordList:
                        newCoord3D = Coordinate3D()
                        newCoord3D.setVolume(tomo)
                        newCoord3D.setX(element[0], constants.BOTTOM_LEFT_CORNER)
                        newCoord3D.setY(element[1], constants.BOTTOM_LEFT_CORNER)
                        newCoord3D.setZ(element[2], constants.BOTTOM_LEFT_CORNER)

                        output.append(newCoord3D)
                        output.update(newCoord3D)
                    output.write()
                    self._store(output)
                else:
                    self.createOutputFailedSet(tomo)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        coords3D = getattr(self, OUTPUT_COORDINATES_3D_NAME, None)
        if coords3D is not None:
            summary.append(f"Input tomograms: {self.getInputSet().getSize()}\n"
                           "Output coordinates 3D: "
                           f"{coords3D.getSize()}")
        return summary

    # --------------------------- UTILS functions -----------------------------
    def getInputSet(self, pointer=False):
        return (self.inputSetOfTomograms.get() if
                not pointer else self.inputSetOfTomograms)
