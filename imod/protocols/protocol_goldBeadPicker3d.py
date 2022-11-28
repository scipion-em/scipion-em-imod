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
import pyworkflow.protocol.params as params
from pyworkflow.utils import path, removeBaseExt
from pyworkflow.object import Set
import tomo.objects as tomoObj
import tomo.constants as constants

from .. import Plugin, utils
from .protocol_base import ProtImodBase


class ProtImodGoldBeadPicker3d(ProtImodBase):
    """
    3-dimensional gold bead picker using the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/findbeads3d.html
    """

    _label = 'Gold bead picker 3D'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTomograms',
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
                      default='18',
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

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self.defineExecutionPararell()

        allOutputId = []

        for ts in self.inputSetOfTomograms.get():
            pickId = self._insertFunctionStep(self.pickGoldBeadsStep,
                                              ts.getObjId(),
                                              prerequisites=[])

            convertId = self._insertFunctionStep(self.convertModelToCoordinatesStep,
                                                 ts.getObjId(),
                                                 prerequisites=[pickId])

            outputID = self._insertFunctionStep(self.createOutputStep,
                                                ts.getObjId(),
                                                prerequisites=[convertId])

            allOutputId.append(outputID)

        self._insertFunctionStep('closeOutputSetStep',
                                 prerequisites=allOutputId)

    # --------------------------- STEPS functions -----------------------------
    def pickGoldBeadsStep(self, tsObjId):
        tomo = self.inputSetOfTomograms.get()[tsObjId]
        fileName = removeBaseExt(tomo.getFileName())
        extraPrefix = self._getExtraPath(fileName)
        path.makePath(extraPrefix)

        """ Run findbeads3d IMOD program """
        paramsFindbeads3d = {
            'inputFile': tomo.getFileName(),
            'outputFile': os.path.join(extraPrefix, "%s.mod" % fileName),
            'beadSize': self.beadDiameter.get(),
            'minRelativeStrength': self.minRelativeStrength.get(),
            'minSpacing': self.minSpacing.get(),
        }

        argsFindbeads3d = "-InputFile %(inputFile)s " \
                          "-OutputFile %(outputFile)s " \
                          "-BeadSize %(beadSize)d " \
                          "-MinRelativeStrength %(minRelativeStrength)f " \
                          "-StorageThreshold 0.0 " \
                          "-MinSpacing %(minSpacing)f "

        if self.beadsColor.get() == 1:
            argsFindbeads3d += "-LightBeads "

        Plugin.runImod(self, 'findbeads3d', argsFindbeads3d % paramsFindbeads3d)

    def convertModelToCoordinatesStep(self, tsObjId):
        tomo = self.inputSetOfTomograms.get()[tsObjId]
        location = tomo.getFileName()
        fileName, _ = os.path.splitext(location)

        extraPrefix = self._getExtraPath(os.path.basename(fileName))

        """ Run model2point IMOD program """
        paramsModel2Point = {
            'inputFile': os.path.join(extraPrefix, "%s.mod" % os.path.basename(fileName)),
            'outputFile': os.path.join(extraPrefix, "%s.xyz" % os.path.basename(fileName)),
        }

        argsModel2Point = "-InputFile %(inputFile)s " \
                          "-OutputFile %(outputFile)s "

        Plugin.runImod(self, 'model2point', argsModel2Point % paramsModel2Point)

    def createOutputStep(self, tsObjId):
        tomo = self.inputSetOfTomograms.get()[tsObjId]
        location = tomo.getFileName()
        fileName, _ = os.path.splitext(location)

        extraPrefix = self._getExtraPath(os.path.basename(fileName))

        """ Create the output set of coordinates 3D from gold beads detected """
        output = self.getOutputSetOfCoordinates3Ds(self.inputSetOfTomograms.get(),
                                                   self.inputSetOfTomograms.get())

        coordFilePath = os.path.join(extraPrefix,
                                     "%s.xyz" % os.path.basename(fileName))

        coordList = utils.formatGoldBead3DCoordinatesList(coordFilePath)

        with self._lock:
            for element in coordList:
                newCoord3D = tomoObj.Coordinate3D()
                newCoord3D.setVolume(tomo)
                newCoord3D.setX(element[0], constants.BOTTOM_LEFT_CORNER)
                newCoord3D.setY(element[1], constants.BOTTOM_LEFT_CORNER)
                newCoord3D.setZ(element[2], constants.BOTTOM_LEFT_CORNER)

                newCoord3D.setVolId(tsObjId)
                output.append(newCoord3D)
                output.update(newCoord3D)

            output.setBoxSize(self.beadDiameter.get())
            output.write()

        self._store()

    def closeOutputSetStep(self):
        self.Coordinates3D.setStreamState(Set.STREAM_CLOSED)
        self.Coordinates3D.write()
        self._store()
