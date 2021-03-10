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

from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pyworkflow.utils import path
from pyworkflow.object import Set
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
import tomo.objects as tomoObj
import tomo.constants as constants
from imod import Plugin
from imod import utils

import os


class ProtImodGoldBeadPicker3d(EMProtocol, ProtTomoBase):
    """
    3-dimensional gold bead picker using the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/findbeads3d.html
    """

    _label = 'Gold bead picker 3D'
    _devStatus = BETA

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        ProtTomoBase.__init__(self)
        self.stepsExecutionMode = STEPS_PARALLEL

# -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTomograms',
                      params.PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Input set of tilt-series',
                      help='Input set of tomograms from which gold beads will be picked.')

        form.addParam('beadDiameter',
                      params.IntParam,
                      label='Fiducial diameter (pixels)',
                      default='10',
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
                      label='Minimum relative Strength',
                      default=0.05,
                      help='Minimum relative peak strength for keeping a peak in the analysis.  The square root of the '
                           'specified value is used for comparing with the square root of peak strength, for '
                           'compatibility with existing command files. The default is 0.05, which corresponds to a '
                           'relative square root peak strength of 0.22.  Too many weak peaks can prevent a dip from '
                           'showing up in the smoothed histogram of strengths.  If the program fails to find a '
                           'histogram dip, one strategy is to try raising this value.')

        form.addParam('minSpacing',
                      params.FloatParam,
                      label='Minimum spacing',
                      default=0.9,
                      help='Minimum spacing between peaks as a fraction of the bead size. When two peaks are closer '
                           'than this distance apart, the weaker one is eliminated unless the -both option is entered. '
                           'The default is 0.9.  A value less than 1 is helpful for picking both beads in a pair.')

        form.addParallelSection(threads=4, mpi=1)

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        allOutputId = []
        for ts in self.inputSetOfTomograms.get():
            pickId = self._insertFunctionStep('pickGoldBeadsStep',
                                              ts.getObjId(),
                                              prerequisites=[])

            convertId = self._insertFunctionStep('convertModelToCoordinatesStep',
                                                 ts.getObjId(),
                                                 prerequisites=[pickId])

            outputID = self._insertFunctionStep('createOutputStep',
                                     ts.getObjId(),
                                     prerequisites=[convertId])

            allOutputId.append(outputID)

        self._insertFunctionStep('closeOutputSetStep',
                                 prerequisites=allOutputId)

    # --------------------------- STEPS functions ----------------------------
    def pickGoldBeadsStep(self, tsObjId):
        tomo = self.inputSetOfTomograms.get()[tsObjId]
        location = tomo.getLocation()[1]
        fileName, _ = os.path.splitext(location)

        extraPrefix = self._getExtraPath(os.path.basename(fileName))
        path.makePath(extraPrefix)

        """ Run findbeads3d IMOD program """
        paramsFindbeads3d = {
            'inputFile': location,
            'outputFile': os.path.join(extraPrefix, "%s.mod" % os.path.basename(fileName)),
            'beadSize': self.beadDiameter.get(),
            'minRelativeStrength': self.minRelativeStrength.get(),
            'minSpacing': self.minSpacing.get(),
        }

        argsFindbeads3d = "-InputFile %(inputFile)s " \
                          "-OutputFile %(outputFile)s " \
                          "-BeadSize %(beadSize)d " \
                          "-MinRelativeStrength %(minRelativeStrength)f " \
                          "-MinSpacing %(minSpacing)d "

        if self.beadsColor.get() == 1:
            argsFindbeads3d += "-LightBeads "

        Plugin.runImod(self, 'findbeads3d', argsFindbeads3d % paramsFindbeads3d)

    def convertModelToCoordinatesStep(self, tsObjId):
        tomo = self.inputSetOfTomograms.get()[tsObjId]
        location = tomo.getLocation()[1]
        fileName, _ = os.path.splitext(location)

        extraPrefix = self._getExtraPath(os.path.basename(fileName))

        """ Run model2point IMOD program """
        paramsModel2Point = {
            'inputFile': os.path.join(extraPrefix, "%s.mod" % os.path.basename(fileName)),
            'outputFile': os.path.join(extraPrefix, "%s.xyz" % os.path.basename(fileName)),
        }

        argsModel2Point = "-InputFile %(inputFile)s " \
                          "-OutputFile %(outputFile)s " \

        Plugin.runImod(self, 'model2point', argsModel2Point % paramsModel2Point)

    def createOutputStep(self, tsObjId):
        tomo = self.inputSetOfTomograms.get()[tsObjId]
        location = tomo.getLocation()[1]
        fileName, _ = os.path.splitext(location)

        extraPrefix = self._getExtraPath(os.path.basename(fileName))

        """ Create the output set of coordinates 3D from gold beads detected """
        outputSetOfCoordinates3D = self.getOutputSetOfCoordinates3Ds()

        coordFilePath = os.path.join(extraPrefix, "%s.xyz" % os.path.basename(fileName))

        coordList = utils.formatGoldBead3DCoordinatesList(coordFilePath)

        for element in coordList:
            newCoord3D = tomoObj.Coordinate3D()
            newCoord3D.setVolume(tomo)
            newCoord3D.setX(element[0], constants.BOTTOM_LEFT_CORNER)
            newCoord3D.setY(element[1], constants.BOTTOM_LEFT_CORNER)
            newCoord3D.setZ(element[2], constants.BOTTOM_LEFT_CORNER)

            newCoord3D.setVolId(tsObjId)
            outputSetOfCoordinates3D.append(newCoord3D)
            outputSetOfCoordinates3D.update(newCoord3D)

        outputSetOfCoordinates3D.write()

        self._store()

    def closeOutputSetStep(self):
        self.getOutputSetOfCoordinates3Ds().setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfCoordinates3Ds(self):
        if hasattr(self, "outputSetOfCoordinates3D"):
            self.outputSetOfCoordinates3D.enableAppend()
        else:
            outputSetOfCoordinates3D = self._createSetOfCoordinates3D(volSet=self.inputSetOfTomograms.get(),
                                                                      suffix='LandmarkModel')
            outputSetOfCoordinates3D.setSamplingRate(self.inputSetOfTomograms.get().getSamplingRate())
            outputSetOfCoordinates3D.setPrecedents(self.inputSetOfTomograms)
            outputSetOfCoordinates3D.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfCoordinates3D=outputSetOfCoordinates3D)
            self._defineSourceRelation(self.inputSetOfTomograms, outputSetOfCoordinates3D)
        return self.outputSetOfCoordinates3D

