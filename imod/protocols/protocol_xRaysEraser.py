# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
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
from pwem.protocols import EMProtocol
import pyworkflow.utils.path as path
import pyworkflow.protocol.params as params
from tomo.protocols import ProtTomoBase
import imod.utils as utils
from imod import Plugin


class ProtImodXraysEraser(EMProtocol, ProtTomoBase):
    """
    Erase X-rays, defects, and fiducial markers from aligned tilt-series based on the IMOD procedure.
    More info:
            https://bio3d.colorado.edu/imod/doc/man/ccderaser.html
    """

    _label = 'x-rays eraser'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

        # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series.')

        form.addParam('inputSetOfLandmarkModels',
                      params.PointerParam,
                      pointerClass='SetOfLandmarkModels',
                      important=True,
                      label='Input set of tilt-series',
                      help='Input set of landmark models containing the location of the gold beads through the series')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('generateFiducialModel', ts.getObjId())
            self._insertFunctionStep('eraseXrays', ts.getObjId())

    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(extraPrefix)

        """Apply the transformation form the input tilt-series"""
        outputTsFileName = os.path.join(tmpPrefix, ts.getFirstItem().parseFileName())
        ts.applyTransform(outputTsFileName)

    def generateFiducialModel(self, tsObjId):
        # TODO: check si es el landmark model correcto
        lm = self.inputSetOfLandmarkModels.get()[tsObjId]
        ts = self.inputSetOfTiltSeries.get()[tsObjId]

        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        landmarkTextFilePath = os.path.join(extraPrefix,
                                            ts.getFirstItem().parseFileName(suffix="_fid", extension=".txt"))
        landmarkModelPath = os.path.join(extraPrefix,
                                         ts.getFirstItem().parseFileName(suffix="_fid", extension=".mod."))

        # Generate the IMOD file containing the information from the landmark model
        utils.generateIMODFiducialTextFile(landmarkModel=lm,
                                           outputFilePath=landmarkTextFilePath)

        # Convert IMOD file into IMOD model
        paramsPoint2Model = {
            'inputFile': landmarkTextFilePath,
            'outputFile': landmarkModelPath,
        }
        argsPoint2Model = "-InputFile %(inputFile)s " \
                          "-OutputFile %(outputFile)s"

        Plugin.runImod(self, 'point2model', argsPoint2Model % paramsPoint2Model)

    def eraseXrays(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]

        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsCcderaser = {
            'input': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName()),
            'output': os.path.join(extraPrefix, ts.getFirstItem().parseFileName()),
            'findPeaks': 1,
            'peakCriterion': 8.0, -
            'diffCriterion': 6.0, -
            'growCriterion': 4,
            'scanCriterion': 3,
            'maximumRadius': 4.2, -
            'giantCriterion': 12,
            'extraLargeRadius': 8,
            'bigDiffCriterion': 19, -
            'annulusWidth': 2.0,
            'xyScanSize': 100,
            'edgeExclusionWidth': 4,
            'pointModel': os.path.join(extraPrefix, ts.getFirstItem().parseFileName(suffix="_fid", extension=".mod.")),
            'borderSize': 2,
            'polynomialOrder': 2,
        }

        argsCcderaser = "-InputFile %()s " \
                        "-OutputFile %()s " \
                        "-FindPeaks %()d " \
                        "-PeakCriterion %()f " \
                        "-DiffCriterion %()f " \
                        "-GrowCriterion %()d " \
                        "-ScanCriterion %()d " \
                        "-MaximumRadius %()f " \
                        "-GiantCriterion %()d " \
                        "-ExtraLargeRadius %()d " \
                        "-BigDiffCriterion %()d " \
                        "-AnnulusWidth %()f " \
                        "-XYScanSize %()d " \
                        "-EdgeExclusionWidth %()d " \
                        "-PointModel %()s " \
                        "-BorderSize %()d " \
                        "-PolynomialOrder %()d "

        Plugin.runImod(self, 'ccderaser', argsCcderaser % paramsCcderaser)
