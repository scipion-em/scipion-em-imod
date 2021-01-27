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
from pyworkflow.object import Set
from tomo.protocols import ProtTomoBase
import tomo.objects as tomoObj
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

        form.addParam('peakCriterion',
                      params.FloatParam,
                      default=8.0,
                      label='Peak criterion',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Criterion # of SDs above local mean for erasing peak based on intensity (the default is 10 '
                           'SDs)')

        form.addParam('diffCriterion',
                      params.FloatParam,
                      default=6.0,
                      label='Difference criterion',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Criterion # of SDs above mean pixel-to-pixel difference for erasing a peak based on '
                           'differences (the default is 10 SDs).')

        form.addParam('maximumRadius',
                      params.FloatParam,
                      default=4.2,
                      label='Maximum radius',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Maximum radius of peak area to erase (the default is 2.1 pixels).')

        form.addParam('bigDiffCriterion',
                      params.Integer,
                      default=19,
                      label='Big difference criterion',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='An extra-large peak will be erased only if the value for the maximum difference between '
                           'adjacent pixels, averaged over the most extreme one-fourth of the pixels in the patch, '
                           'exceeds this criterion, evaluated as the number of SDs above the mean absolute difference '
                           'between adjacent pixels in the scan area.  The default is 19.  This high a value is needed '
                           'to prevent gold erasure on low-noise data sets with small gold particles, and a lower '
                           'value may be needed to make extra-large peak removal useful.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('generateFiducialModelStep', ts.getObjId())
            self._insertFunctionStep('eraseXraysStep', ts.getObjId())
            self._insertFunctionStep('createOutputStep', ts.getObjId())

    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(extraPrefix)

        """Apply the transformation form the input tilt-series"""
        outputTsFileName = os.path.join(tmpPrefix, ts.getFirstItem().parseFileName())
        ts.applyTransform(outputTsFileName)

    def generateFiducialModelStep(self, tsObjId):
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

    def eraseXraysStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]

        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsCcderaser = {
            'input': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName()),
            'output': os.path.join(extraPrefix, ts.getFirstItem().parseFileName()),
            'findPeaks': 1,
            'peakCriterion': self.peakCriterion.get(),
            'diffCriterion': self.diffCriterion.get(),
            'growCriterion': 4,
            'scanCriterion': 3,
            'maximumRadius': self.maximumRadius.get(),
            'giantCriterion': 12,
            'extraLargeRadius': 8,
            'bigDiffCriterion': self.bigDiffCriterion.get(),
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

    def createOutputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]

        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        outputSetOfTiltSeries = self.getOutputSetOfTiltSeries()

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)

        outputSetOfTiltSeries.append(newTs)

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(
                index + 1,
                os.path.join(extraPrefix, ts.getFirstItem().parseFileName())
            )
            newTs.append(newTi)

        newTs.write(properties=False)

        outputSetOfTiltSeries.update(newTs)
        outputSetOfTiltSeries.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfTiltSeries(self):
        if hasattr(self, "outputSetOfTiltSeries"):
            self.outputSetOfTiltSeries.enableAppend()
        else:
            outputSetOfTiltSeries = self._createSetOfTiltSeries()
            outputSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfTiltSeries=outputSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfTiltSeries)
        return self.outputSetOfTiltSeries