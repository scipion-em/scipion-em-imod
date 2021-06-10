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
from pwem.protocols import EMProtocol
from pyworkflow import BETA
import pyworkflow.utils.path as path
import pyworkflow.protocol.params as params
from pyworkflow.object import Set
from tomo.protocols import ProtTomoBase
import tomo.objects as tomoObj
import imod.utils as utils
from imod import Plugin


class ProtImodGoldBeadEraser(EMProtocol, ProtTomoBase):
    """
    Erase fiducial markers from aligned tilt-series based on the IMOD procedure.
    More info:
            https://bio3d.colorado.edu/imod/doc/man/ccderaser.html
    """

    _label = 'gold bead eraser'
    _devStatus = BETA

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
                      label='Input set of landmark models',
                      help='Input set of landmark models containing the location of the gold beads through the series')

        form.addParam('betterRadius',
                      params.IntParam,
                      default=10,
                      label='Bead diameter (pixels)',
                      help="For circle objects, this entry specifies a radius to use for points without an individual "
                           "point size instead of the object's default sphere radius.  This entry is floating point "
                           "and can be used to overcome the limitations of having an integer default sphere radius. If "
                           "there are multiple circle objects, enter one value to apply to all objects or a value for "
                           "each object.")

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('generateFiducialModelStep', ts.getObjId())
            self._insertFunctionStep('eraseXraysStep', ts.getObjId())
            self._insertFunctionStep('createOutputStep', ts.getObjId())
        self._insertFunctionStep('closeOutputStep')

    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
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
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(extraPrefix)
        path.makePath(tmpPrefix)

        landmarkTextFilePath = os.path.join(extraPrefix,
                                            ts.getFirstItem().parseFileName(suffix="_fid", extension=".txt"))
        landmarkModelPath = os.path.join(extraPrefix,
                                         ts.getFirstItem().parseFileName(suffix="_fid", extension=".mod"))

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
            'inputFile': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName()),
            'outputFile': os.path.join(extraPrefix, ts.getFirstItem().parseFileName()),
            'modelFile': os.path.join(extraPrefix, ts.getFirstItem().parseFileName(suffix="_fid", extension=".mod")),
            'betterRadius': self.betterRadius.get(),
            'polynomialOrder': 0,
            'circleObjects': "/"
        }

        argsCcderaser = "-InputFile %(inputFile)s " \
                        "-OutputFile %(outputFile)s " \
                        "-ModelFile %(modelFile)s " \
                        "-BetterRadius %(betterRadius)d " \
                        "-PolynomialOrder %(polynomialOrder)d " \
                        "-CircleObjects %(circleObjects)s " \
                        "-MergePatches " \
                        "-ExcludeAdjacent"

        Plugin.runImod(self, 'ccderaser', argsCcderaser % paramsCcderaser)

    def createOutputStep(self, tsObjId):
        outputXraysErasedSetOfTiltSeries = self.getOutputSetOfXraysErasedTiltSeries()

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputXraysErasedSetOfTiltSeries.append(newTs)

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(index + 1,
                              (os.path.join(extraPrefix, tiltImage.parseFileName())))
            newTs.append(newTi)

        newTs.write(properties=False)
        outputXraysErasedSetOfTiltSeries.update(newTs)
        outputXraysErasedSetOfTiltSeries.write()
        self._store()

    def closeOutputStep(self):
        self.getOutputSetOfXraysErasedTiltSeries().setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfXraysErasedTiltSeries(self):
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
