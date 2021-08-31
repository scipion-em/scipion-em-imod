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
from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pyworkflow.object import Set
from tomo.objects import TiltSeries, TiltImage
import imod.utils as utils
from imod import Plugin
from imod.protocols.protocol_base import ProtImodBase


class ProtImodGoldBeadEraser(ProtImodBase):
    """
    Erase fiducial markers from aligned tilt-series based on the IMOD procedure.
    More info:
            https://bio3d.colorado.edu/imod/doc/man/ccderaser.html
    """

    _label = 'Gold bead eraser'
    _devStatus = BETA

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
                      label='Input set of landmark models (no gaps)',
                      help='Input set of landmark models containing the location of the gold beads to be erased '
                           'through the series.\n'
                           'IMPORTANT: It is highly recommended to use a Landmark Model with no gaps.')

        form.addParam('useTMFromTS',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=0,
                      label='Use alignment from TS',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='If this option is set to yes, the alignment use to remove the gold beads from the fiducial '
                           'model is the one associated to the input set of tilt-series. In case other alignment is to'
                           'be used, input a tilt-series.\n'
                           'IMPORTANT: the alignment that should be applied to the model is the one correcting for the '
                           'tilt axis rotation (no prealignment). This corresponds to the _fid.xf in the fiducial '
                           'alignment or etomo protocol.')

        form.addParam('inputSetOfTiltSeriesTransform',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      condition='useTMFromTS==1',
                      important=True,
                      label='Input set of tilt-series transform.',
                      help='Input set of tilt-series from which the alignment will be used to remove the gold beads '
                           'from the fiducial model')

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
            tsObjId = ts.getObjId()
            self._insertFunctionStep(self.convertInputStep, tsObjId, False, False)
            self._insertFunctionStep(self.generateIntermediateFiles, tsObjId)
            self._insertFunctionStep(self.generateFiducialModelStep, tsObjId)
            self._insertFunctionStep(self.eraseGoldBeadStep, tsObjId)
            self._insertFunctionStep(self.createOutputStep, tsObjId)
        self._insertFunctionStep(self.closeOutputStep)

    def generateIntermediateFiles(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        firstItem = ts.getFirstItem()

        extraPrefix = self._getExtraPath(tsId)

        outputTMPath = os.path.join(extraPrefix, firstItem.parseFileName(extension=".xf"))

        if self.useTMFromTS.get() == 0:
            utils.formatTransformFile(ts, outputTMPath)

        else:
            tsTM = self.getTiltSeriesFromTs(self.inputSetOfTiltSeriesTransform.get(), tsId)

            utils.formatTransformFile(tsTM, outputTMPath)

    def generateFiducialModelStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        firstItem = ts.getFirstItem()

        lm = self.getLandMarkModelFromTs(self.inputSetOfLandmarkModels.get(), tsId)

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        outputTMPath = os.path.join(extraPrefix, firstItem.parseFileName(extension=".xf"))

        # Generate model if it does not exist
        # if lm.getModelName() is not None:
        #     landmarkModelPath = lm.getModelName()
        #
        # else:

        if True:
            landmarkTextFilePath = lm.getFileName()
            landmarkModelPath = os.path.join(extraPrefix,
                                             firstItem.parseFileName(extension=".fid")
                                             )

            # Generate the IMOD file containing the information from the landmark model
            utils.generateIMODFiducialTextFile(landmarkModel=lm,
                                               outputFilePath=landmarkTextFilePath)

            # Convert IMOD file into IMOD model
            paramsPoint2Model = {
                'inputFile': landmarkTextFilePath,
                'outputFile': landmarkModelPath,
                'image': os.path.join(tmpPrefix, firstItem.parseFileName())
            }

            argsPoint2Model = "-InputFile %(inputFile)s " \
                              "-OutputFile %(outputFile)s " \
                              "-ImageForCoordinates %(image)s "

            Plugin.runImod(self, 'point2model', argsPoint2Model % paramsPoint2Model)

        # Generate interpolated model
        paramsImodtrans = {
            'transformFile': outputTMPath,
            'image': os.path.join(tmpPrefix, firstItem.parseFileName()),
            'inputFile': landmarkModelPath,
            'outputFile': os.path.join(extraPrefix, firstItem.parseFileName(suffix="_fid", extension=".mod"))
        }

        argsImodtrans = "-2 %(transformFile)s " \
                        "%(inputFile)s " \
                        "%(outputFile)s "

#                        "-i %(image)s " \

        Plugin.runImod(self, 'imodtrans', argsImodtrans % paramsImodtrans)

    def eraseGoldBeadStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        firstItem = ts.getFirstItem()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsCcderaser = {
            'inputFile': os.path.join(tmpPrefix, firstItem.parseFileName()),
            'outputFile': os.path.join(extraPrefix, firstItem.parseFileName()),
            'modelFile': os.path.join(extraPrefix, firstItem.parseFileName(suffix="_fid", extension=".mod")),
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
        self.getOutputSetOfTiltSeries(self.inputSetOfTiltSeries.get())

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        newTs = TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        self.outputSetOfTiltSeries.append(newTs)

        for index, tiltImage in enumerate(ts):
            newTi = TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(index + 1,
                              (os.path.join(extraPrefix, tiltImage.parseFileName())))
            newTs.append(newTi)

        newTs.write(properties=False)
        self.outputSetOfTiltSeries.update(newTs)
        self.outputSetOfTiltSeries.write()
        self._store()

    def closeOutputStep(self):
        self.outputSetOfTiltSeries.setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        tsNoTransformSet = []
        transform = True

        for ts in self.inputSetOfTiltSeries.get():
            if not ts.getFirstItem().hasTransform():
                transform = False

                tsNoTransformSet.append(ts)

        if not transform and not self.useTMFromTS.get():
            validateMsgs.append("The following tilt-series (tsId) does not have an associated alignment and there is "
                                "no optional source from where to input alignment information: ")

            for tsNoTransform in tsNoTransformSet:
                validateMsgs.append(tsNoTransform.getTsId())

        return validateMsgs
