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
import numpy as np

from pyworkflow import BETA
from pyworkflow.object import Set
from pyworkflow.utils import path
import pyworkflow.protocol.params as params
import pwem.objects as data
from pwem.emlib.image import ImageHandler
from tomo.objects import TiltSeries, TiltImage

from .. import utils
from .protocol_base import ProtImodBase, OUTPUT_TILTSERIES_NAME


class ProtImodImportTransformationMatrix(ProtImodBase):
    """
    Import the transformation matrices assigned to an input set of tilt-series
    """

    _label = 'Import transformation matrix'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        self._defineImportParams(form)

        form.addParam('exclusionWords', params.StringParam,
                      label='Exclusion words:',
                      help="List of words separated by a space that the "
                           "path should not have",
                      expertLevel=params.LEVEL_ADVANCED)

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Set of tilt-series on which transformation matrices '
                           'will be assigned.',
                      label='Input set of tilt-series')

        groupMatchBinning = form.addGroup('Match binning')

        groupMatchBinning.addParam('binningTM',
                                   params.IntParam,
                                   default=1,
                                   label='Transformation matrix binning',
                                   help='Binning of the tilt series at which '
                                        'the transformation matrices were '
                                        'calculated.')

        groupMatchBinning.addParam('binningTS',
                                   params.IntParam,
                                   default=1,
                                   label='Tilt-series binning',
                                   help='Binning of the tilt-series.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self.matchBinningFactor = self.binningTM.get() / self.binningTS.get()

        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep(self.generateTransformFileStep,
                                     ts.getObjId())
            self._insertFunctionStep(self.assignTransformationMatricesStep,
                                     ts.getObjId())

        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def generateTransformFileStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        tsFileName = ts.getFirstItem().parseFileName(extension='')

        extraPrefix = self._getExtraPath(tsId)
        path.makePath(extraPrefix)

        outputTransformFile = os.path.join(extraPrefix,
                                           ts.getFirstItem().parseFileName(extension=".xf"))

        for tmFilePath, _ in self.iterFiles():
            tmFileName = os.path.basename(os.path.splitext(tmFilePath)[0])

            if tsFileName == tmFileName:

                if self.matchBinningFactor != 1:

                    inputTransformMatrixList = utils.formatTransformationMatrix(tmFilePath)

                    extraPrefix = self._getExtraPath(tsId)
                    path.makePath(extraPrefix)

                    # Update shifts from the transformation matrix considering
                    # the matching binning between the input tilt
                    # series and the transformation matrix. We create an empty
                    # tilt-series containing only tilt-images
                    # with transform information.
                    transformMatrixList = []

                    ids = ts.getIdSet()
                    for index in ids:
                        inputTransformMatrix = inputTransformMatrixList[:, :, index-1]

                        outputTransformMatrix = inputTransformMatrix
                        outputTransformMatrix[0][0] = inputTransformMatrix[0][0]
                        outputTransformMatrix[0][1] = inputTransformMatrix[0][1]
                        outputTransformMatrix[0][2] = inputTransformMatrix[0][2] * self.matchBinningFactor
                        outputTransformMatrix[1][0] = inputTransformMatrix[1][0]
                        outputTransformMatrix[1][1] = inputTransformMatrix[1][1]
                        outputTransformMatrix[1][2] = inputTransformMatrix[1][2] * self.matchBinningFactor
                        outputTransformMatrix[2][0] = inputTransformMatrix[2][0]
                        outputTransformMatrix[2][1] = inputTransformMatrix[2][1]
                        outputTransformMatrix[2][2] = inputTransformMatrix[2][2]

                        transformMatrixList.append(outputTransformMatrix)

                    utils.formatTransformFileFromTransformList(transformMatrixList, outputTransformFile)

                else:
                    path.createLink(tmFilePath, outputTransformFile)

    def assignTransformationMatricesStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        outputTransformFile = os.path.join(extraPrefix,
                                           ts.getFirstItem().parseFileName(extension=".xf"))

        output = self.getOutputSetOfTiltSeries(self.inputSetOfTiltSeries.get())

        newTs = TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)

        output.append(newTs)

        alignmentMatrix = utils.formatTransformationMatrix(outputTransformFile)

        for index, tiltImage in enumerate(ts):
            newTi = TiltImage()
            newTi.copyInfo(tiltImage, copyId=True, copyTM=False)
            newTi.setAcquisition(tiltImage.getAcquisition())
            newTi.setLocation(tiltImage.getLocation())

            transform = data.Transform()

            if tiltImage.hasTransform():
                previousTransform = tiltImage.getTransform().getMatrix()
                newTransform = alignmentMatrix[:, :, index]
                previousTransformArray = np.array(previousTransform)
                newTransformArray = np.array(newTransform)
                outputTransformMatrix = np.matmul(previousTransformArray, newTransformArray)
                transform.setMatrix(outputTransformMatrix)
                newTi.setTransform(transform)

            else:
                transform.setMatrix(alignmentMatrix[:, :, index])
                newTi.setTransform(transform)

            newTs.append(newTi)

        ih = ImageHandler()
        x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
        newTs.setDim((x, y, z))

        newTs.write(properties=False)

        output.update(newTs)
        output.write()

        self._store()

    def closeOutputSetsStep(self):

        output = getattr(self, OUTPUT_TILTSERIES_NAME)
        output.setStreamState(Set.STREAM_CLOSED)
        output.write()

        self._store()

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        match = False

        for ts in self.inputSetOfTiltSeries.get():
            tsFileName = ts.getFirstItem().parseFileName(extension='')

            for tmFilePath, _ in self.iterFiles():
                tmFileName = os.path.basename(os.path.splitext(tmFilePath)[0])

                if tsFileName == tmFileName:
                    match = True
                    break

            if not match:
                validateMsgs.append("No xf file found for tilt-series %s: image file is %s and have not found its "
                                    "exact match."  % (ts.getTsId(), tsFileName, ))

            match = False

        return validateMsgs

    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append("Input tilt-series: %d\nTransformation matrices "
                           "assigned: %d"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.TiltSeries.getSize()))
        return summary
