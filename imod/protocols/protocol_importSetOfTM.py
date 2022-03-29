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
import numpy as np
from pyworkflow import BETA
from pyworkflow.utils import path
import pyworkflow.protocol.params as params
import pwem.objects as data
from tomo.convert import getAnglesFromTlt
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries
from pwem.emlib.image import ImageHandler
from imod import utils
from imod.protocols.protocol_base import ProtImodBase


class ProtImodImportTransformationMatrix(ProtImodBase):
    """
    Import the transformation matrices assigned to an input set of tilt-series.
    """

    _label = 'Import transformation matrix'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        self._defineImportParams(form)

        form.addParam('exclusionWords', params.StringParam,
                      label='Exclusion words:',
                      help="List of words separated by a space that the path should not have",
                      expertLevel=params.LEVEL_ADVANCED)

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Set of tilt-series on which transformation matrices will be assigned.',
                      label='Assign transformation set of tilt-series')

        groupMatchBinning = form.addGroup('Match binning')

        groupMatchBinning.addParam('binningTM',
                                   params.IntParam,
                                   default=1,
                                   label='Transformation matrix binning',
                                   help='Binning of the tilt series at which the transformation matrices were '
                                        'calculated.')

        groupMatchBinning.addParam('binningTS',
                                   params.IntParam,
                                   default=1,
                                   label='Tilt-series binning',
                                   help='Binning of the tilt-series to which the ')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self.matchBinningFactor = self.binningTM.get() / self.binningTS.get()

        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep(self.generateTransformFileStep, ts.getObjId())
            self._insertFunctionStep(self.assignTransformationMatricesStep, ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def generateTransformFileStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        tsFileName = ts.getFirstItem().parseFileName(extension='')

        extraPrefix = self._getExtraPath(tsId)
        path.makePath(extraPrefix)

        outputTransformFile = os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".xf"))
        outputTltFile = os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".tlt"))

        for tmFilePath, _ in self.iterFiles():
            tmFileName = os.path.basename(os.path.splitext(tmFilePath)[0])

            if tsFileName == tmFileName:

                if self.matchBinningFactor != 1:

                    inputTransformMatrixList = utils.formatTransformationMatrix(tmFilePath)

                    extraPrefix = self._getExtraPath(tsId)
                    path.makePath(extraPrefix)

                    # Update shits from the transformation matrix considering the matching bining between the input tilt
                    # series and the transformation matrix. We create an empty tilt-series containing only tilt-images
                    # with transform information.
                    transformMatrixList = []

                    for index, ti in enumerate(ts):
                        inputTransformMatrix = inputTransformMatrixList[:, :, index]

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
                path.createLink(tmFilePath.replace('.xf', '.tlt'), outputTltFile)

    def assignTransformationMatricesStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        outputTransformFile = os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".xf"))
        outputTltFile = os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".tlt"))

        self.getOutputSetOfTiltSeries(self.inputSetOfTiltSeries.get())

        newTs = TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)

        self.outputSetOfTiltSeries.append(newTs)

        alignmentMatrix = utils.formatTransformationMatrix(outputTransformFile)
        angles = getAnglesFromTlt(outputTltFile)

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

            newTi.setTiltAngle(angles[index])

            newTs.append(newTi)

        ih = ImageHandler()
        x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
        newTs.setDim((x, y, z))

        newTs.write(properties=False)

        self.outputSetOfTiltSeries.update(newTs)
        self.outputSetOfTiltSeries.updateDim()
        self.outputSetOfTiltSeries.setStreamState(SetOfTiltSeries.STREAM_CLOSED)
        self.outputSetOfTiltSeries.write()

        self._store()

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        match = False

        for tmFilePath, _ in self.iterFiles():
            tmFileName = os.path.basename(os.path.splitext(tmFilePath)[0])

            for ts in self.inputSetOfTiltSeries.get():
                tsFileName = ts.getFirstItem().parseFileName(extension='')

                if tsFileName == tmFileName:
                    match = True

            if not match:
                validateMsgs.append("No matching tilt-series found for file: %s" % tmFilePath)

            match = False

        for ts in self.inputSetOfTiltSeries.get():
            tsFileName = ts.getFirstItem().parseFileName(extension='')

            for tmFilePath, _ in self.iterFiles():
                tmFileName = os.path.basename(os.path.splitext(tmFilePath)[0])

                if tsFileName == tmFileName:
                    match = True

            if not match:
                validateMsgs.append("No matching file found for tilt-series: %s (with tsID %s)"
                                    % (tsFileName, ts.getTsId()))

            match = False

        return validateMsgs

    def _summary(self):
        summary = []
        if hasattr(self, 'outputSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nTransformation matrices assigned: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputSetOfTiltSeries'):
            methods.append("The transformation matrix has been assigned to %d Tilt-series from the input set.\n"
                           % (self.outputSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
