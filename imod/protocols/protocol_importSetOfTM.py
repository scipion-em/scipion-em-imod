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
import pyworkflow.protocol.params as params
import pwem.objects as data
import tomo.objects as tomoObj
from pwem.emlib.image import ImageHandler
from imod import utils
from imod.protocols.protocol_base import ProtImodBase


class ProtImodImportTransformationMatrix(ProtImodBase):
    """
    Import the transformation matrices assigned to an input set of tilt-series.
    """

    _label = 'import transformation matrix'
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

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('assignTransformationMatricesStep')

    # --------------------------- STEPS functions ----------------------------
    def assignTransformationMatricesStep(self):
        self.getOutputSetOfTiltSeries()

        inputSetOfTiltSeries = self.inputSetOfTiltSeries.get()

        for ts in inputSetOfTiltSeries:
            tsId = ts.getTsId()

            tsFileName = ts.getFirstItem().parseFileName(extension='')

            for tmFilePath, _ in self.iterFiles():
                tmFileName = os.path.basename(os.path.splitext(tmFilePath)[0])

                if tsFileName == tmFileName:
                    alignmentMatrix = utils.formatTransformationMatrix(tmFilePath)

                    newTs = tomoObj.TiltSeries(tsId=tsId)
                    newTs.copyInfo(ts)
                    self.outputSetOfTiltSeries.append(newTs)

                    for index, tiltImage in enumerate(ts):
                        newTi = tomoObj.TiltImage()
                        newTi.copyInfo(tiltImage, copyId=True)
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

                    self.outputSetOfTiltSeries.update(newTs)
                    self.outputSetOfTiltSeries.updateDim()
                    self.outputSetOfTiltSeries.write()

                    self._store()

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        match = False

        for tmFilePath, _ in self.iterFiles():
            tmFileName = os.path.basename(os.path.splitext(tmFilePath)[0])
            print(tmFilePath)

            for ts in self.inputSetOfTiltSeries.get():
                tsFileName = ts.getFirstItem().parseFileName(extension='')

                if tsFileName == tmFileName:
                    match = True

            if not match:
                validateMsgs.append("No matching tilt-series found for file: %s" % tmFilePath)

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
