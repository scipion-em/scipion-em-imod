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
from pwem.protocols import EMProtocol
import pwem.objects as data
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
from tomo.protocols.protocol_base import ProtTomoImportFiles
from pwem.emlib.image import ImageHandler
from imod import utils


class ProtImodImportTransformationMatrix(ProtTomoImportFiles, EMProtocol, ProtTomoBase):
    """
    Import the transformation matrices assigned to an input set of tilt-series.
    """

    _label = 'Import transformation matrix'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def __init__(self, **args):
        ProtTomoImportFiles.__init__(self, **args)

    def _defineParams(self, form):
        ProtTomoImportFiles._defineImportParams(self, form)

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
        self.getOutputAssignedTransformSetOfTiltSeries()

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
                    self.outputAssignedTransformSetOfTiltSeries.append(newTs)

                    for index, tiltImage in enumerate(ts):
                        newTi = tomoObj.TiltImage()
                        newTi.copyInfo(tiltImage, copyId=True)
                        newTi.setLocation(tiltImage.getLocation())
                        transform = data.Transform()
                        transform.setMatrix(alignmentMatrix[:, :, index])
                        newTi.setTransform(transform)
                        newTs.append(newTi)

                    ih = ImageHandler()
                    x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
                    newTs.setDim((x, y, z))

                    newTs.write(properties=False)

                    self.outputAssignedTransformSetOfTiltSeries.update(newTs)
                    self.outputAssignedTransformSetOfTiltSeries.updateDim()
                    self.outputAssignedTransformSetOfTiltSeries.write()

                    self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputAssignedTransformSetOfTiltSeries(self):
        if not hasattr(self, "outputAssignedTransformSetOfTiltSeries"):
            outputAssignedTransformSetOfTiltSeries = self._createSetOfTiltSeries(suffix='AssignedTransform')
            outputAssignedTransformSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputAssignedTransformSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())

            self._defineOutputs(outputAssignedTransformSetOfTiltSeries=outputAssignedTransformSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputAssignedTransformSetOfTiltSeries)
        return self.outputAssignedTransformSetOfTiltSeries

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
        if hasattr(self, 'outputAssignedTransformSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nTransformation matrices assigned: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputAssignedTransformSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputAssignedTransformSetOfTiltSeries'):
            methods.append("The transformation matrix has been assigned to %d Tilt-series from the input set.\n"
                           % (self.outputAssignedTransformSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
