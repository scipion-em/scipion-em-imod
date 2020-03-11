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
import numpy as np
import pyworkflow as pw
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
from pwem.emlib.image import ImageHandler
from tomo.protocols import ProtTomoBase
from tomo.convert import writeTiStack
from tomo.objects import Tomogram, TomoAcquisition


class ProtTomoReconstruction(EMProtocol, ProtTomoBase):
    """
    Tomogram reconstruction procedure based on the IMOD procedure.

    More info:
        https://bio3d.colorado.edu/imod/doc/etomoTutorial.html
    """

    _label = 'tomo reconstruction'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries', params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-Series')

        form.addParam('tomoThickness', params.FloatParam,
                      default=100,
                      label='Tomogram thickness', important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Size in pixels of the tomogram in the z axis (beam direction).')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('computeReconstructionStep', ts.getObjId())
            self._insertFunctionStep('createOutputStep', ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)
        outputTsFileName = os.path.join(tmpPrefix, "%s.st" % tsId)

        """Apply the transformation form the input tilt-series"""
        ts.applyTransform(outputTsFileName)
        # self.runJob('clip flipzy',
        #             ts.getFirstItem().getLocation()[1] + " " + os.path.join(tmpPrefix, "%s_flip.st" % tsId))
        # #path.createLink(ts.getFirstItem().getLocation()[1], os.path.join(tmpPrefix, "%s_flip.st" % tsId))
        #
        # newStack = True
        # for index, ti in enumerate(ts):
        #     if ti.hasTransform():
        #         ih = ImageHandler()
        #         if newStack:
        #             ih.createEmptyImage(fnOut=os.path.join(tmpPrefix, "%s_flip_transf.st" % tsId),
        #                                 xDim=ti.getXDim(),
        #                                 yDim=ti.getYDim(),
        #                                 nDim=ts.getSize())
        #             newStack = False
        #         transform = ti.getTransform().getMatrix()
        #         transformArray = np.array(transform)
        #         ih.applyTransform(inputFile=str(index + 1) + '@' + os.path.join(tmpPrefix, "%s_flip.st" % tsId),
        #                           outputFile=str(index + 1) + '@' + os.path.join(tmpPrefix,
        #                                                                          "%s_flip_transf.st" % tsId),
        #                           transformMatrix=transformArray,
        #                           shape=(ti.getXDim(), ti.getYDim()))
        # #self.runJob('clip flipyz',
        # #            os.path.join(tmpPrefix, "%s_flip_transf.st" % tsId) + " " + os.path.join(tmpPrefix,
        # #                                                                                     "%s_transf.st" % tsId))

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)
        ts.generateTltFile(angleFilePath)

    def computeReconstructionStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        paramsTilt = {
            'InputProjections': self._getTmpPath(os.path.join(tsId, "%s.st" % tsId)),
            'OutputFile': self._getExtraPath(os.path.join(tsId, "%s.rec" % tsId)),
            'TiltFile': self._getTmpPath(os.path.join(tsId, "%s.rawtlt" % tsId)),
            'Thickness': self.tomoThickness.get()
        }
        argsTilt = "-InputProjections %(InputProjections)s " \
                   "-OutputFile %(OutputFile)s " \
                   "-TILTFILE %(TiltFile)s " \
                   "-THICKNESS %(Thickness)d"
        self.runJob('tilt', argsTilt % paramsTilt)

        paramsNewstack = {
            'input': self._getExtraPath(os.path.join(tsId, "%s.rec" % tsId)),
            'output': self._getExtraPath(os.path.join(tsId, "%s.mrc" % tsId)),
        }
        argsNewstack = "-input %(input)s " \
                       "-output %(output)s"
        self.runJob('newstack', argsNewstack % paramsNewstack)

    def createOutputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        outputSetOfTomograms = self.getOutputSetOfTomograms()

        newTomogram = Tomogram()
        newTomogram.setLocation(os.path.join(extraPrefix, "%s.mrc" % tsId))
        outputSetOfTomograms.append(newTomogram)
        outputSetOfTomograms.update(newTomogram)
        outputSetOfTomograms.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfTomograms(self):
        if not hasattr(self, "outputSetOfTomograms"):
            outputSetOfTomograms = self._createSetOfTomograms()
            outputSetOfTomograms.copyInfo(self.inputSetOfTiltSeries.get())
            self._defineOutputs(outputSetOfTomograms=outputSetOfTomograms)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfTomograms)
        return self.outputSetOfTomograms

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputSetOfTomograms'):
            summary.append("Input Tilt-Series: %d.\nTomograms reconstructed: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputSetOfTomograms.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputSetOfTomograms'):
            methods.append("The reconstruction has been computed for %d "
                           "Tilt-series using the IMOD procedure.\n"
                           % (self.outputSetOfTomograms.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
