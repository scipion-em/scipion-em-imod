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

    _label = 'tomogram reconstruction'

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
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('computeReconstructionStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            extraPrefix = self._getExtraPath(tsId)
            tmpPrefix = self._getTmpPath(tsId)
            path.makePath(tmpPrefix)
            path.makePath(extraPrefix)
            outputTsFileName = os.path.join(tmpPrefix, "%s.st" % tsId)

            """Apply the transformation form the input tilt-series"""
            ts.applyTransform(outputTsFileName)

            """Generate angle file"""
            angleFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)
            ts.generateTltFile(angleFilePath)

    def computeReconstructionStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            tomoPathOut = self._getExtraPath(os.path.join(tsId, "%s.mrc" % tsId))
            tsPathIn = ts.getFirstItem().getLocation()[1]
            angleFilePath = self._getTmpPath(os.path.join(tsId, "%s.rawtlt" % tsId))

            self.runJob('tilt', '-InputProjections %s -OutputFile %s -TILTFILE %s -THICKNESS %d' %
                        (tsPathIn, tomoPathOut, angleFilePath, self.tomoThickness.get()))

    def createOutputStep(self):
        self.outputSetOfTomograms = self._createSetOfTomograms()
        self.outputSetOfTomograms.setSamplingRate(self.inputSetOfTiltSeries.get().getSamplingRate())
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            newTomogram = Tomogram()
            newTomogram.setLocation(os.path.join(self._getExtraPath(tsId), '%s.mrc' % tsId))
            self.outputSetOfTomograms.append(newTomogram)
        self._defineOutputs(outputTomograms=self.outputSetOfTomograms)
        self._defineSourceRelation(self.inputSetOfTiltSeries, self.outputSetOfTomograms)

        path.moveTree(self._getTmpPath(), self._getExtraPath())

