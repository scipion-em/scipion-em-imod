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
import imod.utils as utils
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pyworkflow.object import Set
from tomo.objects import TiltSeries, TiltImage
from imod import Plugin
from pwem.emlib.image import ImageHandler
from imod.protocols.protocol_base import ProtImodBase


class ProtImodApplyTransformationMatrix(ProtImodBase):
    """
    Compute the interpolated tilt-series from its transform matrix.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/newstack.html
    """

    _label = 'Apply transformation'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        form.addParam('binning', params.FloatParam,
                      default=1.0,
                      label='Binning',
                      help='Binning to be applied to the interpolated tilt-series in IMOD convention. Images will be '
                           'binned by the given factor. Must be an integer bigger than 1')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep(self.generateTransformFileStep, ts.getObjId())
            self._insertFunctionStep(self.computeAlignmentStep, ts.getObjId())
            self._insertFunctionStep(self.generateOutputStackStep, ts.getObjId())
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions ----------------------------
    def generateTransformFileStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        path.makePath(extraPrefix)
        utils.formatTransformFile(ts, os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension="_fid.xf")))

    def computeAlignmentStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        firstItem = ts.getFirstItem()

        paramsAlignment = {
            'input': firstItem.getFileName(),
            'output': os.path.join(extraPrefix, firstItem.parseFileName()),
            'xform': os.path.join(extraPrefix, firstItem.parseFileName(extension="_fid.xf")),
            'bin': int(self.binning.get()),
            'imagebinned': 1.0
        }

        argsAlignment = "-input %(input)s " \
                        "-output %(output)s " \
                        "-xform %(xform)s " \
                        "-bin %(bin)d " \
                        "-imagebinned %(imagebinned)s "

        rotationAngleAvg = utils.calculateRotationAngleFromTM(ts)

        # Check if rotation angle is greater than 45º. If so, swap x and y dimensions to adapt output image sizes to
        # the final sample disposition.
        if rotationAngleAvg > 45 or rotationAngleAvg < -45:
            paramsAlignment.update({
                'size': "%d,%d" % (firstItem.getYDim(), firstItem.getXDim())
            })

            argsAlignment += "-size %(size)s "

        Plugin.runImod(self, 'newstack', argsAlignment % paramsAlignment)

    def generateOutputStackStep(self, tsObjId):
        self.getOutputInterpolatedSetOfTiltSeries(self.inputSetOfTiltSeries.get())

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        newTs = TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        self.outputInterpolatedSetOfTiltSeries.append(newTs)

        if self.binning > 1:
            newTs.setSamplingRate(ts.getSamplingRate() * int(self.binning.get()))

        for index, tiltImage in enumerate(ts):
            newTi = TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setAcquisition(tiltImage.getAcquisition())
            newTi.setLocation(index + 1, (os.path.join(extraPrefix, tiltImage.parseFileName())))
            if self.binning > 1:
                newTi.setSamplingRate(tiltImage.getSamplingRate() * int(self.binning.get()))
            newTs.append(newTi)

        ih = ImageHandler()
        x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
        newTs.setDim((x, y, z))

        newTs.write(properties=False)

        self.outputInterpolatedSetOfTiltSeries.update(newTs)
        self.outputInterpolatedSetOfTiltSeries.updateDim()
        self.outputInterpolatedSetOfTiltSeries.write()
        self._store()

    def closeOutputSetsStep(self):
        self.outputInterpolatedSetOfTiltSeries.setStreamState(Set.STREAM_CLOSED)
        self.outputInterpolatedSetOfTiltSeries.write()
        self._store()

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        for ts in self.inputSetOfTiltSeries.get():
            if not ts.getFirstItem().hasTransform():
                validateMsgs.append("Some tilt-series from the input set of tilt-series is missing from a "
                                    "transformation matrix.")
                break

        return validateMsgs

    def _summary(self):
        summary = []
        if hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nInterpolations applied: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            methods.append("The interpolation has been computed for %d "
                           "Tilt-series using the IMOD newstack program.\n"
                           % (self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods