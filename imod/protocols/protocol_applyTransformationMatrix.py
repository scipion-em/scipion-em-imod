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

from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pyworkflow.object import Set
from pwem.emlib.image import ImageHandler
from tomo.objects import TiltSeries, TiltImage

from .. import Plugin, utils
from .protocol_base import ProtImodBase


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
                      help='Binning to be applied to the interpolated '
                           'tilt-series in IMOD convention. Images will be '
                           'binned by the given factor. Must be an integer '
                           'bigger than 1')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep(self.generateTransformFileStep, ts.getObjId())
            self._insertFunctionStep(self.computeAlignmentStep, ts.getObjId())
            self._insertFunctionStep(self.generateOutputStackStep, ts.getObjId())
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions ------------------------------
    def generateTransformFileStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        path.makePath(extraPrefix)
        utils.formatTransformFile(ts,
                                  os.path.join(extraPrefix,
                                               ts.getFirstItem().parseFileName(extension=".xf")))

    def computeAlignmentStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        firstItem = ts.getFirstItem()
        binning = int(self.binning.get())

        paramsAlignment = {
            'input': firstItem.getFileName(),
            'output': os.path.join(extraPrefix, firstItem.parseFileName()),
            'xform': os.path.join(extraPrefix, firstItem.parseFileName(extension=".xf")),
            'bin': binning,
            'imagebinned': 1.0
        }

        argsAlignment = "-input %(input)s " \
                        "-output %(output)s " \
                        "-xform %(xform)s " \
                        "-bin %(bin)d " \
                        "-antialias -1 " \
                        "-imagebinned %(imagebinned)s "

        rotationAngle = ts.getAcquisition().getTiltAxisAngle()

        # Check if rotation angle is greater than 45º. If so,
        # swap x and y dimensions to adapt output image sizes to
        # the final sample disposition.
        if 45 < abs(rotationAngle) < 135:
            paramsAlignment.update({
                'size': "%d,%d" % (round(firstItem.getYDim() / binning),
                                   round(firstItem.getXDim() / binning))
            })

            argsAlignment += "-size %(size)s "

        excludedViews = ts.getExcludedViewsIndex(caster=str, indexOffset=-1)
        if len(excludedViews):
            argsAlignment += "-exclude %s " % ",".join(excludedViews)

        Plugin.runImod(self, 'newstack', argsAlignment % paramsAlignment)

    def generateOutputStackStep(self, tsObjId):
        output = self.getOutputInterpolatedSetOfTiltSeries(self.inputSetOfTiltSeries.get())

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        binning = int(self.binning.get())

        newTs = TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        newTs.setInterpolated(True)
        acq = newTs.getAcquisition()
        acq.setTiltAxisAngle(0.)  # 0 because TS is aligned
        newTs.setAcquisition(acq)
        output.append(newTs)

        if binning > 1:
            newTs.setSamplingRate(ts.getSamplingRate() * binning)

        index = 1
        for tiltImage in ts:
            if tiltImage.isEnabled():
                newTi = TiltImage()
                newTi.copyInfo(tiltImage, copyId=False, copyTM=False)
                acq = tiltImage.getAcquisition()
                acq.setTiltAxisAngle(0.)
                newTi.setAcquisition(acq)
                newTi.setLocation(index, (os.path.join(extraPrefix, tiltImage.parseFileName())))
                index += 1
                if binning > 1:
                    newTi.setSamplingRate(tiltImage.getSamplingRate() * binning)
                newTs.append(newTi)

        ih = ImageHandler()
        x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
        newTs.setDim((x, y, z))

        newTs.write(properties=False)
        output.update(newTs)
        output.write()
        self._store()

    def closeOutputSetsStep(self):
        self.InterpolatedTiltSeries.setStreamState(Set.STREAM_CLOSED)
        self.InterpolatedTiltSeries.write()
        self._store()

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        for ts in self.inputSetOfTiltSeries.get():
            if not ts.getFirstItem().hasTransform():
                validateMsgs.append("Some tilt-series from the input set "
                                    "are missing a transformation matrix.")
                break

        return validateMsgs

    def _summary(self):
        summary = []
        if self.InterpolatedTiltSeries:
            summary.append("Input tilt-series: %d\nInterpolations applied: %d\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.InterpolatedTiltSeries.getSize()))
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.InterpolatedTiltSeries:
            methods.append("The interpolation has been computed for %d "
                           "tilt-series using the IMOD *newstack* command.\n"
                           % (self.InterpolatedTiltSeries.getSize()))
        return methods
