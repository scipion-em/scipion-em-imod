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
import pyworkflow.em as pyem
import pyworkflow.em.data as data
import pyworkflow.protocol.params as params
import tomo.objects as tomoObj

from tomo.protocols import ProtTomoBase, ProtTomoReconstruct
from tomo.convert import writeTiStack


class ProtImodXcorr(pyem.EMProtocol, ProtTomoBase):
    """
    Tilt-series's cross correlation alignment based on the IMOD procedure.

    More info:
        https://bio3d.colorado.edu/imod/doc/etomoTutorial.html
    """

    _label = 'xcorr prealignment'

    def __init__(self, **kwargs):
        pyem.EMProtocol.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries', params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-Series')

        form.addParam('computeAlignment', params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Generate interpolated tilt-series', important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Generate and save the interpolated tilt-series applying the'
                           'obtained transformation matrices.')

        group = form.addGroup('Interpolated tilt-series',
                      condition='computeAlignment==0')

        group.addParam('binning', params.FloatParam,
                       default=1.0,
                       label='Binning',
                       help='Binning to be applied to the interpolated tilt-series')

        form.addParam('rotationAngle',
                      params.FloatParam,
                      label='Tilt rotation angle (deg)',
                      default='0.0',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Angle from the vertical to the tilt axis in raw images.")

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('computeXcorrStep')
        if self.computeAlignment.get() == 0:
            self._insertFunctionStep('computeInterpolatedStack')
        self._insertFunctionStep('_createOutputStep')

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)
            prefix = os.path.join(workingFolder, tsId)
            pw.utils.makePath(workingFolder)

            tiList = [ti.clone() for ti in ts]
            tiList.sort(key=lambda ti: ti.getTiltAngle())

            writeTiStack(tiList,
                         outputStackFn=prefix + '.st',
                         outputTltFn=prefix + '.rawtlt')

    def computeXcorrStep(self):
        # Compute transformation matrix for each tilt series
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)
            paramsXcorr = {
                'input': '%s.st' % tsId,
                'output': '%s.prexf' % tsId,
                'tiltfile': '%s.rawtlt' % tsId,
                'RotationAngle': self.rotationAngle.get(),
                'FilterSigma1': 0.03,
                'FilterSigma2': 0.05,
                'FilterRadius2': 0.25
            }
            argsXcorr = "-input %(input)s " \
                        "-output %(output)s " \
                        "-tiltfile %(tiltfile)s " \
                        "-RotationAngle %(RotationAngle)f " \
                        "-FilterSigma1 %(FilterSigma1)f " \
                        "-FilterSigma2 %(FilterSigma2)f " \
                        "-FilterRadius2 %(FilterRadius2)f"

            self.runJob('tiltxcorr', argsXcorr % paramsXcorr, cwd=workingFolder)

        self.outputSetOfTiltSeries = self._createSetOfTiltSeries()
        self.outputSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
        self.outputSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())

        # Generate output tilt series
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            alignmentMatrix = self.formatTransformationMatrix(self._getExtraPath('%s/%s.prexf' % (tsId, tsId)))

            # Create new  output tiltSeries
            tsObj = tomoObj.TiltSeries(tsId=tsId)
            # we need this to set mapper before adding any item
            self.outputSetOfTiltSeries.append(tsObj)
            i = 0

            # For each tilt series image in the input
            for tiltImage in ts:

                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)

                # Set the tansformation matrix
                transform = data.Transform(alignmentMatrix[:, :, i])
                newTi.setTransform(transform)
                i += 1
                tsObj.append(newTi)
            self.outputSetOfTiltSeries.update(tsObj)  # update items and size info

    def computeInterpolatedStack(self):
        self.outputInterpolatedSetOfTiltSeries = self._createSetOfTiltSeries()
        self.outputInterpolatedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
        self.outputInterpolatedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())

        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            tsObj = tomoObj.TiltSeries(tsId=tsId)
            workingFolder = self._getExtraPath(tsId)
            paramsAlginment = {
                'input': "%s.st" % tsId,
                'output': '%s.preali' % tsId,
                'xform': "%s.prexf" % tsId,
                'bin': int(self.binning.get())
            }
            argsAlignment = "-input %(input)s " \
                            "-output %(output)s " \
                            "-xform %(xform)s " \
                            "-bin %(bin)d"
            self.runJob('newstack', argsAlignment % paramsAlginment, cwd=workingFolder)

    def _createOutputStep(self):
        if self.computeAlignment.get() != 1:
            samplingRate = self.inputTiltSeries.get().getSamplingRate()
            outTiltSeries.copyInfo(self.inputTiltSeries.get())
            interpolatedPath =os.path.join(self._getExtraPath(), '%s.preali' % tsId)
            outTiltSeries.writeStack(interpolatedPath)
            if self.binning > 1:
                samplingRate *= self.binning.get()
            outTiltSeries.setSamplingRate(samplingRate)
       # for image in outTiltSeries.iterItems():
       #     print(image.getTransform())
        self._defineOutputs(outputTiltSeries=self.outputSetOfTiltSeries)
        self._defineSourceRelation(self.inputSetOfTiltSeries, self.outputSetOfTiltSeries)

    # --------------------------- UTILS functions ----------------------------
    def formatTransformationMatrix(self, matrixFile):
        with open(matrixFile, "r") as matrix:
            lines = matrix.readlines()
        numberLines = len(lines)
        frameMatrix = np.empty([3, 3, numberLines])
        i = 0
        for line in lines:
            values = line.split()
            frameMatrix[0, 0, i] = float(values[0])
            frameMatrix[1, 0, i] = float(values[1])
            frameMatrix[0, 1, i] = float(values[2])
            frameMatrix[1, 1, i] = float(values[3])
            frameMatrix[0, 2, i] = float(values[4])
            frameMatrix[1, 2, i] = float(values[5])
            frameMatrix[2, 0, i] = 0.0
            frameMatrix[2, 1, i] = 0.0
            frameMatrix[2, 2, i] = 1.0
            i += 1

        return frameMatrix
"""
        inputTs = self.inputTiltSeries.get()
        outTomos = self._createSetOfTomograms()
        samplingRate = inputTs.getSamplingRate()


        if self.binning > 1:
            samplingRate *= self.binning.get()

        outTomos.setSamplingRate(samplingRate)

        t = Tomogram(location=tomoFn + ':mrc')
        t.setObjId(inputTs.getObjId())
        t.setTsId(inputTs.getTsId())
        outTomos.append(t)

        self._defineOutputs(outputTomograms=outTomos)
        self._defineSourceRelation(self.inputTiltSeries, outTomos)
        """