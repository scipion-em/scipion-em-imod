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
from tomo.protocols import ProtTomoBase
from tomo.convert import writeTiStack


class ProtImodXcorr(pyem.EMProtocol, ProtTomoBase):
    """
    Tilt-series' cross correlation alignment based on the IMOD procedure.

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
                       help='Binning to be applied to the interpolated tilt-series. '
                            'Must be a integer bigger than 1')

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
            self._insertFunctionStep('computeInterpolatedStackStep')
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
            tiList.reverse()

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

            paramsXftoxg ={
                'input': '%s.prexf' % tsId,
                'goutput': '%s.prexg' % tsId,
            }
            argsXftoxg = "-input %(input)s " \
                        "-goutput %(goutput)s"
            self.runJob('xftoxg', argsXftoxg % paramsXftoxg, cwd=workingFolder)

        # Generate output tilt series
        outputSetOfTiltSeries = self.getOutputSetOfTiltSeries()

        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            alignmentMatrix = self.formatTransformationMatrix(self._getExtraPath('%s/%s.prexg' % (tsId, tsId)))
            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            outputSetOfTiltSeries.append(newTs)
            i = 0

            # For each tilt image in the series, assign its transform matrix
            for tiltImage in ts:
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)
                newTi.setLocation(tiltImage.getLocation())

                # Set the tansformation matrix
                transform = data.Transform()
                transform.setMatrix(alignmentMatrix[:, :, i])
                newTi.setTransform(transform)
                i += 1
                newTs.append(newTi)
            outputSetOfTiltSeries.update(newTs)  # update items and size info
        self._store()

    def computeInterpolatedStackStep(self):
        outputInterpolatedSetOfTiltSeries = self.getOutputInterpolatedSetOfTiltSeries()

        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            outputInterpolatedSetOfTiltSeries.append(newTs)
            workingFolder = self._getExtraPath(tsId)
            paramsAlginment = {
                'input': "%s.st" % tsId,
                'output': '%s_preali.st' % tsId,
                'xform': "%s.prexg" % tsId,
                'bin': int(self.binning.get()),
                'mode': 0,
                'float': 2,
                'imagebinned': 1.0}
            argsAlignment = "-input %(input)s " \
                            "-output %(output)s " \
                            "-xform %(xform)s " \
                            "-bin %(bin)d " \
                            "-mode %(mode)s " \
                            "-float %(float)s " \
                            "-imagebinned %(imagebinned)s"

            self.runJob('newstack', argsAlignment % paramsAlginment, cwd=workingFolder)
            for tiltImage in ts:
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)
                newTi.setLocation(os.path.join(workingFolder, '%s_preali.st' % tsId))
                newTs.append(newTi)
            if self.binning > 1:
                newTs.setSamplingRate(ts.getSamplingRate()*int(self.binning.get()))
            outputInterpolatedSetOfTiltSeries.update(newTs)  # update items and size info
        self._store()

    def _createOutputStep(self):
        pass

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

    def getOutputSetOfTiltSeries(self):
        if not hasattr(self, "outputSetOfTiltSeries"):
            outputSetOfTiltSeries = self._createSetOfTiltSeries()
            outputSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            self._defineOutputs(outputSetOfTiltSeries=outputSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfTiltSeries)
        return self.outputSetOfTiltSeries

    def getOutputInterpolatedSetOfTiltSeries(self):
        if not hasattr(self, "outputInterpolatedSetOfTiltSeries"):
            outputInterpolatedSetOfTiltSeries = self._createSetOfTiltSeries()
            outputInterpolatedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            if self.binning > 1:
                samplingRate = self.inputSetOfTiltSeries.get().getSamplingRate()
                samplingRate *= self.binning.get()
                outputInterpolatedSetOfTiltSeries.setSamplingRate(samplingRate)
            self._defineOutputs(outputInterpolatedSetOfTiltSeries=outputInterpolatedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputInterpolatedSetOfTiltSeries)
        return self.outputInterpolatedSetOfTiltSeries
