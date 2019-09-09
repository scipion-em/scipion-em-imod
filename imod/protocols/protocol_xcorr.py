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
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os

import pyworkflow as pw
import pyworkflow.protocol.params as params
import numpy as np

from tomo.objects import TiltSeriesDict, TiltSeries, Tomogram
from tomo.protocols import ProtTomoReconstruct
from tomo.convert import writeTiStack


class ProtImodXcorr(ProtTomoReconstruct):
    """
    Tilt-series's cross correlation alignment based on the IMOD procedure.

    More info:
        https://bio3d.colorado.edu/imod/doc/etomoTutorial.html
    """

    _label = 'xcorr prealignment'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('computeAlignment', params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Compute alignment', important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Compute and save the aligned tilt-series.')

        form.addParam('inputTiltSeries', params.PointerParam,
                       pointerClass='TiltSeries',
                       important=True,
                       label='Input tilt-Series')

        form.addParam('rotationAngle', params.FloatParam,
                       label='Tilt rotation angle (deg)',
                       help='Angle from the vertical to the tilt axis in raw '
                            'images.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        ts = self.inputTiltSeries.get()
        tsId = ts.getTsId()
        self._insertFunctionStep('convertInputStep', tsId)
        if self.computeAlignment.get == 0:
            self._insertFunctionStep('computeXcorrStep', tsId)
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsId):
        ts = self.inputTiltSeries.get()
        workingFolder = self._getExtraPath(tsId)
        prefix = os.path.join(workingFolder, tsId)
        pw.utils.makePath(workingFolder)

        tiList = [ti.clone() for ti in ts]
        tiList.sort(key=lambda ti: ti.getTiltAngle())

        writeTiStack(tiList,
                     outputStackFn=prefix + '.st',
                     outputTltFn=prefix + '.rawtlt')

    def computeXcorrStep(self, tsId):
        workingFolder = self._getExtraPath(tsId)
        paramsXcorr = {
            'input': '%s.st' % tsId,
            'output': '%s.prexf' % tsId,
            'tiltfile': '%s.rawtlt' % tsId,
            'RotationAngle': 12.5,
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

        self.allocateTransformMatrix(workingFolder + '%s.prexf' % tsId, self.inputTiltSeries.get())

    def computeStackAlignment(self, tsId):
        workingFolder = self._getExtraPath(tsId)
        paramsAlginment = {
            'input': "%s.st" % tsId,
            'output': '%s.preali' % tsId,
            'bin': '%f' % self.binning,
            'xform': "%s.prexf" % tsId,
        }
        argsAlignment = "-input %(input)s " \
                        "-output %(output)s " \
                        "-bin %(bin)f " \
                        "-xform %(xform)s "
        self.runJob('newstack', argsAlignment % paramsAlginment, cwd=workingFolder)

    def _createOutput(self, tomoFn):
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

    # --------------------------- UTILS functions ----------------------------
    def allocateTransformMatrix(self, matrixFile, tiltSeries):
        frameMatrix = np.empty([3, 3])
        for image, line in zip(tiltSeries, matrixFile):
            values = line.split()
            frameMatrix[0, 0] = float(values(0))
            frameMatrix[1, 0] = float(values(1))
            frameMatrix[0, 1] = float(values(2))
            frameMatrix[1, 1] = float(values(3))
            frameMatrix[0, 2] = float(values(4))
            frameMatrix[1, 2] = float(values(5))
            frameMatrix[3, 0] = 0.0
            frameMatrix[3, 1] = 0.0
            frameMatrix[3, 2] = 1.0
            image.setTransform(matrixFile)
