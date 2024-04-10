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
from pyworkflow.object import Set
from pwem.emlib.image import ImageHandler
import tomo.objects as tomoObj

from .. import Plugin
from .protocol_base import ProtImodBase, OUTPUT_TILTSERIES_NAME, EXT_MRCS_TS_ODD_NAME, EXT_MRCS_TS_EVEN_NAME, ODD, \
    MRCS_EXT, EVEN


class ProtImodXraysEraser(ProtImodBase):
    """
    Erase X-rays from aligned tilt-series based on the IMOD procedure.
    More info:
            https://bio3d.colorado.edu/imod/doc/man/ccderaser.html
    """

    _label = 'X-rays eraser'
    _devStatus = BETA
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: tomoObj.SetOfTiltSeries}

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        form.addParam('peakCriterion',
                      params.FloatParam,
                      default=8.0,
                      label='Peak criterion',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Criterion # of SDs above local mean for erasing '
                           'peak based on intensity (the default is 8 SDs)')

        form.addParam('diffCriterion',
                      params.FloatParam,
                      default=6.0,
                      label='Difference criterion',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Criterion # of SDs above mean pixel-to-pixel '
                           'difference for erasing a peak based on '
                           'differences (the default is 6 SDs).')

        form.addParam('maximumRadius',
                      params.FloatParam,
                      default=4.2,
                      label='Maximum radius (px)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Maximum radius of peak area to erase (the '
                           'default is 4.2 pixels).')

        form.addParam('bigDiffCriterion',
                      params.IntParam,
                      default=19,
                      label='Big difference criterion',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='An extra-large peak will be erased only if the '
                           'value for the maximum difference between '
                           'adjacent pixels, averaged over the most extreme '
                           'one-fourth of the pixels in the patch, '
                           'exceeds this criterion, evaluated as the '
                           'number of SDs above the mean absolute difference '
                           'between adjacent pixels in the scan area. The '
                           'default is 19.  This high a value is needed '
                           'to prevent gold erasure on low-noise data sets '
                           'with small gold particles, and a lower value '
                           'may be needed to make extra-large peak removal '
                           'useful.')

        form.addParam('processOddEven',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=True,
                      label='Apply to odd/even',
                      help='If True, the full tilt series and the associated odd/even tilt series will be processed. '
                           'The filter applied to the odd/even tilt series will be exactly the same.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        for tsId in self.tsDict.keys():
            self._insertFunctionStep(self.convertInputStep, tsId)
            self._insertFunctionStep(self.eraseXraysStep, tsId)
            self._insertFunctionStep(self.createOutputStep, tsId)
        self._insertFunctionStep(self.closeOutputStep)

    def _initialize(self):
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self.inputSetOfTiltSeries.get()}

    def convertInputStep(self, tsId, **kwargs):
        oddEvenFlag = self.applyToOddEven(self.inputSetOfTiltSeries.get())
        super().convertInputStep(tsId,
                                 imodInterpolation=None,
                                 generateAngleFile=False,
                                 oddEven=oddEvenFlag)

    def eraseXraysStep(self, tsId):
        ts = self.tsDict[tsId]
        firstItem = ts.getFirstItem()

        paramsCcderaser = {
            'input': self.getTmpOutFile(tsId),
            'output': self.getExtraOutFile(tsId),
            'findPeaks': 1,
            'peakCriterion': self.peakCriterion.get(),
            'diffCriterion': self.diffCriterion.get(),
            'growCriterion': 4.,
            'scanCriterion': 3.,
            'maximumRadius': self.maximumRadius.get(),
            'giantCriterion': 12.,
            'extraLargeRadius': 8.,
            'bigDiffCriterion': self.bigDiffCriterion.get(),
            'annulusWidth': 2.0,
            'xyScanSize': 100,
            'edgeExclusionWidth': 4,
            'pointModel': self.getExtraOutFile(tsId, suffix="fid", ext="mod"),
            'borderSize': 2,
            'polynomialOrder': 2,
        }

        argsCcderaser = "-InputFile %(input)s " \
                        "-OutputFile %(output)s " \
                        "-FindPeaks %(findPeaks)d " \
                        "-PeakCriterion %(peakCriterion).2f " \
                        "-DiffCriterion %(diffCriterion).2f " \
                        "-GrowCriterion %(growCriterion).2f " \
                        "-ScanCriterion %(scanCriterion).2f " \
                        "-MaximumRadius %(maximumRadius).2f " \
                        "-GiantCriterion %(giantCriterion).2f " \
                        "-ExtraLargeRadius %(extraLargeRadius).2f " \
                        "-BigDiffCriterion %(bigDiffCriterion).2f " \
                        "-AnnulusWidth %(annulusWidth).2f " \
                        "-XYScanSize %(xyScanSize)d " \
                        "-EdgeExclusionWidth %(edgeExclusionWidth)d " \
                        "-BorderSize %(borderSize)d " \
                        "-PolynomialOrder %(polynomialOrder)d "

        Plugin.runImod(self, 'ccderaser', argsCcderaser % paramsCcderaser)

        if self.applyToOddEven(ts):
            oddFn = firstItem.getOdd().split('@')[1]
            evenFn = firstItem.getEven().split('@')[1]
            paramsCcderaser['input'] = oddFn
            paramsCcderaser['output'] = self.getExtraOutFile(tsId, suffix=ODD, ext=MRCS_EXT)
            Plugin.runImod(self, 'ccderaser', argsCcderaser % paramsCcderaser)
            paramsCcderaser['input'] = evenFn
            paramsCcderaser['output'] = self.getExtraOutFile(tsId, suffix=EVEN, ext=MRCS_EXT)
            Plugin.runImod(self, 'ccderaser', argsCcderaser % paramsCcderaser)

    def createOutputStep(self, tsId):
        output = self.getOutputSetOfTiltSeries(self.inputSetOfTiltSeries.get())

        ts = self.tsDict[tsId]
        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        output.append(newTs)

        ih = ImageHandler()

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True, copyTM=True)
            newTi.setAcquisition(tiltImage.getAcquisition())
            newTi.setLocation(index + 1, self.getExtraOutFile(tsId))

            if self.applyToOddEven(ts):
                locationOdd = index + 1, self.getExtraOutFile(tsId, suffix=ODD, ext=MRCS_EXT)
                locationEven = index + 1, self.getExtraOutFile(tsId, suffix=EVEN, ext=MRCS_EXT)
                newTi.setOddEven([ih.locationToXmipp(locationOdd), ih.locationToXmipp(locationEven)])
            else:
                newTi.setOddEven([])

            newTs.append(newTi)

        newTs.write(properties=False)
        output.update(newTs)
        output.write()
        self._store()

    def closeOutputStep(self):
        if self.TiltSeries:
            self.TiltSeries.setStreamState(Set.STREAM_CLOSED)
            self.TiltSeries.write()
        self._store()

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append("Input tilt-series: %d\nX-rays erased output "
                           "tilt series: %d"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.TiltSeries.getSize()))
        else:
            summary.append("Outputs are not ready yet.")

        return summary

    def _methods(self):
        methods = []
        if self.TiltSeries:
            methods.append("The x-rays artifacts have been erased for %d "
                           "tilt-series using the IMOD *ccderaser* command.\n"
                           % (self.TiltSeries.getSize()))

        return methods
