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

import pyworkflow.protocol.params as params
from pyworkflow.object import Set
from pwem.emlib.image import ImageHandler as ih
import tomo.objects as tomoObj

from .. import Plugin
from .protocol_base import (ProtImodBase, OUTPUT_TILTSERIES_NAME,
                            ODD, EVEN, MOD_EXT)


class ProtImodXraysEraser(ProtImodBase):
    """
    Erase Xrays from aligned tilt-series based on the IMOD procedure.
    More info:
            https://bio3d.colorado.edu/imod/doc/man/ccderaser.html

    This program replaces deviant pixels with interpolated values from
    surrounding pixels. It is designed to correct defects in electron
    microscope images from CCD cameras. It can use two algorithms to
    automatically remove peaks in intensity caused by X-rays.\n

    The automatic removal of X-rays works by dividing the area of each
    image into patches for scanning. The mean and standard deviation (SD)
    of the pixels in a patch are computed. The patch is then scanned for
    pixels that deviate from the mean by more than a criterion number of
    SDs (the scan criterion, a relatively low number to keep from missing
    peaks). When such a pixel is found, the program searches neighboring
    pixels to find a peak in intensity. It then computes the mean and SD
    of pixels in an annulus around the peak and makes sure that the peak
    deviates from this local mean by more than a criterion number of SDs
    (the peak criterion).  Neighboring pixels inside the inner radius of
    the annulus are added to the list of pixels to be replaced if they
    deviate by a lower criterion (the grow criterion).  The patch of pixels
    is then replaced by fitting a polynomial to adjacent pixels and inter-
    polating from the polynomial.  If the peak does not deviate suffi-
    ciently from this local mean, but is stronger than the mean of the scan
    area by the scan criterion plus 1, then the mean and SD is again com-
    puted in a larger annulus.  If the peak deviates from this mean by a
    number of SDs bigger than another criterion for extra-large peaks, a
    patch of pixels is found, but it is replaced only if enough pixels dif-
    fer from adjacent ones by large enough amounts (see the -big option
    below).  The reason for these two stages is that the inner radius for
    the first stage must be set safely smaller than the radius of gold
    beads to avoid erasing part of the beads, whereas the second stage can
    work with larger areas because it has more stringent criteria that
    reject gold beads.

    After the peaks are found in a scanning patch, the program next finds
    the difference between each pixel and the mean of the eight adjacent
    pixels.  The mean and SD of this difference is computed, then pixels
    are sought that deviate from the mean by yet another criterion, the
    difference criterion.  When such a pixel is found, neighboring pixels
    are examined and added to the patch of pixels to replace if their dif-
    ference exceeds the grow criterion.  If the number of pixels in the
    patch does not exceed a specified maximum, replacement proceeds as
    above; otherwise the patch is ignored.

    Two methods are used because the first method is probably more reliable
    for dealing with strong peaks that extend over several pixels, while
    the second method is definitely better for finding small X-rays.

    After all the patches have been scanned for a section, the program then
    searches for single pixels with large interpixel differences at the
    edges of the image, over the width set by the -border option.  A dif-
    ference between a pixel and the mean of whatever adjacent pixels exist
    is computed and its deviation from the overall mean interpixel differ-
    ence is divided by the maximum SD of interpixel differences over all of
    the scans.  When this value exceeds the difference criterion and the
    interpixel difference is greater than that of its neighbors, the pixel
    is replaced with the mean.  This procedure is iterated up to 4 times to
    catch adjacent extreme pixels.

    Tuning the removal of X-rays would primarily involve adjusting two of
    the criteria.  The peak and difference criteria would be adjusted down
    or up to increase or decrease the number of deviant pixels that are
    found.  The grow criterion could also be adjusted down or up depending
    on whether too few or too many pixels are included in a patch that is
    replaced, but this step is not usually done in practice.  If there are
    strong, large artifacts that are not being removed, the big difference
    criterion for extra-large peaks should be lowered first, then if neces-
    sary, the maximum radius and criterion strength for extra-large peaks
    can be adjusted.

    """

    _label = 'X-rays eraser'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: tomoObj.SetOfTiltSeries}

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt Series')

        form.addParam('peakCriterion',
                      params.FloatParam,
                      default=8.0,
                      label='Peak criterion (in std)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Criterion # of SDs above local mean for erasing '
                           'peak based on intensity (the default is 10 SDs)')

        form.addParam('diffCriterion',
                      params.FloatParam,
                      default=6.0,
                      label='Difference criterion  (in std)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Criterion # of SDs above mean pixel-to-pixel '
                           'difference for erasing a peak based on '
                           'differences (the default is 8 SDs).')

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
                      label='Big difference criterion  (in std)',
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
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self._getSetOfInputTS()}

    def convertInputStep(self, tsId, **kwargs):
        oddEvenFlag = self.applyToOddEven(self._getSetOfInputTS())
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
            'pointModel': self.getExtraOutFile(tsId, suffix="fid", ext=MOD_EXT),
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
            paramsCcderaser['output'] = self.getExtraOutFile(tsId, suffix=ODD)
            Plugin.runImod(self, 'ccderaser', argsCcderaser % paramsCcderaser)
            paramsCcderaser['input'] = evenFn
            paramsCcderaser['output'] = self.getExtraOutFile(tsId, suffix=EVEN)
            Plugin.runImod(self, 'ccderaser', argsCcderaser % paramsCcderaser)

    def createOutputStep(self, tsId):
        output = self.getOutputSetOfTiltSeries(self._getSetOfInputTS())

        ts = self.tsDict[tsId]
        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        output.append(newTs)

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True, copyTM=True)
            newTi.setAcquisition(tiltImage.getAcquisition())
            newTi.setLocation(index + 1, self.getExtraOutFile(tsId))

            if self.applyToOddEven(ts):
                locationOdd = index + 1, self.getExtraOutFile(tsId, suffix=ODD)
                locationEven = index + 1, self.getExtraOutFile(tsId, suffix=EVEN)
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
                           % (self._getSetOfInputTS().getSize(),
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
