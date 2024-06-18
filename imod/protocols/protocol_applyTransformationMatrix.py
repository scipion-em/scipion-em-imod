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
from tomo.objects import TiltSeries, TiltImage
from .. import Plugin, utils
from .protocol_base import ProtImodBase, XF_EXT, ODD, EVEN


class ProtImodApplyTransformationMatrix(ProtImodBase):
    """
    Compute the interpolated tilt-series from its transform matrix.
    The protocol makes use of the IMod command newstack
    More info:
        https://bio3d.colorado.edu/imod/doc/man/newstack.html

    Generally, the tilt series has an associated transformation matrix
    which contains the alignment information. The transformation matrix
    is usually associated but not applied to avoid to accumulate interpolation
    errors during the image processing. This protocol allows to apply
    the transformation matrix to the tilt series
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
                      label='Tilt-series to apply the transformation matrix')

        form.addParam('binning', params.IntParam,
                      default=1,
                      label='Binning for the interpolated',
                      help='Binning to be applied to the interpolated tilt-series in IMOD '
                           'convention. \n'
                           'Binning is an scaling factor given by an integer greater than 1. '
                           'IMOD uses ordinary binning (with antialiasing filter) to reduce images in size by the given factor. '
                           'The value of a binned pixel is the average of pixel values in each block '
                           'of pixels being binned. Binning is applied before all other image '
                           'transformations.')

        form.addParam('taperInside',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=True,
                      label='Taper inwards from the edge?',
                      help='When the image is transformed areas with no information are filled in (e.g. because of rotation).'
                      'Decide whether tapering is done inwards or outwards from the edge.')
        
        form.addParam('linear',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=False,
                      label='Linear interpolation?',
                      help='From newstack man page: Use linear instead of cubic interpolation to transform images. '
                            'Linear interpolation is more suitable when images are very noisy, '
                            'but cubic interpolation will preserve fine detail better when noise is not an issue.')
        
        form.addParam('processOddEven',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=True,
                      label='Apply to odd/even',
                      help='If True, the full tilt series and the associated odd/even tilt series will be processed. '
                           'The transformations applied to the odd/even tilt series will be exactly the same.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        for tsId in self.tsDict.keys():
            self._insertFunctionStep(self.generateTransformFileStep, tsId)
            self._insertFunctionStep(self.computeAlignmentStep, tsId)
            self._insertFunctionStep(self.generateOutputStackStep, tsId)
            self._insertFunctionStep(self.createOutputFailedStep, tsId)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions ------------------------------
    def _initialize(self):
        self._failedTs = []
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self.inputSetOfTiltSeries.get()}

    def generateTransformFileStep(self, tsId):
        ts = self.tsDict[tsId]
        self.genTsPaths(tsId)
        utils.genXfFile(ts, self.getExtraOutFile(tsId, ext=XF_EXT))

    @ProtImodBase.tryExceptDecorator
    def computeAlignmentStep(self, tsId):
        ts = self.tsDict[tsId]
        firstItem = ts.getFirstItem()
        binning = self.binning.get()

        paramsAlignment = {
            'input': firstItem.getFileName(),
            'output': self.getExtraOutFile(tsId),
            'xform': self.getExtraOutFile(tsId, ext=XF_EXT),
            'bin': binning,
            'imagebinned': 1.0,
            'taper': "1,1" if self.taperInside else "1,0"
        }

        argsAlignment = "-input %(input)s " \
                        "-output %(output)s " \
                        "-xform %(xform)s " \
                        "-bin %(bin)d " \
                        "-antialias -1 " \
                        "-imagebinned %(imagebinned)s " \
                        "-taper %(taper)s "
        
        if self.linear.get():
            argsAlignment += "-linear "

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

        if self.applyToOddEven(ts):
            oddFn = firstItem.getOdd().split('@')[1]
            evenFn = firstItem.getEven().split('@')[1]
            paramsAlignment['input'] = oddFn
            paramsAlignment['output'] = self.getExtraOutFile(tsId, suffix=ODD)
            Plugin.runImod(self, 'newstack', argsAlignment % paramsAlignment)
            paramsAlignment['input'] = evenFn
            paramsAlignment['output'] = self.getExtraOutFile(tsId, suffix=EVEN)
            Plugin.runImod(self, 'newstack', argsAlignment % paramsAlignment)

    def generateOutputStackStep(self, tsId):
        ts = self.tsDict[tsId]
        outputLocation = self.getExtraOutFile(tsId)

        if os.path.exists(outputLocation):
            output = self.getOutputInterpolatedSetOfTiltSeries(self.inputSetOfTiltSeries.get())

            binning = self.binning.get()

            newTs = TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            newTs.setInterpolated(True)
            acq = newTs.getAcquisition()
            acq.setTiltAxisAngle(0.)  # 0 because TS is aligned
            newTs.setAcquisition(acq)
            output.append(newTs)

            if binning > 1:
                newTs.setSamplingRate(ts.getSamplingRate() * binning)

            ih = ImageHandler()

            index = 1
            for tiltImage in ts:
                if tiltImage.isEnabled():
                    newTi = TiltImage()
                    newTi.copyInfo(tiltImage, copyId=False, copyTM=False)
                    acq = tiltImage.getAcquisition()
                    acq.setTiltAxisAngle(0.)
                    newTi.setAcquisition(acq)
                    newTi.setLocation(index, outputLocation)
                    if self.applyToOddEven(ts):
                        locationOdd = index, (self.getExtraOutFile(tsId, suffix=ODD))
                        locationEven = index, (self.getExtraOutFile(tsId, suffix=EVEN))
                        newTi.setOddEven([ih.locationToXmipp(locationOdd), ih.locationToXmipp(locationEven)])
                    else:
                        newTi.setOddEven([])

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

    def createOutputFailedStep(self, tsId):
        ts = self.tsDict[tsId]
        super().createOutputFailedSet(ts)

    def closeOutputSetsStep(self):
        for _, output in self.iterOutputAttributes():
            output.setStreamState(Set.STREAM_CLOSED)
            output.write()
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
