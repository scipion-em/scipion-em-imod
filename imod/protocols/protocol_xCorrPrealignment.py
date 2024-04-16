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
import numpy as np

from pyworkflow import BETA
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
from pwem.objects import Transform
from pwem.emlib.image import ImageHandler
import tomo.objects as tomoObj

from .. import Plugin, utils
from .protocol_base import ProtImodBase
from .protocol_forms import CommonIMODforms


class ProtImodXcorrPrealignment(ProtImodBase, CommonIMODforms):
    """
    Tilt-series cross correlation alignment based on the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/tiltxcorr.html

    Tiltxcorr uses cross-correlation to find an initial translational
    alignment between successive images of a tilt series.  For a given pair
    of images, it stretches the image with the larger tilt angle perpendic-
    ular to the tilt axis, by an amount equal to the ratio of the cosines
    of the two tilt angles (cosine stretch).  The stretched image is corre-
    lated with the other image, and the position of the peak of the corre-
    lation indicates the relative shift between the images.


    """

    _label = 'Coarse prealignment'
    _devStatus = BETA
    _possibleOutputs = {"TiltSeries": tomoObj.SetOfTiltSeries}

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-series to be prealigned')

        form.addParam('cumulativeCorr',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Use cumulative correlation?',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='The program will take the image at zero tilt as the first'
                           'reference, and correlate it with the image at the next most negative tilt.'
                           'It will then add the aligned image to the first reference to make the '
                           'reference for the next tilt.  At each tilt, the reference will be the sum of '
                           'images that have already been aligned. When the most negative tilt angle is '
                           'reached, the procedure is repeated from the zero-tilt view to more positive tilt angles.')

        form.addParam('computeAlignment',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Generate interpolated tilt-series?',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Generate and save the interpolated tilt-series applying the obtained transformation '
                           'matrices.\n'
                           'By default, the output of this protocol will be a tilseries that will have associated'
                           'the alignment information as a transformation matrix. When this option is set as Yes, '
                           'then a second output, called interpolated tilt series, is generated. The interpolated tilt '
                           'series should be used for visualization purpose but not for image processing')

        form.addParam('binning',
                      params.IntParam,
                      condition='computeAlignment==0',
                      default=1,
                      label='Binning for the interpolated',
                      help='Binning to be applied to the interpolated  tilt-series in IMOD '
                           'convention. \n'
                           'Binning is an scaling factor given by an integer greater than 1. '
                           'IMOD uses ordinary binning to reduce images in size by the given factor. '
                           'The value of a binned pixel is the average of pixel values in each block '
                           'of pixels being binned. Binning is applied before all other image '
                           'transformations.')

        form.addParam('Trimming parameters', params.LabelParam,
                      label='Tilt axis angle detected from import. In case another value is desired please adjust the '
                            'number below')

        form.addParam('tiltAxisAngle',
                      params.FloatParam,
                      allowsNull=True,
                      label='Tilt axis angle (degrees)',
                      help='The tilt axis angle is the tilt axis rotation relative to the Y axis of the image.'
                           'If it was not properly set in the import of the tilt series, or the imported'
                           'information is not correct you have the chance to correct at in this point.'
                           'Usually, it will be 90 degrees less than the RotationAngle in a '
                           'system with no axis inversions')

        trimming = form.addGroup('Trimming parameters', expertLevel=params.LEVEL_ADVANCED)

        self.trimimgForm(trimming, pxTrimCondition='False', correlationCondition='True', levelType=params.LEVEL_ADVANCED)
        self.filteringParametersForm(form, condition='True', levelType=params.LEVEL_ADVANCED)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep(self.convertInputStep, ts.getObjId())
            self._insertFunctionStep(self.computeXcorrStep, ts.getObjId())
            self._insertFunctionStep(self.generateOutputStackStep,
                                     ts.getObjId())
            if self.computeAlignment.get() == 0:
                self._insertFunctionStep(self.computeInterpolatedStackStep,
                                         ts.getObjId())
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def computeXcorrStep(self, tsObjId):
        """Compute transformation matrix for each tilt series"""
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        tiltAxisAngle = self.getTiltAxisOrientation(ts)

        paramsXcorr = {
            'input': os.path.join(tmpPrefix,
                                  ts.getFirstItem().parseFileName()),
            'output': os.path.join(extraPrefix,
                                   ts.getFirstItem().parseFileName(extension=".prexf")),
            'tiltfile': os.path.join(tmpPrefix,
                                     ts.getFirstItem().parseFileName(extension=".tlt")),
            'rotationAngle': tiltAxisAngle,
            'filterSigma1': self.filterSigma1.get(),
            'filterSigma2': self.filterSigma2.get(),
            'filterRadius1': self.filterRadius1.get(),
            'filterRadius2': self.filterRadius2.get()
        }

        argsXcorr = "-input %(input)s " \
                    "-output %(output)s " \
                    "-tiltfile %(tiltfile)s " \
                    "-RotationAngle %(rotationAngle).2f " \
                    "-FilterSigma1 %(filterSigma1).3f " \
                    "-FilterSigma2 %(filterSigma2).3f " \
                    "-FilterRadius1 %(filterRadius1).3f " \
                    "-FilterRadius2 %(filterRadius2).3f "

        if self.cumulativeCorr == 0:
            argsXcorr += "-CumulativeCorrelation "

        xdim = ts.getDim()[0]
        ydim = ts.getDim()[1]
        xmin = self.xmin.get()
        xmax = self.xmax.get()
        ymin = self.ymin.get()
        ymax = self.ymax.get()
        if (not (xmin is None)) or (not (xmax is None)) or (not (ymin is None)) or (not (ymax is None)):
            if xmin is None:
                xmin = 0
            if xmax is None:
                xmax = xdim - 1
            if ymin is None:
                ymin = 0
            if ymax is None:
                ymax = ydim - 1

            argsXcorr += " -xminmax %i,%i " % (xmin, xmax)
            argsXcorr += " -yminmax %i,%i " % (ymin, ymax)

        # Excluded views
        excludedViews = ts.getExcludedViewsIndex(caster=str)
        if len(excludedViews):
            argsXcorr += f"-SkipViews {','.join(excludedViews)} "

        Plugin.runImod(self, 'tiltxcorr', argsXcorr % paramsXcorr)

        paramsXftoxg = {
            'input': os.path.join(extraPrefix,
                                  ts.getFirstItem().parseFileName(extension=".prexf")),
            'goutput': os.path.join(extraPrefix,
                                    ts.getFirstItem().parseFileName(extension=".prexg")),
        }
        argsXftoxg = "-input %(input)s " \
                     "-NumberToFit 0 " \
                     "-goutput %(goutput)s "
        Plugin.runImod(self, 'xftoxg', argsXftoxg % paramsXftoxg)

    def getTiltAxisOrientation(self, ts):
        if self.tiltAxisAngle.get():
            return self.tiltAxisAngle.get()
        else:
            return ts.getAcquisition().getTiltAxisAngle()

    def generateOutputStackStep(self, tsObjId):
        """ Generate tilt-serie with the associated transform matrix """
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        output = self.getOutputSetOfTiltSeries(self.inputSetOfTiltSeries.get())

        alignmentMatrix = utils.formatTransformationMatrix(
            os.path.join(extraPrefix,
                         ts.getFirstItem().parseFileName(extension=".prexg")))

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        newTs.getAcquisition().setTiltAxisAngle(self.getTiltAxisOrientation(ts))

        output.append(newTs)

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True, copyTM=False)

            if tiltImage.hasTransform():
                transform = Transform()
                previousTransform = tiltImage.getTransform().getMatrix()
                newTransform = alignmentMatrix[:, :, index]
                previousTransformArray = np.array(previousTransform)
                newTransformArray = np.array(newTransform)
                outputTransformMatrix = np.matmul(newTransformArray, previousTransformArray)
                transform.setMatrix(outputTransformMatrix)
                newTi.setTransform(transform)

            else:
                transform = Transform()
                newTransform = alignmentMatrix[:, :, index]
                newTransformArray = np.array(newTransform)
                transform.setMatrix(newTransformArray)
                newTi.setTransform(transform)

            newTi.setAcquisition(tiltImage.getAcquisition())
            newTi.setLocation(tiltImage.getLocation())

            newTs.append(newTi)

        newTs.write(properties=False)

        output.update(newTs)
        output.write()

        self._store()

    def computeInterpolatedStackStep(self, tsObjId):
        output = self.getOutputInterpolatedSetOfTiltSeries(self.inputSetOfTiltSeries.get())

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsAlignment = {
            'input': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName()),
            'output': os.path.join(extraPrefix, ts.getFirstItem().parseFileName()),
            'xform': os.path.join(extraPrefix,
                                  ts.getFirstItem().parseFileName(extension=".prexg")),
            'bin': self.binning.get(),
            'imagebinned': 1.0
        }
        argsAlignment = "-input %(input)s " \
                        "-output %(output)s " \
                        "-xform %(xform)s " \
                        "-bin %(bin)d " \
                        "-imagebinned %(imagebinned)s " \
                        "-antialias -1 " \
                        "-float 2 "

        Plugin.runImod(self, 'newstack', argsAlignment % paramsAlignment)

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        newTs.getAcquisition().setTiltAxisAngle(self.getTiltAxisOrientation(ts))

        newTs.setInterpolated(True)
        output.append(newTs)

        if self.binning > 1:
            newTs.setSamplingRate(ts.getSamplingRate() * self.binning.get())

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(index + 1,
                              (os.path.join(extraPrefix,
                                            tiltImage.parseFileName())))
            if self.binning > 1:
                newTi.setSamplingRate(tiltImage.getSamplingRate() * self.binning.get())
            newTs.append(newTi)

        ih = ImageHandler()
        x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
        newTs.setDim((x, y, z))

        newTs.write(properties=False)

        output.update(newTs)
        output.write()
        self._store()

    def closeOutputSetsStep(self):
        self.TiltSeries.setStreamState(Set.STREAM_CLOSED)
        self.TiltSeries.write()
        if self.computeAlignment.get() == 0:
            self.InterpolatedTiltSeries.setStreamState(Set.STREAM_CLOSED)
            self.InterpolatedTiltSeries.write()

        self._store()

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append("Input tilt-series: %d\nTransformation matrices "
                           "calculated: %d"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.TiltSeries.getSize()))
            if self.InterpolatedTiltSeries:
                summary.append("Interpolated tilt-series: %d"
                               % self.InterpolatedTiltSeries.getSize())
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.TiltSeries:
            methods.append("The transformation matrix has been calculated for %d "
                           "tilt-series using the IMOD *tiltxcorr* command."
                           % (self.TiltSeries.getSize()))
        return methods
