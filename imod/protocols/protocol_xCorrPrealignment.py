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

import pyworkflow.protocol.params as params
from pwem.objects import Transform
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries

from imod import utils
from imod.protocols.protocol_base import ProtImodBase
from imod.constants import TLT_EXT, PREXF_EXT, PREXG_EXT, OUTPUT_TILTSERIES_NAME


class ProtImodXcorrPrealignment(ProtImodBase):
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
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-series to be prealigned')

        form.addParam('cumulativeCorr',
                      params.BooleanParam,
                      default=False,
                      label='Use cumulative correlation?',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='The program will take the image at zero tilt as the first'
                           'reference, and correlate it with the image at the next '
                           'most negative tilt. It will then add the aligned image '
                           'to the first reference to make the reference for the next '
                           'tilt.  At each tilt, the reference will be the sum of '
                           'images that have already been aligned. When the most '
                           'negative tilt angle is reached, the procedure is repeated '
                           'from the zero-tilt view to more positive tilt angles.')

        form.addParam('computeAlignment',
                      params.BooleanParam,
                      default=False,
                      label='Generate interpolated tilt-series?',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Generate and save the interpolated tilt-series applying '
                           'the obtained transformation matrices.\n'
                           'By default, the output of this protocol will be a tilt '
                           'series that will have associated the alignment information '
                           'as a transformation matrix. When this option is set as Yes, '
                           'then a second output, called interpolated tilt series, '
                           'is generated. The interpolated tilt series should be used '
                           'for visualization purpose but not for image processing.')

        form.addParam('binning',
                      params.IntParam,
                      condition='computeAlignment',
                      default=1,
                      label='Binning for the interpolated',
                      help='Binning to be applied to the interpolated tilt-series '
                           'in IMOD convention. \n'
                           'Binning is an scaling factor given by an integer greater '
                           'than 1. IMOD uses ordinary binning to reduce images in '
                           'size by the given factor. The value of a binned pixel is '
                           'the average of pixel values in each block of pixels '
                           'being binned. Binning is applied before all other image '
                           'transformations.')

        form.addParam('Trimming parameters', params.LabelParam,
                      label='Tilt axis angle detected from import. In case another '
                            'value is desired please adjust the number below.')

        form.addParam('tiltAxisAngle',
                      params.FloatParam,
                      allowsNull=True,
                      label='Tilt axis angle (degrees)',
                      help='The tilt axis angle is the tilt axis rotation relative '
                           'to the Y axis of the image. If it was not properly set '
                           'in the import of the tilt series, or the imported'
                           'information is not correct you have the chance to '
                           'correct at in this point. Usually, it will be 90 '
                           'degrees less than the RotationAngle in a system with '
                           'no axis inversions.')

        trimming = form.addGroup('Trimming parameters',
                                 expertLevel=params.LEVEL_ADVANCED)

        self.trimingForm(trimming, pxTrimCondition='False',
                         correlationCondition='True',
                         levelType=params.LEVEL_ADVANCED)
        self.filteringParametersForm(form, condition='True',
                                     levelType=params.LEVEL_ADVANCED)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        binning = self.binning.get()
        for tsId in self.tsDict.keys():
            self._insertFunctionStep(self.convertInputStep, tsId)
            self._insertFunctionStep(self.computeXcorrStep, tsId)
            self._insertFunctionStep(self.generateOutputStackStep, tsId)
            if self.computeAlignment:
                self._insertFunctionStep(self.computeInterpolatedStackStep,
                                         tsId, binning)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.tsDict = {ts.getTsId(): ts.clone() for ts in self.getInputSet()}

    def convertInputStep(self, tsId, **kwargs):
        oddEvenFlag = self.applyToOddEven(self.getInputSet())
        super().convertInputStep(tsId, oddEven=oddEvenFlag)

    def computeXcorrStep(self, tsId):
        """Compute transformation matrix for each tilt series. """
        try:
            ts = self.tsDict[tsId]
            tiltAxisAngle = self.getTiltAxisOrientation(ts)

            paramsXcorr = {
                "-input": self.getTmpOutFile(tsId),
                "-output": self.getExtraOutFile(tsId, ext=PREXF_EXT),
                "-tiltfile": self.getExtraOutFile(tsId, ext=TLT_EXT),
                "-RotationAngle": tiltAxisAngle,
                "-FilterSigma1": self.filterSigma1.get(),
                "-FilterSigma2": self.filterSigma2.get(),
                "-FilterRadius1": self.filterRadius1.get(),
                "-FilterRadius2": self.filterRadius2.get()
            }

            if self.cumulativeCorr:
                paramsXcorr["-CumulativeCorrelation"] = ""

            xdim, ydim, _ = ts.getDim()
            xmin, xmax = self.xmin.get() or 0
            xmax = self.xmax.get() or xdim-1
            ymin = self.ymin.get() or 0
            ymax = self.ymax.get() or ydim-1

            paramsXcorr["-xminmax"] = f"{xmin},{xmax}"
            paramsXcorr["-yminmax"] = f"{ymin},{ymax}"

            # Excluded views
            excludedViews = ts.getExcludedViewsIndex(caster=str)
            if len(excludedViews):
                paramsXcorr["-SkipViews"] = ",".join(excludedViews)

            self.runProgram('tiltxcorr', paramsXcorr)

            paramsXftoxg = {
                "-input": self.getExtraOutFile(tsId, ext=PREXF_EXT),
                "-goutput": self.getExtraOutFile(tsId, ext=PREXG_EXT),
                "-NumberToFit": 0
            }
            self.runProgram('xftoxg', paramsXftoxg)

        except Exception as e:
            self._failedTs.append(tsId)
            self.error(f'tiltxcorr or xftoxg execution failed for tsId {tsId} -> {e}')

    def generateOutputStackStep(self, tsId):
        """ Generate tilt-series with the associated transform matrix """
        ts = self.tsDict[tsId]
        if tsId in self._failedTs:
            self.createOutputFailedSet(ts)
        else:
            outputFn = self.getExtraOutFile(tsId, ext=PREXG_EXT)
            if os.path.exists(outputFn):
                output = self.getOutputSetOfTS(self.getInputSet())
                alignmentMatrix = utils.formatTransformationMatrix(outputFn)

                newTs = TiltSeries(tsId=tsId)
                newTs.copyInfo(ts)
                newTs.getAcquisition().setTiltAxisAngle(self.getTiltAxisOrientation(ts))
                output.append(newTs)

                for index, tiltImage in enumerate(ts):
                    newTi = TiltImage()
                    newTi.copyInfo(tiltImage, copyId=True, copyTM=False)

                    transform = Transform()

                    if tiltImage.hasTransform():
                        previousTransform = tiltImage.getTransform().getMatrix()
                        newTransform = alignmentMatrix[:, :, index]
                        previousTransformArray = np.array(previousTransform)
                        newTransformArray = np.array(newTransform)
                        outputTransformMatrix = np.matmul(newTransformArray, previousTransformArray)
                        transform.setMatrix(outputTransformMatrix)
                        newTi.setTransform(transform)
                    else:
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
                self._store(output)

    def computeInterpolatedStackStep(self, tsId, binning):
        ts = self.tsDict[tsId]
        if tsId in self._failedTs:
            self.createOutputFailedSet(ts)
        else:
            xfFile = self.getExtraOutFile(tsId, ext=PREXG_EXT)
            if os.path.exists(xfFile):
                output = self.getOutputInterpolatedTS(self.getInputSet(), binning)
                firstItem = ts.getFirstItem()
                params = self.getBasicNewstackParams(ts,
                                                     self.getExtraOutFile(tsId),
                                                     inputTsFileName=self.getTmpOutFile(tsId),
                                                     xfFile=xfFile,
                                                     firstItem=firstItem,
                                                     binning=binning,
                                                     doNorm=True)
                self.runProgram('newstack', params)

                newTs = TiltSeries(tsId=tsId)
                newTs.copyInfo(ts)
                newTs.getAcquisition().setTiltAxisAngle(self.getTiltAxisOrientation(ts))
                newTs.setInterpolated(True)
                output.append(newTs)

                if binning > 1:
                    newTs.setSamplingRate(ts.getSamplingRate() * binning)

                for index, tiltImage in enumerate(ts):
                    newTi = TiltImage()
                    newTi.copyInfo(tiltImage, copyId=True)
                    newTi.setLocation(index + 1, self.getExtraOutFile(tsId))
                    if binning > 1:
                        newTi.setSamplingRate(tiltImage.getSamplingRate() * binning)
                    newTs.append(newTi)

                dims = self._getOutputDim(self.getExtraOutFile(tsId))
                newTs.setDim(dims)
                newTs.write(properties=False)

                output.update(newTs)
                output.write()
                self._store(output)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append(f"Input tilt-series: {self.getInputSet().getSize()}\n"
                           "Transformation matrices calculated: "
                           f"{self.TiltSeries.getSize()}")
            if self.InterpolatedTiltSeries:
                summary.append("Interpolated tilt-series: "
                               f"{self.InterpolatedTiltSeries.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.TiltSeries:
            methods.append("The transformation matrix has been calculated for "
                           f"{self.TiltSeries.getSize()} tilt-series using "
                           "the IMOD *tiltxcorr* command.")
        return methods

    # --------------------------- UTILS functions ------------------------------
    def getTiltAxisOrientation(self, ts):
        if self.tiltAxisAngle.hasValue():
            return self.tiltAxisAngle.get()
        else:
            return ts.getAcquisition().getTiltAxisAngle()
