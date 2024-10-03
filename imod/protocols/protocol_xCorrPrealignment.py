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
from imod.protocols.protocol_base import IN_TS_SET, BINNING_FACTOR
from pwem import ALIGN_NONE
from pwem.objects import Transform
from pyworkflow.utils import Message
from tomo.objects import SetOfTiltSeries

from imod import utils
from imod.protocols import ProtImodBase
from imod.constants import (TLT_EXT, PREXF_EXT, PREXG_EXT,
                            OUTPUT_TILTSERIES_NAME,
                            OUTPUT_TS_INTERPOLATED_NAME)


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
        form.addSection(Message.LABEL_INPUT)

        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-series to be prealigned')

        form.addParam('cumulativeCorr',
                      params.BooleanParam,
                      default=False,
                      label='Use cumulative correlation?',
                      help='The program will take the image at zero tilt as the first'
                           'reference, and correlate it with the image at the next '
                           'most negative tilt. It will then add the aligned image '
                           'to the first reference to make the reference for the next '
                           'tilt. At each tilt, the reference will be the sum of '
                           'images that have already been aligned. When the most '
                           'negative tilt angle is reached, the procedure is repeated '
                           'from the zero-tilt view to more positive tilt angles.')

        form.addParam('computeAlignment',
                      params.BooleanParam,
                      default=False,
                      label='Generate interpolated tilt-series?',
                      important=True,
                      help='Generate and save the interpolated tilt-series applying '
                           'the obtained transformation matrices.\n'
                           'By default, the output of this protocol will be a tilt '
                           'series that will have associated the alignment information '
                           'as a transformation matrix. When this option is set as Yes, '
                           'then a second output, called interpolated tilt series, '
                           'is generated. The interpolated tilt series should be used '
                           'for visualization purpose but not for image processing.')

        form.addParam(BINNING_FACTOR,
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

        self.addTrimingParams(trimming,
                              pxTrimCondition=False,
                              correlationCondition=True,
                              levelType=params.LEVEL_ADVANCED)
        self.filteringParametersForm(form,
                                     condition=True,
                                     levelType=params.LEVEL_ADVANCED)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        binning = self.binning.get()
        closeSetStepDeps = []
        for tsId in self.tsDict.keys():
            convId = self._insertFunctionStep(self.convertInputStep, tsId, prerequisites=[])
            compId = self._insertFunctionStep(self.computeXcorrStep, tsId, prerequisites=[convId])
            outId = self._insertFunctionStep(self.generateOutputStackStep, tsId, prerequisites=[compId])
            closeSetStepDeps.append(outId)
            if self.computeAlignment:
                intpId = self._insertFunctionStep(self.computeInterpolatedStackStep, tsId, binning,
                                                  prerequisites=[outId])
                closeSetStepDeps.append(intpId)

        self._insertFunctionStep(self.closeOutputSetsStep, prerequisites=closeSetStepDeps)

    # --------------------------- STEPS functions -----------------------------
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

            doTrim = any([getattr(self, attr).hasValue() for
                          attr in ["xmin", "xmax", "ymin", "ymax"]])
            if doTrim:
                xdim, ydim, _ = ts.getDim()
                xmin, xmax = self.xmin.get() or 0, self.xmax.get() or xdim-1
                ymin, ymax = self.ymin.get() or 0, self.ymax.get() or ydim-1

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
            self._failedItems.append(tsId)
            self.error(f'tiltxcorr or xftoxg execution failed for tsId {tsId} -> {e}')

    def generateOutputStackStep(self, tsId):
        """ Generate tilt-series with the associated transform matrix """
        ts = self.tsDict[tsId]
        if tsId in self._failedItems:
            self.createOutputFailedSet(ts)
        else:
            outputFn = self.getExtraOutFile(tsId, ext=PREXG_EXT)
            if os.path.exists(outputFn):
                tAx = self.tiltAxisAngle.get()
                output = self.getOutputSetOfTS(self.getInputSet(pointer=True),
                                               tiltAxisAngle=tAx)
                alignmentMatrix = utils.formatTransformationMatrix(outputFn)
                self.copyTsItems(output, ts, tsId,
                                 updateTsCallback=self.updateTsNonInterp,
                                 updateTiCallback=self.updateTiNonInterp,
                                 copyDisabledViews=True,
                                 copyId=True,
                                 copyTM=False,
                                 alignmentMatrix=alignmentMatrix,
                                 tiltAxisAngle=tAx)
            else:
                self.createOutputFailedSet(ts)

    def computeInterpolatedStackStep(self, tsId, binning):
        if tsId not in self._failedItems:
            ts = self.tsDict[tsId]
            xfFile = self.getExtraOutFile(tsId, ext=PREXG_EXT)
            if os.path.exists(xfFile):
                output = self.getOutputSetOfTS(self.getInputSet(pointer=True),
                                               binning,
                                               attrName=OUTPUT_TS_INTERPOLATED_NAME,
                                               suffix="Interpolated")
                firstItem = ts.getFirstItem()
                tsExcludedIndices = ts.getExcludedViewsIndex()
                paramsDict = self.getBasicNewstackParams(ts,
                                                         self.getExtraOutFile(tsId),
                                                         inputTsFileName=self.getTmpOutFile(tsId),
                                                         xfFile=xfFile,
                                                         firstItem=firstItem,
                                                         binning=binning,
                                                         tsExcludedIndices=tsExcludedIndices,
                                                         doNorm=True)
                self.runProgram('newstack', paramsDict)

                self.copyTsItems(output, ts, tsId,
                                 updateTsCallback=self.updateTsInterp,
                                 updateTiCallback=self.updateTi,
                                 copyId=True,
                                 copyTM=False,
                                 excludedViews=len(tsExcludedIndices) > 0)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append(f"Input tilt-series: {self.getInputSet().getSize()}\n"
                           "Transformation matrices calculated: "
                           f"{self.TiltSeries.getSize()}")

            interpTS = getattr(self, OUTPUT_TS_INTERPOLATED_NAME, None)
            if interpTS is not None:
                summary.append("Interpolated tilt-series: "
                               f"{interpTS.getSize()}")
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

    def updateTsNonInterp(self, tsId, ts, tsOut, **kwargs):
        tsOut.getAcquisition().setTiltAxisAngle(self.getTiltAxisOrientation(ts))
        tsOut.setAlignment2D()

    @staticmethod
    def updateTiNonInterp(origIndex, index, tsId, ts, ti, tsOut, tiOut, alignmentMatrix=None,
                          tiltAxisAngle=None, **kwargs):
        transform = Transform()
        newTransform = alignmentMatrix[:, :, index]
        newTransformArray = np.array(newTransform)

        if ti.hasTransform():
            previousTransform = ti.getTransform().getMatrix()
            previousTransformArray = np.array(previousTransform)
            outputTransformMatrix = np.matmul(newTransformArray, previousTransformArray)
            transform.setMatrix(outputTransformMatrix)
        else:
            transform.setMatrix(newTransformArray)

        tiOut.setTransform(transform)
        if tiltAxisAngle:
            tiOut.getAcquisition().setTiltAxisAngle(tiltAxisAngle)

    @staticmethod
    def updateTsInterp(tsId, ts, tsOut, **kwargs):
        tsOut.getAcquisition().setTiltAxisAngle(0.)
        tsOut.setAlignment(ALIGN_NONE)
        tsOut.setInterpolated(True)



