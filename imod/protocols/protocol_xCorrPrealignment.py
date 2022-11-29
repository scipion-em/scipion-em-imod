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


class ProtImodXcorrPrealignment(ProtImodBase):
    """
    Tilt-series cross correlation alignment based on the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/tiltxcorr.html
    """

    _label = 'Coarse prealignment'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        form.addParam('cumulativeCorr',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Use cumulative correlation?',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Use this option to add up previously aligned '
                           'pictures to get the reference for the next '
                           'alignment. Alignments will start at low tilt '
                           'and work up to high tilt.')

        form.addParam('computeAlignment',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Generate interpolated tilt-series?',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Generate and save the interpolated tilt-series '
                           'applying the obtained transformation matrices.')

        group = form.addGroup('Interpolated tilt-series',
                              condition='computeAlignment==0')

        group.addParam('binning',
                       params.FloatParam,
                       default=1.0,
                       label='Binning',
                       help='Binning to be applied to the interpolated '
                            'tilt-series in IMOD convention. Images will be '
                            'binned by the given factor. Must be an integer '
                            'bigger than 1')

        form.addParam('filterRadius1',
                      params.FloatParam,
                      label='Filter radius 1',
                      default='0.0',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Low spatial frequencies in the cross-correlation "
                           "will be attenuated by a Gaussian curve that is 1 "
                           "at this cutoff radius and falls off below this "
                           "radius with a standard deviation specified by "
                           "FilterSigma2. Spatial frequency units range from "
                           "0 to 0.5. Use FilterSigma1 instead of this entry "
                           "for more predictable attenuation of low frequencies.")

        form.addParam('filterRadius2',
                      params.FloatParam,
                      label='Filter radius 2',
                      default='0.25',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="High spatial frequencies in the cross-correlation "
                           "will be attenuated by a Gaussian curve that is 1 "
                           "at this cutoff radius and falls off above this "
                           "radius with a standard deviation specified by "
                           "FilterSigma2.")

        form.addParam('filterSigma1',
                      params.FloatParam,
                      label='Filter sigma 1',
                      default='0.03',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Sigma value to filter low frequencies in the "
                           "correlations with a curve that is an inverted "
                           "Gaussian.  This filter is 0 at 0 frequency and "
                           "decays up to 1 with the given sigma value. "
                           "However, if a negative value of radius1 is entered, "
                           "this filter will be zero from 0 to "
                           "|radius1| then decay up to 1.")

        form.addParam('filterSigma2',
                      params.FloatParam,
                      label='Filter sigma 2',
                      default='0.05',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Sigma value for the Gaussian rolloff below and "
                           "above the cutoff frequencies specified by "
                           "FilterRadius1 and FilterRadius2")

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

        paramsXcorr = {
            'input': os.path.join(tmpPrefix,
                                  ts.getFirstItem().parseFileName()),
            'output': os.path.join(extraPrefix,
                                   ts.getFirstItem().parseFileName(extension=".prexf")),
            'tiltfile': os.path.join(tmpPrefix,
                                     ts.getFirstItem().parseFileName(extension=".tlt")),
            'rotationAngle': ts.getAcquisition().getTiltAxisAngle(),
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
            'bin': int(self.binning.get()),
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
        newTs.setInterpolated(True)
        output.append(newTs)

        if self.binning > 1:
            newTs.setSamplingRate(ts.getSamplingRate() * int(self.binning.get()))

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(index + 1,
                              (os.path.join(extraPrefix,
                                            tiltImage.parseFileName())))
            if self.binning > 1:
                newTi.setSamplingRate(tiltImage.getSamplingRate() * int(self.binning.get()))
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
