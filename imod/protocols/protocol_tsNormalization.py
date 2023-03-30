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
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
from pwem.emlib.image import ImageHandler
import tomo.objects as tomoObj

from .. import Plugin
from .protocol_base import ProtImodBase
from ..utils import formatTransformFile


class ProtImodTSNormalization(ProtImodBase):
    """
    Normalize input tilt-series and change its storing formatting.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/newstack.html
    """

    _label = 'Tilt-series preprocess'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        form.addParam('binning',
                      params.FloatParam,
                      default=1,
                      label='Binning',
                      important=True,
                      help='Binning to be applied to the normalized tilt-series '
                           'in IMOD convention. Images will be binned by the '
                           'given factor. Must be an integer bigger than 1')

        form.addParam('applyAlignment',
                      params.BooleanParam,
                      default=False,
                      label='Apply transformation matrix',
                      help='Apply the tilt series transformation matrix if tilt series have them')

        form.addParam('floatDensities',
                      params.EnumParam,
                      choices=['default', '1', '2', '3', '4'],
                      default=0,
                      label='Adjust densities mode',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Adjust densities of sections individually:\n'
                           '-Default: no adjustment performed\n'
                           '-Mode 1: sections fill the data range\n'
                           '-Mode 2: sections scaled to common mean and standard deviation.\n'
                           '-Mode 3: sections shifted to a common mean without scaling\n'
                           '-Mode 4: sections shifted to a common mean and then '
                           'rescale the resulting minimum and maximum densities '
                           'to the Min and Max values specified')

        form.addParam('modeToOutput',
                      params.EnumParam,
                      choices=['default', '4-bit', 'byte', 'signed 16-bit',
                               'unsigned 16-bit', '32-bit float'],
                      default=0,
                      label='Storage data type',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Apply one density scaling to all sections to '
                           'map current min and max to the given Min and '
                           'Max. The storage mode of the output file. The '
                           'default is the mode of the first input file, '
                           'except for a 4-bit input file, where the default '
                           'is to output as bytes')

        form.addParam('scaleRangeToggle',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      condition="floatDensities==0 or floatDensities==1 or floatDensities==3",
                      default=1,
                      label='Set scaling range values?',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='This option will rescale the densities of all '
                           'sections by the same factors so that the original '
                           'minimum and maximum density will be mapped '
                           'to the Min and Max values that are entered')

        form.addParam('scaleRangeMax',
                      params.FloatParam,
                      condition="(floatDensities==0 or floatDensities==1 or floatDensities==3) and scaleRangeToggle==0",
                      default=255,
                      label='Max.',
                      help='Maximum value for the rescaling')

        form.addParam('scaleRangeMin',
                      params.FloatParam,
                      condition="(floatDensities==0 or floatDensities==1 or floatDensities==3) and scaleRangeToggle==0",
                      default=0,
                      label='Min.',
                      help='Minimum value for the rescaling')

        form.addParam('antialias',
                      params.EnumParam,
                      choices=['None', 'Blackman', 'Triangle', 'Mitchell',
                               'Lanczos 2', 'Lanczos 3'],
                      default=5,
                      label='Antialias method:',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Type of antialiasing filter to use when reducing images.\n'
                           'The available types of filters are:\n\n'
                           'None\n'
                           'Blackman - fast but not as good at antialiasing as slower filters\n'
                           'Triangle - fast but smooths more than Blackman\n'
                           'Mitchell - good at antialiasing, smooths a bit\n'
                           'Lanczos 2 lobes - good at antialiasing, less smoothing than Mitchell\n'
                           'Lanczos 3 lobes - slower, even less smoothing but more risk of ringing\n'
                           'The default is Lanczos 3 as of IMOD 4.7. Although '
                           'many people consider Lanczos 2 the best compromise '
                           'among the various factors, that sentiment may be '
                           'based on images of natural scenes where there are '
                           'sharp edges.')

        groupMeanSd = form.addGroup('Mean and SD',
                                    condition='floatDensities==2',
                                    help='Scale all images to the given mean '
                                         'and standard deviation. This option '
                                         'implies -float 2 and is incompatible '
                                         'with all other scaling options. If no '
                                         'values are set, mean=0 and SD=1 by default')

        groupMeanSd.addParam('meanSdToggle',
                             params.EnumParam,
                             choices=['Yes', 'No'],
                             default=1,
                             label='Set mean and SD?',
                             display=params.EnumParam.DISPLAY_HLIST,
                             help='Set mean and SD values')

        groupMeanSd.addParam('scaleMean',
                             params.FloatParam,
                             default=0,
                             label='Mean',
                             help='Mean value for the rescaling')

        groupMeanSd.addParam('scaleSd',
                             params.FloatParam,
                             default=1,
                             label='SD',
                             help='Standard deviation value for the rescaling')

        groupScale = form.addGroup('Scaling values',
                                   condition='floatDensities==4')

        groupScale.addParam('scaleMax',
                            params.FloatParam,
                            default=255,
                            label='Max.',
                            help='Maximum value for the rescaling')

        groupScale.addParam('scaleMin',
                            params.FloatParam,
                            default=0,
                            label='Min.',
                            help='Minimum value for the rescaling')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep(self.convertInputStep, ts.getObjId())
            self._insertFunctionStep(self.generateOutputStackStep, ts.getObjId())
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self, tsObjId):
        # Interpolation will be done in the generateOutputStep
        super().convertInputStep(tsObjId, imodInterpolation=None,
                                 generateAngleFile=False)

    def generateOutputStackStep(self, tsObjId):
        output = self.getOutputSetOfTiltSeries(self.inputSetOfTiltSeries.get(),
                                               self.binning.get())

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        firstItem = ts.getFirstItem()

        xfFile = None

        if self.applyAlignment.get() and ts.hasAlignment():
            xfFile = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".xf"))
            formatTransformFile(ts, xfFile)

        binning = int(self.binning.get())

        argsNewstack, paramsNewstack = self.getBasicNewstackParams(ts,
                                                                   os.path.join(extraPrefix, firstItem.parseFileName()),
                                                                   inputTsFileName=os.path.join(tmpPrefix, firstItem.parseFileName()),
                                                                   xfFile=xfFile,
                                                                   firstItem=firstItem,
                                                                   binning=binning,
                                                                   )
        paramsNewstack.update({
            'bin': binning,
            'imagebinned': 1.0,
            'antialias': self.antialias.get() + 1
        })

        argsNewstack += "-bin %(bin)d " \
                        "-antialias %(antialias)d " \
                        "-imagebinned %(imagebinned)s "

        if self.floatDensities.get() != 0:
            argsNewstack += " -FloatDensities " + str(self.floatDensities.get())

            if self.floatDensities.get() == 2:
                if self.meanSdToggle.get() == 0:
                    argsNewstack += " -MeanAndStandardDeviation " + str(self.scaleMean.get()) + "," + \
                                    str(self.scaleSd.get())

            elif self.floatDensities.get() == 4:
                argsNewstack += " -ScaleMinAndMax " + str(self.scaleMax.get()) + "," + str(self.scaleMin.get())

            else:
                if self.scaleRangeToggle.get() == 0:
                    argsNewstack += " -ScaleMinAndMax " + str(self.scaleRangeMax.get()) + "," + \
                                    str(self.scaleRangeMin.get())

        if self.getModeToOutput() is not None:
            argsNewstack += " -ModeToOutput " + str(self.getModeToOutput())

        Plugin.runImod(self, 'newstack', argsNewstack % paramsNewstack)

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        output.append(newTs)

        if binning > 1:
            newTs.setSamplingRate(ts.getSamplingRate() * int(self.binning.get()))

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True, copyTM=True)

            # Tranformation matrix
            if tiltImage.hasTransform() and not self.applyAlignment.get():
                newTi = self.updateTM(newTi)
            else:
                newTi.setTransform(None)

            newTi.setAcquisition(tiltImage.getAcquisition())
            newTi.setLocation(index + 1,
                              (os.path.join(extraPrefix, tiltImage.parseFileName())))
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

    def closeOutputSetsStep(self):
        self.TiltSeries.setStreamState(Set.STREAM_CLOSED)
        self.TiltSeries.write()
        self._store()

    # --------------------------- UTILS functions -----------------------------
    def getModeToOutput(self):
        parseParamsOutputMode = {
            0: None,
            1: 101,
            2: 0,
            3: 1,
            4: 6,
            5: 2
        }
        return parseParamsOutputMode[self.modeToOutput.get()]

    def updateTM(self, newTi):
        transform = newTi.getTransform()
        matrix = transform.getMatrix()

        matrix[0][2] /= self.binning.get()
        matrix[1][2] /= self.binning.get()

        transform.setMatrix(matrix)
        newTi.setTransform(transform)
        
        return newTi

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        hasAlign = self.inputSetOfTiltSeries.get().getFirstItem().getFirstItem().hasTransform()
        if self.applyAlignment.get() and not hasAlign:
            errors.append("Input tilt-series do not have alignment information")

        return errors

    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append("Input tilt-series: %d\nInterpolations applied: %d"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.TiltSeries.getSize()))
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.TiltSeries:
            methods.append("%d tilt-series have been normalized using the IMOD "
                           "*newstack* command.\n"
                           % (self.TiltSeries.getSize()))
        return methods
