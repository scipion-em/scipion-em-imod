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

import pyworkflow.protocol.params as params
from imod.protocols.protocol_base import IN_TS_SET, PROCESS_ODD_EVEN
from pwem.emlib.image import ImageHandler as ih
from pyworkflow.utils import Message
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries

from imod.protocols import ProtImodBase
from imod.constants import OUTPUT_TILTSERIES_NAME, XF_EXT, ODD, EVEN
from imod.utils import genXfFile


class ProtImodTsNormalization(ProtImodBase):
    """
    Normalize input tilt-series and change its storing formatting.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/newstack.html

    IMOD tilt series preprocess makes use of the Newstack command.
    In particular, three functionalities are possible:\n

    _1 Binning_: The protocol also allows to bin tilt series. This
    means to reduce the dimensions of the tilt series keeping but
    keeping most of the information. The binning factor or simply
    binning is an integer number and represent the scaling factor
    of the images. Binning 2 means that the original images will
    be twice the binned ones.
    _2 Normalization_: This protocol allows to scale the gray values
    of the images, also called normalization, to a common range or
    mean of density. The most used normalization consists in zero
    mean and standard deviation one.\n

    _3 storage format_: IMOD is able to modify the number of bit of
    the stored data in order to reduce the disc occupancy.

    """

    _label = 'Tilt-series preprocess'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        form.addParam('binning',
                      params.IntParam,
                      default=1,
                      label='Binning',
                      important=True,
                      help='Binning is an scaling factor for the output images. '
                           'Must be an integer greater than 1. IMOD uses ordinary'
                           'binning to reduce images in size by the given factor. '
                           'The value of a binned pixel is the average of pixel '
                           'values in each block of pixels being binned. Binning '
                           'is applied before all')

        # TODO: decide if this option should be removed. If not, the out TS must be renamed to interpolated and the non-
        # interpolated should be also generated
        # form.addParam('applyAlignment',
        #               params.BooleanParam,
        #               default=False,
        #               label='Apply transformation matrix',
        #               help='Apply the tilt series transformation matrix if tilt '
        #                    'series have them')

        form.addParam('floatDensities',
                      params.EnumParam,
                      choices=['No adjust',
                               'adjust each section to fill the data range',
                               'scaled to common mean and standard deviation',
                               'shifted to a common mean without scaling',
                               'shifted to mean and rescaled to a min and max'],
                      default=2,
                      label='Adjust densities mode',
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='Adjust densities of sections individually. Modes:\n\n'
                           '*0 - No adjustment performed*\n\n'
                           
                           '*1 - Adjust each section to fill the data range*\n\n'

                           # '-Range between min and max: This option will scale the gray values'
                           # 'to be in a range given by a minimum and a maximum values.'
                           # 'This is the mode 1 in newstack flag -floatDensities.\n'

                           '*2 - Scaled to common mean and standard deviation*:\n\n'
                           'This is the most '
                           'common normalization procedure. The new tilt series will have'
                           'a meand and a standard deviation introduced by the user. Generaly,'
                           'a zero meand and a standard deviation one is a good choice.'
                           'This is the mode 2 in newstack flag -floatDensities.\n\n'

                           '*3 - Shifted to a common mean without scaling*:\n\n'
                           'This option only'
                           'add an offset to the gray values of the images. The offset will'
                           'be calculated such as the new images will present a mean gray value'
                           'introduced by the user. This is the mode 3 in newstack flag '
                           'floatDensities.\n\n'

                           '*4 - shifted to mean and rescaled to a min and max*:\n\nIn this case, an '
                           'offset is added to the images in order to achieve a mean gray value '
                           'then they are rescale the resulting minimum and maximum densities '
                           'to the Min and Max values specified. This is the mode 4 in newstack '
                           'flag -floatDensities.')
        # NEWSTACK - The -scale, -contrast, -multadd, and -float options are mutually exclusive except with -float 4
        # scaleRangeToggleCond = "floatDensities in [0, 4]"
        # -meansd (-mea) OR -MeanAndStandardDeviation   Two floats
        #               Scale all images to the given mean and standard deviation.  This
        #               option implies -float 2 and is incompatible with all other scaling options.
        groupMeanSd = form.addGroup('Mean and SD',
                                    condition='floatDensities==2',
                                    help='Scale all images to the given mean '
                                         'and standard deviation.')
        floatDensMode2Cond = 'floatDensities == 2'
        groupMeanSd.addParam('scaleMean',
                             params.FloatParam,
                             condition=floatDensMode2Cond,
                             default=0,
                             label='Mean',
                             help='Mean value for the rescaling')

        groupMeanSd.addParam('scaleSd',
                             params.FloatParam,
                             condition=floatDensMode2Cond,
                             default=1,
                             label='SD',
                             help='Standard deviation value for the rescaling')

        groupScale = form.addGroup('Scaling values',
                                   condition='floatDensities in [0, 4]')
        msg = 'This option will rescale the densities of all sections by the same factors so that the original ' \
              'minimum and maximum density will be mapped to the Min and Max values that are entered.'
        groupScale.addParam('scaleMax',
                            params.FloatParam,
                            default=255.,
                            label='Max.',
                            help=f'Maximum value for the rescaling. {msg}')

        groupScale.addParam('scaleMin',
                            params.FloatParam,
                            default=0.,
                            label='Min.',
                            help=f'Minimum value for the rescaling. {msg}')

        form.addParam('modeToOutput',
                      params.EnumParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      choices=['default', '4-bit', '8-bit', 'signed 16-bit',
                               'unsigned 16-bit', '32-bit float'],
                      default=0,
                      label='Storage data type',
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='The storage mode of the output file. The '
                           'default is the mode of the first input file, '
                           'except for a 4-bit input file, where the default '
                           'is to output as bytes')

        form.addParam('antialias',
                      params.EnumParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      choices=['None', 'Blackman', 'Triangle', 'Mitchell',
                               'Lanczos 2 lobes', 'Lanczos 3 lobes'],
                      default=5,
                      label='Antialias method:',
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='Type of antialiasing filter to use when reducing images.\n'
                           'The available types of filters are:\n\n'
                           'None - Antialias will not be applied\n'
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

        form.addParam(PROCESS_ODD_EVEN,
                      params.BooleanParam,
                      default=False,
                      label='Apply to odd/even',
                      help='If True, the full tilt series and the associated odd/even '
                           'tilt series will be processed. The transformations applied '
                           'to the odd/even tilt series will be exactly the same.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        binning = self.binning.get()
        for tsId in self.tsDict.keys():
            self._insertFunctionStep(self.convertInputStep, tsId)
            self._insertFunctionStep(self.generateOutputStackStep, tsId, binning)
            self._insertFunctionStep(self.createOutputStep, tsId, binning)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self.getInputSet()}
        self.oddEvenFlag = self.applyToOddEven(self.getInputSet())

    def convertInputStep(self, tsId, **kwargs):
        # Interpolation will be done in the generateOutputStep
        super().convertInputStep(tsId,
                                 imodInterpolation=None,
                                 generateAngleFile=False,
                                 oddEven=self.oddEvenFlag)

    def generateOutputStackStep(self, tsId, binning):
        try:
            ts = self.tsDict[tsId]
            firstItem = ts.getFirstItem()
            xfFile = None
            norm = self.floatDensities.get()
            paramsDict = self.getBasicNewstackParams(ts,
                                                     self.getExtraOutFile(tsId),
                                                     inputTsFileName=self.getTmpOutFile(tsId),
                                                     xfFile=xfFile,
                                                     firstItem=firstItem,
                                                     binning=binning,
                                                     doNorm=norm != 0)
            paramsDict["-antialias"] = self.antialias.get() + 1
            # Float densities
            if norm > 0:
                paramsDict["-FloatDensities"] = norm
                if norm == 2:
                    paramsDict["-MeanAndStandardDeviation"] = f"{self.scaleMean.get()},{self.scaleSd.get()}"
                elif norm == 4:
                    paramsDict["-ScaleMinAndMax"] = f"{self.scaleMax.get()},{self.scaleMin.get()}"

            if self.getModeToOutput() is not None:
                paramsDict["-ModeToOutput"] = self.getModeToOutput()

            self.runProgram("newstack", paramsDict)

            if self.oddEvenFlag:
                paramsDict['-input'] = self.getTmpOutFile(tsId, suffix=ODD)
                paramsDict['-output'] = self.getExtraOutFile(tsId, suffix=ODD)
                self.runProgram("newstack", paramsDict)

                paramsDict['-input'] = self.getTmpOutFile(tsId, suffix=EVEN)
                paramsDict['-output'] = self.getExtraOutFile(tsId, suffix=EVEN)
                self.runProgram("newstack", paramsDict)

        except Exception as e:
            self._failedTs.append(tsId)
            self.error(f'newstack execution failed for tsId {tsId} -> {e}')

    def createOutputStep(self, tsId, binning):
        ts = self.tsDict[tsId]
        if tsId in self._failedTs:
            self.createOutputFailedSet(ts)
        else:
            outputFn = self.getExtraOutFile(tsId)
            if os.path.exists(outputFn):
                output = self.getOutputSetOfTS(self.getInputSet(), binning)
                newTs = TiltSeries(tsId=tsId)
                newTs.copyInfo(ts)
                output.append(newTs)

                if binning > 1:
                    newTs.setSamplingRate(ts.getSamplingRate() * binning)

                for index, tiltImage in enumerate(ts):
                    newTi = TiltImage()
                    newTi.copyInfo(tiltImage, copyId=True, copyTM=True)

                    # Tranformation matrix
                    if tiltImage.hasTransform():
                        newTi = self.updateTM(newTi, binning)

                    newTi.setAcquisition(tiltImage.getAcquisition())
                    if self.oddEvenFlag:
                        locationOdd = index + 1, self.getExtraOutFile(tsId, suffix=ODD)
                        locationEven = index + 1, self.getExtraOutFile(tsId, suffix=EVEN)
                        newTi.setOddEven([ih.locationToXmipp(locationOdd), ih.locationToXmipp(locationEven)])
                    else:
                        newTi.setOddEven([])

                    newTi.setLocation(index + 1, outputFn)

                    if binning > 1:
                        newTi.setSamplingRate(tiltImage.getSamplingRate() * binning)
                    newTs.append(newTi)

                dims = self._getOutputDim(outputFn)
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
                           f"Interpolations applied: {self.TiltSeries.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.TiltSeries:
            methods.append(f"{self.TiltSeries.getSize()} tilt-series have been "
                           "normalized using the IMOD *newstack* command.")
        return methods

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

    @staticmethod
    def updateTM(newTi, binning):
        transform = newTi.getTransform()
        matrix = transform.getMatrix()

        matrix[0][2] /= binning
        matrix[1][2] /= binning

        transform.setMatrix(matrix)
        newTi.setTransform(transform)

        return newTi
