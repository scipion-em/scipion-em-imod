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
from imod.protocols.protocol_base import IN_TS_SET
from pyworkflow.object import String
from pyworkflow.protocol.constants import STEPS_SERIAL
from pyworkflow.utils import Message
from tomo.objects import Tomogram, SetOfTomograms

from imod import Plugin
from imod.protocols import ProtImodBase
from imod.constants import (TLT_EXT, ODD, EVEN, MRC_EXT,
                            OUTPUT_TOMOGRAMS_NAME)


class ProtImodTomoReconstruction(ProtImodBase):
    """
    Tomogram reconstruction procedure based on the IMOD procedure.

    More info:
        https://bio3d.colorado.edu/imod/doc/man/tilt.html

    It makes use of the IMOD tilt program. Tilt is a program for reconstructing
    a three-dimensional object (tomogram) from a series of 2D projections.
    The projections are assumed to arise from rotation about a fixed tilt axis,
    subject to minor variations from this scheme.\n

    The program uses a number of different numerical strategies depending on
    the complexity of the alignments needed to reconstruct the tomogram
    and on whether the computation is being done by the central processing
    unit (CPU) or the graphical processing unit (GPU).  If there are no
    local alignments being applied, then for processing on the CPU, the
    program will do a preliminary stretching of each input line by the
    cosine of the tilt angle.  This stretching speeds up the direct back-
    projection because each stretched input line is in register with the
    lines of the output planes. When computing on the GPU, the program
    does not use cosine stretching, thus avoiding the consequences of
    interpolating the data twice.\n
    """

    _label = 'Tomo reconstruction'
    _possibleOutputs = {OUTPUT_TOMOGRAMS_NAME: SetOfTomograms}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.stepsExecutionMode = STEPS_SERIAL
        self.widthWarnMsg = String()

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)

        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt Series')

        form.addParam('tomoThickness',
                      params.IntParam,
                      default=300,
                      label='Tomogram thickness (voxels)',
                      important=True,
                      help='Size in voxels of the tomogram along the z '
                           'axis (beam direction).')

        form.addParam('tomoWidth',
                      params.IntParam,
                      default=0,
                      label='Tomogram width (voxels)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Number of pixels to cut out in X, centered on the middle in X. '
                           'Leave 0 for default X.')

        lineShift = form.addLine('Tomogram shift (px)',
                                 expertLevel=params.LEVEL_ADVANCED,
                                 help="This entry allows one to shift the reconstructed"
                                      " slice in X or Z before it is output.  If the "
                                      " X shift is positive, the slice will be shifted to "
                                      " the right, and the output will contain the left "
                                      " part of the whole potentially reconstructable area. "
                                      " If the Z shift is positive, the slice is shifted "
                                      " upward. The Z entry is optional and defaults to 0 when "
                                      " omitted.")

        lineShift.addParam('tomoShiftX',
                           params.FloatParam,
                           default=0,
                           label=' in X ')

        lineShift.addParam('tomoShiftZ',
                           params.FloatParam,
                           default=0,
                           label=' in Z ')

        lineoffSet = form.addLine('Offset (deg) of the ',
                                  expertLevel=params.LEVEL_ADVANCED,
                                  help="* Tilt angle offset: pply an angle offset in "
                                       "degrees to all tilt angles. This offset "
                                       "positively rotates the reconstructed sections "
                                       "anticlockwise.\n* Tilt axis offset: Apply an "
                                       "offset to the tilt axis in a stack of full-sized "
                                       "projection images, cutting the X-axis at NX/2 + "
                                       "offset instead of NX/2.")

        lineoffSet.addParam('angleOffset',
                            params.FloatParam,
                            default=0,
                            label='Tilt angles ',
                            help='Apply an angle offset in degrees to all tilt '
                                 'angles. This offset positively rotates the '
                                 'reconstructed sections anticlockwise.')

        lineoffSet.addParam('tiltAxisOffset',
                            params.FloatParam,
                            default=0,
                            label='Tilt axis',
                            help='Apply an offset to the tilt axis in a stack of '
                                 'full-sized projection images, cutting the '
                                 'X-axis at NX/2 + offset instead of NX/2. The '
                                 'DELXX entry is optional and defaults to 0 '
                                 'when omitted.')

        form.addParam('superSampleFactor',
                      params.IntParam,
                      default=2,
                      label='Super-sampling factor',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Compute slices in pixels smaller by this factor to reduce artifacts.'
                           'Super-sampling refers to computing the back projection in a slice '
                           'larger by an integer factor in each dimension, which is done here in '
                           'one of two ways: by interpolating the projection data at smaller '
                           'intervals during backprojection, or by expanding the input lines by '
                           'the factor using sync interpolation (padding in Fourier space).  '
                           'Super-sampling will reduce the rays along the projection angles that '
                           'appear to reflect from the edges of a Fourier transform of an X/Z slice.'
                           'These rays result from back-projecting into discrete pixels and represent'
                           'inappropriate information generated by the transitions between successive '
                           'pixels along a backprojection ray.  Super-sampling by 2 will remove most '
                           'of these rays, especially oblique ones.  The additional benefit (amount '
                           'of change in the image) of going from 2 to 3 is about 10% as large as '
                           'that of super-sampling by 2; it is about 3% going from 3 to 4 and 1.5% '
                           'going from 4 to 8 (the highest allowed value). The super-sampled slice '
                           'is reduced to the original resolution by cropping its Fourier transform. '
                           'Super-sampling alone with the first method does not reduce the rays in '
                           'the corners of the Fourier transform past 0.5 cycle/pixel. These rays '
                           'are  also inappropriate since they originate from  lower-frequency '
                           'information in the projection images, so they are removed from the '
                           'cropped Fourier transform before inverting it. This removal has an '
                           'added benefit about 1/3 as large as the benefit from supersampling by 2. '
                           'The effect of these removals is a subtle change in the noise, and a '
                           'benefit may show up only with subvolume averaging.\n'
                           'The additional effect of expanding the input lines is to avoid attenu-'
                           'ating frequencies past half Nyquist (0.25/pixel) any more than is '
                           'achieved with the falloff of the radial filter.  This will be noticeable '
                           'in a tomogram and would be particularly helpful if setting the radial '
                           'filter cutoff higher than the default for subvolume averaging.\n')

        form.addParam('fakeInteractionsSIRT',
                      params.IntParam,
                      default=0,
                      label='Iterations of a SIRT-like equivalent filter',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Modify the radial filter to produce a '
                           'reconstruction equivalent to the one produced by '
                           'the given number of iterations of SIRT. The '
                           'Gaussian filter is applied at the high-frequency '
                           'end of the filter. The functioning of this filter '
                           'is described in: \n\t'
                           'https://bio3d.colorado.edu/imod/doc/man/tilt.html'
                           'This entry corresponds to the imod parameter – FakeSIRTiterations')

        groupRadialFrequencies = form.addGroup('Radial filtering',
                                               help='This entry controls low-pass filtering with the radial weighting '
                                                    'function. The radial weighting function is linear away from the '
                                                    'origin out to the distance in reciprocal space specified by the '
                                                    'first value, followed by a Gaussian fall-off determined by the '
                                                    'second value.',
                                               expertLevel=params.LEVEL_ADVANCED)

        groupRadialFrequencies.addParam('radialFirstParameter',
                                        params.FloatParam,
                                        default=0.35,
                                        label='Cutoff linear region',
                                        help='Linear region of the filter. This is the first'
                                             'parameter of the -RADIAL flag in IMOD tilt command.'
                                             'If the cutoff is greater than 1 the distances are interpreted '
                                             'as pixels in Fourier space; otherwise they are treated as '
                                             'frequencies in cycles/pixel, which range from 0 to 0.5. '
                                             'Use a cutoff of 0.5 for no low-pass filtering.')

        groupRadialFrequencies.addParam('radialSecondParameter',
                                        params.FloatParam,
                                        default=0.035,
                                        label='Radial fall-off',
                                        help='Gaussian fall-off parameter. The sigma (or standard deviation) '
                                             'of the Gaussian is the second value times 0.707.')

        form.addHidden(params.USE_GPU,
                       params.BooleanParam,
                       default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation. "
                            "Select the one you want to use.")

        form.addHidden(params.GPU_LIST,
                       params.StringParam,
                       default='0',
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="GPU ID. To pick the best available one set 0. "
                            "For a specific GPU set its number ID "
                            "(starting from 1).")

        self.addOddEvenParams(form)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        widthWarnTsIds = []
        self._initialize()
        for tsId, ts in self.tsDict.items():
            xDim = ts.getXDim()
            tomoWidth = self.tomoWidth.get()
            if tomoWidth > xDim:
                tomoWidth = 0
                widthWarnTsIds.append(tsId)
            self._insertFunctionStep(self.convertInputStep, tsId)
            self._insertFunctionStep(self.computeReconstructionStep, tsId, tomoWidth)
            self._insertFunctionStep(self.createOutputStep, tsId)
        self._insertFunctionStep(self.closeOutputSetsStep)

        if widthWarnTsIds:
            self.widthWarnMsg.set(f'\n\n*WARNING!:*'
                                  f'\nThe introduced width is greater than the X dimension of the tomograms: '
                                  f'*{widthWarnTsIds}*.'
                                  f'\nValue 0 was assumed for all of them.')

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self, tsId, **kwargs):
        presentAcqOrders = self.getPresentAcqOrders(self.tsDict[tsId],
                                                    onlyEnabled=True)  # Re-stack excluding views before reconstructing
        super().convertInputStep(tsId,
                                 doSwap=True,
                                 oddEven=self.oddEvenFlag,
                                 presentAcqOrders=presentAcqOrders)

    def computeReconstructionStep(self, tsId, tomoWidth):
        try:
            # run tilt
            paramsTilt = {
                "-InputProjections": self.getTmpOutFile(tsId),
                "-OutputFile": self.getTmpOutFile(tsId, ext=MRC_EXT),
                "-TILTFILE": self.getExtraOutFile(tsId, ext=TLT_EXT),
                "-THICKNESS": self.tomoThickness.get(),
                "-FalloffIsTrueSigma": 1,
                "-RADIAL": f"{self.radialFirstParameter.get()},{self.radialSecondParameter.get()}",
                "-SHIFT": f"{self.tomoShiftX.get()},{self.tomoShiftZ.get()}",
                "-OFFSET": f"{self.angleOffset.get()},{self.tiltAxisOffset.get()}",
                "-MODE": 1,
                "-PERPENDICULAR": "",
                "-AdjustOrigin": "",
                "SuperSampleFactor": self.superSampleFactor.get()
            }

            if self.fakeInteractionsSIRT.get() != 0:
                paramsTilt["-FakeSIRTiterations"] = self.fakeInteractionsSIRT.get()

            # NOTE: the excluded views were  before at newstack level (this is why the lines below are commented)
            # # Excluded views
            # ts = self.tsDict[tsId]
            # excludedViews = ts.getExcludedViewsIndex(caster=str)
            # if len(excludedViews):
            #     paramsTilt["-EXCLUDELIST2"] = ",".join(excludedViews)

            if self.usesGpu():
                paramsTilt["-UseGPU"] = self.getGpuList()[0]
                paramsTilt["-ActionIfGPUFails"] = "2,2"

            self.runProgram('tilt', paramsTilt)

            # run trimvol
            trimVolOpts = "-rx "
            if tomoWidth > 0:
                trimVolOpts += f" -nx {tomoWidth}"

            paramsTrimVol = {
                'input': self.getTmpOutFile(tsId, ext=MRC_EXT),
                'output': self.getExtraOutFile(tsId, ext=MRC_EXT),
                'options': trimVolOpts
            }

            argsTrimvol = "%(options)s %(input)s %(output)s"
            Plugin.runImod(self, 'trimvol', argsTrimvol % paramsTrimVol)

            oddEvenTmp = [[], []]
            if self.oddEvenFlag:
                # Odd
                paramsTilt['-InputProjections'] = self.getTmpOutFile(tsId, suffix=ODD)
                oddEvenTmp[0] = self.getExtraOutFile(tsId, suffix=ODD, ext=MRC_EXT)
                paramsTilt['-OutputFile'] = oddEvenTmp[0]
                self.runProgram('tilt', paramsTilt)

                paramsTrimVol['input'] = oddEvenTmp[0]
                paramsTrimVol['output'] = self.getExtraOutFile(tsId, suffix=ODD, ext=MRC_EXT)
                Plugin.runImod(self, 'trimvol', argsTrimvol % paramsTrimVol)

                # Even
                paramsTilt['-InputProjections'] = self.getTmpOutFile(tsId, suffix=EVEN)
                oddEvenTmp[1] = self.getExtraOutFile(tsId, suffix=EVEN, ext=MRC_EXT)
                paramsTilt['-OutputFile'] = oddEvenTmp[1]
                self.runProgram('tilt', paramsTilt)

                paramsTrimVol['input'] = oddEvenTmp[1]
                paramsTrimVol['output'] = self.getExtraOutFile(tsId, suffix=EVEN, ext=MRC_EXT)
                Plugin.runImod(self, 'trimvol', argsTrimvol % paramsTrimVol)

        except Exception as e:
            self._failedItems.append(tsId)
            self.error(f'tilt or trimvol execution failed for tsId {tsId} -> {e}')

    def createOutputStep(self, tsId):
        ts = self.tsDict[tsId]
        if tsId in self._failedItems:
            self.createOutputFailedSet(ts)
        else:
            tomoLocation = self.getExtraOutFile(tsId, ext=MRC_EXT)
            if os.path.exists(tomoLocation):
                output = self.getOutputSetOfTomograms(self.getInputSet(pointer=True))

                newTomogram = Tomogram(tsId=tsId)
                newTomogram.copyInfo(ts)
                newTomogram.setLocation(tomoLocation)

                if self.oddEvenFlag:
                    halfMapsList = [self.getExtraOutFile(tsId, suffix=ODD, ext=MRC_EXT),
                                    self.getExtraOutFile(tsId, suffix=EVEN, ext=MRC_EXT)]
                    newTomogram.setHalfMaps(halfMapsList)

                # Set default tomogram origin
                newTomogram.setOrigin(newOrigin=None)
                shiftX = self.tomoShiftX.get()
                shiftZ = self.tomoShiftZ.get()
                if shiftX or shiftZ:
                    sRate = newTomogram.getSamplingRate()
                    x, y, z = newTomogram.getShiftsFromOrigin()
                    shiftXang = shiftX * sRate
                    shiftZang = shiftZ * sRate
                    newTomogram.setShiftsInOrigin(x=x - shiftXang, y=y, z=z - shiftZang)

                output.append(newTomogram)
                output.update(newTomogram)
                output.write()
                self._store(output)
            else:
                self.createOutputFailedSet(ts)

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if self.Tomograms:
            summary.append(f"Input tilt-series: {self.getInputSet().getSize()}\n"
                           f"Tomograms reconstructed: {self.Tomograms.getSize()}")
            widthWarnings = self.widthWarnMsg.get()
            if widthWarnings:
                summary.append(widthWarnings)
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.Tomograms:
            methods.append("The reconstruction has been computed for "
                           f"{self.Tomograms.getSize()} tilt-series using "
                           "the IMOD *tilt* command.")
        return methods
