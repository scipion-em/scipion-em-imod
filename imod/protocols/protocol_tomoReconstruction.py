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
from tomo.objects import Tomogram

from .. import Plugin
from .protocol_base import (ProtImodBase, EXT_MRC_ODD_NAME, EXT_MRC_EVEN_NAME,
                            EXT_MRCS_TS_EVEN_NAME, EXT_MRCS_TS_ODD_NAME, TLT_EXT, ODD, MRCS_EXT, EVEN, MRC_EXT, REC_EXT)


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
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt Series')

        form.addParam('tomoThickness',
                      params.FloatParam,
                      default=1000,
                      label='Tomogram thickness (voxels)',
                      important=True,
                      help='Size in voxels of the tomogram along the z '
                           'axis (beam direction).')

        form.addParam('tomoWidth',
                      params.IntParam,
                      default=0,
                      label='Tomogram width (voxels)',
                      help='Number of pixels to cut out in X, centered on the middle in X. Leave 0 for default X.')

        lineShift = form.addLine('Tomogram shift (Å)',
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
                                  help="* Tilt angle offset: pply an angle offset in degrees to all tilt "
                                       "angles. This offset positively rotates the reconstructed sections anticlockwise.\n"
                                       "* Tilt axis offset: Apply an offset to the tilt axis in a stack of full-sized "
                                       "projection images, cutting the X-axis at NX/2 + offset instead of NX/2.")

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

        form.addParam('processOddEven',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=True,
                      label='Reconstruct odd/even?',
                      help='If True, the full tilt series and the associated odd/even tilt series will be reconstructed. '
                           'The alignment applied to the odd/even tilt series will be exactly the same.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        for tsId in self.tsDict.keys():
            self._insertFunctionStep(self.convertInputStep, tsId)
            self._insertFunctionStep(self.computeReconstructionStep, tsId)
            self._insertFunctionStep(self.createOutputStep, tsId)
            self._insertFunctionStep(self.createOutputFailedStep, tsId)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self._failedTs = []
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self.inputSetOfTiltSeries.get()}

    def convertInputStep(self, tsId, **kwargs):
        # Considering swapXY is required to make tilt axis vertical
        oddEvenFlag = self.applyToOddEven(self.inputSetOfTiltSeries.get())
        super().convertInputStep(tsId, doSwap=True, oddEven=oddEvenFlag)

    # @ProtImodBase.tryExceptDecorator
    def computeReconstructionStep(self, tsId):
        ts = self.tsDict[tsId]
        paramsTilt = {
            'InputProjections': self.getTmpOutFile(tsId),
            'OutputFile': self.getTmpOutFile(tsId, ext=REC_EXT),
            'TiltFile': self.getExtraOutFile(tsId, ext=TLT_EXT),
            'Thickness': self.tomoThickness.get(),
            'FalloffIsTrueSigma': 1,
            'Radial': str(self.radialFirstParameter.get()) + "," + str(self.radialSecondParameter.get()),
            'Shift': str(self.tomoShiftX.get()) + "," + str(self.tomoShiftZ.get()),
            'Offset': str(self.angleOffset.get()) + "," + str(self.tiltAxisOffset.get()),
            'SuperSampleFactor': self.superSampleFactor.get(),
        }

        argsTilt = "-InputProjections %(InputProjections)s " \
                   "-OutputFile %(OutputFile)s " \
                   "-TILTFILE %(TiltFile)s " \
                   "-THICKNESS %(Thickness)d " \
                   "-FalloffIsTrueSigma %(FalloffIsTrueSigma)d " \
                   "-RADIAL %(Radial)s " \
                   "-SHIFT %(Shift)s " \
                   "-OFFSET %(Offset)s " \
                   "-MODE 1 " \
                   "-PERPENDICULAR " \
                   "-SuperSampleFactor %(SuperSampleFactor)d " \
                   "-AdjustOrigin "

        if self.fakeInteractionsSIRT.get() != 0:
            paramsTilt.update({
                'FakeSIRTInteractions': self.fakeInteractionsSIRT.get()
            })
            argsTilt += "-FakeSIRTiterations %(FakeSIRTInteractions)d "

        # Excluded views
        excludedViews = ts.getExcludedViewsIndex(caster=str)
        if len(excludedViews):
            argsTilt += f"-EXCLUDELIST2 {','.join(excludedViews)} "

        if self.usesGpu():
            argsTilt += f"-UseGPU {self.getGpuList()[0]} " \
                        "-ActionIfGPUFails 2,2 "

        Plugin.runImod(self, 'tilt', argsTilt % paramsTilt)

        def getArgs():
            args = "-rx "

            if self.tomoWidth.get():
                args += " -nx %s" % self.tomoWidth.get()

            return args

        oddEvenTmp = [[], []]

        if self.applyToOddEven(ts):
            paramsTilt['InputProjections'] = self.getExtraOutFile(tsId, suffix=ODD, ext=MRCS_EXT)
            oddEvenTmp[0] = self.getExtraOutFile(tsId, suffix=ODD, ext=REC_EXT)
            paramsTilt['OutputFile'] = oddEvenTmp[0]
            Plugin.runImod(self, 'tilt', argsTilt % paramsTilt)

            paramsTilt['InputProjections'] = self.getExtraOutFile(tsId, suffix=EVEN, ext=MRCS_EXT)
            oddEvenTmp[1] = self.getExtraOutFile(tsId, suffix=EVEN, ext=REC_EXT)
            paramsTilt['OutputFile'] = oddEvenTmp[1]
            Plugin.runImod(self, 'tilt', argsTilt % paramsTilt)

        paramsTrimVol = {
            'input': self.getTmpOutFile(tsId, ext=REC_EXT),
            'output': self.getExtraOutFile(tsId, ext=MRC_EXT),
            'options': getArgs()
        }

        argsTrimvol = "%(options)s " \
                      "%(input)s " \
                      "%(output)s "

        Plugin.runImod(self, 'trimvol', argsTrimvol % paramsTrimVol)

        if self.applyToOddEven(ts):
            paramsTrimVol['input'] = oddEvenTmp[0]
            paramsTrimVol['output'] = self.getExtraOutFile(tsId, suffix=ODD, ext=MRCS_EXT)
            Plugin.runImod(self, 'trimvol', argsTrimvol % paramsTrimVol)

            paramsTrimVol['input'] = oddEvenTmp[1]
            paramsTrimVol['output'] = self.getExtraOutFile(tsId, suffix=EVEN, ext=MRCS_EXT)
            Plugin.runImod(self, 'trimvol', argsTrimvol % paramsTrimVol)

    def createOutputStep(self, tsId):
        ts = self.tsDict[tsId]
        tomoLocation = self.getExtraOutFile(tsId, ext=MRC_EXT)

        if os.path.exists(tomoLocation):
            output = self.getOutputSetOfTomograms(self.inputSetOfTiltSeries.get())

            newTomogram = Tomogram()
            newTomogram.setLocation(tomoLocation)

            if self.applyToOddEven(ts):
                halfMapsList = [self.getExtraOutFile(tsId, suffix=ODD, ext=MRCS_EXT),
                                self.getExtraOutFile(tsId, suffix=EVEN, ext=MRCS_EXT)]
                newTomogram.setHalfMaps(halfMapsList)

            newTomogram.setTsId(tsId)
            newTomogram.setSamplingRate(ts.getSamplingRate())

            # Set default tomogram origin
            newTomogram.setOrigin(newOrigin=None)
            if self.tomoShiftZ.get():
                x, y, z = newTomogram.getShiftsFromOrigin()
                shiftZang = self.tomoShiftZ.get() * newTomogram.getSamplingRate()
                newTomogram.setShiftsInOrigin(x=x, y=y, z=z + shiftZang)

            newTomogram.setAcquisition(ts.getAcquisition())

            output.append(newTomogram)
            output.update(newTomogram)
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

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if self.Tomograms:
            summary.append("Input tilt-series: %d\nTomograms reconstructed: %d"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.Tomograms.getSize()))
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.Tomograms:
            methods.append("The reconstruction has been computed for %d "
                           "tilt-series using the IMOD *tilt* command.\n"
                           % (self.Tomograms.getSize()))
        return methods
