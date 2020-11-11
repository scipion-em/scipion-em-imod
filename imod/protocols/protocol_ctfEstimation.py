# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
# **************************************************************************

import os
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
from imod import Plugin


class ProtImodCtfEstimation(EMProtocol, ProtTomoBase):
    """
    CTF estimation of a set of input tilt-series using the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/ctfplotter.html
    """

    _label = 'CTF estimation'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series',
                      help='This should be a raw stack, not an aligned stack, because the interpolation used to make '
                           'an aligned stack attenuates high frequencies and the noise power spectra would no longer '
                           'match.')

        form.addParam('defocusTol',
                      params.FloatParam,
                      label='Defocus tolerance',
                      default=200,
                      important=True,
                      help='Defocus tolerance in nanometers defining the center strips. The center strips are taken '
                           'from the central region of a view that has defocus difference less than this tolerance. '
                           'These kind of center strips from all views within AngleRange are considered to have a '
                           'constant defocus and are used to compute the initial CTF after being further tessellated '
                           'into tiles.')

        form.addParam('expectedDefocusOrigin',
                      params.EnumParam,
                      choices=['Value', 'List'],
                      default=0,
                      label='Input expected defocus as',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Run the protocol through the interactive GUI. If run in auto mode defocus values are saved '
                           'to file and exit after autofitting. The program will not ask for confirmation before '
                           'removing existing entries in the defocus table. If run in interactive mode defocus values'
                           'MUST BE SAVED manually by the user.')

        form.addParam('expectedDefocusValue',
                      params.FloatParam,
                      default=6000,
                      label='Expected defocus value',
                      important=True,
                      condition="expectedDefocusOrigin == 0",
                      help='This value will be applied as the expected defocus in nanometers for every tilt-series '
                           'from the set.')

        form.addParam('expectedDefocusFile',
                      params.PathParam,
                      label='Expected defocus file',
                      important=True,
                      condition="expectedDefocusOrigin == 1",
                      help='File containing a list of expected defoci in nanometes for each tilt-series of the set. '
                           'This file must contain two columns. The first column must be the filename of the '
                           'tilt-series and the second the expected defocus.')

        form.addParam('axisAngle',
                      params.FloatParam,
                      default=0.0,
                      label='Axis angle',
                      important=True,
                      help='Specifies how much the tilt axis deviates from vertical (Y axis). This angle is in degrees.'
                           ' It follows the right hand  rule and counter-clockwise is positive.')

        form.addParam('interactiveMode',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Run interactive GUI',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Run the protocol through the interactive GUI. If run in auto mode defocus values are saved '
                           'to file and exit after autofitting. The program will not ask for confirmation before '
                           'removing existing entries in the defocus table. If run in interactive mode defocus values'
                           'MUST BE SAVED manually by the user.')

        form.addParam('leftDefTol',
                      params.FloatParam,
                      label='Left defocus tolerance',
                      default=2000,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Defocus tolerance in nanometers for strips to the left of the center strip.")

        form.addParam('rightDefTol',
                      params.FloatParam,
                      label='Right defocus tolerance',
                      default=2000,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Defocus tolerance in nanometers for strips to the right of the center strip.")

        form.addParam('tileSize',
                      params.IntParam,
                      label='Tile size',
                      default=256,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="The tile size each strip will be tessellated into. The size is in pixels and the tiles are "
                           "square.  Each view is first divided into strips that are considered to have constant "
                           "defocus.")

        groupAngleRange = form.addGroup('Autorefinement angle settings',
                                        help='This entry sets the range of angles in each fit and the step size between'
                                             'angles for the initial autofitting of the whole tilt-series.')

        groupAngleRange.addParam('angleStep',
                                 params.FloatParam,
                                 default=2,
                                 label='Angle step',
                                 help='Step size between ranges. A value of zero for the step will make it fit to each '
                                      'single image separately, regardless of the value for the range.')

        groupAngleRange.addParam('angleRange',
                                 params.FloatParam,
                                 condition="angleStep != 0",
                                 default=120,
                                 label='Angle range',
                                 help='Size of the angle range in which the CTF is estimated.')

        groupFrequencyRange = form.addGroup('Autorefinement frequency range',
                                            expertLevel=params.LEVEL_ADVANCED,
                                            help='Starting and ending frequencies of range to fit in power spectrum. '
                                                 'The two values will be used to set the "X1 starts" and "X2 ends" '
                                                 'fields in the fitting dialog.')

        groupFrequencyRange.addParam('startFreq',
                                     params.FloatParam,
                                     default=0.0,
                                     label='Start',
                                     help='Starting frequency. "X1 starts". "X2 ends".')

        groupFrequencyRange.addParam('endFreq',
                                     params.FloatParam,
                                     default=0.0,
                                     label='End',
                                     help='Ending frequency.')

        form.addParam('extraZerosToFit',
                      params.FloatParam,
                      label='Extra zeros to fit',
                      default=0.0,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="By default, the ending frequency of the fitting range is set to the expected location of "
                           "the second zero.  With this entry, the range will be extended by the given multiple of the "
                           "interval between first and seconds zeros.  For example, entries of 1 and 2 will fit "
                           "approximately to the third and fourth zeros, respectively.  An entry of more than 0.5 will "
                           "trigger fitting to two exponentials, which is important for fitting multiple peaks between "
                           "zeros.")

        form.addParam('skipAstigmaticViews',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Skip astigmatic phase views',
                      expertLevel=params.LEVEL_ADVANCED,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Skip or break views only when finding astigmatism or phase shift')

        groupAstigmatism = form.addGroup('Astigmatism settings',
                                         expertLevel=params.LEVEL_ADVANCED,
                                         help='Parameters for astigmatism analysis')

        groupAstigmatism.addParam('searchAstigmatism',
                                  params.EnumParam,
                                  choices=['Yes', 'No'],
                                  default=1,
                                  label='Search astigmatism',
                                  display=params.EnumParam.DISPLAY_HLIST,
                                  help='Search for astigmatism when fitting')

        groupAstigmatism.addParam('findAstigPhaseCutonToggle',
                                  params.EnumParam,
                                  choices=['Yes', 'No'],
                                  default=1,
                                  label='Find astigmatism, phase shift, and cut-on frequency?',
                                  display=params.EnumParam.DISPLAY_HLIST,
                                  condition='searchAstigmatism==0',
                                  help='Find astigmatism, phase shift, and cut-on frequency')

        groupAstigmatism.addParam('phaseShiftAstigmatism',
                                  params.IntParam,
                                  label='Phase shift',
                                  condition='searchAstigmatism==0 and findAstigPhaseCutonToggle==0',
                                  help='Phase shift for astigmatism analysis.')

        groupAstigmatism.addParam('cutOnFrequencyAstigmatism',
                                  params.IntParam,
                                  label='Cut-on frequency',
                                  condition='searchAstigmatism==0 and findAstigPhaseCutonToggle==0',
                                  help='Cut-on frequency for astigmatism analysis.')

        groupAstigmatism.addParam('minimumViewsAstigmatism',
                                  params.IntParam,
                                  default=3,
                                  label='Minimum views astigmatism',
                                  condition='searchAstigmatism==0 and findAstigPhaseCutonToggle==0',
                                  help='Minimum views for finding astigmatism.')

        groupAstigmatism.addParam('minimumViewsPhaseShift',
                                  params.IntParam,
                                  default=1,
                                  label='Minimum views phase shift',
                                  condition='searchAstigmatism==0 and findAstigPhaseCutonToggle==0',
                                  help='Minimum views for finding phase shift.')

        groupAstigmatism.addParam('numberSectorsAstigmatism',
                                  params.IntParam,
                                  default=36,
                                  label='Number of sectors',
                                  condition='searchAstigmatism==0',
                                  help='Number of sectors for astigmatism analysis.')

        groupAstigmatism.addParam('maximumAstigmatism',
                                  params.FloatParam,
                                  default=1.2,
                                  label='Maximum astigmatism',
                                  condition='searchAstigmatism==0',
                                  help='Maximum astigmatism in microns.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('ctfEstimation', ts.getObjId())
            self._insertFunctionStep('createOutputStep', ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        """Apply the transformation form the input tilt-series"""
        outputTsFileName = os.path.join(tmpPrefix, ts.getFirstItem().parseFileName())
        ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, ts.getFirstItem().parseFileName(extension=".tlt"))
        ts.generateTltFile(angleFilePath)

    def ctfEstimation(self, tsObjId):
        """Run ctfplotter IMOD program"""
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsCtfPlotter = {
            'inputStack': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName()),
            'angleFile': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName(extension=".tlt")),
            'defocusFile': os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".defocus")),
            'axisAngle': self.axisAngle.get(),
            'pixelSize': self.inputSetOfTiltSeries.get().getSamplingRate() / 10,
            'expectedDefocus': self.getExpectedDefocus(tsId),
            'autoFitRangeAndStep': str(self.angleRange.get()) + "," + str(self.angleStep.get()),
            'voltage': self.inputSetOfTiltSeries.get().getAcquisition().getVoltage(),
            'sphericalAberration': self.inputSetOfTiltSeries.get().getAcquisition().getSphericalAberration(),
            'amplitudeContrast': self.inputSetOfTiltSeries.get().getAcquisition().getAmplitudeContrast(),
            'defocusTol': self.defocusTol.get(),
            'psResolution': 101,
            'leftDefTol': self.leftDefTol.get(),
            'rightDefTol': self.rightDefTol.get(),
            'tileSize': self.tileSize.get(),
        }

        argsCtfPlotter = "-InputStack %(inputStack)s " \
                         "-AngleFile %(angleFile)s " \
                         "-DefocusFile %(defocusFile)s " \
                         "-AxisAngle %(axisAngle)f " \
                         "-PixelSize %(pixelSize)f " \
                         "-ExpectedDefocus %(expectedDefocus)f " \
                         "-AutoFitRangeAndStep %(autoFitRangeAndStep)s " \
                         "-Voltage %(voltage)d " \
                         "-SphericalAberration %(sphericalAberration)f " \
                         "-AmplitudeContrast %(amplitudeContrast)f " \
                         "-DefocusTol %(defocusTol)d " \
                         "-PSResolution %(psResolution)d " \
                         "-LeftDefTol %(leftDefTol)f " \
                         "-RightDefTol %(rightDefTol)f " \
                         "-tileSize %(tileSize)d " \

        if self.startFreq.get() != 0 or self.endFreq.get() != 0:
            paramsCtfPlotter.update({
                'startFreq': self.startFreq.get(),
                'endFreq': self.endFreq.get(),
            })

            argsCtfPlotter += "-FrequencyRangeToFit %(startFreq)f,%(endFreq)f "

        if self.extraZerosToFit.get() != 0:
            paramsCtfPlotter.update({
                'extraZerosToFit': self.extraZerosToFit.get(),
            })

            argsCtfPlotter += "-ExtraZerosToFit %(extraZerosToFit)f "

        if self.skipAstigmaticViews.get() == 0:
            argsCtfPlotter += "-SkipOnlyForAstigPhase "

        if self.searchAstigmatism.get() == 0:

            if self.findAstigPhaseCutonToggle.get() == 0:
                paramsCtfPlotter.update({
                    'phaseShiftAstigmatism': self.phaseShiftAstigmatism.get(),
                    'cutOnFrequencyAstigmatism': self.cutOnFrequencyAstigmatism.get(),
                    'minimumViewsAstigmatism': self.minimumViewsAstigmatism.get(),
                    'minimumViewsPhaseShift': self.minimumViewsPhaseShift.get(),
                })

                argsCtfPlotter += "-FindAstigPhaseCuton 1,%(phaseShiftAstigmatism)d,%(cutOnFrequencyAstigmatism)d " \
                                  "-MinViewsAstigAndPhase %(minimumViewsAstigmatism)d,%(minimumViewsPhaseShift)d "

            else:
                argsCtfPlotter += "-SearchAstigmatism " \

            paramsCtfPlotter.update({
                'numberSectorsAstigmatism': self.numberSectorsAstigmatism.get(),
                'maximumAstigmatism': self.maximumAstigmatism.get(),
            })

            argsCtfPlotter += "-NumberOfSectors %(numberSectorsAstigmatism)d " \
                              "-MaximumAstigmatism %(maximumAstigmatism)f " \

        if self.interactiveMode.get() == 1:
            argsCtfPlotter += "-SaveAndExit "

        Plugin.runImod(self, 'ctfplotter', argsCtfPlotter % paramsCtfPlotter)

    def createOutputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        if os.path.exists(os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".defocus"))):
            outputCtfEstimatedSetOfTiltSeries = self.getOutputCtfEstimatedSetOfTiltSeries()

            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            outputCtfEstimatedSetOfTiltSeries.append(newTs)

            for index, tiltImage in enumerate(ts):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)
                newTi.setLocation(tiltImage.getLocation())
                if tiltImage.hasTransform():
                    newTi.setTransform(tiltImage.getTransform())
                newTs.append(newTi)

            newTs.write(properties=False)
            outputCtfEstimatedSetOfTiltSeries.update(newTs)
            outputCtfEstimatedSetOfTiltSeries.write()
            self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputCtfEstimatedSetOfTiltSeries(self):
        if hasattr(self, "outputCtfEstimatedSetOfTiltSeries"):
            self.outputCtfEstimatedSetOfTiltSeries.enableAppend()
        else:
            outputCtfEstimatedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='CtfEstimated')
            outputCtfEstimatedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputCtfEstimatedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            outputCtfEstimatedSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputCtfEstimatedSetOfTiltSeries=outputCtfEstimatedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputCtfEstimatedSetOfTiltSeries)
        return self.outputCtfEstimatedSetOfTiltSeries

    def getExpectedDefocus(self, tsId):
        if self.expectedDefocusOrigin.get() == 0:
            return self.expectedDefocusValue.get()
        else:
            with open(self.expectedDefocusFile.get()) as f:
                lines = f.readlines()
                for line in lines:
                    defocusTuple = line.split()
                    if tsId in defocusTuple[0]:  # Look for the filename that contains the tsId
                        return float(defocusTuple[1])
                raise Exception("ERROR: tilt-series with tsId %s has not been found in %s" %
                                (tsId, (self.expectedDefocusFile.get())))

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputCtfEstimatedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nnumber of CTF estimated: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputCtfEstimatedSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputCtfEstimatedSetOfTiltSeries'):
            methods.append("%d Tilt-series CTF have been estimated using the IMOD ctfplotter software.\n"
                           % (self.outputCtfEstimatedSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
