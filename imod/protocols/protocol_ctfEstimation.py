# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
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
from pyworkflow import BETA
from pyworkflow.object import Set
import pyworkflow.object as pwobj
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
from imod import Plugin
from imod import utils


class ProtImodCtfEstimation(EMProtocol, ProtTomoBase):
    """
    CTF estimation of a set of input tilt-series using the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/ctfplotter.html
    """

    _label = 'CTF estimation'
    _devStatus = BETA

    defocusUTolerance = 20
    defocusVTolerance = 20

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSet',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries, SetOfCTFTomoSeries',
                      label='Input set of tilt-series',
                      help='This should be a raw stack, not an aligned stack, because the interpolation used to make '
                           'an aligned stack attenuates high frequencies and the noise power spectra would no longer '
                           'match.')

        form.addParam('defocusTol',
                      params.FloatParam,
                      label='Defocus tolerance (nm)',
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
                      label='Expected defocus value (nm)',
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
                      label='Left defocus tolerance (nm)',
                      default=2000,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Defocus tolerance in nanometers for strips to the left of the center strip.")

        form.addParam('rightDefTol',
                      params.FloatParam,
                      label='Right defocus tolerance (nm)',
                      default=2000,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Defocus tolerance in nanometers for strips to the right of the center strip.")

        form.addParam('tileSize',
                      params.IntParam,
                      label='Tile size (pixels)',
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
                                 default=16,
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
                                         help='Parameters for astigmatism analysis')

        groupAstigmatism.addParam('searchAstigmatism',
                                  params.EnumParam,
                                  choices=['Yes', 'No'],
                                  default=1,
                                  label='Search astigmatism',
                                  display=params.EnumParam.DISPLAY_HLIST,
                                  help='Search for astigmatism when fitting.')

        groupAstigmatism.addParam('maximumAstigmatism',
                                  params.FloatParam,
                                  default=1.2,
                                  label='Maximum astigmatism (um)',
                                  condition='searchAstigmatism==0',
                                  help='Maximum astigmatism, in microns.  During the fitting to wedge spectra, the '
                                       'defocus is allowed to vary from the global value by more than half of this '
                                       'amount.')

        groupAstigmatism.addParam('numberSectorsAstigmatism',
                                  params.IntParam,
                                  default=36,
                                  label='Number of sectors',
                                  condition='searchAstigmatism==0',
                                  help='Number of sectors for astigmatism analysis.  A power spectrum is stored '
                                       'separately for each sector; spectra can then be computed fairly quickly for '
                                       'wedges of any size that is a multiple of the sector size.  The default is 36, '
                                       'giving 5 degree sectors.')

        groupAstigmatism.addParam('minimumViewsAstigmatism',
                                  params.IntParam,
                                  default=3,
                                  label='Minimum views astigmatism',
                                  condition='searchAstigmatism==0',
                                  help='Minimum number of views for finding astigmatism.')

        groupPhaseShift = form.addGroup('Phase shift settings',
                                        help='Parameters for phase shift analysis')

        groupPhaseShift.addParam('searchPhaseShift',
                                 params.EnumParam,
                                 choices=['Yes', 'No'],
                                 default=1,
                                 label='Search phase shift',
                                 display=params.EnumParam.DISPLAY_HLIST,
                                 help='Search for phase shift when fitting.')

        groupPhaseShift.addParam('minimumViewsPhaseShift',
                                 params.IntParam,
                                 default=1,
                                 label='Minimum views phase shift',
                                 condition='searchPhaseShift==0',
                                 help='Minimum number of views for finding phase shift.')

        groupCutOnFreq = form.addGroup('Cut-on frequency settings')

        groupCutOnFreq.addParam('searchCutOnFreq',
                                params.EnumParam,
                                choices=['Yes', 'No'],
                                default=1,
                                label='Search cut-on frequency',
                                display=params.EnumParam.DISPLAY_HLIST,
                                help='Search for cut-on frequency when finding phase shift.')

        groupCutOnFreq.addParam('maximumCutOnFreq',
                                params.FloatParam,
                                default=-1,
                                label='Maximum astigmatism (um)',
                                condition='searchCutOnFreq==0',
                                help='Maximum frequency to test when searching for cut-on frequency, in reciprocal '
                                     'nanometers.  The default is the frequency of the first zero at the expected '
                                     'defocus and phase shift. To use the default value set box to -1.')

    # -------------------------- INSERT steps functions ---------------------
    def _getSetOfTiltSeries(self):
        if isinstance(self.inputSet.get(), tomoObj.SetOfCTFTomoSeries):
            return self.inputSet.get().getSetOfTiltSeries()
        return self.inputSet.get()

    def _getTiltSeries(self, itemId):
        ts = self.inputSet.get()[itemId]
        if isinstance(ts, tomoObj.CTFTomoSeries):
            ts = ts.getTiltSeries()
        return ts

    def _insertAllSteps(self):
        for item in self.inputSet.get():
            self._insertFunctionStep('convertInputStep', item.getObjId())
            self._insertFunctionStep('ctfEstimation', item.getObjId())
            self._insertFunctionStep('createOutputStep', item.getObjId())
        self._insertFunctionStep('closeOutputSetsStep')

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self._getTiltSeries(tsObjId)
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
        ts = self._getTiltSeries(tsObjId)
        tsSet = self._getSetOfTiltSeries()
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsCtfPlotter = {
            'inputStack': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName()),
            'angleFile': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName(extension=".tlt")),
            'defocusFile': os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".defocus")),
            'axisAngle': self.axisAngle.get(),
            'pixelSize': tsSet.getSamplingRate() / 10,
            'expectedDefocus': self.getExpectedDefocus(tsId),
            'autoFitRangeAndStep': str(self.angleRange.get()) + "," + str(self.angleStep.get()),
            'voltage': tsSet.getAcquisition().getVoltage(),
            'sphericalAberration': tsSet.getAcquisition().getSphericalAberration(),
            'amplitudeContrast': tsSet.getAcquisition().getAmplitudeContrast(),
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
                         "-tileSize %(tileSize)d "

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
            paramsCtfPlotter.update({
                'maximumAstigmatism': self.maximumAstigmatism.get(),
                'numberOfSectors': self.numberSectorsAstigmatism.get(),
                'minimumViewsAstigmatism': self.minimumViewsAstigmatism.get()
            })

            argsCtfPlotter += "-SearchAstigmatism " \
                              "-MaximumAstigmatism %(maximumAstigmatism)f " \
                              "-NumberOfSectors %(numberOfSectors)d "

        if self.searchPhaseShift.get() == 0:
            paramsCtfPlotter.update({
                'minimumViewsPhaseShift': self.minimumViewsPhaseShift.get()
            })

            argsCtfPlotter += "-SearchPhaseShift "

        if self.searchAstigmatism.get() == 0 and self.searchPhaseShift.get() == 0:
            argsCtfPlotter += "-MinViewsAstigAndPhase %(minimumViewsAstigmatism)d,%(minimumViewsPhaseShift)d "
        elif self.searchAstigmatism.get() == 0:
            argsCtfPlotter += "-MinViewsAstigAndPhase %(minimumViewsAstigmatism)d,0 "
        elif self.searchPhaseShift.get() == 0:
            argsCtfPlotter += "-MinViewsAstigAndPhase 0,%(minimumViewsPhaseShift)d "

        if self.searchCutOnFreq.get() == 0:

            if self.maximumCutOnFreq.get() == -1:
                argsCtfPlotter += "-SearchCutonFrequency "

            else:
                paramsCtfPlotter.update({
                    'maximumCutOnFreq': self.maximumCutOnFreq.get(),
                })

                argsCtfPlotter += "-SearchCutonFrequency " \
                                  "-MaxCutOnToSearch %(maximumCutOnFreq)f "

        if self.interactiveMode.get() == 1:
            argsCtfPlotter += "-SaveAndExit "

        Plugin.runImod(self, 'ctfplotter', argsCtfPlotter % paramsCtfPlotter)

    def createOutputStep(self, tsObjId):
        ts = self._getTiltSeries(tsObjId)
        tsId = ts.getTsId()
        objId = ts.getObjId()

        extraPrefix = self._getExtraPath(tsId)

        defocusFilePath = os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".defocus"))

        if os.path.exists(defocusFilePath):

            self.getOutputSetOfCTFTomoSeries()

            defocusFileFlag = utils.getDefocusFileFlag(defocusFilePath)

            newCTFTomoSeries = tomoObj.CTFTomoSeries()
            newCTFTomoSeries.copyInfo(ts)
            newCTFTomoSeries.setTiltSeries(ts)
            newCTFTomoSeries.setTsId(tsId)
            newCTFTomoSeries.setObjId(objId)
            newCTFTomoSeries.setIMODDefocusFileFlag(defocusFileFlag)

            # We need to create now all the attributes of this object in order
            # to append it to the set and be able to update it posteriorly. "

            newCTFTomoSeries.setNumberOfEstimationsInRange(None)
            self.outputSetOfCTFTomoSeries.append(newCTFTomoSeries)

            if defocusFileFlag == 0:
                " Plain estimation "
                defocusUDict = utils.readCTFEstimationInfoFile(defocusFilePath,
                                                               flag=defocusFileFlag)

            elif defocusFileFlag == 1:
                " Astigmatism estimation "
                defocusUDict, defocusVDict, defocusAngleDict = utils.readCTFEstimationInfoFile(defocusFilePath,
                                                                                               flag=defocusFileFlag)

            elif defocusFileFlag == 4:
                " Phase-shift information "
                defocusUDict, phaseShiftDict = utils.readCTFEstimationInfoFile(defocusFilePath,
                                                                               flag=defocusFileFlag)

            elif defocusFileFlag == 5:
                " Astigmatism and phase shift estimation "
                defocusUDict, defocusVDict, defocusAngleDict, phaseShiftDict = \
                    utils.readCTFEstimationInfoFile(defocusFilePath,
                                                    flag=defocusFileFlag)

            elif defocusFileFlag == 37:
                " Astigmatism, phase shift and cut-on frequency estimation "
                defocusUDict, defocusVDict, defocusAngleDict, phaseShiftDict, cutOnFreqDict = \
                    utils.readCTFEstimationInfoFile(defocusFilePath,
                                                    flag=defocusFileFlag)

            else:
                raise Exception("Defocus file flag do not supported. Only supported formats corresponding to flags 0, "
                                "1, 4, 5, and 37.")

            for index, _ in enumerate(ts):
                newCTFTomo = tomoObj.CTFTomo()
                newCTFTomo.setIndex(pwobj.Integer(index + 1))

                if defocusFileFlag == 0:
                    " Plain estimation "
                    newCTFTomo._defocusUList = pwobj.CsvList(pType=float)
                    newCTFTomo.setDefocusUList(defocusUDict[index + 1])

                elif defocusFileFlag == 1:
                    " Astigmatism estimation "
                    newCTFTomo._defocusUList = pwobj.CsvList(pType=float)
                    newCTFTomo.setDefocusUList(defocusUDict[index + 1])

                    newCTFTomo._defocusVList = pwobj.CsvList(pType=float)
                    newCTFTomo.setDefocusVList(defocusVDict[index + 1])

                    newCTFTomo._defocusAngleList = pwobj.CsvList(pType=float)
                    newCTFTomo.setDefocusAngleList(defocusAngleDict[index + 1])

                elif defocusFileFlag == 4:
                    " Phase-shift information "
                    newCTFTomo._defocusUList = pwobj.CsvList(pType=float)
                    newCTFTomo.setDefocusUList(defocusUDict[index + 1])

                    newCTFTomo._phaseShiftList = pwobj.CsvList(pType=float)
                    newCTFTomo.setPhaseShiftList(phaseShiftDict[index + 1])

                elif defocusFileFlag == 5:
                    " Astigmatism and phase shift estimation "
                    newCTFTomo._defocusUList = pwobj.CsvList(pType=float)
                    newCTFTomo.setDefocusUList(defocusUDict[index + 1])

                    newCTFTomo._defocusVList = pwobj.CsvList(pType=float)
                    newCTFTomo.setDefocusVList(defocusVDict[index + 1])

                    newCTFTomo._defocusAngleList = pwobj.CsvList(pType=float)
                    newCTFTomo.setDefocusAngleList(defocusAngleDict[index + 1])

                    newCTFTomo._phaseShiftList = pwobj.CsvList(pType=float)
                    newCTFTomo.setPhaseShiftList(phaseShiftDict[index + 1])

                elif defocusFileFlag == 37:
                    " Astigmatism, phase shift and cut-on frequency estimation "
                    newCTFTomo._defocusUList = pwobj.CsvList(pType=float)
                    newCTFTomo.setDefocusUList(defocusUDict[index + 1])

                    newCTFTomo._defocusVList = pwobj.CsvList(pType=float)
                    newCTFTomo.setDefocusVList(defocusVDict[index + 1])

                    newCTFTomo._defocusAngleList = pwobj.CsvList(pType=float)
                    newCTFTomo.setDefocusAngleList(defocusAngleDict[index + 1])

                    newCTFTomo._phaseShiftList = pwobj.CsvList(pType=float)
                    newCTFTomo.setPhaseShiftList(phaseShiftDict[index + 1])

                    newCTFTomo._cutOnFreqList = pwobj.CsvList(pType=float)
                    newCTFTomo.setCutOnFreqList(cutOnFreqDict[index + 1])

                    defocusUDict, defocusVDict, defocusAngleDict, phaseShiftDict, cutOnFreqDict = \
                        utils.readCTFEstimationInfoFile(defocusFilePath,
                                                        flag=defocusFileFlag)

                newCTFTomo.completeInfoFromList()
                newCTFTomoSeries.append(newCTFTomo)

            newCTFTomoSeries.setNumberOfEstimationsInRangeFromDefocusList()

            newCTFTomoSeries.calculateDefocusUDeviation(defocusUTolerance=self.defocusUTolerance)
            newCTFTomoSeries.calculateDefocusVDeviation(defocusVTolerance=self.defocusVTolerance)

            if not (newCTFTomoSeries.getIsDefocusUDeviationInRange() and
                    newCTFTomoSeries.getIsDefocusVDeviationInRange()):
                newCTFTomoSeries.setEnabled(False)

            newCTFTomoSeries.write(properties=False)

            self.outputSetOfCTFTomoSeries.update(newCTFTomoSeries)
            self.outputSetOfCTFTomoSeries.write()

            self._store()

    def closeOutputSetsStep(self):
        self.outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfCTFTomoSeries(self):
        if hasattr(self, "outputSetOfCTFTomoSeries"):
            self.outputSetOfCTFTomoSeries.enableAppend()
        else:
            outputSetOfCTFTomoSeries = tomoObj.SetOfCTFTomoSeries.create(self._getPath(),
                                                                         template='CTFmodels%s.sqlite')
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(self._getSetOfTiltSeries())
            outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfCTFTomoSeries=outputSetOfCTFTomoSeries)
        return self.outputSetOfCTFTomoSeries

    def getExpectedDefocus(self, tsId):
        if self.expectedDefocusOrigin.get() == 0:
            return self.expectedDefocusValue.get()
        else:
            with open(self.expectedDefocusFile.get()) as f:
                lines = f.readlines()
                for line in lines:
                    defocusTuple = line.split()
                    if tsId in defocusTuple[0]:
                        " Look for the filename that contains the tsId "
                        return float(defocusTuple[1])
                raise Exception("ERROR: tilt-series with tsId %s has not been found in %s" %
                                (tsId, (self.expectedDefocusFile.get())))

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputSetOfCTFTomoSeries'):
            summary.append("Input Tilt-Series: %d.\nnumber of CTF estimated: %d.\n"
                           % (self._getSetOfTiltSeries().getSize(),
                              self.outputSetOfCTFTomoSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputSetOfCTFTomoSeries'):
            methods.append("%d Tilt-series CTF have been estimated using the IMOD ctfplotter software.\n"
                           % (self.outputSetOfCTFTomoSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
