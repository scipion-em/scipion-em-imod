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
from pyworkflow.utils import Message
from tomo.objects import CTFTomoSeries, SetOfCTFTomoSeries

from imod.protocols import ProtImodBase
from imod.constants import OUTPUT_CTF_SERIE, TLT_EXT, DEFOCUS_EXT


class ProtImodAutomaticCtfEstimation(ProtImodBase):
    """
    CTF estimation of a set of input tilt-series using the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/ctfplotter.html
    """

    _label = 'CTF estimation (auto)'
    _possibleOutputs = {OUTPUT_CTF_SERIE: SetOfCTFTomoSeries}
    _interactiveMode = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.sRate = None
        self.acq = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      label="Input tilt-series",
                      pointerClass='SetOfTiltSeries',
                      help='This should be a *raw stack*, not an aligned stack, '
                           'because the interpolation used to make '
                           'an aligned stack attenuates high frequencies and '
                           'the noise power spectra would no longer match.')

        form.addParam('defocusTol',
                      params.IntParam,
                      label='Defocus tolerance (nm)',
                      default=200,
                      important=True,
                      help='Defocus tolerance in nanometers defining the '
                           'center strips. The center strips are taken '
                           'from the central region of a view that has defocus '
                           'difference less than this tolerance. '
                           'These kind of center strips from all views within '
                           'AngleRange are considered to have a constant '
                           'defocus and are used to compute the initial CTF '
                           'after being further tessellated into tiles.')

        form.addParam('expectedDefocusOrigin',
                      params.EnumParam,
                      choices=['Value', 'List'],
                      default=0,
                      label='Input expected defocus as:',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST)

        form.addParam('expectedDefocusValue',
                      params.FloatParam,
                      default=6000.,
                      label='Expected defocus value (nm)',
                      important=True,
                      condition="expectedDefocusOrigin == 0",
                      help='This value will be applied as the expected '
                           'defocus in nanometers for every tilt-series '
                           'from the set.')

        form.addParam('expectedDefocusFile',
                      params.PathParam,
                      label='Expected defocus file',
                      important=True,
                      condition="expectedDefocusOrigin == 1",
                      help='File containing the expected defocus in nanometers for each '
                           'tilt-series belonging to the set.\n\n'
                           'The format of the text file must be two columns, '
                           'the first one being the tilt series ID '
                           'and the second the defocus value.\n\n'
                           'An example of this file comes as follows:\n'
                           'TS_01 4000\n'
                           'TS_02 1500\n'
                           '...')

        form.addParam('leftDefTol',
                      params.FloatParam,
                      label='Left defocus tolerance (nm)',
                      default=2000.,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Defocus tolerance in nanometers for strips "
                           "to the left of the center strip.")

        form.addParam('rightDefTol',
                      params.FloatParam,
                      label='Right defocus tolerance (nm)',
                      default=2000.,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Defocus tolerance in nanometers for strips "
                           "to the right of the center strip.")

        form.addParam('tileSize',
                      params.IntParam,
                      label='Tile size (px)',
                      default=256,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="The tile size each strip will be tessellated "
                           "into. The size is in pixels and the tiles are "
                           "square. Each view is first divided into strips "
                           "that are considered to have constant defocus.")

        if not self._interactiveMode:
            groupAngleRange = form.addGroup('Autorefinement angle settings',
                                            help='This entry sets the range of '
                                                 'angles for each fit and the step '
                                                 'size between angles for the '
                                                 'initial autofitting of the whole '
                                                 'tilt-series.')

            groupAngleRange.addParam('angleStep',
                                     params.FloatParam,
                                     default=2.0,
                                     label='Angle step',
                                     help='Step size between ranges. A value of '
                                          'zero for the step will make it fit to '
                                          'each single image separately, '
                                          'regardless of the value for the range.')

            groupAngleRange.addParam('angleRange',
                                     params.FloatParam,
                                     condition="angleStep != 0",
                                     default=16.0,
                                     label='Angle range',
                                     help='Size of the angle range for which the '
                                          'CTF is estimated.')

            groupFrequencyRange = form.addGroup('Autorefinement frequency range',
                                                expertLevel=params.LEVEL_ADVANCED,
                                                help='Starting and ending frequencies '
                                                     'of range to fit in power spectrum. '
                                                     'The two values will be used to '
                                                     'set the "X1 starts" and "X2 ends" '
                                                     'fields in the fitting dialog.')

            groupFrequencyRange.addParam('startFreq',
                                         params.FloatParam,
                                         default=0.0,
                                         label='Start',
                                         help='Starting frequency (X1 starts)')

            groupFrequencyRange.addParam('endFreq',
                                         params.FloatParam,
                                         default=0.0,
                                         label='End',
                                         help='Ending frequency (X2 ends)')

            form.addParam('extraZerosToFit',
                          params.FloatParam,
                          label='Extra zeros to fit',
                          default=0.0,
                          expertLevel=params.LEVEL_ADVANCED,
                          help="By default, the ending frequency of the fitting "
                               "range is set to the expected location of "
                               "the second zero. With this entry, the range will "
                               "be extended by the given multiple of the interval "
                               "between first and seconds zeros. For example, "
                               "entries of 1 and 2 will fit approximately to the "
                               "third and fourth zeros, respectively. An entry of "
                               "more than 0.5 will trigger fitting to two "
                               "exponentials, which is important for fitting "
                               "multiple peaks between zeros.")

            form.addParam('skipAstigmaticViews',
                          params.BooleanParam,
                          default=False,
                          label='Skip astigmatic phase views?',
                          expertLevel=params.LEVEL_ADVANCED,
                          help='Skip or break views only when finding astigmatism '
                               'or phase shift')

            groupAstigmatism = form.addGroup('Astigmatism settings',
                                             help='Parameters for astigmatism analysis')

            groupAstigmatism.addParam('searchAstigmatism',
                                      params.BooleanParam,
                                      default=False,
                                      label='Search astigmatism?',
                                      help='Search for astigmatism when fitting.')

            groupAstigmatism.addParam('maximumAstigmatism',
                                      params.FloatParam,
                                      default=1.2,
                                      label='Maximum astigmatism (um)',
                                      condition='searchAstigmatism',
                                      help='Maximum astigmatism, in microns. '
                                           'During the fitting to wedge spectra, '
                                           'the defocus is allowed to vary from '
                                           'the global value by more than half '
                                           'of this amount.')

            groupAstigmatism.addParam('numberSectorsAstigmatism',
                                      params.IntParam,
                                      default=36,
                                      label='Number of sectors',
                                      condition='searchAstigmatism',
                                      help='Number of sectors for astigmatism '
                                           'analysis.  A power spectrum is stored '
                                           'separately for each sector; spectra '
                                           'can then be computed fairly quickly '
                                           'for wedges of any size that is a '
                                           'multiple of the sector size. The '
                                           'default is 36, giving 5 degree sectors.')

            groupAstigmatism.addParam('minimumViewsAstigmatism',
                                      params.IntParam,
                                      default=3,
                                      label='Minimum views astigmatism',
                                      condition='searchAstigmatism',
                                      help='Minimum number of views for '
                                           'finding astigmatism.')

            groupPhaseShift = form.addGroup('Phase shift settings',
                                            help='Parameters for phase shift analysis')

            groupPhaseShift.addParam('searchPhaseShift',
                                     params.BooleanParam,
                                     default=False,
                                     label='Search phase shift?',
                                     help='Search for phase shift when fitting.')

            groupPhaseShift.addParam('minimumViewsPhaseShift',
                                     params.IntParam,
                                     default=1,
                                     label='Minimum views phase shift',
                                     condition='searchPhaseShift',
                                     help='Minimum number of views for '
                                          'finding phase shift.')

            groupCutOnFreq = form.addGroup('Cut-on frequency settings')

            groupCutOnFreq.addParam('searchCutOnFreq',
                                    params.BooleanParam,
                                    default=False,
                                    label='Search cut-on frequency?',
                                    help='Search for cut-on frequency when '
                                         'finding phase shift.')

            groupCutOnFreq.addParam('maximumCutOnFreq',
                                    params.FloatParam,
                                    default=-1,
                                    label='Maximum astigmatism (1/nm)',
                                    condition='searchCutOnFreq',
                                    help='Maximum frequency to test when searching '
                                         'for cut-on frequency, in reciprocal '
                                         'nanometers.  The default is the '
                                         'frequency of the first zero at the '
                                         'expected defocus and phase shift. '
                                         'To use the default value set box to -1.')

            form.addParallelSection(threads=4, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        expDefoci = self.getExpectedDefocus()

        closeSetStepDeps = []
        for tsId in self.tsDict.keys():
            convId = self._insertFunctionStep(self.convertInputStep, tsId,
                                              prerequisites=[])
            compId = self._insertFunctionStep(self.ctfEstimation, tsId,
                                              expDefoci, prerequisites=[convId])
            outId = self._insertFunctionStep(self.createOutputStep, tsId,
                                             prerequisites=[compId])
            closeSetStepDeps.append(outId)

        self._insertFunctionStep(self.closeOutputSetsStep,
                                 prerequisites=closeSetStepDeps)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        tsSet = self.getInputSet()
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in tsSet}
        self.sRate = tsSet.getSamplingRate()
        self.acq = tsSet.getAcquisition()

    def convertInputStep(self, tsId, **kwargs):
        """ Implement the convertStep to cancel interpolation of the tilt series."""
        super().convertInputStep(tsId, imodInterpolation=False)

    def ctfEstimation(self, tsId, expDefoci):
        """Run ctfplotter IMOD program"""
        try:
            ts = self.tsDict[tsId]

            paramsCtfPlotter = {
                "-InputStack": self.getTmpOutFile(tsId),
                "-AngleFile": self.getExtraOutFile(tsId, ext=TLT_EXT),
                "-DefocusFile": self.getExtraOutFile(tsId, ext=DEFOCUS_EXT),
                "-AxisAngle": ts.getAcquisition().getTiltAxisAngle(),
                "-PixelSize": self.sRate / 10,  # nm
                "-Voltage": int(self.acq.getVoltage()),
                "-SphericalAberration": self.acq.getSphericalAberration(),
                "-AmplitudeContrast": self.acq.getAmplitudeContrast(),
                "-DefocusTol": self.defocusTol.get(),
                "-PSResolution": 101,
                "-LeftDefTol": self.leftDefTol.get(),
                "-RightDefTol": self.rightDefTol.get(),
                "-tileSize": self.tileSize.get(),
            }

            if self.expectedDefocusOrigin.get() == 0:
                paramsCtfPlotter["-ExpectedDefocus"] = self.expectedDefocusValue.get()
            else:
                self.debug(f"Expected defoci: {expDefoci}")
                defocus = expDefoci.get(tsId, None)
                if defocus is None:
                    raise ValueError(f"{tsId} not found in the provided defocus file.")

                paramsCtfPlotter["-ExpectedDefocus"] = float(defocus)

            # Excluded views
            excludedViews = ts.getExcludedViewsIndex(caster=str)
            if excludedViews:
                paramsCtfPlotter["-ViewsToSkip"] = ",".join(excludedViews)

            if self._interactiveMode:
                paramsCtfPlotter["-AngleRange"] = "-20.0,20.0"
            else:
                if self.extraZerosToFit.get() != 0:
                    paramsCtfPlotter["-ExtraZerosToFit"] = self.extraZerosToFit.get()

                if self.skipAstigmaticViews:
                    paramsCtfPlotter["-SkipOnlyForAstigPhase"] = ""

                if self.searchAstigmatism:
                    paramsCtfPlotter.update({
                        "-SearchAstigmatism": "",
                        "-MaximumAstigmatism": self.maximumAstigmatism.get(),
                        "-NumberOfSectors": self.numberSectorsAstigmatism.get()
                    })

                if self.searchPhaseShift:
                    paramsCtfPlotter["-SearchPhaseShift"]: ""

                if self.searchAstigmatism and self.searchPhaseShift:
                    paramsCtfPlotter["-MinViewsAstigAndPhase"] = f"{self.minimumViewsAstigmatism.get()},"\
                                                                 f"{self.minimumViewsPhaseShift.get()}"
                elif self.searchAstigmatism:
                    paramsCtfPlotter["-MinViewsAstigAndPhase"] = f"{self.minimumViewsAstigmatism.get()},0"
                elif self.searchPhaseShift:
                    paramsCtfPlotter["-MinViewsAstigAndPhase"] = f"0,{self.minimumViewsPhaseShift.get()}"

                if self.searchCutOnFreq:
                    paramsCtfPlotter["-SearchCutonFrequency"] = ""

                    if self.maximumCutOnFreq.get() != -1.0:
                        paramsCtfPlotter["-MaxCutOnToSearch"] = self.maximumCutOnFreq.get()

                paramsCtfPlotter["-AutoFitRangeAndStep"] = f"{self.angleRange.get()},{self.angleStep.get()}"
                paramsCtfPlotter["-SaveAndExit"] = ""

                if self.startFreq.get() != 0 or self.endFreq.get() != 0:
                    paramsCtfPlotter["-FrequencyRangeToFit"] = f"{self.startFreq.get()},{self.endFreq.get()}"

            self.runProgram('ctfplotter', paramsCtfPlotter)

        except Exception as e:
            self._failedItems.append(tsId)
            self.error(f'ctfplotter execution failed for tsId {tsId} -> {e}')

    def createOutputStep(self, tsId, outputSetName=OUTPUT_CTF_SERIE):
        ts = self.tsDict[tsId]

        with self._lock:
            if tsId in self._failedItems:
                self.createOutputFailedSet(ts)
            else:
                defocusFilePath = self.getExtraOutFile(tsId, ext=DEFOCUS_EXT)
                if os.path.exists(defocusFilePath):
                    output = self.getOutputSetOfCTFTomoSeries(self.getInputSet(pointer=True),
                                                              outputSetName)
                    newCTFTomoSeries = CTFTomoSeries(tsId=tsId)
                    newCTFTomoSeries.copyInfo(ts)
                    newCTFTomoSeries.setTiltSeries(ts)

                    # flags below will be updated in parseTSDefocusFile
                    newCTFTomoSeries.setIMODDefocusFileFlag(1)
                    newCTFTomoSeries.setNumberOfEstimationsInRange(0)
                    output.append(newCTFTomoSeries)

                    self.parseTSDefocusFile(ts, defocusFilePath, newCTFTomoSeries)

                    output.update(newCTFTomoSeries)
                    output.write()
                    self._store(output)
                else:
                    self.createOutputFailedSet(ts)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        ctfSeries = getattr(self, OUTPUT_CTF_SERIE, None)
        if ctfSeries is not None:
            summary.append(f"Input tilt-series: {self.getInputSet().getSize()}\n"
                           f"Number of CTF estimated: {ctfSeries.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []

        ctfSeries = getattr(self, OUTPUT_CTF_SERIE, None)
        if ctfSeries is not None:
            methods.append(f"{ctfSeries.getSize()} tilt-series CTF "
                           "have been estimated using the IMOD *ctfplotter* "
                           "command.")
        return methods

    # --------------------------- UTILS functions -----------------------------
    def allowsDelete(self, obj):
        return True

    def getExpectedDefocus(self):
        if self.expectedDefocusOrigin.get() == 1:
            with open(self.expectedDefocusFile.get()) as f:
                lines = f.readlines()
                lines = filter(lambda x: x.strip(), lines)
            result = {line.split()[0]: line.split()[1] for line in lines}
            return result
        else:
            return None

    def getInputSet(self, pointer=False):
        """ Reimplemented from the base class for CTF case. """
        inputSet = getattr(self, IN_TS_SET)
        return inputSet.get() if not pointer else inputSet
