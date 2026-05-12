# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
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
import logging
import traceback
from collections import Counter
from os.path import exists
from typing import Union
import pyworkflow.protocol.params as params
from imod.convert.dataimport import ImodCtfParser
from pyworkflow.protocol import STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.utils import Message, cyanStr, redStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import CTFTomoSeries, SetOfCTFTomoSeries, TiltSeries
from imod.protocols import ProtImodBase
from imod.constants import OUTPUT_CTF_SERIE, TLT_EXT, DEFOCUS_EXT, CTFPLOTTER_PROGRAM

logger = logging.getLogger(__name__)


class ProtImodAutomaticCtfEstimation(ProtImodBase, ProtStreamingBase):
    """
    Performs automatic Contrast Transfer Function estimation for cryo-electron
    tomography tilt-series using the IMOD ctfplotter workflow. The protocol is
    designed to analyze the optical behavior of tilt images and determine
    defocus-related parameters required for accurate tomographic reconstruction,
    downstream correction procedures, and high-quality structural interpretation.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/ctfplotter.html

    AI Generated:

    Automatic CTF Estimation (ProtImodAutomaticCtfEstimation) — User Manual
        Overview

        The Automatic CTF Estimation protocol provides an automated strategy
        for estimating the Contrast Transfer Function of tilt-series acquired
        for cryo-electron tomography experiments. Its primary objective is to
        determine defocus and related optical parameters across the complete
        tilt-series so that subsequent reconstruction and refinement steps can
        be performed with improved accuracy and reliability.

        In cryo-electron tomography, each tilt image is acquired under varying
        geometric conditions that produce changes in apparent defocus across
        the field of view. Correctly characterizing these effects is essential
        for preserving high-resolution information and for minimizing artifacts
        during tomogram reconstruction. The protocol is particularly useful in
        workflows where large collections of tilt-series must be processed in
        an automated or streaming manner.

        Biological Motivation

        Accurate CTF estimation is one of the most important preprocessing
        stages in cryo-electron tomography because it directly affects the
        interpretability of the reconstructed density. Poorly estimated
        defocus values can lead to reduced contrast recovery, incorrect phase
        correction, and loss of structural detail. For biological users, this
        can compromise the interpretation of macromolecular assemblies,
        membrane organization, or in situ cellular architecture.

        The protocol is intended to support both routine data collection and
        demanding high-resolution studies. It can be applied to datasets
        acquired with conventional defocus contrast approaches as well as
        modern phase-shift imaging strategies. The automated workflow helps
        reduce manual intervention while maintaining consistency across large
        experimental datasets.

        Inputs and Experimental Requirements

        The protocol requires a set of aligned or imported tilt-series
        together with the acquisition metadata describing microscope
        parameters such as voltage, spherical aberration, amplitude contrast,
        and tilt geometry. Reliable metadata are particularly important
        because CTF estimation strongly depends on accurate knowledge of the
        imaging conditions.

        Users may define an expected defocus either as a single global value
        or as a per-tilt-series list. Providing realistic expected defocus
        estimates improves fitting stability, especially for noisy datasets or
        strongly tilted images where Thon rings may be difficult to identify.
        In practice, values close to the microscope acquisition settings
        generally provide the best starting point.

        Defocus Tolerance and Spatial Stability

        One of the key concepts in this protocol is the use of defocus
        tolerance ranges to identify image regions that can be treated as
        having approximately constant defocus. Because tilted specimens
        introduce systematic defocus gradients across the image, the protocol
        analyzes subsets of the field of view that remain sufficiently stable
        for reliable frequency analysis.

        From a biological perspective, choosing appropriate tolerance values
        is important because excessively restrictive settings may reduce the
        amount of usable signal, whereas overly permissive values may combine
        regions with substantially different optical conditions. For most
        routine datasets, moderate default tolerances provide a good balance
        between robustness and sensitivity.

        Angle Range Refinement

        The protocol supports automated fitting across angular ranges of the
        tilt-series. This strategy improves estimation stability by combining
        information from neighboring tilts with related imaging conditions.
        The angular step and angular range parameters determine how broadly
        information is grouped during the estimation process.

        Broad angular ranges can improve robustness for noisy datasets,
        particularly in low-dose tomography experiments. Narrower ranges are
        often preferable when imaging conditions vary rapidly with tilt angle
        or when high precision is required. For datasets with highly stable
        acquisition conditions, fitting individual views may provide more
        localized estimates.

        Frequency Range Selection

        Advanced users can define the frequency region used during fitting.
        This determines which spatial frequencies contribute most strongly to
        the CTF estimation process. Lower frequencies generally provide more
        robust signal, while higher frequencies contain information relevant
        to high-resolution refinement.

        Extending the fitting region toward additional CTF zeros may improve
        accuracy when the signal quality supports it. However, aggressive
        fitting at high frequencies can become unstable in noisy tomographic
        data. Biological users should therefore interpret very high-resolution
        estimates cautiously, particularly in thick specimens or strongly
        tilted images.

        Astigmatism and Phase Shift Analysis

        The protocol optionally supports astigmatism estimation and phase
        shift analysis. These advanced analyses are especially important for
        modern cryo-electron microscopy experiments using phase plates or
        datasets acquired under imperfect optical alignment conditions.

        Astigmatism estimation helps identify anisotropic defocus variations
        that may reduce image quality or distort structural interpretation.
        For biological projects targeting high resolution, monitoring
        astigmatism can significantly improve downstream correction and
        reconstruction quality.

        Phase shift estimation is particularly relevant for phase-plate data,
        where additional optical phase changes modify the appearance of the
        power spectrum. Correctly estimating these parameters improves
        consistency between tilt images and contributes to more accurate
        tomographic reconstructions.

        Streaming and Automated Processing

        The protocol is designed to operate efficiently in streaming
        environments, allowing incoming tilt-series to be processed as they
        become available during acquisition. This capability is especially
        useful in automated facility pipelines and high-throughput cryo-ET
        experiments where rapid quality assessment is essential.

        Continuous processing enables users to identify problematic datasets
        early during data collection. Examples include excessive drift,
        unstable defocus behavior, missing metadata, or poor signal quality.
        Early detection helps optimize microscope usage and reduces the risk
        of investing computational resources in unsuitable datasets.

        Outputs and Interpretation

        The protocol produces a set of CTF estimation results associated with
        each tilt-series. These outputs can subsequently be used for CTF
        correction, tomographic reconstruction, subtomogram averaging, or
        quality control analyses.

        The estimated parameters should always be interpreted in the context
        of the experimental data quality. Smooth and physically reasonable
        defocus trends across tilt angles generally indicate stable
        acquisition conditions. Abrupt discontinuities or highly variable
        estimates may reflect specimen thickness, low signal, charging
        effects, or acquisition instabilities.

        Practical Recommendations

        For routine cryo-electron tomography workflows, it is generally
        advisable to begin with default parameters and visually inspect the
        resulting estimates. If the dataset contains very noisy images or
        strong tilt-induced gradients, adjusting the defocus tolerance and
        angular fitting parameters may improve stability.

        When processing phase-plate datasets or experiments targeting
        near-atomic resolution, enabling astigmatism and phase shift analysis
        is often beneficial. However, these advanced analyses typically
        require higher signal quality and may increase computational cost.

        Providing accurate expected defocus values greatly improves robustness,
        particularly in automated pipelines processing heterogeneous datasets.
        Biological users should also ensure that acquisition metadata are
        complete and consistent before starting large-scale processing.

        Final Perspective

        Automatic CTF estimation is a foundational step in cryo-electron
        tomography because it defines how accurately microscope-induced image
        distortions can be characterized and corrected. Reliable estimation
        improves contrast recovery, enhances structural interpretability, and
        supports higher-quality tomographic reconstructions. Careful parameter
        selection combined with realistic expectations about data quality is
        essential for obtaining biologically meaningful results.
    """

    _label = 'CTF estimation (auto)'
    _possibleOutputs = {OUTPUT_CTF_SERIE: SetOfCTFTomoSeries}
    _interactiveMode = False
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.sRate = None
        self.acq = None

    @classmethod
    def worksInStreaming(cls):
        return True
    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        super().addInTsSetFormParam(form)
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

            form.addParallelSection(threads=3, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def stepsGeneratorStep(self) -> None:
        self._initialize()
        closeSetStepDeps = []
        expDefoci = self.getExpectedDefocus()
        inTsSet = self.getInputTsSet()
        outCtfSet = getattr(self, OUTPUT_CTF_SERIE, None)
        self.readingOutput(outCtfSet)

        while True:
            with self._lock:
                inTsIds = set(inTsSet.getTSIds())

            if not inTsSet.isStreamOpen() and Counter(self.tsIdReadList) == Counter(inTsIds):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetsStep,
                                         OUTPUT_CTF_SERIE,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break

            nonProcessedTsIds = inTsIds - set(self.tsIdReadList)
            tsToProcessDict = {tsId: ts.clone() for ts in inTsSet.iterItems()
                               if (tsId := ts.getTsId()) in nonProcessedTsIds  # Only not processed tsIds
                               and ts.getSize() > 0}  # Avoid processing empty TS
            for tsId, ts in tsToProcessDict.items():
                convId = self._insertFunctionStep(self.convertInputStep,
                                                  ts,
                                                  prerequisites=[],
                                                  needsGPU=False)
                compId = self._insertFunctionStep(self.ctfEstimation,
                                                  ts,
                                                  expDefoci,
                                                  prerequisites=convId,
                                                  needsGPU=False)
                outId = self._insertFunctionStep(self.createOutputStep,
                                                 ts,
                                                 prerequisites=compId,
                                                 needsGPU=False)
                closeSetStepDeps.append(outId)
                logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                self.tsIdReadList.append(tsId)

            self.refreshStreaming(inTsSet)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        tsSet = self.getInputTsSet()
        self.sRate = tsSet.getSamplingRate()
        self.acq = tsSet.getAcquisition()

    def ctfEstimation(self,
                      ts: TiltSeries,
                      expDefoci: dict = None):
        """Run ctfplotter IMOD program"""
        tsId = ts.getTsId()
        if tsId in self.failedItems:
            return
        try:
            logger.info(cyanStr(f'tsId = {tsId} -> Estimating the CTF...'))
            paramsCtfPlotter = self._genCtfPlotterParams(ts, expDefoci)
            self.runProgram(CTFPLOTTER_PROGRAM, paramsCtfPlotter)
        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {CTFPLOTTER_PROGRAM} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def createOutputStep(self, ts: TiltSeries, outputSetName=OUTPUT_CTF_SERIE):
        tsId = ts.getTsId()
        if tsId in self.failedItems:
            self.addToOutFailedSet(ts)
            return

        try:
            defocusFilePath = self.getExtraOutFile(tsId, ext=DEFOCUS_EXT)
            if exists(defocusFilePath):
                self._registerOutput(ts, defocusFilePath, outputSetName=outputSetName)
            else:
                logger.error(redStr(f'tsId = {tsId} -> Output file {defocusFilePath} was not generated. Skipping... '))

        except Exception as e:
           logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
           logger.error(traceback.format_exc())

    @retry_on_sqlite_lock(log=logger)
    def _registerOutput(self,
                        ts: TiltSeries,
                        defocusFilePath: str,
                        outputSetName: str):
        tsId = ts.getTsId()
        with self._lock:
            outputCtfSet = self.getOutputSetOfCTFTomoSeries(self.getInputTsSet(pointer=True),
                                                            outputSetName)
            outCtfSeries = CTFTomoSeries(tsId=tsId)
            outCtfSeries.copyInfo(ts)
            outCtfSeries.setTiltSeries(ts)

            # flags below will be updated in parseTSDefocusFile
            outCtfSeries.setIMODDefocusFileFlag(1)
            outCtfSeries.setNumberOfEstimationsInRange(0)
            outputCtfSet.append(outCtfSeries)
            ctfParser = ImodCtfParser(self)
            ctfParser.parseTSDefocusFile(ts, defocusFilePath, outCtfSeries)

            outCtfSeries.write()
            outputCtfSet.update(outCtfSeries)
            outputCtfSet.write()
            self._store(outputCtfSet)
            # Close explicitly the outputs (for streaming)
            self.closeOutputsForStreaming()

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        ctfSeries = getattr(self, OUTPUT_CTF_SERIE, None)
        if ctfSeries is not None:
            summary.append(f"Input tilt-series: {self.getInputTsSet().getSize()}\n"
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

    def getExpectedDefocus(self) -> Union[dict, None]:
        if self.expectedDefocusOrigin.get() == 1:
            with open(self.expectedDefocusFile.get()) as f:
                lines = f.readlines()
                lines = filter(lambda x: x.strip(), lines)
            result = {line.split()[0]: line.split()[1] for line in lines}
            return result
        else:
            return None

    def _genCtfPlotterParams(self,
                             ts: TiltSeries,
                             expDefoci: dict = None) -> dict:
        tsId = ts.getTsId()
        ctfPlotterParams = {
            "-InputStack": self.getTmpOutFile(tsId),
            "-AngleFile": self.getExtraOutFile(tsId, ext=TLT_EXT),
            "-DefocusFile": self.getExtraOutFile(tsId, ext=DEFOCUS_EXT),
            "-AxisAngle": ts.getAcquisition().getTiltAxisAngle(),
            "-PixelSize": self.sRate / 10,
            "-Voltage": int(self.acq.getVoltage()),
            "-SphericalAberration": self.acq.getSphericalAberration(),
            "-AmplitudeContrast": self.acq.getAmplitudeContrast(),
            "-DefocusTol": self.defocusTol.get(),
            "-PSResolution": 101,
            "-LeftDefTol": self.leftDefTol.get(),
            "-RightDefTol": self.rightDefTol.get(),
            "-tileSize": self.tileSize.get()
        }
        self._addExpectedDefocus(tsId, ctfPlotterParams, expDefoci)
        self._addExcludedViews(ts, ctfPlotterParams)
        self._addAdvancedParams(ctfPlotterParams)
        return ctfPlotterParams

    def _addExpectedDefocus(self,
                            tsId: str,
                            paramsDict: dict,
                            expDefoci: dict) -> None:
        if self.expectedDefocusOrigin.get() == 0:
            paramsDict["-ExpectedDefocus"] = self.expectedDefocusValue.get()
        else:
            defocus = expDefoci.get(tsId)
            if defocus is None:
                raise ValueError(f"{tsId} not found in the provided defocus file.")
            paramsDict["-ExpectedDefocus"] = float(defocus)

    @staticmethod
    def _addExcludedViews(ts: TiltSeries,
                          paramsDict: dict) -> None:
        excludedViews = ts.getTsExcludedViewsIndices(ts.getTsPresentAcqOrders())
        if excludedViews:
            logger.info(cyanStr(f'tsId = {ts.getTsId()} -> Excluded views detected {excludedViews}'))
            paramsDict["-ViewsToSkip"] = ",".join(map(str, excludedViews))

    def _addAdvancedParams(self, paramsDict: dict) -> None:
        if self._interactiveMode:
            paramsDict["-AngleRange"] = "-20.0,20.0"
        else:
            if self.extraZerosToFit.get() != 0:
                paramsDict["-ExtraZerosToFit"] = self.extraZerosToFit.get()
            if self.skipAstigmaticViews:
                paramsDict["-SkipOnlyForAstigPhase"] = ""
            if self.searchAstigmatism:
                paramsDict.update({
                    "-SearchAstigmatism": "",
                    "-MaximumAstigmatism": self.maximumAstigmatism.get(),
                    "-NumberOfSectors": self.numberSectorsAstigmatism.get()
                })
            if self.searchPhaseShift:
                paramsDict["-SearchPhaseShift"] = ""
            self._addMinViewsAstigAndPhase(paramsDict)
            if self.searchCutOnFreq:
                paramsDict["-SearchCutonFrequency"] = ""
                if self.maximumCutOnFreq.get() != -1.0:
                    paramsDict["-MaxCutOnToSearch"] = self.maximumCutOnFreq.get()
            paramsDict["-AutoFitRangeAndStep"] = f"{self.angleRange.get()},{self.angleStep.get()}"
            if self.startFreq.get() != 0 or self.endFreq.get() != 0:
                paramsDict["-FrequencyRangeToFit"] = f"{self.startFreq.get()},{self.endFreq.get()}"
            paramsDict["-SaveAndExit"] = ""

    def _addMinViewsAstigAndPhase(self, paramsDict: dict) -> None:
        astig = self.searchAstigmatism
        phase = self.searchPhaseShift
        if astig and phase:
            paramsDict["-MinViewsAstigAndPhase"] = (f"{self.minimumViewsAstigmatism.get()},"
                                                    f"{self.minimumViewsPhaseShift.get()}")
        elif astig:
            paramsDict["-MinViewsAstigAndPhase"] = f"{self.minimumViewsAstigmatism.get()},0"
        elif phase:
            paramsDict["-MinViewsAstigAndPhase"] = f"0,{self.minimumViewsPhaseShift.get()}"
