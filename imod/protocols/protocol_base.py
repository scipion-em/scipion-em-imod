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
import time
import traceback
import typing
from typing import Union, Tuple, List
from imod.constants import NEWSTACK_PROGRAM, ODD, EVEN, XF_EXT, TLT_EXT, OUTPUT_TILTSERIES_NAME, \
    OUTPUT_FIDUCIAL_NO_GAPS_NAME, OUTPUT_TOMOGRAMS_NAME, OUTPUT_TS_FAILED_NAME, OUTPUT_TOMOS_FAILED_NAME, MRC_EXT, \
    MRCS_EXT
from imod.convert.convert import genXfFile
from pyworkflow.object import Set, Boolean, Pointer
from pyworkflow.protocol import params
from pyworkflow.utils import path, cyanStr, redStr, yellowStr
from pwem.protocols import EMProtocol
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.protocols.protocol_base import ProtTomoBase
from tomo.objects import (SetOfTiltSeries, SetOfTomograms, SetOfCTFTomoSeries,
                          TiltSeries, TiltImage, SetOfLandmarkModels, Tomogram)
from imod import Plugin



logger = logging.getLogger(__name__)
IN_TS_SET = 'inputSetOfTiltSeries'
IN_TOMO_SET = 'inputSetOfTomograms'
IN_CTF_TOMO_SET = 'inputSetOfCtfTomoSeries'
PROCESS_ODD_EVEN = 'processOddEven'
BINNING_FACTOR = 'binning'


class ProtImodBase(EMProtocol, ProtTomoBase):
    """
    Base protocol for IMOD tomography workflows. It provides the common
    infrastructure required to prepare, organize, process, and export
    tilt-series and tomographic datasets within Scipion environments.

    AI Generated:

    IMOD Base Protocol (ProtImodBase) - User Manual
        Overview

        The IMOD Base protocol serves as the foundational framework for
        tomography-oriented processing workflows based on IMOD utilities.
        Its purpose is to standardize how tilt-series, tomograms,
        fiducial models, and CTF information are prepared and managed
        throughout cryo-electron tomography pipelines. Rather than
        representing a standalone biological analysis procedure, this
        protocol establishes the operational environment that allows
        higher-level tomography methods to function consistently and
        reliably.

        In practical cryo-ET workflows, preprocessing and data handling
        are often as important as the reconstruction or alignment stages
        themselves. Datasets may contain excluded views, prior alignment
        information, varying acquisition orders, or associated odd/even
        datasets used for resolution estimation and validation. This
        protocol provides a unified strategy for handling all of these
        situations while maintaining compatibility with IMOD processing
        conventions.

        General Role in Tomography Pipelines

        The protocol is designed to support workflows involving tilt-series
        alignment, tomogram reconstruction, fiducial tracking, and CTF
        estimation. It manages the preparation of input data before
        execution of downstream IMOD programs and ensures that outputs are
        consistently registered into Scipion project structures.

        From a biological perspective, this organizational layer is
        essential because tomography datasets are rarely uniform. Different
        acquisition sessions may include missing tilts, interrupted angular
        coverage, partially aligned data, or heterogeneous preprocessing
        histories. Proper management of these conditions is necessary to
        avoid reconstruction artifacts, geometric inconsistencies, or loss
        of biologically meaningful structural information.

        Tilt-Series Preparation and Alignment Handling

        One of the central responsibilities of the protocol is the
        preparation of tilt-series for downstream IMOD operations.
        Depending on the characteristics of the dataset, the protocol can
        work with already aligned tilt-series, unaligned data, or datasets
        containing excluded projections.

        In workflows where alignment information already exists, the
        protocol preserves and propagates the geometric transformations so
        that downstream reconstruction or refinement steps operate in a
        consistent spatial framework. This is particularly important when
        combining processing stages originating from different software
        environments or previous Scipion executions.

        When excluded views are present, the protocol ensures that the
        angular metadata and associated image transformations remain
        synchronized. Biologically, maintaining correct angular ordering is
        critical because tomographic reconstruction accuracy strongly
        depends on the consistency between projection geometry and image
        content. Mishandling excluded projections can introduce distortions
        that compromise downstream interpretation of cellular structures or
        macromolecular assemblies.

        Support for Odd and Even Data

        The protocol also supports workflows involving odd and even
        datasets. These half datasets are commonly used in cryo-ET to
        estimate resolution, validate reconstructions, and assess
        overfitting. When enabled, the same processing operations applied
        to the full dataset are propagated consistently to the odd and even
        subsets.

        From a biological interpretation standpoint, preserving identical
        geometric treatment between half datasets is essential for reliable
        Fourier shell correlation measurements and for obtaining meaningful
        quality estimates in subtomogram averaging or tomographic
        reconstruction workflows.

        Streaming and Incremental Processing

        The protocol includes infrastructure for streaming-oriented
        processing. In high-throughput cryo-ET facilities, tilt-series may
        arrive progressively during data acquisition rather than being
        available as a complete dataset from the start. The protocol
        supports incremental handling of these datasets while maintaining
        stable bookkeeping of processed and pending items.

        This capability becomes especially valuable during automated
        acquisition campaigns or facility-scale processing pipelines where
        rapid feedback is needed to evaluate sample quality, alignment
        stability, or acquisition performance before the microscopy session
        has finished.

        Management of Multiple Output Types

        The protocol standardizes the generation of several important
        tomography outputs, including processed tilt-series, reconstructed
        tomograms, fiducial models, failed datasets, and CTF estimation
        series. Each output type is organized into dedicated collections so
        that downstream protocols can consume them transparently.

        Failed datasets are also tracked explicitly. In large cryo-ET
        experiments, some tilt-series may fail due to acquisition issues,
        excessive drift, poor contrast, or inconsistent metadata. Keeping
        these datasets separated from successful outputs allows users to
        continue processing valid data while preserving traceability and
        reproducibility.

        Integration with IMOD Programs

        The protocol provides a unified interface for executing IMOD
        utilities, particularly workflows based on newstack and related
        geometric operations. These operations include interpolation,
        restacking, binning, tapering, normalization, and handling of
        excluded projections.

        Binning support is especially relevant in tomography because large
        datasets can become computationally demanding. Reduced-resolution
        datasets are often used during exploratory alignment stages or
        rapid quality-control analyses before final high-resolution
        reconstruction. Care must be taken, however, because aggressive
        binning may remove biologically important fine structural details.

        Spatial Consistency and Coordinate Preservation

        Throughout all processing stages, the protocol emphasizes
        preservation of geometric consistency between projections,
        transformations, angular metadata, and reconstructed outputs. This
        is biologically important because cryo-ET interpretation depends
        directly on accurate spatial relationships between projections and
        reconstructed density.

        Errors in coordinate propagation can lead to blurred tomograms,
        distorted fiducial trajectories, inaccurate particle localization,
        or unreliable averaging results. By maintaining coherent metadata
        management across processing stages, the protocol reduces the risk
        of introducing subtle geometric inconsistencies.

        Practical Recommendations

        In routine cryo-ET workflows, users should ensure that input
        tilt-series contain consistent metadata before processing. Angular
        ordering, excluded view annotations, and alignment information
        should be verified carefully, particularly when importing data from
        external software packages.

        Odd/even processing should generally be enabled when downstream
        resolution validation or subtomogram averaging is anticipated.
        Similarly, binning should initially be kept moderate in order to
        preserve structural interpretability during exploratory analyses.

        For datasets containing excluded projections or incomplete tilt
        ranges, users should inspect reconstruction quality carefully
        because missing angular information may introduce anisotropic
        resolution or elongation artifacts.

        Final Perspective

        The IMOD Base protocol provides the structural backbone required
        for reliable tomography processing within Scipion. Although much of
        its functionality is organizational and infrastructural, these
        operations directly influence the biological quality and
        reproducibility of downstream cryo-ET analyses.

        Proper management of tilt geometry, alignment consistency,
        odd/even datasets, and reconstruction metadata is essential for
        generating tomograms that can support trustworthy structural and
        cellular interpretation.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsDict = None
        self.tomoDict = None
        self.failedItems = []
        self.doOddEven = False

        # Streaming
        self.tsIdReadList = []

    # -------------------------- DEFINE param functions -----------------------
    @staticmethod
    def addInTsSetFormParam(form):
        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt Series')

    @staticmethod
    def addInTomoSetFormParam(form):
        form.addParam(IN_TOMO_SET,
                      params.PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Input set of tomograms')

    @staticmethod
    def addOddEvenParams(form, isTomogram=False):
        objStr = 'tomograms' if isTomogram else 'tilt-series'
        form.addParam(PROCESS_ODD_EVEN,
                      params.BooleanParam,
                      default=False,
                      label='Apply to odd/even',
                      help=f'If True, the full {objStr} and the associated odd/even '
                           f'{objStr} will be processed. The transformations applied '
                           f'to the odd/even {objStr} will be exactly the same.')

    @classmethod
    def worksInStreaming(cls):
        """ So far none of them work in streaming. """
        return False

    def allowsDelete(self, obj):
        return True

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.doOddEven = self.applyToOddEven(self.getInputTsSet())

    def refreshStreaming(self,
                         inSet: Union[SetOfTiltSeries, SetOfLandmarkModels, SetOfCTFTomoSeries]):
        # Refresh status for the streaming
        time.sleep(10)
        if inSet.isStreamOpen():
            with self._lock:
                inSet.loadAllProperties()  # refresh status for the streaming

    def closeOutputsForStreaming(self):
        # Close explicitly the outputs (for streaming)
        for outputName in self._possibleOutputs.keys():
            output = getattr(self, outputName, None)
            if output:
                output.close()

    def linkTsStep(self, ts: TiltSeries) -> None:
        try:
            self._linkTs(ts)
        except Exception as e:
            logger.error(redStr(f'tsId = {ts.getTsId()} -> input conversion failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def _linkTs(self, ts: TiltSeries):
        tsId = ts.getTsId()
        self.genTsPaths(tsId)
        outTsFn = self.getTmpOutFile(tsId)
        firstTi = ts.getFirstEnabledItem()
        # Make the link using the tsId instead of the original name prevent IMOD from
        # failing in case of strange characters or even numeric names
        tsFn = firstTi.getFileName()
        logger.info(cyanStr(f"tsId = {tsId}: link TS: {outTsFn} -> {tsFn}"))
        self.linkTs(tsFn, outTsFn)
        if self.doOddEven:
            # ODD
            inTsOddFn = firstTi.getOdd()
            outTsFnOdd = self.getTmpOutFile(tsId, suffix=ODD)
            logger.info(cyanStr(f"tsId = {tsId}: link TS ODD: {outTsFnOdd} -> {inTsOddFn}"))
            self.linkTs(inTsOddFn, outTsFnOdd)
            # Even
            inTsEvenFn = firstTi.getEven()
            outTsFnEven = self.getTmpOutFile(tsId, suffix=EVEN)
            self.linkTs(inTsEvenFn, outTsFnEven)
            logger.info(cyanStr(f"tsId = {tsId}: link TS EVEN: {outTsFnEven} -> {inTsEvenFn}"))

    def convertInputStep(self,
                         ts: TiltSeries,
                         presentAcqOrders: typing.Optional[typing.Set[int]] = None):
        tsId = ts.getTsId()
        try:
            self.genTsPaths(tsId)
            if presentAcqOrders:
                self.convertInputForNonEvProgram(ts, presentAcqOrders)
            else:
                self.convertInputForEvProgram(ts)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> input conversion failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def convertInputForEvProgram(self, ts: TiltSeries) -> None:
        """The input converters of the protocols that use main IMOD programs (excluding newstack,
        considered here as an auxiliary program) that can manage the excluded views must behave as
        described below to work as expected:

            -> Excluded views must be ignored in the generation of the tlt and, if the tilt-series
               has alignment, in the generation of the xf. The program newstack is only executed in
               case of alignment for interpolation.

        These programs are the ones used in the coarse pre-alignment, fiducial model and fiducial alignment.
        """
        tsId = ts.getTsId()
        if ts.hasAlignment():
            logger.info(f"tsId = {tsId}: alignment will be applied with {NEWSTACK_PROGRAM}")
            xfFile = self.getExtraOutFile(tsId, ext=XF_EXT)
            # The xf file must contain all thw views to interpolate and re-stack
            genXfFile(ts, xfFile)
            self.runNewStackBasic(ts, xfFile=xfFile)
        else:
            # Link it, so the input file expected is in the same place in both sides of the "if"
            self._linkTs(ts)

        # Generate the tlt file
        tltFile = self.getExtraOutFile(tsId, ext=TLT_EXT)
        ts.generateTltFile(tltFile)

    def convertInputForNonEvProgram(self,
                                    ts: TiltSeries,
                                    presentAcqOrders: typing.Set[int]) -> None:
        """The input converters of the protocols that use main IMOD programs (excluding newstack,
          considered here as an auxiliary program) that can manage the excluded views must behave as
          described below to work as expected:

              -> Excluded views must be considered in the generation of the tlt.
              -> If the tilt-series has alignment, a first xf file must be generated containing all
                 the views and the program newstack must be executed with the whole xf file and passing
                 the present acquisition orders. This will generate an interpolated and re-stacked
                 tilt-series. After that, a new xf file will be generated containing only the
                 transformations that correspond to the non-excluded views to be used by the protocol's
                 main program.
              -> If the tilt-series does not have alignment, the acquisition orders will be used for
                 re-stacking the tilt-series.
          """
        tsId = ts.getTsId()
        hasExcludedViews = ts.hasExcludedViews()
        hasAlignment = ts.hasAlignment()
        if not hasAlignment and not hasExcludedViews:
            # Link it, so the input file expected is in the same place in both sides of the "if"
            self._linkTs(ts)
        else:
            if hasAlignment:
                xfFile = self.getTmpOutFile(ts.getTsId(), ext=XF_EXT)
                try:
                    logger.info(f"tsId = {tsId}: alignment will be applied with {NEWSTACK_PROGRAM}")
                    # The xf file must contain all the views to make newstack interpolate and
                    # re-stack using its own excluded views feature
                    genXfFile(ts, xfFile)
                    self.runNewStackBasic(ts,
                                          xfFile=xfFile,
                                          presentAcqOrders=presentAcqOrders)
                except Exception as e:
                    # In some cases, newstack may fail (e.g. the assigning the transformation matrix from one
                    # tilt-series with smaller number of elements to a bigger one). In that case, newstack fails
                    # because it is not prepared to manage a tilt-series binary file with more tilt-images than
                    # lines in the alignment file, but Scipion can manage that case
                    logger.info(yellowStr(f'tsId = {tsId} - program {NEWSTACK_PROGRAM} failed with the exception '
                                          f'{e}'))
                    logger.info(cyanStr(f'Trying with Scipion...'))
                    outTsFn, _, _ = self.getTmpFileNames(tsId)
                    ts.applyTransform(outTsFn)

                # After that, for the following programs, a new xfFile without the
                # excluded views must be generated to be used by the protocol main program
                xfFile = self.getExtraOutFile(ts.getTsId(), ext=XF_EXT)
                genXfFile(ts, xfFile, presentAcqOrders=presentAcqOrders)
            else:
                # Re-stack
                logger.info(f"tsId = {tsId}: tilt-series re-stacking will be carried out with {NEWSTACK_PROGRAM}")
                self.runNewStackBasic(ts, presentAcqOrders=presentAcqOrders)

        # Generate the tlt file without the excluded views must be generated to be
        # used by the protocol main program
        tltFile = self.getExtraOutFile(tsId, ext=TLT_EXT)
        ts.generateTltFile(tltFile, presentAcqOrders=presentAcqOrders)

    def closeOutputSetsStep(self, attrib: Union[List[str], str]):
        self._closeOutputSet()
        attribList = [attrib] if type(attrib) is str else attrib
        failedOutputList = []
        for attr in attribList:
            output = getattr(self, attr, None)
            if not output or (output and len(output) == 0):
                failedOutputList.append(attr)
        if failedOutputList:
            raise Exception(f'No output/s {failedOutputList} were generated. Please check the '
                            f'Output Log > run.stdout and run.stderr')

    # --------------------------- OUTPUT functions ----------------------------
    def getOutputSetOfTS(self,
                         inputPtr,
                         binning=1,
                         attrName=OUTPUT_TILTSERIES_NAME,
                         tiltAxisAngle=None,
                         suffix="") -> SetOfTiltSeries:
        """ Method to generate output of set of tilt-series.
        :param inputPtr: input set pointer (TS or tomograms)
        :param binning: binning factor
        :param attrName: output attr name
        :param tiltAxisAngle: Only applies to TS. If not None, the corresponding value of the
        set acquisition will be updated (xCorr prot)
        :param suffix: output set suffix
        """
        inputSet = inputPtr.get()
        outputSet = getattr(self, attrName, None)
        if outputSet:
            outputSet.enableAppend()
        else:
            outputSet = self._createSetOfTiltSeries(suffix=suffix)

            if isinstance(inputSet, SetOfTiltSeries):
                outputSet.copyInfo(inputSet)
                if tiltAxisAngle:
                    outputSet.getAcquisition().setTiltAxisAngle(tiltAxisAngle)

            elif isinstance(inputSet, SetOfTomograms):
                # _anglesCount does not exist for SetOfTomograms, so we can't use copyInfo
                outputSet.setAcquisition(inputSet.getAcquisition())
                outputSet.setSamplingRate(inputSet.getSamplingRate())

            if binning > 1:
                samplingRate = inputSet.getSamplingRate()
                samplingRate *= binning
                outputSet.setSamplingRate(samplingRate)

            outputSet.setStreamState(Set.STREAM_OPEN)

            # Write set properties, otherwise it may expose the set (sqlite) without properties.
            outputSet.write()

            self._defineOutputs(**{attrName: outputSet})
            self._defineSourceRelation(inputPtr, outputSet)

        return outputSet

    def getOutputSetOfCTFTomoSeries(self,
                                    inputPtr: Pointer,
                                    outputSetName: str) -> SetOfCTFTomoSeries:

        outputSetOfCTFTomoSeries = getattr(self, outputSetName, None)

        if outputSetOfCTFTomoSeries is not None:
            outputSetOfCTFTomoSeries.enableAppend()
        else:
            outputSetOfCTFTomoSeries = SetOfCTFTomoSeries.create(self._getPath(),
                                                                 template='CTFmodels%s.sqlite')
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(inputPtr)
            outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputSetName: outputSetOfCTFTomoSeries})
            self._defineCtfRelation(inputPtr, outputSetOfCTFTomoSeries)

        return outputSetOfCTFTomoSeries

    def getOutputFiducialModel(self,
                               inputPtr: Pointer,
                               attrName: str = OUTPUT_FIDUCIAL_NO_GAPS_NAME,
                               suffix: str = "NoGaps",
                               forceNew=False) -> SetOfLandmarkModels:
        """ Method to generate output of set fiducial models.
                :param inputPtr: input TS set pointer
                :param attrName: output attr name
                :param suffix: output set suffix
                :param forceNew: Forces to have a new output even if exists
        """
        if not inputPtr.isPointer():
            logger.warning("FOR DEVELOPERS: inputSet must be a pointer!")
            inputSet = inputPtr
        else:
            inputSet = inputPtr.get()

        fidModel = getattr(self, attrName, None)

        if fidModel is not None and not forceNew:
            fidModel.enableAppend()
        else:
            fidModel = self._createSetOfLandmarkModels(suffix=suffix)
            fidModel.copyInfo(inputSet)
            fidModel.setSetOfTiltSeries(inputPtr)
            fidModel.setHasResidualInfo(True)
            fidModel.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{attrName: fidModel})
            self._defineSourceRelation(inputPtr, fidModel)

        return fidModel

    def getOutputSetOfTomograms(self,
                                inputPtr: Pointer,
                                binning: int = 1) -> SetOfTomograms:
        inputSet = inputPtr.get()
        outputSet = getattr(self, OUTPUT_TOMOGRAMS_NAME, None)

        if outputSet:
            outputSet.enableAppend()
        else:
            outputSet = self._createSetOfTomograms()
            outputSet.copyInfo(inputSet)

            if binning > 1:
                samplingRate = inputSet.getSamplingRate() * binning
                outputSet.setSamplingRate(samplingRate)

            outputSet.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_TOMOGRAMS_NAME: outputSet})
            self._defineSourceRelation(inputPtr, outputSet)

        return outputSet

    def getOutputFailedSet(self,
                           inputPtr: Pointer,
                           inputsAreTs: bool = True) -> Union[SetOfTiltSeries, SetOfTomograms]:
        """ Create output set for failed TS or tomograms. """
        inputSet = inputPtr.get()
        if inputsAreTs:
            failedTs = getattr(self, OUTPUT_TS_FAILED_NAME, None)

            if failedTs:
                failedTs.enableAppend()
            else:
                logger.info(cyanStr('Create the set of failed TS'))
                failedTs = self._createSetOfTiltSeries(suffix='Failed')
                failedTs.copyInfo(inputSet)
                failedTs.setStreamState(Set.STREAM_OPEN)
                self._defineOutputs(**{OUTPUT_TS_FAILED_NAME: failedTs})
                self._defineSourceRelation(inputPtr, failedTs)

            return failedTs

        else:
            failedTomos = getattr(self, OUTPUT_TOMOS_FAILED_NAME, None)
            if failedTomos:
                failedTomos.enableAppend()
            else:
                logger.info(cyanStr('Create the set of failed tomograms'))
                failedTomos = self._createSetOfTomograms(suffix='Failed')
                failedTomos.copyInfo(inputSet)
                failedTomos.setStreamState(Set.STREAM_OPEN)
                self._defineOutputs(**{OUTPUT_TOMOS_FAILED_NAME: failedTomos})
                self._defineSourceRelation(inputPtr, failedTomos)

            return failedTomos

    @retry_on_sqlite_lock(log=logger)
    def addToOutFailedSet(self,
                          item: Union[TiltSeries, Tomogram]) -> None:
        """ Just copy input item to the failed output set. """
        tsId = item.getTsId()
        logger.info(cyanStr(f'Failed TS ---> {tsId}'))
        try:
            inputsAreTs = True if isinstance(item, TiltSeries) else False
            with self._lock:
                inputSet = self.getInputTsSet(pointer=True) if inputsAreTs else self.getInputTomoSet(pointer=True)
                output = self.getOutputFailedSet(inputSet, inputsAreTs=inputsAreTs)
                newItem = item.clone()
                newItem.copyInfo(item)
                output.append(newItem)

                if inputsAreTs:
                    newItem.copyItems(item)
                    newItem.write()

                output.update(newItem)
                output.write()
                self._store(output)
                # Close explicitly the outputs (for streaming)
                output.close()
        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the failed output with '
                                f'exception {e}. Skipping... '))

    def setTsOddEven(self, tsId: str, outTi: TiltImage, binGenerated: bool = False) -> None:
        if self.doOddEven:
            if binGenerated:
                outTi.setOddEven([self.getExtraOutFile(tsId, suffix=ODD),
                                  self.getExtraOutFile(tsId, suffix=EVEN)])
        else:
            outTi.setOddEven([])  # the input may have odd/even but the user may have decided not
            # to consider them in the current execution, so they should be set to empty to avoid
            # next protocols be confused about having them.

    def setTomoOddEven(self, tsId: str, outTomo: Tomogram) -> None:
        if self.doOddEven:
            halfMapsList = [self.getExtraOutFile(tsId, suffix=ODD, ext=MRC_EXT),
                            self.getExtraOutFile(tsId, suffix=EVEN, ext=MRC_EXT)]
            outTomo.setHalfMaps(halfMapsList)
        else:
            outTomo.setHalfMaps([])

    # --------------------------- UTILS functions -----------------------------
    def getInputTsSet(self, pointer: bool = False) -> Union[Pointer, SetOfTiltSeries]:
        tsSetPointer = getattr(self, IN_TS_SET)
        return tsSetPointer if pointer else tsSetPointer.get()

    def getInputTomoSet(self, pointer: bool = False) -> Union[Pointer, SetOfTomograms]:
        tomoSetPointer = getattr(self, IN_TOMO_SET)
        return tomoSetPointer if pointer else tomoSetPointer.get()

    def getInputCtfSet(self, pointer: bool = False) -> Union[Pointer, SetOfCTFTomoSeries]:
        tomoSetPointer = getattr(self, IN_CTF_TOMO_SET)
        return tomoSetPointer if pointer else tomoSetPointer.get()

    def genTsPaths(self, tsId):
        """Generate the subdirectories corresponding to the
        current tilt-series in tmp and extra"""
        path.makePath(*[self._getExtraPath(tsId), self._getTmpPath(tsId)])

    def readingOutput(self,
                      outSet: Union[SetOfTiltSeries, SetOfTomograms,
                      SetOfLandmarkModels, SetOfCTFTomoSeries],
                      tsIdListName: str = None) -> None:
        if outSet:
            if tsIdListName:
                tsIdList = getattr(self, tsIdListName)
            else:
                tsIdList = self.tsIdReadList
            for item in outSet:
                tsIdList.append(item.getTsId())
            self.info(cyanStr(f'Item processed {tsIdList}'))
        else:
            self.info(cyanStr('No items have been processed yet'))

    @staticmethod
    def getOutTsFileName(tsId, suffix=None, ext=MRCS_EXT):
        return f'{tsId}_{suffix}.{ext}' if suffix else f'{tsId}.{ext}'

    def getTmpOutFile(self, tsId, suffix=None, ext=MRCS_EXT):
        return self._getTmpPath(tsId,
                                self.getOutTsFileName(tsId, suffix=suffix, ext=ext))

    def getExtraOutFile(self, tsId, suffix=None, ext=MRCS_EXT):
        return self._getExtraPath(tsId,
                                  self.getOutTsFileName(tsId, suffix=suffix, ext=ext))

    def applyToOddEven(self, setOfTs):
        return (hasattr(self, "processOddEven") and
                self.processOddEven.get() and
                setOfTs.hasOddEven())

    def warningOddEven(self, inSet: Union[SetOfTiltSeries, SetOfTomograms], warnMsgList: list):
        if getattr(self, PROCESS_ODD_EVEN, Boolean(False).get()) and not inSet.hasOddEven():
            warnMsgList.append('The even/odd tilt-series or tomograms were not found in the introduced tilt-series or '
                               'tomograns metadata. Thus, only the full tilt-series or tomograms will be processed.')

    def runProgram(self, program, params, cwd=None):
        """ Shortcut method to run IMOD's command given input params dict. """
        args = ' '.join(['%s %s' % (k, str(v)) for k, v in params.items()])
        Plugin.runImod(self, program, args, cwd)

    @staticmethod
    def getBasicNewstackParams(ts: Union[TiltSeries, Tomogram],
                               inFileName: str,
                               outFileName: str,
                               xfFile: str = None,
                               binning: int = 1,
                               doSwap: bool = False,
                               tsExcludedIndices: set = None,
                               doTaper: bool = False,
                               doNorm: bool = False):
        """ Returns basic newstack arguments

        :param ts: Title Series or Tomogram object
        :param inFileName: Input tilt-series file name.
        :param outFileName:  Output tilt-series file name.
        :param xfFile: xf file name, if passed, alignment will be generated and used
        :param binning: Default to 1. to apply to output size
        :param doSwap: Default False.
        :param tsExcludedIndices: List of indices to be excluded in the tilt-series, starting from 1
        :param doTaper: optionally taper the tilt-series
        :param doNorm: optionally normalize the tilt-series
        """
        logger.info(cyanStr(f'tsId = {ts.getTsId()} -> Executing {NEWSTACK_PROGRAM}:.'))
        # Apply interpolation
        paramsNs = {
            '-input': inFileName,
            '-output': outFileName,
            '-bin': binning,
            '-antialias': -1,
            '-imagebinned': 1.0
        }
        if doTaper:
            paramsNs["-taper"] = "1,1"
        if doNorm:
            paramsNs["-FloatDensities"] = 2

        if xfFile is not None:
            paramsNs['-xform'] = xfFile
            logger.info(cyanStr(f'\t--> xf file detected. The tilt-series will be interpolated with {xfFile}'))

            if doSwap:
                dimX, dimY, _ = ts.getFirstEnabledItem().getDim()
                paramsNs["-size"] = f"{round(dimY / binning)}," \
                                    f"{round(dimX / binning)}"

        if tsExcludedIndices:
            paramsNs["-exclude"] = ",".join(map(str, tsExcludedIndices))
            # From IMOD's newstack doc: "sections are numbered from 0 unless -fromone is entered"
            paramsNs["-fromone"] = ""
            logger.info(cyanStr(f'\t--> Excluded views detected. ==> {tsExcludedIndices}. '
                                f'The tilt-series passed to the main program/s of the protocol '
                                f'will be re-stacked '))

        if binning > 1:
            logger.info(cyanStr(f'\t--> The tilt-series will be re-scaled to bin {binning}.'))

        return paramsNs

    @staticmethod
    def getNewstackDoSwap(ti: TiltImage, xfFile: str) -> bool:
        doSwap = False
        if xfFile:  # No xf file means that newstack will be used only
            # for scaling or re-stacking, but not for interpolation
            rotationAngle = ti.getRotationAngle()
            doSwap = 45 < abs(rotationAngle) < 135
        return doSwap

    def runNewStackBasic(self,
                         ts: TiltSeries,
                         xfFile: str = None,
                         binning: int = 1,
                         presentAcqOrders: typing.Set[int] = None) -> None:

        tsExcludedIndices = None
        tsId = ts.getTsId()
        outTsFn, outTsOddFn, outTsEvenFn = self.getTmpFileNames(tsId)
        firstTi = ts.getFirstEnabledItem()
        doSwap = self.getNewstackDoSwap(firstTi, xfFile)
        if presentAcqOrders:
            tsExcludedIndices = ts.getTsExcludedViewsIndices(presentAcqOrders)

        logger.info(cyanStr(f'tsId = {tsId}: running {NEWSTACK_PROGRAM}...'))
        param = self.getBasicNewstackParams(ts,
                                            firstTi.getFileName(),
                                            outTsFn,
                                            xfFile=xfFile,
                                            doSwap=doSwap,
                                            tsExcludedIndices=tsExcludedIndices,
                                            binning=binning)
        self.runProgram(NEWSTACK_PROGRAM, param)

        if self.doOddEven:
            # ODD
            logger.info(cyanStr(f'tsId = {tsId} ODD: running {NEWSTACK_PROGRAM}...'))
            param = self.getBasicNewstackParams(ts,
                                                ts.getOddFileName(),
                                                outTsOddFn,
                                                xfFile=xfFile,
                                                doSwap=doSwap,
                                                tsExcludedIndices=tsExcludedIndices,
                                                binning=binning)
            self.runProgram(NEWSTACK_PROGRAM, param)
            # EVEN
            logger.info(cyanStr(f'tsId = {tsId} EVEN: running {NEWSTACK_PROGRAM}...'))
            param = self.getBasicNewstackParams(ts,
                                                ts.getEvenFileName(),
                                                outTsEvenFn,
                                                xfFile=xfFile,
                                                doSwap=doSwap,
                                                tsExcludedIndices=tsExcludedIndices,
                                                binning=binning)
            self.runProgram(NEWSTACK_PROGRAM, param)

    def getTmpFileNames(self, tsId: str) -> Tuple:
        tsFn = self.getTmpOutFile(tsId)
        if self.doOddEven:
            tsFnOdd = self.getTmpOutFile(tsId, suffix=ODD)
            tsFnEven = self.getTmpOutFile(tsId, suffix=EVEN)
            return tsFn, tsFnOdd, tsFnEven
        else:
            return tsFn, None, None

    @staticmethod
    def linkTs(inFileName: str, outFileName: str):
        logger.info(cyanStr("\t--> Tilt-series file linked."))
        path.createAbsLink(inFileName, outFileName)

    # --------------------------- INFO functions ------------------------------
    def _warnings(self):
        warnMsgList = []
        if getattr(self, IN_TS_SET, None):
            inSet = self.getInputTsSet()
            self.warningOddEven(inSet, warnMsgList)
        elif getattr(self, IN_TOMO_SET, None):
            inSet = self.getInputTomoSet()
            self.warningOddEven(inSet, warnMsgList)
        return warnMsgList
