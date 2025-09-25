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
from subprocess import CalledProcessError
from typing import Union, Tuple, List
from imod.convert.convert import genXfFile
from pyworkflow.object import Set, Boolean, Pointer
from pyworkflow.protocol import params
from pyworkflow.utils import path, cyanStr, redStr, yellowStr
from pwem.protocols import EMProtocol
from tomo.protocols.protocol_base import ProtTomoBase
from tomo.objects import (SetOfTiltSeries, SetOfTomograms, SetOfCTFTomoSeries,
                          TiltSeries, TiltImage, CTFTomoSeries,
                          SetOfLandmarkModels, Tomogram)
from imod import Plugin
from imod.constants import *

logger = logging.getLogger(__name__)
IN_TS_SET = 'inputSetOfTiltSeries'
IN_TOMO_SET = 'inputSetOfTomograms'
IN_CTF_TOMO_SET = 'inputSetOfCtfTomoSeries'
PROCESS_ODD_EVEN = 'processOddEven'
BINNING_FACTOR = 'binning'


class ProtImodBase(EMProtocol, ProtTomoBase):
    """ Base class with methods used in the rest of the imod protocols. """

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

    def refreshStreaming(self, inSet: Union[SetOfTiltSeries, SetOfLandmarkModels]):
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

    def linkTsStep(self, tsId: str):
        # Link the tils-series to tmp using the tsId as basename, preventing problematic
        # filenames, such as the ones starting by a number or containing special characters. That
        # won't happen with the tsId as it is normalized
        try:
            self.genTsPaths(tsId)
            with self._lock:
                ts = self.getCurrentTs(tsId)
            self.linkTs(ts.getFirstItem().getFileName(), self.getTmpOutFile(tsId))
        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> input conversion failed with the exception -> {e}'))

    def convertInputStep(self,
                         tsId: str,
                         presentAcqOrders: typing.Optional[typing.Set[int]] = None):
        try:
            self.genTsPaths(tsId)
            if presentAcqOrders:
                self.convertInputForNonEvProgram(tsId, presentAcqOrders)
            else:
                self.convertInputForEvProgram(tsId)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> input conversion failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def convertInputForEvProgram(self, tsId: str) -> None:
        """The input converters of the protocols that use main IMOD programs (excluding newstack,
        considered here as an auxiliary program) that can manage the excluded views must behave as
        described below to work as expected:

            -> Excluded views must be ignored in the generation of the tlt and, if the tilt-series
               has alignment, in the generation of the xf. The program newstack is only executed in
               case of alignment for interpolation.

        These programs are the ones used in the coarse pre-alignment, fiducial model and fiducial alignment.
        """
        with self._lock:
            ts = self.getCurrentTs(tsId)
            firstTi = ts.getFirstItem()

        hasAlignment = firstTi.hasTransform()
        if hasAlignment:
            xfFile = self.getExtraOutFile(ts.getTsId(), ext=XF_EXT)
            # The xf file must contain all thw views to interpolate and re-stack
            genXfFile(ts, xfFile)
            self.runNewStackBasic(ts, xfFile=xfFile)
        else:
            # Link it, so the input file expected is in the same place in both sides of the "if"
            outTsFn, _, _ = self.getTmpFileNames(ts)
            self.linkTs(firstTi.getFileName(), outTsFn)

        # Generate the tlt file
        tltFile = self.getExtraOutFile(tsId, ext=TLT_EXT)
        ts.generateTltFile(tltFile)

    def convertInputForNonEvProgram(self,
                                    tsId: str,
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
        with self._lock:
            ts = self.getCurrentTs(tsId)

        firstTi = ts.getFirstItem()
        hasExcludedViews = ts.hasExcludedViews()
        hasAlignment = firstTi.hasTransform()
        if not hasAlignment and not hasExcludedViews:
            # Link it, so the input file expected is in the same place in both sides of the "if"
            outTsFn, _, _ = self.getTmpFileNames(ts)
            self.linkTs(firstTi.getFileName(), outTsFn)
        else:
            if hasAlignment:
                xfFile = self.getExtraOutFile(ts.getTsId(), ext=XF_EXT)
                try:
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
                    outTsFn, _, _ = self.getTmpFileNames(ts)
                    ts.applyTransform(outTsFn)

                # After that, for the following programs, a new xfFile without the
                # excluded views must be generated to be used by the protocol main program
                genXfFile(ts, xfFile, presentAcqOrders=presentAcqOrders)
            else:
                # Re-stack
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

    def addToOutFailedSet(self,
                          tsId: str,
                          inputsAreTs: bool = True) -> None:
        """ Just copy input item to the failed output set. """
        logger.info(cyanStr(f'Failed TS ---> {tsId}'))
        try:
            with self._lock:
                inputSet = self.getInputTsSet(pointer=True) if inputsAreTs else self.getInputTomoSet(pointer=True)
                output = self.getOutputFailedSet(inputSet, inputsAreTs=inputsAreTs)
                item = self.getCurrentTs(tsId) if inputsAreTs else self.getCurrentTomo(tsId)
                newItem = item.clone()
                newItem.copyInfo(item)
                output.append(newItem)

                if isinstance(item, TiltSeries):
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

    # --------------------------- UTILS functions -----------------------------
    def getInputTsSet(self, pointer: bool = False) -> Union[Pointer, SetOfTiltSeries]:
        tsSetPointer = getattr(self, IN_TS_SET)
        return tsSetPointer if pointer else tsSetPointer.get()

    def getCurrentTs(self, tsId: str) -> TiltSeries:
        return self.getInputTsSet().getItem(TiltSeries.TS_ID_FIELD, tsId)

    def getInputTomoSet(self, pointer: bool = False) -> Union[Pointer, SetOfTomograms]:
        tomoSetPointer = getattr(self, IN_TOMO_SET)
        return tomoSetPointer if pointer else tomoSetPointer.get()

    def getCurrentTomo(self, tsId: str) -> Tomogram:
        return self.getInputTomoSet().getItem(Tomogram.TS_ID_FIELD, tsId)

    def getInputCtfSet(self, pointer: bool = False) -> Union[Pointer, SetOfCTFTomoSeries]:
        tomoSetPointer = getattr(self, IN_CTF_TOMO_SET)
        return tomoSetPointer if pointer else tomoSetPointer.get()

    def getCurrentCtf(self, tsId: str) -> CTFTomoSeries:
        return self.getInputCtfSet().getItem(CTFTomoSeries.TS_ID_FIELD, tsId)

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
            self.info(cyanStr(f'Item processed {self.tsIdReadList}'))
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
                dimX, dimY, _ = ts.getFirstItem().getDim()
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
        outTsFn, outTsOddFn, outTsEvenFn = self.getTmpFileNames(ts)
        with self._lock:
            firstTi = ts.getFirstItem()
        doSwap = self.getNewstackDoSwap(firstTi, xfFile)
        if presentAcqOrders:
            tsExcludedIndices = ts.getTsExcludedViewsIndices(presentAcqOrders)

        param = self.getBasicNewstackParams(ts,
                                            firstTi.getFileName(),
                                            outTsFn,
                                            xfFile=xfFile,
                                            doSwap=doSwap,
                                            tsExcludedIndices=tsExcludedIndices,
                                            binning=binning)
        self.runProgram(NEWSTACK_PROGRAM, param)

        if self.doOddEven:
            logger.info(cyanStr(f'\t--> running {NEWSTACK_PROGRAM} with the ODD tilt-series'))
            param = self.getBasicNewstackParams(ts,
                                                ts.getOddFileName(),
                                                outTsOddFn,
                                                xfFile=xfFile,
                                                doSwap=doSwap,
                                                tsExcludedIndices=tsExcludedIndices,
                                                binning=binning)
            self.runProgram(NEWSTACK_PROGRAM, param)

            logger.info(cyanStr(f'\t--> running {NEWSTACK_PROGRAM} with the EVEN tilt-series'))
            param = self.getBasicNewstackParams(ts,
                                                ts.getEvenFileName(),
                                                outTsEvenFn,
                                                xfFile=xfFile,
                                                doSwap=doSwap,
                                                tsExcludedIndices=tsExcludedIndices,
                                                binning=binning)
            self.runProgram(NEWSTACK_PROGRAM, param)

    def getTmpFileNames(self, ts: TiltSeries) -> Tuple:
        tsId = ts.getTsId()
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
