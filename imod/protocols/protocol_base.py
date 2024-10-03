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
import logging
from typing import Union

from pyworkflow.object import Set, CsvList, Boolean
from pyworkflow.protocol import params
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.utils import path
from pwem.emlib.image import ImageHandler as ih
from pwem.protocols import EMProtocol

from tomo.protocols.protocol_base import ProtTomoBase
from tomo.objects import (SetOfTiltSeries, SetOfTomograms, SetOfCTFTomoSeries,
                          CTFTomo, SetOfTiltSeriesCoordinates, TiltSeries,
                          TiltImage, CTFTomoSeries)

from imod import Plugin, utils
from imod.constants import *
from tomo.utils import getCommonTsAndCtfElements

logger = logging.getLogger(__name__)
IN_TS_SET = 'inputSetOfTiltSeries'
IN_TOMO_SET = 'inputSetOfTomograms'
IN_CTF_TOMO_SET = 'inputSetOfCtfTomoSeries'
PROCESS_ODD_EVEN = 'processOddEven'
BINNING_FACTOR = 'binning'


class ProtImodBase(EMProtocol, ProtTomoBase):
    """ Base class with methods used in the rest of the imod protocols. """
    _label = None

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsDict = None
        self.tomoDict = None
        self._failedItems = []
        self.oddEvenFlag = False

        # Possible outputs (synchronize these names with the constants)
        self.TiltSeries = None
        self.Tomograms = None

        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
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

    @staticmethod
    def addTrimingParams(form, pxTrimCondition=False, correlationCondition=True,
                         levelType=params.LEVEL_ADVANCED):
        """
        Generally, this form will be integrated in a groupForm,
        the group form argument is form. A set of flags
        control what elements are shown
        """
        form.addParam('pxTrim',
                      params.NumericListParam,
                      condition=pxTrimCondition,
                      label='Pixels to trim (x y without coma separator)',
                      default="40 40",
                      help='Pixels to trim off each side in X and Y.\n'
                           'Some trimming should be used for patch tracking',
                      expertLevel=levelType)

        xtrimming = form.addLine('Pixels to do correlation along X-axis',
                                 expertLevel=levelType,
                                 condition=correlationCondition,
                                 help="Starting and ending X coordinates of "
                                      "a region to correlate, based on the "
                                      "position of the region at zero tilt.")

        xtrimming.addParam('xmin',
                           params.IntParam,
                           label='X axis min (left)',
                           allowsNull=True,
                           expertLevel=levelType)

        xtrimming.addParam('xmax',
                           params.IntParam,
                           label='X axis max (right)',
                           allowsNull=True,
                           expertLevel=levelType)

        ytrimming = form.addLine('Pixels to do correlation along Y-axis',
                                 expertLevel=levelType,
                                 condition=correlationCondition,
                                 help="Starting and ending Y coordinates "
                                      "of a region to correlate, based on "
                                      "the position of the region at zero tilt.")

        ytrimming.addParam('ymin',
                           params.IntParam,
                           label='Y axis min (top)',
                           allowsNull=True,
                           expertLevel=levelType)

        ytrimming.addParam('ymax',
                           params.IntParam,
                           label='Y axis max (botton)',
                           allowsNull=True,
                           expertLevel=levelType)

    @staticmethod
    def filteringParametersForm(form, condition, levelType=params.LEVEL_NORMAL):
        filtering = form.addGroup('Filtering parameters',
                                  condition=condition,
                                  expertLevel=levelType)

        line1 = filtering.addLine('High pass filter',
                                  expertLevel=levelType,
                                  help="Some high pass filtering, using a small value of Sigma1 such "
                                       "as 0.03, may be needed to keep the program from being misled by very "
                                       "large scale features in the images.  If the images are noisy, some low "
                                       "pass filtering with Sigma2 and Radius2 is appropriate (e.g. 0.05 for "
                                       " Sigma2, 0.25 for Radius2).  If the images are binned, these values "
                                       "specify frequencies in the binned image, so a higher cutoff (less filtering) "
                                       "might be appropriate.\n\n"
                                       "*FilterRadius1*: Low spatial frequencies in the cross-correlation "
                                       "will be attenuated by a Gaussian curve that is 1 "
                                       "at this cutoff radius and falls off below this "
                                       "radius with a standard deviation specified by "
                                       "FilterSigma2. Spatial frequency units range from "
                                       "0 to 0.5.\n"
                                       "*Filter sigma 1*: Sigma value to filter low frequencies in the "
                                       "correlations with a curve that is an inverted "
                                       "Gaussian.  This filter is 0 at 0 frequency and "
                                       "decays up to 1 with the given sigma value. "
                                       "However, if a negative value of radius1 is entered, "
                                       "this filter will be zero from 0 to "
                                       "|radius1| then decay up to 1.")

        line1.addParam('filterRadius1',
                       params.FloatParam,
                       label='Filter radius 1',
                       default=0.0,
                       expertLevel=levelType)

        line1.addParam('filterSigma1',
                       params.FloatParam,
                       label='Filter sigma 1',
                       default=0.03,
                       expertLevel=levelType)

        line2 = filtering.addLine('Low pass filter',
                                  expertLevel=levelType,
                                  help="If the images are noisy, some low "
                                       "pass filtering with Sigma2 and Radius2 is appropriate (e.g. 0.05 for "
                                       " Sigma2, 0.25 for Radius2).  If the images are binned, these values "
                                       "specify frequencies in the binned image, so a higher cutoff (less filtering) "
                                       "might be appropriate.\n\n"
                                       "*Filter radius 2*: High spatial frequencies in the cross-correlation "
                                       "will be attenuated by a Gaussian curve that is 1 "
                                       "at this cutoff radius and falls off above this "
                                       "radius with a standard deviation specified by "
                                       "FilterSigma2.\n"
                                       "*Filter sigma 2*: Sigma value for the Gaussian rolloff below and "
                                       "above the cutoff frequencies specified by "
                                       "FilterRadius1 and FilterRadius2")

        line2.addParam('filterRadius2',
                       params.FloatParam,
                       label='Filter radius 2',
                       default=0.25,
                       expertLevel=levelType)

        line2.addParam('filterSigma2',
                       params.FloatParam,
                       label='Filter sigma 2',
                       default=0.05,
                       expertLevel=levelType)

    @classmethod
    def worksInStreaming(cls):
        """ So far none of them work in streaming. """
        return False

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[])
                       for ts in self.getInputSet()}
        self.oddEvenFlag = self.applyToOddEven(self.getInputSet())

    @staticmethod
    def getPresentAcqOrders(ts: Union[TiltSeries, None] = None,
                            ctf: Union[CTFTomoSeries, None] = None,
                            onlyEnabled: bool = True):
        if not ts or not ctf:
            obj = ts if ts else ctf
            if onlyEnabled:
                return {ti.getAcquisitionOrder() for ti in obj if ti.isEnabled()}
            else:
                return {ti.getAcquisitionOrder() for ti in obj}
        else:
            return getCommonTsAndCtfElements(ts, ctf, onlyEnabled=onlyEnabled)

    def convertInputStep(self, tsId, generateAngleFile=True,
                         imodInterpolation=True, doSwap=False,
                         oddEven=False, presentAcqOrders=None):
        """
        :param tsId: Tilt-series identifier
        :param generateAngleFile:  Boolean(True) to generate IMOD angle file
        :param imodInterpolation: Boolean (True) to interpolate the tilt series with
                                  imod in case there is a TM.
                                  Pass None to cancel interpolation.
        :param doSwap: if applying alignment, consider swapping X/Y
        :param oddEven: process odd/even sets
        :param presentAcqOrders: set containing the present acq orders in both
        the given TS and CTFTomoSeries. Used to generate the xf file, the tlt file,
        and the interpolated TS with IMOD's newstack program.
        """
        ts = self.tsDict[tsId]
        self.genTsPaths(tsId)
        self.genAlignmentFiles(ts, generateAngleFile=generateAngleFile,
                               imodInterpolation=imodInterpolation,
                               doSwap=doSwap,
                               oddEven=oddEven,
                               presentAcqOrders=presentAcqOrders)

    def closeOutputSetsStep(self):
        self._closeOutputSet()

    def applyNewStackBasic(self, ts, outputTsFileName, inputTsFileName,
                           xfFile=None, doSwap=None, tsExcludedIndices=None):

        params = self.getBasicNewstackParams(ts,
                                             outputTsFileName,
                                             inputTsFileName=inputTsFileName,
                                             xfFile=xfFile,
                                             firstItem=ts.getFirstItem(),
                                             doSwap=doSwap,
                                             tsExcludedIndices=tsExcludedIndices)

        self.runProgram("newstack", params)

    def genAlignmentFiles(self, ts: TiltSeries,
                          generateAngleFile: bool = True,
                          imodInterpolation: bool = True,
                          doSwap: bool = False,
                          oddEven: bool = False,
                          presentAcqOrders: Union[set, None] = None):
        """
        :param ts: Tilt-series
        :param generateAngleFile:  Boolean(True) to generate IMOD angle file
        :param imodInterpolation: Boolean (True) to interpolate the tilt series with
                                  imod in case there is a TM.
                                  Pass None to cancel interpolation.
        :param doSwap: if applying alignment, consider swapping X/Y
        :param oddEven: process odd/even sets
        :param presentAcqOrders: set containing the present acq orders in both the
        given TS and/or the CTFTomoSeries. Used to generate the xf file, the tlt file,
        and the interpolated TS with IMOD's newstack program.
        """
        def _linkTs():
            logger.info(f"TS [{tsId}] linked")
            path.createLink(firstTi.getFileName(), outputTsFileName)

            if oddEven:
                path.createLink(fnOdd, outputOddTsFileName)
                path.createLink(fnEven, outputEvenTsFileName)

        def _applyNewStackBasic():
            logger.info(f"TS [{tsId}] re-stacked with IMOD")
            self.applyNewStackBasic(ts, outputTsFileName, inTsFileName,
                                    xfFile=xfFile,
                                    doSwap=doSwap,
                                    tsExcludedIndices=tsExcludedIndices)
            if oddEven:
                self.applyNewStackBasic(ts, outputOddTsFileName, fnOdd,
                                        xfFile=xfFile,
                                        doSwap=doSwap,
                                        tsExcludedIndices=tsExcludedIndices)
                self.applyNewStackBasic(ts, outputEvenTsFileName, fnEven,
                                        xfFile=xfFile,
                                        doSwap=doSwap,
                                        tsExcludedIndices=tsExcludedIndices)

        # Initialization
        tsId = ts.getTsId()
        tsExcludedIndices = None
        firstTi = ts.getFirstItem()
        inTsFileName = firstTi.getFileName()
        outputTsFileName = self.getTmpOutFile(tsId)
        fnOdd = None
        fnEven = None
        outputOddTsFileName = None
        outputEvenTsFileName = None
        if oddEven:
            fnOdd = ts.getOddFileName()
            fnEven = ts.getEvenFileName()
            outputOddTsFileName = self.getTmpOutFile(tsId, suffix=ODD)
            outputEvenTsFileName = self.getTmpOutFile(tsId, suffix=EVEN)

        # Interpolation
        if imodInterpolation is None:
            _linkTs()

        elif imodInterpolation:
            xfFile = self.getExtraOutFile(tsId, ext=XF_EXT)
            # Use IMOD newstack interpolation
            if firstTi.hasTransform():
                # Generate transformation matrices file (xf)
                utils.genXfFile(ts, xfFile)

                # Generate the interpolated TS with IMOD's newstack program
                logger.info(f"TS [{tsId}] will be interpolated with IMOD")
                if presentAcqOrders and len(ts) != len(presentAcqOrders):
                    tsExcludedIndices = [ti.getIndex() for ti in ts if not ti.getAcquisitionOrder() in presentAcqOrders]
                # else:
                #     tsExcludedIndices = ts.getExcludedViewsIndex()
                _applyNewStackBasic()

                # If some views were excluded to generate the new stack,
                # a new xfFile containing them should be generated
                if presentAcqOrders and len(ts) != len(presentAcqOrders):
                    utils.genXfFile(ts, xfFile, presentAcqOrders=presentAcqOrders)

            else:
                # The given TS is interpolated
                logger.info(f"TS [{tsId}] is interpolated or not aligned")
                if presentAcqOrders:
                    if len(presentAcqOrders) == len(ts):
                        _linkTs()
                    else:
                        xfFile = None
                        _applyNewStackBasic()
                else:
                    _linkTs()

        # Use Xmipp interpolation via Scipion
        else:
            logger.info(f"TS [{tsId}] interpolated with emlib")
            ts.applyTransform(outputTsFileName)

        logger.info(f"TS [{tsId}] available for processing at {outputTsFileName}")

        # Generate the tlt file
        if generateAngleFile:
            logger.info(f"Generate angle file for the TS [{tsId}]")
            angleFilePath = self.getExtraOutFile(tsId, ext=TLT_EXT)
            ts.generateTltFile(angleFilePath, presentAcqOrders=presentAcqOrders)

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

            self._defineOutputs(**{attrName: outputSet})
            self._defineSourceRelation(inputPtr, outputSet)

        return outputSet

    def getOutputFiducialModel(self, inputPtr,
                               attrName=OUTPUT_FIDUCIAL_NO_GAPS_NAME,
                               suffix="NoGaps"):
        """ Method to generate output of set fiducial models.
                :param inputPtr: input TS set pointer
                :param attrName: output attr name
                :param suffix: output set suffix
        """
        if not inputPtr.isPointer():
            logger.warning("FOR DEVELOPERS: inputSet must be a pointer!")
            inputSet = inputPtr
        else:
            inputSet = inputPtr.get()

        fidModel = getattr(self, attrName, None)
        if fidModel is not None:
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

    def getOutputSetOfTiltSeriesCoordinates(self, inputPtr):
        tsCoords = getattr(self, OUTPUT_TS_COORDINATES_NAME, None)
        if tsCoords is not None:
            tsCoords.enableAppend()
        else:
            tsCoords = SetOfTiltSeriesCoordinates.create(self._getPath(),
                                                         suffix='Fiducials3D')
            tsCoords.setSetOfTiltSeries(inputPtr)
            tsCoords.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_TS_COORDINATES_NAME: tsCoords})
            self._defineSourceRelation(inputPtr, tsCoords)

        return tsCoords

    def getOutputSetOfCoordinates3Ds(self, inputPtr, outputSet):
        coords3D = getattr(self, OUTPUT_COORDINATES_3D_NAME, None)
        if coords3D is not None:
            coords3D.enableAppend()
        else:
            coords3D = self._createSetOfCoordinates3D(volSet=outputSet,
                                                      suffix='Fiducials3D')
            coords3D.setSamplingRate(outputSet.getSamplingRate())
            coords3D.setPrecedents(outputSet)
            coords3D.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_COORDINATES_3D_NAME: coords3D})
            self._defineSourceRelation(inputPtr, coords3D)

        return coords3D

    def getOutputSetOfTomograms(self, inputPtr, binning=1):
        inputSet = inputPtr.get()

        if self.Tomograms:
            getattr(self, OUTPUT_TOMOGRAMS_NAME).enableAppend()
        else:
            outputSetOfTomograms = self._createSetOfTomograms()
            outputSetOfTomograms.copyInfo(inputSet)

            if binning > 1:
                samplingRate = inputSet.getSamplingRate() * binning
                outputSetOfTomograms.setSamplingRate(samplingRate)

            outputSetOfTomograms.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_TOMOGRAMS_NAME: outputSetOfTomograms})
            self._defineSourceRelation(inputPtr, outputSetOfTomograms)

        return self.Tomograms

    def getOutputSetOfCTFTomoSeries(self, inputPtr, outputSetName):
        inputSet = inputPtr.get()

        outputSetOfCTFTomoSeries = getattr(self, outputSetName, None)

        if outputSetOfCTFTomoSeries is not None:
            outputSetOfCTFTomoSeries.enableAppend()
        else:
            outputSetOfCTFTomoSeries = SetOfCTFTomoSeries.create(self._getPath(),
                                                                 template='CTFmodels%s.sqlite')
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(inputPtr)
            outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputSetName: outputSetOfCTFTomoSeries})
            self._defineCtfRelation(inputSet, outputSetOfCTFTomoSeries)

        return outputSetOfCTFTomoSeries

    def getOutputFailedSet(self, inputPtr):
        """ Create output set for failed TS or tomograms. """
        inputSet = inputPtr.get()
        if isinstance(inputSet, SetOfTiltSeries):
            failedTs = getattr(self, OUTPUT_TS_FAILED_NAME, None)

            if failedTs:
                failedTs.enableAppend()
            else:
                logger.info('Create the set of failed TS')
                failedTs = self._createSetOfTiltSeries(suffix='Failed')
                failedTs.copyInfo(inputSet)
                failedTs.setStreamState(Set.STREAM_OPEN)
                self._defineOutputs(**{OUTPUT_TS_FAILED_NAME: failedTs})
                self._defineSourceRelation(inputPtr, failedTs)

            return failedTs

        elif isinstance(inputSet, SetOfTomograms):
            failedTomos = getattr(self, OUTPUT_TOMOS_FAILED_NAME, None)
            if failedTomos:
                failedTomos.enableAppend()
            else:
                logger.info('Create the set of failed tomograms')
                failedTomos = self._createSetOfTomograms(suffix='Failed')
                failedTomos.copyInfo(inputSet)
                failedTomos.setStreamState(Set.STREAM_OPEN)
                self._defineOutputs(**{OUTPUT_TOMOS_FAILED_NAME: failedTomos})
                self._defineSourceRelation(inputPtr, failedTomos)

            return failedTomos

    def createOutputFailedSet(self, item):
        """ Just copy input item to the failed output set. """
        logger.info(f'Failed TS ---> {item.getTsId()}')
        inputSet = self.getInputSet(pointer=True)
        output = self.getOutputFailedSet(inputSet)
        newItem = item.clone()
        newItem.copyInfo(item)
        output.append(newItem)

        if isinstance(item, TiltSeries):
            newItem.copyItems(item)
            newItem.write(properties=False)

        output.update(newItem)
        output.write()
        self._store(output)

    # --------------------------- UTILS functions -----------------------------
    def genTsPaths(self, tsId):
        """Generate the subdirectories corresponding to the
        current tilt-series in tmp and extra"""
        path.makePath(*[self._getExtraPath(tsId), self._getTmpPath(tsId)])

    @staticmethod
    def getOutTsFileName(tsId, suffix=None, ext=MRCS_EXT):
        return f'{tsId}_{suffix}.{ext}' if suffix else f'{tsId}.{ext}'

    def getTmpOutFile(self, tsId, suffix=None, ext=MRCS_EXT):
        return self._getTmpPath(tsId,
                                self.getOutTsFileName(tsId, suffix=suffix, ext=ext))

    def getExtraOutFile(self, tsId, suffix=None, ext=MRCS_EXT):
        return self._getExtraPath(tsId,
                                  self.getOutTsFileName(tsId, suffix=suffix, ext=ext))

    def getInputSet(self, pointer=False):
        return self.inputSetOfTiltSeries.get() if not pointer else self.inputSetOfTiltSeries

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
    def getBasicNewstackParams(ts, outputTsFileName, inputTsFileName=None,
                               xfFile=None, firstItem=None, binning=1,
                               doSwap=False, tsExcludedIndices=None,
                               doTaper=False, doNorm=False):
        """ Returns basic newstack arguments

        :param ts: Title Series object
        :param outputTsFileName: tilt series output file name after newstack
        :param inputTsFileName: Input tilt series file name. Default to firstItem.getFilename()
        :param xfFile: xf file name, if passed, alignment will be generated and used
        :param firstItem: Optional, otherwise it will be taken from ts
        :param binning: Default to 1. to apply to output size
        :param doSwap: Default False.
        :param tsExcludedIndices: List of indices to be excluded in the tilt-series, starting from 1
        :param doTaper: optionally taper the tilt-series
        :param doNorm: optionally normalize the tilt-series
        """

        if firstItem is None:
            firstItem = ts.getFirstItem()

        if inputTsFileName is None:
            inputTsFileName = firstItem.getFileName()

        # Apply interpolation
        params = {
            '-input': inputTsFileName,
            '-output': outputTsFileName,
            '-bin': binning,
            '-antialias': -1,
            '-imagebinned': 1.0
        }
        if doTaper:
            params["-taper"] = "1,1"
        if doNorm:
            params["-FloatDensities"] = 2

        if xfFile is not None:
            params['-xform'] = xfFile

            if doSwap:
                rotationAngle = ts.getAcquisition().getTiltAxisAngle()
                # Check if rotation angle is greater than 45º. If so,
                # swap x and y dimensions to adapt output image sizes to
                # the final sample disposition.
                if 45 < abs(rotationAngle) < 135:
                    dimX, dimY, _ = firstItem.getDim()
                    params["-size"] = f"{round(dimY / binning)}," \
                                      f"{round(dimX / binning)}"

        if tsExcludedIndices:
            params["-exclude"] = ",".join(map(str, tsExcludedIndices))
            # From IMOD's newstack doc: "sections are numbered from 0 unless -fromone is entered"
            params["-fromone"] = ""

        return params

    @staticmethod
    def parseTSDefocusFile(inputTs, defocusFilePath, newCTFTomoSeries):
        """ Parse tilt-series ctf estimation file.
        :param inputTs: input tilt-series
        :param defocusFilePath: input *.defocus file to be parsed
        :param newCTFTomoSeries: output CTFTomoSeries
        """
        defocusFileFlag = utils.getDefocusFileFlag(defocusFilePath)

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
            raise ValueError(
                f"Defocus file flag {defocusFileFlag} is not supported. Only supported formats "
                "correspond to flags 0, 1, 4, 5, and 37.")

        for i, ti in enumerate(inputTs):
            tiObjId = ti.getObjId()
            newCTFTomo = CTFTomo()
            " Plain estimation (any defocus flag)"
            newCTFTomo._defocusUList = CsvList(pType=float)
            newCTFTomo.setDefocusUList(defocusUDict.get(tiObjId, [0.]))

            if ti.isEnabled():
                # newCTFTomo.setAcquisitionOrder(ti.getAcquisitionOrder())
                # newCTFTomo.setIndex(ti.getIndex())

                # if tiObjId not in defocusUDict.keys() and not ti.isEnabled():
                #     raise IndexError("ERROR IN TILT-SERIES %s: NO CTF ESTIMATED FOR VIEW %d, TILT ANGLE %f" % (
                #         inputTs.getTsId(), tiObjId, inputTs[tiObjId].getTiltAngle()))

                if defocusFileFlag == 1:
                    " Astigmatism estimation "
                    newCTFTomo._defocusVList = CsvList(pType=float)
                    newCTFTomo.setDefocusVList(defocusVDict.get(tiObjId, [0.]))

                    newCTFTomo._defocusAngleList = CsvList(pType=float)
                    newCTFTomo.setDefocusAngleList(defocusAngleDict.get(tiObjId, [0.]))

                elif defocusFileFlag == 4:
                    " Phase-shift information "
                    newCTFTomo._phaseShiftList = CsvList(pType=float)
                    newCTFTomo.setPhaseShiftList(phaseShiftDict.get(tiObjId, [0.]))

                elif defocusFileFlag == 5:
                    " Astigmatism and phase shift estimation "
                    newCTFTomo._defocusVList = CsvList(pType=float)
                    newCTFTomo.setDefocusVList(defocusVDict.get(tiObjId, [0.]))

                    newCTFTomo._defocusAngleList = CsvList(pType=float)
                    newCTFTomo.setDefocusAngleList(defocusAngleDict.get(tiObjId, [0.]))

                    newCTFTomo._phaseShiftList = CsvList(pType=float)
                    newCTFTomo.setPhaseShiftList(phaseShiftDict.get(tiObjId, [0.]))

                elif defocusFileFlag == 37:
                    " Astigmatism, phase shift and cut-on frequency estimation "
                    newCTFTomo._defocusVList = CsvList(pType=float)
                    newCTFTomo.setDefocusVList(defocusVDict.get(tiObjId, [0.]))

                    newCTFTomo._defocusAngleList = CsvList(pType=float)
                    newCTFTomo.setDefocusAngleList(defocusAngleDict.get(tiObjId, [0.]))

                    newCTFTomo._phaseShiftList = CsvList(pType=float)
                    newCTFTomo.setPhaseShiftList(phaseShiftDict.get(tiObjId, [0.]))

                    newCTFTomo._cutOnFreqList = CsvList(pType=float)
                    newCTFTomo.setCutOnFreqList(cutOnFreqDict.get(tiObjId, [0.]))

                newCTFTomo.completeInfoFromList()
            else:
                newCTFTomo.setWrongDefocus()
                newCTFTomo.setEnabled(False)

            newCTFTomo.setIndex(i + 1)
            newCTFTomo.setAcquisitionOrder(ti.getAcquisitionOrder())
            newCTFTomoSeries.append(newCTFTomo)

        newCTFTomoSeries.setIMODDefocusFileFlag(defocusFileFlag)
        newCTFTomoSeries.setNumberOfEstimationsInRangeFromDefocusList()

    def copyTsItems(self, outputTsSet, ts, tsId,
                    updateTsCallback=None,
                    updateTiCallback=None,
                    copyDisabledViews=False,
                    copyId=False,
                    copyTM=True,
                    excludedViews=None,
                    isSemiStreamified=True,
                    **kwargs):
        """ Re-implemented function from tomo.objects. Works on a single TS object.
        Params:
            outputSet: output set of tilt series.
            ts: input TiltSeries.
            tsId: can be used by other methods
            updateTsCallback: optional callback after TiltSeries is created
            updateTiCallback: optional callback after TiltImage is created
            copyDisabled: if True, also copy disabled views.
            copyId: copy ObjId.
            copyTM: copy transformation matrix
            excludedViews: list of excluded views, starting from 1
            isSemiStreamified: boolean used to indicate if the protocol is semiStreamified or not. If True, the outputs
            will be generated updated and stored in the execution of the createOutputStep for each batch of steps
            generated and executed.
        """
        tsOut = TiltSeries(tsId=tsId)
        tsOut.copyInfo(ts, copyId=copyId)
        if updateTsCallback:
            updateTsCallback(tsId, ts, tsOut, **kwargs)
        outputTsSet.append(tsOut)

        angleMin = 999
        angleMax = -999
        accumDose = 0
        initialDose = 999
        tiList = []
        for index, ti in enumerate(ts.iterItems()):
            if ti.isEnabled() or (not ti.isEnabled() and copyDisabledViews):
                tiOut = TiltImage(tsId=tsId)
                tiOut.copyInfo(ti, copyId=copyId, copyTM=copyTM, copyStatus=True)
                if updateTiCallback:
                    originIndex = ti.getIndex()
                    updateTiCallback(originIndex, index, tsId, ts, ti, tsOut, tiOut, **kwargs)
                # Update the acquisition of the TS. The accumDose, angle min and angle max for the re-stacked TS, as
                # these values may change if the removed tilt-images are the first or the last, for example.
                tiAngle = ti.getTiltAngle()
                angleMin = min(tiAngle, angleMin)
                angleMax = max(tiAngle, angleMax)
                accumDose = max(ti.getAcquisition().getAccumDose(), accumDose)
                initialDose = min(ti.getAcquisition().getDoseInitial(), initialDose)

                tiList.append(tiOut)

        if excludedViews:
            # Update the acquisition minAngle and maxAngle values of the tilt-series
            acq = tsOut.getAcquisition()
            acq.setAngleMin(angleMin)
            acq.setAngleMax(angleMax)
            acq.setAccumDose(accumDose)
            acq.setDoseInitial(initialDose)
            tsOut.setAcquisition(acq)
            # Update the acquisition minAngle and maxAngle values of each tilt-image acq while preserving their
            # specific accum and initial dose values
            for tiOut in tiList:
                tiAcq = tiOut.getAcquisition()
                tiAcq.setAngleMin(angleMin)
                tiAcq.setAngleMax(angleMax)
                tiOut.setAcquisition(tiAcq)
                tsOut.append(tiOut)
            tsOut.setAnglesCount(len(tsOut))
        else:
            for tiOut in tiList:
                tsOut.append(tiOut)

        outputTsSet.update(tsOut)
        if isSemiStreamified:
            self._store(outputTsSet)

    def updateTi(self, origIndex, index, tsId, ts, ti, tsOut, tiOut, **kwargs):
        outputLocation = self.getExtraOutFile(tsId)
        tiOut.setLocation(index + 1, outputLocation)

        if self.oddEvenFlag:
            locationOdd = index + 1, self.getExtraOutFile(tsId, suffix=ODD)
            locationEven = index + 1, self.getExtraOutFile(tsId, suffix=EVEN)
            tiOut.setOddEven([ih.locationToXmipp(locationOdd),
                              ih.locationToXmipp(locationEven)])
        else:
            tiOut.setOddEven([])

    # --------------------------- INFO functions ------------------------------
    def _warnings(self):
        warnMsgList = []
        self.warningOddEven(self.getInputSet(), warnMsgList)
        return warnMsgList
