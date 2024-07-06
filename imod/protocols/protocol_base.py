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

from pyworkflow.object import Set, CsvList, Pointer
from pyworkflow.protocol import params
from pyworkflow.utils import path
from pwem.emlib.image import ImageHandler as ih
from pwem.protocols import EMProtocol

from tomo.protocols.protocol_base import ProtTomoBase
from tomo.objects import (SetOfTiltSeries, SetOfTomograms, SetOfCTFTomoSeries,
                          CTFTomo, SetOfTiltSeriesCoordinates, TiltSeries)

from imod import Plugin, utils
from imod.constants import *

logger = logging.getLogger(__name__)


class ProtImodBase(EMProtocol, ProtTomoBase):
    """ Base class with methods used in the rest of the imod protocols. """
    _label = None
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}

    def __init__(self, **kwargs):
        self.tsDict = None
        self.tomoDict = None
        self._failedTs = []
        self._failedTomos = []

        # Possible outputs (synchronize these names with the constants)
        self.TiltSeriesCoordinates = None
        self.FiducialModelNoGaps = None
        self.FiducialModelGaps = None
        self.TiltSeries = None
        self.InterpolatedTiltSeries = None
        self.CTFTomoSeries = None
        self.Tomograms = None
        self.Coordinates3D = None

        self.FailedTiltSeries = None
        self.FailedTomograms = None

        EMProtocol.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions -----------------------
    @staticmethod
    def trimingForm(form, pxTrimCondition=False, correlationCondition=True,
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

    # --------------------------- CALCULUS functions ---------------------------
    def convertInputStep(self, tsObjId, generateAngleFile=True,
                         imodInterpolation=True, doSwap=False,
                         oddEven=False, presentAcqOrders=None):
        """
        :param tsObjId: Tilt-series identifier
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
        if type(tsObjId) is str:
            ts = self.tsDict[tsObjId]
        else:
            tsSet = self.inputSetOfTiltSeries,
            ts = tsSet.get()[tsObjId] if isinstance(tsSet, Pointer) else tsSet[tsObjId]

        self.genTsPaths(ts.getTsId())
        self.genAlignmentFiles(ts, generateAngleFile=generateAngleFile,
                               imodInterpolation=imodInterpolation,
                               doSwap=doSwap, oddEven=oddEven,
                               presentAcqOrders=presentAcqOrders)

    def closeOutputSetsStep(self):
        self._closeOutputSet()

    def runProgram(self, program, params, cwd=None):
        """ Shortcut method to run IMOD's command given input params dict. """
        args = ' '.join(['%s %s' % (k, str(v)) for k, v in params.items()])
        Plugin.runImod(self, program, args, cwd)

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

    def genAlignmentFiles(self, ts, generateAngleFile=True,
                          imodInterpolation=True, doSwap=False,
                          oddEven=False, presentAcqOrders=None):
        """
        :param ts: Tilt-series
        :param generateAngleFile:  Boolean(True) to generate IMOD angle file
        :param imodInterpolation: Boolean (True) to interpolate the tilt series with
                                  imod in case there is a TM.
                                  Pass None to cancel interpolation.
        :param doSwap: if applying alignment, consider swapping X/Y
        :param oddEven: process odd/even sets
        :param presentAcqOrders: set containing the present acq orders in both the
        given TS and CTFTomoSeries. Used to generate the xf file, the tlt file,
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
            logger.info(f"TS [{tsId}] linked")
            path.createLink(inTsFileName, outputTsFileName)

        elif imodInterpolation:
            xfFile = self.getExtraOutFile(tsId, ext=XF_EXT)
            # Use IMOD newstack interpolation
            if firstTi.hasTransform():
                # Generate transformation matrices file (xf)
                utils.genXfFile(ts, xfFile)

                # Generate the interpolated TS with IMOD's newstack program
                logger.info(f"TS [{tsId}] will be interpolated with IMOD")
                if presentAcqOrders:
                    tsExcludedIndices = [ti.getIndex() for ti in ts if not ti.getAcquisitionOrder() in presentAcqOrders]
                else:
                    tsExcludedIndices = [ti.getIndex() for ti in ts if not ti.getAcquisitionOrder()]
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
                    params["-size"] = f"{round(firstItem.getYDim() / binning)},"\
                                      f"{round(firstItem.getXDim() / binning)}"

        if tsExcludedIndices:
            params["-exclude"] = ",".join(map(str, tsExcludedIndices))
            # From IMOD's newstack doc: "sections are numbered from 0 unless -fromone is entered"
            params["-fromone"] = ""

        return params

    # --------------------------- OUTPUT functions ----------------------------
    def getOutputSetOfTS(self, inputSet, binning=1) -> SetOfTiltSeries:
        """ Method to generate output classes of set of tilt-series"""

        outputSetOfTiltSeries = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if outputSetOfTiltSeries:
            outputSetOfTiltSeries.enableAppend()
        else:
            outputSetOfTiltSeries = self._createSetOfTiltSeries()

            if isinstance(inputSet, SetOfTiltSeries):
                outputSetOfTiltSeries.copyInfo(inputSet)
                outputSetOfTiltSeries.setDim(inputSet.getDim())

            elif isinstance(inputSet, SetOfTomograms):
                outputSetOfTiltSeries.setAcquisition(inputSet.getAcquisition())
                outputSetOfTiltSeries.setSamplingRate(inputSet.getSamplingRate())
                outputSetOfTiltSeries.setDim(inputSet.getDim())

            if binning > 1:
                samplingRate = inputSet.getSamplingRate()
                samplingRate *= binning
                outputSetOfTiltSeries.setSamplingRate(samplingRate)

            outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_TILTSERIES_NAME: outputSetOfTiltSeries})
            self._defineSourceRelation(inputSet, outputSetOfTiltSeries)

        return outputSetOfTiltSeries

    def getOutputInterpolatedTS(self, inputSet, binning=1):
        """ Method to generate output interpolated classes of set of tilt-series"""

        if self.InterpolatedTiltSeries:
            self.InterpolatedTiltSeries.enableAppend()
        else:
            interpTS = self._createSetOfTiltSeries(suffix='Interpolated')

            if isinstance(inputSet, SetOfTiltSeries):
                interpTS.copyInfo(inputSet)
                interpTS.setDim(inputSet.getDim())

            elif isinstance(inputSet, SetOfTomograms):
                interpTS.setAcquisition(inputSet.getAcquisition())
                interpTS.setDim(inputSet.getDim())

            if binning > 1:
                samplingRate = inputSet.getSamplingRate() * binning
                interpTS.setSamplingRate(samplingRate)

            interpTS.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_TS_INTERPOLATED_NAME: interpTS})
            self._defineSourceRelation(self.inputSetOfTiltSeries, interpTS)

        return self.InterpolatedTiltSeries

    def getOutputFailedSet(self, inputSet):
        """ Create output set for failed TS or tomograms. """
        if isinstance(inputSet, SetOfTiltSeries):
            if self.FailedTiltSeries:
                self.FailedTiltSeries.enableAppend()
            else:
                failedSet = self._createSetOfTiltSeries(suffix='Failed')
                failedSet.copyInfo(inputSet)
                failedSet.setStreamState(Set.STREAM_OPEN)
                self._defineOutputs(**{OUTPUT_TS_FAILED_NAME: failedSet})
                self._defineSourceRelation(inputSet, failedSet)

            return self.FailedTiltSeries

        elif isinstance(inputSet, SetOfTomograms):
            if self.FailedTomograms:
                self.FailedTomograms.enableAppend()
            else:
                failedSet = self._createSetOfTomograms(suffix='Failed')
                failedSet.copyInfo(inputSet)
                failedSet.setStreamState(Set.STREAM_OPEN)
                self._defineOutputs(**{OUTPUT_TOMOS_FAILED_NAME: failedSet})
                self._defineSourceRelation(inputSet, failedSet)

            return self.FailedTomograms

    def getOutputFiducialModelNoGaps(self, tiltSeries=None):
        if self.FiducialModelNoGaps:
            self.FiducialModelNoGaps.enableAppend()
        else:
            if tiltSeries is None:
                tiltSeriesPointer = self.inputSetOfTiltSeries
                tiltSeries = tiltSeriesPointer.get()
            else:
                tiltSeriesPointer = tiltSeries

            fidModelNoGaps = self._createSetOfLandmarkModels(suffix='NoGaps')

            fidModelNoGaps.copyInfo(tiltSeries)
            fidModelNoGaps.setSetOfTiltSeries(tiltSeriesPointer)
            fidModelNoGaps.setHasResidualInfo(True)
            fidModelNoGaps.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_FIDUCIAL_NO_GAPS_NAME: fidModelNoGaps})
            self._defineSourceRelation(tiltSeriesPointer, fidModelNoGaps)

        return self.FiducialModelNoGaps

    def getOutputFiducialModelGaps(self):
        if self.FiducialModelGaps:
            self.FiducialModelGaps.enableAppend()
        else:
            fidModelGaps = self._createSetOfLandmarkModels(suffix='Gaps')

            fidModelGaps.copyInfo(self.getInputSet())
            fidModelGaps.setSetOfTiltSeries(self.inputSetOfTiltSeries)
            fidModelGaps.setHasResidualInfo(False)
            fidModelGaps.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_FIDUCIAL_GAPS_NAME: fidModelGaps})
            self._defineSourceRelation(self.inputSetOfTiltSeries, fidModelGaps)

        return self.FiducialModelGaps

    def getOutputSetOfTiltSeriesCoordinates(self, setOfTiltSeries=None):
        if self.TiltSeriesCoordinates:
            self.TiltSeriesCoordinates.enableAppend()
        else:
            coords3D = SetOfTiltSeriesCoordinates.create(self._getPath(),
                                                         suffix='Fiducials3D')
            coords3D.setSetOfTiltSeries(setOfTiltSeries)
            coords3D.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_TS_COORDINATES_NAME: coords3D})
            self._defineSourceRelation(setOfTiltSeries, coords3D)

        return self.TiltSeriesCoordinates

    def getOutputSetOfCoordinates3Ds(self, inputSet=None, outputSet=None):
        if self.Coordinates3D:
            self.Coordinates3D.enableAppend()
        else:
            coords3D = self._createSetOfCoordinates3D(volSet=outputSet,
                                                      suffix='Fiducials3D')
            coords3D.setSamplingRate(outputSet.getSamplingRate())
            coords3D.setPrecedents(outputSet)
            coords3D.setBoxSize(32)
            coords3D.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_COORDINATES_3D_NAME: coords3D})
            self._defineSourceRelation(inputSet, coords3D)

        return self.Coordinates3D

    def getOutputSetOfTomograms(self, inputSet, binning=1):
        if self.Tomograms:
            getattr(self, OUTPUT_TOMOGRAMS_NAME).enableAppend()
        else:
            outputSetOfTomograms = self._createSetOfTomograms()

            if isinstance(inputSet, SetOfTomograms):
                outputSetOfTomograms.copyInfo(inputSet)

            elif isinstance(inputSet, SetOfTiltSeries):
                outputSetOfTomograms.setAcquisition(inputSet.getAcquisition())
                outputSetOfTomograms.setSamplingRate(inputSet.getSamplingRate())

            if binning > 1:
                samplingRate = inputSet.getSamplingRate()
                samplingRate *= self.binning.get()
                outputSetOfTomograms.setSamplingRate(samplingRate)

            outputSetOfTomograms.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_TOMOGRAMS_NAME: outputSetOfTomograms})
            self._defineSourceRelation(inputSet, outputSetOfTomograms)

        return self.Tomograms

    def getOutputSetOfCTFTomoSeries(self, outputSetName):
        outputSetOfCTFTomoSeries = getattr(self, outputSetName, None)

        if outputSetOfCTFTomoSeries:
            outputSetOfCTFTomoSeries.enableAppend()
        else:
            outputSetOfCTFTomoSeries = SetOfCTFTomoSeries.create(self._getPath(),
                                                                 template='CTFmodels%s.sqlite')
            ts = self.getInputSet(pointer=True)
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(ts)
            outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputSetName: outputSetOfCTFTomoSeries})
            self._defineCtfRelation(outputSetOfCTFTomoSeries, ts.get())

        return outputSetOfCTFTomoSeries

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

        for ti in inputTs:
            tiObjId = ti.getObjId()
            newCTFTomo = CTFTomo()
            newCTFTomo.setAcquisitionOrder(ti.getAcquisitionOrder())
            newCTFTomo.setIndex(ti.getIndex())

            if tiObjId not in defocusUDict.keys() and not ti.isEnabled():
                raise IndexError("ERROR IN TILT-SERIES %s: NO CTF ESTIMATED FOR VIEW %d, TILT ANGLE %f" % (
                    inputTs.getTsId(), tiObjId, inputTs[tiObjId].getTiltAngle()))

            " Plain estimation (any defocus flag)"
            newCTFTomo._defocusUList = CsvList(pType=float)
            newCTFTomo.setDefocusUList(defocusUDict.get(tiObjId, [0.]))

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
            newCTFTomoSeries.append(newCTFTomo)

        newCTFTomoSeries.setIMODDefocusFileFlag(defocusFileFlag)
        newCTFTomoSeries.setNumberOfEstimationsInRangeFromDefocusList()
        newCTFTomoSeries.calculateDefocusUDeviation(defocusUTolerance=20)
        newCTFTomoSeries.calculateDefocusVDeviation(defocusVTolerance=20)

    def createOutputFailedSet(self, item, presentAcqOrders=None):
        """ Just copy input item to the failed output set. """
        inputSet = self.getInputSet()
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

    def getTsFromTsId(self, tsId):
        tsSet = self.getInputSet()
        return tsSet.getItem(TiltSeries.TS_ID_FIELD, tsId)

    def applyToOddEven(self, ts):
        return (hasattr(self, "processOddEven") and
                self.processOddEven and
                ts.hasOddEven())

    @staticmethod
    def _getOutputDim(fn: str):
        x, y, z, _ = ih.getDimensions(fn)
        return x, y, z
