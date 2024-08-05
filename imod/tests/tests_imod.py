# *
# * Authors:     Scipion Team
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *  e-mail address 'scipion-users@lists.sourceforge.net'
# *
# **************************************************************************
import math

import numpy as np

from imod.constants import OUTPUT_TILTSERIES_NAME, SCIPION_IMPORT, FIXED_DOSE, OUTPUT_TS_INTERPOLATED_NAME
from imod.protocols import ProtImodXraysEraser, ProtImodDoseFilter, ProtImodTsNormalization, \
    ProtImodApplyTransformationMatrix, ProtImodImportTransformationMatrix, ProtImodXcorrPrealignment
from pwem import ALIGN_NONE, ALIGN_2D
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from tomo.protocols import ProtImportTs
from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer

# 5 TS with no. tilt-images:
#   - TS_01 = 40
#   - TS_03 = 40
#   - TS_43 = 41
#   - TS_45 = 41
#   - TS_54 = 41
#
# 5 tomograms with a thickness of (px):
#   - TS_01 = 340
#   - TS_03 = 280
#   - TS_43 = 300
#   - TS_45 = 300
#   - TS_54 = 280
TS_01 = 'TS_01'
TS_03 = 'TS_03'
TS_43 = 'TS_43'
TS_45 = 'TS_45'
TS_54 = 'TS_54'
unbinnedTiDims = [3710, 3838]
tsOriginAngst = - DataSetRe4STATuto.unbinnedPixSize.value * np.array(unbinnedTiDims) / 2
tomoDimsThk300 = [928, 928, 300]
tomoDimsThk280 = [928, 928, 280]
tomoDimsThk340 = [928, 928, 340]


class TestImodBase(TestBaseCentralizedLayer):
    importedTs = None
    unbinnedSRate = DataSetRe4STATuto.unbinnedPixSize.value
    expectedTsSetSize = 2
    testAcqObjDict = {
        TS_03: DataSetRe4STATuto.testAcq03.value,
        TS_54: DataSetRe4STATuto.testAcq54.value,
    }
    testInterpAcqObjDict = {
        TS_03: DataSetRe4STATuto.testAcq03Interp.value,
        TS_54: DataSetRe4STATuto.testAcq54Interp.value,
    }

    anglesCountDict = {
        TS_03: 40,
        TS_54: 41,
    }

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
        cls._runPreviousProtocols()

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTs = cls._runImportTs()

    @staticmethod
    def _getExpectedDimsDict(binningFactor=1, swapXY=False):
        dims = []
        for iDim in unbinnedTiDims:
            newDim = math.ceil(iDim / binningFactor)
            # Imod always generates images with even dimensions
            if newDim % 2 != 0:
                newDim += 1
            dims.append(newDim)

        if swapXY:
            dims.reverse()
        expectedDimensions = {
            TS_03: dims + [40],
            TS_54: dims + [41]
        }
        return expectedDimensions

    @classmethod
    def _runImportTs(cls, exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value):
        print(magentaStr("\n==> Importing the tilt series:"))
        protImportTs = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                       filesPattern=DataSetRe4STATuto.tsPattern.value,
                                       exclusionWords=exclusionWords,
                                       anglesFrom=2,  # From tlt file
                                       voltage=DataSetRe4STATuto.voltage.value,
                                       magnification=DataSetRe4STATuto.magnification.value,
                                       sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                       amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value,
                                       samplingRate=cls.unbinnedSRate,
                                       doseInitial=DataSetRe4STATuto.initialDose.value,
                                       dosePerFrame=DataSetRe4STATuto.dosePerTiltImgWithTltFile.value,
                                       tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value)

        cls.launchProtocol(protImportTs)
        tsImported = getattr(protImportTs, 'outputTiltSeries', None)
        return tsImported

    @classmethod
    def _runXRayEraser(cls, inTsSet):
        print(magentaStr("\n==> Running the X-Ray eraser:"))

        protXRayEraser = cls.newProtocol(ProtImodXraysEraser, inputSetOfTiltSeries=inTsSet)
        cls.launchProtocol(protXRayEraser)
        tsXRayErased = getattr(protXRayEraser, OUTPUT_TILTSERIES_NAME, None)
        return tsXRayErased

    @classmethod
    def _runDoseFilter(cls, inTsSet, fixedDose=False, fixedDoseValue=0):
        if fixedDose:
            doseMsg = 'fixed dose'
            doseType = FIXED_DOSE
        else:
            doseMsg = 'Scipion import dose'
            doseType = SCIPION_IMPORT
        print(magentaStr(f"\n==> Running the dose filter with {doseMsg}:"))
        protDoseFilter = cls.newProtocol(ProtImodDoseFilter,
                                         inputSetOfTiltSeries=inTsSet,
                                         inputDoseType=doseType)
        if fixedDose:
            protDoseFilter.fixedImageDose.set(fixedDoseValue)
        protDoseFilter.setObjLabel(doseMsg)
        cls.launchProtocol(protDoseFilter)
        tsDoseFiltered = getattr(protDoseFilter, OUTPUT_TILTSERIES_NAME, None)
        return tsDoseFiltered

    @classmethod
    def _runImportTrMatrix(cls, inTsSet, binningTM=1, binningTS=1):
        print(magentaStr("\n==> Importing the TS' transformation matrices with IMOD:"
                         f"\n\t- Transformation matrix binning = {binningTM}"
                         f"\n\t- TS binning = {binningTS}"))
        protImportTrMatrix = cls.newProtocol(ProtImodImportTransformationMatrix,
                                             filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                             filesPattern=DataSetRe4STATuto.transformPattern.value,
                                             inputSetOfTiltSeries=inTsSet,
                                             binningTM=binningTM,
                                             binningTS=binningTS)
        protImportTrMatrix.setObjLabel(f'trMat_b{binningTM} ts_b{binningTS}')
        cls.launchProtocol(protImportTrMatrix)
        outTsSet = getattr(protImportTrMatrix, OUTPUT_TILTSERIES_NAME, None)
        return outTsSet

    @classmethod
    def _runApplytTrMatrix(cls, inTsSet, binning=1, taperInside=False, linearInterp=False):
        interpMsg = 'linear' if linearInterp else 'cubic'
        print(magentaStr(f"\n==> Applying the TS' transformation matrices with IMOD:"
                         f"\n\t- Binning = {binning}"
                         f"\n\t- Tapper inside = {taperInside}"
                         f"\n\t- Interpolation type = {interpMsg}"))
        protApplyTrMat = cls.newProtocol(ProtImodApplyTransformationMatrix,
                                         inputSetOfTiltSeries=inTsSet,
                                         binning=binning,
                                         taperInside=taperInside,
                                         linear=linearInterp)
        protApplyTrMat.setObjLabel(f'b_{binning} tapIns_{taperInside} interp_{interpMsg}')
        cls.launchProtocol(protApplyTrMat)
        outTsSet = getattr(protApplyTrMat, OUTPUT_TS_INTERPOLATED_NAME, None)
        return outTsSet

    @classmethod
    def _runTsPreprocess(cls, inTsSet, binning=1, applyAli=False, densAdjustMode=2, **kwargs):
        choices = ['No adjust',
                   'range between min and max',
                   'scaled to common mean and standard deviation',
                   'shifted to a common mean without scaling',
                   'shifted to mean and rescaled to a min and max']
        print(magentaStr(f"\n==> Running the TS preprocessing:"
                         f"\n\t- Binning factor = {binning}"
                         f"\n\t- Apply transformation matrix = {applyAli}"
                         f"\n\t- Adjust densities mode = {choices[densAdjustMode]}"))
        protTsNorm = cls.newProtocol(ProtImodTsNormalization,
                                     inputSetOfTiltSeries=inTsSet,
                                     binning=binning,
                                     applyAlignment=applyAli,
                                     floatDensities=densAdjustMode,
                                     **kwargs)
        protTsNorm.setObjLabel(f'Bin_{binning} Ali_{applyAli} Mode_{densAdjustMode}')
        cls.launchProtocol(protTsNorm)
        tsPreprocessed = getattr(protTsNorm, OUTPUT_TILTSERIES_NAME, None)
        return tsPreprocessed

    @classmethod
    def _runXcorrAli(cls, inTsSet, genInterp=False, cumulativeCorr=False, interpBinning=1, tiltAxisAngle=None):
        tAxMsg = 'manually introduced' if tiltAxisAngle else 'from Scipion metadata'
        print(magentaStr(f"\n==> Running the TS xCorr pre-alignment:"
                         f"\n\t- Generate the interpolated TS = {genInterp}"
                         f"\n\t- Interpolated binning = {interpBinning}"
                         f"\n\t- Tilt axis angle {tAxMsg}"))
        protXcorr = cls.newProtocol(ProtImodXcorrPrealignment,
                                    inputSetOfTiltSeries=inTsSet,
                                    computeAlignment=genInterp,
                                    binning=interpBinning,
                                    cumulativeCorr=cumulativeCorr)
        if tiltAxisAngle:
            protXcorr.tiltAxisAngle.set(tiltAxisAngle)

        protXcorr.setObjLabel(f'GenInterp_{genInterp} cumCor_{cumulativeCorr} ib_{interpBinning} tAx_{tiltAxisAngle}')
        cls.launchProtocol(protXcorr)
        tsXcorr = getattr(protXcorr, OUTPUT_TILTSERIES_NAME, None)
        tsXcorrInterp = getattr(protXcorr, OUTPUT_TS_INTERPOLATED_NAME, None)
        return tsXcorr, tsXcorrInterp


class TestImodXRayEraser(TestImodBase):

    def testXRayEraser(self):
        tsXRayErased = self._runXRayEraser(self.importedTs)
        self.checkTiltSeries(tsXRayErased,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate,
                             imported=True,
                             expectedDimensions=self._getExpectedDimsDict(),
                             testAcqObj=self.testAcqObjDict,
                             anglesCount=self.anglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst)


class TestImodDoseFilter(TestImodBase):

    def _checkTiltSeries(self, inTsSet):
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate,
                             imported=True,
                             expectedDimensions=self._getExpectedDimsDict(),
                             testAcqObj=self.testAcqObjDict,
                             anglesCount=self.anglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst)

    def testDoseFilter01(self):
        tsDoseFilterred = self._runDoseFilter(self.importedTs, fixedDose=SCIPION_IMPORT)
        self._checkTiltSeries(tsDoseFilterred)

    def testDoseFilter02(self):
        tsDoseFilterred = self._runDoseFilter(self.importedTs,
                                              fixedDose=FIXED_DOSE,
                                              fixedDoseValue=DataSetRe4STATuto.dosePerTiltImgWithTltFile.value)
        self._checkTiltSeries(tsDoseFilterred)


class TestImodTsPreprocess(TestImodBase):

    def _checkTiltSeries(self, inTsSet, binningFactor=1, imported=True, hasAlignment=False, alignment=ALIGN_NONE):
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             imported=imported,
                             hasAlignment=hasAlignment,
                             alignment=alignment,
                             expectedDimensions=self._getExpectedDimsDict(binningFactor),
                             testAcqObj=self.testAcqObjDict,
                             anglesCount=self.anglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst)

    def testTsPreprocess00(self):
        binningFactor = 4
        tsPreprocessed = self._runTsPreprocess(self.importedTs,
                                               binning=binningFactor,
                                               applyAli=False,
                                               densAdjustMode=0)  # No adjust
        self._checkTiltSeries(tsPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess01(self):
        binningFactor = 1
        tsPreprocessed = self._runTsPreprocess(self.importedTs,
                                               binning=binningFactor,
                                               applyAli=False,
                                               densAdjustMode=0,  # No adjust
                                               scaleMax=200,
                                               scaleMin=20)
        self._checkTiltSeries(tsPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess02(self):
        binningFactor = 8
        tsPreprocessed = self._runTsPreprocess(self.importedTs,
                                               binning=binningFactor,
                                               applyAli=False,
                                               densAdjustMode=1)  # range between min and max
        self._checkTiltSeries(tsPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess03(self):
        binningFactor = 2
        tsPreprocessed = self._runTsPreprocess(self.importedTs,
                                               binning=binningFactor,
                                               applyAli=False,
                                               densAdjustMode=2,  # scaled to common mean and standard deviation
                                               scaleMean=0,
                                               scaleSd=1)
        self._checkTiltSeries(tsPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess04(self):
        binningFactor = 4
        tsPreprocessed = self._runTsPreprocess(self.importedTs,
                                               binning=binningFactor,
                                               applyAli=False,
                                               densAdjustMode=2,  # scaled to common mean and standard deviation
                                               meanSdToggle=False)
        self._checkTiltSeries(tsPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess05(self):
        binningFactor = 3
        tsPreprocessed = self._runTsPreprocess(self.importedTs,
                                               binning=binningFactor,
                                               applyAli=False,
                                               densAdjustMode=3,  # shifted to a common mean without scaling
                                               meanSdToggle=False)
        self._checkTiltSeries(tsPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess06(self):
        binningFactor = 6
        tsPreprocessed = self._runTsPreprocess(self.importedTs,
                                               binning=binningFactor,
                                               applyAli=False,
                                               densAdjustMode=4,  # shifted to mean and rescaled to a min and max
                                               scaleMax=200,
                                               scaleMin=20)
        self._checkTiltSeries(tsPreprocessed, binningFactor=binningFactor)

    # def testTsPreprocess07(self):
    #     binningFactor = 4
    #     tsImportedMat = self._runImportTrMatrix(self.importedTs)
    #     tsPreprocessed = self._runTsPreprocess(tsImportedMat,
    #                                            binning=binningFactor,
    #                                            applyAli=True,
    #                                            densAdjustMode=4)  # shifted to mean and rescaled to a min and max
    #     self._checkTiltSeries(tsPreprocessed,
    #                           binningFactor=binningFactor,
    #                           imported=False,
    #                           hasAlignment=True,
    #                           alignment=ALIGN_2D)


class TestImodImportTrMatrix(TestImodBase):

    def _checkTiltSeries(self, inTsSet, binningFactor=1):
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             hasAlignment=True,
                             alignment=ALIGN_2D,
                             expectedDimensions=self._getExpectedDimsDict(binningFactor),
                             testAcqObj=self.testAcqObjDict,
                             anglesCount=self.anglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst)

    def testImportTrMatrix01(self):
        tsImportedTrMat = self._runImportTrMatrix(self.importedTs)
        self._checkTiltSeries(tsImportedTrMat)

    def testImportTrMatrix02(self):
        binningFactor = 4
        tsPreprocessed = self._runTsPreprocess(self.importedTs, binning=binningFactor)
        tsImportedTrMat = self._runImportTrMatrix(tsPreprocessed, binningTS=binningFactor)
        self._checkTiltSeries(tsImportedTrMat, binningFactor=binningFactor)


class TestImodApplyTrMatrix(TestImodBase):

    def _checkTiltSeries(self, inTsSet, binningFactor=1):
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             isInterpolated=True,
                             expectedDimensions=self._getExpectedDimsDict(binningFactor, swapXY=True),
                             testAcqObj=self.testInterpAcqObjDict,
                             anglesCount=self.anglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst)

    def testApplyTrMatrix01(self):
        binningFactor = 4
        tsImportedTrMat = self._runImportTrMatrix(self.importedTs)
        tsTrMatrixApplied = self._runApplytTrMatrix(tsImportedTrMat, binning=binningFactor)
        self._checkTiltSeries(tsTrMatrixApplied, binningFactor=binningFactor)

    def testApplyTrMatrix02(self):
        binningFactor = 8
        tsImportedTrMat = self._runImportTrMatrix(self.importedTs)
        tsTrMatrixApplied = self._runApplytTrMatrix(tsImportedTrMat,
                                                    binning=binningFactor,
                                                    taperInside=True,
                                                    linearInterp=True)
        self._checkTiltSeries(tsTrMatrixApplied, binningFactor=binningFactor)


class TestXcorrAlignment(TestImodBase):

    def _checkTiltSeries(self, inTsSet, binningFactor=1):
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             expectedDimensions=self._getExpectedDimsDict(binningFactor),
                             hasAlignment=True,
                             alignment=ALIGN_2D,
                             testAcqObj=self.testAcqObjDict,
                             anglesCount=self.anglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst)

    def _checkInterpTiltSeries(self, inTsSet, binningFactor=1):
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             isInterpolated=True,
                             expectedDimensions=self._getExpectedDimsDict(binningFactor),  # No swap, only translations
                             testAcqObj=self.testInterpAcqObjDict,
                             anglesCount=self.anglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst)

    def testXcorAli01(self):
        xCorrTs, xCorrTsInterp = self._runXcorrAli(self.importedTs, genInterp=False)
        # Check the TS
        self._checkTiltSeries(xCorrTs)
        # Check the interpolated TS
        self.assertIsNone(xCorrTsInterp)

    def testXcorAli02(self):
        interptBinningFactor = 4
        xCorrTs, xCorrTsInterp = self._runXcorrAli(self.importedTs,
                                                   genInterp=True,
                                                   cumulativeCorr=True,
                                                   interpBinning=interptBinningFactor)
        # Check the TS
        self._checkTiltSeries(xCorrTs)
        # Check the interpolated TS
        self._checkInterpTiltSeries(xCorrTsInterp, binningFactor=interptBinningFactor)

    def testXcorAli03(self):
        xCorrTs, xCorrTsInterp = self._runXcorrAli(self.importedTs, genInterp=False)
        # Check the TS
        self._checkTiltSeries(xCorrTs)
        # Check the interpolated TS
        self.assertIsNone(xCorrTsInterp)

# genInterp=False, cumulativeCorr=False, interpBinning=1, tiltAxisAngle=None):