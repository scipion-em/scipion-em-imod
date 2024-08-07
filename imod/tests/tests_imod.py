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
from os.path import exists

import numpy as np

from imod.constants import OUTPUT_TILTSERIES_NAME, SCIPION_IMPORT, FIXED_DOSE, OUTPUT_TS_INTERPOLATED_NAME, \
    FIDUCIAL_MODEL, PT_FRACTIONAL_OVERLAP, OUTPUT_FIDUCIAL_GAPS_NAME, PATCH_TRACKING, PT_NUM_PATCHES, \
    OUTPUT_FIDUCIAL_NO_GAPS_NAME, OUTPUT_TOMOGRAMS_NAME
from imod.protocols import ProtImodXraysEraser, ProtImodDoseFilter, ProtImodTsNormalization, \
    ProtImodApplyTransformationMatrix, ProtImodImportTransformationMatrix, ProtImodXcorrPrealignment, \
    ProtImodFiducialModel, ProtImodTomoReconstruction
from imod.protocols.protocol_fiducialAlignment import GROUP_ROTATIONS, GROUP_TILTS, DIST_DISABLED, \
    ROT_SOLUTION_CHOICES, MAG_SOLUTION_CHOICES, TILT_SOLUTION_CHOICES, DISTORTION_SOLUTION_CHOICES, \
    ProtImodFiducialAlignment, ONE_ROTATION, FIXED_MAG, DIST_FULL_SOLUTION, ALL_ROTATIONS, ALL_EXCEPT_MIN, \
    DIST_SKEW_ONLY
from imod.protocols.protocol_tsPreprocess import FLOAT_DENSITIES_CHOICES
from pwem import ALIGN_NONE, ALIGN_2D
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr, cyanStr
from tomo.protocols import ProtImportTs, ProtImportTomograms
from tomo.protocols.protocol_import_tomograms import OUTPUT_NAME
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
        cls.runPrevProtocols()

    @classmethod
    def runPrevProtocols(cls):
        print(cyanStr('--------------------------------- RUNNING PREVIOUS PROTOCOLS ---------------------------------'))
        cls._runPreviousProtocols()
        print(cyanStr('\n-------------------------------- PREVIOUS PROTOCOLS FINISHED --------------------------------'))


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
    def _runImportTomograms(cls):
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomos = cls.newProtocol(ProtImportTomograms,
                                          filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                          filesPattern=DataSetRe4STATuto.tomosPattern.value,
                                          samplingRate=DataSetRe4STATuto.sRateBin4.value)  # Bin 4
        cls.launchProtocol(protImportTomos)
        outTomos = getattr(protImportTomos, OUTPUT_NAME, None)
        return outTomos

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
    def _runTsPreprocess(cls, inTsSet, binning=1, densAdjustMode=2, **kwargs):
        print(magentaStr(f"\n==> Running the TS preprocessing:"
                         f"\n\t- Binning factor = {binning}"
                         f"\n\t- Adjust densities mode = {FLOAT_DENSITIES_CHOICES[densAdjustMode]}"))
        protTsNorm = cls.newProtocol(ProtImodTsNormalization,
                                     inputSetOfTiltSeries=inTsSet,
                                     binning=binning,
                                     floatDensities=densAdjustMode,
                                     **kwargs)
        protTsNorm.setObjLabel(f'Bin_{binning} Mode_{densAdjustMode}')
        cls.launchProtocol(protTsNorm)
        tsPreprocessed = getattr(protTsNorm, OUTPUT_TILTSERIES_NAME, None)
        return tsPreprocessed

    @classmethod
    def _runXcorrAli(cls, inTsSet, genInterp=False, cumulativeCorr=False, interpBinning=1, tiltAxisAngle=None):
        tAxMsg = f'manually introduced of {tiltAxisAngle} deg.' if tiltAxisAngle else 'from Scipion metadata'
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

    @classmethod
    def _genFiducialModel(cls, inTsSet, modelType=FIDUCIAL_MODEL, bothSurfaces=False, trackWithModel=True,
                          sizeOfPatches='680 680', patchLayout=PT_FRACTIONAL_OVERLAP, iterationsSubpixel=1,
                          overlapPatches='0.33 0.33', numberOfPatches='2 2', objLabel=None):
        if modelType == FIDUCIAL_MODEL:
            modelTypeStr = 'Make seed and Track'
            displayMsg = (f'\n\t- Find beads on two surfaces = {bothSurfaces}'
                          f'\n\t- Track with fiducial model as seed = {trackWithModel}')
        else:
            if patchLayout == PT_FRACTIONAL_OVERLAP:
                patchLayoutMsg = 'Fractional overlap'
                patchModeInfoMsg = f'\n\t- Fractional overlap value = {overlapPatches}'
            else:
                patchLayoutMsg = 'Number of patches'
                patchModeInfoMsg = f'\n\t- Number of patches = {numberOfPatches}'
            modelTypeStr = 'Patch Tracking'
            displayMsg = (f'\n\t- Size of the patches (X,Y) = {sizeOfPatches}'
                          f'\n\t- Iterations subpixel = {iterationsSubpixel}'
                          f'\n\t- Patch layout = {patchLayoutMsg}'
                          f'{patchModeInfoMsg}')

        print(magentaStr(f"\n==> Generating the TS fiducial model:"
                         f"\n\t- Model Type = {modelTypeStr},"
                         f"{displayMsg}"))
        argsDict = {
            'inputSetOfTiltSeries': inTsSet,
            'typeOfModel': modelType}
        if modelType == FIDUCIAL_MODEL:
            argsDict['twoSurfaces'] = bothSurfaces,
            argsDict['doTrackWithModel'] = trackWithModel
        else:
            argsDict['sizeOfPatches'] = sizeOfPatches
            argsDict['patchLayout'] = patchLayout
            argsDict['iterationsSubpixel'] = iterationsSubpixel
            if patchLayout == PT_FRACTIONAL_OVERLAP:
                argsDict['overlapPatches'] = overlapPatches
            else:
                argsDict['numberOfPatches'] = numberOfPatches

        protFiduAli = cls.newProtocol(ProtImodFiducialModel, **argsDict)
        if objLabel:
            protFiduAli.setObjLabel(objLabel)
        cls.launchProtocol(protFiduAli)
        fiducialModels = getattr(protFiduAli, OUTPUT_FIDUCIAL_GAPS_NAME, None)
        return fiducialModels

    def _checkFiducialModels(self, inFiducialsSet, expectedSetSize=2, expectedFiduSizeAngs=100,
                             presentTsIds=(TS_03, TS_54)):
        self.assertSetSize(inFiducialsSet, expectedSetSize)
        for fiducialModel in inFiducialsSet:
            self.assertTrue(fiducialModel.getTsId() in presentTsIds)
            self.assertTrue(exists(fiducialModel.getFileName()))
            self.assertTrue(exists(fiducialModel.getModelName()))
            self.assertEqual(expectedFiduSizeAngs, fiducialModel.getSize())
            self.assertGreater(fiducialModel.getCount(), 0)

    @classmethod
    def _runFiducialAli(cls, inFiduModels, bothSurfaces=False, genInterp=False, interpBinFactor=-1,
                        rotationType=ONE_ROTATION, magnifType=FIXED_MAG, tiltAngleType=GROUP_TILTS,
                        distortionType=DIST_DISABLED, objLabel=None):
        msg = (f"\n==> Running the TS alignment:"
               f"\n\t- Beads on two surfaces = {bothSurfaces}"
               f"\n\t- Generate the interpolated TS = {genInterp}")
        if genInterp:
            msg += f"\n\t- Interpolated TS binning factor = {interpBinFactor}"
        msg += (f"\n\t- Rotation solution type = {ROT_SOLUTION_CHOICES[rotationType]}"
                f"\n\t- Magnification solution type = {MAG_SOLUTION_CHOICES[magnifType]}"
                f"\n\t- Tilt angle solution type = {TILT_SOLUTION_CHOICES[tiltAngleType]}"
                f"\n\t- Distortion solution type = {DISTORTION_SOLUTION_CHOICES[distortionType]}")
        print(magentaStr(msg))

        protFiduAli = cls.newProtocol(ProtImodFiducialAlignment,
                                      inputSetOfLandmarkModels=inFiduModels,
                                      twoSurfaces=bothSurfaces,
                                      computeAlignment=genInterp,
                                      binning=interpBinFactor,
                                      rotationSolutionType=rotationType,
                                      magnificationSolutionType=magnifType,
                                      tiltAngleSolutionType=tiltAngleType,
                                      distortionSolutionType=distortionType)
        if objLabel:
            protFiduAli.setObjLabel(objLabel)
        cls.launchProtocol(protFiduAli)
        tsAli = getattr(protFiduAli, OUTPUT_TILTSERIES_NAME, None)
        tsInterp = getattr(protFiduAli, OUTPUT_TS_INTERPOLATED_NAME, None)
        fiducialModels = getattr(protFiduAli, OUTPUT_FIDUCIAL_NO_GAPS_NAME, None)
        return tsAli, tsInterp, fiducialModels

    @classmethod
    def _runTomoRec(cls, inTsSet, tomoThickness=-1, tomoWidth=0, tomoShiftX=0, tomoShiftZ=0, superSampleFactor=2,
                    angleOffset=0, tiltAxisOffset=0, fakeInteractionsSIRT=0, objLabel=None):
        print(magentaStr(f"\n==> Reconstructing the tomogram:"
                         f"\n\t- Tomogram thickness = {tomoThickness}"
                         f"\n\t- Tomogram width = {tomoWidth}"
                         f"\n\t- Tomogram shift [x, z] = [{tomoShiftX}, {tomoShiftZ}]"
                         f"\n\t- Tilt angles offset (deg) = {angleOffset}"
                         f"\n\t- Tilt axis angle offset (deg) = {tiltAxisOffset}"
                         f"\n\t- Super-sampling factor = {superSampleFactor}"
                         f"\n\t- Iter. SIRT equivalent filter = {fakeInteractionsSIRT}"))
        protTomoRec = cls.newProtocol(ProtImodTomoReconstruction,
                                      inputSetOfTiltSeries=inTsSet,
                                      tomoThickness=tomoThickness,
                                      tomoWidth=tomoWidth,
                                      tomoShiftX=tomoShiftX,
                                      tomoShiftZ=tomoShiftZ,
                                      angleOffset=angleOffset,
                                      tiltAxisOffset=tiltAxisOffset,
                                      superSampleFactor=superSampleFactor,
                                      fakeInteractionsSIRT=fakeInteractionsSIRT)
        if objLabel:
            protTomoRec.setObjLabel(objLabel)
        cls.launchProtocol(protTomoRec)
        tomograms = getattr(protTomoRec, OUTPUT_TOMOGRAMS_NAME, None)
        return tomograms


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
                                               densAdjustMode=0)  # No adjust
        self._checkTiltSeries(tsPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess01(self):
        binningFactor = 1
        tsPreprocessed = self._runTsPreprocess(self.importedTs,
                                               binning=binningFactor,
                                               densAdjustMode=0,  # No adjust
                                               scaleMax=200,
                                               scaleMin=20)
        self._checkTiltSeries(tsPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess02(self):
        binningFactor = 8
        tsPreprocessed = self._runTsPreprocess(self.importedTs,
                                               binning=binningFactor,
                                               densAdjustMode=1)  # range between min and max
        self._checkTiltSeries(tsPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess03(self):
        binningFactor = 2
        tsPreprocessed = self._runTsPreprocess(self.importedTs,
                                               binning=binningFactor,
                                               densAdjustMode=2,  # scaled to common mean and standard deviation
                                               scaleMean=0,
                                               scaleSd=1)
        self._checkTiltSeries(tsPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess04(self):
        binningFactor = 4
        tsPreprocessed = self._runTsPreprocess(self.importedTs,
                                               binning=binningFactor,
                                               densAdjustMode=2,  # scaled to common mean and standard deviation
                                               meanSdToggle=False)
        self._checkTiltSeries(tsPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess05(self):
        binningFactor = 3
        tsPreprocessed = self._runTsPreprocess(self.importedTs,
                                               binning=binningFactor,
                                               densAdjustMode=3,  # shifted to a common mean without scaling
                                               meanSdToggle=False)
        self._checkTiltSeries(tsPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess06(self):
        binningFactor = 6
        tsPreprocessed = self._runTsPreprocess(self.importedTs,
                                               binning=binningFactor,
                                               densAdjustMode=4,  # shifted to mean and rescaled to a min and max
                                               scaleMax=200,
                                               scaleMin=20)
        self._checkTiltSeries(tsPreprocessed, binningFactor=binningFactor)


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


class TestImodXcorrAlignment(TestImodBase):

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

    def testXcorAli04(self):
        interptBinningFactor = 8
        tiltAxisAngle = 89.1
        xCorrTs, xCorrTsInterp = self._runXcorrAli(self.importedTs,
                                                   genInterp=True,
                                                   interpBinning=interptBinningFactor,
                                                   tiltAxisAngle=tiltAxisAngle)
        # Check the TS
        for tsId, acq in self.testAcqObjDict.items():
            # Update the expected acquisition with the tilt axis angle value introduced manually
            acq.setTiltAxisAngle(tiltAxisAngle)
            self.testAcqObjDict[tsId] = acq
        self._checkTiltSeries(xCorrTs)
        # Check the interpolated TS
        self._checkInterpTiltSeries(xCorrTsInterp, binningFactor=interptBinningFactor)


class TestImodGenFiducialModel(TestImodBase):
    binningFactor = 4

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTs = cls._runImportTs()
        cls.tsPreprocessed = cls._runTsPreprocess(cls.importedTs, binning=cls.binningFactor)
        cls.preAliTsSet, _ = cls._runXcorrAli(cls.tsPreprocessed, genInterp=False)

    def testFiducialModel01(self):
        fiducialModels = self._genFiducialModel(self.preAliTsSet, objLabel='testFiducialModel01')
        # Check the fiducial models
        self._checkFiducialModels(fiducialModels)

    def testFiducialModel02(self):
        fiducialModels = self._genFiducialModel(self.preAliTsSet,
                                                objLabel='testFiducialModel02',
                                                bothSurfaces=True)
        # Check the fiducial models
        self._checkFiducialModels(fiducialModels)

    def testFiducialModel03(self):
        fiducialModels = self._genFiducialModel(self.preAliTsSet,
                                                objLabel='testFiducialModel03',
                                                trackWithModel=False)
        # Check the fiducial models
        self._checkFiducialModels(fiducialModels)

    def testFiducialModel04(self):
        fiducialModels = self._genFiducialModel(self.preAliTsSet,
                                                objLabel='testFiducialModel04',
                                                modelType=PATCH_TRACKING)
        # Check the fiducial models
        self._checkFiducialModels(fiducialModels)

    def testFiducialModel05(self):
        fiducialModels = self._genFiducialModel(self.preAliTsSet,
                                                objLabel='testFiducialModel05',
                                                modelType=PATCH_TRACKING,
                                                sizeOfPatches='420 400')
        # Check the fiducial models
        self._checkFiducialModels(fiducialModels)

    def testFiducialModel06(self):
        fiducialModels = self._genFiducialModel(self.preAliTsSet,
                                                objLabel='testFiducialModel06',
                                                modelType=PATCH_TRACKING,
                                                overlapPatches='0.25 0.33',
                                                iterationsSubpixel=2)
        # Check the fiducial models
        self._checkFiducialModels(fiducialModels)

    def testFiducialModel07(self):
        fiducialModels = self._genFiducialModel(self.preAliTsSet,
                                                objLabel='testFiducialModel07',
                                                modelType=PATCH_TRACKING,
                                                patchLayout=PT_NUM_PATCHES)
        # Check the fiducial models
        self._checkFiducialModels(fiducialModels)


class TestImodTsAlignment(TestImodBase):
    binningFactor = 4

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTs = cls._runImportTs()
        cls.preAliTsSet, _ = cls._runXcorrAli(cls.importedTs, genInterp=False)
        cls.fiducialModels = cls._genFiducialModel(cls.preAliTsSet)

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
                             expectedDimensions=self._getExpectedDimsDict(binningFactor, swapXY=True),
                             testAcqObj=self.testInterpAcqObjDict,
                             anglesCount=self.anglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst)

    def testFiducialAli01(self):
        tsAli, tsInterp, fiducialModels = self._runFiducialAli(self.fiducialModels, objLabel='testFiducialAli01')
        # Check the generated TS
        self._checkTiltSeries(tsAli)
        # Check the interpolated TS
        self.assertIsNone(tsInterp)
        # Check the fiducial models
        self._checkFiducialModels(fiducialModels)

    def testFiducialAli02(self):
        tsAli, tsInterp, fiducialModels = self._runFiducialAli(self.fiducialModels,
                                                               objLabel='testFiducialAli02',
                                                               bothSurfaces=True,
                                                               genInterp=True,
                                                               interpBinFactor=self.binningFactor)
        # Check the generated TS
        self._checkTiltSeries(tsAli)
        # Check the interpolated TS
        self._checkInterpTiltSeries(tsInterp, binningFactor=self.binningFactor)
        # Check the fiducial models
        self._checkFiducialModels(fiducialModels)

    def testFiducialAli03(self):
        tsAli, tsInterp, fiducialModels = self._runFiducialAli(self.fiducialModels,
                                                               objLabel='testFiducialAli03',
                                                               rotationType=GROUP_ROTATIONS,
                                                               distortionType=DIST_FULL_SOLUTION)
        # Check the generated TS
        self._checkTiltSeries(tsAli)
        # Check the interpolated TS
        self.assertIsNone(tsInterp)
        # Check the fiducial models
        self._checkFiducialModels(fiducialModels)

    def testFiducialAli04(self):
        tsAli, tsInterp, fiducialModels = self._runFiducialAli(self.fiducialModels,
                                                               objLabel='testFiducialAli04',
                                                               genInterp=True,
                                                               interpBinFactor=self.binningFactor,
                                                               rotationType=ALL_ROTATIONS,
                                                               tiltAngleType=ALL_EXCEPT_MIN,
                                                               distortionType=DIST_SKEW_ONLY)
        # Check the generated TS
        self._checkTiltSeries(tsAli)
        # Check the interpolated TS
        self._checkInterpTiltSeries(tsInterp, binningFactor=self.binningFactor)
        # Check the fiducial models
        self._checkFiducialModels(fiducialModels)


class TestImodTomoReconstruction(TestImodBase):
    binningFactor = 4

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTs = cls._runImportTs()
        cls.tsPreprocessed = cls._runTsPreprocess(cls.importedTs, binning=cls.binningFactor)
        cls.preAliTsSet, _ = cls._runXcorrAli(cls.tsPreprocessed, genInterp=False)
        cls.fiducialModels = cls._genFiducialModel(cls.preAliTsSet)
        cls.tsAli, _, _ = cls._runFiducialAli(cls.fiducialModels)

    def _checkTomos(self, inTomos, expectedTomoDims=None, expectedOriginShifts=None):
        binnedSRate = self.unbinnedSRate * self.binningFactor
        testOriginShifts = expectedOriginShifts
        if expectedOriginShifts:
            testOriginShifts = - np.array(expectedTomoDims) * binnedSRate / 2
            testOriginShifts[0] -= expectedOriginShifts[0] * binnedSRate
            testOriginShifts[2] -= expectedOriginShifts[1] * binnedSRate
        self.checkTomograms(inTomos,
                            expectedSetSize=len(self.testAcqObjDict),
                            expectedSRate=binnedSRate,
                            expectedDimensions=expectedTomoDims,
                            expectedOriginShifts=testOriginShifts,
                            isHeterogeneousSet=False,
                            testAcqObj=self.testAcqObjDict)

    def testTomoRec01(self):
        tomoThk = 300
        tomograms = self._runTomoRec(self.tsAli,
                                     objLabel='testTomoRec01',
                                     tomoThickness=tomoThk)
        # Check the tomograms
        self._checkTomos(tomograms, expectedTomoDims=[960, 928, tomoThk])

    def testTomoRec02(self):
        tomoThk = 250
        tomoShiftXAngst = 120
        tomoShiftZAngst = 30
        tomograms = self._runTomoRec(self.tsAli,
                                     objLabel='testTomoRec02',
                                     tomoThickness=tomoThk,
                                     tomoShiftX=tomoShiftXAngst,
                                     tomoShiftZ=tomoShiftZAngst)
        # Check the tomograms
        self._checkTomos(tomograms,
                         expectedTomoDims=[960, 928, tomoThk],
                         expectedOriginShifts=[tomoShiftXAngst, tomoShiftZAngst])

    def testTomoRec03(self):
        tomoThk = 320
        tomoWidth = 900
        tomograms = self._runTomoRec(self.tsAli,
                                     objLabel='testTomoRec03',
                                     tomoThickness=tomoThk,
                                     tomoWidth=tomoWidth,
                                     superSampleFactor=4,
                                     )
        # Check the tomograms
        self._checkTomos(tomograms, expectedTomoDims=[tomoWidth, 928, tomoThk])

    def testTomoRec04(self):
        tomoThk = 280
        tomograms = self._runTomoRec(self.tsAli,
                                     objLabel='testTomoRec04',
                                     tomoThickness=tomoThk,
                                     angleOffset=2,
                                     tiltAxisOffset=3)
        # Check the tomograms
        self._checkTomos(tomograms, expectedTomoDims=[960, 928, tomoThk])

    def testTomoRec05(self):
        tomoThk = 220
        tomograms = self._runTomoRec(self.tsAli,
                                     objLabel='testTomoRec05',
                                     tomoThickness=tomoThk,
                                     fakeInteractionsSIRT=5)
        # Check the tomograms
        self._checkTomos(tomograms, expectedTomoDims=[960, 928, tomoThk])






