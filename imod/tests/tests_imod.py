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
import copy
import math
from os.path import exists
from typing import Union
import numpy as np

from pwem import ALIGN_NONE, ALIGN_2D
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr, cyanStr
from tomo.objects import TomoAcquisition, SetOfTiltSeries, SetOfCTFTomoSeries
from tomo.protocols import ProtImportTs, ProtImportTomograms, ProtImportTsCTF
from tomo.protocols.protocol_base import ProtTomoImportAcquisition
from tomo.protocols.protocol_import_ctf import ImportChoice
from tomo.protocols.protocol_import_tomograms import OUTPUT_NAME
from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto, TS_03, TS_54, TS_01, TS_43, TS_45
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer

from cistem.protocols import CistemProtTsCtffind
from imod.constants import *
from imod.protocols import *
from imod.protocols.protocol_base_preprocess import FLOAT_DENSITIES_CHOICES
from imod.protocols.protocol_fiducialAlignment import (GROUP_ROTATIONS, GROUP_TILTS, DIST_DISABLED,
                                                       ROT_SOLUTION_CHOICES, MAG_SOLUTION_CHOICES,
                                                       TILT_SOLUTION_CHOICES, DISTORTION_SOLUTION_CHOICES,
                                                       ONE_ROTATION, FIXED_MAG, DIST_FULL_SOLUTION,
                                                       ALL_ROTATIONS, ALL_EXCEPT_MIN, DIST_SKEW_ONLY)

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
unbinnedTiDims = [3710, 3838]
tsOriginAngst = - DataSetRe4STATuto.unbinnedPixSize.value * np.array(unbinnedTiDims) / 2
tomoDimsThk300 = DataSetRe4STATuto.tomoDimsThk300.value
tomoDimsThk280 = DataSetRe4STATuto.tomoDimsThk280.value
tomoDimsThk340 = DataSetRe4STATuto.tomoDimsThk340.value


class TestImodBase(TestBaseCentralizedLayer):
    importedTs = None
    importedTomos = None
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
    # Excluded views stuff
    excludedViewsDict = {
        TS_03: [0, 38, 39],
        TS_54: [0, 1, 38, 39, 40]
    }
    anglesCountDictExcluded = {
        TS_03: 37,
        TS_54: 36,
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
        print(
            cyanStr('\n-------------------------------- PREVIOUS PROTOCOLS FINISHED ---------------------------------'))

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTs = cls._runImportTs()

    @classmethod
    def _gentestAcqObjDictReStacked(cls, isInterp=False):
        # The accumDose, angle min and angle max for the re-stacked TS, as these values may change if the
        # removed tilt-images are the first or the last, for example.
        testAcqObjDictReStacked = {}
        acq_TS_03 = cls.testAcqObjDict[TS_03].clone()
        acq_TS_03.setAccumDose(111)
        acq_TS_03.setAngleMin(-54)
        acq_TS_03.setAngleMax(54)
        if isInterp:
            acq_TS_03.setTiltAxisAngle(0)
        testAcqObjDictReStacked[TS_03] = acq_TS_03

        acq_TS_54 = cls.testAcqObjDict[TS_54].clone()
        acq_TS_54.setAccumDose(108)
        acq_TS_54.setAngleMin(-54)
        acq_TS_54.setAngleMax(51)
        if isInterp:
            acq_TS_54.setTiltAxisAngle(0)
        testAcqObjDictReStacked[TS_54] = acq_TS_54
        return testAcqObjDictReStacked

    @classmethod
    def _getExpectedDimsDict(cls, unbinnedXYDims=None, nImgsDict=None, binningFactor=1, swapXY=False):
        if not nImgsDict:
            nImgsDict = cls.anglesCountDict
        if not unbinnedXYDims:
            unbinnedXYDims = unbinnedTiDims
        dims = []
        for iDim in unbinnedXYDims:
            newDim = math.ceil(iDim / binningFactor)
            # Imod always generates images with even dimensions, at least for the TS
            if newDim % 2 != 0:
                newDim += 1
            dims.append(newDim)

        if swapXY:
            dims.reverse()
        expectedDimensions = {tsId: dims + [nImgs] for tsId, nImgs in nImgsDict.items()}
        return expectedDimensions

    @classmethod
    def _runImportTs(cls, filesPattern=DataSetRe4STATuto.tsPattern.value,
                     exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value):
        print(magentaStr("\n==> Importing the tilt series:"))
        protImportTs = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                       filesPattern=filesPattern,
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
    def _runImportCtf(cls, isTsSet):
        print(magentaStr("\n==> Importing the CTFs:"))
        protImportCtf = cls.newProtocol(ProtImportTsCTF,
                                        filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                        filesPattern=DataSetRe4STATuto.ctfPattern.value,
                                        importFrom=ImportChoice.CTFFIND.value,
                                        inputSetOfTiltSeries=isTsSet)
        cls.launchProtocol(protImportCtf)
        importedCtfs = getattr(protImportCtf, protImportCtf._possibleOutputs.CTFs.name, None)
        return importedCtfs

    @classmethod
    def _runImportTomograms(cls, filesPattern=DataSetRe4STATuto.tomosPattern.value):
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomos = cls.newProtocol(ProtImportTomograms,
                                          filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                          filesPattern=filesPattern,
                                          samplingRate=DataSetRe4STATuto.sRateBin4.value,  # Bin 4
                                          importAcquisitionFrom=ProtTomoImportAcquisition.MANUAL_IMPORT,
                                          oltage=DataSetRe4STATuto.voltage.value,
                                          sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                          amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value)
        cls.launchProtocol(protImportTomos)
        outTomos = getattr(protImportTomos, OUTPUT_NAME, None)
        return outTomos

    @classmethod
    def _runXRayEraser(cls, inTsSet, excludedViews=None):
        excludedViewsMsg = 'eV' if excludedViews else ''
        print(magentaStr(f"\n==> Running the X-Ray eraser:"
                         f"\n\t- Excluded views = {excludedViews is not None}"))

        protXRayEraser = cls.newProtocol(ProtImodXraysEraser, inputSetOfTiltSeries=inTsSet)
        protXRayEraser.setObjLabel(f'X-Ray eraser {excludedViewsMsg}')
        cls.launchProtocol(protXRayEraser)
        tsXRayErased = getattr(protXRayEraser, OUTPUT_TILTSERIES_NAME, None)
        return tsXRayErased

    @classmethod
    def _runDoseFilter(cls, inTsSet, fixedDose=False, fixedDoseValue=0, excludedViews=None):
        if fixedDose:
            doseMsg = 'fixed dose'
            doseType = FIXED_DOSE
        else:
            doseMsg = 'Scipion import dose'
            doseType = SCIPION_IMPORT
        excludedViewsMsg = 'eV' if excludedViews else None
        print(magentaStr(f"\n==> Running the dose filter with {doseMsg}:"
                         f"\n\tExcluded views = {excludedViews is not None}"))
        protDoseFilter = cls.newProtocol(ProtImodDoseFilter,
                                         inputSetOfTiltSeries=inTsSet,
                                         inputDoseType=doseType)
        if fixedDose:
            protDoseFilter.fixedImageDose.set(fixedDoseValue)
        protDoseFilter.setObjLabel(f'{doseMsg}, {excludedViewsMsg}')
        cls.launchProtocol(protDoseFilter)
        tsDoseFiltered = getattr(protDoseFilter, OUTPUT_TILTSERIES_NAME, None)
        return tsDoseFiltered

    @classmethod
    def _runImportTrMatrix(cls, inTsSet, binningTM=1, binningTS=1,
                           filesPattern=DataSetRe4STATuto.transformPattern.value,
                           exclusionWords=None):
        print(magentaStr("\n==> Importing the TS' transformation matrices with IMOD:"
                         f"\n\t- Files pattern = {filesPattern}"
                         f"\n\t- Excluded words = {exclusionWords}"
                         f"\n\t- Transformation matrix binning = {binningTM}"
                         f"\n\t- TS binning = {binningTS}"))
        protImportTrMatrix = cls.newProtocol(ProtImodImportTransformationMatrix,
                                             filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                             filesPattern=filesPattern,
                                             exclusionWords=exclusionWords,
                                             inputSetOfTiltSeries=inTsSet,
                                             binningTM=binningTM,
                                             binningTS=binningTS)
        protImportTrMatrix.setObjLabel(f'trMat_b{binningTM} ts_b{binningTS}')
        cls.launchProtocol(protImportTrMatrix)
        outTsSet = getattr(protImportTrMatrix, OUTPUT_TILTSERIES_NAME, None)
        return outTsSet

    @classmethod
    def _runApplytTrMatrix(cls, inTsSet, binning=1, taperInside=False, linearInterp=False, excludeViews=False):
        interpMsg = 'linear' if linearInterp else 'cubic'
        excludeViewsMsg = ''
        print(magentaStr(f"\n==> Applying the TS' transformation matrices with IMOD:"
                         f"\n\t- Binning = {binning}"
                         f"\n\t- Tapper inside = {taperInside}"
                         f"\n\t- Interpolation type = {interpMsg}"))
        protApplyTrMat = cls.newProtocol(ProtImodApplyTransformationMatrix,
                                         inputSetOfTiltSeries=inTsSet,
                                         binning=binning,
                                         taperInside=taperInside,
                                         linear=linearInterp)
        if excludeViews:
            # Some views were excluded in the input at metadata level
            cls._excludeSetViews(inTsSet)
            excludeViewsMsg = 'eV'
        protApplyTrMat.setObjLabel(f'b_{binning} tapIns_{taperInside} interp_{interpMsg} {excludeViewsMsg}')
        cls.launchProtocol(protApplyTrMat)
        outTsSet = getattr(protApplyTrMat, OUTPUT_TS_INTERPOLATED_NAME, None)
        return outTsSet

    @classmethod
    def _runTsPreprocess(cls, inTsSet, binning=1, densAdjustMode=2, excludedViews=False, **kwargs):
        print(magentaStr(f"\n==> Running the TS preprocessing:"
                         f"\n\t- Binning factor = {binning}"
                         f"\n\t- Excluded views = {excludedViews}"
                         f"\n\t- Adjust densities mode = {FLOAT_DENSITIES_CHOICES[densAdjustMode]}"))
        excludeViewsMsg = 'eV' if excludedViews else ''
        protTsNorm = cls.newProtocol(ProtImodTsNormalization,
                                     inputSetOfTiltSeries=inTsSet,
                                     binning=binning,
                                     floatDensities=densAdjustMode,
                                     **kwargs)
        # if excludedViews:
            # # # Exclude some views at metadata level
            # # cls._excludeSetViews(inTsSet)
            # excludeViewsMsg = 'eV'
        protTsNorm.setObjLabel(f'Bin_{binning} Mode_{densAdjustMode} {excludeViewsMsg}')
        cls.launchProtocol(protTsNorm)
        tsPreprocessed = getattr(protTsNorm, OUTPUT_TILTSERIES_NAME, None)
        return tsPreprocessed

    @classmethod
    def _runXcorrAli(cls, inTsSet, genInterp=False, cumulativeCorr=False, interpBinning=1, tiltAxisAngle=None,
                     excludeViews=False):
        tAxMsg = f'manually introduced of {tiltAxisAngle} deg.' if tiltAxisAngle else 'from Scipion metadata'
        excludeViewsMsg = ''
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
        if excludeViews:
            # Some views were excluded in the input at metadata level
            cls._excludeSetViews(inTsSet)
            excludeViewsMsg = 'eV'
        protXcorr.setObjLabel(f'GenInterp_{genInterp} cumCor_{cumulativeCorr} ib_{interpBinning} tAx_{tiltAxisAngle} '
                              f'{excludeViewsMsg}')
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
                    angleOffset=0, tiltAxisOffset=0, fakeInteractionsSIRT=0, excludedViews=False, objLabel=None):
        print(magentaStr(f"\n==> Reconstructing the tomogram:"
                         f"\n\t- Tomogram thickness = {tomoThickness}"
                         f"\n\t- Tomogram width = {tomoWidth}"
                         f"\n\t- Tomogram shift [x, z] = [{tomoShiftX}, {tomoShiftZ}]"
                         f"\n\t- Excluded views = {excludedViews}"
                         f"\n\t- Tilt angles offset (deg) = {angleOffset}"
                         f"\n\t- Tilt axis angle offset (deg) = {tiltAxisOffset}"
                         f"\n\t- Super-sampling factor = {superSampleFactor}"
                         f"\n\t- Iter. SIRT equivalent filter = {fakeInteractionsSIRT}"))
        if excludedViews:
            # Exclude some views at metadata level
            cls._excludeSetViews(inTsSet)
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

    @classmethod
    def _runTomogramsPreprocess(cls, inTomoSet, binning=1, densAdjustMode=2, **kwargs):
        print(magentaStr(f"\n==> Running the tomograms preprocessing:"
                         f"\n\t- Binning factor = {binning}"
                         f"\n\t- Adjust densities mode = {FLOAT_DENSITIES_CHOICES[densAdjustMode]}"))
        protTomoPreprocess = cls.newProtocol(ProtImodTomoNormalization,
                                             inputSetOfTomograms=inTomoSet,
                                             binning=binning,
                                             floatDensities=densAdjustMode,
                                             **kwargs)
        protTomoPreprocess.setObjLabel(f'Bin_{binning} Mode_{densAdjustMode}')
        cls.launchProtocol(protTomoPreprocess)
        tomoPreprocessed = getattr(protTomoPreprocess, OUTPUT_TOMOGRAMS_NAME, None)
        return tomoPreprocessed

    @classmethod
    def _runTomoProjection(cls, inTomoSet, minAngle=-60, maxAngle=60, angleStep=10,
                           rotationAxis=ProtImodTomoProjection.AXIS_Y, objLabel=None):
        rotAxisChoices = ['X', 'Y', 'Z']
        print(magentaStr(f"\n==> Running the tomograms projection:"
                         f"\n\t- [AngleMin, AngleMax, Step] deg = [{minAngle}, {maxAngle}, {angleStep}]"
                         f"\n\t- Rotation axis = {rotAxisChoices[rotationAxis]}"))
        protTomoProj = cls.newProtocol(ProtImodTomoProjection,
                                       inputSetOfTomograms=inTomoSet,
                                       minAngle=minAngle,
                                       maxAngle=maxAngle,
                                       stepAngle=angleStep,
                                       rotationAxis=rotationAxis)
        if objLabel:
            protTomoProj.setObjLabel(objLabel)
        cls.launchProtocol(protTomoProj)
        tsProjected = getattr(protTomoProj, OUTPUT_TILTSERIES_NAME, None)
        return tsProjected

    @classmethod
    def _excludeSetViews(cls, inSet: Union[SetOfTiltSeries, SetOfCTFTomoSeries],
                         excludedViewsDict: Union[dict, None] = None) -> None:
        if not excludedViewsDict:
            excludedViewsDict = cls.excludedViewsDict
        objList = [obj.clone(ignoreAttrs=[]) for obj in inSet]
        for obj in objList:
            cls._excIntermediateSetViews(inSet, obj, excludedViewsDict[obj.getTsId()])

    @staticmethod
    def _excIntermediateSetViews(inSet, obj, excludedViewsList):
        tiList = [ti.clone() for ti in obj]
        for i, ti in enumerate(tiList):
            if i in excludedViewsList:
                ti._objEnabled = False
                obj.update(ti)
        obj.write()
        inSet.update(obj)
        inSet.write()

    @classmethod
    def _runExcludeViewsProt(cls, inTsSet, objLabel=None):
        print(magentaStr("\n==> Running the TS exclusion of views:"))
        protExcViews = cls.newProtocol(ProtImodExcludeViews, inputSetOfTiltSeries=inTsSet)
        if objLabel:
            protExcViews.setObjLabel(objLabel)
        cls.launchProtocol(protExcViews)
        outTsSet = getattr(protExcViews, OUTPUT_TILTSERIES_NAME, None)
        return outTsSet

    @classmethod
    def _runEstimateCtf(cls, inTsSet, expectedDefocusValue=1500, defocusTol=200, angleStep=0, angleRange=16,
                        objLabel=None):
        print(magentaStr(f"\n==> Running the CTF estimation:"
                         f"\n\t- Expected defocus value (nm) = {expectedDefocusValue}"
                         f"\n\t- Defocus tol (nm) = {defocusTol}"
                         f"\n\t- Angle step (deg) = {angleStep}"
                         f"\n\t- Angle range (deg) = {angleRange}"))
        protEstimateCtf = cls.newProtocol(ProtImodAutomaticCtfEstimation,
                                          inputSetOfTiltSeries=inTsSet,
                                          expectedDefocusValue=expectedDefocusValue,
                                          defocusTol=defocusTol,
                                          angleStep=angleStep,
                                          angleRange=angleRange)
        if objLabel:
            protEstimateCtf.setObjLabel(objLabel)
        cls.launchProtocol(protEstimateCtf)
        outTsSet = getattr(protEstimateCtf, OUTPUT_CTF_SERIE, None)
        return outTsSet

    @classmethod
    def _runCtfCorrection(cls, inTsSet, inCtfSet, tsSetMsg, ctfSetMsg, defocusTol=200, interpWidth=15):
        print(magentaStr(f"\n==> Running the CTF correction:"
                         f"\n\t- Tilt-series = {tsSetMsg}"
                         f"\n\t- CTFs: {ctfSetMsg}"
                         f"\n\t- Defocus tolerance (nm) = {defocusTol}"
                         f"\n\t- Interpolation width (px) = {interpWidth}"))
        protCtfCorr = cls.newProtocol(ProtImodCtfCorrection,
                                      inputSetOfTiltSeries=inTsSet,
                                      inputSetOfCtfTomoSeries=inCtfSet,
                                      defocusTol=defocusTol,
                                      interpolationWidth=interpWidth)
        objLabel = f'ts {tsSetMsg}, ctf {ctfSetMsg}'
        protCtfCorr.setObjLabel(objLabel)
        cls.launchProtocol(protCtfCorr)
        outTsSet = getattr(protCtfCorr, OUTPUT_TILTSERIES_NAME, None)
        return outTsSet

    @classmethod
    def _runGoldBeadPicker(cls, inTomoSet, beadDiameter=10, objLabel=None):
        print(magentaStr(f"\n==> Running the gold bead picker:"
                         f"\n\t- Bead diameter (nm) = {beadDiameter}"))
        protGbp = cls.newProtocol(ProtImodGoldBeadPicker3d,
                                  inputSetOfTomograms=inTomoSet,
                                  beadDiameter=beadDiameter)
        if objLabel:
            protGbp.setObjLabel(objLabel)
        cls.launchProtocol(protGbp)
        outFiduCoords = getattr(protGbp, OUTPUT_COORDINATES_3D_NAME, None)
        return outFiduCoords


class TestImodXRayEraser(TestImodBase):

    def _checkTiltSeries(self, inTsSet, excludedViewsDict=None):
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate,
                             imported=True,
                             expectedDimensions=self._getExpectedDimsDict(),
                             testAcqObj=self.testAcqObjDict,
                             anglesCount=self.anglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst,
                             excludedViewsDict=excludedViewsDict)

    def testXRayEraser01(self):
        tsXRayErased = self._runXRayEraser(self.importedTs)
        self._checkTiltSeries(tsXRayErased)

    def testXRayEraser02(self):
        importedTs = self._runImportTs()
        # Exclude some views at metadata level
        self._excludeSetViews(importedTs)
        # Run the protocol
        tsXRayErased = self._runXRayEraser(importedTs, excludedViews=True)
        # Check the results
        self._checkTiltSeries(tsXRayErased, excludedViewsDict=self.excludedViewsDict)


class TestImodDoseFilter(TestImodBase):

    def _checkTiltSeries(self, inTsSet, excludedViewsDict=None):
        testAcqDict = {}
        for tsId, acq in self.testAcqObjDict.items():
            # Set initial dose and accum dose to 0 as the output ts are dose-weighted
            newAcq = acq.clone()
            newAcq.setAccumDose(0)
            newAcq.setDoseInitial(0)
            testAcqDict[tsId] = newAcq

        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate,
                             imported=True,
                             expectedDimensions=self._getExpectedDimsDict(),
                             testAcqObj=testAcqDict,
                             anglesCount=self.anglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst,
                             excludedViewsDict=excludedViewsDict)

    def testDoseFilter01(self):
        tsDoseFilterred = self._runDoseFilter(self.importedTs, fixedDose=SCIPION_IMPORT)
        self._checkTiltSeries(tsDoseFilterred)

    def testDoseFilter02(self):
        tsDoseFilterred = self._runDoseFilter(self.importedTs,
                                              fixedDose=FIXED_DOSE,
                                              fixedDoseValue=DataSetRe4STATuto.dosePerTiltImgWithTltFile.value)
        self._checkTiltSeries(tsDoseFilterred)

    def testDoseFilter03(self):
        importedTs = self._runImportTs()
        # Exclude some views at metadata level
        self._excludeSetViews(importedTs)
        # Run the protocol
        tsDoseFilterred = self._runDoseFilter(importedTs, fixedDose=SCIPION_IMPORT)
        # Check the results
        self._checkTiltSeries(tsDoseFilterred, excludedViewsDict=self.excludedViewsDict)


class TestImodTsPreprocess(TestImodBase):

    def _checkTiltSeries(self, inTsSet, binningFactor=1, imported=True, hasAlignment=False, alignment=ALIGN_NONE,
                         testAcqObjDict=None, anglesCountDict=None, excludedViewsDict=None):
        if not testAcqObjDict:
            testAcqObjDict = self.testAcqObjDict
        if not anglesCountDict:
            anglesCountDict = self.anglesCountDict
        expectedDimensions = self._getExpectedDimsDict(binningFactor=binningFactor, nImgsDict=anglesCountDict)
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             imported=imported,
                             hasAlignment=hasAlignment,
                             alignment=alignment,
                             expectedDimensions=expectedDimensions,
                             testAcqObj=testAcqObjDict,
                             anglesCount=anglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst,
                             excludedViewsDict=excludedViewsDict)

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

    def testTsPreprocess07(self):
        binningFactor = 4
        importedTs = self._runImportTs()
        # Exclude some views at metadata level
        self._excludeSetViews(importedTs)
        # Run the protocol
        tsPreprocessed = self._runTsPreprocess(importedTs,
                                               binning=binningFactor,
                                               densAdjustMode=0)  # No adjust
        self._checkTiltSeries(tsPreprocessed,
                              binningFactor=binningFactor,
                              excludedViewsDict=self.excludedViewsDict)


class TestImodImportTrMatrix(TestImodBase):

    def _checkTiltSeries(self, inTsSet, binningFactor=1):
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             hasAlignment=True,
                             alignment=ALIGN_2D,
                             expectedDimensions=self._getExpectedDimsDict(binningFactor=binningFactor),
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


class TestImodImportTrMatrixWithPattern(TestImodBase):

    @classmethod
    def _runPreviousProtocols(cls):
        # Import the 5 tilt-series
        cls.importedTs = cls._runImportTs(filesPattern='*/{TS}.mrc', exclusionWords='output')

    def _checkTiltSeries(self, inTsSet, binningFactor=1, testAcqObjDict=None, expectedDimensionsDict=None,
                         anglesCountDict=None):
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=len(testAcqObjDict),
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             hasAlignment=True,
                             alignment=ALIGN_2D,
                             expectedDimensions=expectedDimensionsDict,
                             testAcqObj=testAcqObjDict,
                             anglesCount=anglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst)

    def testImportTrMatrixWpat01(self):
        # exclusionWords = '01 43 45'  # Imported CTF should be TS_03 and TS_54
        tsImportedTrMat = self._runImportTrMatrix(self.importedTs, filesPattern='*/{TS}.xf')
        testAcqObjDict, expectedDimensionsDict, anglesCountDict = DataSetRe4STATuto.genTestTsDicts()
        self._checkTiltSeries(tsImportedTrMat,
                              testAcqObjDict=testAcqObjDict,
                              expectedDimensionsDict=expectedDimensionsDict,
                              anglesCountDict=anglesCountDict)

    def testImportTrMatrixWpat02(self):
        exclusionWords = '01 43 45'  # Imported CTF should be TS_03 and TS_54
        tsImportedTrMat = self._runImportTrMatrix(self.importedTs,
                                                  filesPattern='*/{TS}.xf',
                                                  exclusionWords=exclusionWords)
        tsIdList = (TS_03, TS_54)
        testAcqObjDict, expectedDimensionsDict, anglesCountDict = DataSetRe4STATuto.genTestTsDicts(tsIdList=tsIdList)
        self._checkTiltSeries(tsImportedTrMat,
                              testAcqObjDict=testAcqObjDict,
                              expectedDimensionsDict=expectedDimensionsDict,
                              anglesCountDict=anglesCountDict)

    def testImportTrMatrixWpat03(self):
        exclusionWords = '03 54'  # Imported CTF should be TS_01, TS_43, TS_45
        tsImportedTrMat = self._runImportTrMatrix(self.importedTs,
                                                  filesPattern='*/{TS}.xf',
                                                  exclusionWords=exclusionWords)
        tsIdList = (TS_01, TS_43, TS_45)
        testAcqObjDict, expectedDimensionsDict, anglesCountDict = DataSetRe4STATuto.genTestTsDicts(tsIdList=tsIdList)
        self._checkTiltSeries(tsImportedTrMat,
                              testAcqObjDict=testAcqObjDict,
                              expectedDimensionsDict=expectedDimensionsDict,
                              anglesCountDict=anglesCountDict)


class TestImodApplyTrMatrix(TestImodBase):

    def _checkTiltSeries(self, inTsSet, binningFactor=1, testAcqObjDict=None, testAnglesCountDict=None):
        if not testAcqObjDict:
            testAcqObjDict = self.testAcqObjDict
        if not testAnglesCountDict:
            testAnglesCountDict = self.anglesCountDict
        expectedDimensions = self._getExpectedDimsDict(binningFactor=binningFactor,
                                                       swapXY=True,
                                                       nImgsDict=testAnglesCountDict)
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             isInterpolated=True,
                             expectedDimensions=expectedDimensions,
                             testAcqObj=testAcqObjDict,
                             anglesCount=testAnglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst)

    def testApplyTrMatrix01(self):
        binningFactor = 4
        tsImportedTrMat = self._runImportTrMatrix(self.importedTs)
        # Run the protocol
        tsTrMatrixApplied = self._runApplytTrMatrix(tsImportedTrMat, binning=binningFactor)
        # Check the results
        self._checkTiltSeries(tsTrMatrixApplied,
                              binningFactor=binningFactor,
                              testAcqObjDict=self.testInterpAcqObjDict)

    def testApplyTrMatrix02(self):
        binningFactor = 8
        tsImportedTrMat = self._runImportTrMatrix(self.importedTs)
        # Run the protocol
        tsTrMatrixApplied = self._runApplytTrMatrix(tsImportedTrMat,
                                                    binning=binningFactor,
                                                    taperInside=True,
                                                    linearInterp=True)
        # Check the results
        self._checkTiltSeries(tsTrMatrixApplied,
                              binningFactor=binningFactor,
                              testAcqObjDict=self.testInterpAcqObjDict)

    def testApplyTrMatrix03(self):
        binningFactor = 4
        tsImportedTrMat = self._runImportTrMatrix(self.importedTs)
        # Exclude some views at metadata level
        self._excludeSetViews(tsImportedTrMat)
        # Run the protocol
        tsTrMatrixApplied = self._runApplytTrMatrix(tsImportedTrMat,
                                                    binning=binningFactor,
                                                    excludeViews=True)
        # Check the results
        self._checkTiltSeries(tsTrMatrixApplied,
                              binningFactor=binningFactor,
                              testAcqObjDict=self._gentestAcqObjDictReStacked(isInterp=True),
                              testAnglesCountDict=self.anglesCountDictExcluded)


class TestImodXcorrAlignment(TestImodBase):

    def _checkTiltSeries(self, inTsSet, testAcqObjDict=None, binningFactor=1, excludedViewsDict=None):
        if not testAcqObjDict:
            testAcqObjDict = self.testAcqObjDict
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             expectedDimensions=self._getExpectedDimsDict(binningFactor=binningFactor),
                             hasAlignment=True,
                             alignment=ALIGN_2D,
                             testAcqObj=testAcqObjDict,
                             anglesCount=self.anglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst,
                             excludedViewsDict=excludedViewsDict)

    def _checkInterpTiltSeries(self, inTsSet, testAcqObjDict, testAnglesCountDict, binningFactor=1):
        expectedDimensions = self._getExpectedDimsDict(binningFactor=binningFactor,
                                                       swapXY=False,  # No swap, only translations
                                                       nImgsDict=testAnglesCountDict)
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             isInterpolated=True,
                             expectedDimensions=expectedDimensions,
                             testAcqObj=testAcqObjDict,
                             anglesCount=testAnglesCountDict,
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
        self._checkInterpTiltSeries(xCorrTsInterp,
                                    self.testInterpAcqObjDict,
                                    self.anglesCountDict,
                                    binningFactor=interptBinningFactor)

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
        testAcqDict = copy.deepcopy(self.testAcqObjDict)
        for tsId, acq in testAcqDict.items():
            # Update the expected acquisition with the tilt axis angle value introduced manually
            acq.setTiltAxisAngle(tiltAxisAngle)
            testAcqDict[tsId] = acq
        self._checkTiltSeries(xCorrTs, testAcqObjDict=testAcqDict)
        # Check the interpolated TS
        self._checkInterpTiltSeries(xCorrTsInterp,
                                    self.testInterpAcqObjDict,
                                    self.anglesCountDict,
                                    binningFactor=interptBinningFactor)

    def testXcorAli05(self):
        interptBinningFactor = 4
        importedTs = self._runImportTs()
        # Exclude some views at metadata level
        self._excludeSetViews(importedTs)
        # Run the protocol
        xCorrTs, xCorrTsInterp = self._runXcorrAli(importedTs,
                                                   genInterp=True,
                                                   cumulativeCorr=True,
                                                   interpBinning=interptBinningFactor)
        # Check the TS
        self._checkTiltSeries(xCorrTs, excludedViewsDict=self.excludedViewsDict)  # Excluded at metadata level
        # Check the interpolated TS
        self._checkInterpTiltSeries(xCorrTsInterp,
                                    self._gentestAcqObjDictReStacked(isInterp=True),
                                    self.anglesCountDictExcluded,
                                    binningFactor=interptBinningFactor)


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

    def testFiducialModel08(self):
        preAliTsSet, _ = self._runXcorrAli(self.tsPreprocessed, genInterp=False)
        # Exclude some views at metadata level
        self._excludeSetViews(preAliTsSet)
        # Run the protocol
        fiducialModels = self._genFiducialModel(preAliTsSet, objLabel='testFiducialModel08, eV')
        # Check the fiducial models
        self._checkFiducialModels(fiducialModels)


class TestImodTsAlignment(TestImodBase):
    binningFactor = 4

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTs = cls._runImportTs()
        cls.preAliTsSet, _ = cls._runXcorrAli(cls.importedTs, genInterp=False)
        cls.fiducialModels = cls._genFiducialModel(cls.preAliTsSet)

    def _checkTiltSeries(self, inTsSet, binningFactor=1, excludedViewsDict=None):
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             expectedDimensions=self._getExpectedDimsDict(binningFactor=binningFactor),
                             hasAlignment=True,
                             alignment=ALIGN_2D,
                             testAcqObj=self.testAcqObjDict,
                             anglesCount=self.anglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst,
                             excludedViewsDict=excludedViewsDict)

    def _checkInterpTiltSeries(self, inTsSet, testAcqObjDict=None, testAnglesCountDict=None, binningFactor=1):
        if not testAcqObjDict:
            testAcqObjDict = self.testInterpAcqObjDict
        if not testAnglesCountDict:
            testAnglesCountDict = self.anglesCountDict
        expectedDimensions = self._getExpectedDimsDict(binningFactor=binningFactor,
                                                       swapXY=True,
                                                       nImgsDict=testAnglesCountDict)
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             isInterpolated=True,
                             expectedDimensions=expectedDimensions,
                             testAcqObj=testAcqObjDict,
                             anglesCount=testAnglesCountDict,
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

    def testFiducialAli05(self):
        preAliTsSet, _ = self._runXcorrAli(self.importedTs, genInterp=False)
        # Exclude some views at metadata level and commpute the fiducial models using them
        self._excludeSetViews(preAliTsSet)
        fiducialModels = self._genFiducialModel(preAliTsSet)
        # Run the protocol
        tsAli, tsInterp, fiducialModels = self._runFiducialAli(fiducialModels,
                                                               objLabel='testFiducialAli05, eV',
                                                               bothSurfaces=True,
                                                               genInterp=True,
                                                               interpBinFactor=self.binningFactor)
        # Check the generated TS
        self._checkTiltSeries(tsAli, excludedViewsDict=self.excludedViewsDict)  # Excluded at metadata level
        # Check the interpolated TS
        self._checkInterpTiltSeries(tsInterp,
                                    self._gentestAcqObjDictReStacked(isInterp=True),
                                    self.anglesCountDictExcluded,
                                    binningFactor=self.binningFactor)
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
                                     tomoWidth=1000,  # The protocol will ignore this as it's greater than the tomo Xdim
                                     fakeInteractionsSIRT=5)
        # Check the tomograms
        self._checkTomos(tomograms, expectedTomoDims=[960, 928, tomoThk])

    def testTomoRec06(self):
        tomoThk = 320
        tomoWidth = 900
        tsAli, _, _ = self._runFiducialAli(self.fiducialModels)
        # Exclude some views at metadata level
        self._excludeSetViews(tsAli)
        # Run the protocol
        tomograms = self._runTomoRec(tsAli,
                                     objLabel='testTomoRec06, eV',
                                     tomoThickness=tomoThk,
                                     tomoWidth=tomoWidth,
                                     superSampleFactor=4,
                                     excludedViews=True)
        # Check the tomograms
        self._checkTomos(tomograms, expectedTomoDims=[tomoWidth, 928, tomoThk])


class TestImodTomogramPreprocess(TestImodBase):
    binningFactor = 4

    def _checkTomos(self, inTomos, binningFactor=1):
        # The input tomograms are at binning 4, so it will be the reference binning
        binnedSRate = self.unbinnedSRate * self.binningFactor * binningFactor
        testAcqObjDict = {
            TS_01: DataSetRe4STATuto.testAcq01.value,
            TS_03: DataSetRe4STATuto.testAcq03.value,
            TS_43: DataSetRe4STATuto.testAcq43.value,
            TS_45: DataSetRe4STATuto.testAcq45.value,
            TS_54: DataSetRe4STATuto.testAcq54.value,
        }
        tomoDimsThk280b = [round(iDim / binningFactor) for iDim in tomoDimsThk280]
        tomoDimsThk300b = [round(iDim / binningFactor) for iDim in tomoDimsThk300]
        tomoDimsThk340b = [round(iDim / binningFactor) for iDim in tomoDimsThk340]
        expectedDimensionsDict = {
            TS_01: tomoDimsThk340b,
            TS_03: tomoDimsThk280b,
            TS_43: tomoDimsThk300b,
            TS_45: tomoDimsThk300b,
            TS_54: tomoDimsThk280b,
        }
        self.checkTomograms(inTomos,
                            expectedSetSize=len(testAcqObjDict),
                            expectedSRate=binnedSRate,
                            expectedDimensions=expectedDimensionsDict,
                            isHeterogeneousSet=False,
                            testAcqObj=testAcqObjDict)

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTomos = cls._runImportTomograms()

    def testTsPreprocess00(self):
        tomosPreprocessed = self._runTomogramsPreprocess(self.importedTomos,
                                                         densAdjustMode=0)  # No adjust
        self._checkTomos(tomosPreprocessed)

    def testTsPreprocess01(self):
        binningFactor = 1
        tomosPreprocessed = self._runTomogramsPreprocess(self.importedTomos,
                                                         binning=binningFactor,
                                                         densAdjustMode=0,  # No adjust
                                                         scaleMax=200,
                                                         scaleMin=20)
        self._checkTomos(tomosPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess02(self):
        binningFactor = 2
        tomosPreprocessed = self._runTomogramsPreprocess(self.importedTomos,
                                                         binning=binningFactor,
                                                         densAdjustMode=1)  # range between min and max
        self._checkTomos(tomosPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess03(self):
        binningFactor = 2
        tomosPreprocessed = self._runTomogramsPreprocess(self.importedTomos,
                                                         binning=binningFactor,
                                                         densAdjustMode=2,
                                                         # scaled to common mean and standard deviation
                                                         scaleMean=0,
                                                         scaleSd=1)
        self._checkTomos(tomosPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess04(self):
        tomosPreprocessed = self._runTomogramsPreprocess(self.importedTomos,
                                                         densAdjustMode=2,
                                                         # scaled to common mean and standard deviation
                                                         meanSdToggle=False)
        self._checkTomos(tomosPreprocessed)

    def testTsPreprocess05(self):
        binningFactor = 3
        tomosPreprocessed = self._runTomogramsPreprocess(self.importedTomos,
                                                         binning=binningFactor,
                                                         densAdjustMode=3,  # shifted to a common mean without scaling
                                                         meanSdToggle=False)
        self._checkTomos(tomosPreprocessed, binningFactor=binningFactor)

    def testTsPreprocess06(self):
        tomosPreprocessed = self._runTomogramsPreprocess(self.importedTomos,
                                                         densAdjustMode=4,
                                                         # shifted to mean and rescaled to a min and max
                                                         scaleMax=200,
                                                         scaleMin=20)
        self._checkTomos(tomosPreprocessed)


class TestImodTomoProjection(TestImodBase):

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTomos = cls._runImportTomograms(filesPattern='*3.mrc')  # TS_03 and TS_43

    def _checkTiltSeries(self, inTsSet, testAcqObjDict, anglesCountDict, binningFactor=4):  # Binned 4 tomograms used
        expectedDimensions = self._getExpectedDimsDict(nImgsDict=anglesCountDict,
                                                       unbinnedXYDims=[3712, 3712],  # Square tomograms
                                                       binningFactor=binningFactor)
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             expectedDimensions=expectedDimensions,
                             imported=True,
                             testAcqObj=testAcqObjDict,
                             anglesCount=anglesCountDict,
                             expectedOrigin=tsOriginAngst)

    def checkTomoAcquisition(self, testAcq: Union[TomoAcquisition, dict], currentAcq: TomoAcquisition,
                             tsId: Union[str, None] = None,
                             isTomogramAcq: bool = False) -> None:
        # As the TS' are generated from imported tomograms, some acquisition params won't be present as in the "normal"
        # TS acquisition. Thus, we make this override to custom the required checks for this specific batch of tests

        testAcq = testAcq[tsId] if type(testAcq) is dict else testAcq
        self.assertAlmostEqual(testAcq.getVoltage(), currentAcq.getVoltage(), delta=1)
        self.assertAlmostEqual(testAcq.getSphericalAberration(), currentAcq.getSphericalAberration(), delta=0.01)
        self.assertAlmostEqual(testAcq.getAmplitudeContrast(), currentAcq.getAmplitudeContrast(), delta=0.01)
        self.assertAlmostEqual(testAcq.getTiltAxisAngle(), currentAcq.getTiltAxisAngle(), delta=0.5)
        self.assertAlmostEqual(testAcq.getAngleMin(), currentAcq.getAngleMin(), delta=0.01)
        self.assertAlmostEqual(testAcq.getAngleMax(), currentAcq.getAngleMax(), delta=0.01)
        self.assertAlmostEqual(testAcq.getStep(), currentAcq.getStep(), delta=0.1)

    @staticmethod
    def _genTestData(minAngle, maxAngle, angleStep):
        nAngles = len(range(minAngle, maxAngle + 1, angleStep))
        nImgsDict = {
            TS_03: nAngles,
            TS_43: nAngles
        }
        acqDict = {
            TS_03: DataSetRe4STATuto.testAcq03.value.clone(),
            TS_43: DataSetRe4STATuto.testAcq43.value.clone()
        }
        for tsId, acq in acqDict.items():
            acq.setAngleMin(minAngle)
            acq.setAngleMax(maxAngle)
            acq.setStep(angleStep)
            acq.setMagnification(None)
            acq.setTiltAxisAngle(0)
        return acqDict, nImgsDict

    def testTomoProj01(self):
        minAngle = -60
        maxAngle = 60
        angleStep = 10
        # Launch the protocol
        projTs = self._runTomoProjection(inTomoSet=self.importedTomos,
                                         minAngle=minAngle,
                                         maxAngle=maxAngle,
                                         angleStep=angleStep,
                                         objLabel='testTomoProj01')
        # Check results
        acqDict, nImgsDict = self._genTestData(minAngle, maxAngle, angleStep)
        self._checkTiltSeries(projTs, acqDict, nImgsDict)

    def testTomoProj02(self):
        minAngle = -54
        maxAngle = 51
        angleStep = 6
        # Launch the protocol
        projTs = self._runTomoProjection(inTomoSet=self.importedTomos,
                                         minAngle=minAngle,
                                         maxAngle=maxAngle,
                                         angleStep=angleStep,
                                         rotationAxis=ProtImodTomoProjection.AXIS_X,
                                         objLabel='testTomoProj02')
        # Check results
        acqDict, nImgsDict = self._genTestData(minAngle, maxAngle, angleStep)
        self._checkTiltSeries(projTs, acqDict, nImgsDict)


class TestImodExcludeViews(TestImodBase):

    @classmethod
    def _runPreviousProtocols(cls):
        pass

    def _checkTiltSeries(self, inTsSet, testAcqObjDict, anglesCountDict, binningFactor=1):
        expectedDimensions = self._getExpectedDimsDict(nImgsDict=anglesCountDict, binningFactor=binningFactor)
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             expectedDimensions=expectedDimensions,
                             imported=True,
                             testAcqObj=testAcqObjDict,
                             anglesCount=anglesCountDict,
                             isHeterogeneousSet=True,
                             expectedOrigin=tsOriginAngst)

    def testExcludeViews01(self):
        importedTs = self._runImportTs()
        # There are no excluded views at metadata level, so the protocol should do nothing
        # Run the protocol
        outTsSet = self._runExcludeViewsProt(importedTs, objLabel='testExcludeViews01')
        # Check the results
        self._checkTiltSeries(outTsSet, self.testAcqObjDict, self.anglesCountDict)

    def testExcludeViews02(self):
        importedTs = self._runImportTs()
        # Exclude some views at metadata level
        self._excludeSetViews(importedTs)
        # Run the protocol
        outTsSet = self._runExcludeViewsProt(importedTs, objLabel='testExcludeViews02')
        # Check the results
        self._checkTiltSeries(outTsSet,
                              testAcqObjDict=self._gentestAcqObjDictReStacked(),
                              anglesCountDict=self.anglesCountDictExcluded)

    def testExcludeViews03(self):
        # Other views excluded respecting the previous tests
        excludedViewsDict = {
            TS_03: [0, 1, 38, 39],
            TS_54: [0, 39, 40]
        }
        anglesCountDictExcluded = {
            TS_03: 36,
            TS_54: 38,
        }
        testAcqObjDictReStacked = {}
        acq_TS_03 = self.testAcqObjDict[TS_03].clone()
        acq_TS_03.setAccumDose(111)
        acq_TS_03.setAngleMin(-51)
        acq_TS_03.setAngleMax(54)
        testAcqObjDictReStacked[TS_03] = acq_TS_03

        acq_TS_54 = self.testAcqObjDict[TS_54].clone()
        acq_TS_54.setAccumDose(117)
        acq_TS_54.setAngleMin(-57)
        acq_TS_54.setAngleMax(54)
        testAcqObjDictReStacked[TS_54] = acq_TS_54

        importedTs = self._runImportTs()
        # Exclude some views at metadata level
        self._excludeSetViews(importedTs, excludedViewsDict=excludedViewsDict)
        # Run the protocol
        outTsSet = self._runExcludeViewsProt(importedTs, objLabel='testExcludeViews03')
        # Check the results
        self._checkTiltSeries(outTsSet,
                              testAcqObjDict=testAcqObjDictReStacked,
                              anglesCountDict=anglesCountDictExcluded)


class TestImodEstimateCtf(TestImodBase):

    def _checkCtfs(self, inCtfSet, excludedViewsDict=None):
        expectedSetSize = 2  # TS_03 and TS_54
        self.checkCTFs(inCtfSet,
                       expectedSetSize=expectedSetSize,
                       excludedViewsDict=excludedViewsDict)

    def testEstimateCtf01(self):
        ctfs = self._runEstimateCtf(self.importedTs, objLabel='testEstimateCtf01')
        self._checkCtfs(ctfs)

    def testEstimateCtf02(self):
        ctfs = self._runEstimateCtf(self.importedTs,
                                    angleStep=2,
                                    objLabel='testEstimateCtf02')
        self._checkCtfs(ctfs)

    def testEstimateCtf03(self):
        ctfs = self._runEstimateCtf(self.importedTs,
                                    angleStep=2,
                                    angleRange=10,
                                    objLabel='testEstimateCtf03')
        self._checkCtfs(ctfs)

    def testEstimateCtf04(self):
        importedTs = self._runImportTs()
        self._excludeSetViews(importedTs)  # Excluded some views in the TS at metadata level
        # Run the protocol
        ctfs = self._runEstimateCtf(importedTs, objLabel='testEstimateCtf04, eV')
        # Check the results
        self._checkCtfs(ctfs, excludedViewsDict=self.excludedViewsDict)


class TestImodCtfCorrection(TestImodBase):
    UNMODIFIED = 'unmodified'
    EXC_VIEWS = 'exc. views'
    RE_STACKED = 're-stacked'
    ctfExcludedViewsDict = {
        TS_03: [0, 1, 38, 39],
        TS_54: [0, 39, 40]
    }
    ctfAnglesCountDictExcluded = {
        TS_03: 36,
        TS_54: 38,
    }
    intersectAnglesCountDictExcluded = {
        TS_03: 36,
        TS_54: 36,
    }

    # excludedViewsDict = {
    #     TS_03: [0, 38, 39],
    #     TS_54: [0, 1, 38, 39, 40]
    # }
    # anglesCountDictExcluded = {
    #     TS_03: 37,
    #     TS_54: 36,
    # }

    @classmethod
    def _runPrevProts(cls, importCtf=True):
        importedCtfs = None
        if importCtf:
            importedCtfs = cls._runImportCtf(cls.importedTs)
        tsWithAlignment = cls._runImportTrMatrix(cls.importedTs)
        tsWithAliBin4 = cls._runTsPreprocess(tsWithAlignment, binning=4)
        return importedCtfs, tsWithAliBin4

    @classmethod
    def _gentestAcqObjDictReStacked(cls, ts03MinAngle=-999, ts03MaxAngle=999, ts54MinAngle=-999, ts54MaxAngle=999):
        # The angle min and angle max for the re-stacked TS, as these values may change if the
        # removed tilt-images are the first or the last, for example. FGor the CTF correction, the tilt
        # axis angle and the initial and accum dose will be always 0.
        testAcqObjDictReStacked = {}
        acq_TS_03 = cls.testAcqObjDict[TS_03].clone()
        acq_TS_03.setDoseInitial(0)
        acq_TS_03.setAccumDose(0)
        acq_TS_03.setTiltAxisAngle(0)
        acq_TS_03.setAngleMin(ts03MinAngle)
        acq_TS_03.setAngleMax(ts03MaxAngle)
        testAcqObjDictReStacked[TS_03] = acq_TS_03

        acq_TS_54 = cls.testAcqObjDict[TS_54].clone()
        acq_TS_54.setDoseInitial(0)
        acq_TS_54.setAccumDose(0)
        acq_TS_54.setTiltAxisAngle(0)
        acq_TS_54.setAngleMin(ts54MinAngle)
        acq_TS_54.setAngleMax(ts54MaxAngle)
        testAcqObjDictReStacked[TS_54] = acq_TS_54
        return testAcqObjDictReStacked

    @classmethod
    def _runCistemEstimateCtf(cls, inTsSet):
        print(magentaStr("\n==> Estimating the CTF with Cistem:"))
        protEstimateCtf = cls.newProtocol(CistemProtTsCtffind,
                                          inputTiltSeries=inTsSet,
                                          lowRes=50,
                                          highRes=5,
                                          minDefocus=5000,
                                          maxDefocus=50000)
        cls.launchProtocol(protEstimateCtf)
        ctfs = getattr(protEstimateCtf, CistemProtTsCtffind._possibleOutputs.CTFs.name, None)
        return ctfs

    @classmethod
    def _genReStackedCtf(cls):
        importedTs = cls._runImportTs()
        # Exclude some views from the TS at metadata level
        cls._excludeSetViews(importedTs, excludedViewsDict=cls.ctfExcludedViewsDict)
        # Re-stack that TS
        reStackedTsSet = cls._runExcludeViewsProt(importedTs)
        # Estimate the CTF using the re-stacked TS
        return cls._runCistemEstimateCtf(reStackedTsSet)

    def _checkInterpTiltSeries(self, inTsSet, testAcqObjDict, anglesCountDict, binningFactor=4):
        expectedDimensions = self._getExpectedDimsDict(nImgsDict=anglesCountDict,
                                                       binningFactor=binningFactor,
                                                       swapXY=True)
        self.checkTiltSeries(inTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             isInterpolated=True,
                             expectedDimensions=expectedDimensions,
                             testAcqObj=testAcqObjDict,
                             anglesCount=anglesCountDict,
                             isHeterogeneousSet=True,
                             hasCtfCorrected=True,
                             expectedOrigin=tsOriginAngst)

    def testCtfCorrection01(self):
        importedCtfs, tsWithAliBin4 = self._runPrevProts()
        tsSetCtfCorr = self._runCtfCorrection(tsWithAliBin4, importedCtfs,
                                              tsSetMsg=self.UNMODIFIED,
                                              ctfSetMsg=self.UNMODIFIED)
        self._checkInterpTiltSeries(tsSetCtfCorr,
                                    testAcqObjDict=self.testInterpAcqObjDict,
                                    anglesCountDict=self.anglesCountDict)

    def testCtfCorrection02(self):
        importedCtfs, tsWithAliBin4 = self._runPrevProts()
        self._excludeSetViews(importedCtfs,
                              excludedViewsDict=self.ctfExcludedViewsDict)  # Excluded some views in the CTF at metadata level
        tsSetCtfCorr = self._runCtfCorrection(tsWithAliBin4, importedCtfs,
                                              tsSetMsg=self.UNMODIFIED,
                                              ctfSetMsg=self.EXC_VIEWS)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-51,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-57,
                                                          ts54MaxAngle=54)
        self._checkInterpTiltSeries(tsSetCtfCorr,
                                    testAcqObjDict=testAcqObjDict,
                                    anglesCountDict=self.ctfAnglesCountDictExcluded)

    def testCtfCorrection03(self):
        importedCtfs, tsWithAliBin4 = self._runPrevProts()
        self._excludeSetViews(tsWithAliBin4)  # Excluded some views in the TS at metadata level
        tsSetCtfCorr = self._runCtfCorrection(tsWithAliBin4, importedCtfs,
                                              tsSetMsg=self.EXC_VIEWS,
                                              ctfSetMsg=self.UNMODIFIED)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-54,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-54,
                                                          ts54MaxAngle=51)
        self._checkInterpTiltSeries(tsSetCtfCorr,
                                    testAcqObjDict=testAcqObjDict,
                                    anglesCountDict=self.anglesCountDictExcluded)

    def testCtfCorrection04(self):
        importedCtfs, tsWithAliBin4 = self._runPrevProts()
        self._excludeSetViews(tsWithAliBin4)  # Excluded some views in the TS at metadata level
        tsSetReStacked = self._runExcludeViewsProt(tsWithAliBin4)  # Re-stack the TS
        tsSetCtfCorr = self._runCtfCorrection(tsSetReStacked, importedCtfs,
                                              tsSetMsg=self.RE_STACKED,
                                              ctfSetMsg=self.UNMODIFIED)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-54,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-54,
                                                          ts54MaxAngle=51)
        self._checkInterpTiltSeries(tsSetCtfCorr,
                                    testAcqObjDict=testAcqObjDict,
                                    anglesCountDict=self.anglesCountDictExcluded)

    def testCtfCorrection05(self):
        importedCtfs, tsWithAliBin4 = self._runPrevProts()
        self._excludeSetViews(tsWithAliBin4)  # Excluded some views in the TS at metadata level
        self._excludeSetViews(tsWithAliBin4,
                              excludedViewsDict=self.ctfExcludedViewsDict)  # Excluded some views in the CTF at metadata level
        tsSetCtfCorr = self._runCtfCorrection(tsWithAliBin4, importedCtfs,
                                              tsSetMsg=self.EXC_VIEWS,
                                              ctfSetMsg=self.EXC_VIEWS)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-51,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-54,
                                                          ts54MaxAngle=51)
        self._checkInterpTiltSeries(tsSetCtfCorr,
                                    testAcqObjDict=testAcqObjDict,
                                    anglesCountDict=self.intersectAnglesCountDictExcluded)

    def testCtfCorrection06(self):
        importedCtfs, tsWithAliBin4 = self._runPrevProts()
        self._excludeSetViews(tsWithAliBin4)  # Excluded some views in the TS at metadata level
        tsSetReStacked = self._runExcludeViewsProt(tsWithAliBin4)  # Re-stack the TS
        self._excludeSetViews(importedCtfs,
                              excludedViewsDict=self.ctfExcludedViewsDict)  # Excluded some views in the CTF at metadata level
        tsSetCtfCorr = self._runCtfCorrection(tsSetReStacked, importedCtfs,
                                              tsSetMsg=self.RE_STACKED,
                                              ctfSetMsg=self.EXC_VIEWS)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-51,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-54,
                                                          ts54MaxAngle=51)
        self._checkInterpTiltSeries(tsSetCtfCorr,
                                    testAcqObjDict=testAcqObjDict,
                                    anglesCountDict=self.intersectAnglesCountDictExcluded)

    def testCtfCorrection07(self):
        _, tsWithAliBin4 = self._runPrevProts(importCtf=False)
        ctfSetReStacked = self._genReStackedCtf()  # Gen a CTF estimated on a re-stacked TS
        tsSetCtfCorr = self._runCtfCorrection(tsWithAliBin4, ctfSetReStacked,
                                              tsSetMsg=self.UNMODIFIED,
                                              ctfSetMsg=self.RE_STACKED)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-51,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-57,
                                                          ts54MaxAngle=54)
        self._checkInterpTiltSeries(tsSetCtfCorr,
                                    testAcqObjDict=testAcqObjDict,
                                    anglesCountDict=self.ctfAnglesCountDictExcluded)

    def testCtfCorrection08(self):
        _, tsWithAliBin4 = self._runPrevProts(importCtf=False)
        self._excludeSetViews(tsWithAliBin4)  # Excluded some views in the TS at metadata level
        ctfSetReStacked = self._genReStackedCtf()  # Gen a CTF estimated on a re-stacked TS
        tsSetCtfCorr = self._runCtfCorrection(tsWithAliBin4, ctfSetReStacked,
                                              tsSetMsg=self.EXC_VIEWS,
                                              ctfSetMsg=self.RE_STACKED)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-51,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-54,
                                                          ts54MaxAngle=51)
        self._checkInterpTiltSeries(tsSetCtfCorr,
                                    testAcqObjDict=testAcqObjDict,
                                    anglesCountDict=self.intersectAnglesCountDictExcluded)

    def testCtfCorrection09(self):
        _, tsWithAliBin4 = self._runPrevProts(importCtf=False)
        self._excludeSetViews(tsWithAliBin4)  # Excluded some views in the TS at metadata level
        tsSetReStacked = self._runExcludeViewsProt(tsWithAliBin4)  # Re-stack the TS
        ctfSetReStacked = self._genReStackedCtf()  # Gen a CTF estimated on a re-stacked TS
        tsSetCtfCorr = self._runCtfCorrection(tsSetReStacked, ctfSetReStacked,
                                              tsSetMsg=self.RE_STACKED,
                                              ctfSetMsg=self.RE_STACKED)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-51,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-54,
                                                          ts54MaxAngle=51)
        self._checkInterpTiltSeries(tsSetCtfCorr,
                                    testAcqObjDict=testAcqObjDict,
                                    anglesCountDict=self.intersectAnglesCountDictExcluded)


class TestImodGoldBeadPicker(TestImodBase):

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTomos = cls._runImportTomograms(filesPattern='*4*.mrc')  # TS_43, TS_45, TS_54

    def testGoldBeadPicker01(self):
        beadDiameter = 12
        fiduCoords = self._runGoldBeadPicker(self.importedTomos, beadDiameter=beadDiameter)
        # Check the results
        tomoSRate = self.importedTomos.getSamplingRate()

        self.checkCoordinates(fiduCoords,
                              expectedSetSize=120,
                              expectedBoxSize=beadDiameter,
                              expectedSRate=tomoSRate)
