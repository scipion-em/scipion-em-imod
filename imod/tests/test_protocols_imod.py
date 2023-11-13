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

from pyworkflow.tests import *
from pyworkflow.utils import path
from pwem.emlib.image import ImageHandler
import tomo

from ..protocols import *


class TestImodBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @classmethod
    def _runImportTiltSeries(cls, filesPath, pattern, voltage, magnification,
                             sphericalAberration, amplitudeContrast,
                             samplingRate, doseInitial, dosePerFrame,
                             anglesFrom=0, minAngle=0.0, maxAngle=0.0,
                             stepAngle=1.0, tiltAxisAngle=-12.5):
        cls.protImportTS = cls.newProtocol(tomo.protocols.ProtImportTs,
                                           filesPath=filesPath,
                                           filesPattern=pattern,
                                           voltage=voltage,
                                           anglesFrom=anglesFrom,
                                           magnification=magnification,
                                           sphericalAberration=sphericalAberration,
                                           amplitudeContrast=amplitudeContrast,
                                           samplingRate=samplingRate,
                                           doseInitial=doseInitial,
                                           dosePerFrame=dosePerFrame,
                                           minAngle=minAngle,
                                           maxAngle=maxAngle,
                                           stepAngle=stepAngle,
                                           tiltAxisAngle=tiltAxisAngle)
        cls.launchProtocol(cls.protImportTS)
        return cls.protImportTS

    @classmethod
    def _runImportTransformationMatrix(cls, filesPath, pattern,
                                       exclusionWords, inputSetOfTiltSeries):
        cls.protImportTM = cls.newProtocol(ProtImodImportTransformationMatrix,
                                           filesPath=filesPath,
                                           filesPattern=pattern,
                                           exclusionWords=exclusionWords,
                                           inputSetOfTiltSeries=inputSetOfTiltSeries, )
        cls.launchProtocol(cls.protImportTM)
        return cls.protImportTM

    @classmethod
    def _runXRaysEraser(cls, inputSoTS, peakCriterion, diffCriterion,
                        maximumRadius, bigDiffCriterion):
        cls.protXRaysEraser = cls.newProtocol(ProtImodXraysEraser,
                                              inputSetOfTiltSeries=inputSoTS,
                                              peakCriterion=peakCriterion,
                                              diffCriterion=diffCriterion,
                                              maximumRadius=maximumRadius,
                                              bigDiffCriterion=bigDiffCriterion)
        cls.launchProtocol(cls.protXRaysEraser)
        return cls.protXRaysEraser

    @classmethod
    def _runDoseFilter(cls, inputSoTS, initialDose,
                       inputDoseType, fixedImageDose):
        cls.protDoseFilter = cls.newProtocol(ProtImodDoseFilter,
                                             inputSetOfTiltSeries=inputSoTS,
                                             initialDose=initialDose,
                                             inputDoseType=inputDoseType,
                                             fixedImageDose=fixedImageDose)
        cls.launchProtocol(cls.protDoseFilter)
        return cls.protDoseFilter

    @classmethod
    def _runExcludeViews(cls, inputSoTS, excludeViewsFile):
        cls.protExcludeViews = cls.newProtocol(ProtImodExcludeViews,
                                               inputSetOfTiltSeries=inputSoTS,
                                               excludeViewsFile=excludeViewsFile)
        cls.launchProtocol(cls.protExcludeViews)
        return cls.protExcludeViews

    @classmethod
    def _runTSNormalization(cls, inputSoTS, binning, floatDensities,
                            modeToOutput, scaleRangeToggle, scaleRangeMax,
                            scaleRangeMin, meanSdToggle, scaleMean,
                            scaleSd, scaleMax, scaleMin):
        cls.protTSNormalization = cls.newProtocol(ProtImodTSNormalization,
                                                  inputSetOfTiltSeries=inputSoTS,
                                                  binning=binning,
                                                  floatDensities=floatDensities,
                                                  modeToOutput=modeToOutput,
                                                  scaleRangeToggle=scaleRangeToggle,
                                                  scaleRangeMax=scaleRangeMax,
                                                  scaleRangeMin=scaleRangeMin,
                                                  meanSdToggle=meanSdToggle,
                                                  scaleMean=scaleMean,
                                                  scaleSd=scaleSd,
                                                  scaleMax=scaleMax,
                                                  scaleMin=scaleMin)
        cls.launchProtocol(cls.protTSNormalization)
        return cls.protTSNormalization

    @classmethod
    def _runXcorrPrealignment(cls, inputSoTS, computeAlignmentToggle,
                              binning, rotationAngle, xmin, ymin):
        cls.protXcorr = cls.newProtocol(ProtImodXcorrPrealignment,
                                        inputSetOfTiltSeries=inputSoTS,
                                        computeAlignment=computeAlignmentToggle,
                                        binning=binning,
                                        tiltAxisAngle=rotationAngle,
                                        xmin=xmin,
                                        ymin=ymin)
        cls.launchProtocol(cls.protXcorr)
        return cls.protXcorr

    @classmethod
    def _runFiducialModels(cls, inputSoTS, twoSurfaces, fiducialRadius,
                           numberFiducial, rotationAngle,
                           shiftsNearZeroFraction) -> ProtImodFiducialModel:
        cls.protFiducialAlignment = cls.newProtocol(ProtImodFiducialModel,
                                                    inputSetOfTiltSeries=inputSoTS,
                                                    twoSurfaces=twoSurfaces,
                                                    fiducialRadius=fiducialRadius,
                                                    numberFiducial=numberFiducial,
                                                    rotationAngle=rotationAngle,
                                                    shiftsNearZeroFraction=shiftsNearZeroFraction)
        cls.launchProtocol(cls.protFiducialAlignment)
        return cls.protFiducialAlignment

    @classmethod
    def _runFiducialAlignemnt(cls, inputSoLM, twoSurfaces, rotationAngle,
                              computeAlignment, binning) -> ProtImodFiducialAlignment:
        cls.protFiducialAlignment = cls.newProtocol(ProtImodFiducialAlignment,
                                                    inputSetOfLandmarkModels=inputSoLM,
                                                    twoSurfaces=twoSurfaces,
                                                    rotationAngle=rotationAngle,
                                                    computeAlignment=computeAlignment,
                                                    binning=binning)
        cls.launchProtocol(cls.protFiducialAlignment)
        return cls.protFiducialAlignment

    @classmethod
    def _runApplyTransformationMatrix(cls, inputSoTS, binning):
        cls.protApplyTransformationMatrix = cls.newProtocol(ProtImodApplyTransformationMatrix,
                                                            inputSetOfTiltSeries=inputSoTS,
                                                            binning=binning)
        cls.launchProtocol(cls.protApplyTransformationMatrix)
        return cls.protApplyTransformationMatrix

    @classmethod
    def _runTomoReconstruction(cls, inputSoTS, tomoThickness, tomoShiftX,
                               tomoShiftZ, angleOffset, tiltAxisOffset,
                               fakeInteractionsSIRT, radialFirstParameter,
                               radialSecondParameter):
        cls.protTomoReconstruction = cls.newProtocol(ProtImodTomoReconstruction,
                                                     inputSetOfTiltSeries=inputSoTS,
                                                     tomoThickness=tomoThickness,
                                                     tomoShiftX=tomoShiftX,
                                                     tomoShiftZ=tomoShiftZ,
                                                     angleOffset=angleOffset,
                                                     tiltAxisOffset=tiltAxisOffset,
                                                     fakeInteractionsSIRT=fakeInteractionsSIRT,
                                                     radialFirstParameter=radialFirstParameter,
                                                     radialSecondParameter=radialSecondParameter)
        cls.launchProtocol(cls.protTomoReconstruction)
        return cls.protTomoReconstruction

    @classmethod
    def _runTomoNormalization(cls,
                              inputSetOfTomograms, binning, floatDensities,
                              modeToOutput, scaleRangeToggle,
                              scaleRangeMax, scaleRangeMin, meanSdToggle,
                              scaleMean, scaleSd, scaleMax, scaleMin) -> ProtImodTomoNormalization:

        cls.protTomoNormalization = cls.newProtocol(ProtImodTomoNormalization,
                                                    inputSetOfTomograms=inputSetOfTomograms,
                                                    binning=binning,
                                                    floatDensities=floatDensities,
                                                    modeToOutput=modeToOutput,
                                                    scaleRangeToggle=scaleRangeToggle,
                                                    scaleRangeMax=scaleRangeMax,
                                                    scaleRangeMin=scaleRangeMin,
                                                    meanSdToggle=meanSdToggle,
                                                    scaleMean=scaleMean,
                                                    scaleSd=scaleSd,
                                                    scaleMax=scaleMax,
                                                    scaleMin=scaleMin)
        cls.launchProtocol(cls.protTomoNormalization)
        return cls.protTomoNormalization

    @classmethod
    def _runGoldBeadPiker3D(cls, inputSetOfTomograms, beadDiameter,
                            beadsColor, minRelativeStrength, minSpacing) -> ProtImodGoldBeadPicker3d:
        cls.protGoldBeadPicker3d = cls.newProtocol(ProtImodGoldBeadPicker3d,
                                                   inputSetOfTomograms=inputSetOfTomograms,
                                                   beadDiameter=beadDiameter,
                                                   beadsColor=beadsColor,
                                                   minRelativeStrength=minRelativeStrength,
                                                   minSpacing=minSpacing)
        cls.launchProtocol(cls.protGoldBeadPicker3d)
        return cls.protGoldBeadPicker3d

    @classmethod
    def _runTomoProjection(cls, inputSetOfTomograms, minAngle,
                           maxAngle, stepAngle, rotationAxis):
        cls.protTomoProjection = cls.newProtocol(ProtImodTomoProjection,
                                                 inputSetOfTomograms=inputSetOfTomograms,
                                                 minAngle=minAngle,
                                                 maxAngle=maxAngle,
                                                 stepAngle=stepAngle,
                                                 rotationAxis=rotationAxis)
        cls.launchProtocol(cls.protTomoProjection)
        return cls.protTomoProjection

    @classmethod
    def _runImportSetOfCtfSeries(cls, filesPath, filesPattern, inputSetOfTiltSeries):
        cls.protImportSetOfCtfSeries = cls.newProtocol(ProtImodImportSetOfCtfTomoSeries,
                                                       filesPath=filesPath,
                                                       filesPattern=filesPattern,
                                                       inputSetOfTiltSeries=inputSetOfTiltSeries)
        cls.launchProtocol(cls.protImportSetOfCtfSeries)
        return cls.protImportSetOfCtfSeries

    @classmethod
    def _runCTFEstimation(cls, inputSoTS, defocusTol, expectedDefocusOrigin,
                          expectedDefocusValue, expectedDefocusFile,
                          axisAngle, interactiveMode, leftDefTol, rightDefTol,
                          tileSize, angleStep, angleRange, startFreq, endFreq,
                          extraZerosToFit, skipAstigmaticViews, searchAstigmatism,
                          findAstigPhaseCutonToggle, phaseShiftAstigmatism,
                          cutOnFrequencyAstigmatism, minimumViewsAstigmatism,
                          minimumViewsPhaseShift, numberSectorsAstigmatism,
                          maximumAstigmatism):
        cls.protCTFEstimation = cls.newProtocol(ProtImodAutomaticCtfEstimation,
                                                inputSet=inputSoTS,
                                                defocusTol=defocusTol,
                                                expectedDefocusOrigin=expectedDefocusOrigin,
                                                expectedDefocusValue=expectedDefocusValue,
                                                expectedDefocusFile=expectedDefocusFile,
                                                axisAngle=axisAngle,
                                                interactiveMode=interactiveMode,
                                                leftDefTol=leftDefTol,
                                                rightDefTol=rightDefTol,
                                                tileSize=tileSize,
                                                angleStep=angleStep,
                                                angleRange=angleRange,
                                                startFreq=startFreq,
                                                endFreq=endFreq,
                                                extraZerosToFit=extraZerosToFit,
                                                skipAstigmaticViews=skipAstigmaticViews,
                                                searchAstigmatism=searchAstigmatism,
                                                findAstigPhaseCutonToggle=findAstigPhaseCutonToggle,
                                                phaseShiftAstigmatism=phaseShiftAstigmatism,
                                                cutOnFrequencyAstigmatism=cutOnFrequencyAstigmatism,
                                                minimumViewsAstigmatism=minimumViewsAstigmatism,
                                                minimumViewsPhaseShift=minimumViewsPhaseShift,
                                                numberSectorsAstigmatism=numberSectorsAstigmatism,
                                                maximumAstigmatism=maximumAstigmatism)
        cls.launchProtocol(cls.protCTFEstimation)
        return cls.protCTFEstimation

    @classmethod
    def _runCTFCorrection(cls, inputSetOfTiltSeries, inputSetOfCtfTomoSeries,
                          defocusTol, interpolationWidth):
        cls.protCTFCorrection = cls.newProtocol(ProtImodCtfCorrection,
                                                inputSetOfTiltSeries=inputSetOfTiltSeries,
                                                inputSetOfCtfTomoSeries=inputSetOfCtfTomoSeries,
                                                defocusTol=defocusTol,
                                                interpolationWidth=interpolationWidth)
        cls.launchProtocol(cls.protCTFCorrection)
        return cls.protCTFCorrection


class TestImodReconstructionWorkflow(TestImodBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.inputDataSet = DataSet.getDataSet('tomo-em')
        cls.inputSoTS = cls.inputDataSet.getFile('ts1')

        cls.excludeViewsFile = cls.inputDataSet.getFile('excludeViewsFile')

        cls.inputTMFolder = os.path.split(cls.inputDataSet.getFile('tm1'))[0]

        cls.excludeViewsOutputSizes = {'a': 57, 'b': 56}

        cls.binningTsNormalization = 2

        cls.binningPrealignment = 2

        cls.binningFiducialAlignment = 2

        cls.binningApplyTransformMatrix = 2

        cls.thicknessTomo = 100

        cls.binningTomoNormalization = 2

        cls.protImportTS = cls._runImportTiltSeries(filesPath=os.path.split(cls.inputSoTS)[0],
                                                    pattern="BB{TS}.st",
                                                    anglesFrom=0,
                                                    voltage=300,
                                                    magnification=105000,
                                                    sphericalAberration=2.7,
                                                    amplitudeContrast=0.1,
                                                    samplingRate=20.2,
                                                    doseInitial=0,
                                                    dosePerFrame=3.0,
                                                    minAngle=-55,
                                                    maxAngle=65.0,
                                                    stepAngle=2.0)

        cls.protImportTM = cls._runImportTransformationMatrix(filesPath=cls.inputTMFolder,
                                                              pattern="BB*.prexg",
                                                              exclusionWords='',
                                                              inputSetOfTiltSeries=cls.protImportTS.outputTiltSeries)

        cls.protXRaysEraser = cls._runXRaysEraser(inputSoTS=cls.protImportTS.outputTiltSeries,
                                                  peakCriterion=8.0,
                                                  diffCriterion=6.0,
                                                  maximumRadius=4.2,
                                                  bigDiffCriterion=19)

        cls.protDoseFilter = cls._runDoseFilter(inputSoTS=cls.protXRaysEraser.TiltSeries,
                                                initialDose=0,
                                                inputDoseType=1,
                                                fixedImageDose=2.0)

        cls.protExcludeViews = cls._runExcludeViews(inputSoTS=cls.protXRaysEraser.TiltSeries,
                                                    excludeViewsFile=cls.excludeViewsFile)

        cls.protTSNormalization = cls._runTSNormalization(inputSoTS=cls.protDoseFilter.TiltSeries,
                                                          binning=cls.binningTsNormalization,
                                                          floatDensities=0,
                                                          modeToOutput=0,
                                                          scaleRangeToggle=1,
                                                          scaleRangeMax=255,
                                                          scaleRangeMin=0,
                                                          meanSdToggle=1,
                                                          scaleMean=0,
                                                          scaleSd=1,
                                                          scaleMax=255,
                                                          scaleMin=0)

        cls.protXcorr = cls._runXcorrPrealignment(inputSoTS=cls.protDoseFilter.TiltSeries,
                                                  computeAlignmentToggle=0,
                                                  binning=cls.binningPrealignment,
                                                  rotationAngle=-12.5,
                                                  xmin=10,
                                                  ymin=10)

        cls.protFiducialModels = cls._runFiducialModels(inputSoTS=cls.protXcorr.TiltSeries,
                                                        twoSurfaces=0,
                                                        fiducialRadius=4.95,
                                                        numberFiducial=25,
                                                        rotationAngle=-12.5,
                                                        shiftsNearZeroFraction=0.2)

        cls.protFiducialAlignment = cls._runFiducialAlignemnt(inputSoLM=cls.protFiducialModels.FiducialModelGaps,
                                                              twoSurfaces=0,
                                                              rotationAngle=-12.5,
                                                              computeAlignment=0,
                                                              binning=cls.binningFiducialAlignment)

        cls.protApplyTransformationMatrix = \
            cls._runApplyTransformationMatrix(inputSoTS=cls.protFiducialAlignment.TiltSeries,
                                              binning=cls.binningApplyTransformMatrix)

        cls.protTomoReconstruction = \
            cls._runTomoReconstruction(inputSoTS=cls.protFiducialAlignment.TiltSeries,
                                       tomoThickness=cls.thicknessTomo,
                                       tomoShiftX=0.0,
                                       tomoShiftZ=0.0,
                                       angleOffset=0.0,
                                       tiltAxisOffset=0.0,
                                       fakeInteractionsSIRT=0,
                                       radialFirstParameter=0.35,
                                       radialSecondParameter=0.035)

        cls.protTomoNormalization = \
            cls._runTomoNormalization(inputSetOfTomograms=cls.protTomoReconstruction.Tomograms,
                                      binning=cls.binningTomoNormalization,
                                      floatDensities=0,
                                      modeToOutput=0,
                                      scaleRangeToggle=1,
                                      scaleRangeMax=255,
                                      scaleRangeMin=0,
                                      meanSdToggle=1,
                                      scaleMean=0,
                                      scaleSd=1,
                                      scaleMax=255,
                                      scaleMin=0)

        cls.protGoldBeadPicker3D = \
            cls._runGoldBeadPiker3D(inputSetOfTomograms=cls.protTomoReconstruction.Tomograms,
                                    beadDiameter=10,
                                    beadsColor=0,
                                    minRelativeStrength=0.05,
                                    minSpacing=0.9)

        cls.protTomoProjection = \
            cls._runTomoProjection(inputSetOfTomograms=cls.protTomoNormalization.Tomograms,
                                   minAngle=-60.0,
                                   maxAngle=60.0,
                                   stepAngle=2.0,
                                   rotationAxis=1)

    def test_importTMOutput(self):

        tseries = self.protImportTM.TiltSeries
        self.assertSetSize(tseries, size=2)
        self.assertTrue(tseries.hasAlignment(), "Tilt series does not have alignment flag")
        for ts in tseries:
            self.assertTrue(ts.getFirstItem().hasTransform())

    def test_doseFilterOutputTS(self):
        ts = self.protDoseFilter.TiltSeries
        self.assertSetSize(ts)

        tsId = ts.getFirstItem().getTsId()

        self.assertTrue(os.path.exists(os.path.join(self.protDoseFilter._getExtraPath(tsId),
                                                    "BB" + tsId + ".st")))

    def test_xRaysEraserOutputTS(self):

        ts = self.protXRaysEraser.TiltSeries
        self.assertSetSize(ts)

        tsId = ts.getFirstItem().getTsId()

        self.assertTrue(os.path.exists(os.path.join(self.protXRaysEraser._getExtraPath(tsId),
                                                    "BB" + tsId + ".st")))

    def test_excludeViewsOutputTS(self):

        ts = self.protExcludeViews.TiltSeries
        self.assertSetSize(ts)

        tsId = ts.getFirstItem().getTsId()

        self.assertTrue(os.path.exists(os.path.join(self.protExcludeViews._getExtraPath(tsId),
                                                    "BB" + tsId + ".st")))

        for index, tsOut in enumerate(ts):
            self.assertEqual(tsOut.getSize(), self.excludeViewsOutputSizes[tsOut.getTsId()])

    def test_normalizationOutputTS(self):

        ts = self.protTSNormalization.TiltSeries
        self.assertSetSize(ts)

        tsId = ts.getFirstItem().getTsId()

        self.assertTrue(os.path.exists(os.path.join(self.protTSNormalization._getExtraPath(tsId),
                                                    "BB" + tsId + ".st")))

        inSamplingRate = self.protTSNormalization.inputSetOfTiltSeries.get().getSamplingRate()
        outSamplingRate = ts.getSamplingRate()

        self.assertEqual(inSamplingRate * self.binningTsNormalization, outSamplingRate)

    def test_prealignmentOutputTS(self):

        ts = self.protXcorr.TiltSeries
        self.assertSetSize(ts)

        tsId = ts.getFirstItem().getTsId()
        outputLocation = os.path.join(self.protXcorr._getExtraPath(tsId),
                                      "BB" + tsId + ".st")

        self.assertTrue(os.path.exists(outputLocation))

        self.assertIsNotNone(ts.getFirstItem().getFirstItem().getTransform())

    def test_prealignmentOutputInterpolatedTS(self):

        ts = self.protXcorr.InterpolatedTiltSeries
        self.assertSetSize(ts)
        self.assertFalse(ts.hasAlignment(), "Tilt series does not have alignment flag canceled")

        tsId = ts.getFirstItem().getTsId()
        outputLocation = os.path.join(self.protXcorr._getExtraPath(tsId),
                                      "BB" + tsId + ".st")

        self.assertTrue(os.path.exists(outputLocation))

        inSamplingRate = self.protXcorr.inputSetOfTiltSeries.get().getSamplingRate()
        outSamplingRate = ts.getSamplingRate()

        self.assertEqual(inSamplingRate * self.binningPrealignment, outSamplingRate)

    def test_fiducialModelstOutputFiducialModelGaps(self):

        output = self.protFiducialModels.FiducialModelGaps
        self.assertSetSize(output, size=2)

        tsId = output.getFirstItem().getTsId()
        outputLocationImod = os.path.join(self.protFiducialModels._getExtraPath(tsId),
                                          "BB" + tsId + "_gaps.fid")
        outputLocationScipion = os.path.join(self.protFiducialModels._getExtraPath(tsId),
                                             "BB" + tsId + "_gaps.sfid")

        self.assertTrue(os.path.exists(outputLocationImod))
        self.assertTrue(os.path.exists(outputLocationScipion))

    def test_fiducialAlignmentOutputTS(self):

        output = self.protFiducialAlignment.TiltSeries
        self.assertSetSize(output, size=2)
        self.assertTrue(output.hasAlignment(), "Fiducial alignment Tilt series does not have alignment flag")

        tsId = output.getFirstItem().getTsId()
        outputLocation = os.path.join(self.protFiducialAlignment._getExtraPath(tsId),
                                      "BB" + tsId + ".st")

        self.assertTrue(os.path.exists(outputLocation))

        self.assertIsNotNone(output.getFirstItem().getFirstItem().getTransform())

        # Interpolated
        output = self.protFiducialAlignment.InterpolatedTiltSeries
        self.assertSetSize(output, size=2)
        self.assertFalse(output.hasAlignment(), "Interpolated Tilt series does have alignment flag")

        tsId = output.getFirstItem().getTsId()
        outputLocation = os.path.join(self.protFiducialAlignment._getExtraPath(tsId),
                                      "BB" + tsId + ".st")

        self.assertTrue(os.path.exists(outputLocation))

        inputSoTS = self.protFiducialAlignment.inputSetOfLandmarkModels.get().getSetOfTiltSeries()

        inSamplingRate = inputSoTS.getSamplingRate()
        outSamplingRate = output.getSamplingRate()

        self.assertEqual(inSamplingRate * self.binningFiducialAlignment, outSamplingRate)

        output = self.protFiducialAlignment.FiducialModelNoGaps
        self.assertSetSize(output, size=2)

        tsId = output.getFirstItem().getTsId()
        outputLocation = os.path.join(self.protFiducialAlignment._getExtraPath(tsId),
                                      "BB" + tsId + "_noGaps.sfid")

        self.assertTrue(os.path.exists(outputLocation))

        # Output coordinates
        output = self.protFiducialAlignment.TiltSeriesCoordinates
        self.assertSetSize(output)

        tsId = output.getFirstItem().getTsId()
        outputLocation = os.path.join(self.protFiducialAlignment._getExtraPath(tsId),
                                      "BB" + tsId + "_fid.xyz")

        self.assertTrue(os.path.exists(outputLocation))

        tolerance = 20
        expectedSize = 35

        self.assertTrue(abs(output.getSize() - expectedSize) <= tolerance)

    def test_applyTransformationMatrixOutputInterpolatedTS(self):

        output = self.protApplyTransformationMatrix.InterpolatedTiltSeries
        self.assertSetSize(output, size=2)
        self.assertFalse(output.hasAlignment(), "Tilt series interpolated does have alignment flag")

        tsId = output.getFirstItem().getTsId()
        outputLocation = os.path.join(self.protApplyTransformationMatrix._getExtraPath(tsId),
                                      "BB" + tsId + ".st")

        self.assertTrue(os.path.exists(outputLocation))

        inSamplingRate = self.protApplyTransformationMatrix.inputSetOfTiltSeries.get().getSamplingRate()
        outSamplingRate = output.getSamplingRate()

        self.assertEqual(inSamplingRate * self.binningApplyTransformMatrix, outSamplingRate)

    def test_tomoReconstructionOutputTomogram(self):

        output = self.protTomoReconstruction.Tomograms
        self.assertSetSize(output, size=2)

        tomoId = self.protTomoReconstruction.inputSetOfTiltSeries.get().getFirstItem().getTsId()
        outputLocation = os.path.join(self.protTomoReconstruction._getExtraPath(tomoId),
                                      "BB" + tomoId + ".mrc")

        self.assertTrue(os.path.exists(outputLocation))

        # Dimensions
        ih = ImageHandler()
        outputDimensions = ih.getDimensions(output.getFirstItem())

        self.assertEqual(outputDimensions, (512, 512, self.thicknessTomo, 1))

    def test_goldBeadPeaker3DOutput(self):

        output = self.protGoldBeadPicker3D.Coordinates3D
        self.assertSetSize(output, size=52, diffDelta=30)

        tomoId = self.protGoldBeadPicker3D.inputSetOfTomograms.get().getFirstItem().getTsId()
        location = self.protGoldBeadPicker3D._getExtraPath("BB" + tomoId)

        self.assertTrue(os.path.exists(os.path.join(location,
                                                    "BB" + tomoId + ".mod")))

    def test_tomoNormalizationOutput(self):

        output = self.protTomoNormalization.Tomograms
        self.assertSetSize(output, size=2)

        location = self.protTomoNormalization.inputSetOfTomograms.get().getFirstItem().getFileName()
        tomoId = pwutils.removeBaseExt(location)

        self.assertTrue(os.path.exists(os.path.join(self.protTomoNormalization._getExtraPath(tomoId),
                                                    tomoId + ".mrc")))

        inSamplingRate = self.protTomoNormalization.inputSetOfTomograms.get().getSamplingRate()
        outSamplingRate = output.getSamplingRate()

        self.assertEqual(inSamplingRate * self.binningTomoNormalization, outSamplingRate)

    def test_tomoProjectionOutputTiltSeriesSize(self):

        output = self.protTomoProjection.TiltSeries
        self.assertSetSize(output, size=2)

        # Dimensions
        ih = ImageHandler()
        outputDimensions = ih.getDimensions(output.getFirstItem().getFirstItem().getFileName())

        self.assertEqual(outputDimensions, (256, 256, 61, 1))

        # Sampling rate
        inSamplingRate = self.protTomoProjection.inputSetOfTomograms.get().getSamplingRate()
        outSamplingRate = output.getSamplingRate()

        self.assertEqual(inSamplingRate, outSamplingRate)


class TestImodCTFCorrectionWorkflow(TestImodBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.inputDataSet = DataSet.getDataSet('tutorialDataImodCTF')
        cls.inputSoTS = cls.inputDataSet.getFile('tsCtf1')

        cls.inputCtfFile = cls.inputDataSet.getFile('inputCtfFile')

        # Create links to the input tilt-series and its associated mdoc file to test the protocols with a set of two
        # elements to make the tests more robust
        linkTs = os.path.join(os.path.split(cls.inputSoTS)[0],
                              "WTI042413_1series4_copy.st")

        if not os.path.exists(linkTs):
            path.createLink(cls.inputSoTS, linkTs)

        cls.protImportTS = cls._runImportTiltSeries(filesPath=os.path.split(cls.inputSoTS)[0],
                                                    pattern="*.mdoc",
                                                    anglesFrom=0,
                                                    voltage=300,
                                                    magnification=50000,
                                                    sphericalAberration=2.7,
                                                    amplitudeContrast=0.07,
                                                    samplingRate=6.73981,
                                                    doseInitial=0,
                                                    dosePerFrame=3.0)

        cls.protImportSetOfCtfSeries = \
            cls._runImportSetOfCtfSeries(filesPath=os.path.split(cls.inputCtfFile)[0],
                                         filesPattern='WTI042413_1series4.defocus',
                                         inputSetOfTiltSeries=cls.protImportTS.outputTiltSeries)

        cls.protCTFEstimation = cls._runCTFEstimation(inputSoTS=cls.protImportTS.outputTiltSeries,
                                                      defocusTol=200.0,
                                                      expectedDefocusOrigin=0,
                                                      expectedDefocusValue=6000,
                                                      expectedDefocusFile="",
                                                      axisAngle=0.0,
                                                      interactiveMode=1,
                                                      leftDefTol=2000.0,
                                                      rightDefTol=2000.0,
                                                      tileSize=256,
                                                      angleStep=2.0,
                                                      angleRange=20.0,
                                                      startFreq=0.0,
                                                      endFreq=0.0,
                                                      extraZerosToFit=0.0,
                                                      skipAstigmaticViews=1,
                                                      searchAstigmatism=1,
                                                      findAstigPhaseCutonToggle=1,
                                                      phaseShiftAstigmatism=0,
                                                      cutOnFrequencyAstigmatism=0,
                                                      minimumViewsAstigmatism=3,
                                                      minimumViewsPhaseShift=1,
                                                      numberSectorsAstigmatism=36,
                                                      maximumAstigmatism=1.2)

        cls.protCTFCorrection = cls._runCTFCorrection(inputSetOfTiltSeries=cls.protImportTS.outputTiltSeries,
                                                      inputSetOfCtfTomoSeries=cls.protCTFEstimation.CTFTomoSeries,
                                                      defocusTol=200,
                                                      interpolationWidth=15)

    def test_importCtfTomoSeriesOutput(self):
        self.assertSetSize(self.protImportSetOfCtfSeries.CTFTomoSeries, size=1)

    def test_ctfEstimationOutputSize(self):
        self.assertSetSize(self.protCTFEstimation.CTFTomoSeries, size=2)

    def test_ctfEstimationOutputDefocusFile(self):
        for ts in self.protCTFEstimation.inputSet.get():
            tsId = ts.getTsId()
            defocusFile = os.path.join(self.protCTFEstimation._getExtraPath(tsId),
                                       '%s.defocus' % tsId)

            self.assertTrue(os.path.exists(defocusFile))

    def test_ctfCorrectionOutput(self):

        output = self.protCTFCorrection.TiltSeries
        self.assertSetSize(output, size=2)

        for ts in output:
            tsId = ts.getTsId()
            outputLocation = os.path.join(self.protCTFCorrection._getExtraPath(tsId),
                                          '%s.st' % tsId)

            self.assertTrue(os.path.exists(outputLocation))

        self.assertTrue(output.interpolated(), "Tilt series does not have interpolated flag")
