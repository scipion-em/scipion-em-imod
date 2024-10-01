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
from pwem.emlib.image import ImageHandler as ih
from tomo.protocols import ProtImportTs, ProtImportTsCTF

from imod.protocols import *


class TestImodBase(BaseTest):
    atsId = 'BBa'
    btsId = 'BBb'

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @classmethod
    def _runImportTiltSeries(cls, minAngle=0.0, maxAngle=0.0, stepAngle=1.0,
                             tiltAxisAngle=-12.5, **kwargs):
        cls.protImportTS = cls.newProtocol(ProtImportTs,
                                           minAngle=minAngle,
                                           maxAngle=maxAngle,
                                           stepAngle=stepAngle,
                                           tiltAxisAngle=tiltAxisAngle,
                                           **kwargs)
        cls.launchProtocol(cls.protImportTS)
        return cls.protImportTS

    @classmethod
    def _runImportTransformationMatrix(cls, **kwargs):
        cls.protImportTM = cls.newProtocol(ProtImodImportTransformationMatrix,
                                           **kwargs)
        cls.launchProtocol(cls.protImportTM)
        return cls.protImportTM

    @classmethod
    def _runXRaysEraser(cls, **kwargs):
        cls.protXRaysEraser = cls.newProtocol(ProtImodXraysEraser, **kwargs)
        cls.launchProtocol(cls.protXRaysEraser)
        return cls.protXRaysEraser

    @classmethod
    def _runDoseFilter(cls, **kwargs):
        cls.protDoseFilter = cls.newProtocol(ProtImodDoseFilter, **kwargs)
        cls.launchProtocol(cls.protDoseFilter)
        return cls.protDoseFilter

    @classmethod
    def _runTSNormalization(cls, **kwargs):
        cls.protTSNormalization = cls.newProtocol(ProtImodTsNormalization, **kwargs)
        cls.launchProtocol(cls.protTSNormalization)
        return cls.protTSNormalization

    @classmethod
    def _runXcorrPrealignment(cls, **kwargs):
        cls.protXcorr = cls.newProtocol(ProtImodXcorrPrealignment, **kwargs)
        cls.launchProtocol(cls.protXcorr)
        return cls.protXcorr

    @classmethod
    def _runFiducialModels(cls, **kwargs):
        cls.protFiducialAlignment = cls.newProtocol(ProtImodFiducialModel,
                                                    typeOfModel=0, **kwargs)
        cls.launchProtocol(cls.protFiducialAlignment)
        return cls.protFiducialAlignment

    @classmethod
    def _runFiducialModelsPT(cls, **kwargs):
        cls.protFiducialAlignment = cls.newProtocol(ProtImodFiducialModel,
                                                    typeOfModel=1, **kwargs)
        cls.launchProtocol(cls.protFiducialAlignment)
        return cls.protFiducialAlignment

    @classmethod
    def _runFiducialAlignemnt(cls, **kwargs):
        cls.protFiducialAlignment = cls.newProtocol(ProtImodFiducialAlignment,
                                                    **kwargs)
        cls.launchProtocol(cls.protFiducialAlignment)
        return cls.protFiducialAlignment

    @classmethod
    def _runApplyTransformationMatrix(cls, **kwargs):
        cls.protApplyTransformationMatrix = cls.newProtocol(ProtImodApplyTransformationMatrix,
                                                            **kwargs)
        cls.launchProtocol(cls.protApplyTransformationMatrix)
        return cls.protApplyTransformationMatrix

    @classmethod
    def _runTomoReconstruction(cls, **kwargs):
        cls.protTomoReconstruction = cls.newProtocol(ProtImodTomoReconstruction,
                                                     **kwargs)
        cls.launchProtocol(cls.protTomoReconstruction)
        return cls.protTomoReconstruction

    @classmethod
    def _runTomoNormalization(cls, **kwargs):
        cls.protTomoNormalization = cls.newProtocol(ProtImodTomoNormalization,
                                                    **kwargs)
        cls.launchProtocol(cls.protTomoNormalization)
        return cls.protTomoNormalization

    @classmethod
    def _runGoldBeadPiker3D(cls, **kwargs):
        cls.protGoldBeadPicker3d = cls.newProtocol(ProtImodGoldBeadPicker3d,
                                                   **kwargs)
        cls.launchProtocol(cls.protGoldBeadPicker3d)
        return cls.protGoldBeadPicker3d

    @classmethod
    def _runTomoProjection(cls, **kwargs):
        cls.protTomoProjection = cls.newProtocol(ProtImodTomoProjection,
                                                 **kwargs)
        cls.launchProtocol(cls.protTomoProjection)
        return cls.protTomoProjection

    @classmethod
    def _runImportSetOfCtfSeries(cls, **kwargs):
        cls.protImportSetOfCtfSeries = cls.newProtocol(ProtImportTsCTF,
                                                       importFrom=1,  # imod
                                                       **kwargs)
        cls.launchProtocol(cls.protImportSetOfCtfSeries)
        return cls.protImportSetOfCtfSeries

    @classmethod
    def _runCTFEstimation(cls, **kwargs):
        cls.protCTFEstimation = cls.newProtocol(ProtImodAutomaticCtfEstimation,
                                                **kwargs)
        cls.launchProtocol(cls.protCTFEstimation)
        return cls.protCTFEstimation

    @classmethod
    def _runCTFCorrection(cls, **kwargs):
        cls.protCTFCorrection = cls.newProtocol(ProtImodCtfCorrection,
                                                **kwargs)
        cls.launchProtocol(cls.protCTFCorrection)
        return cls.protCTFCorrection


class TestImodReconstructionWorkflow(TestImodBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.inputDataSet = DataSet.getDataSet('tomo-em')
        cls.inputSoTS = cls.inputDataSet.getFile('ts1')
        cls.inputTMFolder = os.path.split(cls.inputDataSet.getFile('tm1'))[0]
        cls.binningTsNormalization = 2
        cls.binningPrealignment = 2
        cls.binningFiducialAlignment = 2
        cls.binningApplyTransformMatrix = 2
        cls.thicknessTomo = 100
        cls.binningTomoNormalization = 2

        cls.protImportTS = cls._runImportTiltSeries(filesPath=os.path.split(cls.inputSoTS)[0],
                                                    filesPattern="{TS}.st",
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
                                                              filesPattern="BB*.prexg",
                                                              inputSetOfTiltSeries=cls.protImportTS.outputTiltSeries)

        cls.protXRaysEraser = cls._runXRaysEraser(inputSetOfTiltSeries=cls.protImportTS.outputTiltSeries,
                                                  peakCriterion=8.0,
                                                  diffCriterion=6.0,
                                                  maximumRadius=4.2,
                                                  bigDiffCriterion=19)

        cls.protDoseFilter = cls._runDoseFilter(inputSetOfTiltSeries=cls.protXRaysEraser.TiltSeries,
                                                initialDose=0,
                                                inputDoseType=1,
                                                fixedImageDose=2.0)

        cls.protTSNormalization = cls._runTSNormalization(inputSetOfTiltSeries=cls.protDoseFilter.TiltSeries,
                                                          binning=cls.binningTsNormalization,
                                                          floatDensities=2,
                                                          modeToOutput=0)

        cls.protXcorr = cls._runXcorrPrealignment(inputSetOfTiltSeries=cls.protDoseFilter.TiltSeries,
                                                  computeAlignment=True,
                                                  binning=cls.binningPrealignment,
                                                  rotationAngle=-12.5,
                                                  xmin=10,
                                                  ymin=10)

        cls.protFiducialModels = cls._runFiducialModels(inputSetOfTiltSeries=cls.protXcorr.TiltSeries,
                                                        twoSurfaces=False,
                                                        fiducialDiameter=10,
                                                        numberFiducial=25,
                                                        rotationAngle=-12.5)

        cls.protFiducialModelsPT = cls._runFiducialModelsPT(inputSetOfTiltSeries=cls.protXcorr.TiltSeries,
                                                            sizeOfPatches="100 100")

        cls.protFiducialAlignment = cls._runFiducialAlignemnt(inputSetOfLandmarkModels=cls.protFiducialModels.FiducialModelGaps,
                                                              twoSurfaces=False,
                                                              rotationAngle=-12.5,
                                                              computeAlignment=True,
                                                              binning=cls.binningFiducialAlignment)

        cls.protApplyTransformationMatrix = \
            cls._runApplyTransformationMatrix(inputSetOfTiltSeries=cls.protFiducialAlignment.TiltSeries,
                                              binning=cls.binningApplyTransformMatrix)

        cls.protTomoReconstruction = \
            cls._runTomoReconstruction(inputSetOfTiltSeries=cls.protFiducialAlignment.TiltSeries,
                                       tomoThickness=cls.thicknessTomo)

        cls.protTomoNormalization = \
            cls._runTomoNormalization(inputSetOfTomograms=cls.protTomoReconstruction.Tomograms,
                                      binning=cls.binningTomoNormalization,
                                      floatDensities=0,
                                      modeToOutput=0)

        cls.protGoldBeadPicker3D = \
            cls._runGoldBeadPiker3D(inputSetOfTomograms=cls.protTomoReconstruction.Tomograms,
                                    beadDiameter=10,
                                    beadsColor=0)

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
            self.assertTrue(ts.hasAlignment())

    def test_doseFilterOutputTS(self):
        ts = self.protDoseFilter.TiltSeries
        self.assertSetSize(ts)

        tsId = ts.getFirstItem().getTsId()

        self.assertTrue(os.path.exists(self.protDoseFilter._getExtraPath(tsId,
                                                                         tsId + ".mrcs")))

    def test_xRaysEraserOutputTS(self):
        ts = self.protXRaysEraser.TiltSeries
        self.assertSetSize(ts)

        tsId = ts.getFirstItem().getTsId()

        self.assertTrue(os.path.exists(self.protXRaysEraser._getExtraPath(tsId,
                                                                          tsId + ".mrcs")))

    def test_normalizationOutputTS(self):
        ts = self.protTSNormalization.TiltSeries
        self.assertSetSize(ts)

        tsId = ts.getFirstItem().getTsId()

        self.assertTrue(os.path.exists(self.protTSNormalization._getExtraPath(tsId,
                                                                              tsId + ".mrcs")))

        inSamplingRate = self.protTSNormalization.inputSetOfTiltSeries.get().getSamplingRate()
        outSamplingRate = ts.getSamplingRate()

        self.assertEqual(inSamplingRate * self.binningTsNormalization, outSamplingRate)

    def test_prealignmentOutputTS(self):
        ts = self.protXcorr.TiltSeries
        self.assertSetSize(ts)

        self.assertTrue(ts.hasAlignment(), "Tilt series does not have alignment flag")

    def test_prealignmentOutputInterpolatedTS(self):
        ts = self.protXcorr.InterpolatedTiltSeries
        self.assertSetSize(ts)
        self.assertFalse(ts.hasAlignment(), "Tilt series has alignment flag")

        tsId = ts.getFirstItem().getTsId()
        outputLocation = self.protXcorr._getExtraPath(tsId, tsId + ".mrcs")

        self.assertTrue(os.path.exists(outputLocation))

        inSamplingRate = self.protXcorr.inputSetOfTiltSeries.get().getSamplingRate()
        outSamplingRate = ts.getSamplingRate()

        self.assertEqual(inSamplingRate * self.binningPrealignment, outSamplingRate)

    def test_fiducialModelstOutputFiducialModelGaps(self):
        output = self.protFiducialModels.FiducialModelGaps
        self.assertSetSize(output, size=2)

        tsId = output.getFirstItem().getTsId()
        outputLocationImod = self.protFiducialModels._getExtraPath(tsId,
                                                                   tsId + "_gaps.fid")
        outputLocationScipion = self.protFiducialModels._getExtraPath(tsId,
                                                                      tsId + "_gaps.sfid")

        self.assertTrue(os.path.exists(outputLocationImod))
        self.assertTrue(os.path.exists(outputLocationScipion))

    def test_fiducialAlignmentOutputTS(self):
        output = self.protFiducialAlignment.TiltSeries
        self.assertSetSize(output, size=2)
        self.assertTrue(output.hasAlignment(),
                        "Fiducial alignment Tilt series does not have alignment flag")

        # Interpolated
        output = self.protFiducialAlignment.InterpolatedTiltSeries
        self.assertSetSize(output, size=2)
        self.assertFalse(output.hasAlignment(), "Interpolated Tilt series does have alignment flag")

        tsId = output.getFirstItem().getTsId()
        outputLocation = self.protFiducialAlignment._getExtraPath(tsId,
                                                                  tsId + ".mrcs")

        self.assertTrue(os.path.exists(outputLocation))

        inputSoTS = self.protFiducialAlignment.inputSetOfLandmarkModels.get().getSetOfTiltSeries()

        inSamplingRate = inputSoTS.getSamplingRate()
        outSamplingRate = output.getSamplingRate()

        self.assertEqual(inSamplingRate * self.binningFiducialAlignment, outSamplingRate)

        output = self.protFiducialAlignment.FiducialModelNoGaps
        self.assertSetSize(output, size=2)

        tsId = output.getFirstItem().getTsId()
        outputLocation = self.protFiducialAlignment._getExtraPath(tsId,
                                                                  tsId + "_noGaps.sfid")

        self.assertTrue(os.path.exists(outputLocation))

        # Output coordinates
        output = self.protFiducialAlignment.TiltSeriesCoordinates
        self.assertSetSize(output)

        tsId = output.getFirstItem().getTsId()
        outputLocation = self.protFiducialAlignment._getExtraPath(tsId,
                                                                  tsId + "_fid.xyz")

        self.assertTrue(os.path.exists(outputLocation))

        tolerance = 20
        expectedSize = 35

        self.assertTrue(abs(output.getSize() - expectedSize) <= tolerance)

    def test_applyTransformationMatrixOutputInterpolatedTS(self):
        output = self.protApplyTransformationMatrix.InterpolatedTiltSeries
        self.assertSetSize(output, size=2)
        self.assertFalse(output.hasAlignment(), "Tilt series interpolated has alignment flag")

        tsId = output.getFirstItem().getTsId()
        outputLocation = self.protApplyTransformationMatrix._getExtraPath(tsId,
                                                                          tsId + ".mrcs")

        self.assertTrue(os.path.exists(outputLocation))

        inSamplingRate = self.protApplyTransformationMatrix.inputSetOfTiltSeries.get().getSamplingRate()
        outSamplingRate = output.getSamplingRate()

        self.assertEqual(inSamplingRate * self.binningApplyTransformMatrix, outSamplingRate)

    def test_tomoReconstructionOutputTomogram(self):
        output = self.protTomoReconstruction.Tomograms
        self.assertSetSize(output, size=2)

        tomoId = self.protTomoReconstruction.inputSetOfTiltSeries.get().getFirstItem().getTsId()
        outputLocation = self.protTomoReconstruction._getExtraPath(tomoId,
                                                                   tomoId + ".mrc")

        self.assertTrue(os.path.exists(outputLocation))

        # Dimensions
        outputDimensions = ih.getDimensions(output.getFirstItem())

        self.assertEqual(outputDimensions, (512, 512, self.thicknessTomo, 1))

    def test_goldBeadPeaker3DOutput(self):
        output = self.protGoldBeadPicker3D.Coordinates3D
        self.assertSetSize(output, size=52, diffDelta=30)

        tomoId = self.protGoldBeadPicker3D.inputSetOfTomograms.get().getFirstItem().getTsId()
        location = self.protGoldBeadPicker3D._getExtraPath(tomoId)

        self.assertTrue(os.path.exists(os.path.join(location, tomoId + ".mod")))

    def test_tomoNormalizationOutput(self):
        output = self.protTomoNormalization.Tomograms
        self.assertSetSize(output, size=2)

        location = self.protTomoNormalization.inputSetOfTomograms.get().getFirstItem().getFileName()
        tomoId = pwutils.removeBaseExt(location)

        self.assertTrue(os.path.exists(self.protTomoNormalization._getExtraPath(tomoId,
                                                                                tomoId + ".mrc")))

        inSamplingRate = self.protTomoNormalization.inputSetOfTomograms.get().getSamplingRate()
        outSamplingRate = output.getSamplingRate()

        self.assertEqual(inSamplingRate * self.binningTomoNormalization, outSamplingRate)

    def test_tomoProjectionOutputTiltSeriesSize(self):

        output = self.protTomoProjection.TiltSeries
        self.assertSetSize(output, size=2)

        # Dimensions
        outputDimensions = ih.getDimensions(output.getFirstItem().getFirstItem().getFileName())

        self.assertEqual(outputDimensions, (256, 256, 1, 61))

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
                              "WTI042413_1series4_copy.mrc")

        if not os.path.exists(linkTs):
            path.createLink(cls.inputSoTS, linkTs)

        cls.protImportTS = cls._runImportTiltSeries(filesPath=os.path.split(cls.inputSoTS)[0],
                                                    filesPattern="*.mdoc",
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

        cls.protCTFEstimation = cls._runCTFEstimation(inputSetOfTiltSeries=cls.protImportTS.outputTiltSeries,
                                                      defocusTol=200,
                                                      expectedDefocusOrigin=0,
                                                      expectedDefocusValue=6000,
                                                      angleStep=2.0,
                                                      angleRange=20.0,
                                                      searchAstigmatism=True)

        cls.protCTFCorrection = cls._runCTFCorrection(inputSetOfTiltSeries=cls.protImportTS.outputTiltSeries,
                                                      inputSetOfCtfTomoSeries=cls.protCTFEstimation.CTFTomoSeries,
                                                      defocusTol=200,
                                                      interpolationWidth=15)

    def test_importCtfTomoSeriesOutput(self):
        self.assertSetSize(self.protImportSetOfCtfSeries.CTFs, size=2)

    def test_ctfEstimationOutputSize(self):
        self.assertSetSize(self.protCTFEstimation.CTFTomoSeries, size=2)

    def test_ctfEstimationOutputDefocusFile(self):
        for ts in self.protCTFEstimation.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            defocusFile = self.protCTFEstimation._getExtraPath(tsId, tsId + '.defocus')

            self.assertTrue(os.path.exists(defocusFile))

    def test_ctfCorrectionOutput(self):
        output = self.protCTFCorrection.TiltSeries
        self.assertSetSize(output, size=2)

        for ts in output:
            tsId = ts.getTsId()
            outputLocation = self.protCTFCorrection._getExtraPath(tsId, tsId + '.mrcs')

            self.assertTrue(os.path.exists(outputLocation))

        self.assertTrue(output.interpolated(),
                        "Tilt series does not have interpolated flag")
