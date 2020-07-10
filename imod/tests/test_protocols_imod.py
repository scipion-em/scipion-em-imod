# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.tests import *
from imod.protocols import *
import tomo


class TestImodBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @classmethod
    def _runImportTiltSeries(cls, filesPath, pattern, voltage, magnification, sphericalAberration, amplitudeContrast,
                             samplingRate, doseInitial, dosePerFrame, minAngle, maxAngle, stepAngle):
        cls.protImportTS = cls.newProtocol(tomo.protocols.ProtImportTs,
                                           filesPath=filesPath,
                                           filesPattern=pattern,
                                           voltage=voltage,
                                           magnification=magnification,
                                           sphericalAberration=sphericalAberration,
                                           amplitudeContrast=amplitudeContrast,
                                           samplingRate=samplingRate,
                                           doseInitial=doseInitial,
                                           dosePerFrame=dosePerFrame,
                                           minAngle=minAngle,
                                           maxAngle=maxAngle,
                                           stepAngle=stepAngle)
        cls.launchProtocol(cls.protImportTS)
        return cls.protImportTS

    @classmethod
    def _runXcorrPrealignment(cls, inputSoTS, computeAlignmentToggle, binning, rotationAngle):
        cls.protXcorr = cls.newProtocol(ProtImodXcorrPrealignment,
                                        inputSetOfTiltSeries=inputSoTS,
                                        computeAlignment=computeAlignmentToggle,
                                        binning=binning,
                                        rotationAngle=rotationAngle)
        cls.launchProtocol(cls.protXcorr)
        return cls.protXcorr

    @classmethod
    def _runTSNormalization(cls, inputSoTS, binning, floatDensities, modeToOutput, scaleRangeToggle, scaleRangeMax,
                            scaleRangeMin, meanSdToggle, scaleMean, scaleSd, scaleMax, scaleMin):
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


class TestImodXcorrPrealignment(TestImodBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.inputDataSet = DataSet.getDataSet('tomo-em')
        cls.inputSoTS = cls.inputDataSet.getFile('ts1')
        cls.binning = 2

        cls.protImportTS = cls._runImportTiltSeries(filesPath=os.path.split(cls.inputSoTS)[0],
                                                    pattern="BB{TS}.st",
                                                    voltage=300,
                                                    magnification=105000,
                                                    sphericalAberration=2.7,
                                                    amplitudeContrast=0.1,
                                                    samplingRate=20.2,
                                                    doseInitial=0,
                                                    dosePerFrame=0.3,
                                                    minAngle=-55,
                                                    maxAngle=65.0,
                                                    stepAngle=2.0)

        cls.protXcorr = cls._runXcorrPrealignment(inputSoTS=cls.protImportTS.outputTiltSeries,
                                                  computeAlignmentToggle=0,
                                                  binning=cls.binning,
                                                  rotationAngle=-12.5)

    def test_outputTS(self):
        self.assertIsNotNone(self.protXcorr.outputSetOfTiltSeries)

    def test_outputInterpolatedTS(self):
        self.assertIsNotNone(self.protXcorr.outputInterpolatedSetOfTiltSeries)

    def test_outputSamplingRate(self):
        inSamplingRate = self.protXcorr.inputSetOfTiltSeries.get().getSamplingRate()
        outSamplingRate = self.protXcorr.outputInterpolatedSetOfTiltSeries.getSamplingRate()
        self.assertTrue(inSamplingRate * self.binning == outSamplingRate)


class TestTiltSeriesNormalization(TestImodBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.protImportTS = TestImodXcorrPrealignment.protImportTS

        cls.protNormalizeTiltSeries = cls._runTSNormalization(inputSoTS=cls.protImportTS.outputTiltSeries,
                                                              binning=2,
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

    def test_outputNormalizedTS(self):
        pass
