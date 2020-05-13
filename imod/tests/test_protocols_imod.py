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
        cls.inputDataSet = DataSet.getDataSet('tomo-em')
        cls.inputTS = cls.inputDataSet.getFile('ts1')


class TestImodXcorrPrealignment(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.inputDataSet = DataSet.getDataSet('tomo-em')
        cls.inputTS = cls.inputDataSet.getFile('ts1')

    def _runImportTiltSeries(self):
        protImportTS = self.newProtocol(tomo.protocols.ProtImportTs,
                                        filesPath=self.inputTS,
                                        filesPattern=self.inputTS,
                                        voltage=300,
                                        magnification=105000,
                                        sphericalAberration=2.7,
                                        amplitudeContrast=0.1,
                                        samplingRate=20.2,
                                        doseInitial=0,
                                        dosePerFrame=0.3,
                                        minAngle=-55.0,
                                        maxAngle=65.0,
                                        stepAngle=2.0)
        self.launchProtocol(protImportTS)
        return protImportTS

    def _runXcorrPrealignment(self, inputTS, computeAlignmentToggle,
                              binning, rotationAngle):
        protXcorr = self.newProtocol(ProtImodXcorrPrealignment,
                                     inputSetOfTiltSeries=inputTS,
                                     computeAlignment=computeAlignmentToggle,
                                     binning=binning,
                                     rotationAngle=rotationAngle)
        self.launchProtocol(protXcorr)
        return protXcorr

    def test_outputTS(self):
        protImportTS = self._runImportTiltSeries()
        inputTS = protImportTS.outputTiltSeries
        protXcorr = self._runXcorrPrealignment(inputTS=inputTS,
                                               computeAlignmentToggle=1,
                                               binning=1,
                                               rotationAngle=-12.5)
        self.assertIsNotNone(protXcorr.outputSetOfTiltSeries)

    def test_outputInterpolatedTS(self):
        protImportTS = self._runImportTiltSeries()
        inputTS = protImportTS.outputTiltSeries
        protXcorr = self._runXcorrPrealignment(inputTS=inputTS,
                                               computeAlignmentToggle=0,
                                               binning=1,
                                               rotationAngle=-12.5)
        self.assertIsNotNone(protXcorr.outputInterpolatedSetOfTiltSeries)

    def test_outputSamplingRate(self):
        inputBinning=2
        protImportTS = self._runImportTiltSeries()
        inputTS = protImportTS.outputTiltSeries
        protXcorr = self._runXcorrPrealignment(inputTS=inputTS,
                                               computeAlignmentToggle=0,
                                               binning=inputBinning,
                                               rotationAngle=-12.5)
        inSamplingRate = protXcorr.inputSetOfTiltSeries.get().getSamplingRate()
        outSamplingRate = protXcorr.outputInterpolatedSetOfTiltSeries.getSamplingRate()
        self.assertTrue(inSamplingRate * inputBinning == outSamplingRate)
