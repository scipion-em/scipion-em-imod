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
from imod.constants import OUTPUT_TILTSERIES_NAME
from imod.protocols import ProtImodXraysEraser
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
tsDims40 = [3710, 3838, 40]
tsDims41 = [3710, 3838, 41]
tomoDimsThk300 = [928, 928, 300]
tomoDimsThk280 = [928, 928, 280]
tomoDimsThk340 = [928, 928, 340]


class TestImodBase(TestBaseCentralizedLayer):
    importedTs = None
    unbinnedSRate = DataSetRe4STATuto.unbinnedPixSize.value
    expectedTsSetSize = 2
    expectedDimensions = {
        TS_03: tsDims40,
        TS_54: tsDims41,
    }
    testAcqObjDict = {
        TS_03: DataSetRe4STATuto.testAcq03.value,
        TS_54: DataSetRe4STATuto.testAcq54.value,
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
        pass

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


class TestXRayEraser(TestImodBase):

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTs = cls._runImportTs()

    def testXRayEraser(self):
        tsXRayErased = self._runXRayEraser(self.importedTs)
        self.checkTiltSeries(tsXRayErased,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate,
                             imported=True,
                             expectedDimensions=self.expectedDimensions,
                             testAcqObj=self.testAcqObjDict,
                             anglesCount=self.anglesCountDict,
                             isHetereogeneousSet=True)

