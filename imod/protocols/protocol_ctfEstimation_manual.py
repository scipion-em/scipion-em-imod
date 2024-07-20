# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# **************************************************************************
from pyworkflow.protocol.constants import STEPS_SERIAL

from tomo.objects import SetOfCTFTomoSeries

from imod.protocols import ProtImodAutomaticCtfEstimation
from imod.constants import OUTPUT_CTF_SERIE


class ProtImodManualCtfEstimation(ProtImodAutomaticCtfEstimation):
    """
    CTF estimation of a set of input tilt-series using the IMOD procedure.
    Runs the protocol through the interactive GUI. The resulting defocus values
    MUST BE SAVED manually by the user.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/ctfplotter.html

    This GUI program will plot the logarithm of a rotationally averaged
    power spectrum of an input tilt series after subtracting the noise
    floor.  The method is based on periodogram averaging; namely, averaging
    of spectra from small, overlapping areas of the images, referred to as
    tiles.  The user can interactively choose which projection views are
    included in the averaging.  It is also possible to run the program non-
    interactively for automatic fitting.
    """

    _label = 'CTF estimation (manual)'
    _interactiveMode = True
    _possibleOutputs = {OUTPUT_CTF_SERIE: SetOfCTFTomoSeries}

    def __init__(self, **args):
        ProtImodAutomaticCtfEstimation.__init__(self, **args)
        self.stepsExecutionMode = STEPS_SERIAL
        self.OUTPUT_PREFIX = OUTPUT_CTF_SERIE

    def _insertAllSteps(self):
        self.inputTiltSeries = None
        self._insertFunctionStep(self.runCTFEtimationStep, interactive=True)

    # --------------------------- STEPS functions ----------------------------

    def runCTFEtimationStep(self):
        from imod.viewers import ImodGenericView
        self.inputSetOfTiltSeries = self.getInputSet()
        view = ImodGenericView(None, self, self.inputSetOfTiltSeries,
                               createSetButton=True,
                               isInteractive=True)
        view.show()
        self.createOutput()

    def runAllSteps(self, obj):
        tsId = obj.getTsId()
        self.convertInputStep(tsId)
        expDefoci = self.getExpectedDefocus()
        self.ctfEstimation(tsId, expDefoci)

    def createOutput(self):
        suffix = self._getOutputSuffix(SetOfCTFTomoSeries)
        outputSetName = self.OUTPUT_PREFIX + str(suffix)
        tsIdList = self.inputTiltSeries.getTSIds()
        for tsId in tsIdList:
            self.createOutputStep(tsId, outputSetName)

        self.closeOutputSetsStep()

    def _summary(self):
        summary = []
        return summary
