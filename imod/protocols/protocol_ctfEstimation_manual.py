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

from pyworkflow import BETA
from tomo.objects import SetOfCTFTomoSeries

from .protocol_ctfEstimation_automatic import ProtImodAutomaticCtfEstimation
from .protocol_base import OUTPUT_CTF_SERIE


class ProtImodManualCtfEstimation(ProtImodAutomaticCtfEstimation):
    """
    CTF estimation of a set of input tilt-series using the IMOD procedure.
    Runs the protocol through the interactive GUI. The resulting defocus values
    MUST BE SAVED manually by the user.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/ctfplotter.html

    """

    _label = 'Manual CTF estimation'
    _devStatus = BETA
    _interactiveMode = True

    def __init__(self, **args):
        ProtImodAutomaticCtfEstimation.__init__(self, **args)
        self.OUTPUT_PREFIX = OUTPUT_CTF_SERIE

    def _insertAllSteps(self):
        self.inputTiltSeries = None
        self._insertFunctionStep(self.runCTFEtimationStep, interactive=True)

    # --------------------------- STEPS functions ----------------------------

    def runCTFEtimationStep(self):
        from imod.viewers import ImodGenericViewer
        self.inputSetOfTiltSeries = self._getSetOfTiltSeries()
        view = ImodGenericViewer(None, self, self.inputSetOfTiltSeries,
                                 createSetButton=True,
                                 isInteractive=True,
                                 itemDoubleClick=True)
        view.show()
        self.createOutput()

    def runAllSteps(self, obj):
        objId = obj.getObjId()
        # DO not use this to avoid interpolation. CTF will be done on the input TS (hopefully non interpolated).
        self.convertInputStep(objId)
        expDefoci = self.getExpectedDefocus()
        self.ctfEstimation(objId, expDefoci)

    def createOutput(self):
        suffix = self._getOutputSuffix(SetOfCTFTomoSeries)
        outputSetName = self.OUTPUT_PREFIX + str(suffix)
        setOfTiltseries = self._getSetOfTiltSeries()
        for item in setOfTiltseries.iterItems(iterate=False):
            self.createOutputStep(item.getObjId(), outputSetName)
        self.closeOutputSetsStep()

    def _summary(self):
        summary = []
        return summary
