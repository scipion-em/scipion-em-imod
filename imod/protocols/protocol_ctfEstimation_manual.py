# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os

from pyworkflow import BETA
from tomo.objects import SetOfCTFTomoSeries
from .protocol_ctfEstimation_automatic import ProtImodAutomaticCtfEstimation


class ProtImodManualCtfEstimation(ProtImodAutomaticCtfEstimation):
    """
    CTF estimation of a set of input tilt-series using the IMOD procedure.
    Run the protocol through the interactive GUI. Defocus values are saved
    to file and exit after autofitting. The program will not ask for confirmation before
    removing existing entries in the defocus table. If run in interactive mode defocus values
    MUST BE SAVED manually by the user.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/ctfplotter.html

    """

    _label = 'manual CTF estimation (step 2)'
    _devStatus = BETA
    _interactiveMode = True

    def __init__(self, **args):
        ProtImodAutomaticCtfEstimation.__init__(self, **args)

    def _insertAllSteps(self):
        self.inputTiltSeries = None
        self._insertFunctionStep(self.runCTFEtimationStep, interactive=True)

    # --------------------------- STEPS functions ----------------------------

    def runCTFEtimationStep(self):
        from imod.viewers import ImodGenericViewer
        self.inputSetOfTiltSeries = self._getSetOfTiltSeries()
        view = ImodGenericViewer(None, self, self.inputSetOfTiltSeries,
                                 displayAllButton=False, createSetButton=True,
                                 isInteractive=True,
                                 itemDoubleClick=True)
        view.show()

    def runAllSteps(self, obj):
        objId = obj.getObjId()
        self.convertInputStep(objId)
        self.ctfEstimation(objId)

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

    # ------------------ UTILS METHODS ------------------------------------

    def getFilePath(self, ts, suffix="", extension=""):
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        return os.path.join(extraPrefix, ts.getFirstItem().parseFileName(suffix=suffix,
                                                                         extension=extension))
