# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
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

import os
from pyworkflow import BETA
import pyworkflow.object as pwobj
import pyworkflow.protocol.params as params
from pyworkflow import BETA
from pyworkflow.object import Set
import tomo.objects as tomoObj
from imod import utils
from imod.protocols.protocol_base import ProtImodBase


class ProtImodImportSetOfCtfTomoSeries(ProtImodBase):
    """
    Protocol to import estimations of CTF series from tilt-series into Scipion.
    """

    _label = 'import tomo CTFs'
    _devStatus = BETA

    defocusUTolerance = 20
    defocusVTolerance = 20

    def _defineParams(self, form):
        self._defineImportParams(form)

        form.addParam('exclusionWords', params.StringParam,
                      label='Exclusion words:',
                      help="List of words separated by a space that the path should not have",
                      expertLevel=params.LEVEL_ADVANCED)

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Input tilt-series',
                      help='Select the tilt-series to which the imported estimation of CTF will be paired. The file '
                           'names of the file and the defocus file must be the same (except the extension).')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importSetOfCtfTomoSeries')
        self._insertFunctionStep('closeOutputSetsStep')

    # --------------------------- STEPS functions ----------------------------
    def importSetOfCtfTomoSeries(self):

        inputSetOfTiltSeries = self.inputSetOfTiltSeries.get()

        self.getOutputSetOfCTFTomoSeries()

        for ts in inputSetOfTiltSeries:
            tsId = ts.getTsId()

            tsFileName = ts.getFirstItem().parseFileName(extension='')

            for ctfFile, _ in self.iterFiles():
                defocusFilePath = ctfFile
                defocusFileName = os.path.basename(os.path.splitext(ctfFile)[0])

                if tsFileName == defocusFileName:

                    print("Parsing file: " + defocusFilePath)

                    self.addCTFTomoSeriesToSetFromDefocusFile(ts, defocusFilePath)

    def closeOutputSetsStep(self):
        self.outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputSetOfCTFTomoSeries'):
            summary.append("Imported CTF tomo series: %d.\n"
                           % (self.outputSetOfCTFTomoSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputSetOfCTFTomoSeries'):
            methods.append("%d CTF tomo series have been imported into the scipion framework.\n"
                           % (self.outputSetOfCTFTomoSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods

