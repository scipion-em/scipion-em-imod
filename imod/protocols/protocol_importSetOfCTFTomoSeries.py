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

import os

import pyworkflow.protocol.params as params
from pyworkflow import BETA
from pyworkflow.object import Set

from .protocol_base import ProtImodBase, OUTPUT_CTF_SERIE


class ProtImodImportSetOfCtfTomoSeries(ProtImodBase):
    """
    Protocol to import estimations of CTF series from tilt-series into Scipion.
    """

    _label = 'Import tomo CTFs'
    _devStatus = BETA

    defocusUTolerance = 20
    defocusVTolerance = 20

    def _defineParams(self, form):
        self._defineImportParams(form)

        form.addParam('exclusionWords', params.StringParam,
                      label='Exclusion words:',
                      help="List of words separated by a space that the "
                           "path should not have",
                      expertLevel=params.LEVEL_ADVANCED)

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Input tilt-series',
                      help='Select the tilt-series to which the imported '
                           'estimation of CTF will be paired. The file '
                           'names of the file and the defocus file must '
                           'be the same (except the extension).')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):

        self._insertFunctionStep(self.importSetOfCtfTomoSeries)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def importSetOfCtfTomoSeries(self):

        inputSetOfTiltSeries = self.inputSetOfTiltSeries.get()

        output = self.getOutputSetOfCTFTomoSeries(OUTPUT_CTF_SERIE)

        for ts in inputSetOfTiltSeries:
            tsFileName = ts.getFirstItem().parseFileName(extension='')

            for ctfFile, _ in self.iterFiles():
                defocusFilePath = ctfFile
                defocusFileName = os.path.basename(os.path.splitext(ctfFile)[0])

                if tsFileName == defocusFileName:

                    self.info("Parsing file: " + defocusFilePath)

                    self.addCTFTomoSeriesToSetFromDefocusFile(ts, defocusFilePath,
                                                              output)

    def closeOutputSetsStep(self):
        self.CTFTomoSeries.setStreamState(Set.STREAM_CLOSED)
        self.CTFTomoSeries.write()
        self._store()

    # --------------------------- UTILS functions -----------------------------
    def _getSetOfTiltSeries(self, pointer=False):
        return self.inputSetOfTiltSeries.get() if not pointer else self.inputSetOfTiltSeries

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.CTFTomoSeries:
            summary.append("Imported CTF tomo series: %d"
                           % (self.CTFTomoSeries.getSize()))

        return summary

    def _validate(self):
        errorMsg = []
        if not self.getMatchFiles():
            errorMsg.append('Unable to find the files provided:\n\n'
                            '\t-filePath = %s\n'
                            '\t-pattern = %s\n' % (self.filesPath.get(),
                                                   self.filesPattern.get()))

        return errorMsg

