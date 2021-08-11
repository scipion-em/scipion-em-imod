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

from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from tomo.protocols.protocol_base import ProtTomoImportFiles
from pyworkflow.object import Set


class ProtImodBase(ProtTomoImportFiles, EMProtocol, ProtTomoBase):
    """
    Base class with methods used in the rest of the imod protocols
    """

    def __init__(self, **args):
        ProtTomoImportFiles.__init__(self, **args)

    def _defineImportParams(self, form):
        """ Method to define import params in protocol form """
        ProtTomoImportFiles._defineImportParams(self, form)

    # --------------------------- OUTPUT functions ----------------------------
    def getOutputSetOfTiltSeries(self):
        """ Method to generate output classes of set of tilt-series"""

        if hasattr(self, "outputSetOfTiltSeries"):
            self.outputSetOfTiltSeries.enableAppend()

        else:
            outputSetOfTiltSeries = self._createSetOfTiltSeries()
            outputSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfTiltSeries=outputSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfTiltSeries)

        return self.outputSetOfTiltSeries

    def getOutputInterpolatedSetOfTiltSeries(self):
        """ Method to generate output interpolated classes of set of tilt-series"""

        if hasattr(self, "outputInterpolatedSetOfTiltSeries"):
            self.outputInterpolatedSetOfTiltSeries.enableAppend()

        else:
            outputInterpolatedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Interpolated')
            outputInterpolatedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputInterpolatedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())

            if self.binning > 1:
                samplingRate = self.inputSetOfTiltSeries.get().getSamplingRate()
                samplingRate *= self.binning.get()
                outputInterpolatedSetOfTiltSeries.setSamplingRate(samplingRate)

            outputInterpolatedSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputInterpolatedSetOfTiltSeries=outputInterpolatedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputInterpolatedSetOfTiltSeries)

        return self.outputInterpolatedSetOfTiltSeries

    # --------------------------- UTILS functions ----------------------------
    def iterFiles(self):
        """ Iterate through the files matched with the pattern.
        Provide the fileName and fileId.
        """
        filePaths = self.getMatchFiles()

        filePaths = self._excludeByWords(filePaths)

        for fileName in filePaths:
            if self._idRegex:
                # Try to match the file id from filename
                # this is set by the user by using #### format in the pattern
                match = self._idRegex.match(fileName)
                if match is None:
                    raise Exception("File '%s' doesn't match the pattern '%s'"
                                    % (fileName, self.getPattern()))

                fileId = int(match.group(1))

            else:
                fileId = None

            yield fileName, fileId

    def _excludeByWords(self, files):
        exclusionWords = self.exclusionWords.get()

        if exclusionWords is None:
            return files

        exclusionWordList = exclusionWords.split()

        allowedFiles = []

        for file in files:
            if any(bannedWord in file for bannedWord in exclusionWordList):
                print("%s excluded. Contains any of %s" % (file, exclusionWords))
                continue
            allowedFiles.append(file)

        return allowedFiles
