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

import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from imod.protocols.protocol_base import IN_TS_SET
from pyworkflow.utils import Message
from tomo.objects import SetOfTiltSeries
from imod.protocols import ProtImodBase
from imod.constants import OUTPUT_TILTSERIES_NAME


class ProtImodExcludeViews(ProtImodBase):
    """
    excludeviews - Reversibly remove views from a tilt series stack

    By default, the protocol will remove disabled tilt images from the input TS.
    Alternatively, you can provide a text file with a list of tilts to exclude.

    If you use this protocol, make sure this output tilt series is use for everything else
    CTF estimation, per particle per tilt, tomogram reconstruction....
    More info:
        https://bio3d.colorado.edu/imod/doc/man/excludeviews.html
    """

    _label = 'Exclude views'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        self.addOddEvenParams(form)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []
        for tsId in self.tsDict.keys():
            exclStepId = self._insertFunctionStep(self.excludeViewsStep,
                                                  tsId,
                                                  prerequisites=[])
            outStepId = self._insertFunctionStep(self.createOutputStep,
                                                 tsId,
                                                 prerequisites=[exclStepId])
            closeSetStepDeps.append(outStepId)
        self._insertFunctionStep(self.closeOutputSetsStep,
                                 prerequisites=closeSetStepDeps)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[])
                       for ts in self.getInputSet()}

    def excludeViewsStep(self, tsId):
        ts = self.tsDict[tsId]
        tsFileName = ts.getFirstItem().getFileName()
        self.genTsPaths(tsId)

        outputFileName = self.getExtraOutFile(tsId)
        excludedViews = self.getExcludedViews(ts)

        if excludedViews:
            path.copyFile(tsFileName, outputFileName)
            params = {
                '-StackName': outputFileName,
                '-ViewsToExclude': ",".join(map(str, excludedViews)),
            }
            self.runProgram("excludeviews", params)
        else:
            # Just create the link
            self.info(f"No views to exclude for {tsId}")
            path.createLink(tsFileName, outputFileName)

    def createOutputStep(self, tsId):
        ts = self.tsDict[tsId]
        excludedViews = self.getExcludedViews(ts)

        with self._lock:
            output = self.getOutputSetOfTS(self.getInputSet(pointer=True))
            if excludedViews:
                self.copyTsItems(output, ts, tsId,
                                 updateTiCallback=self.updateTi,
                                 copyId=False, copyTM=True,
                                 excludedViews=excludedViews)
            else:
                newTs = ts.clone(ignoreAttrs=[])
                output.append(newTs)
                for ti in ts:
                    newTi = ti.clone()
                    newTs.append(newTi)

                output.update(newTs)
                self._store(output)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append("Excluded views:\n")

            tsInSet = self.getInputSet().iterItems(orderBy='_tsId')
            tsOutSet = self.TiltSeries.iterItems(orderBy='_tsId')
            for tsIn, tsOut in zip(tsInSet, tsOutSet):
                summary.append(f"Tilt-series: {tsIn.getTsId()}; "
                               f"Size: {tsIn.getSize()} ---> {tsOut.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")

        return summary

    # --------------------------- UTILS functions -----------------------------
    @staticmethod
    def getExcludedViews(ts):
        """ Returns the indexes of the tilt to exclude for a
        specific tilt series starting from 1. """
        return ts.getExcludedViewsIndex()
