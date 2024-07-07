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
from pyworkflow.protocol.constants import STEPS_PARALLEL
import pyworkflow.utils.path as path
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries

from imod import utils
from imod.protocols.protocol_base import ProtImodBase
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

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
        self.excludedViewsFromFile = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        form.addParam('excludeViewsFile',
                      params.FileParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Exclude views file',
                      help='File containing the views to be excluded for each '
                           'tilt-series belonging to the set.\n\n'
                           'The format of the text file must be two columns, '
                           'the first one being the tilt series ID of '
                           'the series from which the views will be excluded '
                           'and the second the views to exclude, numbered from '
                           '1. The syntax for this exclude list is a comma '
                           'separated list of ranges with no spaces between '
                           'them (e.g., 1,4-5,60-70). \n\n'
                           'An example of this file comes as follows:\n'
                           'TS_01 1,4-6,8,44-47\n'
                           'TS_02 3,10-12,24\n'
                           '...')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        if self.excludeViewsFile.get():
            self.excludedViewsFromFile = utils.readExcludeViewsFile(self.excludeViewsFile.get())

        closeSetStepDeps = []
        for tsId in self.tsDict.keys():
            exclStepId = self._insertFunctionStep(self.excludeViewsStep,
                                                  tsId,
                                                  prerequisites=[])
            outStepId = self._insertFunctionStep(self.generateOutputStackStep,
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
        firstItem = ts.getFirstItem()
        self.genTsPaths(tsId)

        outputFileName = self.getExtraOutFile(tsId)
        excludedViews = self.getExcludedViews(ts)

        if excludedViews:
            path.copyFile(firstItem.getFileName(), outputFileName)
            params = {
                '-StackName': outputFileName,
                '-ViewsToExclude': ",".join(map(str, excludedViews)),
            }
            self.runProgram("excludeviews", params)
        else:
            # Just create the link
            self.info(f"No views to exclude for {tsId}")
            path.createLink(firstItem.getFileName(), outputFileName)

    def generateOutputStackStep(self, tsId):
        ts = self.tsDict[tsId]
        output = self.getOutputSetOfTS(self.getInputSet())

        newTs = TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        output.append(newTs)

        excludedViews = self.getExcludedViews(ts)

        i = 1
        for index, tiltImage in enumerate(ts):
            if (index + 1) not in excludedViews:
                newTi = TiltImage()
                newTi.copyInfo(tiltImage, copyId=False, copyTM=True)
                newTi.setAcquisition(tiltImage.getAcquisition())
                newTi.setLocation(i, self.getExtraOutFile(tsId))
                newTs.append(newTi)
                i += 1

        newTs.write(properties=False)
        output.update(newTs)
        output.write()
        self._store(output)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append("Excluded views:\n")

            for tsIn, tsOut in zip(self.getInputSet(), self.TiltSeries):
                summary.append(f"Tilt-series: {tsIn.getTsId()}; "
                               f"Size: {tsIn.getSize()} ---> {tsOut.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")

        return summary

    # --------------------------- UTILS functions -----------------------------
    def getExcludedViews(self, ts):
        """ Returns the indexes of the tilt to exclude for a
        specific tilt series"""
        if self.excludeViewsFile.get():
            matrix = self.excludedViewsFromFile
            return matrix.get(ts.getTsId(), [])
        else:
            return ts._getExcludedViewsIndex()
