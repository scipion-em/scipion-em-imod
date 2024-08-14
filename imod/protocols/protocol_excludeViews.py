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
from imod.protocols.protocol_base import IN_TS_SET
from pyworkflow.protocol.constants import STEPS_PARALLEL
import pyworkflow.utils.path as path
from pyworkflow.utils import Message
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries
from imod.protocols import ProtImodBase
from imod.constants import OUTPUT_TILTSERIES_NAME, EVEN, ODD


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
        # if self.excludeViewsFile.get():
        #     self.excludedViewsFromFile = utils.readExcludeViewsFile(self.excludeViewsFile.get())

        closeSetStepDeps = []
        for tsId in self.tsDict.keys():
            exclStepId = self._insertFunctionStep(self.excludeViewsStep, tsId, prerequisites=[])
            outStepId = self._insertFunctionStep(self.generateOutputStackStep, tsId, prerequisites=[exclStepId])
            closeSetStepDeps.append(outStepId)
        self._insertFunctionStep(self.closeOutputSetsStep, prerequisites=closeSetStepDeps)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self.getInputSet()}
        self.oddEvenFlag = self.applyToOddEven(self.getInputSet())

    def excludeViewsStep(self, tsId):
        ts = self.tsDict[tsId]
        tsFileName = ts.getFirstItem().getFileName()
        self.genTsPaths(tsId)

        outputFileName = self.getExtraOutFile(tsId)
        excludedViews = self.getExcludedViews(ts)

        if excludedViews:
            tsFileNames = [tsFileName]
            outputFileNames = [outputFileName]
            if self.oddEvenFlag:
                tsFileNames.extend([ts.getEven(), ts.getOdd()])
                outputFileNames.extend([self.getExtraOutFile(tsId, suffix=EVEN),
                                        self.getExtraOutFile(tsId, suffix=ODD)])

            for inFile, outFile in zip(tsFileNames, outputFileNames):
                path.copyFile(inFile, outFile)
                paramsDict = {
                    '-StackName': outFile,
                    '-ViewsToExclude': ",".join(map(str, excludedViews)),
                }
                self.runProgram("excludeviews", paramsDict)
        else:
            # Just create the link
            self.info(f"No views to exclude for {tsId}")
            path.createLink(tsFileName, outputFileName)

    def generateOutputStackStep(self, tsId):
        ts = self.tsDict[tsId]
        output = self.getOutputSetOfTS(self.getInputSet())
        excludedViews = self.getExcludedViews(ts)
        if excludedViews:
            newTs = TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            oddEvenFNames = [self.getExtraOutFile(tsId, suffix=EVEN),
                             self.getExtraOutFile(tsId, suffix=ODD)]
            output.append(newTs)

            angleList = []
            doseList = []
            counter = 1
            for ti in ts:
                if ti.getIndex() not in excludedViews:
                    newTi = TiltImage()
                    newTi.copyInfo(ti, copyId=False, copyTM=True)
                    newTi.setAcquisition(ti.getAcquisition())
                    newTi.setLocation(counter, self.getExtraOutFile(tsId))
                    angleList.append(ti.getTiltAngle())
                    doseList.append(ti.getAcquisition().getAccumDose())
                    if self.oddEvenFlag:
                        newTi.setOddEven(oddEvenFNames)
                    newTs.append(newTi)
                    counter += 1
            # Update the acquisition of the TS. The accumDose, angle min and angle max for the re-stacked TS, as
            # these values may change if the removed tilt-images are the first or the last, for example.
            acq = newTs.getAcquisition()
            acq.setAngleMin(min(angleList))
            acq.setAngleMax(max(angleList))
            acq.setAccumDose(max(doseList))
            newTs.setAcquisition(acq)
            newTs.setAnglesCount(len(newTs))
        else:
            newTs = ts.clone(ignoreAttrs=[])
            output.append(newTs)
            for ti in ts:
                newTi = ti.clone()
                newTs.append(newTi)

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
    @staticmethod
    def getExcludedViews(ts):
        """ Returns the indexes of the tilt to exclude for a
        specific tilt series"""
        # if self.excludeViewsFile.get():
        #     matrix = self.excludedViewsFromFile
        #     return matrix.get(ts.getTsId(), [])
        # else:
        return ts.getExcludedViewsIndex()

