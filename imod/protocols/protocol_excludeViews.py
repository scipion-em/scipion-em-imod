# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * National Center of Biotechnology, CSIC, Spain
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
import logging

import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from imod.protocols.protocol_base import IN_TS_SET
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.utils import Message, cyanStr, redStr
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from imod.protocols import ProtImodBase
from imod.constants import OUTPUT_TILTSERIES_NAME, ODD, EVEN, EXCLUDE_VIEWS_PROGRAM

logger = logging.getLogger(__name__)


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
    stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        super().addInTsSetFormParam(form)
        self.addOddEvenParams(form)
        form.addParallelSection(threads=2, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []
        for ts in self.getInputTsSet():
            tsId = ts.getTsId()
            exclStepId = self._insertFunctionStep(self.excludeViewsStep,
                                                  tsId,
                                                  prerequisites=[],
                                                  needsGPU=False)
            outStepId = self._insertFunctionStep(self.createOutputStep,
                                                 tsId,
                                                 prerequisites=exclStepId,
                                                 needsGPU=False)
            closeSetStepDeps.append(outStepId)
        self._insertFunctionStep(self.closeOutputSetsStep,
                                 OUTPUT_TILTSERIES_NAME,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def excludeViewsStep(self, tsId: str):
        self.genTsPaths(tsId)
        try:
            with self._lock:
                ts = self.getCurrentTs(tsId)
                firstTi = ts.getFirstItem()

            tsFileName = firstTi.getFileName()
            outputFileName = self.getExtraOutFile(tsId)

            if ts.hasExcludedViews():
                logger.info(cyanStr(f'tsId = {tsId}: excluding the disabled views...'))
                path.copyFile(tsFileName, outputFileName)
                excludeViewsInd = ts.getTsExcludedViewsIndices(ts.getTsPresentAcqOrders())
                eVParams = {
                    '-StackName': outputFileName,
                    '-ViewsToExclude': ",".join(map(str, excludeViewsInd)),
                }
                self.runProgram(EXCLUDE_VIEWS_PROGRAM, eVParams)
            else:
                # Just create the link
                logger.info(cyanStr(f"tsId = {tsId} -> No views to exclude."))
                path.createLink(tsFileName, outputFileName)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {EXCLUDE_VIEWS_PROGRAM} execution failed with the exception -> {e}'))

    def createOutputStep(self, tsId):
        if tsId in self.failedItems:
            self.addToOutFailedSet(tsId)
            return

        try:
            with self._lock:
                ts = self.getCurrentTs(tsId)
                outFn = self.getExtraOutFile(tsId)
                outTsSet = self.getOutputSetOfTS(self.getInputTsSet(pointer=True))
                outTs = TiltSeries()
                outTs.copyInfo(ts)
                outTsSet.append(outTs)

                angleMin = 999
                angleMax = -999
                accumDose = 0
                initialDose = 999
                tiList = []
                for ti in ts.iterItems():
                    if ti.isEnabled():
                        # Update the acquisition of the TS. The accumDose, angle min and angle max for the re-stacked TS, as
                        # these values may change if the removed tilt-images are the first or the last, for example.
                        tiAngle = ti.getTiltAngle()
                        angleMin = min(tiAngle, angleMin)
                        angleMax = max(tiAngle, angleMax)
                        accumDose = max(ti.getAcquisition().getAccumDose(), accumDose)
                        initialDose = min(ti.getAcquisition().getDoseInitial(), initialDose)

                        outTi = TiltImage()
                        outTi.copyInfo(ti)
                        outTi.setFileName(outFn)
                        if self.doOddEven:
                            outTi.setOddEven([self.getExtraOutFile(tsId, suffix=ODD),
                                              self.getExtraOutFile(tsId, suffix=EVEN)])
                        else:
                            outTi.setOddEven([])  # the input may have odd/even but the user may have decided not
                            # to consider them in the current execution, so they should be set to empty to avoid
                            # next protocols be confused about having them.
                        tiList.append(outTi)

                if ts.hasExcludedViews():
                    # Update the acquisition minAngle and maxAngle values of the tilt-series
                    acq = outTs.getAcquisition()
                    acq.setAngleMin(angleMin)
                    acq.setAngleMax(angleMax)
                    acq.setAccumDose(accumDose)
                    acq.setDoseInitial(initialDose)
                    outTs.setAcquisition(acq)
                    # Update the acquisition minAngle and maxAngle values of each tilt-image acq while preserving their
                    # specific accum and initial dose values
                    for tiOut in tiList:
                        tiAcq = tiOut.getAcquisition()
                        tiAcq.setAngleMin(angleMin)
                        tiAcq.setAngleMax(angleMax)
                        tiOut.setAcquisition(tiAcq)
                        outTs.append(tiOut)
                    outTs.setAnglesCount(len(outTs))
                else:
                    for tiOut in tiList:
                        outTs.append(tiOut)

                outTs.write()
                outTsSet.update(outTs)
                outTsSet.write()
                self._store(outTsSet)

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            summary.append("Excluded views:\n")

            tsInSet = self.getInputTsSet().iterItems(orderBy='_tsId')
            tsOutSet = output.iterItems(orderBy='_tsId')
            for tsIn, tsOut in zip(tsInSet, tsOutSet):
                summary.append(f"Tilt-series: {tsIn.getTsId()}; "
                               f"Size: {tsIn.getSize()} ---> {tsOut.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")

        return summary

    # --------------------------- UTILS functions -----------------------------

