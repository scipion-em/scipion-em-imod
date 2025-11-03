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
import subprocess
import traceback
import typing
from collections import Counter

from imod import Plugin
from imod.constants import OUTPUT_TILTSERIES_NAME, TLT_EXT, PATCH_TRACKING, FIDUCIAL_MODEL, \
    BRT_ENV_NAME
from imod.protocols.protocol_base import IN_TS_SET
from imod.protocols.protocol_base_ts_align import ProtImodBaseTsAlign
from imod.protocols.protocol_fiducialAlignment import TILT_ALIGN_PROGRAM
from pyworkflow.constants import BETA
from pyworkflow.protocol import PointerParam, EnumParam, IntParam, GT, STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.utils import Message, cyanStr, redStr
from tomo.objects import SetOfTiltSeries, TiltSeries

logger = logging.getLogger(__name__)


class ProtImodBRT(ProtImodBaseTsAlign, ProtStreamingBase):
    """Automatic tilt-series alignment using IMOD's batchruntomo
    (https://bio3d.colorado.edu/imod/doc/man/batchruntomo.html) wrapper made by Team Tomo
    (yet-another-imod-wrapper https://teamtomo.org/teamtomo-site-archive/).
    """

    _label = "teamtomo/batchruntomo"
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}
    stepsExecutionMode = STEPS_PARALLEL
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsReadList = []

    @classmethod
    def worksInStreaming(cls):
        return True

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        super().addInTsSetFormParam(form)
        form.addParam('alignMode', EnumParam,
                      important=True,
                      choices=['Fiducial alignment', 'Patch tracking'],
                      default=FIDUCIAL_MODEL,
                      label='Alignment mode')
        form.addParam('fidSize', IntParam,
                      condition=f'alignMode == {FIDUCIAL_MODEL}',
                      default=10,
                      validators=[GT(0)],
                      label='Fiducial diameter (nm)')
        patchTrackingCond = f'alignMode == {PATCH_TRACKING}'
        form.addParam('patchSize', IntParam,
                      condition=patchTrackingCond,
                      default=500,
                      validators=[GT(0)],
                      label='Patch side-length (A)')
        form.addParam('patchOverlapPercent', IntParam,
                      condition=patchTrackingCond,
                      default=80,
                      validators=[GT(0)],
                      label='Patch overlap percent',
                      help='Percentage of tile-length to overlap on each side.')
        form.addParallelSection(threads=3, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def stepsGeneratorStep(self) -> None:
        """
        This step should be implemented by any streaming protocol.
        It should check its input and when ready conditions are met
        call the self._insertFunctionStep method.
        """
        closeSetStepDeps = []
        inTsSet = self.getInputTsSet()
        self.readingOutput(getattr(self, OUTPUT_TILTSERIES_NAME, None))

        while True:
            listTSInput = inTsSet.getTSIds()
            if not inTsSet.isStreamOpen() and Counter(self.tsReadList) == Counter(listTSInput):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetsStep,
                                         OUTPUT_TILTSERIES_NAME,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break
            closeSetStepDeps = []
            for ts in inTsSet.iterItems():
                tsId = ts.getTsId()
                if tsId not in self.tsReadList and ts.getSize() > 0:  # Avoid processing empty TS (before the Tis are added)
                        cInputId = self._insertFunctionStep(self.convertInStep, tsId,
                                                            prerequisites=[],
                                                            needsGPU=False)
                        predFidId = self._insertFunctionStep(self.runBRT, tsId,
                                                             prerequisites=cInputId,
                                                             needsGPU=False)
                        cOutId = self._insertFunctionStep(self.createOutputStep, tsId,
                                                          prerequisites=predFidId,
                                                          needsGPU=False)
                        closeSetStepDeps.append(cOutId)
                        logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                        self.tsReadList.append(tsId)

            self.refreshStreaming(inTsSet)

    # --------------------------- STEPS functions -----------------------------
    def convertInStep(self, tsId: str):
        with self._lock:
            ts = self.getCurrentTs(tsId)
        presentAcqOrders = ts.getTsPresentAcqOrders()
        super().convertInputStep(tsId, presentAcqOrders=presentAcqOrders)

    def runBRT(self, tsId: str):
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'tsId = {tsId}: aligning...'))
                with self._lock:
                    ts = self.getCurrentTs(tsId)
                self.genTsPaths(tsId)
                args = self._getFiducialAliCmd(ts) if self.alignMode.get() == FIDUCIAL_MODEL \
                    else self._getPatchTrackingCmd(ts)
                Plugin.runBRT(self, args)
            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> {TILT_ALIGN_PROGRAM} execution failed with the exception -> {e}'))
                logger.error(traceback.format_exc())

    def createOutputStep(self, tsId: str):
        if tsId in self.failedItems:
            self.addToOutFailedSet(tsId)
            return
        try:
            with self._lock:
                ts = self.getCurrentTs(tsId)
                self.createOutTs(ts, self.getInputTsSet(pointer=True))

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    # --------------------------- INFO functions ------------------------------
    def _validate(self) -> typing.List[str]:
        # Check if the environment required by the BRT is installed (for git pulls in
        # devel mode mainly)
        errorMsg = []
        result = subprocess.run(['conda', 'env', 'list'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True)
        if BRT_ENV_NAME not in result.stdout:
            errorMsg.append('Unable to run the program batchruntomo. Please check if the '
                            'configuration variable BRT_ENV_ACTIVATION is pointing to an '
                            'existing conda environment. Otherwise, reinstall the latest '
                            'version of the plugin scipion-em-imod.')
        return errorMsg

    # --------------------------- UTILS functions -----------------------------
    def _getCommonCmd(self, ts: TiltSeries):
        tsId = ts.getTsId()
        cmd = [
            f'--tilt-series {self.getTmpOutFile(tsId)}',
            f'--tilt-angles {self.getExtraOutFile(tsId, ext=TLT_EXT)}',
            f'--output-directory {self._getExtraPath(tsId)}',
            f'--pixel-size {ts.getSamplingRate():.3f}',
            f'--nominal-rotation-angle {ts.getAcquisition().getTiltAxisAngle():.2f}',
            f'--basename {tsId}'
        ]
        return ' '.join(cmd)

    def _getFiducialAliCmd(self, ts: TiltSeries):
        cmd = [
            'fiducials',
            self._getCommonCmd(ts),
            f'--fiducial-size {self.fidSize.get()}',
        ]
        return ' '.join(cmd)

    def _getPatchTrackingCmd(self, ts: TiltSeries):
        cmd = [
            'patch-tracking',
            self._getCommonCmd(ts),
            f'--patch-size {self.patchSize.get()}',
            f'--patch-overlap-percentage {self.patchOverlapPercent.get()}'
        ]
        return ' '.join(cmd)

    def getTltFilePath(self, tsId):
        return self.getExtraOutFile(tsId, suffix="fid", ext=TLT_EXT)
