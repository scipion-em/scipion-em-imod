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
import typing

from imod import Plugin
from imod.constants import OUTPUT_TILTSERIES_NAME, TLT_EXT, XF_EXT
from imod.protocols import ProtImodFiducialAlignment
from imod.protocols.protocol_base import IN_TS_SET
from imod.utils import formatTransformationMatrix, formatAngleList
from pyworkflow.constants import BETA
from pyworkflow.object import Pointer
from pyworkflow.protocol import PointerParam, STEPS_PARALLEL, EnumParam, IntParam, GT
from pyworkflow.utils import Message, cyanStr
from tomo.objects import SetOfTiltSeries, TiltSeries

logger = logging.getLogger(__name__)

# Alignment modes
FIDUCIAL_ALIGNMENT = 0
PATCH_TRACKING = 1


class ProtImodBRT(ProtImodFiducialAlignment):
    """Automatic tilt-series alignment using IMOD's batchruntomo
    (https://bio3d.colorado.edu/imod/doc/man/batchruntomo.html) wrapper made by Team Tomo
    (https://teamtomo.org/teamtomo-site-archive/).
    """

    _label = "Tilt-series alignment (Team Tomo's wrapper)"
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}
    _devStatus = BETA
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._failedItems = []

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        form.addParam(IN_TS_SET, PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-series')
        form.addParam('alignMode', EnumParam,
                      important=True,
                      choices=['Fiducial alignment', 'Patch tracking'],
                      default=FIDUCIAL_ALIGNMENT,
                      label='Alignment mode')
        form.addParam('fidSize', IntParam,
                      condition=f'alignMode == {FIDUCIAL_ALIGNMENT}',
                      default=10,
                      validators=[GT(0)],
                      label='Fiducial diameter (nm)')
        patchTrackingCond = f'alignMode == {PATCH_TRACKING}'
        form.addParam('patchSize', IntParam,
                      condition=patchTrackingCond,
                      default=500,
                      validators=[GT(0)],
                      label='Patch sidelength (A)')
        form.addParam('patchOverlapPercent', IntParam,
                      condition=patchTrackingCond,
                      default=33,
                      validators=[GT(0)],
                      label='Patch overlap percent',
                      help='Percentage of tile-length to overlap on each side.')
        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        closeSetStepDeps = []
        for ts in self.getInputSet().iterItems():
            tsId = ts.getTsId()
            cInputId = self._insertFunctionStep(self.convertInputStep, tsId,
                                                prerequisites=[],
                                                needsGPU=False)
            predFidId = self._insertFunctionStep(self.runBRT, tsId,
                                                 prerequisites=cInputId,
                                                 needsGPU=True)
            cOutId = self._insertFunctionStep(self.createOutputStep, tsId,
                                              prerequisites=predFidId,
                                              needsGPU=False)
            closeSetStepDeps.append(cOutId)
        self._insertFunctionStep(self._closeOutputSet,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def runBRT(self, tsId: str):
        logger.info(cyanStr(f'===> tsId = {tsId}: aligning...'))
        try:
            ts = self.getCurrentItem(tsId)
            self.genTsPaths(tsId)
            args = self._getFiducialAliCmd(ts) if self.alignMode.get() == FIDUCIAL_ALIGNMENT \
                else self._getPatchTrackingCmd(ts)
            Plugin.runBRT(self, args)
        except Exception as e:
            self._failedItems.append(tsId)
            logger.error(f'tiltalign execution failed for tsId {tsId} -> {e}')

    def createOutputStep(self, tsId: str):
        with self._lock:
            ts = self.getCurrentItem(tsId)
            if tsId not in self._failedItems:
                xfFile = self.getExtraOutFile(tsId, suffix="fid", ext=XF_EXT)
                tltFile = self.getExtraOutFile(tsId, suffix="fid", ext=TLT_EXT)
                aliMatrix = formatTransformationMatrix(xfFile)
                tiltAngles = formatAngleList(tltFile)
                output = self.getOutputSetOfTS(self.getInputSet(pointer=True))
                self.copyTsItems(output, ts, tsId,
                                 updateTsCallback=self.updateTsNonInterp,
                                 updateTiCallback=self.updateTiNonInterp,
                                 copyDisabledViews=True,
                                 copyId=True,
                                 copyTM=False,
                                 alignmentMatrix=aliMatrix,
                                 tltList=tiltAngles)
            else:
                self.createOutputFailedSet(ts)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        pass

    def _methods(self):
        pass

    # --------------------------- UTILS functions -----------------------------
    def getInputSet(self, pointer: bool = False) -> typing.Union[Pointer, SetOfTiltSeries]:
        tsPointer = getattr(self, IN_TS_SET)
        return tsPointer if pointer else tsPointer.get()

    def getCurrentItem(self, tsId: str) -> TiltSeries:
        return self.getInputSet().getItem(TiltSeries.TS_ID_FIELD, tsId)

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