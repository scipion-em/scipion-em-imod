# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
# *
# * [1] National Center of Biotechnology, CSIC, Spain
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
from os import stat
from os.path import exists
from typing import Tuple, List
import numpy as np
from imod.constants import XF_EXT
from imod.convert.convert import readXfFile
from imod.protocols import ProtImodBase
from imod.utils import formatAngleList
from pwem.objects import Transform
from pyworkflow.object import Pointer
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import TiltSeries, TiltImage

logger = logging.getLogger(__name__)


IDENTITY_MATRIX = np.eye(3)  # Store in memory instead of multiple creation


class ProtImodBaseTsAlign(ProtImodBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- STEPS functions -----------------------------

    # --------------------------- UTILS functions -----------------------------
    def getTltFilePath(self, tsId: str):
        # To be defined by the child classes
        pass

    @retry_on_sqlite_lock(log=logger)
    def createOutTs(self,
                    ts: TiltSeries,
                    inTsSetPointer: Pointer) -> None:
        tsId = ts.getTsId()
        xfFile = self.getExtraOutFile(tsId, suffix="fid", ext=XF_EXT)
        if exists(xfFile) and stat(xfFile).st_size != 0:
            tltFile = self.getTltFilePath(tsId)
            aliMatrix = readXfFile(xfFile)
            tiltAngles = formatAngleList(tltFile)
            # Set of tilt-series
            with self._lock:
                outTsSet = self.getOutputSetOfTS(inTsSetPointer)
                # Tilt-series
                outTs = TiltSeries()
                outTs.copyInfo(ts)
                outTs.setAlignment2D()
                outTsSet.append(outTs)
                # Tilt-images
                stackIndex = 0
                for ti in ts.iterItems(orderBy=TiltImage.INDEX_FIELD):
                    outTi = TiltImage()
                    outTi.copyInfo(ti)
                    if ti.isEnabled():
                        tiltAngle, newTransformArray = self._getTrDataEnabled(stackIndex,
                                                                              aliMatrix,
                                                                              tiltAngles)
                        stackIndex += 1
                    else:
                        tiltAngle, newTransformArray = self._getTrDataDisabled(ti)
                    self._updateTiltImage(ti, outTi, newTransformArray, tiltAngle)
                    self.setTsOddEven(tsId, outTi, binGenerated=False)
                    outTs.append(outTi)
                # Data persistence
                outTs.write()
                outTsSet.update(outTs)
                outTsSet.write()
                self._store(outTsSet)
                # Close explicitly the outputs (for streaming)
                self.closeOutputsForStreaming()
        else:
            logger.error(f'tsId = {tsId} -> Output file {xfFile} was not generated or is empty. Skipping... ')

    @staticmethod
    def _getTrDataEnabled(stackIndex: int,
                          alignmentMatrix: np.ndarray,
                          tiltAngleList: List[float]) -> Tuple[float, np.ndarray]:
        newTransform = alignmentMatrix[:, :, stackIndex]
        newTransformArray = np.array(newTransform)
        tiltAngle = float(tiltAngleList[stackIndex])
        return tiltAngle, newTransformArray

    @staticmethod
    def _getTrDataDisabled(ti: TiltImage) -> Tuple[float, np.ndarray]:
        tiltAngle = ti.getTiltAngle()
        if ti.hasTransform():
            newTransformArray = ti.getTransform().getMatrix()
        else:
            newTransformArray = IDENTITY_MATRIX
        return tiltAngle, newTransformArray

    @staticmethod
    def _updateTiltImage(ti: TiltImage,
                         outTi: TiltImage,
                         newTransformArray: np.ndarray,
                         tiltAngle: float) -> None:
        transform = Transform()
        if ti.hasTransform():
            previousTransform = ti.getTransform().getMatrix()
            previousTransformArray = np.array(previousTransform)
            outputTransformMatrix = np.matmul(newTransformArray, previousTransformArray)
            transform.setMatrix(outputTransformMatrix)
        else:
            transform.setMatrix(newTransformArray)

        outTi.setTransform(transform)
        outTi.setTiltAngle(tiltAngle)
