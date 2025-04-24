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
from os import stat
from os.path import exists

import numpy as np
from imod.constants import XF_EXT, OUTPUT_TS_INTERPOLATED_NAME, TLT_EXT
from imod.protocols import ProtImodBase
from imod.protocols.protocol_base import BINNING_FACTOR
from imod.utils import formatAngleList, formatTransformationMatrix
from pwem import ALIGN_NONE
from pwem.objects import Transform
from pyworkflow.protocol import params
from pyworkflow.utils import cyanStr

logger = logging.getLogger(__name__)


class ProtImodBaseTsAlign(ProtImodBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # To be filled in the child classes
        self.inTsSetPointer = None
        self.isStreamified = None
        self.isSemiStreamified = None

    # -------------------------- DEFINE param functions -----------------------
    @staticmethod
    def _insertInterpTsParams(form):
        form.addParam('computeAlignment',
                      params.BooleanParam,
                      default=False,
                      label='Generate interpolated tilt-series?',
                      important=True,
                      help='Generate and save the interpolated tilt-series applying the obtained transformation '
                           'matrices.\n'
                           'By default, the output of this protocol will be a tilseries that will have associated'
                           'the alignment information as a transformation matrix. When this option is set as Yes, '
                           'then a second output, called interpolated tilt series, is generated. The interpolated tilt '
                           'series should be used for visualization purpose but not for image processing')

        form.addParam(BINNING_FACTOR,
                      params.IntParam,
                      default=1,
                      condition='computeAlignment',
                      label='Binning for the interpolated',
                      help='Binning to be applied to the interpolated  tilt-series in IMOD '
                           'convention. \n'
                           'Binning is an scaling factor given by an integer greater than 1. '
                           'IMOD uses ordinary binning to reduce images in size by the given factor. '
                           'The value of a binned pixel is the average of pixel values in each block '
                           'of pixels being binned. Binning is applied before all other image '
                           'transformations.')

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self, tsId: str, **kwargs):
        with self._lock:
            ts = self.getCurrentItem(tsId)
            presentAcqOrders = self.getPresentAcqOrders(ts, onlyEnabled=True)  # Re-stack excluding views before reconstructing
            super().convertInputStep(tsId,
                                     oddEven=self.oddEvenFlag,
                                     presentAcqOrders=presentAcqOrders)

    def createOutTs(self, tsId, isSemiStreamified, isStreamified):
        ts = self.getCurrentItem(tsId)
        if tsId not in self.failedItems:
            xfFile = self.getExtraOutFile(tsId, suffix="fid", ext=XF_EXT)
            if exists(xfFile) and stat(xfFile).st_size != 0:
                tltFile = self.getTltFilePath(tsId)
                aliMatrix = formatTransformationMatrix(xfFile)
                tiltAngles = formatAngleList(tltFile)
                output = self.getOutputSetOfTS(self.inTsSetPointer)
                self.copyTsItems(output, ts, tsId,
                                 updateTsCallback=self.updateTsNonInterp,
                                 updateTiCallback=self.updateTiNonInterp,
                                 copyDisabledViews=True,
                                 copyId=True,
                                 copyTM=False,
                                 alignmentMatrix=aliMatrix,
                                 tltList=tiltAngles,
                                 isStreamified=isStreamified,
                                 isSemiStreamified=isSemiStreamified)
            else:
                self.createOutputFailedSet(ts)

    def createOutInterpTs(self, tsId, isSemiStreamified, isStreamified):
        if self.computeAlignment and tsId not in self.failedItems:
            tltFilePath = self.getTltFilePath(tsId)
            if exists(tltFilePath) and stat(tltFilePath).st_size != 0:
                binning = self.binning.get()
                output = self.getOutputSetOfTS(self.inTsSetPointer, binning,
                                               attrName=OUTPUT_TS_INTERPOLATED_NAME,
                                               suffix="Interpolated")
                ts = self.getCurrentItem(tsId)
                tsExcludedIndices = ts.getExcludedViewsIndex()
                tltList = formatAngleList(tltFilePath)
                self.copyTsItems(output, ts, tsId,
                                 updateTsCallback=self.updateTsInterp,
                                 updateTiCallback=self.updateTiInterp,
                                 copyId=True,
                                 copyTM=False,
                                 excludedViews=len(tsExcludedIndices) > 0,
                                 tltList=tltList,
                                 isStreamified=isStreamified,
                                 isSemiStreamified=isSemiStreamified)

    def computeInterpTsStep(self, tsId):
        """ Generate interpolated stack. """
        if self.computeAlignment and tsId not in self.failedItems:
            logger.info(cyanStr(f"tsId = {tsId}: calculating the interpolated tilt-series"))
            binning = self.binning.get()
            if tsId not in self.failedItems:
                tmpFileName = self.getExtraOutFile(tsId, suffix="fid", ext=XF_EXT)
                if exists(tmpFileName) and stat(tmpFileName).st_size != 0:
                    ts = self.getCurrentItem(tsId)
                    firstItem = ts.getFirstItem()
                    paramsDict = self.getBasicNewstackParams(
                        ts,
                        self.getExtraOutFile(tsId),
                        inputTsFileName=self.getTmpOutFile(tsId),
                        xfFile=tmpFileName,
                        firstItem=firstItem,
                        binning=binning,
                        doSwap=True,
                        doTaper=True,
                    )
                    self.runProgram('newstack', paramsDict)

    # --------------------------- UTILS functions -----------------------------
    def getTltFilePath(self, tsId):
        # To be defined by the child classes
        pass

    @staticmethod
    def updateTsNonInterp(tsId, ts, tsOut, **kwargs):
        tsOut.setAlignment2D()

    @staticmethod
    def updateTiNonInterp(origIndex, index, tsId, ts, ti, tsOut, tiOut, alignmentMatrix=None, tltList=None, **kwargs):
        transform = Transform()
        identityMatrix = np.eye(3)
        if ti.isEnabled():
            newTransform = alignmentMatrix[:, :, index]
            newTransformArray = np.array(newTransform)
            tiltAngle = float(tltList[index])  # It may have been refined
        else:
            tiltAngle = ti.getTiltAngle()
            if ti.hasTransform():
                newTransformArray = ti.getTransform().getMatrix()
            else:
                newTransformArray = identityMatrix

        if ti.hasTransform():
            previousTransform = ti.getTransform().getMatrix()
            previousTransformArray = np.array(previousTransform)
            outputTransformMatrix = np.matmul(newTransformArray, previousTransformArray)
            transform.setMatrix(outputTransformMatrix)
        else:
            transform.setMatrix(newTransformArray)

        tiOut.setTransform(transform)
        tiOut.setTiltAngle(tiltAngle)

    @staticmethod
    def updateTsInterp(tsId, ts, tsOut, **kwargs):
        tsOut.getAcquisition().setTiltAxisAngle(0.)
        tsOut.setAlignment(ALIGN_NONE)
        tsOut.setInterpolated(True)

    def updateTiInterp(self, origIndex, index, tsId, ts, ti, tsOut, tiOut, tltList=None, **kwargs):
        super().updateTi(origIndex, index, tsId, ts, ti, tsOut, tiOut, **kwargs)
        tiOut.setTiltAngle(float(tltList[index]))
