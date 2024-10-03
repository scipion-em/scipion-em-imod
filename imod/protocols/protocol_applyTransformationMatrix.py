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
from imod.protocols.protocol_base import IN_TS_SET, BINNING_FACTOR
from pwem import ALIGN_NONE
from pyworkflow.utils import Message
from tomo.objects import SetOfTiltSeries

from imod import utils
from imod.protocols import ProtImodBase
from imod.constants import (XF_EXT, ODD, EVEN, OUTPUT_TILTSERIES_NAME,
                            OUTPUT_TS_INTERPOLATED_NAME)


class ProtImodApplyTransformationMatrix(ProtImodBase):
    """
    Compute the interpolated tilt-series from its transform matrix.
    The protocol makes use of the IMod command newstack
    More info:
        https://bio3d.colorado.edu/imod/doc/man/newstack.html

    Generally, the tilt series has an associated transformation matrix
    which contains the alignment information. The transformation matrix
    is usually associated but not applied to avoid to accumulate interpolation
    errors during the image processing. This protocol allows to apply
    the transformation matrix to the tilt series
    """

    _label = 'Apply transformation'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-series to apply the transformation matrix')

        form.addParam(BINNING_FACTOR, params.IntParam,
                      default=1,
                      label='Binning for the interpolated',
                      help='Binning to be applied to the interpolated tilt-series '
                           'in IMOD convention.\nBinning is a scaling factor '
                           'given by an integer greater than 1. IMOD uses ordinary '
                           'binning (with antialiasing filter) to reduce images in '
                           'size by the given factor. The value of a binned pixel '
                           'is the average of pixel values in each block of pixels '
                           'being binned. Binning is applied before all other image '
                           'transformations.')

        form.addParam('taperInside',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=True,
                      label='Taper inwards from the edge?',
                      help='When the image is transformed areas with no information '
                           'are filled in (e.g. because of rotation). '
                           'Decide whether tapering is done inwards or outwards '
                           'from the edge.')

        form.addParam('linear',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=False,
                      label='Linear interpolation?',
                      help='From newstack man page: Use linear instead of cubic '
                           'interpolation to transform images. Linear interpolation '
                           'is more suitable when images are very noisy, but cubic '
                           'interpolation will preserve fine detail better when '
                           'noise is not an issue.')

        self.addOddEvenParams(form)

        form.addParallelSection(threads=4, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []
        for tsId in self.tsDict.keys():
            compId = self._insertFunctionStep(self.computeAlignmentStep,
                                              tsId,
                                              prerequisites=[])
            outId = self._insertFunctionStep(self.createOutputStep, tsId,
                                             prerequisites=[compId])
            closeSetStepDeps.append(outId)
        self._insertFunctionStep(self.closeOutputSetsStep,
                                 prerequisites=closeSetStepDeps)

    # --------------------------- STEPS functions ------------------------------
    def computeAlignmentStep(self, tsId):
        try:
            ts = self.tsDict[tsId]
            firstItem = ts.getFirstItem()
            self.genTsPaths(tsId)
            utils.genXfFile(ts, self.getExtraOutFile(tsId, ext=XF_EXT))

            paramsDict = self.getBasicNewstackParams(ts,
                                                     self.getExtraOutFile(tsId),
                                                     firstItem=firstItem,
                                                     xfFile=self.getExtraOutFile(tsId, ext=XF_EXT),
                                                     binning=self.binning.get(),
                                                     doSwap=True,
                                                     tsExcludedIndices=ts.getExcludedViewsIndex(),
                                                     doTaper=True)
            paramsDict["-taper"] = "1,1" if self.taperInside else "1,0"

            if self.linear:
                paramsDict["-linear"] = ""

            self.runProgram("newstack", paramsDict)

            if self.oddEvenFlag:
                paramsDict['-input'] = ts.getOddFileName()
                paramsDict['-output'] = self.getExtraOutFile(tsId, suffix=ODD)
                self.runProgram("newstack", paramsDict)

                paramsDict['-input'] = ts.getEvenFileName()
                paramsDict['-output'] = self.getExtraOutFile(tsId, suffix=EVEN)
                self.runProgram("newstack", paramsDict)

        except Exception as e:
            self._failedItems.append(tsId)
            self.error(f'Newstack execution failed for tsId {tsId} -> {e}')

    def createOutputStep(self, tsId):
        ts = self.tsDict[tsId]
        with self._lock:
            if tsId in self._failedItems:
                self.createOutputFailedSet(ts)
            else:
                outputLocation = self.getExtraOutFile(tsId)
                if os.path.exists(outputLocation):
                    output = self.getOutputSetOfTS(self.getInputSet(pointer=True),
                                                   binning=self.binning.get(),
                                                   attrName=OUTPUT_TS_INTERPOLATED_NAME,
                                                   suffix="Interpolated",)

                    self.copyTsItems(output, ts, tsId,
                                     updateTsCallback=self.updateTs,
                                     updateTiCallback=self.updateTi,
                                     copyId=False,
                                     copyTM=False,
                                     excludedViews=ts.getExcludedViewsIndex())
                else:
                    self.createOutputFailedSet(ts)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        for ts in self.getInputSet():
            if not ts.hasAlignment():
                validateMsgs.append("Some tilt-series from the input set "
                                    "are missing a transformation matrix.")
                break

        return validateMsgs

    def _summary(self):
        summary = []
        output = getattr(self, OUTPUT_TS_INTERPOLATED_NAME, None)
        if output is not None:
            summary.append(f"Input tilt-series: {self.getInputSet().getSize()}\n"
                           f"Interpolations applied: {output.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        output = getattr(self, OUTPUT_TS_INTERPOLATED_NAME, None)
        if output is not None:
            methods.append("The interpolation has been computed for "
                           f"{output.getSize()} "
                           "tilt-series using the IMOD *newstack* command.")
        return methods

    # --------------------------- UTILS functions -----------------------------
    @staticmethod
    def updateTs(tsId, ts, tsOut, **kwargs):
        tsOut.setInterpolated(True)
        tsOut.setAlignment(ALIGN_NONE)
        tsOut.getAcquisition().setTiltAxisAngle(0.)  # 0 because TS is aligned

    def updateTi(self, origIndex, index, tsId, ts, ti, tsOut, tiOut, **kwargs):
        super().updateTi(origIndex, index, tsId, ts, ti, tsOut, tiOut, **kwargs)
        tiOut.getAcquisition().setTiltAxisAngle(0.)
