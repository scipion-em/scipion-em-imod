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
import logging
import os
import time
from os.path import exists

import pyworkflow.protocol.params as params
from imod.convert import genXfFile
from imod.protocols.protocol_base import IN_TS_SET, NEWSTACK_PROGRAM
from imod.protocols.protocol_base_preprocess import ProtImodBasePreprocess
from pyworkflow.protocol import STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.utils import Message, cyanStr, redStr
from tomo.objects import SetOfTiltSeries, TiltImage, TiltSeries
from imod.constants import OUTPUT_TILTSERIES_NAME, ODD, EVEN, NO_TS_PROCESSED_MSG, OUTPUT_TS_FAILED_NAME, XF_EXT

logger = logging.getLogger(__name__)


class ProtImodTsNormalization(ProtImodBasePreprocess, ProtStreamingBase):
    """
    Normalize input tilt-series and change its storing formatting.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/newstack.html

    IMOD tilt series preprocess makes use of the Newstack command.
    In particular, three functionalities are possible:\n

    _1 Binning_: The protocol also allows to bin tilt series. This
    means to reduce the dimensions of the tilt series keeping but
    keeping most of the information. The binning factor or simply
    binning is an integer number and represent the scaling factor
    of the images. Binning 2 means that the original images will
    be twice the binned ones.
    _2 Normalization_: This protocol allows to scale the gray values
    of the images, also called normalization, to a common range or
    mean of density. The most used normalization consists in zero
    mean and standard deviation one.\n

    _3 storage format_: IMOD is able to modify the number of bit of
    the stored data in order to reduce the disc occupancy.

    """

    _label = 'Tilt-series preprocess'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsReadList = []

    @classmethod
    def worksInStreaming(cls):
        return True

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form, *args):
        form.addSection(Message.LABEL_INPUT)
        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')
        super()._defineParams(form)

    # -------------------------- INSERT steps functions -----------------------
    def stepsGeneratorStep(self) -> None:
        """
        This step should be implemented by any streaming protocol.
        It should check its input and when ready conditions are met
        call the self._insertFunctionStep method.
        """
        self._initialize()
        binning = self.binning.get()
        inTsSet = self.getInputSet()
        self.readingOutput(inTsSet)
        closeSetStepDeps = []

        while True:
            listTSInput = inTsSet.getTSIds()
            if not inTsSet.isStreamOpen() and self.tsReadList == listTSInput:
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetsStep,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break
            closeSetStepDeps = []
            for ts in inTsSet.iterItems():
                tsId = ts.getTsId()
                if tsId not in self.tsReadList and ts.getSize() > 0:  # Avoid processing empty TS (before the Tis are added)
                    convId = self._insertFunctionStep(self.convertInputStep,
                                                      tsId,
                                                      prerequisites=[],
                                                      needsGPU=False)
                    compId = self._insertFunctionStep(self.generateOutputStackStep,
                                                      tsId,
                                                      binning,
                                                      prerequisites=convId,
                                                      needsGPU=False)
                    outId = self._insertFunctionStep(self.createOutputStep,
                                                     tsId,
                                                     binning,
                                                     prerequisites=compId,
                                                     needsGPU=False)
                    closeSetStepDeps.append(outId)
                    logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                    self.tsReadList.append(tsId)
            time.sleep(10)
            if inTsSet.isStreamOpen():
                with self._lock:
                    inTsSet.loadAllProperties()  # refresh status for the streaming

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self, tsId: str):
        self.genTsPaths(tsId)
        with self._lock:
            ts = self.getCurrentItem(tsId)
            firstTi = ts.getFirstItem()
        # Generate the xf alignment file
        if firstTi.hasTransform():
            xfFile = self.getExtraOutFile(ts.getTsId(), ext=XF_EXT)
            genXfFile(ts,
                      xfFile,
                      ignoreExcludedViews=True)  # The newstack program can exclude them itself if necessary

    def generateOutputStackStep(self, tsId: str, binning: int):
        logger.info(cyanStr(f'===> tsId = {tsId}: preprocessing...'))
        try:
            with self._lock:
                ts = self.getCurrentItem(tsId)
            xfFile = self.getExtraOutFile(ts.getTsId(), ext=XF_EXT)
            xfFile = xfFile if exists(xfFile) else None
            norm = self.floatDensities.get()
            paramsDict = self.getBasicNewstackParams(ts,
                                                     ts.getFirstItem().getFileName(),
                                                     self.getExtraOutFile(tsId),
                                                     xfFile=xfFile,
                                                     binning=binning,
                                                     doNorm=norm != 0)

            paramsDict["-antialias"] = self.antialias.get() + 1
            # Float densities
            if norm > 0:
                paramsDict["-FloatDensities"] = norm
                if norm == 2:
                    paramsDict["-MeanAndStandardDeviation"] = f"{self.scaleMean.get()},{self.scaleSd.get()}"
                elif norm == 4:
                    paramsDict["-ScaleMinAndMax"] = f"{self.scaleMax.get()},{self.scaleMin.get()}"

            if self.getModeToOutput() is not None:
                paramsDict["-ModeToOutput"] = self.getModeToOutput()

            self.runProgram(NEWSTACK_PROGRAM, paramsDict)

            if self.doOddEven:
                paramsDict['-input'] = self.getTmpOutFile(tsId, suffix=ODD)
                paramsDict['-output'] = self.getExtraOutFile(tsId, suffix=ODD)
                self.runProgram(NEWSTACK_PROGRAM, paramsDict)

                paramsDict['-input'] = self.getTmpOutFile(tsId, suffix=EVEN)
                paramsDict['-output'] = self.getExtraOutFile(tsId, suffix=EVEN)
                self.runProgram(NEWSTACK_PROGRAM, paramsDict)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(f'tsId = {tsId} -> {NEWSTACK_PROGRAM} execution '
                         f'failed with the exception -> {e}')

    def createOutputStep(self, tsId: str, binning: int):
        if tsId in self.failedItems:
            self.addToOutFailedSet(tsId)
        else:
            try:
                outputFn = self.getExtraOutFile(tsId)
                if exists(outputFn):
                    with self._lock:
                        ts = self.getCurrentItem(tsId)
                        outTsSet = self.getOutputSetOfTS(self.getInputSet(pointer=True), binning)
                        outTs = TiltSeries()
                        outTs.copyInfo(ts)
                        outTsSet.append(outTs)
                        for ti in ts.iterItems(orderBy=TiltImage.INDEX_FIELD):
                            outTi = TiltImage()
                            outTi.copyInfo(ti)
                            outTi.setFileName(outputFn)
                            self.updateTransformMatrix(outTi, binning=binning)
                            if self.doOddEven:
                                outTi.setOddEven([self.getExtraOutFile(tsId, suffix=ODD),
                                                  self.getExtraOutFile(tsId, suffix=EVEN)])
                            outTs.append(outTi)
                        outTs.write()
                        outTsSet.update(outTs)
                        outTsSet.write()
                        self._store(outTsSet)
                else:
                    logger.error(f'tsId = {tsId} -> Output file {outputFn} was not generated. Skipping... ')
            except Exception as e:
                logger.error(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... ')

    def closeOutputSetsStep(self):
        if not getattr(self, OUTPUT_TILTSERIES_NAME, None):
            raise Exception('No tilt-series was generated. Please '
                            'check the Output Log > run.stdout and run.stderr')
        self._closeOutputSet()

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            summary.append(f"Input tilt-series: {self.getInputSet().getSize()}\n"
                           f"Interpolations applied: {output.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            methods.append(f"{output.getSize()} tilt-series have been "
                           "normalized using the IMOD *newstack* command.")
        return methods

    # --------------------------- UTILS functions -----------------------------
    def getModeToOutput(self):
        parseParamsOutputMode = {
            0: None,
            1: 101,
            2: 0,
            3: 1,
            4: 6,
            5: 2
        }
        return parseParamsOutputMode[self.modeToOutput.get()]

    @staticmethod
    def updateTransformMatrix(ti: TiltImage, binning: int = 1) -> None:
        if ti.hasTransform() and binning != 1:
            transform = ti.getTransform()
            matrix = transform.getMatrix()

            matrix[0][2] /= binning
            matrix[1][2] /= binning

            transform.setMatrix(matrix)
            ti.setTransform(transform)


