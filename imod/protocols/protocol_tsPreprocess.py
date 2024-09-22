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
from imod.protocols.protocol_base import IN_TS_SET
from imod.protocols.protocol_base_preprocess import ProtImodBasePreprocess
from pyworkflow.utils import Message
from tomo.objects import SetOfTiltSeries
from imod.constants import OUTPUT_TILTSERIES_NAME, ODD, EVEN


class ProtImodTsNormalization(ProtImodBasePreprocess):
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

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')
        super()._defineParams(form)
        form.addParallelSection(threads=4, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        binning = self.binning.get()
        closeSetStepDeps = []

        for tsId in self.tsDict.keys():
            convId = self._insertFunctionStep(self.convertInputStep,
                                              tsId,
                                              prerequisites=[])
            compId = self._insertFunctionStep(self.generateOutputStackStep,
                                              tsId,
                                              binning,
                                              prerequisites=[convId])
            outId = self._insertFunctionStep(self.createOutputStep,
                                             tsId,
                                             binning,
                                             prerequisites=[compId])
            closeSetStepDeps.append(outId)

        self._insertFunctionStep(self.closeOutputSetsStep,
                                 prerequisites=closeSetStepDeps)

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self, tsId, **kwargs):
        super().convertInputStep(tsId,
                                 imodInterpolation=None,
                                 generateAngleFile=False,
                                 oddEven=self.oddEvenFlag)

    def generateOutputStackStep(self, tsId, binning):
        try:
            ts = self.tsDict[tsId]
            firstItem = ts.getFirstItem()
            xfFile = None
            norm = self.floatDensities.get()
            paramsDict = self.getBasicNewstackParams(ts,
                                                     self.getExtraOutFile(tsId),
                                                     inputTsFileName=self.getTmpOutFile(tsId),
                                                     xfFile=xfFile,
                                                     firstItem=firstItem,
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

            self.runProgram("newstack", paramsDict)

            if self.oddEvenFlag:
                paramsDict['-input'] = self.getTmpOutFile(tsId, suffix=ODD)
                paramsDict['-output'] = self.getExtraOutFile(tsId, suffix=ODD)
                self.runProgram("newstack", paramsDict)

                paramsDict['-input'] = self.getTmpOutFile(tsId, suffix=EVEN)
                paramsDict['-output'] = self.getExtraOutFile(tsId, suffix=EVEN)
                self.runProgram("newstack", paramsDict)

        except Exception as e:
            self._failedItems.append(tsId)
            self.error(f'newstack execution failed for tsId {tsId} -> {e}')

    def createOutputStep(self, tsId, binning):
        ts = self.tsDict[tsId]
        with self._lock:
            if tsId in self._failedItems:
                self.createOutputFailedSet(ts)
            else:
                outputFn = self.getExtraOutFile(tsId)
                if os.path.exists(outputFn):
                    output = self.getOutputSetOfTS(self.getInputSet(pointer=True), binning)
                    self.copyTsItems(output, ts, tsId,
                                     updateTiCallback=self.updateTi,
                                     copyDisabledViews=True,
                                     copyId=True,
                                     copyTM=True,
                                     binning=binning)
                else:
                    self.createOutputFailedSet(ts)

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
    def updateTM(newTi, binning=1):
        if binning != 1:
            transform = newTi.getTransform()
            matrix = transform.getMatrix()

            matrix[0][2] /= binning
            matrix[1][2] /= binning

            transform.setMatrix(matrix)
            newTi.setTransform(transform)

    def updateTi(self, origIndex, index, tsId, ts, ti, tsOut, tiOut, **kwargs):
        super().updateTi(origIndex, index, tsId, ts, ti, tsOut, tiOut, **kwargs)
        # Transformation matrix
        if ti.hasTransform():
            self.updateTM(tiOut, binning=self.binning.get())

