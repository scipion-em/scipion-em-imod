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

import os.path

import pyworkflow.protocol.params as params
from imod.protocols.protocol_base import IN_TS_SET
from pyworkflow.utils import Message
from tomo.objects import SetOfTiltSeries

from imod import utils
from imod.protocols import ProtImodBase
from imod.constants import (ODD, EVEN, SCIPION_IMPORT, FIXED_DOSE,
                            OUTPUT_TILTSERIES_NAME)


class ProtImodDoseFilter(ProtImodBase):
    """
    Tilt-series dose filtering based on the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/mtffilter.html

    A specialized filter can be applied to perform dose weight-filtering of
    cryoEM images, particularly ones from tilt series.  The filter is as
    described in Grant and Grigorieff, 2015 (DOI: 10.7554/eLife.06980) and
    the implementation follows that in their "unblur" program.  At any fre-
    quency, the filter follows an exponential decay with dose, where the
    exponential is of the dose divided by 2 times a "critical dose" for
    that frequency.  This critical dose was empirically found to be approx-
    imated by a * k^b + c, where k is frequency; the values of a, b, c in
    that paper are used by default.
    """

    _label = 'Dose filter'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)

        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt Series',
                      help='This input tilt-series will be low pass '
                           'filtered according to their accumulated dose.')

        form.addParam('initialDose',
                      params.FloatParam,
                      default=0.0,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Initial dose (e/Å^2)',
                      help='Dose applied before any of the images in the '
                           'input file were taken; this value will be '
                           'added to all the dose values.')

        # TODO: add more options for inputting the dose information
        form.addParam('inputDoseType',
                      params.EnumParam,
                      choices=['Scipion import', 'Fixed dose'],
                      default=SCIPION_IMPORT,
                      label='Input dose source',
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='Where to find the dose information:\n'
                           '- Scipion import: use the dose provided '
                           'during import of the tilt-series\n'
                           '- Fixed dose: manually input fixed dose '
                           'for each image of the input file, '
                           'in electrons/Å^2.')

        form.addParam('fixedImageDose',
                      params.FloatParam,
                      default=FIXED_DOSE,
                      label='Fixed dose (e/Å^2)',
                      condition='inputDoseType == %i' % FIXED_DOSE,
                      help='Fixed dose for each image of the input file, '
                           'in electrons/square Ångstrom.')

        self.addOddEvenParams(form)
        form.addParallelSection(threads=4, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []
        for tsId in self.tsDict.keys():
            compId = self._insertFunctionStep(self.doseFilterStep,
                                              tsId,
                                              prerequisites=[])
            outId = self._insertFunctionStep(self.createOutputStep, tsId,
                                             prerequisites=[compId])
            closeSetStepDeps.append(outId)

        self._insertFunctionStep(self.closeOutputSetsStep,
                                 prerequisites=closeSetStepDeps)

    # --------------------------- STEPS functions -----------------------------
    def doseFilterStep(self, tsId):
        """Apply the dose filter to every tilt series"""
        try:
            ts = self.tsDict[tsId]
            firstItem = ts.getFirstItem()
            self.genTsPaths(tsId)

            params = {
                '-input': firstItem.getFileName(),
                '-output': self.getExtraOutFile(tsId),
                '-PixelSize': ts.getSamplingRate(),
                '-Voltage': int(ts.getAcquisition().getVoltage()),
            }

            if self.initialDose.get() != 0.0:
                params["-InitialDose"] = self.initialDose.get()

            if self.inputDoseType.get() == SCIPION_IMPORT:
                outputDoseFilePath = self.getExtraOutFile(tsId, ext="dose")
                utils.generateDoseFile(ts, outputDoseFilePath)
                params["-TypeOfDoseFile"] = 2
                params["-DoseWeightingFile"] = outputDoseFilePath

            elif self.inputDoseType.get() == FIXED_DOSE:
                params["-FixedImageDose"] = self.fixedImageDose.get()

            self.runProgram("mtffilter", params)

            if self.oddEvenFlag:
                params['-input'] = ts.getOddFileName()
                params['-output'] = self.getExtraOutFile(tsId, suffix=ODD)
                self.runProgram("mtffilter", params)

                params['-input'] = ts.getEvenFileName()
                params['-output'] = self.getExtraOutFile(tsId, suffix=EVEN)
                self.runProgram("mtffilter", params)

        except Exception as e:
            self._failedItems.append(tsId)
            self.error(f'Mtffilter execution failed for tsId {tsId} -> {e}')

    def createOutputStep(self, tsId):
        """Generate output filtered tilt series"""
        ts = self.tsDict[tsId]
        with self._lock:
            if tsId in self._failedItems:
                self.createOutputFailedSet(ts)
            else:
                outputLocation = self.getExtraOutFile(tsId)
                if os.path.exists(outputLocation):
                    output = self.getOutputSetOfTS(self.getInputSet(pointer=True))

                    self.copyTsItems(output, ts, tsId,
                                     updateTsCallback=self.updateTs,
                                     updateTiCallback=self.updateTi,
                                     copyDisabledViews=True,
                                     copyId=True,
                                     copyTM=True)
                else:
                    self.createOutputFailedSet(ts)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        if self.inputDoseType.get() == SCIPION_IMPORT:
            for ts in self.getInputSet():
                if ts.getFirstItem().getAcquisition().getDosePerFrame() is None:
                    validateMsgs.append(f"{ts.getTsId()} has no dose information stored "
                                        "in Scipion Metadata. To solve this, re-import "
                                        "tilt-series using the mdoc option.")
                    break

        return validateMsgs

    def _summary(self):
        summary = []

        if self.TiltSeries:
            summary.append(f"Input tilt-series: {self.getInputSet().getSize()}\n"
                           "Dose weighting applied: "
                           f"{self.TiltSeries.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")

        return summary

    def _methods(self):
        methods = []
        if self.TiltSeries:
            methods.append("The dose-weighting has been applied to "
                           f"{self.TiltSeries.getSize()} "
                           "tilt-series using the IMOD *mtffilter* command.")
        return methods

    # --------------------------- UTILS functions -----------------------------
    def updateTi(self, origIndex, index, tsId, ts, ti, tsOut, tiOut, **kwargs):
        super().updateTi(origIndex, index, tsId, ts, ti, tsOut, tiOut, **kwargs)
        # output is dose-weighted
        acq = ti.getAcquisition()
        acq.setDoseInitial(0.)
        acq.setAccumDose(0.)
        tiOut.setAcquisition(acq)

    @staticmethod
    def updateTs(tsId, ts, tsOut, **kwargs):
        acq = tsOut.getAcquisition()
        acq.setAccumDose(0.)
        acq.setDoseInitial(0.)
        tsOut.setAcquisition(acq)
