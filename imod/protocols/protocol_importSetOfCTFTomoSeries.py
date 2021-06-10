# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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

import os
from pwem.protocols import EMProtocol
from pyworkflow import BETA
import pyworkflow.object as pwobj
import pyworkflow.protocol.params as params
from pyworkflow import BETA
from pyworkflow.object import Set
from tomo.protocols import ProtTomoBase
from tomo.protocols.protocol_base import ProtTomoImportFiles
import tomo.objects as tomoObj
from imod import utils


class ProtImodImportSetOfCtfTomoSeries(ProtTomoImportFiles, EMProtocol, ProtTomoBase):
    """
    Protocol to import estimations of CTF series from tilt-series into Scipion.
    """

    _label = 'import tomo CTFs'
    _devStatus = BETA

    def __init__(self, **args):
        ProtTomoImportFiles.__init__(self, **args)

    def _defineParams(self, form):
        ProtTomoImportFiles._defineImportParams(self, form)

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Input tilt-series',
                      help='Select the tilt-series to which the imported estimation of CTF will be paired. The file '
                           'names of the file and the defocus file must be the same (except the extension).')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importSetOfCtfTomoSeries')
        self._insertFunctionStep('closeOutputSetsStep')

    # --------------------------- STEPS functions ----------------------------
    def importSetOfCtfTomoSeries(self):

        inputSetOfTiltSeries = self.inputSetOfTiltSeries.get()

        self.getOutputSetOfCTFTomoSeries()

        for ts in inputSetOfTiltSeries:
            tsId = ts.getTsId()

            tsFileName = ts.getFirstItem().parseFileName(extension='')

            for ctfFile, _ in self.iterFiles():
                defocusFilePath = ctfFile
                defocusFileName = os.path.basename(os.path.splitext(ctfFile)[0])

                if tsFileName == defocusFileName:

                    print("Parsing file: " + defocusFilePath)

                    defocusFileFlag = utils.getDefocusFileFlag(defocusFilePath)

                    newCTFTomoSeries = tomoObj.CTFTomoSeries()
                    newCTFTomoSeries.copyInfo(ts)
                    newCTFTomoSeries.setTiltSeries(ts)
                    newCTFTomoSeries.setTsId(tsId)
                    newCTFTomoSeries.setIMODDefocusFileFlag(defocusFileFlag)

                    # We need to create now all the attributes of this object in order to append it to the set and be
                    # able to update it posteriorly.

                    newCTFTomoSeries.setNumberOfEstimationsInRange(None)
                    self.outputSetOfCTFTomoSeries.append(newCTFTomoSeries)

                    if defocusFileFlag == 0:
                        " Plain estimation "
                        defocusUDict = utils.readCTFEstimationInfoFile(defocusFilePath,
                                                                       flag=defocusFileFlag)

                    elif defocusFileFlag == 1:
                        " Astigmatism estimation "
                        defocusUDict, defocusVDict, defocusAngleDict = utils.readCTFEstimationInfoFile(defocusFilePath,
                                                                                                       flag=defocusFileFlag)

                    elif defocusFileFlag == 4:
                        " Phase-shift information "
                        defocusUDict, phaseShiftDict = utils.readCTFEstimationInfoFile(defocusFilePath,
                                                                                       flag=defocusFileFlag)

                    elif defocusFileFlag == 5:
                        " Astigmatism and phase shift estimation "
                        defocusUDict, defocusVDict, defocusAngleDict, phaseShiftDict = \
                            utils.readCTFEstimationInfoFile(defocusFilePath,
                                                            flag=defocusFileFlag)

                    elif defocusFileFlag == 37:
                        " Astigmatism, phase shift and cut-on frequency estimation "
                        defocusUDict, defocusVDict, defocusAngleDict, phaseShiftDict, cutOnFreqDict = \
                            utils.readCTFEstimationInfoFile(defocusFilePath,
                                                            flag=defocusFileFlag)

                    else:
                        raise Exception(
                            "Defocus file flag do not supported. Only supported formats corresponding to flags 0, "
                            "1, 4, 5, and 37.")

                    for index, _ in enumerate(ts):
                        newCTFTomo = tomoObj.CTFTomo()
                        newCTFTomo.setIndex(pwobj.Integer(index + 1))

                        if defocusFileFlag == 0:
                            " Plain estimation "
                            newCTFTomo._defocusUList = pwobj.CsvList(pType=float)
                            newCTFTomo.setDefocusUList(defocusUDict[index + 1])

                        elif defocusFileFlag == 1:
                            " Astigmatism estimation "
                            newCTFTomo._defocusUList = pwobj.CsvList(pType=float)
                            newCTFTomo.setDefocusUList(defocusUDict[index + 1])

                            newCTFTomo._defocusVList = pwobj.CsvList(pType=float)
                            newCTFTomo.setDefocusVList(defocusVDict[index + 1])

                            newCTFTomo._defocusAngleList = pwobj.CsvList(pType=float)
                            newCTFTomo.setDefocusAngleList(defocusAngleDict[index + 1])

                        elif defocusFileFlag == 4:
                            " Phase-shift information "
                            newCTFTomo._defocusUList = pwobj.CsvList(pType=float)
                            newCTFTomo.setDefocusUList(defocusUDict[index + 1])

                            newCTFTomo._phaseShiftList = pwobj.CsvList(pType=float)
                            newCTFTomo.setPhaseShiftList(phaseShiftDict[index + 1])

                        elif defocusFileFlag == 5:
                            " Astigmatism and phase shift estimation "
                            newCTFTomo._defocusUList = pwobj.CsvList(pType=float)
                            newCTFTomo.setDefocusUList(defocusUDict[index + 1])

                            newCTFTomo._defocusVList = pwobj.CsvList(pType=float)
                            newCTFTomo.setDefocusVList(defocusVDict[index + 1])

                            newCTFTomo._defocusAngleList = pwobj.CsvList(pType=float)
                            newCTFTomo.setDefocusAngleList(defocusAngleDict[index + 1])

                            newCTFTomo._phaseShiftList = pwobj.CsvList(pType=float)
                            newCTFTomo.setPhaseShiftList(phaseShiftDict[index + 1])

                        elif defocusFileFlag == 37:
                            " Astigmatism, phase shift and cut-on frequency estimation "
                            newCTFTomo._defocusUList = pwobj.CsvList(pType=float)
                            newCTFTomo.setDefocusUList(defocusUDict[index + 1])

                            newCTFTomo._defocusVList = pwobj.CsvList(pType=float)
                            newCTFTomo.setDefocusVList(defocusVDict[index + 1])

                            newCTFTomo._defocusAngleList = pwobj.CsvList(pType=float)
                            newCTFTomo.setDefocusAngleList(defocusAngleDict[index + 1])

                            newCTFTomo._phaseShiftList = pwobj.CsvList(pType=float)
                            newCTFTomo.setPhaseShiftList(phaseShiftDict[index + 1])

                            newCTFTomo._cutOnFreqList = pwobj.CsvList(pType=float)
                            newCTFTomo.setCutOnFreqList(cutOnFreqDict[index + 1])

                            defocusUDict, defocusVDict, defocusAngleDict, phaseShiftDict, cutOnFreqDict = \
                                utils.readCTFEstimationInfoFile(defocusFilePath,
                                                                flag=defocusFileFlag)

                        newCTFTomo.completeInfoFromList()

                        newCTFTomoSeries.append(newCTFTomo)

                    newCTFTomoSeries.setNumberOfEstimationsInRangeFromDefocusList()

                    newCTFTomoSeries.write(properties=False)

                    self.outputSetOfCTFTomoSeries.update(newCTFTomoSeries)
                    self.outputSetOfCTFTomoSeries.write()

                    self._store()

    def closeOutputSetsStep(self):
        self.outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfCTFTomoSeries(self):
        if hasattr(self, "outputSetOfCTFTomoSeries"):
            self.outputSetOfCTFTomoSeries.enableAppend()
        else:
            outputSetOfCTFTomoSeries = tomoObj.SetOfCTFTomoSeries.create(self._getPath(),
                                                                         template='CTFmodels%s.sqlite')

            outputSetOfCTFTomoSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(self.inputSetOfTiltSeries.get())
            outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfCTFTomoSeries=outputSetOfCTFTomoSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfCTFTomoSeries)
        return self.outputSetOfCTFTomoSeries

