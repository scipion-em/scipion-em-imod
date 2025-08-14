# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
# *
# * [1] National Center of Biotechnology, CSIC, Spain
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
# **************************************************************************
from pyworkflow.object import CsvList, Float
import pyworkflow.utils as pwutils

with pwutils.weakImport('tomo'):
    from tomo.objects import CTFTomo


class ImodCtfParser:

    def __init__(self, protocol):
        self.protocol = protocol

    def parseTSDefocusFile(self, inputTs, defocusFilePath, newCTFTomoSeries):
        """ Parse tilt-series ctf estimation file.
        :param inputTs: input tilt-series
        :param defocusFilePath: input *.defocus file to be parsed
        :param newCTFTomoSeries: output CTFTomoSeries
        """
        defocusFileFlag = self.getDefocusFileFlag(defocusFilePath)

        if defocusFileFlag == 0:
            " Plain estimation "
            defocusUDict = self.readCTFEstimationInfoFile(defocusFilePath,
                                                          flag=defocusFileFlag)

        elif defocusFileFlag == 1:
            " Astigmatism estimation "
            defocusUDict, defocusVDict, defocusAngleDict = self.readCTFEstimationInfoFile(defocusFilePath,
                                                                                          flag=defocusFileFlag)

        elif defocusFileFlag == 4:
            " Phase-shift information "
            defocusUDict, phaseShiftDict = self.readCTFEstimationInfoFile(defocusFilePath,
                                                                          flag=defocusFileFlag)

        elif defocusFileFlag == 5:
            " Astigmatism and phase shift estimation "
            defocusUDict, defocusVDict, defocusAngleDict, phaseShiftDict = \
                self.readCTFEstimationInfoFile(defocusFilePath,
                                               flag=defocusFileFlag)

        elif defocusFileFlag == 37:
            " Astigmatism, phase shift and cut-on frequency estimation "
            defocusUDict, defocusVDict, defocusAngleDict, phaseShiftDict, cutOnFreqDict = \
                self.readCTFEstimationInfoFile(defocusFilePath,
                                               flag=defocusFileFlag)

        else:
            raise ValueError(
                f"Defocus file flag {defocusFileFlag} is not supported. Only supported formats "
                "correspond to flags 0, 1, 4, 5, and 37.")

        for i, ti in enumerate(inputTs):
            tiObjId = ti.getObjId()
            newCTFTomo = CTFTomo()
            " Plain estimation (any defocus flag)"
            newCTFTomo._defocusUList = CsvList(pType=float)
            newCTFTomo.setDefocusUList(defocusUDict.get(tiObjId, [0.]))

            if ti.isEnabled():

                if defocusFileFlag == 1:
                    " Astigmatism estimation "
                    newCTFTomo._defocusVList = CsvList(pType=float)
                    newCTFTomo.setDefocusVList(defocusVDict.get(tiObjId, [0.]))

                    newCTFTomo._defocusAngleList = CsvList(pType=float)
                    newCTFTomo.setDefocusAngleList(defocusAngleDict.get(tiObjId, [0.]))

                elif defocusFileFlag == 4:
                    " Phase-shift information "
                    newCTFTomo._phaseShiftList = CsvList(pType=float)
                    newCTFTomo.setPhaseShiftList(phaseShiftDict.get(tiObjId, [0.]))

                elif defocusFileFlag == 5:
                    " Astigmatism and phase shift estimation "
                    newCTFTomo._defocusVList = CsvList(pType=float)
                    newCTFTomo.setDefocusVList(defocusVDict.get(tiObjId, [0.]))

                    newCTFTomo._defocusAngleList = CsvList(pType=float)
                    newCTFTomo.setDefocusAngleList(defocusAngleDict.get(tiObjId, [0.]))

                    newCTFTomo._phaseShiftList = CsvList(pType=float)
                    newCTFTomo.setPhaseShiftList(phaseShiftDict.get(tiObjId, [0.]))

                elif defocusFileFlag == 37:
                    " Astigmatism, phase shift and cut-on frequency estimation "
                    newCTFTomo._defocusVList = CsvList(pType=float)
                    newCTFTomo.setDefocusVList(defocusVDict.get(tiObjId, [0.]))

                    newCTFTomo._defocusAngleList = CsvList(pType=float)
                    newCTFTomo.setDefocusAngleList(defocusAngleDict.get(tiObjId, [0.]))

                    newCTFTomo._phaseShiftList = CsvList(pType=float)
                    newCTFTomo.setPhaseShiftList(phaseShiftDict.get(tiObjId, [0.]))

                    newCTFTomo._cutOnFreqList = CsvList(pType=float)
                    newCTFTomo.setCutOnFreqList(cutOnFreqDict.get(tiObjId, [0.]))

                newCTFTomo.completeInfoFromList()
            else:
                newCTFTomo.setWrongDefocus()
                newCTFTomo.setEnabled(False)

            newCTFTomo.setIndex(i + 1)
            newCTFTomo.setAcquisitionOrder(ti.getAcquisitionOrder())
            newCTFTomoSeries.append(newCTFTomo)

        newCTFTomoSeries.setIMODDefocusFileFlag(defocusFileFlag)
        newCTFTomoSeries.setNumberOfEstimationsInRangeFromDefocusList()

    @staticmethod
    def getDefocusFileFlag(defocusFilePath: str):
        """ This method returns the flag that indicate the
        information contained in an IMOD defocus file. The flag
        value "is the sum of:
              1 if the file has astigmatism values
              2 if the astigmatism axis angle is in radians, not degrees
              4 if the file has phase shifts
              8 if the phase shifts are in radians, not degrees
             16 if tilt angles need to be inverted to match what the
                 program expects (what Ctfplotter would produce)
                 with the -invert option
             32 if the file has cut-on frequencies attenuating the phase
                 at low frequencies"

                 from https://bio3d.colorado.edu/imod/doc/man/ctfphaseflip.html """

        with open(defocusFilePath) as f:
            lines = f.readlines()

        " File contains only defocus information (no astigmatism, no phase shift, no cut-on frequency) "
        if len(lines) == 1:
            return 0
        elif len(lines[1].split()) == 5:
            return 0

        else:
            " File contains more information apart "
            return int(lines[0].split()[0])

    def readCTFEstimationInfoFile(self, defocusFilePath, flag):
        """ This method takes an IMOD-based file path containing the
        information associated to a CTF estimation and produces a set
        of dictionaries containing the information of each parameter
        for each tilt-image belonging to the tilt-series. These
        dictionaries are readable information for Scipion, useful
        to generate the corresponding output CTFTomoSeries object. """

        # Read info as table
        ctfInfoIMODTable = self.readDefocusFileAsTable(defocusFilePath)

        if flag == 0:
            # Plain estimation
            return self.refactorCTFDefocusEstimationInfo(ctfInfoIMODTable)

        elif flag == 1:
            # Astigmatism estimation
            return self.refactorCTFDesfocusAstigmatismEstimationInfo(ctfInfoIMODTable)

        elif flag == 4:
            # Phase-shift estimation
            return self.refactorCTFDefocusPhaseShiftEstimationInfo(ctfInfoIMODTable)

        elif flag == 5:
            # Astigmatism and phase shift estimation
            return self.refactorCTFDefocusAstigmatismPhaseShiftEstimationInfo(ctfInfoIMODTable)

        elif flag == 37:
            # Astigmatism, phase shift and cut-on frequency estimation
            return self.refactorCTFDefocusAstigmatismPhaseShiftCutOnFreqEstimationInfo(ctfInfoIMODTable)

        else:
            raise ValueError("Defocus file flag do not supported. Only supported formats corresponding to flags 0, "
                             "1, 4, 5, and 37.")

    @staticmethod
    def readDefocusFileAsTable(defocusFilePath):
        """ This method takes an IMOD-based ctf estimation file path
        and returns a table containing the CTF estimation
        information of each tilt-image (per line) belonging to the tilt-series. """
        defocusTable = []

        with open(defocusFilePath) as f:
            lines = f.readlines()

            for index, line in enumerate(lines):
                vector = line.split()
                [float(i) for i in vector]

                if index == 0 and len(lines) == 1:
                    "CTF estimation is plain (no astigmatism, no phase shift, no cut-on frequency) and is the first line."
                    "Remove last element from the first line (it contains the mode of the estimation run). This case "
                    "considers that the estimation file only has one line. "
                    vector.pop()
                    defocusTable.append(vector)

                elif index == 0 and len(lines[1].split()) == 5:
                    "CTF estimation is plain (no astigmatism, no phase shift, no cut-on frequency) and is the first line."
                    "Remove last element from the first line (it contains the mode of the estimation run). "
                    vector.pop()
                    defocusTable.append(vector)

                elif index == 0 and len(lines[1].split()) != 5:
                    "CTF estimation is not plain and is the first line."
                    "Do not add this line to the table. Only contains flag and format info."
                    pass

                else:
                    "Any posterior line that is not the first one is added to the table ."
                    defocusTable.append(vector)

            return defocusTable
    
    @staticmethod
    def appendToDict(dictionary: dict,
                     index: int,
                     value: float) -> None:
        #Python Dictionary setdefault() returns the value of a key (if the key is in dictionary).
        # Else, it inserts a key with the default value to the dictionary.
        dictionary.setdefault(index, []).append(Float(value))

    def refactorCTFDefocusEstimationInfo(self, ctfInfoIMODTable):
        """ This method takes a table containing the information of
        an IMOD-based CTF estimation containing only defocus
        information (5 columns) and produces a new dictionary
        containing the same information in a format readable for
        Scipion. Flag 0. """

        if len(ctfInfoIMODTable[0]) != 5:
            raise RuntimeError("Misleading file format, CTF estimation with no astigmatism should be 5 columns long")

        defocusUDict = {}

        for element in ctfInfoIMODTable:
            start, end = int(element[0]), int(element[1])
            defocus = float(element[4]) * 10

            for index in range(start, end + 1):
                self.appendToDict(defocusUDict, index, defocus)

        return defocusUDict

    def refactorCTFDesfocusAstigmatismEstimationInfo(self, ctfInfoIMODTable):
        """ This method takes a table containing the information of an
        IMOD-based CTF estimation containing defocus and
        astigmatism information (7 columns) and produces a set
        of dictionaries table containing the same information in a
        format readable for Scipion. Flag 1. """

        if len(ctfInfoIMODTable[0]) != 7:
            raise RuntimeError("Misleading file format, CTF estimation "
                               "with astigmatism should be 7 columns long")
        defocusUDict = {}
        defocusVDict = {}
        defocusAngleDict = {}

        for element in ctfInfoIMODTable:
            start, end = int(element[0]), int(element[1])
            defocusU = float(element[4]) * 10  # From nm to angstroms
            defocusV = float(element[5]) * 10  # From nm to angstroms
            defocusAngle = float(element[6])

            # Segregate information from range
            for index in range(start, end + 1):
                self.appendToDict(defocusUDict, index, defocusU)
                self.appendToDict(defocusVDict, index, defocusV)
                self.appendToDict(defocusAngleDict, index, defocusAngle)

        return defocusUDict, defocusVDict, defocusAngleDict

    @staticmethod
    def refactorCTFDefocusPhaseShiftEstimationInfo(self, ctfInfoIMODTable):
        """ This method takes a table containing the information of
        an IMOD-based CTF estimation containing defocus, and phase
        shift information (6 columns) and produces a new set of
        dictionaries containing the same information in a format
        readable for Scipion. Flag 4. """

        if len(ctfInfoIMODTable[0]) != 6:
            raise RuntimeError(
                "Misleading file format, CTF estimation with defocus and phase shift should be 6 columns long")

        defocusUDict = {}
        phaseShiftDict = {}

        for element in ctfInfoIMODTable:
            start, end = int(element[0]), int(element[1])
            defocusU = float(element[4]) * 10
            phaseShift = float(element[5])

            for index in range(start, end + 1):
                self.appendToDict(defocusUDict, index, defocusU)
                self.appendToDict(phaseShiftDict, index, phaseShift)

        return defocusUDict, phaseShiftDict

    def refactorCTFDefocusAstigmatismPhaseShiftEstimationInfo(self, ctfInfoIMODTable):
        """ This method takes a table containing the information of
        an IMOD-based CTF estimation containing defocus, astigmatism
        and phase shift information (8 columns) and produces a new
        set of dictionaries containing the same information in a
        format readable for Scipion. Flag 5. """

        if len(ctfInfoIMODTable[0]) != 8:
            raise RuntimeError(
                "Misleading file format, CTF estimation with astigmatism and phase shift should be 8 columns long")

        defocusUDict = {}
        defocusVDict = {}
        defocusAngleDict = {}
        phaseShiftDict = {}

        for element in ctfInfoIMODTable:
            start, end = int(element[0]), int(element[1])
            defocusU = float(element[4]) * 10
            defocusV = float(element[5]) * 10
            angle = float(element[6])
            phaseShift = float(element[7])

            for index in range(start, end + 1):
                self.appendToDict(defocusUDict, index, defocusU)
                self.appendToDict(defocusVDict, index, defocusV)
                self.appendToDict(defocusAngleDict, index, angle)
                self.appendToDict(phaseShiftDict, index, phaseShift)

        return defocusUDict, defocusVDict, defocusAngleDict, phaseShiftDict

    def refactorCTFDefocusAstigmatismPhaseShiftCutOnFreqEstimationInfo(self, ctfInfoIMODTable):
        """ This method takes a table containing the information of an
        IMOD-based CTF estimation containing defocus, astigmatism, phase
        shift information and cut-on frequency (8 columns) and produces a
        new set of dictionaries containing the same information in a
        format readable for Scipion. Flag 37. """

        if len(ctfInfoIMODTable[0]) != 9:
            raise RuntimeError(
                "Misleading file format, CTF estimation with astigmatism, phase shift and cut-on frequency should be 9 columns long")

        defocusUDict = {}
        defocusVDict = {}
        defocusAngleDict = {}
        phaseShiftDict = {}
        cutOnFreqDict = {}

        for element in ctfInfoIMODTable:
            start, end = int(element[0]), int(element[1])
            defocusU = float(element[4]) * 10
            defocusV = float(element[5]) * 10
            angle = float(element[6])
            phaseShift = float(element[7])
            cutOnFreq = float(element[8])

            for index in range(start, end + 1):
                self.appendToDict(defocusUDict, index, defocusU)
                self.appendToDict(defocusVDict, index, defocusV)
                self.appendToDict(defocusAngleDict, index, angle)
                self.appendToDict(phaseShiftDict, index, phaseShift)
                self.appendToDict(cutOnFreqDict, index, cutOnFreq)

        return defocusUDict, defocusVDict, defocusAngleDict, phaseShiftDict, cutOnFreqDict

