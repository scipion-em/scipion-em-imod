# *****************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *              Scipion Team (scipion@cnb.csic.es) [1]
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

from imod.protocols.protocol_base import IN_TS_SET, IN_CTF_TOMO_SET
from pwem import ALIGN_NONE
from pyworkflow.object import String
import pyworkflow.protocol.params as params
from pwem.emlib.image import ImageHandler as ih
from pyworkflow.utils import Message
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries
from tomo.utils import getCommonTsAndCtfElements

from imod import utils
from imod.protocols import ProtImodBase
from imod.constants import (DEFOCUS_EXT, TLT_EXT, XF_EXT, ODD,
                            EVEN, OUTPUT_TILTSERIES_NAME)


class ProtImodCtfCorrection(ProtImodBase):
    """
    CTF correction of a set of input tilt-series using the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/ctfphaseflip.html

    This program will correct the CTF of an input tilt series by phase
    flipping, with an option to attenuate frequencies near the zeros of the
    CTF.

    Ctfphaseflip corrects each view strip by strip.  A strip is defined as
    an image region whose defocus difference is less than a user specified
    value, the defocus tolerance.  Normally, the strips are vertically ori-
    ented and defocus is assumed to be the same along a vertical line.
    Thus, the tilt series must be aligned so that the tilt axis is vertical
    before applying this correction.  The original thinking was that an
    image region with defocus difference less than the tolerance could be
    considered to have constant defocus and could be corrected as one
    strip.  However, the validity of the correction at the center of the
    strip probably does not depend on whether it contains material beyond
    this focus range, since only vertical lines near or at the center are
    used in the corrected image.  The program may limit the width further
    to reduce computation time, or expand it to retain enough resolution
    between successive zeros in the X direction of frequency space.

    Through most of the image, each strip is corrected based on the defocus
    at the center of the strip.  However, the strips at the left and right
    edges of the image may be corrected repeatedly, at different defocus
    values, in order to extend the correction close enough to the edges of
    the image.


    """

    _label = 'CTF correction'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.matchingMsg = String()
        self.ctfDict = None
        self.presentTsIds = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)

        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      label="Input tilt-series",
                      pointerClass='SetOfTiltSeries',
                      help='Select the set of tilt-series to be '
                           'CTF corrected. Usually this will be the '
                           'tilt-series with alignment information.')

        form.addParam(IN_CTF_TOMO_SET,
                      params.PointerParam,
                      label="Input CTF estimation",
                      pointerClass='SetOfCTFTomoSeries',
                      help='Select the CTF estimation for the set '
                           'of tilt-series.')

        form.addParam('defocusTol',
                      params.IntParam,
                      label='Defocus tolerance (nm)',
                      default=200,
                      important=True,
                      help='The value introduced must be the same used for '
                           'CTF estimation with IMOD.\n\n'
                           'Defocus tolerance in nanometers defining the '
                           'center strips. The center strips are taken '
                           'from the central region of a view that has defocus '
                           'difference less than this tolerance. '
                           'These kind of center strips from all views within '
                           'AngleRange are considered to have a '
                           'constant defocus and are used to compute the '
                           'initial CTF after being further tessellated '
                           'into tiles.')

        form.addParam('interpolationWidth',
                      params.IntParam,
                      label='Interpolation width (px)',
                      default=15,
                      important=True,
                      help="The distance in pixels between the center lines "
                           "of two consecutive strips. A pixel inside the "
                           "region between those two center lines resides in "
                           "both strips. As the two strips are corrected "
                           "separately, that pixel will have 2 corrected "
                           "values. The final value for that pixel is a "
                           "linear interpolation of the 2 corrected "
                           "values. If a value of 1 is entered, there is "
                           "no such interpolation. For a value greater "
                           "than one, the entered value will be used "
                           "whenever the strip width is less than 256 "
                           "(i.e., at high tilt), and the value will be "
                           "scaled proportional to the strip width for widths "
                           "above 256.  This scaling keeps the computational "
                           "time down and is reasonable because the defocus "
                           "difference between adjacent wide strips at "
                           "wider intervals is still less than that between "
                           "the narrower strips at high tilt. However, strips "
                           "at constant spacing can still be obtained by "
                           "entering the negative of the desired spacing, "
                           "which disables the scaling of the spacing.")

        form.addHidden(params.USE_GPU,
                       params.BooleanParam,
                       default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation."
                            "Select the one you want to use.")

        form.addHidden(params.GPU_LIST,
                       params.StringParam,
                       default='0',
                       label="Choose GPU IDs",
                       help="GPU ID. To pick the best available one set 0. "
                            "For a specific GPU set its number ID "
                            "(starting from 1).")

        self.addOddEvenParams(form)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        pIdList = []
        for tsId in self.presentTsIds:
            presentAcqOrders = getCommonTsAndCtfElements(self.tsDict[tsId], self.ctfDict[tsId])
            pidConvert = self._insertFunctionStep(self.convertInputsStep, tsId, presentAcqOrders, prerequisites=[])
            pidProcess = self._insertFunctionStep(self.ctfCorrection, tsId, prerequisites=pidConvert)
            pidCreateOutput = self._insertFunctionStep(self.createOutputStep, tsId, presentAcqOrders,
                                                       prerequisites=pidProcess)
            pIdList.append(pidCreateOutput)
        self._insertFunctionStep(self.closeOutputSetsStep, prerequisites=pIdList)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        tsSet = self.getInputSet()
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in tsSet}
        self.ctfDict = {ctf.getTsId(): ctf.clone(ignoreAttrs=[]) for ctf in self.inputSetOfCtfTomoSeries.get()}
        # Manage the present and not present tsIds
        tsIds = list(self.tsDict.keys())
        ctfTsIds = list(self.ctfDict.keys())
        self.presentTsIds = set(tsIds) & set(ctfTsIds)
        allTsIds = set(tsIds + ctfTsIds)
        nonMatchingTsIds = [tsId for tsId in allTsIds if tsId not in self.presentTsIds]
        # Update the msg for each protocol execution to avoid duplicities in the summary
        if nonMatchingTsIds:
            self.matchingMsg.set(f'WARNING! No CTFTomoSeries found for the tilt-series: {nonMatchingTsIds}')
            self._store(self.matchingMsg)
        self.sRate = tsSet.getSamplingRate()
        self.acq = tsSet.getAcquisition()
        self.oddEvenFlag = self.applyToOddEven(self.getInputSet())

    def convertInputsStep(self, tsId, presentAcqOrders):
        # Generate the alignment-related files: xf, tlt, and a possible mrc
        super().convertInputStep(tsId,  # Considering swapXY is required to make tilt axis vertical
                                 doSwap=True,
                                 oddEven=self.oddEvenFlag,
                                 presentAcqOrders=presentAcqOrders)
        # Generate the defocus file
        self.generateDefocusFile(tsId, presentAcqOrders=presentAcqOrders)

    def ctfCorrection(self, tsId):
        try:
            ts = self.tsDict[tsId]

            paramsCtfPhaseFlip = {
                "-InputStack": self.getTmpOutFile(tsId),
                "-AngleFile": self.getExtraOutFile(tsId, ext=TLT_EXT),
                "-OutputFileName": self.getExtraOutFile(tsId),
                "-DefocusFile": self.getExtraOutFile(tsId, ext=DEFOCUS_EXT),
                "-Voltage": int(self.acq.getVoltage()),
                "-SphericalAberration": self.acq.getSphericalAberration(),
                "-DefocusTol": self.defocusTol.get(),
                "-PixelSize": self.sRate / 10,  # nm
                "-AmplitudeContrast": self.acq.getAmplitudeContrast(),
                "-InterpolationWidth": self.interpolationWidth.get()
            }

            if self.usesGpu():
                paramsCtfPhaseFlip["-UseGPU"] = self.getGpuList()[0]
                paramsCtfPhaseFlip["-ActionIfGPUFails"] = "2,2"

            if ts.hasAlignment():
                paramsCtfPhaseFlip["-TransformFile"] = self.getExtraOutFile(tsId, ext=XF_EXT)

            self.runProgram('ctfphaseflip', paramsCtfPhaseFlip)

            if self.oddEvenFlag:
                # ODD
                paramsCtfPhaseFlip["-InputStack"] = self.getTmpOutFile(tsId, suffix=ODD)
                paramsCtfPhaseFlip["-OutputFileName"] = self.getExtraOutFile(tsId, suffix=ODD)
                self.runProgram('ctfphaseflip', paramsCtfPhaseFlip)

                # EVEN
                paramsCtfPhaseFlip["-InputStack"] = self.getTmpOutFile(tsId, suffix=EVEN)
                paramsCtfPhaseFlip["-OutputFileName"] = self.getExtraOutFile(tsId, suffix=EVEN)
                self.runProgram('ctfphaseflip', paramsCtfPhaseFlip)

        except Exception as e:
            self._failedItems.append(tsId)
            self.error(f"ctfphaseflip execution failed for tsId {tsId} -> {e}")

    def createOutputStep(self, tsId, presentAcqOrders):
        ts = self.tsDict[tsId]
        if tsId in self._failedItems:
            self.createOutputFailedSet(ts)
        else:
            outputFn = self.getExtraOutFile(tsId)
            if os.path.exists(outputFn):
                inTsSet = self.getInputSet(pointer=True)
                outputSetOfTs = self.getOutputSetOfTS(inTsSet)
                newTs = TiltSeries(tsId=tsId)
                ts = self.tsDict[tsId]
                newTs.copyInfo(ts)
                newTs.setAlignment(ALIGN_NONE)
                newTs.setAnglesCount(len(presentAcqOrders))
                newTs.setCtfCorrected(True)
                newTs.setInterpolated(True)
                newTs.getAcquisition().setTiltAxisAngle(0.)  # 0 because TS is aligned
                outputSetOfTs.append(newTs)

                angleMin = 999
                angleMax = -999
                tiList = []
                for index, inTi in enumerate(ts.iterItems()):
                    if inTi.getAcquisitionOrder() in presentAcqOrders:
                        newTi = TiltImage(tsId=tsId)
                        newTi.copyInfo(inTi, copyId=True, copyTM=False)
                        acq = inTi.getAcquisition()
                        acq.setTiltAxisAngle(0.)  # Is interpolated
                        newTi.setAcquisition(acq)
                        newTi.setLocation(index + 1, outputFn)
                        if self.oddEvenFlag:
                            locationOdd = index + 1, self.getExtraOutFile(tsId, suffix=ODD)
                            locationEven = index + 1, self.getExtraOutFile(tsId, suffix=EVEN)
                            newTi.setOddEven([ih.locationToXmipp(locationOdd),
                                              ih.locationToXmipp(locationEven)])
                        else:
                            newTi.setOddEven([])
                        # Update the acquisition of the TS. The accumDose, angle min and angle max for the re-stacked TS, as
                        # these values may change if the removed tilt-images are the first or the last, for example.
                        tiAngle = newTi.getTiltAngle()
                        angleMin = min(tiAngle, angleMin)
                        angleMax = max(tiAngle, angleMax)

                        tiList.append(newTi)
                    
                if len(presentAcqOrders) != max(len(ts), len(self.ctfDict[tsId])):
                    # Update the acquisition minAngle and maxAngle values of the tilt-series
                    acq = newTs.getAcquisition()
                    acq.setAngleMin(angleMin)
                    acq.setAngleMax(angleMax)
                    acq.setAccumDose(0)
                    acq.setDoseInitial(0)
                    newTs.setAcquisition(acq)
                    # Update the acquisition minAngle and maxAngle values of each tilt-image acq while preserving their
                    # specific accum and initial dose values
                    for tiOut in tiList:
                        tiAcq = tiOut.getAcquisition()
                        tiAcq.setAngleMin(angleMin)
                        tiAcq.setAngleMax(angleMax)
                        tiAcq.setAccumDose(0)
                        tiAcq.setDoseInitial(0)
                        newTs.append(tiOut)
                    newTs.setAnglesCount(len(newTs))
                else:
                    for tiOut in tiList:
                        newTs.append(tiOut)

                outputSetOfTs.update(newTs)
                self._store(outputSetOfTs)

    # --------------------------- UTILS functions -----------------------------
    def generateDefocusFile(self, tsId, presentAcqOrders=None):
        ts = self.tsDict[tsId]
        ctfTomoSeries = self.ctfDict[tsId]

        self.debug(f"Generating defocus file for {tsId}")
        defocusFilePath = self.getExtraOutFile(tsId, ext=DEFOCUS_EXT)
        utils.generateDefocusIMODFileFromObject(ctfTomoSeries, defocusFilePath,
                                                inputTiltSeries=ts,
                                                presentAcqOrders=presentAcqOrders)

    # --------------------------- INFO functions ------------------------------
    def _warnings(self):
        warnings = []
        for ts in self.getInputSet():
            if not ts.hasAlignment():
                warnings.append(f"Input tilt-series {ts.getTsId()} does not have "
                                "alignment information! The recommended workflow is to "
                                "estimate CTF on raw tilt-series and then here "
                                "provide tilt-series with alignment "
                                "(non-interpolated).")
                break

        return warnings

    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append(f"Input tilt-series: {self.getInputSet().getSize()}\n"
                           f"CTF corrections applied: {self.TiltSeries.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        if self.matchingMsg.get():
            summary.append(self.matchingMsg.get())
        return summary

    def _methods(self):
        methods = []
        if self.TiltSeries:
            methods.append(f"{self.TiltSeries.getSize()} tilt-series have been "
                           "CTF corrected using the IMOD *ctfphaseflip* program.")
        return methods
