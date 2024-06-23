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

from pyworkflow.object import Set, String
import pyworkflow.protocol.params as params
from pwem.emlib.image import ImageHandler as ih
from tomo.objects import TiltSeries, TiltImage
from tomo.utils import getCommonTsAndCtfElements
from .. import Plugin, utils
from .protocol_base import ProtImodBase, DEFOCUS_EXT, TLT_EXT, XF_EXT, ODD, EVEN


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

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.matchingMsg = String()
        self.tsDict = None
        self.ctfDict = None
        self.presentTsIds = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      label="Input tilt-series",
                      pointerClass='SetOfTiltSeries',
                      help='Select the set of tilt-series to be '
                           'CTF corrected. Usually this will be the '
                           'tilt-series with alignment information.')

        form.addParam('inputSetOfCtfTomoSeries',
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
                      default='15',
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

        form.addParam('processOddEven',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=True,
                      label='Correct odd/even',
                      help='If True, the full tilt series and the associated odd/even tilt series will be processed. '
                           'The CTF correction applied to the odd/even tilt series will be exactly the same.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        for tsId in self.presentTsIds:  # Stores the steps serializing the tsId instead of the whole ts object
            presentAcqOrders = getCommonTsAndCtfElements(self.tsDict[tsId], self.ctfDict[tsId])
            self._insertFunctionStep(self.convertInputsStep, tsId, presentAcqOrders)
            self._insertFunctionStep(self.ctfCorrection, tsId)
            self._insertFunctionStep(self.createOutputStep, tsId, presentAcqOrders)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self._getInputSetOfTS()}
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

    def convertInputsStep(self, tsId, presentAcqOrders):
        # Generate the alignment-related files: xf, tlt, and a possible mrc
        super().convertInputStep(tsId,  # Considering swapXY is required to make tilt axis vertical
                                 doSwap=True,
                                 oddEven=self.applyToOddEven(self._getInputSetOfTS()),
                                 presentAcqOrders=presentAcqOrders)
        # Generate the defocus file
        self.generateDefocusFile(tsId, presentAcqOrders=presentAcqOrders)

    def ctfCorrection(self, tsId):
        try:
            ts = self.tsDict[tsId]
            acq = ts.getAcquisition()

            """Run ctfphaseflip IMOD program"""
            paramsCtfPhaseFlip = {
                'inputStack': self.getTmpOutFile(tsId),
                'angleFile': self.getExtraOutFile(tsId, ext=TLT_EXT),
                'outputFileName': self.getExtraOutFile(tsId),
                'defocusFile': self.getExtraOutFile(tsId, ext=DEFOCUS_EXT),
                'voltage': acq.getVoltage(),
                'sphericalAberration': acq.getSphericalAberration(),
                'defocusTol': self.defocusTol.get(),
                'pixelSize': ts.getSamplingRate() / 10,  # nm/px
                'amplitudeContrast': acq.getAmplitudeContrast(),
                'interpolationWidth': self.interpolationWidth.get()
            }

            argsCtfPhaseFlip = "-InputStack %(inputStack)s " \
                               "-AngleFile %(angleFile)s " \
                               "-OutputFileName %(outputFileName)s " \
                               "-DefocusFile %(defocusFile)s " \
                               "-Voltage %(voltage)d " \
                               "-SphericalAberration %(sphericalAberration)f " \
                               "-DefocusTol %(defocusTol)d " \
                               "-PixelSize %(pixelSize)f " \
                               "-AmplitudeContrast %(amplitudeContrast)f " \
                               "-InterpolationWidth %(interpolationWidth)d "

            if self.usesGpu():
                argsCtfPhaseFlip += f"-UseGPU {self.getGpuList()[0]} " \
                                    "-ActionIfGPUFails 2,2 "

            if ts.hasAlignment():
                paramsCtfPhaseFlip['xformFile'] = self.getExtraOutFile(tsId, ext=XF_EXT)
                argsCtfPhaseFlip += "-TransformFile %(xformFile)s "

            Plugin.runImod(self, 'ctfphaseflip', argsCtfPhaseFlip % paramsCtfPhaseFlip)

            if self.applyToOddEven(ts):
                # ODD
                paramsCtfPhaseFlip['inputStack'] = self.getTmpOutFile(tsId, suffix=ODD)
                paramsCtfPhaseFlip['outputFileName'] = self.getExtraOutFile(tsId, suffix=ODD)
                Plugin.runImod(self, 'ctfphaseflip', argsCtfPhaseFlip % paramsCtfPhaseFlip)

                # EVEN
                paramsCtfPhaseFlip['inputStack'] = self.getTmpOutFile(tsId, suffix=EVEN)
                paramsCtfPhaseFlip['outputFileName'] = self.getExtraOutFile(tsId, suffix=EVEN)
                Plugin.runImod(self, 'ctfphaseflip', argsCtfPhaseFlip % paramsCtfPhaseFlip)

        except Exception as e:
            self._failedTs.append(tsId)
            self.error('Ctf correction execution failed for tsId %s -> %s' % (tsId, e))

    def createOutputStep(self, tsId, presentAcqOrders):
        inTsSet = self._getInputSetOfTS()
        if tsId in self._failedTs:
            outputSetOfFailedTs = self.getOutputFailedSetOfTiltSeries(inTsSet)
            ts = self.tsDict[tsId]
            newTs = ts.clone()
            newTs.copyInfo(ts)
            outputSetOfFailedTs.append(newTs)
            newTs.copyItems(ts)
            newTs.write(properties=False)
            outputSetOfFailedTs.update(newTs)
            outputSetOfFailedTs.write()
            self._store(outputSetOfFailedTs)
        else:
            outputSetOfTs = self.getOutputSetOfTS(inTsSet)
            newTs = TiltSeries(tsId=tsId)
            ts = self.tsDict[tsId]
            newTs.copyInfo(ts)
            newTs.setCtfCorrected(True)
            newTs.setInterpolated(True)
            acq = newTs.getAcquisition()
            acq.setTiltAxisAngle(0.)  # 0 because TS is aligned
            newTs.setAcquisition(acq)
            outputSetOfTs.append(newTs)

            for index, inTi in enumerate(ts):
                if inTi.getAcquisitionOrder() in presentAcqOrders:
                    newTi = TiltImage()
                    newTi.copyInfo(inTi, copyId=True, copyTM=False)
                    acq = inTi.getAcquisition()
                    acq.setTiltAxisAngle(0.)  # Is interpolated
                    newTi.setAcquisition(acq)
                    newTi.setLocation(index + 1, self.getExtraOutFile(tsId))
                    if self.applyToOddEven(ts):
                        locationOdd = index + 1, self.getExtraOutFile(tsId, suffix=ODD)
                        locationEven = index + 1, self.getExtraOutFile(tsId, suffix=EVEN)
                        newTi.setOddEven([ih.locationToXmipp(locationOdd), ih.locationToXmipp(locationEven)])
                    else:
                        newTi.setOddEven([])
                    newTs.append(newTi)

            newTs.write(properties=False)
            outputSetOfTs.update(newTs)
            outputSetOfTs.write()
            self._store(outputSetOfTs)

    def createOutputFailedStep(self, tsId, presentAcqOrders):
        ts = self.tsDict[tsId]
        super().createOutputFailedSet(ts, presentAcqOrders=presentAcqOrders)

    def closeOutputSetsStep(self):
        for _, output in self.iterOutputAttributes():
            output.setStreamState(Set.STREAM_CLOSED)
            output.write()
        self._store()

    # --------------------------- UTILS functions -----------------------------
    def generateDefocusFile(self, tsId, presentAcqOrders=None):
        ts = self.tsDict[tsId]
        ctfTomoSeries = self.ctfDict[tsId]

        self.debug(f"Generating defocus file for {tsId} (ObjId), {tsId} (TsId)")
        # Compose the defocus file path
        defocusFilePath = self.getExtraOutFile(tsId, ext=DEFOCUS_EXT)
        """Generate defocus file"""
        utils.generateDefocusIMODFileFromObject(ctfTomoSeries, defocusFilePath,
                                                inputTiltSeries=ts,
                                                presentAcqOrders=presentAcqOrders)

    # --------------------------- INFO functions ------------------------------
    def _warnings(self):
        warnings = []
        for ts in self._getInputSetOfTS():
            if not ts.hasAlignment():
                warnings.append("Input tilt-series do not have alignment "
                                "information! The recommended workflow is to "
                                "estimate CTF on raw tilt-series and then here "
                                "provide tilt-series with alignment "
                                "(non-interpolated).")
                break

        return warnings

    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append("Input tilt-series: %d\nCTF corrections applied: %d"
                           % (self.inputSetOfCtfTomoSeries.get().getSize(),
                              self.TiltSeries.getSize()))
        else:
            summary.append("Outputs are not ready yet.")
        if self.matchingMsg.get():
            summary.append(self.matchingMsg.get())
        return summary

    def _methods(self):
        methods = []
        if self.TiltSeries:
            methods.append("%d tilt-series have been CTF corrected "
                           "using the IMOD *ctfphaseflip* program."
                           % (self.TiltSeries.getSize()))
        return methods
