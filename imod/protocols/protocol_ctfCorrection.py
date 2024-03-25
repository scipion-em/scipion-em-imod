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
from os.path import join
from pyworkflow import BETA
from pyworkflow.object import Set, String
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.emlib.image import ImageHandler
from tomo.objects import TiltSeries, TiltImage
from .. import Plugin, utils
from .protocol_base import ProtImodBase, EXT_MRCS_TS_ODD_NAME, EXT_MRCS_TS_EVEN_NAME


class ProtImodCtfCorrection(ProtImodBase):
    """
    CTF correction of a set of input tilt-series using the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/ctfphaseflip.html
    """

    _label = 'CTF correction'
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsDict = None
        self.ctfDict = None

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
        nonMatchingTsIds = []
        self._initialize()
        for tsId in self.tsDict.keys():  # Stores the steps serializing the tsId instead of the whole ts object
            matchingCtfTomoSeries = self.ctfDict.get(tsId, None)
            if matchingCtfTomoSeries:
                self._insertFunctionStep(self.convertInputStep, tsId)
                self._insertFunctionStep(self.generateDefocusFile, tsId)
                self._insertFunctionStep(self.ctfCorrection, tsId)
                self._insertFunctionStep(self.createOutputStep, tsId)
                self._insertFunctionStep(self.createOutputFailedSet, tsId)
            else:
                nonMatchingTsIds.append(tsId)
        self._insertFunctionStep(self.closeOutputSetsStep)
        if nonMatchingTsIds:
            self.matchingMsg.set(f'WARNING! No CTFTomoSeries found for the tilt-series {nonMatchingTsIds}')
            self._store(self.matchingMsg)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.matchingMsg = String()  # Update the msg for each protocol execution to avoid duplicities in the summary
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self.inputSetOfTiltSeries.get()}
        self.ctfDict = {ctf.getTsId(): ctf.clone(ignoreAttrs=[]) for ctf in self.inputSetOfCtfTomoSeries.get()}

    def convertInputStep(self, tsId, **kwargs):
        oddEvenFlag = self.applyToOddEven(self.inputSetOfTiltSeries.get())

        # Considering swapXY is required to make tilt axis vertical
        super().convertInputStep(tsId, doSwap=True, oddEven=oddEvenFlag, onlyEnabled=True)

    def generateDefocusFile(self, tsId):
        ts = self.tsDict[tsId]
        ctfTomoSeries = self.ctfDict[tsId]

        self.debug(f"Generating defocus file for {tsId} (ObjId), {tsId} (TsId)")
        # # At this step, a tlt file has been generated considering the excluded views in the current TS. If the no.
        # # tilt images is different from the no. CTFTomos, the attribute _isEnabled of the elements that are not present
        # # in the TS or are disabled tilt-images will be set to False and considered when generating the defocus IMOD
        # # files
        # if len(ts) != len(ctfTomoSeries):
        #     ctfTomoSeries.updateCtfTomoEnableFromTs(ts)
        # Compose the defocus file path
        defocusFilePath = self.getDefocusFileName(ts)
        """Generate defocus file"""
        utils.generateDefocusIMODFileFromObject(ctfTomoSeries, defocusFilePath, inputTiltSeries=ts)

    @ProtImodBase.tryExceptDecorator
    def ctfCorrection(self, tsId):
        ts = self.tsDict[tsId]
        acq = ts.getAcquisition()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        firstItem = ts.getFirstItem()

        """Run ctfphaseflip IMOD program"""
        paramsCtfPhaseFlip = {
            'inputStack': join(tmpPrefix, firstItem.parseFileName()),
            'angleFile': join(tmpPrefix, firstItem.parseFileName(extension=".tlt")),
            'outputFileName': join(extraPrefix, firstItem.parseFileName()),
            'defocusFile': self.getDefocusFileName(ts),
            'voltage': acq.getVoltage(),
            'sphericalAberration': ts.getAcquisition().getSphericalAberration(),
            'defocusTol': self.defocusTol.get(),
            'pixelSize': ts.getSamplingRate() / 10,
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

        if ts.getFirstItem().hasTransform():
            paramsCtfPhaseFlip['xformFile'] = join(tmpPrefix, firstItem.parseFileName(extension=".xf"))
            argsCtfPhaseFlip += "-TransformFile %(xformFile)s "

        Plugin.runImod(self, 'ctfphaseflip', argsCtfPhaseFlip % paramsCtfPhaseFlip)

        if self.applyToOddEven(ts):
            # ODD
            oddFn = self.getTmpTSFile(tsId, tmpPrefix=tmpPrefix, suffix=EXT_MRCS_TS_ODD_NAME)
            paramsCtfPhaseFlip['inputStack'] = oddFn
            paramsCtfPhaseFlip['outputFileName'] = join(extraPrefix, tsId+EXT_MRCS_TS_ODD_NAME)
            Plugin.runImod(self, 'ctfphaseflip', argsCtfPhaseFlip % paramsCtfPhaseFlip)

            # EVEN
            evenFn = self.getTmpTSFile(tsId, tmpPrefix=tmpPrefix, suffix=EXT_MRCS_TS_EVEN_NAME)
            paramsCtfPhaseFlip['inputStack'] = evenFn
            paramsCtfPhaseFlip['outputFileName'] = join(extraPrefix, tsId+EXT_MRCS_TS_EVEN_NAME)
            Plugin.runImod(self, 'ctfphaseflip', argsCtfPhaseFlip % paramsCtfPhaseFlip)

    def createOutputStep(self, tsId):
        if tsId not in self._failedTs:
            inTsSet = self.inputSetOfTiltSeries.get()
            outputSetOfTs = self.getOutputSetOfTiltSeries(inTsSet)
            extraPrefix = self._getExtraPath(tsId)

            newTs = TiltSeries(tsId=tsId)
            ts = self.tsDict[tsId]
            newTs.copyInfo(ts)
            newTs.setCtfCorrected(True)
            newTs.setInterpolated(True)
            acq = newTs.getAcquisition()
            acq.setTiltAxisAngle(0.)  # 0 because TS is aligned
            newTs.setAcquisition(acq)
            outputSetOfTs.append(newTs)

            ih = ImageHandler()

            for index, tiltImage in enumerate(ts):
                newTi = TiltImage()
                newTi.copyInfo(tiltImage, copyId=True, copyTM=False)
                acq = tiltImage.getAcquisition()
                acq.setTiltAxisAngle(0.)  # Is interpolated
                newTi.setAcquisition(acq)
                newTi.setLocation(index + 1, (join(extraPrefix, tiltImage.parseFileName())))
                if self.applyToOddEven(ts):
                    locationOdd = index + 1, (join(extraPrefix, tsId + EXT_MRCS_TS_ODD_NAME))
                    locationEven = index + 1, (join(extraPrefix, tsId + EXT_MRCS_TS_EVEN_NAME))
                    newTi.setOddEven([ih.locationToXmipp(locationOdd), ih.locationToXmipp(locationEven)])
                else:
                    newTi.setOddEven([])
                newTs.append(newTi)

            newTs.write(properties=False)
            outputSetOfTs.update(newTs)
            outputSetOfTs.write()
            self._store(outputSetOfTs)

    def createOutputFailedSet(self, tsId):
        ts = self.tsDict[tsId]
        super().createOutputFailedSet(ts)

    def closeOutputSetsStep(self):
        for _, output in self.iterOutputAttributes():
            output.setStreamState(Set.STREAM_CLOSED)
            output.write()
        self._store()

    # --------------------------- UTILS functions -----------------------------
    def getDefocusFileName(self, ts):
        """ Returns the path of the defocus filename based on
         the tilt series and creates the folder/s"""

        tmpPrefix = self._getTmpPath(ts.getTsId())
        path.makePath(tmpPrefix)
        defocusFn = ts.getFirstItem().parseFileName(extension=".defocus")
        defocusFilePath = join(tmpPrefix, defocusFn)
        return defocusFilePath

    # --------------------------- INFO functions ------------------------------
    def _warnings(self):
        warnings = []
        ts = self.inputSetOfTiltSeries.get()
        if not ts.getFirstItem().getFirstItem().hasTransform():
            warnings.append("Input tilt-series do not have alignment "
                            "information! The recommended workflow is to "
                            "estimate CTF on raw tilt-series and then here "
                            "provide tilt-series with alignment "
                            "(non-interpolated).")

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
