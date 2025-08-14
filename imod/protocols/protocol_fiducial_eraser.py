# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
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
from os.path import exists
import pyworkflow.protocol.params as params
from imod.protocols.protocol_base import IN_TS_SET
from pyworkflow.protocol import STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.utils import Message, cyanStr, redStr
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from imod.protocols import ProtImodBase
from imod.constants import OUTPUT_TILTSERIES_NAME, IMODFINDBEADS_PROGRAM, CCDERASER_PROGRAM, ODD, EVEN

logger = logging.getLogger(__name__)


class ProtImodFiducialEraser(ProtImodBase, ProtStreamingBase):
    """
    This protocol will erase the fiducial gold beads present in the tilt
     series images using IMOD procedures.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/imodfindbeads.html
        https://bio3d.colorado.edu/imod/doc/man/ccderaser.html

    IMODFINDBEADS finds gold particles or other circular densities (beads)
    in images by a combination of cross-correlation and other methods.  It
    starts by correlating with a model of a spherical bead of a specified
    size, then forms an average out of the most strongly-correlating subset
    and repeats the procedure by correlating with the average.  It analyzes
    the distribution of correlation strengths to find the strength that
    best separates the particles of interest from similar densities.  The
    positions of the beads are stored in an IMOD model along with strengths
    for each.  The points can then be visualized in 3dmod, and with the
    help of the Bead Fixer module, the threshold can be adjusted and points
    below threshold can be deleted.

    Rather than cross-correlating with a model or averaged bead, the pro-
    gram applies an edge-detecting filter (Sobel, by default) to both the
    images and the reference, and correlates the filtered images.  This
    method improves the detectability of the beads and may improve the
    accuracy of the center positions.  However, it only works well for
    beads in a certain size range, so the program first scales the images
    to bring the beads to a specified size (8, by default).  The peaks in
    this correlation are the set of candidate positions for the beads.

    At each position, the program then computes an integral of the bead
    density relative to the background in an annulus around the bead.  The
    program can then work with three measures of peak strength.  One is the
    strength of the Sobel-filter correlation (which includes a component
    based on the density of the bead, a factor lost when using a normalized
    correlation coefficient).  The second is the integrated density, and
    the third is the geometric mean of the first two.  Whichever measure is
    chosen, it is scaled so that the maximum value is 1.

    Correlation with a simple item like a bead always produces many more
    peaks than actual beads, but a histogram of peak strength generally
    shows a dip between actual beads and spurious peaks.  The program thus
    computes a histogram and smooths it with kernel smoothing, whereby a
    narrow distribution function instead of a single point is added into
    the histogram at every peak position.  The width of this function is
    the kernel width, referred to as H in program output.  The program
    tries a series of widths, from 0.2 downward, until it finds a dip in
    the histogram; then it computes it again with a kernel width of 0.05 in
    order to locate the dip more accurately.

    After the initial correlation with a model bead, the program uses the
    histogram analysis to select beads to average as a template for the
    second round correlation.  If the analysis fails, it is possible to
    bypass it by entering a relative peak strength to use as a criterion
    for selecting beads.  After the second round of correlations, the loca-
    tion of the dip is used to determine which points to output in the
    model.  The default is to output a number of points below the dip so
    that the user can check and adjust the threshold if necessary.  How-
    ever, with the -store option, you can output just the points above the
    dip or a certain fraction of the strongest peaks above the dip.  Or, if
    the histogram analysis fails, this option can be used to bypass it and
    specify the actual peak strength to use as the criterion for output.

    CCDERASER replaces deviant pixels with interpolated values from sur-
    rounding pixels.  It is designed to correct defects in electron micro-
    scope images from CCD cameras.  It can use two algorithms to automati-
    cally remove peaks in intensity caused by X-rays.  It can also take an
    IMOD model file with specifications of regions to be replaced; in this
    mode it can be used to erase gold fiducial markers.  With a model, the
    program can replace a group of adjacent pixels with interpolated val-
    ues, or all of the pixels along a line.  It can do this on only a spe-
    cific image, or on all of the sections in the file. The program can
    operate in trial mode, without making an output file, and it can output
    a model file with points at the pixels to be replaced.


    """

    _label = 'fiducial eraser'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @classmethod
    def worksInStreaming(cls):
        return True

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)

        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-series')

        form.addParam('fidDiameter',
                      params.FloatParam,
                      important=True,
                      default= 10.0,
                      label='Fiducial diameter (nm)')

        form.addParam('erasedDiameter',
                      params.FloatParam,
                      important=True,
                      default=12.0,
                      label='Diameter to erase (nm)')

        form.addParallelSection(threads=2, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def stepsGeneratorStep(self) -> None:
        self._initialize()
        closeSetStepDeps = []
        inTsSet = self.getInputTsSet()
        outTsSet = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        self.readingOutput(outTsSet)

        while True:
            listInTsIds = inTsSet.getTSIds()
            if not inTsSet.isStreamOpen() and self.tsIdReadList == listInTsIds:
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetsStep,
                                         OUTPUT_TILTSERIES_NAME,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break

            for ts in inTsSet.iterItems():
                tsId = ts.getTsId()
                if tsId not in self.tsIdReadList and ts.getSize() > 0:  # Avoid processing empty TS
                    cInId = self._insertFunctionStep(self.linkTsStep,
                                                     tsId,
                                                     prerequisites=[],
                                                     needsGPU=False)
                    beadsId = self._insertFunctionStep(self.imodfindbeadsStep,
                                                       tsId,
                                                       prerequisites=cInId,
                                                       needsGPU=False)
                    eraserId = self._insertFunctionStep(self.ccderaserStep,
                                                        tsId,
                                                        prerequisites=beadsId,
                                                        needsGPU=False)
                    createOutputId = self._insertFunctionStep(self.createOutputStep,
                                                              tsId,
                                                              prerequisites=eraserId,
                                                              needsGPU=False)
                    closeSetStepDeps.append(createOutputId)
                    logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                    self.tsIdReadList.append(tsId)

                self.refreshStreaming(inTsSet)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        super()._initialize()

    def imodfindbeadsStep(self, tsId: str):
        """This step creates a fiducial model"""
        try:
            logger.info(cyanStr(f'tsId = {tsId} -> Finding the fiducials...'))
            with self._lock:
                ts = self.getCurrentTs(tsId)
                firstTi = ts.getFirstItem()

            paramsImodFindBeads = {
                "-inp": firstTi.getFileName(),
                "-o": self._getModelFileName(tsId),
                "-size": self._getFiducialDiameterPx(ts)
            }
            self.runProgram(IMODFINDBEADS_PROGRAM, paramsImodFindBeads)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {IMODFINDBEADS_PROGRAM} execution '
                                f'failed with the exception -> {e}'))

    def ccderaserStep(self, tsId: str):
        """This step erase the gold beads from the fiducial model"""
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'tsId = {tsId} -> Erasing the gold beads...'))
                with self._lock:
                    ts = self.getCurrentTs(tsId)
                    firstTi = ts.getFirstItem()

                paramsCCDeraser = {
                    "-ExpandCircleIterations": 3,
                    "-BetterRadius": self._getFiducialDiameterPx(ts) / 2,
                    "-input": firstTi.getFileName(),
                    "-output": self.getExtraOutFile(tsId),
                    "-ModelFile": self._getModelFileName(tsId),
                    "-MergePatches": '',
                    "-ExcludeAdjacent": '',
                    "-CircleObjects": '/',
                    "-SkipTurnedOffPoints": 1,
                    "-PolynomialOrder": -1,
                }

                self.runProgram(CCDERASER_PROGRAM, paramsCCDeraser)

                if self.doOddEven:
                    logger.info(cyanStr(f'tsId = {tsId} -> Erasing the gold beads (ODD Tilt-series) ...'))
                    paramsCCDeraser['-input'] = ts.getOddFileName(),
                    paramsCCDeraser['-output'] = self.getExtraOutFile(tsId, suffix=ODD)
                    self.runProgram(CCDERASER_PROGRAM, paramsCCDeraser)

                    logger.info(cyanStr(f'tsId = {tsId} -> Erasing the gold beads (EVEN Tilt-series) ...'))
                    paramsCCDeraser['-input'] = ts.getEvenFileName(),
                    paramsCCDeraser['-output'] = self.getExtraOutFile(tsId, suffix=EVEN)
                    self.runProgram(CCDERASER_PROGRAM, paramsCCDeraser)

            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> {CCDERASER_PROGRAM} execution failed'
                                    f' with the exception -> {e}'))

    def createOutputStep(self, tsId: str):
        if tsId in self.failedItems:
            self.addToOutFailedSet(tsId)
            return

        try:
            outTsFile = self.getExtraOutFile(tsId)
            if exists(outTsFile):
                with self._lock:
                    ts = self.getCurrentTs(tsId)
                    # Set of tilt-series
                    outTsSet = self.getOutputSetOfTS(self.getInputTsSet(pointer=True))
                    # Tilt-series
                    outTs = TiltSeries()
                    outTs.copyInfo(ts)
                    outTsSet.append(outTs)
                    # Tilt-images
                    for ti in ts.iterItems():
                        outTi = TiltImage()
                        outTi.copyInfo(ti)
                        outTi.setFileName(outTsFile)
                        if self.doOddEven:
                            outTi.setOddEven([self.getExtraOutFile(tsId, suffix=ODD),
                                              self.getExtraOutFile(tsId, suffix=EVEN)])
                        else:
                            outTi.setOddEven([])  # the input may have odd/even but the user may have decided not
                            # to consider them in the current execution, so they should be set to empty to avoid
                            # next protocols be confused about having them.
                        outTs.append(outTi)
                    # Data persistence
                    outTs.write()
                    outTsSet.update(outTs)
                    outTsSet.write()
                    self._store(outTsSet)
                    # Close explicitly the outputs (for streaming)
                    for outputName in self._possibleOutputs.keys():
                        output = getattr(self, outputName, None)
                        if output:
                            output.close()
            else:
                logger.error(redStr(f'tsId = {tsId} -> Output file {outTsFile} was not generated. Skipping... '))

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append(f"Input tilt-series: {self.getInputTsSet().getSize()}\n"
                           "Transformation matrices calculated: "
                           f"{self.TiltSeries.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.TiltSeries:
            methods.append("The transformation matrix has been calculated for "
                           f"{self.TiltSeries.getSize()} tilt-series using "
                           "the IMOD *tiltxcorr* command.")
        return methods

    # --------------------------- UTILS functions -----------------------------
    def _getFiducialDiameterPx(self, ts: TiltSeries) -> float:
        return 10 * self.fidDiameter.get() / ts.getSamplingRate()

    def _getModelFileName(self, tsId: str) -> str:
        return self._getExtraPath(tsId, f'{tsId}_model')