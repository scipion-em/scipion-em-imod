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
import traceback
from os.path import exists
from typing import List, Tuple

from imod import Plugin
from imod.constants import OUTPUT_TOMOGRAMS_NAME, MRC_EXT
from imod.protocols import ProtImodBase
from pwem.convert.headers import setMRCSamplingRate
from pwem.objects import Transform
from pyworkflow import BETA
from pyworkflow.object import Integer
from pyworkflow.protocol import IntParam, STEPS_PARALLEL
from pyworkflow.utils import Message, cyanStr, redStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import SetOfTomograms, Tomogram

logger = logging.getLogger(__name__)
X_AXIS = 'x'
Y_AXIS = 'y'
Z_AXIS = 'z'
EVEN_SUFFIX = '_even'
ODD_SUFFIX = '_odd'


class ProtImodCropTomograms(ProtImodBase):
    """
    Trims and crops tomographic volumes to isolate a selected region of
interest within a larger tomogram.

    AI Generated:

    Crop Tomograms (ProtImodCropTomograms) — User Manual
        Overview

        The Crop Tomograms protocol extracts a defined subvolume from one or
        more tomograms in order to focus subsequent analysis on a biologically
        relevant region. This operation is commonly used in cryo-electron
        tomography workflows to reduce dataset size, improve visualization,
        isolate structures of interest, or prepare data for downstream tasks
        such as segmentation, subtomogram averaging, particle picking, or
        manual interpretation.

        In practical biological workflows, tomograms frequently contain large
        regions that are not relevant for the current study. These may include
        empty solvent regions, support film, fiducial markers, neighboring
        cells, or peripheral areas outside the target structure. Cropping
        allows researchers to retain only the meaningful portion of the volume,
        simplifying data handling and reducing computational demands in later
        processing stages.

        Inputs and Biological Context

        The protocol operates on reconstructed tomograms that already exist as
        volumetric datasets. These tomograms may correspond to cellular
        environments, macromolecular assemblies, membrane regions, organelles,
        or purified complexes reconstructed through cryo-electron tomography.

        The selected crop region is defined independently along the X, Y, and Z
        axes. Users can therefore isolate compact structures, elongated
        assemblies, membrane slabs, or localized intracellular regions
        depending on the biological question. The resulting cropped tomograms
        preserve the original voxel sampling and remain compatible with
        downstream tomography workflows.

        In many biological studies, cropping is not merely a convenience step
        but an important strategy for improving interpretability. Removing large
        irrelevant regions often improves visualization clarity and reduces the
        influence of unrelated densities during later analyses.

        Coordinate Selection and Volume Definition

        The protocol defines the cropped region through start and end
        coordinates along each spatial axis. These coordinates determine the
        final dimensions of the extracted volume. Because the selected limits
        are inclusive, users should carefully choose the upper boundaries when
        exact output dimensions are required.

        From a biological perspective, it is usually preferable to include a
        small safety margin around the target structure instead of cropping too
        tightly. Excessively restrictive cropping may truncate flexible domains,
        membrane extensions, or neighboring contextual features that later
        become biologically important.

        Cropping only selected axes is also possible. This flexibility is
        particularly useful when the structure occupies the full extent of one
        dimension but only a restricted region in others, such as membrane
        sheets or filamentous assemblies.

        Tomogram Geometry and Spatial Interpretation

        An important aspect of tomogram cropping is the preservation of spatial
        consistency. The protocol updates the tomogram spatial origin after
        cropping so that downstream visualization and analysis tools continue
        to interpret the coordinates correctly.

        This behavior is biologically important because cryo-electron tomography
        studies often rely on positional interpretation. Examples include
        measuring distances to membranes, locating macromolecular complexes
        inside cells, or correlating tomograms with fluorescence microscopy or
        segmentation data. Maintaining accurate geometric information helps
        preserve biological meaning after trimming the volume.

        Odd and Even Tomogram Support

        The protocol can also operate on odd and even tomogram variants when
        these datasets are available. Such datasets are commonly used in
        resolution estimation procedures or validation workflows. By applying
        the same cropping operation consistently to all related tomograms, the
        protocol preserves compatibility between the different data partitions.

        This is especially useful in workflows involving subtomogram averaging,
        local resolution estimation, or independent half-map validation, where
        maintaining strict spatial correspondence between datasets is essential.

        Outputs and Their Interpretation

        The protocol produces a new set of cropped tomograms corresponding to
        the selected spatial region. The resulting tomograms preserve the
        original biological information within the retained area while reducing
        unnecessary volume outside the region of interest.

        These cropped tomograms are typically easier to visualize, transfer,
        archive, and process computationally. Smaller datasets can
        substantially accelerate downstream analyses, particularly in workflows
        involving iterative refinement, segmentation, denoising, or machine
        learning approaches.

        Practical Recommendations

        In routine cryo-electron tomography workflows, it is generally advisable
        to inspect the tomogram carefully before defining crop coordinates.
        Including sufficient contextual information around the target structure
        often improves interpretation and prevents accidental exclusion of
        biologically relevant features.

        When preparing tomograms for subtomogram averaging or segmentation,
        moderate cropping is usually preferable to extremely tight trimming.
        Preserving nearby membrane boundaries, neighboring complexes, or local
        cellular context can remain valuable even if these features are not the
        primary target of the analysis.

        Cropping is also highly beneficial when handling very large cellular
        tomograms, where reducing the volume size can dramatically decrease
        storage requirements and processing time.

        Final Perspective

        For many cryo-electron tomography projects, tomogram cropping is a
        fundamental preparation step that improves both computational efficiency
        and biological focus. Careful selection of the cropped region allows
        researchers to concentrate analyses on the most informative part of the
        dataset while preserving the spatial and structural context necessary
        for meaningful biological interpretation.
    """

    _label = 'crop tomograms'
    _possibleOutputs = {OUTPUT_TOMOGRAMS_NAME: SetOfTomograms}
    stepsExecutionMode = STEPS_PARALLEL
    _devStatus = BETA
    program = 'trimvol'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.xyzParams = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        self.addInTomoSetFormParam(form)
        self.addOddEvenParams(form, isTomogram=True)
        g = form.addGroup('Crop values (px)')
        g.addLine('Final dim = (p1 - p0) + 1',
                  help='Example:\n\n'
                       'With x0 = 270 and x1 = 500 -> X axis: (500 - 270) + 1 = 231 px.\n\n'
                       'This is because of the cut coordinates are inclusive. To get even '
                       'dimensions, consider to add or subtract 1 to the upper limit of each '
                       'pair of coordinates.')
        l = g.addLine('Crop X',
                      help='Starting and ending X coordinate of region to cut out. '
                           'Numbers lower than 0 means to ignore the crop in that axis.')
        l.addParam(f'{X_AXIS}0', IntParam, label='x start', default=-1)
        l.addParam(f'{X_AXIS}1', IntParam, label='x end', default=-1)
        l = g.addLine('Crop Y',
                      help='Starting and ending Y coordinate of region to cut out. '
                           'Numbers lower than 0 means to ignore the crop in that axis.')
        l.addParam(f'{Y_AXIS}0', IntParam, label='y start', default=-1)
        l.addParam(f'{Y_AXIS}1', IntParam, label='y end', default=-1)
        l = g.addLine('Crop Z',
                      help='Starting and ending Z coordinate of region to cut out. '
                           'Numbers lower than 0 means to ignore the crop in that axis.')
        l.addParam(f'{Z_AXIS}0', IntParam, label='z start', default=-1)
        l.addParam(f'{Z_AXIS}1', IntParam, label='z end', default=-1)
        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []
        for tsId in self.tomoDict.keys():
            compId = self._insertFunctionStep(self.cropTomogramStep,
                                              tsId,
                                              prerequisites=[],
                                              needsGPU=False)
            outId = self._insertFunctionStep(self.generateOutputStep,
                                             tsId,
                                             prerequisites=compId,
                                             needsGPU=False)
            closeSetStepDeps.append(outId)

        self._insertFunctionStep(self.closeOutputSetsStep,
                                 OUTPUT_TOMOGRAMS_NAME,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        inTomosSet = self.getInputTomoSet()
        self.doOddEven = self.applyToOddEven(inTomosSet)
        self.tomoDict = {tomo.getTsId(): tomo.clone() for tomo in inTomosSet.iterItems()}
        params = ''
        if self._doCropInAxis(X_AXIS):
            x0, x1 = self._getAxisCoords(X_AXIS)
            params += f'-{X_AXIS} {x0},{x1} '
        if self._doCropInAxis(Y_AXIS):
            y0, y1 = self._getAxisCoords(Y_AXIS)
            params += f'-{Y_AXIS} {y0},{y1} '
        if self._doCropInAxis(Z_AXIS):
            z0, z1 = self._getAxisCoords(Z_AXIS)
            params += f'-{Z_AXIS} {z0},{z1} '
        self.xyzParams = params

    def cropTomogramStep(self, tsId: str):
        try:
            self._runCropTomo(tsId)
            if self.doOddEven:
                self._runCropTomo(tsId, suffix=EVEN_SUFFIX)
                self._runCropTomo(tsId, suffix=ODD_SUFFIX)
        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {self.program} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def generateOutputStep(self, tsId: str):
        tomo = self.tomoDict[tsId]
        if tsId in self.failedItems:
            self.addToOutFailedSet(tomo)
            return

        try:
            outputFn = self._getCroppedTomoFn(tsId)
            if exists(outputFn):
                self._registerOutput(tomo, outputFn)
            else:
                logger.error(redStr(f'tsId = {tsId} -> Output file {outputFn} '
                                    f'was not generated. Skipping... '))

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output '
                                f'with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    @retry_on_sqlite_lock(log=logger)
    def _registerOutput(self, tomo: Tomogram, outputFn: str):
        tsId = tomo.getTsId()
        with self._lock:
            # Set of tomograms
            outTomoSet = self.getOutputSetOfTomograms(self.getInputTomoSet(pointer=True))
            setMRCSamplingRate(outputFn, outTomoSet.getSamplingRate())  # Update the apix value in file header
            # Tomogram
            outTomo = Tomogram(tsId=tsId)
            outTomo.copyInfo(tomo)
            outTomo.setFileName(outputFn)
            self._manageTomoOrigin(outTomo)
            self.setTomoOddEven(tsId, outTomo)
            # Data persistence
            outTomoSet.append(outTomo)
            outTomoSet.updateDim()
            outTomoSet.update(outTomo)
            outTomoSet.write()
            self._store(outTomoSet)

    # --------------------------- INFO functions ------------------------------
    def _validate(self) -> List[str]:
        errorMsg = []

        def _validateVals(axisName: str) -> str:
            val0, val1 = self._getAxisCoords(axisName)
            if val0 > val1:
                return f'{axisName}1 must be greater than {axisName}0'
            if val0 == val1 and val0 > 0 and val1 > 0:
                return f'{axisName}0 and {axisName}1 cannot be equal and greater than 0.'
            return ''

        if (not self._doCropInAxis(X_AXIS) and
                not self._doCropInAxis(Y_AXIS) and
                not self._doCropInAxis(Z_AXIS)):
            errorMsg.append('No cropping will be carried out with the current x, y and z values.')
        xMsg = _validateVals(X_AXIS)
        if xMsg:
            errorMsg.append(xMsg)
        yMsg = _validateVals(Y_AXIS)
        if yMsg:
            errorMsg.append(yMsg)
        zMsg = _validateVals(Z_AXIS)
        if zMsg:
            errorMsg.append(zMsg)
        return errorMsg

    # --------------------------- UTILS functions -----------------------------
    def _getAxisCoords(self, axisName: str) -> Tuple[int, int]:
        val0 = getattr(self, f'{axisName}0', Integer(-1)).get()
        val1 = getattr(self, f'{axisName}1', Integer(-1)).get()
        return val0, val1

    def _doCropInAxis(self, axisName: str) -> bool:
        val0, val1 = self._getAxisCoords(axisName)
        return True if val0 >= 0 and val1 >= 0 else False

    def _getCroppedTomoFn(self, tsId: str, suffix: str = '') -> str:
        return self._getExtraPath(f'{tsId}{suffix}.mrc')

    def _runCropTomo(self, tsId: str, suffix: str = '') -> None:
        logger.info(cyanStr(f'tsId = {tsId}: Cropping the tomogram {suffix}...'))
        params = self.xyzParams
        params += f'{self.tomoDict[tsId].getFileName()} '  # input file
        params += f'{self._getCroppedTomoFn(tsId, suffix=suffix)} '  # output file
        Plugin.runImod(self, self.program, params)

    def _manageTomoOrigin(self, tomo: Tomogram) -> None:
        prevOriginShifts = tomo.getShiftsFromOrigin()
        sxPrev, syPrev, szPrev = prevOriginShifts
        sr = tomo.getSamplingRate()
        x0, x1 = self._getAxisCoords(X_AXIS)
        newSx = - sr * (x0 + (x1 - x0) / 2) if self._doCropInAxis(X_AXIS) else sxPrev
        y0, y1 = self._getAxisCoords(Y_AXIS)
        newSy = - sr * (y0 + (y1 - y0) / 2) if self._doCropInAxis(Y_AXIS) else syPrev
        z0, z1 = self._getAxisCoords(Z_AXIS)
        newSz = - sr * (z0 + (z1 - z0) / 2) if self._doCropInAxis(Z_AXIS) else szPrev
        newOrigin = Transform()
        newOrigin.setShifts(newSx, newSy, newSz)
        tomo.setOrigin(newOrigin)
