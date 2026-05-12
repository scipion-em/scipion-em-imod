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
import os
import numpy as np
import pyworkflow as pw
import pyworkflow.protocol.params as params
from imod.convert.convert import readXfFile, fiducialModel2List, fidResidualModel2List
from pyworkflow.protocol.constants import STEPS_SERIAL
from pyworkflow.object import Set
import pyworkflow.utils as pwutils
from pwem.emlib.image import ImageHandler as ih
from pwem.objects import Transform
import tomo.objects as tomoObj
from imod import Plugin, utils
from imod.protocols import ProtImodBase
from imod.protocols.protocol_base import IN_TS_SET
from imod.constants import *
from tomo.objects import TiltSeries

logger = logging.getLogger(__name__)


class ProtImodEtomo(ProtImodBase):
    """
    Provides an interactive environment for manual tomographic reconstruction
    workflows using the IMOD eTomo package. The protocol is designed to guide
    users through the complete reconstruction process for single-axis tilt
    series, allowing manual supervision and refinement of alignment,
    fiducial-based tracking, preprocessing, and tomogram generation. More info:
    https://bio3d.colorado.edu/imod/doc/etomoTutorial.html

    AI Generated:

    Etomo Interactive Reconstruction (ProtImodEtomo) - User Manual
        Overview

        The Etomo Interactive protocol provides a bridge between Scipion and the
        IMOD eTomo environment for users who require detailed manual control over
        tomographic reconstruction workflows. Unlike fully automated pipelines,
        this protocol is intended for situations where visual inspection,
        iterative adjustment, and expert supervision are essential for obtaining
        optimal reconstructions.

        In cryo-electron tomography, reconstruction quality depends strongly on
        accurate tilt-series alignment, reliable fiducial tracking, proper
        handling of excluded views, and appropriate reconstruction parameters.
        Biological samples often contain uneven contrast, contamination,
        radiation damage, or structural flexibility that may limit the
        performance of fully automated methods. This protocol allows users to
        interact directly with the reconstruction process and refine each stage
        according to the characteristics of the specimen.

        The protocol is particularly valuable for difficult datasets, including
        thick cellular tomograms, low signal-to-noise acquisitions, samples with
        sparse fiducials, or tilt series affected by beam-induced motion and
        imaging artifacts. In these scenarios, manual intervention frequently
        improves alignment stability and final tomogram quality.

        Inputs and Reconstruction Workflow

        The protocol requires a set of tilt series as input. Each tilt series is
        prepared for processing within the eTomo environment, where users can
        interactively execute alignment, tracking, reconstruction, and
        postprocessing tasks.

        Fiducial marker information plays a central role in the workflow. Users
        typically define the approximate diameter of gold beads present in the
        sample so that eTomo can identify and track fiducial markers throughout
        the tilt range. Accurate fiducial detection is critical because the
        final alignment quality directly influences tomogram resolution and
        interpretability.

        The protocol also supports workflows in which previous alignment
        transformations already exist. In these cases, users may choose whether
        the reconstruction process should preserve and refine the incoming
        geometric alignment information or start from the original unaligned
        tilt series.

        Biological Importance of Manual Alignment

        Tomographic reconstruction is highly sensitive to alignment accuracy.
        Small errors in tilt geometry or fiducial positioning can propagate into
        severe reconstruction artifacts, including elongation effects, blurred
        densities, or loss of structural detail.

        Manual alignment inspection becomes especially important for biological
        specimens with heterogeneous thickness, uneven fiducial distribution, or
        local deformations. Users can visually evaluate tracking quality,
        correct problematic fiducials, remove unreliable tilt images, and adjust
        alignment parameters iteratively until the reconstructed volume becomes
        biologically meaningful.

        This flexibility is one of the primary advantages of the eTomo workflow.
        Instead of relying exclusively on automated optimization, researchers
        can apply domain knowledge about specimen geometry, expected structural
        features, and experimental artifacts.

        Excluded Views and Dataset Quality

        During tomography acquisition, some tilt images may become unusable due
        to contamination, charging, severe drift, ice thickness variation, or
        beam damage. Retaining these problematic views may significantly degrade
        alignment accuracy and reconstruction quality.

        The protocol supports workflows in which selected views are excluded
        during reconstruction while preserving the integrity of the remaining
        tilt series. From a biological perspective, removing a limited number of
        low-quality images often improves the final tomogram more than retaining
        all acquired data.

        However, excessive exclusion may reduce angular coverage and worsen the
        missing wedge problem. Users should therefore balance data quality
        against angular completeness when deciding which views to discard.

        Fiducial Models and Landmark Coordinates

        The protocol can generate fiducial landmark models and three-dimensional
        coordinate information associated with aligned tilt series. These outputs
        are biologically useful for evaluating alignment precision and for
        downstream analyses involving particle localization or spatial
        interpretation of intracellular structures.

        Landmark residuals provide an indirect estimate of alignment quality.
        Large residual errors may indicate unstable fiducial tracking, specimen
        deformation, or insufficient marker coverage. In practical workflows,
        users often inspect these residuals before accepting a reconstruction as
        suitable for interpretation or subtomogram analysis.

        Reconstruction Outputs

        Multiple reconstruction products may be generated depending on the
        stages completed within eTomo. These outputs can include prealigned tilt
        series, fully aligned tilt series, original tilt series carrying updated
        alignment information, raw reconstructed tomograms, and postprocessed
        tomograms.

        Prealigned stacks are typically useful for evaluating coarse alignment
        quality and fiducial tracking consistency. Fully aligned tilt series are
        usually employed for reconstruction and downstream subtomogram
        processing. Tomograms generated after postprocessing may contain
        filtering or refinement steps intended to improve visualization and
        interpretability.

        Biological users should remember that postprocessing can alter apparent
        contrast and density continuity. Quantitative interpretation should
        therefore consider the reconstruction history and any applied filtering
        operations.

        Interactive Workflow Philosophy

        This protocol is intentionally designed around interactive execution.
        Rather than enforcing a rigid automated pipeline, it provides an
        environment in which reconstruction decisions remain under user control.
        This philosophy reflects the reality that many tomographic datasets
        require case-by-case optimization.

        Experienced users often iterate several times through alignment,
        fiducial refinement, exclusion of problematic views, and reconstruction
        parameter tuning before obtaining satisfactory results. The protocol
        supports this exploratory workflow while preserving integration with the
        broader Scipion ecosystem.

        Practical Recommendations

        For most biological applications, users should begin by carefully
        inspecting fiducial visibility and tilt-series quality before launching
        reconstruction. Ensuring a consistent fiducial distribution across the
        field of view significantly improves alignment robustness.

        When reconstruction artifacts appear, common corrective actions include
        excluding damaged tilt images, refining fiducial tracking, adjusting
        alignment parameters, or reviewing the tilt-axis orientation. Thick or
        crowded cellular samples may require additional manual supervision due
        to overlapping densities and reduced fiducial contrast.

        It is also advisable to compare reconstructed tomograms before and after
        postprocessing. Although filtering often improves visual appearance,
        excessive processing can obscure subtle biological features or introduce
        misleading interpretations.

        Final Perspective

        Interactive tomographic reconstruction remains one of the most important
        stages in cryo-electron tomography because reconstruction quality
        directly determines the interpretability of biological structures.
        Careful supervision of alignment, fiducial tracking, and reconstruction
        parameters can substantially improve the quality of final tomograms,
        especially for challenging experimental datasets.

        The Etomo Interactive protocol provides users with the flexibility and
        control necessary to adapt reconstruction strategies to the biological
        complexity of real specimens while maintaining compatibility with
        integrated Scipion tomography workflows.
    """

    _label = 'Etomo interactive'
    _possibleOutputs = {
        OUTPUT_TILTSERIES_NAME: tomoObj.TiltSeries,
        OUTPUT_PREALI_TILTSERIES_NAME: tomoObj.TiltSeries,
        OUTPUT_ALI_TILTSERIES_NAME: tomoObj.TiltSeries,
        OUTPUT_TS_COORDINATES_NAME: tomoObj.SetOfTiltSeriesCoordinates,
        OUTPUT_FIDUCIAL_NO_GAPS_NAME: tomoObj.SetOfLandmarkModels,
        "FullTomograms": tomoObj.SetOfTomograms,
        "PostProcessTomograms": tomoObj.SetOfTomograms
    }
    stepsExecutionMode = STEPS_SERIAL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series',
                      help='Input set of tilt-series to be processed with eTomo.')

        form.addParam('markersDiameter',
                      params.FloatParam,
                      default=10,
                      label='Fiducial markers diameter (nm)',
                      help='Diameter of gold beads in nanometers. Note that fiducials are'
                           'small gold beads that are used as marker elements in the images.'
                           'They can be used as reference points to align the tilt series')

        form.addParam('applyAlignment',
                      params.BooleanParam,
                      default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Apply transformation matrix?',
                      help='Apply the transformation matrix if input'
                           'tilt series have it.')

        form.addParallelSection(threads=3, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    # Overwrite the following function to prevent streaming from base class
    def _stepsCheck(self):
        pass

    def _insertAllSteps(self):
        self.inputTiltSeries = None
        self._insertFunctionStep(self.runEtomoStep,
                                 interactive=True,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def runEtomoStep(self):
        from imod.viewers import ImodGenericView
        setOftiltSeries = self.getInputTsSet()
        view = ImodGenericView(None, self, setOftiltSeries,
                               isInteractive=True)
        view.show()
        self.createOutput()

    def runAllSteps(self, obj):
        """ Used by the viewer!!"""
        for item in self.getInputTsSet():
            if item.getTsId() == obj.getTsId():
                self.runEtomo(item)
                break

    def convertInputStep(self, ts, **kwargs):
        tsId = ts.getTsId()
        acq = ts.getAcquisition()
        pixSize = ts.getSamplingRate() / 10.  # nm
        self.genTsPaths(tsId)
        firstItem = ts.getFirstItem()
        inputTsFileName = firstItem.getFileName()

        outputTsFileName = self.getExtraOutFile(tsId, ext=MRC_EXT)

        """Apply the transformation from the input tilt-series"""
        if self.applyAlignment:
            ts.applyTransform(outputTsFileName)
        else:
            pwutils.createAbsLink(os.path.abspath(inputTsFileName),
                                  outputTsFileName)

        """Generate angle file"""
        angleFilePath = self.getExtraOutFile(tsId, ext=RAWTLT_EXT)
        ts.generateTltFile(angleFilePath)

        """Generate etomo config file"""
        copytomoParams = {
            '-name': pwutils.removeBaseExt(outputTsFileName),
            '-gold': self.markersDiameter,
            '-pixel': pixSize,
            '-rotation': acq.getTiltAxisAngle(),
            '-userawtlt': "",
            '-fei': 1,
            '-change': Plugin.getHome("SystemTemplate/cryoSample.adoc"),
            '-NamingStyle': 1,  # mrc style,
            '-StackExtension': 'mrc',
            '-binning': 1.0,
            '-Cs': acq.getSphericalAberration(),
            '-voltage': int(acq.getVoltage())
        }

        if ts.getExcludedViewsIndex():
            copytomoParams["-ViewsToSkip"] = ",".join(ts.getExcludedViewsIndex(caster=str))

        self.runProgram('copytomocoms', copytomoParams,
                        cwd=self._getExtraPath(tsId))

        edfFn = self.getExtraOutFile(tsId, ext=EDF_EXT)
        minTilt = min(utils.formatAngleList(angleFilePath))
        self._writeEtomoEdf(edfFn,
                            {
                                'date': pw.utils.prettyTime(),
                                'name': pwutils.removeBaseExt(outputTsFileName),
                                'pixelSize': pixSize,
                                'version': pw.__version__,
                                'minTilt': minTilt,
                                'markerDiameter': self.markersDiameter,
                                'rotationAngle': acq.getTiltAxisAngle(),
                                'imodDir': Plugin.getHome(),
                                'useCpu': self.numberOfThreads > 1
                            })

    def runEtomo(self, ts):
        tsId = ts.getTsId()
        edfFilePath = self.getExtraOutFile(tsId, ext=EDF_EXT)

        if not os.path.exists(edfFilePath):
            self.convertInputStep(ts)

        if ts is not None:
            params = {"--fg": self.getOutTsFileName(tsId, ext=EDF_EXT)}
            self.runProgram('etomo', params, cwd=self._getExtraPath(tsId))

    def createOutput(self):
        inputTS = self.getInputTsSet()
        outputTs = None # getattr(self, OUTPUT_TILTSERIES_NAME, None)  # original TS with new alignment
        outputPreAliTs = None  # getattr(self, OUTPUT_PREALI_TILTSERIES_NAME, None)
        outputAliTs = None  # getattr(self, OUTPUT_ALI_TILTSERIES_NAME, None)
        outputTSCoords = None  # getattr(self, OUTPUT_TS_COORDINATES_NAME, None)
        outputTomos = None  # getattr(self, "FullTomograms", None)
        outputPostProcessTomos = None  # getattr(self, "PostProcessTomograms", None)
        outputSetOfLandmarkModelsNoGaps = None
        for ts in inputTS:
            self.inputTiltSeries = ts
            tsId = ts.getTsId()

            self.debug(f"Registering output for {tsId}.")

            """Prealigned tilt-series"""
            prealiFilePath = self.getExtraOutFile(tsId, suffix="preali", ext=MRC_EXT)
            if os.path.exists(prealiFilePath):
                xPrealiDims, newPixSize = self.getNewPixAndDim(prealiFilePath)
                self.debug(f"{prealiFilePath}: pix = {newPixSize}, dims = {xPrealiDims}")
                if outputPreAliTs is None:
                    outputPreAliTs = self._createSetOfTiltSeries(suffix='_prealigned')
                    outputPreAliTs.copyInfo(inputTS)
                    outputPreAliTs.setSamplingRate(newPixSize)
                    outputPreAliTs.setStreamState(Set.STREAM_OPEN)
                    self._defineOutputs(**{OUTPUT_PREALI_TILTSERIES_NAME: outputPreAliTs})
                    self._defineSourceRelation(self.getInputTsSet(pointer=True),
                                           outputPreAliTs)
                # else:
                #     outputPreAliTs.enableAppend()

                newTs = TiltSeries()
                newTs.copyInfo(ts)
                newTs.setInterpolated(True)
                outputPreAliTs.append(newTs)

                # Getting the excluded views in order to disable the
                # prealigned tilt-series
                xcorrFn = self._getExtraPath(tsId, 'xcorr.com')
                excludedViewList = self.getExcludedViewList(xcorrFn,
                                                            reservedWord='SkipViews')
                self.debug(f"Excluding views for xcorr: {excludedViewList}")

                for tiltImage in ts.iterItems(iterate=False):
                    newTi = tiltImage.clone()
                    newTi.copyInfo(tiltImage, copyId=True, copyTM=False)
                    newTi.setAcquisition(tiltImage.getAcquisition())
                    newTi.setOddEven([])
                    sliceIndex = newTi.getIndex()
                    newTi.setLocation(sliceIndex, prealiFilePath)
                    if sliceIndex in excludedViewList:
                        newTi.setEnabled(False)
                    newTs.append(newTi)

                newTs.setDim(xPrealiDims)

                outputPreAliTs.update(newTs)
                outputPreAliTs.write()

            """Aligned tilt-series"""
            aligFilePath = self.getExtraOutFile(tsId, suffix="ali", ext=MRC_EXT)
            if os.path.exists(aligFilePath):
                aliDims, newPixSize = self.getNewPixAndDim(aligFilePath)
                self.debug(f"{aligFilePath}: pix = {newPixSize}, dims = {aliDims}")

                if outputAliTs is None:
                    outputAliTs = self._createSetOfTiltSeries(suffix='_aligned')
                    outputAliTs.copyInfo(inputTS)
                    outputAliTs.setSamplingRate(newPixSize)
                    outputAliTs.setStreamState(Set.STREAM_OPEN)
                    outputAliTs.write()
                    self._defineOutputs(**{OUTPUT_ALI_TILTSERIES_NAME: outputAliTs})
                    self._defineSourceRelation(self.getInputTsSet(pointer=True),
                                           outputAliTs)
                # else:
                #     outputAliTs.enableAppend()

                newTs = TiltSeries()
                newTs.copyInfo(ts)
                newTs.setInterpolated(True)
                outputAliTs.append(newTs)

                tltList = self.getNewTiltAngles(tsId)

                # Getting the excluded views in order to disable in the
                # aligned tilt-series
                alignFn = self._getExtraPath(tsId, 'align.com')
                excludedViewList = self.getExcludedViewList(alignFn,
                                                            reservedWord='ExcludeList')
                self.debug(f"Excluding views for align: {excludedViewList}")

                for tiltImage in ts.iterItems(iterate=False):
                    newTi = tiltImage.clone()
                    newTi.copyInfo(tiltImage, copyId=True, copyTM=False)

                    acq = tiltImage.getAcquisition()
                    acq.setTiltAxisAngle(0.)
                    newTi.setAcquisition(acq)

                    newTi.setOddEven([])
                    sliceIndex = newTi.getIndex()
                    newTi.setLocation(sliceIndex, aligFilePath)
                    if tltList is not None:
                        newTi.setTiltAngle(float(tltList[sliceIndex - 1]))
                    if sliceIndex in excludedViewList:
                        newTi.setEnabled(False)
                    newTs.append(newTi)

                acq = newTs.getAcquisition()
                acq.setTiltAxisAngle(0.)  # 0 because TS is aligned
                newTs.setAcquisition(acq)
                newTs.setDim(aliDims)

                outputAliTs.update(newTs)
                outputAliTs.write()

            """Original tilt-series with alignment information"""
            inFilePath = self.getExtraOutFile(tsId, ext=MRC_EXT)
            tmFilePath = self.getExtraOutFile(tsId, suffix="fid", ext=XF_EXT)

            # input TS were renamed by etomo when "using fixed stack"
            if os.path.exists(inFilePath) and os.path.exists(tmFilePath):

                if outputTs is None:
                    outputTs = self._createSetOfTiltSeries(suffix='_orig')
                    outputTs.copyInfo(inputTS)
                    outputTs.setStreamState(Set.STREAM_OPEN)
                    self._defineOutputs(**{OUTPUT_TILTSERIES_NAME: outputTs})
                    self._defineSourceRelation(self.getInputTsSet(pointer=True),
                                               outputTs)
                # else:
                #     outputTs.enableAppend()

                newTs = TiltSeries()
                newTs.copyInfo(ts)
                outputTs.append(newTs)

                tltList = self.getNewTiltAngles(tsId)

                # Getting the excluded views in order to disable in the
                # aligned tilt-series
                alignFn = self._getExtraPath(tsId, 'align.com')
                excludedViewList = self.getExcludedViewList(alignFn,
                                                            reservedWord='ExcludeList')
                self.debug(f"Excluding views for align: {excludedViewList}")

                newTransformationMatricesList = readXfFile(tmFilePath)

                if self.applyAlignment:
                    # input TS were interpolated during the convertInputStep,
                    # etomo was run on interpolated TS. So we have to save both alignments
                    doMatrixMultiplication = True
                else:
                    # input TS might have alignment, but it was ignored.
                    # we have to save only etomo alignment, replacing any previous one
                    doMatrixMultiplication = False

                for index, tiltImage in enumerate(ts):
                    newTi = tiltImage.clone()
                    newTi.copyInfo(tiltImage, copyId=True, copyTM=False)
                    acq = tiltImage.getAcquisition()
                    newTi.setAcquisition(acq)
                    sliceIndex = newTi.getIndex()
                    newTi.setLocation(sliceIndex, inFilePath)
                    if tltList is not None:
                        newTi.setTiltAngle(float(tltList[sliceIndex - 1]))
                    if sliceIndex in excludedViewList:
                        newTi.setEnabled(False)

                    transform = Transform()
                    newTransform = newTransformationMatricesList[:, :, index]
                    newTransformArray = np.array(newTransform)

                    if doMatrixMultiplication:
                        previousTransform = tiltImage.getTransform().getMatrix()
                        previousTransformArray = np.array(previousTransform)
                        outputTransformMatrix = np.matmul(newTransformArray, previousTransformArray)
                        transform.setMatrix(outputTransformMatrix)
                    else:
                        transform.setMatrix(newTransformArray)

                    newTi.setTransform(transform)
                    newTs.append(newTi)

                outputTs.update(newTs)
                outputTs.write()

            """Output set of coordinates 3D (associated to the aligned tilt-series)"""
            coordFilePath = self.getExtraOutFile(tsId, suffix='fid', ext=XYZ_EXT)
            coordFilePath = coordFilePath.replace("_fid", "fid")  # due to etomo bug

            if os.path.exists(coordFilePath) and outputAliTs is not None:
                if outputTSCoords is None:
                    outputTSCoords = tomoObj.SetOfTiltSeriesCoordinates.create(self._getPath(),
                                                                               suffix='Fiducials3D')
                    outputTSCoords.setSetOfTiltSeries(outputAliTs)
                    outputTSCoords.setStreamState(Set.STREAM_OPEN)
                    self._defineOutputs(**{OUTPUT_TS_COORDINATES_NAME: outputTSCoords})
                    self._defineSourceRelation(self.getInputTsSet(pointer=True),
                                               outputTSCoords)
                # else:
                #     outputTSCoords.enableAppend()

                coordList, _, _ = self.format3DCoordinatesList(coordFilePath)

                for element in coordList:
                    newCoord3D = tomoObj.TiltSeriesCoordinate()
                    newCoord3D.setTsId(tsId)
                    self.debug(f"Setting tilt series coordinate x, y, z: {element}")
                    newCoord3D.setX(element[0])
                    newCoord3D.setY(element[1])
                    newCoord3D.setZ(element[2])

                    outputTSCoords.append(newCoord3D)
                    outputTSCoords.update(newCoord3D)

                outputTSCoords.write()

            """Landmark models with no gaps"""
            modelFilePath = self.getExtraOutFile(tsId, suffix="nogaps", ext=FID_EXT)
            residFilePath = self.getExtraOutFile(tsId, ext=RESID_EXT)

            if os.path.exists(modelFilePath) and os.path.exists(residFilePath):
                modelFilePathTxt = self.getExtraOutFile(tsId, suffix="nogaps_fid",
                                                        ext=TXT_EXT)

                paramsModel2Point = {
                    '-InputFile': modelFilePath,
                    '-OutputFile': modelFilePathTxt
                }
                self.runProgram('model2point', paramsModel2Point)

                outputSetOfLandmarkModelsNoGaps = self.getOutputFiducialModel(outputPreAliTs, forceNew=outputSetOfLandmarkModelsNoGaps is None)

                fiducialNoGapList = fiducialModel2List(modelFilePathTxt)

                landmarkModelNoGapsFilePath = self.getExtraOutFile(tsId, suffix="nogaps",
                                                                   ext=SFID_EXT)
                fiducialNoGapsResidList = fidResidualModel2List(residFilePath)
                landmarkModelNoGaps = tomoObj.LandmarkModel(tsId=tsId,
                                                            fileName=landmarkModelNoGapsFilePath,
                                                            modelName=modelFilePath,
                                                            size=self.markersDiameter.get() * 10 / ts.getSamplingRate())

                prevTiltIm = 0
                chainId = 0
                indexFake = 0
                for fiducial in fiducialNoGapList:
                    if int(float(fiducial[2])) <= prevTiltIm:
                        chainId += 1
                    prevTiltIm = int(float(fiducial[2]))
                    if (indexFake < len(fiducialNoGapsResidList) and
                            fiducial[2] == fiducialNoGapsResidList[indexFake][2]):
                        landmarkModelNoGaps.addLandmark(xCoor=fiducial[0],
                                                        yCoor=fiducial[1],
                                                        tiltIm=fiducial[2],
                                                        chainId=chainId,
                                                        xResid=fiducialNoGapsResidList[indexFake][3],
                                                        yResid=fiducialNoGapsResidList[indexFake][4])
                        indexFake += 1
                    else:
                        landmarkModelNoGaps.addLandmark(xCoor=fiducial[0],
                                                        yCoor=fiducial[1],
                                                        tiltIm=fiducial[2],
                                                        chainId=chainId,
                                                        xResid='0',
                                                        yResid='0')

                outputSetOfLandmarkModelsNoGaps.append(landmarkModelNoGaps)
                outputSetOfLandmarkModelsNoGaps.write()

            """Full reconstructed tomogram"""
            reconstructTomoFilePath = self.getExtraOutFile(tsId, suffix="full_rec",
                                                           ext=MRC_EXT)
            if os.path.exists(reconstructTomoFilePath):
                tomoDims, newPixSize = self.getNewPixAndDim(reconstructTomoFilePath)
                self.debug(f"{reconstructTomoFilePath}: pix = {newPixSize}, dims = {tomoDims}")
                if outputTomos is None:
                    outputTomos = self._createSetOfTomograms(suffix='_raw')
                    outputTomos.copyInfo(inputTS)
                    outputTomos.setSamplingRate(newPixSize)
                    outputTomos.setStreamState(Set.STREAM_OPEN)
                    self._defineOutputs(FullTomograms=outputTomos)
                    self._defineSourceRelation(self.getInputTsSet(pointer=True),
                                               outputTomos)
                # else:
                #     outputTomos.enableAppend()

                newTomogram = tomoObj.Tomogram()
                newTomogram.setLocation(reconstructTomoFilePath)
                newTomogram.setTsId(tsId)
                newTomogram.setSamplingRate(newPixSize)
                # Set default tomogram origin
                newTomogram.setOrigin(newOrigin=None)

                outputTomos.append(newTomogram)
                outputTomos.write()

            """Post-processed reconstructed tomogram"""
            posprocessedRecTomoFilePath = self.getExtraOutFile(tsId, suffix="rec",
                                                               ext=MRC_EXT)
            if os.path.exists(posprocessedRecTomoFilePath):
                tomoDims, newPixSize = self.getNewPixAndDim(posprocessedRecTomoFilePath)
                self.debug(f"{posprocessedRecTomoFilePath}: pix = {newPixSize}, dims = {tomoDims}")
                if outputPostProcessTomos is None:
                    outputPostProcessTomos = self._createSetOfTomograms()
                    outputPostProcessTomos.copyInfo(inputTS)
                    outputPostProcessTomos.setSamplingRate(newPixSize)
                    outputPostProcessTomos.setStreamState(Set.STREAM_OPEN)
                    self._defineOutputs(PostProcessTomograms=outputPostProcessTomos)
                    self._defineSourceRelation(self.getInputTsSet(pointer=True),
                                               outputPostProcessTomos)
                # else:
                #     outputPostProcessTomos.enableAppend()

                newTomogram = tomoObj.Tomogram()
                newTomogram.setLocation(posprocessedRecTomoFilePath)
                newTomogram.setTsId(tsId)
                newTomogram.setSamplingRate(newPixSize)

                # Set default tomogram origin
                newTomogram.setOrigin(newOrigin=None)

                outputPostProcessTomos.append(newTomogram)
                outputPostProcessTomos.write()

        self.inputTiltSeries = None  # This temp variable is persisted,even keeping the ID received when saving it as part of the set.
        self._closeOutputSet()

    # --------------------------- UTILS functions -----------------------------
    @staticmethod
    def _writeEtomoEdf(fn, paramsDict):
        template = """
#%(date)s  // Generated by Scipion, version: %(version)s
Setup.DataSource=CCD
ReconstructionState.InvalidEdgeFunctionsA=no result
Setup.BackupDirectory=
ReconstructionState.InvalidEdgeFunctionsB=no result
ProcessTrack.PostProcessing=Not started
Setup.Combine.ManualCleanup=false
Setup.AxisA.TiltAngle.Type=File
Setup.FiducialessAlignmentB=false
Setup.FiducialessAlignmentA=false
ProcessTrack.FinalAlignedStack-A=Not started
ProcessTrack.FinalAlignedStack-B=Not started
Setup.tiltalign.TargetPatchSizeXandY=700,700
Setup.AxisB.TiltAngle.RangeMin=%(minTilt)f
Setup.Combine.TempDirectory=
Setup.Combine.UseList=
Setup.Combine.PatchBoundaryYMax=0
Setup.WholeTomogramSampleA=true
Setup.DatasetName=%(name)s
Setup.FiducialDiameter=%(markerDiameter)f
Setup.WholeTomogramSampleB=true
Setup.SetFEIPixelSize=true
ReconstructionState.A.AdjustOrigin=true
Setup.Setup.OrigImageStackExt=.st
ProcessTrack.RevisionNumber=2.0
Setup.ProjectLog.FrameSize.Width=683
Setup.Combine.PatchBoundaryYMin=0
Setup.Combine.MaxPatchBoundaryZMax=0
ProcessTrack.CoarseAlignment-A=Not started
Setup.ViewType=Single View
ProcessTrack.CoarseAlignment-B=Not started
ReconstructionState.TrimvolFlipped=no result
Setup.Track.B.TrackMethod=Seed
Setup.Track.A.TrackMethod=Seed
Setup.A.SizeToOutputInXandY=/
ProcessTrack.TomogramGeneration-A=Not started
ProcessTrack.TomogramGeneration-B=Not started
ProcessTrack.CleanUp=Not started
Setup.AxisA.TiltAngle.RangeMin=%(minTilt)f
Setup.AxisB.TiltAngle.TiltAngleFilename=
Setup.Pos.A.NewDialog=true
Setup.Squeezevol.LinearInterpolation=false
Setup.Stack.B.Is.Twodir=false
Setup.Combine.RevisionNumber=1.2
Setup.Stack.B.CTF.AutoFit.RangeAndStep=-Infinity,-Infinity
Setup.UseLocalAlignmentsB=true
Setup.UseLocalAlignmentsA=true
Setup.AxisB.TiltAngle.RangeStep=1.0
Setup.Combine.FiducialMatchListA=
Setup.Combine.FiducialMatchListB=
Setup.B.SizeToOutputInXandY=/
Setup.Binning=1
Setup.Combine.ModelBased=false
ReconstructionState.MadeZFactorsB=no result
ReconstructionState.MadeZFactorsA=no result
ReconstructionState.SqueezevolFlipped=no result
ProcessTrack.FiducialModel-B=Not started
ProcessTrack.FiducialModel-A=Not started
Setup.Combine.PatchBoundaryXMax=0
ProcessTrack.Setup=Complete
ProcessTrack.PreProcessing-A=Not started
ProcessTrack.PreProcessing-B=Not started
ProcessTrack.FineAlignment-B=Not started
ProcessTrack.FineAlignment-A=Not started
Setup.AxisA.TiltAngle.TiltAngleFilename=
Setup.AxisA.TiltAngle.RangeStep=1.0
Setup=-Infinity,-Infinity
Setup.Combine.PatchBoundaryXMin=0
Setup.MagGradientFile=
Setup.RevisionNumber=1.12
Setup.Track.B.SeedModel.Transfer=true
Setup.Track.A.Raptor.UseRawStack=false
Setup.Orig.SystemTemplate=%(imodDir)s/SystemTemplate/cryoSample.adoc
Setup.PixelSize=%(pixelSize)f
ReconstructionState.B.AdjustOrigin=true
ReconstructionState.Combine.ScriptsCreated=no result
Setup.Combine.Transfer=true
Setup.DistortionFile=
ReconstructionState.UsedLocalAlignmentsA=no result
ReconstructionState.UsedLocalAlignmentsB=no result
Setup.ProjectLog.FrameSize.Height=230
Setup.Stack.B.Twodir=0.0
Setup.Combine.FiducialMatch=BothSides
Setup.AxisType=Single Axis
Setup.ImageRotationA=%(rotationAngle)f
ProcessTrack.TomogramPositioning-A=Not started
Setup.ImageRotationB=
Setup.Stack.A.Is.Twodir=false
Setup.Pos.B.NewDialog=true
ProcessTrack.TomogramPositioning-B=Not started
Setup.Combine.PatchBoundaryZMax=0
Setup.DefaultParallel=%(useCpu)s
Setup.DefaultGpuProcessing=false
Setup.Track.A.SeedModel.Auto=true
Setup.Combine.PatchSize=M
Setup.AxisB.TiltAngle.Type=Extract
Setup.Combine.PatchBoundaryZMin=0
Setup.ProjectLog.FrameLocation.Y=55
Setup.ProjectLog.FrameLocation.X=95
Setup.ProjectLog.Visible=true
Setup.tiltalign.NumberOfLocalPatchesXandY=5,5
Setup.Combine.PatchRegionModel=
ReconstructionState.NewstFiducialessAlignmentA=no result
ReconstructionState.NewstFiducialessAlignmentB=no result
Setup.Stack.A.Twodir=0.0
ProcessTrack.TomogramCombination=Not started
        """
        with open(fn, 'w') as f:
            f.write(template % paramsDict)

    def getNewPixAndDim(self, fn):
        dims = ih.getDimensions(fn)
        dims = dims[:-1]
        origDimX, origDimY = self._getInputDims()
        originalDim = max(origDimX, origDimY)
        outputDim = max(dim for dim in dims[:2])
        newPixSize = self.inputTiltSeries.getSamplingRate() * round(originalDim / outputDim)

        return dims, newPixSize

    def _getInputDims(self):
        """ Return XY size of the input TS. """
        x, y, _, _ = ih.getDimensions(self.inputTiltSeries.getFirstItem().getFileName())
        return x, y

    @staticmethod
    def getExcludedViewList(fn, reservedWord="ExcludeList"):
        excludedViewList = []
        with open(fn) as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith(reservedWord):
                    excludedStr = line.strip().split("\t")[1]
                    excludedViewList = pwutils.getListFromRangeString(excludedStr)
                    break

        return excludedViewList

    def getNewTiltAngles(self, tsId):
        tltList = None
        tltFilePath = self.getExtraOutFile(tsId, suffix='fid', ext=TLT_EXT)
        if os.path.exists(tltFilePath):
            tltList = utils.formatAngleList(tltFilePath)
            self.debug(f"{tltFilePath} read: {tltList}")

        return tltList

    @staticmethod
    def format3DCoordinatesList(coordFilePath):
        """ This method takes an IMOD-based fiducial coordinates
        file path and returns a list containing each coordinate
        for each fiducial belonging to the tilt-series. """

        coorList = []

        with open(coordFilePath) as f:
            coorText = f.read().splitlines()

            for i, line in enumerate(coorText):
                if line != '':

                    logger.debug("Fiducial coordinate line is: %s" % line)
                    vector = line.replace('-', ' -').split()

                    logger.debug("Fiducial vector is: %s" % vector)
                    if i == 0:
                        xDim = int(vector[-2])
                        yDim = int(vector[-1])

                    coorList.append([float(vector[1]), float(vector[2]), float(vector[3])])

        return coorList, xDim, yDim

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = ["The following outputs have been generated from the "
                   "operations performed on the input tilt-series:"]

        preAliTs = getattr(self, OUTPUT_PREALI_TILTSERIES_NAME, None)
        aliTs = getattr(self, OUTPUT_ALI_TILTSERIES_NAME, None)
        outputTsWithAli = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        fidModelNoGaps = getattr(self, OUTPUT_FIDUCIAL_NO_GAPS_NAME, None)
        fidModelGaps = getattr(self, OUTPUT_FIDUCIAL_GAPS_NAME, None)
        tsCoords = getattr(self, OUTPUT_TS_COORDINATES_NAME, None)
        tomos = getattr(self, "FullTomograms", None)
        postProcessTomos = getattr(self, "PostProcessTomograms", None)

        if preAliTs is not None:
            summary.append("- Pre-aligned tilt-series")

        if aliTs is not None:
            summary.append("- Aligned tilt-series")

        if outputTsWithAli is not None:
            summary.append("- Tilt-series with alignment")

        if tsCoords is not None:
            summary.append("- Landmark 3D coordinates")

        if fidModelGaps is not None:
            summary.append("- Landmark model with gaps")

        if fidModelNoGaps is not None:
            summary.append("- Landmark model without gaps")

        if tomos is not None:
            summary.append("- Raw reconstructed tomogram")

        if postProcessTomos is not None:
            summary.append("- Post-processed tomogram")

        if len(summary) == 1:
            summary = ["Outputs are not ready yet"]

        return summary
