# *****************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Centro Nacional de Biotecnologia, CSIC, Spain
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
import numpy as np

import pyworkflow as pw
import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import STEPS_SERIAL
import pyworkflow.utils as pwutils
from pwem.emlib.image import ImageHandler as ih
from pwem.objects import Transform
import tomo.objects as tomoObj

from imod import Plugin, utils
from imod.protocols import ProtImodBase
from imod.constants import (OUTPUT_PREALI_TILTSERIES_NAME, OUTPUT_ALI_TILTSERIES_NAME,
                            OUTPUT_TILTSERIES_NAME, OUTPUT_TS_COORDINATES_NAME,
                            OUTPUT_FIDUCIAL_NO_GAPS_NAME, RAWTLT_EXT, TLT_EXT, XF_EXT,
                            EDF_EXT, MRC_EXT, SFID_EXT, XYZ_EXT, FID_EXT,
                            RESID_EXT, TXT_EXT)


class ProtImodEtomo(ProtImodBase):
    """
    Simple wrapper around etomo to manually reconstruct a Tomogram.

    More info:
        https://bio3d.colorado.edu/imod/doc/etomoTutorial.html

    Etomo is software tool for assisting users in the tomographic reconstruction
    process of both single and dual axis tilt series. Throughout this procedure,
    eTomo executes numerous program commands and frequently initiates 3dmod
    and Midas to enable users to make precise adjustments. Some of the main features
    are:\n
    - Xray eraser\n
    - dose filtering\n
    - Tilt series alignment\n
    - Gold beads detection and eraser\n
    - Tomogram reconstruction\n
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

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.stepsExecutionMode = STEPS_SERIAL
        self.PrealignedTiltSeries = None
        self.AlignedTiltSeries = None
        self.FullTomograms = None
        self.PostProcessTomograms = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
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

        form.addParallelSection(threads=4, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    # Overwrite the following function to prevent streaming from base class
    def _stepsCheck(self):
        pass

    def _insertAllSteps(self):
        self.inputTiltSeries = None
        self._insertFunctionStep(self.runEtomoStep, interactive=True)

    # --------------------------- STEPS functions -----------------------------
    def runEtomoStep(self):
        from imod.viewers import ImodGenericView
        setOftiltSeries = self.getInputSet()
        view = ImodGenericView(None, self, setOftiltSeries,
                               isInteractive=True)
        view.show()
        self.createOutput()

    def runAllSteps(self, obj):
        for item in self.getInputSet():  # FIXME: why?
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
            '-name': pwutils.removeBaseExt(inputTsFileName),
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
                                'name': pwutils.removeBaseExt(inputTsFileName),
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
        outputSetOfTiltSeries = None  # original TS with new alignment
        outputPrealiSetOfTiltSeries = None
        outputAliSetOfTiltSeries = None
        setOfTSCoords = None
        outputSetOfFullTomograms = None
        outputSetOfPostProcessTomograms = None
        inputTS = self.getInputSet()

        for ts in inputTS:
            self.inputTiltSeries = ts
            tsId = ts.getTsId()

            """Prealigned tilt-series"""
            prealiFilePath = self.getExtraOutFile(tsId, suffix="preali", ext=MRC_EXT)
            if os.path.exists(prealiFilePath):
                xPrealiDims, newPixSize = self.getNewPixAndDim(prealiFilePath)
                self.debug(f"{prealiFilePath}: pix = {newPixSize}, dims = {xPrealiDims}")
                if outputPrealiSetOfTiltSeries is None:
                    outputPrealiSetOfTiltSeries = self._createSetOfTiltSeries(suffix='_prealigned')
                    outputPrealiSetOfTiltSeries.copyInfo(inputTS)
                    outputPrealiSetOfTiltSeries.setSamplingRate(newPixSize)
                    self._defineOutputs(**{OUTPUT_PREALI_TILTSERIES_NAME: outputPrealiSetOfTiltSeries})
                    self._defineSourceRelation(self.getInputSet(pointer=True),
                                               outputPrealiSetOfTiltSeries)
                else:
                    outputPrealiSetOfTiltSeries.enableAppend()

                newTs = ts.clone()
                newTs.copyInfo(ts)
                newTs.setInterpolated(True)
                outputPrealiSetOfTiltSeries.append(newTs)

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

                outputPrealiSetOfTiltSeries.update(newTs)

            """Aligned tilt-series"""
            aligFilePath = self.getExtraOutFile(tsId, suffix="ali", ext=MRC_EXT)
            if os.path.exists(aligFilePath):
                aliDims, newPixSize = self.getNewPixAndDim(aligFilePath)
                self.debug(f"{aligFilePath}: pix = {newPixSize}, dims = {aliDims}")

                if outputAliSetOfTiltSeries is None:
                    outputAliSetOfTiltSeries = self._createSetOfTiltSeries(suffix='_aligned')
                    outputAliSetOfTiltSeries.copyInfo(inputTS)
                    outputAliSetOfTiltSeries.setSamplingRate(newPixSize)
                    self._defineOutputs(**{OUTPUT_ALI_TILTSERIES_NAME: outputAliSetOfTiltSeries})
                    self._defineSourceRelation(self.getInputSet(pointer=True),
                                               outputAliSetOfTiltSeries)
                else:
                    outputAliSetOfTiltSeries.enableAppend()

                newTs = ts.clone()
                newTs.copyInfo(ts)
                newTs.setInterpolated(True)
                outputAliSetOfTiltSeries.append(newTs)

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

                outputAliSetOfTiltSeries.update(newTs)

            """Original tilt-series with alignment information"""
            inFilePath = self.getExtraOutFile(tsId, suffix="orig", ext=MRC_EXT)
            # input TS were renamed by etomo when "using fixed stack"
            if os.path.exists(inFilePath):

                if outputSetOfTiltSeries is None:
                    outputSetOfTiltSeries = self._createSetOfTiltSeries(suffix='_orig')
                    outputSetOfTiltSeries.copyInfo(inputTS)
                    self._defineOutputs(**{OUTPUT_TILTSERIES_NAME: outputSetOfTiltSeries})
                    self._defineSourceRelation(self.getInputSet(pointer=True),
                                               outputSetOfTiltSeries)
                else:
                    outputSetOfTiltSeries.enableAppend()

                newTs = ts.clone()
                newTs.copyInfo(ts)
                outputSetOfTiltSeries.append(newTs)

                tltList = self.getNewTiltAngles(tsId)

                # Getting the excluded views in order to disable in the
                # aligned tilt-series
                alignFn = self._getExtraPath(tsId, 'align.com')
                excludedViewList = self.getExcludedViewList(alignFn,
                                                            reservedWord='ExcludeList')
                self.debug(f"Excluding views for align: {excludedViewList}")

                tmFilePath = self.getExtraOutFile(tsId, suffix="fid", ext=XF_EXT)
                newTransformationMatricesList = utils.formatTransformationMatrix(tmFilePath)

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

                    if doMatrixMultiplication:
                        previousTransform = tiltImage.getTransform().getMatrix()
                        newTransform = newTransformationMatricesList[:, :, index]
                        previousTransformArray = np.array(previousTransform)
                        newTransformArray = np.array(newTransform)
                        outputTransformMatrix = np.matmul(newTransformArray, previousTransformArray)
                        transform.setMatrix(outputTransformMatrix)
                        newTi.setTransform(transform)
                    else:
                        newTransform = newTransformationMatricesList[:, :, index]
                        newTransformArray = np.array(newTransform)
                        transform.setMatrix(newTransformArray)
                        newTi.setTransform(transform)

                    newTs.append(newTi)

                #acq = newTs.getAcquisition()
                #newTs.setAcquisition(acq)
                newTs.write(properties=False)

                outputSetOfTiltSeries.update(newTs)
                outputSetOfTiltSeries.write()
                self._store(outputSetOfTiltSeries)

            """Output set of coordinates 3D (associated to the aligned tilt-series)"""
            coordFilePath = self.getExtraOutFile(tsId, suffix='fid', ext=XYZ_EXT)
            coordFilePath = coordFilePath.replace("_fid", "fid")  # due to etomo bug

            if os.path.exists(coordFilePath) and outputAliSetOfTiltSeries is not None:
                if setOfTSCoords is None:
                    setOfTSCoords = tomoObj.SetOfTiltSeriesCoordinates.create(self._getPath(),
                                                                              suffix='Fiducials3D')
                    setOfTSCoords.setSetOfTiltSeries(outputAliSetOfTiltSeries)
                    self._defineOutputs(**{OUTPUT_TS_COORDINATES_NAME: setOfTSCoords})
                    self._defineSourceRelation(self.getInputSet(pointer=True),
                                               setOfTSCoords)
                else:
                    setOfTSCoords.enableAppend()

                coordList, _, _ = utils.format3DCoordinatesList(coordFilePath)

                for element in coordList:
                    newCoord3D = tomoObj.TiltSeriesCoordinate()
                    newCoord3D.setTsId(tsId)
                    self.debug(f"Setting tilt series coordinate x, y, z: {element}")
                    newCoord3D.setX(element[0])
                    newCoord3D.setY(element[1])
                    newCoord3D.setZ(element[2])

                    setOfTSCoords.append(newCoord3D)
                    setOfTSCoords.update(newCoord3D)

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

                outputSetOfLandmarkModelsNoGaps = self.getOutputFiducialModelNoGaps(outputPrealiSetOfTiltSeries)

                fiducialNoGapList = utils.formatFiducialList(modelFilePathTxt)

                landmarkModelNoGapsFilePath = self.getExtraOutFile(tsId, suffix="nogaps",
                                                                   ext=SFID_EXT)
                fiducialNoGapsResidList = utils.formatFiducialResidList(residFilePath)
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
                outputSetOfLandmarkModelsNoGaps.update(landmarkModelNoGaps)

            """Full reconstructed tomogram"""
            reconstructTomoFilePath = self.getExtraOutFile(tsId, suffix="full_rec",
                                                           ext=MRC_EXT)
            if os.path.exists(reconstructTomoFilePath):
                tomoDims, newPixSize = self.getNewPixAndDim(reconstructTomoFilePath)
                self.debug(f"{reconstructTomoFilePath}: pix = {newPixSize}, dims = {tomoDims}")
                if outputSetOfFullTomograms is None:
                    outputSetOfFullTomograms = self._createSetOfTomograms(suffix='_raw')
                    outputSetOfFullTomograms.copyInfo(inputTS)
                    outputSetOfFullTomograms.setSamplingRate(newPixSize)
                    self._defineOutputs(FullTomograms=outputSetOfFullTomograms)
                    self._defineSourceRelation(self.getInputSet(pointer=True),
                                               outputSetOfFullTomograms)
                else:
                    outputSetOfFullTomograms.enableAppend()

                newTomogram = tomoObj.Tomogram()
                newTomogram.setLocation(reconstructTomoFilePath)
                newTomogram.setTsId(tsId)
                newTomogram.setSamplingRate(newPixSize)
                # Set default tomogram origin
                newTomogram.setOrigin(newOrigin=None)

                outputSetOfFullTomograms.append(newTomogram)
                outputSetOfFullTomograms.update(newTomogram)

            """Post-processed reconstructed tomogram"""
            posprocessedRecTomoFilePath = self.getExtraOutFile(tsId, suffix="rec",
                                                               ext=MRC_EXT)
            if os.path.exists(posprocessedRecTomoFilePath):
                tomoDims, newPixSize = self.getNewPixAndDim(posprocessedRecTomoFilePath)
                self.debug(f"{posprocessedRecTomoFilePath}: pix = {newPixSize}, dims = {tomoDims}")
                if outputSetOfPostProcessTomograms is None:
                    outputSetOfPostProcessTomograms = self._createSetOfTomograms()
                    outputSetOfPostProcessTomograms.copyInfo(inputTS)
                    outputSetOfPostProcessTomograms.setSamplingRate(newPixSize)
                    self._defineOutputs(PostProcessTomograms=outputSetOfPostProcessTomograms)
                    self._defineSourceRelation(self.getInputSet(pointer=True),
                                               outputSetOfPostProcessTomograms)
                else:
                    outputSetOfPostProcessTomograms.enableAppend()

                newTomogram = tomoObj.Tomogram()
                newTomogram.setLocation(posprocessedRecTomoFilePath)
                newTomogram.setTsId(tsId)
                newTomogram.setSamplingRate(newPixSize)

                # Set default tomogram origin
                newTomogram.setOrigin(newOrigin=None)

                outputSetOfPostProcessTomograms.append(newTomogram)
                outputSetOfPostProcessTomograms.update(newTomogram)

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

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = ["The following outputs have been generated from the "
                   "operations performed on the input tilt-series:"]

        if self.PreAlignedTiltSeries:
            summary.append("- Pre-aligned tilt-series")

        if self.AlignedTiltSeries:
            summary.append("- Aligned tilt-series")

        if self.TiltSeriesCoordinates:
            summary.append("- Landmark 3D coordinates")

        if self.FiducialModelGaps:
            summary.append("- Landmark model with gaps")

        if self.FiducialModelNoGaps:
            summary.append("- Landmark model without gaps")

        if self.FullTomograms:
            summary.append("- Raw reconstructed tomogram")

        if self.PostProcessTomograms:
            summary.append("- Post-processed tomogram")

        if len(summary) == 1:
            summary = ["Outputs are not ready yet"]

        return summary
