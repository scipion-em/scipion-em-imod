# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Centro Nacional de Biotecnologia, CSIC, Spain
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

import pyworkflow as pw
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
from tomo.convert import writeTiStack
from imod import Plugin
from imod import utils
from pwem.emlib.image import ImageHandler


class ProtImodEtomo(EMProtocol, ProtTomoBase):
    """
    Simple wrapper around etomo to manually reconstruct a Tomogram.

    More info:
        https://bio3d.colorado.edu/imod/doc/etomoTutorial.html
    """

    _label = 'etomo interactive'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputTiltSeries',
                      params.PointerParam,
                      pointerClass='TiltSeries',
                      important=True,
                      label='Input Tilt-Series',
                      help='Input tilt-series to be processed with etomo.')

        form.addParam('excludeList',
                      params.StringParam,
                      default='',
                      label='Exclusion list',
                      help='Provide tilt images IDs (usually starting at 1) that you want to exclude from the '
                           'processing separated by blank spaces.')

        form.addParam('binning',
                      params.IntParam,
                      default=2,
                      label='Bin the input images',
                      help='Binning of the input images.')

        form.addParam('markersDiameter',
                      params.FloatParam,
                      default=10,
                      label='Fiducial markers diameter (nm)',
                      help='Diameter of gold beads in nanometers.')

        form.addParam('rotationAngle',
                      params.FloatParam,
                      default=0.0,
                      label='Tilt rotation angle in degrees',
                      help='Angle from the vertical to the tilt axis in raw images.')

    # -------------------------- INSERT steps functions ---------------------
    # Overwrite the following function to prevent streaming from base class
    def _stepsCheck(self):
        pass

    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runEtomoStep', interactive=True)

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self):
        ts = self.inputTiltSeries.get()
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        outputTsFileName = os.path.join(extraPrefix, ts.getFirstItem().parseFileName())

        """Apply transformation matrices and remove excluded views"""
        if self.excludeList.get() == '':

            """Apply the transformation form the input tilt-series"""
            ts.applyTransform(outputTsFileName)

            """Generate angle file"""
            angleFilePath = os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".rawtlt"))
            ts.generateTltFile(angleFilePath)

        else:
            interpolatedTsFileName = os.path.join(tmpPrefix, ts.getFirstItem().parseFileName())
            angleFilePath = os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".rawtlt"))

            """Apply the transformation form the input tilt-series and generate a new ts object"""
            ts.applyTransform(interpolatedTsFileName)

            interpolatedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Interpolated')
            interpolatedSetOfTiltSeries.copyInfo(self.inputTiltSeries.get())
            interpolatedSetOfTiltSeries.setDim(self.inputTiltSeries.get().getDim())

            interpolatedTs = tomoObj.TiltSeries(tsId=tsId)
            interpolatedTs.copyInfo(ts)

            interpolatedSetOfTiltSeries.append(interpolatedTs)

            for index, tiltImage in enumerate(ts):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)
                newTi.setLocation(index + 1, interpolatedTsFileName)
                interpolatedTs.append(newTi)
            interpolatedTs.write()

            """Write a new stack discarding excluded tilts"""
            excludeList = [int(i) for i in self.excludeList.get().split()]
            tiList = [ti.clone() for ti in interpolatedTs]
            tiList.sort(key=lambda ti: ti.getTiltAngle())

            writeTiStack(tiList,
                         outputStackFn=outputTsFileName,
                         outputTltFn=angleFilePath,
                         excludeList=excludeList)

        """Generate etomo config file"""
        args = '-name %s ' % ts.getFirstItem().parseFileName(extension="")
        args += '-gold %0.3f ' % self.markersDiameter

        # Imod use the pixel size in NM
        pixelSizeNm = ts.getSamplingRate() / 10.

        args += '-pixel %0.3f ' % pixelSizeNm
        args += '-binning %d ' % self.binning
        args += '-rotation %0.3f ' % self.rotationAngle
        args += '-userawtlt'

        Plugin.runImod(self, 'copytomocoms', args, cwd=extraPrefix)

        edfFn = os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".edf"))
        minTilt = min(utils.formatAngleList(os.path.join(extraPrefix,
                                                         ts.getFirstItem().parseFileName(extension=".rawtlt"))))
        self._writeEtomoEdf(edfFn,
                            {
                                'date': pw.utils.prettyTime(),
                                'name': ts.getFirstItem().parseFileName(extension=''),
                                'pixelSize': pixelSizeNm,
                                'version': pw.__version__,
                                'minTilt': minTilt,
                                'binning': self.binning,
                                'markerDiameter': self.markersDiameter,
                                'rotationAngle': self.rotationAngle
                            })

    def runEtomoStep(self):
        ts = self.inputTiltSeries.get()
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        args = '--fg '
        args += ts.getFirstItem().parseFileName(extension=".edf")
        Plugin.runImod(self, 'etomo', args, cwd=extraPrefix)
        self.createOutputStep()

    def createOutputStep(self):
        ts = self.inputTiltSeries.get()
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        """Prealigned tilt-series"""
        if os.path.exists(os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".preali"))):
            outputPrealiSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Preali')
            outputPrealiSetOfTiltSeries.copyInfo(self.inputTiltSeries.get())
            outputPrealiSetOfTiltSeries.setDim(self.inputTiltSeries.get().getDim())
            self._defineOutputs(outputPrealignedSetOfTiltSeries=outputPrealiSetOfTiltSeries)
            self._defineSourceRelation(self.inputTiltSeries, outputPrealiSetOfTiltSeries)

            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            outputPrealiSetOfTiltSeries.append(newTs)

            ih = ImageHandler()

            for index, tiltImage in enumerate(ts):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)
                newTi.setLocation(
                    index + 1,
                    (os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".preali")))
                )
                xPreali, _, _, _ = ih.getDimensions(newTi.getFileName()+":mrc")
                newTi.setSamplingRate(self.getPixSizeFromDimensions(xPreali))
                newTs.append(newTi)

            xPreali, yPreali, zPreali, _ = ih.getDimensions(newTs.getFirstItem().getFileName()+":mrc")
            newTs.setDim((xPreali, yPreali, zPreali))

            newTs.write(properties=False)

            outputPrealiSetOfTiltSeries.setSamplingRate(self.getPixSizeFromDimensions(xPreali))
            outputPrealiSetOfTiltSeries.update(newTs)
            outputPrealiSetOfTiltSeries.write()
            self._store()

        """Aligned tilt-series"""
        if os.path.exists(os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".ali"))):
            outputAliSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Ali')
            outputAliSetOfTiltSeries.copyInfo(self.inputTiltSeries.get())
            outputAliSetOfTiltSeries.setDim(self.inputTiltSeries.get().getDim())
            self._defineOutputs(outputAlignedSetOfTiltSeries=outputAliSetOfTiltSeries)
            self._defineSourceRelation(self.inputTiltSeries, outputAliSetOfTiltSeries)

            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            outputAliSetOfTiltSeries.append(newTs)

            ih = ImageHandler()

            tltFilePath = os.path.join(extraPrefix, ts.getFirstItem().parseFileName(suffix='_fid', extension=".tlt"))
            tltList = utils.formatAngleList(tltFilePath)

            for index, tiltImage in enumerate(ts):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)
                newTi.setLocation(
                    index + 1,
                    (os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".ali")))
                )
                newTi.setTiltAngle(float(tltList[index]))
                xAli, _, _, _ = ih.getDimensions(newTi.getFileName()+":mrc")
                newTi.setSamplingRate(self.getPixSizeFromDimensions(xAli))
                newTs.append(newTi)

            xAli, yAli, zAli, _ = ih.getDimensions(newTs.getFirstItem().getFileName() + ":mrc")
            newTs.setDim((xAli, yAli, zAli))

            newTs.write(properties=False)

            outputAliSetOfTiltSeries.setSamplingRate(self.getPixSizeFromDimensions(xAli))
            outputAliSetOfTiltSeries.update(newTs)
            outputAliSetOfTiltSeries.write()
            self._store()

            """Output set of coordinates 3D (associated to the aligned tilt-series)"""
            if os.path.exists(os.path.join(extraPrefix,
                                           ts.getFirstItem().parseFileName(suffix='fid', extension=".xyz"))):
                outputSetOfCoordinates3D = self._createSetOfCoordinates3D(volSet=outputAliSetOfTiltSeries,
                                                                          suffix='LandmarkModel')
                outputSetOfCoordinates3D.setSamplingRate(self.inputTiltSeries.get().getSamplingRate())
                outputSetOfCoordinates3D.setPrecedents(outputAliSetOfTiltSeries)
                self._defineOutputs(outputSetOfCoordinates3D=outputSetOfCoordinates3D)
                self._defineSourceRelation(self.inputTiltSeries, outputSetOfCoordinates3D)

                coordFilePath = os.path.join(extraPrefix,
                                             ts.getFirstItem().parseFileName(suffix='fid', extension=".xyz"))
                coordList = utils.format3DCoordinatesList(coordFilePath, xAli, yAli)
                for element in coordList:
                    newCoord3D = tomoObj.Coordinate3D(x=element[0],
                                                      y=element[1],
                                                      z=element[2])
                    newCoord3D.setVolume(ts)
                    newCoord3D.setVolId(ts.getObjId())
                    outputSetOfCoordinates3D.append(newCoord3D)
                    outputSetOfCoordinates3D.update(newCoord3D)
                outputSetOfCoordinates3D.write()
                self._store()

        """Landmark models with gaps"""
        if os.path.exists(os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".fid"))) and \
                os.path.exists(os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".resid"))):
            paramsGapPoint2Model = {
                'inputFile': os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".fid")),
                'outputFile': os.path.join(extraPrefix,
                                           ts.getFirstItem().parseFileName(suffix="_fid", extension=".txt"))
            }
            argsGapPoint2Model = "-InputFile %(inputFile)s " \
                                 "-OutputFile %(outputFile)s"
            Plugin.runImod(self, 'model2point', argsGapPoint2Model % paramsGapPoint2Model)

            outputSetOfLandmarkModelsGaps = self._createSetOfLandmarkModels(suffix='Gaps')
            outputSetOfLandmarkModelsGaps.copyInfo(self.inputTiltSeries.get())
            self._defineOutputs(outputSetOfLandmarkModelsGaps=outputSetOfLandmarkModelsGaps)
            self._defineSourceRelation(self.inputTiltSeries, outputSetOfLandmarkModelsGaps)

            fiducialModelGapPath = os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".fid"))

            landmarkModelGapsFilePath = os.path.join(extraPrefix,
                                                     ts.getFirstItem().parseFileName(suffix="_gaps", extension=".sfid"))

            landmarkModelGapsResidPath = os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".resid"))
            fiducialGapResidList = utils.formatFiducialResidList(landmarkModelGapsResidPath)

            landmarkModelGaps = tomoObj.LandmarkModel(tsId, landmarkModelGapsFilePath, fiducialModelGapPath)

            prevTiltIm = 0
            chainId = 0
            for index, fiducial in enumerate(fiducialGapResidList):
                if int(fiducial[2]) <= prevTiltIm:
                    chainId += 1
                prevTiltIm = int(fiducial[2])
                landmarkModelGaps.addLandmark(xCoor=fiducial[0],
                                              yCoor=fiducial[1],
                                              tiltIm=fiducial[2],
                                              chainId=chainId,
                                              xResid=fiducial[3],
                                              yResid=fiducial[4])

            outputSetOfLandmarkModelsGaps.append(landmarkModelGaps)
            outputSetOfLandmarkModelsGaps.update(landmarkModelGaps)
            outputSetOfLandmarkModelsGaps.write()
            self._store()

        """Landmark models with no gaps"""
        if os.path.exists(os.path.join(extraPrefix,
                                       ts.getFirstItem().parseFileName(suffix="_nogaps", extension=".fid"))) and \
                os.path.exists(os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".resid"))):

            paramsNoGapPoint2Model = {
                'inputFile': os.path.join(
                    extraPrefix,
                    ts.getFirstItem().parseFileName(suffix="_nogaps", extension=".fid")
                ),
                'outputFile': os.path.join(
                    extraPrefix,
                    ts.getFirstItem().parseFileName(suffix="_nogaps_fid", extension=".txt")
                )
            }

            argsNoGapPoint2Model = "-InputFile %(inputFile)s " \
                                   "-OutputFile %(outputFile)s"

            Plugin.runImod(self, 'model2point', argsNoGapPoint2Model % paramsNoGapPoint2Model)

            outputSetOfLandmarkModelsNoGaps = self._createSetOfLandmarkModels(suffix='NoGaps')
            outputSetOfLandmarkModelsNoGaps.copyInfo(self.inputTiltSeries.get())
            self._defineOutputs(outputSetOfLandmarkModelsNoGaps=outputSetOfLandmarkModelsNoGaps)
            self._defineSourceRelation(self.inputTiltSeries, outputSetOfLandmarkModelsNoGaps)

            fiducialNoGapFilePath = os.path.join(
                extraPrefix,
                ts.getFirstItem().parseFileName(suffix="_nogaps_fid", extension=".txt")
            )

            fiducialNoGapList = utils.formatFiducialList(fiducialNoGapFilePath)

            fiducialModelNoGapPath = os.path.join(
                extraPrefix,
                ts.getFirstItem().parseFileName(suffix="_nogaps", extension=".fid")
            )

            landmarkModelNoGapsFilePath = os.path.join(
                extraPrefix,
                ts.getFirstItem().parseFileName(suffix="_nogaps", extension=".sfid")
            )

            landmarkModelNoGapsResidPath = os.path.join(
                extraPrefix,
                ts.getFirstItem().parseFileName(extension=".resid")

            )

            fiducialNoGapsResidList = utils.formatFiducialResidList(landmarkModelNoGapsResidPath)

            landmarkModelNoGaps = tomoObj.LandmarkModel(tsId, landmarkModelNoGapsFilePath, fiducialModelNoGapPath)

            prevTiltIm = 0
            chainId = 0
            indexFake = 0
            for fiducial in fiducialNoGapList:
                if int(float(fiducial[2])) <= prevTiltIm:
                    chainId += 1
                prevTiltIm = int(float(fiducial[2]))
                if indexFake < len(fiducialNoGapsResidList) and \
                        fiducial[2] == fiducialNoGapsResidList[indexFake][2]:
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
            outputSetOfLandmarkModelsNoGaps.write()
            self._store()

        """Full reconstructed tomogram"""
        if os.path.exists(os.path.join(extraPrefix, ts.getFirstItem().parseFileName(suffix="_full", extension=".rec"))):
            outputSetOfFullTomograms = self._createSetOfTomograms(suffix='Full')
            outputSetOfFullTomograms.copyInfo(self.inputTiltSeries.get())
            self._defineOutputs(outputSetOfFullTomograms=outputSetOfFullTomograms)
            self._defineSourceRelation(self.inputTiltSeries, outputSetOfFullTomograms)

            newTomogram = tomoObj.Tomogram()
            newTomogram.setLocation(os.path.join(
                extraPrefix,
                ts.getFirstItem().parseFileName(suffix="_full", extension=".rec"))
            )
            outputSetOfFullTomograms.append(newTomogram)
            outputSetOfFullTomograms.update(newTomogram)
            outputSetOfFullTomograms.write()
            self._store()

        """Post-processed reconstructed tomogram"""
        if os.path.exists(os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".rec"))):
            outputSetOfPostProcessTomograms = self._createSetOfTomograms()
            outputSetOfPostProcessTomograms.copyInfo(self.inputTiltSeries.get())
            self._defineOutputs(outputSetOfPostProcessTomograms=outputSetOfPostProcessTomograms)
            self._defineSourceRelation(self.inputTiltSeries, outputSetOfPostProcessTomograms)

            newTomogram = tomoObj.Tomogram()
            newTomogram.setLocation(os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".rec")))
            outputSetOfPostProcessTomograms.append(newTomogram)
            outputSetOfPostProcessTomograms.update(newTomogram)
            outputSetOfPostProcessTomograms.write()
            self._store()
        self.closeMappers()

    # --------------------------- UTILS functions ----------------------------
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
Setup.WholeTomogramSampleA=false
Setup.DatasetName=%(name)s
Setup.FiducialDiameter=%(markerDiameter)f
Setup.WholeTomogramSampleB=false
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
Setup.Binning=%(binning)s
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
Setup.DefaultGpuProcessing=false
Setup.Track.A.SeedModel.Auto=true
Setup.Combine.PatchSize=M
Setup.AxisB.TiltAngle.Type=Extract
Setup.Combine.PatchBoundaryZMin=0
Setup.ProjectLog.FrameLocation.Y=55
Setup.ProjectLog.FrameLocation.X=95
Setup.DefaultParallel=false
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

    def getPixSizeFromDimensions(self, outputDim):
        ih = ImageHandler()
        originalDim, _, _, _ = ih.getDimensions(self.inputTiltSeries.get().getFirstItem().getFileName())
        return self.inputTiltSeries.get().getSamplingRate() * round(originalDim/outputDim)

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = ["The following outputs have been generated from the "
                   "operations performed over the input tilt-series:"]

        if hasattr(self, 'outputPrealignedSetOfTiltSeries'):
            summary.append("- Tilt-series prealignment.")

        if hasattr(self, 'outputAlignedSetOfTiltSeries'):
            summary.append("- Tilt-series alignment.")

        if hasattr(self, 'outputSetOfCoordinates3D'):
            summary.append("- Landmark 3D coordinates have been extracted.")

        if hasattr(self, 'outputSetOfLandmarkModelsGaps'):
            summary.append("- Landmark model with gaps has been generated.")

        if hasattr(self, 'outputSetOfLandmarkModelsNoGaps'):
            summary.append("- Landmark model without gaps has been generated.")

        if hasattr(self, 'outputSetOfFullTomograms'):
            summary.append("- Full raw reconstructed tomogram.")

        if hasattr(self, 'outputSetOfPostProcessTomograms'):
            summary.append("- Post processed reconstructed tomogram.")

        if summary == ["The following operations has been performed over the input tilt-series:"]:
            summary = ["Output classes not ready yet."]

        return summary

    def _methods(self):
        methods = ["The following outputs have been generated from the "
                   "operations performed over the input tilt-series:"]

        if hasattr(self, 'outputPrealignedSetOfTiltSeries'):
            methods.append("- Tilt-series prealignment.")

        if hasattr(self, 'outputAlignedSetOfTiltSeries'):
            methods.append("- Tilt-series alignment.")

        if hasattr(self, 'outputSetOfCoordinates3D'):
            methods.append("- Landmark 3D coordinates have been extracted.")

        if hasattr(self, 'outputSetOfLandmarkModelsGaps'):
            methods.append("- Landmark model with gaps has been generated.")

        if hasattr(self, 'outputSetOfLandmarkModelsNoGaps'):
            methods.append("- Landmark model without gaps has been generated.")

        if hasattr(self, 'outputSetOfFullTomograms'):
            methods.append("- Full raw reconstructed tomogram.")

        if hasattr(self, 'outputSetOfPostProcessTomograms'):
            methods.append("- Post processed reconstructed tomogram.")

        if methods == ["The following operations has been performed over the input tilt-series:"]:
            methods = ["Output classes not ready yet."]

        return methods
