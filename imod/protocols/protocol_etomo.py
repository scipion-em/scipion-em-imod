# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [2]
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
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
from tomo.convert import writeTiStack
import tomo.constants as constants
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
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of Tilt-Series',
                      help='Input set of tilt-series to be processed with etomo.')

        form.addParam('excludeList',
                      params.StringParam,
                      default='',
                      label='Exclusion list',
                      help='Provide tilt images IDs (usually starting at 1) that you want to exclude from the '
                           'processing separated by blank spaces.')

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
        self.inputTiltSeries = None
        self._insertFunctionStep(self.runEtomoStep, interactive=True)

    # --------------------------- STEPS functions ----------------------------
    def runEtomoStep(self):
        from imod.viewers import ImodGenericViewer
        setOftiltSeries = self.inputSetOfTiltSeries.get()
        view = ImodGenericViewer(None, self, setOftiltSeries,
                                 displayAllButton=False, isInteractive=True,
                                 itemDoubleClick=True)
        view.show()
        self.createOutput()

    def runAllSteps(self, obj):
        for item in self.inputSetOfTiltSeries.get():
            if item.getTsId() == obj.getTsId():
                self.runEtomo(item)
                break

    def getFilePath(self, ts, suffix="", extension=""):
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        return os.path.join(extraPrefix, ts.getFirstItem().parseFileName(suffix=suffix,
                                                                         extension=extension))

    def convertInputStep(self, ts):
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)
        firstItem = ts.getFirstItem()

        outputTsFileName = self.getFilePath(ts, extension=".st")

        """Apply transformation matrices and remove excluded views"""
        if self.excludeList.get() == '':

            """Apply the transformation form the input tilt-series"""
            ts.applyTransform(outputTsFileName)

            """Generate angle file"""
            angleFilePath = self.getFilePath(ts, extension=".rawtlt")
            ts.generateTltFile(angleFilePath)

        else:
            interpolatedTsFileName = os.path.join(tmpPrefix, firstItem.parseFileName())
            angleFilePath = self.getFilePath(ts, extension=".rawtlt")

            """Apply the transformation form the input tilt-series and generate a new ts object"""
            ts.applyTransform(interpolatedTsFileName)

            interpolatedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Interpolated')
            interpolatedSetOfTiltSeries.copyInfo(ts)
            interpolatedSetOfTiltSeries.setDim(ts.getDim())

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
        args = '-name %s ' % firstItem.parseFileName(extension="")
        args += '-gold %0.3f ' % self.markersDiameter

        # Imod use the pixel size in NM
        pixelSizeNm = ts.getSamplingRate() / 10.

        args += '-pixel %0.3f ' % pixelSizeNm
        args += '-rotation %0.3f ' % self.rotationAngle
        args += '-userawtlt '

        # 0 for output image files to have descriptive extensions like ".preali", 1 for extension ".mrc", or 2 for
        # extension ".hdf". In the latter two cases the usual descriptive text is put before the extension, and command
        # files will contain an environment variable setting to make programs generate files of the corresponding type.
        # From: https://bio3d.colorado.edu/imod/doc/man/copytomocoms.html
        args += '-NamingStyle 0 '


        # Extension of raw stack excluding the period.  If this is not specified, the program will assume the extension
        # ".st" unless the -style option is entered.  With a -style option and no specified stack extension, it will
        # look for ".st", ".mrc", ".hdf",".tif", and ".tiff" and require that only one of those types is present. With
        # this entry, which could in principle be arbitrary, it will not care if files with other extensions are
        # present.
        # From: https://bio3d.colorado.edu/imod/doc/man/copytomocoms.html
        args += '-StackExtension ""'

        Plugin.runImod(self, 'copytomocoms', args, cwd=extraPrefix)

        edfFn = self.getFilePath(ts, extension=".edf")
        minTilt = min(utils.formatAngleList(self.getFilePath(ts, extension=".rawtlt")))
        self._writeEtomoEdf(edfFn,
                            {
                                'date': pw.utils.prettyTime(),
                                'name': firstItem.parseFileName(extension=''),
                                'pixelSize': pixelSizeNm,
                                'version': pw.__version__,
                                'minTilt': minTilt,
                                'markerDiameter': self.markersDiameter,
                                'rotationAngle': self.rotationAngle
                            })

    def runEtomo(self, ts):
        self.convertInputStep(ts)
        if ts is not None:
            tsId = ts.getTsId()
            extraPrefix = self._getExtraPath(tsId)
            args = '--fg '
            args += ts.getFirstItem().parseFileName(extension=".edf")
            Plugin.runImod(self, 'etomo', args, cwd=extraPrefix)

    def createOutput(self):
        outputPrealiSetOfTiltSeries = None
        outputAliSetOfTiltSeries = None
        outputSetOfLandmarkModelsNoGaps = None
        outputSetOfCoordinates3D = None
        outputSetOfFullTomograms = None
        outputSetOfPostProcessTomograms = None

        for ts in self.inputSetOfTiltSeries.get():
            self.inputTiltSeries = ts
            tsId = ts.getTsId()
            """Prealigned tilt-series"""
            prealiFilePath = self.getFilePath(ts, extension=".preali")
            if os.path.exists(prealiFilePath):
                if outputPrealiSetOfTiltSeries is None:
                    outputPrealiSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Preali')
                    outputPrealiSetOfTiltSeries.copyInfo(self.inputTiltSeries)
                    outputPrealiSetOfTiltSeries.setDim(self.inputTiltSeries.getDim())
                    self._defineOutputs(outputPrealignedSetOfTiltSeries=outputPrealiSetOfTiltSeries)
                    self._defineSourceRelation(self.inputTiltSeries,
                                               outputPrealiSetOfTiltSeries)

                newTs = ts.clone()
                newTs.copyInfo(ts)
                outputPrealiSetOfTiltSeries.append(newTs)

                ih = ImageHandler()
                index = 0
                for tiltImage in ts.iterItems(iterate=False):
                    newTi = tiltImage.clone()
                    newTi.copyInfo(tiltImage, copyId=True)
                    newTi.setLocation(index + 1, prealiFilePath)
                    index += 1
                    xPreali, _, _, _ = ih.getDimensions(newTi.getFileName()+":mrc")
                    newTi.setSamplingRate(self.getPixSizeFromDimensions(xPreali))
                    newTs.append(newTi)

                xPreali, yPreali, zPreali, _ = ih.getDimensions(newTs.getFirstItem().getFileName()+":mrc")
                newTs.setDim((xPreali, yPreali, zPreali))

                newTs.write(properties=False)

                outputPrealiSetOfTiltSeries.setSamplingRate(self.getPixSizeFromDimensions(xPreali))
                outputPrealiSetOfTiltSeries.write()
                self._store(outputPrealiSetOfTiltSeries)

            """Aligned tilt-series"""
            aligFilePath = self.getFilePath(ts, extension=".ali")
            if os.path.exists(aligFilePath):
                if outputAliSetOfTiltSeries is None:
                    outputAliSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Ali')
                    outputAliSetOfTiltSeries.copyInfo(self.inputTiltSeries)
                    outputAliSetOfTiltSeries.setDim(self.inputTiltSeries.getDim())
                    self._defineOutputs(outputAlignedSetOfTiltSeries=outputAliSetOfTiltSeries)
                    self._defineSourceRelation(self.inputSetOfTiltSeries,
                                               outputAliSetOfTiltSeries)

                newTs = ts.clone()
                newTs.copyInfo(ts)
                outputAliSetOfTiltSeries.append(newTs)

                ih = ImageHandler()

                tltFilePath = self.getFilePath(ts, suffix='_fid',
                                               extension=".tlt")
                if os.path.exists(tltFilePath):
                    tltList = utils.formatAngleList(tltFilePath)
                else:
                    tltList = None
                index = 0
                for tiltImage in ts.iterItems(iterate=False):
                    newTi = tiltImage.clone()
                    newTi.copyInfo(tiltImage, copyId=True)
                    newTi.setLocation(index + 1,   self.getFilePath(ts, extension=".ali"))
                    if tltList is not None:
                        newTi.setTiltAngle(float(tltList[index]))
                    xAli, _, _, _ = ih.getDimensions(newTi.getFileName()+":mrc")
                    newTi.setSamplingRate(self.getPixSizeFromDimensions(xAli))
                    newTs.append(newTi)

                xAli, yAli, zAli, _ = ih.getDimensions(newTs.getFirstItem().getFileName() + ":mrc")
                newTs.setDim((xAli, yAli, zAli))

                newTs.write(properties=False)

                outputAliSetOfTiltSeries.setSamplingRate(self.getPixSizeFromDimensions(xAli))
                outputAliSetOfTiltSeries.write()
                self._store(outputAliSetOfTiltSeries)

            """Output set of coordinates 3D (associated to the aligned tilt-series)"""
            coordFilePath = self.getFilePath(ts, suffix='fid',
                                             extension=".xyz")
            if os.path.exists(coordFilePath):
                if outputSetOfCoordinates3D is None:
                    outputSetOfCoordinates3D = self._createSetOfCoordinates3D(volSet=outputAliSetOfTiltSeries,
                                                                              suffix='LandmarkModel')
                    outputSetOfCoordinates3D.setSamplingRate(self.inputTiltSeries.getSamplingRate())
                    outputSetOfCoordinates3D.setPrecedents(outputAliSetOfTiltSeries)
                    self._defineOutputs(outputSetOfCoordinates3D=outputSetOfCoordinates3D)
                    self._defineSourceRelation(self.inputSetOfTiltSeries,
                                               outputSetOfCoordinates3D)

                coordList = utils.format3DCoordinatesList(coordFilePath, xAli, yAli)
                for element in coordList:
                    newCoord3D = tomoObj.Coordinate3D(x=element[0],
                                                      y=element[1],
                                                      z=element[2])
                    newCoord3D.setVolume(ts)
                    newCoord3D.setVolId(ts.getObjId())
                    newCoord3D.write(properties=False)
                    outputSetOfCoordinates3D.append(newCoord3D)
                outputSetOfCoordinates3D.write()
                self._store(outputSetOfCoordinates3D)

            """Landmark models with no gaps"""
            if (os.path.exists(self.getFilePath(ts, suffix="_nogaps",
                                                extension=".fid")) and
                    os.path.exists(self.getFilePath(ts, extension=".resid"))):

                paramsNoGapPoint2Model = {
                    'inputFile': self.getFilePath(ts, suffix="_nogaps",
                                                  extension=".fid"),
                    'outputFile': self.getFilePath(ts, suffix="_nogaps_fid",
                                                   extension=".txt")
                }

                argsNoGapPoint2Model = "-InputFile %(inputFile)s " \
                                       "-OutputFile %(outputFile)s"

                Plugin.runImod(self, 'model2point', argsNoGapPoint2Model % paramsNoGapPoint2Model)

                if outputSetOfLandmarkModelsNoGaps is None:
                    outputSetOfLandmarkModelsNoGaps = self._createSetOfLandmarkModels(suffix='NoGaps')
                    outputSetOfLandmarkModelsNoGaps.copyInfo(self.inputTiltSeries)
                    self._defineOutputs(outputSetOfLandmarkModelsNoGaps=outputSetOfLandmarkModelsNoGaps)
                    self._defineSourceRelation(self.inputSetOfTiltSeries,
                                               outputSetOfLandmarkModelsNoGaps)

                fiducialNoGapFilePath = self.getFilePath(ts, suffix="_nogaps_fid",
                                                         extension=".txt")

                fiducialNoGapList = utils.formatFiducialList(fiducialNoGapFilePath)

                fiducialModelNoGapPath = self.getFilePath(ts, suffix="_nogaps", extension=".fid")
                landmarkModelNoGapsFilePath = self.getFilePath(ts, suffix="_nogaps", extension=".sfid")
                landmarkModelNoGapsResidPath = self.getFilePath(ts, extension=".resid")
                fiducialNoGapsResidList = utils.formatFiducialResidList(landmarkModelNoGapsResidPath)
                landmarkModelNoGaps = tomoObj.LandmarkModel(tsId,
                                                            landmarkModelNoGapsFilePath,
                                                            fiducialModelNoGapPath)

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
                self._store(outputSetOfLandmarkModelsNoGaps)

            """Full reconstructed tomogram"""
            reconstructTomoFilePath = self.getFilePath(ts, suffix="_full",
                                                       extension=".rec")
            if os.path.exists(reconstructTomoFilePath):
                if outputSetOfFullTomograms is None:
                    outputSetOfFullTomograms = self._createSetOfTomograms(suffix='Full')
                    outputSetOfFullTomograms.copyInfo(self.inputTiltSeries)
                    self._defineOutputs(outputSetOfFullTomograms=outputSetOfFullTomograms)
                    self._defineSourceRelation(self.inputSetOfTiltSeries,
                                               outputSetOfFullTomograms)

                newTomogram = tomoObj.Tomogram()
                newTomogram.setLocation(reconstructTomoFilePath)
                newTomogram.setSamplingRate(ts.getSamplingRate())
                outputSetOfFullTomograms.append(newTomogram)
                outputSetOfFullTomograms.write()
                self._store(outputSetOfFullTomograms)

            """Post-processed reconstructed tomogram"""
            posprocessedRecTomoFilePath = self.getFilePath(ts, extension=".rec")
            if os.path.exists(posprocessedRecTomoFilePath):
                if outputSetOfPostProcessTomograms is None:
                    outputSetOfPostProcessTomograms = self._createSetOfTomograms()
                    outputSetOfPostProcessTomograms.copyInfo(self.inputTiltSeries)
                    self._defineOutputs(outputSetOfPostProcessTomograms=outputSetOfPostProcessTomograms)
                    self._defineSourceRelation(self.inputSetOfTiltSeries,
                                               outputSetOfPostProcessTomograms)

                newTomogram = tomoObj.Tomogram()
                newTomogram.setLocation(posprocessedRecTomoFilePath)
                outputSetOfPostProcessTomograms.append(newTomogram)
                outputSetOfPostProcessTomograms.write()
                self._store(outputSetOfPostProcessTomograms)
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
        originalDim, _, _, _ = ih.getDimensions(self.inputTiltSeries.getFirstItem().getFileName())
        return self.inputTiltSeries.getSamplingRate() * round(originalDim/outputDim)


    def getResizeFactorFromDimensions(self, outputDim):
        ih = ImageHandler()
        originalDim, _, _, _ = ih.getDimensions(self.inputTiltSeries.get().getFirstItem().getFileName())
        return  round(outputDim / originalDim)

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
