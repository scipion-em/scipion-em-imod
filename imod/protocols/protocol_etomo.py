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
from imod.protocols.protocol_base import OUTPUT_TILTSERIES_NAME, OUTPUT_COORDINATES_3D_NAME, ProtImodBase, \
    OUTPUT_FIDUCIAL_NO_GAPS_NAME
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


class ProtImodEtomo(ProtImodBase):
    """
    Simple wrapper around etomo to manually reconstruct a Tomogram.

    More info:
        https://bio3d.colorado.edu/imod/doc/etomoTutorial.html
    """

    _label = 'Etomo interactive'
    _devStatus = BETA

    def __init__(self, **kwargs):

        super().__init__(**kwargs)
        self.PrealignedTiltSeries = None
        self.FullTomograms = None
        self.PostProcessedTomograms = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of Tilt-Series',
                      help='Input set of tilt-series to be processed with etomo.')

        form.addParam('markersDiameter',
                      params.FloatParam,
                      default=10,
                      label='Fiducial markers diameter (nm)',
                      help='Diameter of gold beads in nanometers.')

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

        """Apply the transformation form the input tilt-series"""
        ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = self.getFilePath(ts, extension=".rawtlt")
        ts.generateTltFile(angleFilePath)

        """Generate etomo config file"""
        args = '-name %s ' % firstItem.parseFileName(extension="")
        args += '-gold %0.3f ' % self.markersDiameter

        # Imod use the pixel size in NM
        pixelSizeNm = ts.getSamplingRate() / 10.

        args += '-pixel %0.3f ' % pixelSizeNm
        args += '-rotation %0.3f ' % ts.getAcquisition().getTiltAxisAngle()
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
                                'rotationAngle': ts.getAcquisition().getTiltAxisAngle()
                            })

    def runEtomo(self, ts):
        tsId = ts.getTsId()
        edfFilePath = self._getExtraPath(os.path.join(tsId, ts.getFirstItem().parseFileName(extension=".edf")))

        if not os.path.exists(edfFilePath):
            self.convertInputStep(ts)

        if ts is not None:
            extraPrefix = self._getExtraPath(tsId)
            args = '--fg '
            args += ts.getFirstItem().parseFileName(extension=".edf")
            Plugin.runImod(self, 'etomo', args, cwd=extraPrefix)

    def createOutput(self):

        outputPrealiSetOfTiltSeries = None
        outputAliSetOfTiltSeries = None
        self.FiducialModelNoGaps = None # This will reset the output. Is this what we want?
        setOfTSCoords = None
        outputSetOfFullTomograms = None
        outputSetOfPostProcessTomograms = None
        setOfTiltSeries = self.inputSetOfTiltSeries.get()

        for ts in setOfTiltSeries:
            self.inputTiltSeries = ts
            tsId = ts.getTsId()

            """Prealigned tilt-series"""
            prealiFilePath = self.getFilePath(ts, extension=".preali")
            if os.path.exists(prealiFilePath):
                if outputPrealiSetOfTiltSeries is None:
                    outputPrealiSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Preali')
                    outputPrealiSetOfTiltSeries.copyInfo(setOfTiltSeries)
                    outputPrealiSetOfTiltSeries.setDim(setOfTiltSeries.getDim())
                    self._defineOutputs(PrealignedTiltSeries=outputPrealiSetOfTiltSeries)
                    self._defineSourceRelation(self.inputTiltSeries,
                                               outputPrealiSetOfTiltSeries)
                else:
                    outputPrealiSetOfTiltSeries.enableAppend()

                newTs = ts.clone()
                newTs.copyInfo(ts)
                outputPrealiSetOfTiltSeries.append(newTs)

                ih = ImageHandler()

                # Getting the excluded views in order to disable the
                # prealigned tilt-series
                xcorrFn = self._getExtraPath(tsId, 'xcorr.com')
                excludedViewList = self.getExcludedViewList(xcorrFn,
                                                            reservedWord='SkipViews')

                for tiltImage in ts.iterItems(iterate=False):
                    newTi = tiltImage.clone()
                    newTi.copyInfo(tiltImage, copyId=True, copyTM=False)
                    newTi.setAcquisition(tiltImage.getAcquisition())
                    sliceIndex = newTi.getIndex()
                    newTi.setLocation(sliceIndex, prealiFilePath)
                    xPreali, _, _, _ = ih.getDimensions(newTi.getFileName() + ":mrc")
                    newTi.setSamplingRate(self.getPixSizeFromDimensions(xPreali))
                    if sliceIndex in excludedViewList:
                        newTi.setEnabled(False)
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
                    outputAliSetOfTiltSeries.copyInfo(setOfTiltSeries)
                    outputAliSetOfTiltSeries.setDim(setOfTiltSeries.getDim())
                    self._defineOutputs(**{OUTPUT_TILTSERIES_NAME:outputAliSetOfTiltSeries})
                    self._defineSourceRelation(self.inputSetOfTiltSeries,
                                               outputAliSetOfTiltSeries)
                else:
                    outputAliSetOfTiltSeries.enableAppend()

                newTs = ts.clone()
                newTs.copyInfo(ts)
                outputAliSetOfTiltSeries.append(newTs)

                ih = ImageHandler()

                tltFilePath = self.getFilePath(ts, suffix='_fid',
                                               extension=".tlt")
                if os.path.exists(tltFilePath):
                    tltList = utils.formatAngleList(tltFilePath)
                    self.debug("%s read: %s" % (tltFilePath, tltList))
                else:
                    tltList = None

                # Getting the excluded views in order to disable in the
                # aligned tiltserie
                alignFn = self._getExtraPath(tsId, 'align.com')
                excludedViewList = self.getExcludedViewList(alignFn,
                                                            reservedWord='ExcludeList')

                for tiltImage in ts.iterItems(iterate=False):
                    newTi = tiltImage.clone()
                    newTi.copyInfo(tiltImage, copyId=True, copyTM=False)
                    newTi.setAcquisition(tiltImage.getAcquisition())
                    sliceIndex = newTi.getIndex()
                    self.debug("Slice index is %s" % sliceIndex)
                    newTi.setLocation(sliceIndex, aligFilePath)
                    if tltList is not None:
                        newTi.setTiltAngle(float(tltList[sliceIndex-1]))
                    xAli, _, _, _ = ih.getDimensions(newTi.getFileName() + ":mrc")
                    newTi.setSamplingRate(self.getPixSizeFromDimensions(xAli))
                    if sliceIndex in excludedViewList:
                        newTi.setEnabled(False)
                    newTs.append(newTi)

                xAli, yAli, zAli, _ = ih.getDimensions(newTs.getFirstItem().getFileName() + ":mrc")
                newTs.setDim((xAli, yAli, zAli))

                newTs.write(properties=False)

                outputAliSetOfTiltSeries.setSamplingRate(self.getPixSizeFromDimensions(xAli))
                outputAliSetOfTiltSeries.write()
                self._store(outputAliSetOfTiltSeries)

                """Output set of coordinates 3D (associated to the aligned tilt-series)"""
                coordFilePath = self.getFilePath(ts, suffix='fid', extension=".xyz")

                if os.path.exists(coordFilePath):
                    if setOfTSCoords is None:
                        setOfTSCoords = tomoObj.SetOfTiltSeriesCoordinates.create(self._getPath(), suffix='LandmarkModel')
                        setOfTSCoords.setSetOfTiltSeries(outputAliSetOfTiltSeries)
                        self._defineOutputs(**{OUTPUT_TILTSERIES_NAME:setOfTSCoords})
                        self._defineSourceRelation(self.inputSetOfTiltSeries,
                                                   setOfTSCoords)
                    else:
                        setOfTSCoords.enableAppend()

                    coordList, xDim, yDim = utils.format3DCoordinatesList(coordFilePath)

                    for element in coordList:
                        newCoord3D = tomoObj.TiltSeriesCoordinate()
                        newCoord3D.setTsId(ts.getTsId())
                        self.debug("Setting tilt series coordinate x, y, z: %s, %s, %s." % (element[0], element[1], element[2]))
                        newCoord3D.setX(element[0])
                        newCoord3D.setY(element[1])
                        newCoord3D.setZ(element[2])

                        setOfTSCoords.append(newCoord3D)
                        # setOfTSCoords.update(newCoord3D)

                    setOfTSCoords.write()
                    self._store(setOfTSCoords)

            """Landmark models with no gaps"""
            modelFilePath = self.getFilePath(ts, suffix="_nogaps", extension=".fid")
            residFilePath = self.getFilePath(ts, extension=".resid")

            if os.path.exists(modelFilePath) and os.path.exists(residFilePath):
                modelFilePathTxt = self.getFilePath(ts, suffix="_nogaps_fid", extension=".txt")

                paramsNoGapPoint2Model = {
                    'inputFile': modelFilePath,
                    'outputFile': modelFilePathTxt
                }

                argsNoGapPoint2Model = "-InputFile %(inputFile)s " \
                                       "-OutputFile %(outputFile)s"

                Plugin.runImod(self, 'model2point', argsNoGapPoint2Model % paramsNoGapPoint2Model)

                outputSetOfLandmarkModelsNoGaps = self.getOutputFiducialModelNoGaps()

                fiducialNoGapList = utils.formatFiducialList(modelFilePathTxt)

                landmarkModelNoGapsFilePath = self.getFilePath(ts, suffix="_nogaps", extension=".sfid")
                fiducialNoGapsResidList = utils.formatFiducialResidList(residFilePath)
                landmarkModelNoGaps = tomoObj.LandmarkModel(tsId,
                                                            landmarkModelNoGapsFilePath,
                                                            modelFilePath)

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
                    outputSetOfFullTomograms.copyInfo(setOfTiltSeries)
                    self._defineOutputs(FullTomograms=outputSetOfFullTomograms)
                    self._defineSourceRelation(self.inputSetOfTiltSeries,
                                               outputSetOfFullTomograms)
                else:
                    outputSetOfFullTomograms.enableAppend()

                newTomogram = tomoObj.Tomogram()

                newTomogram.setLocation(reconstructTomoFilePath)
                newTomogram.setTsId(tsId)
                newTomogram.setSamplingRate(ts.getSamplingRate())

                # Set default tomogram origin
                newTomogram.setOrigin(newOrigin=False)

                outputSetOfFullTomograms.append(newTomogram)
                outputSetOfFullTomograms.write()
                self._store(outputSetOfFullTomograms)

            """Post-processed reconstructed tomogram"""
            posprocessedRecTomoFilePath = self.getFilePath(ts, extension=".rec")
            if os.path.exists(posprocessedRecTomoFilePath):
                if outputSetOfPostProcessTomograms is None:
                    outputSetOfPostProcessTomograms = self._createSetOfTomograms()
                    outputSetOfPostProcessTomograms.copyInfo(setOfTiltSeries)
                    self._defineOutputs(PostProcessTomograms=outputSetOfPostProcessTomograms)
                    self._defineSourceRelation(self.inputSetOfTiltSeries,
                                               outputSetOfPostProcessTomograms)
                else:
                    outputSetOfPostProcessTomograms.enableAppend()

                newTomogram = tomoObj.Tomogram()

                ih = ImageHandler()
                outputDim, _, _, _ = ih.getDimensions(posprocessedRecTomoFilePath)

                newTomogram.setLocation(posprocessedRecTomoFilePath)
                newTomogram.setTsId(tsId)
                newTomogram.setSamplingRate(self.getPixSizeFromDimensions(outputDim))

                # Set default tomogram origin
                newTomogram.setOrigin(newOrigin=False)

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
        return round(outputDim / originalDim)

    def getExcludedViewList(self, fn, reservedWord="ExcludeList"):
        with open(fn, 'r') as f:
            data = f.readlines()
        excludedViewList = []
        for line in data:
            # Skip comments
            if line.startswith("#"):
                continue

            if line.startswith(reservedWord):
                excludedViewList = line.strip().replace('\t', ' ').replace(',', ' ')
                excludedViewList = excludedViewList.split(' ')[1:]
                excludedViewList = [int(sliceIndex) for sliceIndex in excludedViewList]
                break
        return excludedViewList


    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = ["The following outputs have been generated from the "
                   "operations performed over the input tilt-series:"]

        if self.PrealignedTiltSeries:
            summary.append("- Tilt-series prealignment.")

        if self.TiltSeries:
            summary.append("- Tilt-series alignment.")

        if self.Coordinates3D:
            summary.append("- Landmark 3D coordinates have been extracted.")

        if self.FiducialModelGaps:
            summary.append("- Landmark model with gaps has been generated.")

        if self.FiducialModelNoGaps:
            summary.append("- Landmark model without gaps has been generated.")

        if self.FullTomograms:
            summary.append("- Full raw reconstructed tomogram.")

        if self.PostProcessedTomograms:
            summary.append("- Post processed reconstructed tomogram.")

        if summary == ["The following operations has been performed over the input tilt-series:"]:
            summary = ["Output classes not ready yet."]

        return summary

    def _methods(self):
        methods = ["The following outputs have been generated from the "
                   "operations performed over the input tilt-series:"]

        if self.PrealignedTiltSeries:
            methods.append("- Tilt-series prealignment.")

        if self.TiltSeries:
            methods.append("- Tilt-series alignment.")

        if self.Coordinates3D:
            methods.append("- Landmark 3D coordinates have been extracted.")

        if self.FiducialModelGaps:
            methods.append("- Landmark model with gaps has been generated.")

        if self.FiducialModelNoGaps:
            methods.append("- Landmark model without gaps has been generated.")

        if self.FullTomograms:
            methods.append("- Full raw reconstructed tomogram.")

        if self.PostProcessTomograms:
            methods.append("- Post processed reconstructed tomogram.")

        if methods == ["The following operations has been performed over the input tilt-series:"]:
            methods = ["Output classes not ready yet."]

        return methods
