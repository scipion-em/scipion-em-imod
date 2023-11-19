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

import pyworkflow as pw
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.emlib.image import ImageHandler
import tomo.objects as tomoObj

from .. import Plugin, utils
from .protocol_base import (OUTPUT_TILTSERIES_NAME, ProtImodBase,
                            OUTPUT_TS_COORDINATES_NAME)


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
                      label='Input set of tilt-series',
                      help='Input set of tilt-series to be processed with eTomo.')

        form.addParam('markersDiameter',
                      params.FloatParam,
                      default=10,
                      label='Fiducial markers diameter (nm)',
                      help='Diameter of gold beads in nanometers.')

        form.addParam('applyAlignment',
                      params.BooleanParam,
                      default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Apply transformation matrix?',
                      help='Apply the transformation matrix if input'
                           'tilt series have it.')

        form.addParallelSection(threads=1, mpi=0)

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
        setOftiltSeries = self.inputSetOfTiltSeries.get()
        view = ImodGenericView(None, self, setOftiltSeries,
                               isInteractive=True)
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

    def convertInputStep(self, ts, **kwargs):
        tsId = ts.getTsId()
        acq = ts.getAcquisition()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)
        firstItem = ts.getFirstItem()

        outputTsFileName = self.getFilePath(ts, extension=".mrc")

        """Apply the transformation from the input tilt-series"""
        if self.applyAlignment:
            ts.applyTransform(outputTsFileName)
        else:
            path.createAbsLink(os.path.abspath(firstItem.getFileName()),
                               outputTsFileName)

        """Generate angle file"""
        angleFilePath = self.getFilePath(ts, extension=".rawtlt")
        ts.generateTltFile(angleFilePath)

        """Generate etomo config file"""
        args = '-name %s ' % firstItem.parseFileName(extension="")
        args += '-gold %0.3f ' % self.markersDiameter

        # Imod use the pixel size in NM
        pixelSizeNm = ts.getSamplingRate() / 10.

        args += '-pixel %0.3f ' % pixelSizeNm
        args += '-rotation %0.3f ' % acq.getTiltAxisAngle()
        args += '-userawtlt -fei 1 -change "%s/SystemTemplate/cryoSample.adoc" ' % Plugin.getHome()

        # 0 for output image files to have descriptive extensions like ".preali", 1 for extension ".mrc", or 2 for
        # extension ".hdf". In the latter two cases the usual descriptive text is put before the extension, and command
        # files will contain an environment variable setting to make programs generate files of the corresponding type.
        # From: https://bio3d.colorado.edu/imod/doc/man/copytomocoms.html
        args += '-NamingStyle 1 '

        # Extension of raw stack excluding the period.  If this is not specified, the program will assume the extension
        # ".st" unless the -style option is entered.  With a -style option and no specified stack extension, it will
        # look for ".st", ".mrc", ".hdf",".tif", and ".tiff" and require that only one of those types is present. With
        # this entry, which could in principle be arbitrary, it will not care if files with other extensions are
        # present.
        # From: https://bio3d.colorado.edu/imod/doc/man/copytomocoms.html
        args += '-StackExtension mrc '

        args += f'-binning 1.0 -Cs {acq.getSphericalAberration()} -voltage {int(acq.getVoltage())} '

        if ts.getExcludedViewsIndex():
            args += f'-ViewsToSkip {",".join(ts.getExcludedViewsIndex())} '

        Plugin.runImod(self, 'copytomocoms', args, cwd=extraPrefix)

        edfFn = self.getFilePath(ts, extension=".edf")
        minTilt = min(utils.formatAngleList(angleFilePath))
        self._writeEtomoEdf(edfFn,
                            {
                                'date': pw.utils.prettyTime(),
                                'name': firstItem.parseFileName(extension=''),
                                'pixelSize': pixelSizeNm,
                                'version': pw.__version__,
                                'minTilt': minTilt,
                                'markerDiameter': self.markersDiameter,
                                'rotationAngle': ts.getAcquisition().getTiltAxisAngle(),
                                'imodDir': Plugin.getHome(),
                                'useCpu': self.numberOfThreads > 1
                            })

    def runEtomo(self, ts):
        tsId = ts.getTsId()
        edfFilePath = self._getExtraPath(os.path.join(tsId,
                                                      ts.getFirstItem().parseFileName(extension=".edf")))

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
        self.FiducialModelNoGaps = None  # This will reset the output. Is this what we want?
        setOfTSCoords = None
        outputSetOfFullTomograms = None
        outputSetOfPostProcessTomograms = None
        setOfTiltSeries = self.inputSetOfTiltSeries.get()
        ih = ImageHandler()

        for ts in setOfTiltSeries:
            self.inputTiltSeries = ts
            tsId = ts.getTsId()

            """Prealigned tilt-series"""
            prealiFilePath = self.getFilePath(ts, suffix="_preali", extension=".mrc")
            if os.path.exists(prealiFilePath):
                xPrealiDims, newPixSize = self.getNewPixAndDim(ih, prealiFilePath)
                self.debug(f"{prealiFilePath}: pix = {newPixSize}, dims = {xPrealiDims}")
                if outputPrealiSetOfTiltSeries is None:
                    outputPrealiSetOfTiltSeries = self._createSetOfTiltSeries(suffix='_prealigned')
                    outputPrealiSetOfTiltSeries.copyInfo(setOfTiltSeries)
                    outputPrealiSetOfTiltSeries.setSamplingRate(newPixSize)
                    self._defineOutputs(PrealignedTiltSeries=outputPrealiSetOfTiltSeries)
                    self._defineSourceRelation(self.inputTiltSeries,
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
                    sliceIndex = newTi.getIndex()
                    newTi.setLocation(sliceIndex, prealiFilePath)
                    if sliceIndex in excludedViewList:
                        newTi.setEnabled(False)
                    newTs.append(newTi)

                newTs.setDim(xPrealiDims)
                newTs.write(properties=False)

                outputPrealiSetOfTiltSeries.update(newTs)
                outputPrealiSetOfTiltSeries.write()
                self._store(outputPrealiSetOfTiltSeries)

            """Aligned tilt-series"""
            aligFilePath = self.getFilePath(ts, suffix="_ali", extension=".mrc")
            if os.path.exists(aligFilePath):
                aliDims, newPixSize = self.getNewPixAndDim(ih, aligFilePath)
                self.debug(f"{aligFilePath}: pix = {newPixSize}, dims = {aliDims}")

                if outputAliSetOfTiltSeries is None:
                    outputAliSetOfTiltSeries = self._createSetOfTiltSeries(suffix='_aligned')
                    outputAliSetOfTiltSeries.copyInfo(setOfTiltSeries)
                    outputAliSetOfTiltSeries.setSamplingRate(newPixSize)
                    self._defineOutputs(**{OUTPUT_TILTSERIES_NAME: outputAliSetOfTiltSeries})
                    self._defineSourceRelation(self.inputSetOfTiltSeries,
                                               outputAliSetOfTiltSeries)
                else:
                    outputAliSetOfTiltSeries.enableAppend()

                newTs = ts.clone()
                newTs.copyInfo(ts)
                newTs.setInterpolated(True)
                outputAliSetOfTiltSeries.append(newTs)

                tltFilePath = self.getFilePath(ts, suffix='_fid',
                                               extension=".tlt")
                if os.path.exists(tltFilePath):
                    tltList = utils.formatAngleList(tltFilePath)
                    self.debug("%s read: %s" % (tltFilePath, tltList))
                else:
                    tltList = None

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
                    sliceIndex = newTi.getIndex()
                    self.debug("Slice index is %s" % sliceIndex)
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
                newTs.write(properties=False)

                outputAliSetOfTiltSeries.update(newTs)
                outputAliSetOfTiltSeries.write()
                self._store(outputAliSetOfTiltSeries)

            """Output set of coordinates 3D (associated to the aligned tilt-series)"""
            coordFilePath = self.getFilePath(ts, suffix='fid', extension=".xyz")

            if os.path.exists(coordFilePath) and outputAliSetOfTiltSeries is not None:
                if setOfTSCoords is None:
                    setOfTSCoords = tomoObj.SetOfTiltSeriesCoordinates.create(self._getPath(),
                                                                              suffix='Fiducials3D')
                    setOfTSCoords.setSetOfTiltSeries(outputAliSetOfTiltSeries)
                    self._defineOutputs(**{OUTPUT_TS_COORDINATES_NAME: setOfTSCoords})
                    self._defineSourceRelation(self.inputSetOfTiltSeries,
                                               setOfTSCoords)
                else:
                    setOfTSCoords.enableAppend()

                coordList, _, _ = utils.format3DCoordinatesList(coordFilePath)

                for element in coordList:
                    newCoord3D = tomoObj.TiltSeriesCoordinate()
                    newCoord3D.setTsId(ts.getTsId())
                    self.debug("Setting tilt series coordinate x, y, z: %s, %s, %s." % (
                        element[0], element[1], element[2]))
                    newCoord3D.setX(element[0])
                    newCoord3D.setY(element[1])
                    newCoord3D.setZ(element[2])

                    setOfTSCoords.append(newCoord3D)
                    setOfTSCoords.update(newCoord3D)

                setOfTSCoords.write()
                self._store(setOfTSCoords)

            """Landmark models with no gaps"""
            modelFilePath = self.getFilePath(ts, suffix="_nogaps", extension=".fid")
            residFilePath = self.getFilePath(ts, extension=".resid")

            if os.path.exists(modelFilePath) and os.path.exists(residFilePath):
                modelFilePathTxt = self.getFilePath(ts, suffix="_nogaps_fid",
                                                    extension=".txt")

                paramsNoGapPoint2Model = {
                    'inputFile': modelFilePath,
                    'outputFile': modelFilePathTxt
                }

                argsNoGapPoint2Model = "-InputFile %(inputFile)s " \
                                       "-OutputFile %(outputFile)s"

                Plugin.runImod(self, 'model2point', argsNoGapPoint2Model % paramsNoGapPoint2Model)

                outputSetOfLandmarkModelsNoGaps = self.getOutputFiducialModelNoGaps(outputPrealiSetOfTiltSeries)

                fiducialNoGapList = utils.formatFiducialList(modelFilePathTxt)

                landmarkModelNoGapsFilePath = self.getFilePath(ts, suffix="_nogaps",
                                                               extension=".sfid")
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
                outputSetOfLandmarkModelsNoGaps.write()
                self._store(outputSetOfLandmarkModelsNoGaps)

            """Full reconstructed tomogram"""
            reconstructTomoFilePath = self.getFilePath(ts, suffix="_full_rec",
                                                       extension=".mrc")
            if os.path.exists(reconstructTomoFilePath):
                tomoDims, newPixSize = self.getNewPixAndDim(ih, reconstructTomoFilePath)
                self.debug(f"{reconstructTomoFilePath}: pix = {newPixSize}, dims = {tomoDims}")
                if outputSetOfFullTomograms is None:
                    outputSetOfFullTomograms = self._createSetOfTomograms(suffix='_raw')
                    outputSetOfFullTomograms.copyInfo(setOfTiltSeries)
                    outputSetOfFullTomograms.setSamplingRate(newPixSize)
                    self._defineOutputs(FullTomograms=outputSetOfFullTomograms)
                    self._defineSourceRelation(self.inputSetOfTiltSeries,
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
                outputSetOfFullTomograms.write()
                self._store(outputSetOfFullTomograms)

            """Post-processed reconstructed tomogram"""
            posprocessedRecTomoFilePath = self.getFilePath(ts, suffix="_rec",
                                                           extension=".mrc")
            if os.path.exists(posprocessedRecTomoFilePath):
                tomoDims, newPixSize = self.getNewPixAndDim(ih, posprocessedRecTomoFilePath)
                self.debug(f"{posprocessedRecTomoFilePath}: pix = {newPixSize}, dims = {tomoDims}")
                if outputSetOfPostProcessTomograms is None:
                    outputSetOfPostProcessTomograms = self._createSetOfTomograms()
                    outputSetOfPostProcessTomograms.copyInfo(setOfTiltSeries)
                    outputSetOfPostProcessTomograms.setSamplingRate(newPixSize)
                    self._defineOutputs(PostProcessTomograms=outputSetOfPostProcessTomograms)
                    self._defineSourceRelation(self.inputSetOfTiltSeries,
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
                outputSetOfPostProcessTomograms.write()
                self._store(outputSetOfPostProcessTomograms)

        self.closeMappers()

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

    def getNewPixAndDim(self, ih, fn):
        dims = ih.getDimensions(fn)
        dims = dims[:-1]
        origDimX, origDimY, _, _ = ih.getDimensions(self.inputTiltSeries.getFirstItem().getFileName())
        originalDim = max(origDimX, origDimY)
        outputDim = max(dim for dim in dims[:2])
        newPixSize = self.inputTiltSeries.getSamplingRate() * round(originalDim / outputDim)

        return dims, newPixSize

    def getExcludedViewList(self, fn, reservedWord="ExcludeList"):
        with open(fn) as f:
            data = f.readlines()
        excludedViewList = []

        for line in data:
            if line.startswith(reservedWord):
                excludedRange = line.strip().split("\t")[1]
                for part in excludedRange.split(','):
                    if '-' in part:
                        a, b = map(int, part.split('-'))
                        excludedViewList.extend(range(a, b + 1))
                    else:
                        a = int(part)
                        excludedViewList.append(a)
                break
        return excludedViewList

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = ["The following outputs have been generated from the "
                   "operations performed on the input tilt-series:"]

        if self.PrealignedTiltSeries:
            summary.append("- Pre-aligned tilt-series")

        if self.TiltSeries:
            summary.append("- Aligned tilt-series")

        if self.TiltSeriesCoordinates:
            summary.append("- Landmark 3D coordinates")

        if self.FiducialModelGaps:
            summary.append("- Landmark model with gaps")

        if self.FiducialModelNoGaps:
            summary.append("- Landmark model without gaps")

        if self.FullTomograms:
            summary.append("- Raw reconstructed tomogram")

        if self.PostProcessedTomograms:
            summary.append("- Post-processed tomogram")

        if len(summary) == 1:
            summary = ["Outputs are not ready yet"]

        return summary
