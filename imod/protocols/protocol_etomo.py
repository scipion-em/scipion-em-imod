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
from tomo.objects import Tomogram
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
                           'processing.')

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

        """Apply transformation matrices and remove excluded views"""
        if self.excludeList.get() == '':
            outputTsFileName = os.path.join(extraPrefix, "%s.st" % tsId)
            angleFilePath = os.path.join(extraPrefix, "%s.rawtlt" % tsId)

            """Apply the transformation form the input tilt-series"""
            ts.applyTransform(outputTsFileName)

            """Generate angle file"""
            ts.generateTltFile(angleFilePath)

        else:
            interpolatedTsFileName = os.path.join(tmpPrefix, "%s.st" % tsId)
            outputTsFileName = os.path.join(extraPrefix, "%s.st" % tsId)
            angleFilePath = os.path.join(extraPrefix, "%s.rawtlt" % tsId)

            """Apply the transformation form the input tilt-series and generate a new ts object"""
            ts.applyTransform(interpolatedTsFileName)

            interpolatedTs = tomoObj.TiltSeries(tsId=tsId)
            interpolatedTs.copyInfo(ts)
            outputSetOfTiltSeries.append(interpolatedTs)

            for index, tiltImage in enumerate(ts):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)
                newTi.setLocation(index + 1, interpolatedTsFileName)
                interpolatedTs.append(newTi)
            interpolatedTs.write()

            """Write a new stack discarding excluded tilts"""
            excludeList = map(int, self.excludeList.get().split())
            tiList = [ti.clone() for ti in interpolatedTs]
            tiList.sort(key=lambda ti: ti.getTiltAngle())

            writeTiStack(tiList,
                         outputStackFn=outputTsFileName,
                         outputTltFn=angleFilePath,
                         excludeList=excludeList)

        """Generate etomo config file"""
        workingFolder = self._getExtraPath(tsId)

        args = '-name %s ' % tsId
        args += '-gold %0.3f ' % self.markersDiameter

        # Imod use the pixel size in NM
        pixelSizeNm = ts.getSamplingRate() / 10.

        args += '-pixel %0.3f ' % pixelSizeNm
        args += '-binning %d ' % self.binning
        args += '-rotation %0.3f ' % self.rotationAngle
        args += '-userawtlt'

        Plugin.runImod(self, 'copytomocoms', args, cwd=workingFolder)

        edfFn = os.path.join(workingFolder, '%s.edf' % tsId)
        minTilt = min(utils.formatAngleList(os.path.join(extraPrefix, "%s.rawtlt" % tsId)))
        self._writeEtomoEdf(edfFn,
                            {
                                'date': pw.utils.prettyTime(),
                                'name': tsId,
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
        Plugin.runImod(self, 'etomo', '--fg %s.edf' % tsId, cwd=extraPrefix)
        self.createOutputStep()

    def createOutputStep(self):
        ts = self.inputTiltSeries.get()
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        """Prealigned tilt-series"""
        if os.path.exists(os.path.join(extraPrefix, "%s.preali" % tsId)):
            outputPrealiSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Preali')
            outputPrealiSetOfTiltSeries.copyInfo(self.inputTiltSeries.get())
            outputPrealiSetOfTiltSeries.setDim(self.inputTiltSeries.get().getDim())
            self._defineOutputs(outputPrealignedSetOfTiltSeries=outputPrealiSetOfTiltSeries)
            self._defineSourceRelation(self.inputTiltSeries, outputPrealiSetOfTiltSeries)

            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            outputPrealiSetOfTiltSeries.append(newTs)

            for index, tiltImage in enumerate(ts):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)
                newTi.setLocation(index + 1, (os.path.join(extraPrefix, '%s.preali' % tsId)))
                newTs.append(newTi)

            ih = ImageHandler()
            xPreali, yPreali, _, _ = ih.getDimensions(newTs.getFirstItem().getFileName()+":mrc")
            newTs.setDim((xPreali, yPreali, _))
            newTs.write()

            outputPrealiSetOfTiltSeries.setSamplingRate(self.getPixSizeFromDimensions(xPreali))
            outputPrealiSetOfTiltSeries.update(newTs)
            outputPrealiSetOfTiltSeries.updateDim()
            outputPrealiSetOfTiltSeries.write()
            self._store()

            """Aligned tilt-series"""
            if os.path.exists(os.path.join(extraPrefix, "%s.ali" % tsId)):
                outputAliSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Ali')
                outputAliSetOfTiltSeries.copyInfo(self.inputTiltSeries.get())
                outputAliSetOfTiltSeries.setDim(self.inputTiltSeries.get().getDim())
                self._defineOutputs(outputAlignedSetOfTiltSeries=outputAliSetOfTiltSeries)
                self._defineSourceRelation(self.inputTiltSeries, outputAliSetOfTiltSeries)

                newTs = tomoObj.TiltSeries(tsId=tsId)
                newTs.copyInfo(ts)
                outputAliSetOfTiltSeries.append(newTs)

                tltFilePath = os.path.join(extraPrefix, "%s.ali" % tsId)
                tltList = utils.formatAngleList(tltFilePath)

                for index, tiltImage in enumerate(ts):
                    newTi = tomoObj.TiltImage()
                    newTi.copyInfo(tiltImage, copyId=True)
                    newTi.setLocation(index + 1, (os.path.join(extraPrefix, '%s.ali' % tsId)))
                    newTi.setTiltAngle(float(tltList[index]))
                    newTs.append(newTi)

                ih = ImageHandler()
                xAli, yAli, _, _ = ih.getDimensions(newTs.getFirstItem().getFileName() + ":mrc")
                newTs.setDim((xAli, yAli, _))
                newTs.write()

                outputAliSetOfTiltSeries.setSamplingRate(self.getPixSizeFromDimensions(xAli))
                outputAliSetOfTiltSeries.update(newTs)
                outputAliSetOfTiltSeries.updateDim()
                outputAliSetOfTiltSeries.write()
                self._store()

                """Output set of coordinates 3D (associated to the aligned tilt-series)"""
                if os.path.exists(os.path.join(extraPrefix, "%sfid.xyz" % tsId)):
                    outputSetOfCoordinates3D = self._createSetOfCoordinates3D(volSet=outputAliSetOfTiltSeries,
                                                                              suffix='LandmarkModel')
                    outputSetOfCoordinates3D.copyInfo(self.inputSetOfTiltSeries.get())
                    self._defineOutputs(outputSetOfCoordinates3D=outputSetOfCoordinates3D)
                    self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfCoordinates3D)

                    coordFilePath = os.path.join(extraPrefix, "%sfid.xyz" % tsId)
                    coordList = utils.format3DCoordinatesList(coordFilePath, xAli, yAli)
                    for element in coordList:
                        newCoord3D = tomoObj.Coordinate3D(x=element[0],
                                                          y=element[1],
                                                          z=element[2])
                        newCoord3D.setVolId(tsObjId + 1)
                        newCoord3D.setVolName(tsId)
                        outputSetOfCoordinates3D.append(newCoord3D)
                        outputSetOfCoordinates3D.update(newCoord3D)
                    outputSetOfCoordinates3D.write()
                    self._store()

            """Create the output set of landmark models with gaps"""
            if os.path.exists(os.path.join(extraPrefix, "%s.fid" % tsId)):
                paramsGapPoint2Model = {
                    'inputFile': os.path.join(extraPrefix, '%s.fid' % tsId),
                    'outputFile': os.path.join(extraPrefix, '%s_fid.txt' % tsId)
                }
                argsGapPoint2Model = "-InputFile %(inputFile)s " \
                                     "-OutputFile %(outputFile)s"
                Plugin.runImod(self, 'model2point', argsGapPoint2Model % paramsGapPoint2Model)

                outputFiducialModelGaps = self._createSetOfLandmarkModels(suffix='Gaps')
                outputFiducialModelGaps.copyInfo(self.inputSetOfTiltSeries.get())
                self._defineOutputs(outputFiducialModelGaps=outputFiducialModelGaps)
                self._defineSourceRelation(self.inputSetOfTiltSeries, outputFiducialModelGaps)

                fiducialModelGapPath = os.path.join(extraPrefix, "%s.fid" % tsId)

                landmarkModelGapsFilePath = os.path.join(extraPrefix, "%s_gaps.sfid" % tsId)

                landmarkModelGapsResidPath = os.path.join(extraPrefix, '%s.resid' % tsId)
                fiducialGapResidList = utils.formatFiducialResidList(landmarkModelGapsResidPath)

                landmarkModelGaps = LandmarkModel(tsId, landmarkModelGapsFilePath, fiducialModelGapPath)

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

            """Create the output set of landmark models with no gaps"""
            if os.path.exists(os.path.join(extraPrefix, "%s_nogaps.fid" % tsId)):
                paramsNoGapPoint2Model = {
                    'inputFile': os.path.join(extraPrefix, '%s_nogaps.fid' % tsId),
                    'outputFile': os.path.join(extraPrefix, '%s_nogaps_fid.txt' % tsId)
                }
                argsNoGapPoint2Model = "-InputFile %(inputFile)s " \
                                       "-OutputFile %(outputFile)s"
                Plugin.runImod(self, 'model2point', argsNoGapPoint2Model % paramsNoGapPoint2Model)

                outputFiducialModelNoGaps = self._createSetOfLandmarkModels(suffix='NoGaps')
                outputFiducialModelNoGaps.copyInfo(self.inputSetOfTiltSeries.get())
                self._defineOutputs(outputFiducialModelNoGaps=outputFiducialModelNoGaps)
                self._defineSourceRelation(self.inputSetOfTiltSeries, outputFiducialModelNoGaps)

                fiducialNoGapFilePath = os.path.join(extraPrefix, tsId + "_nogaps_fid.txt")

                fiducialNoGapList = utils.formatFiducialList(fiducialNoGapFilePath)

                fiducialModelNoGapPath = os.path.join(extraPrefix, tsId + "_nogaps.fid")

                landmarkModelNoGapsFilePath = os.path.join(extraPrefix, tsId + "_nogaps.sfid")

                landmarkModelNoGapsResidPath = os.path.join(extraPrefix, '%s.resid' % tsId)
                fiducialNoGapsResidList = utils.formatFiducialResidList(landmarkModelNoGapsResidPath)

                landmarkModelNoGaps = LandmarkModel(tsId, landmarkModelNoGapsFilePath, fiducialModelNoGapPath)

                prevTiltIm = 0
                chainId = 0
                indexFake = 0
                for fiducial in fiducialNoGapList:
                    if int(float(fiducial[2])) <= prevTiltIm:
                        chainId += 1
                    prevTiltIm = int(float(fiducial[2]))
                    if indexFake < len(fiducialNoGapsResidList) and fiducial[2] == fiducialNoGapsResidList[indexFake][2]:
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

    # --------------------------- UTILS functions ----------------------------
    def _writeEtomoEdf(self, fn, paramsDict):
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

    def _getWorkingPath(self, *paths):
        ts = self.inputTiltSeries.get()
        return self._getExtraPath(ts.getTsId(), *paths)

    def _registerTs(self, outputName, outputPath):
        """ Register a tilt-series. """

    def _registerReconsTomo(self):
        outputName = 'outputTomogram'

        if hasattr(self, outputName):
            raise Exception("The output tomogram has already been registered. ")

        ts = self.inputTiltSeries.get()
        tsId = ts.getTsId()
        tomoFn = self._getWorkingPath('%s_full.rec' % tsId)

        if not os.path.exists(tomoFn):
            raise Exception('Output file %s does not exists. ' % tomoFn)

        samplingRate = ts.getSamplingRate()
        if self.binning > 1:
            samplingRate *= self.binning.get()

        t = Tomogram(location=tomoFn + ':mrc')
        t.setSamplingRate(samplingRate)
        t.setObjId(ts.getObjId())
        t.setTsId(tsId)

        return outputName, t

    def registerOutput(self, outputKey):
        """ Method used to register output selected by the user. """
        outputDict = {
            'saveTsPreAli': (self._registerTs, ['.preali']),
            'saveTsAli': (self._registerTs, ['.ali']),
            'saveReconsTomo': (self._registerReconsTomo, [])
        }

        if outputKey not in outputDict:
            raise Exception("Invalid key '%s' for registerOutput." % outputKey)

        # Call the specific function to create the output for this case
        func, args = outputDict[outputKey]
        outputName, output = func(*args)

        self._defineOutputs(**{outputName: output})
        self._defineSourceRelation(self.inputTiltSeries, output)
        self._store()

    def getPixSizeFromDimensions(self, outputDim):
        ih = ImageHandler()
        originalDim, _, _, _ = ih.getDimensions(self.inputTiltSeries.get().getFirstItem().getFileName())
        return self.inputTiltSeries.get().getSamplingRate() * round(originalDim/outputDim)
