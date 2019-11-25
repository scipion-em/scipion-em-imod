# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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
import pyworkflow.em as pyem
import pyworkflow.protocol.params as params
import tomo.objects as tomoObj

from tomo.objects import Landmark
from tomo.protocols import ProtTomoBase
from tomo.convert import writeTiStack


class ProtFiducialModel(pyem.EMProtocol, ProtTomoBase):
    """
    Construction of a fiducial model based on the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/etomoTutorial.html
    """

    _label = 'fiducial model'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-Series')

        form.addParam('twoSurfaces',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=0,
                      label='Find on two surfaces',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help="Track fiducials differentiating in which side of the sample are located.")

        form.addParam('numberFiducial',
                      params.IntParam,
                      label='Number of fiducials',
                      default='25',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Number of fiducials to be tracked for alignment.")

        form.addParam('rotationAngle',
                      params.FloatParam,
                      label='Tilt rotation angle (deg)',
                      default='0.0',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Angle from the vertical to the tilt axis in raw images.")

        form.addParam('importTrackFile',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Import track file',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Import a customized track file for the execution of the program.")

        groupTrackFile = form.addGroup('Track file',
                                       expertLevel=params.LEVEL_ADVANCED,
                                       condition='importTrackFile==0')

        groupTrackFile.addParam('trackFilePath',
                                params.PathParam,
                                label='File path',
                                expertLevel=params.LEVEL_ADVANCED,
                                help="Customized file path location.")

        form.addParam('computeAlignment',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Generate interpolated tilt-series',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Generate and save the interpolated tilt-series applying the'
                           'obtained transformation matrices.')

        groupInterpolation = form.addGroup('Interpolated tilt-series',
                                           condition='computeAlignment==0')

        groupInterpolation.addParam('binning', params.FloatParam,
                                    default=1.0,
                                    label='Binning',
                                    help='Binning to be applied to the interpolated tilt-series. '
                                         'Must be a integer bigger than 1')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('generateTrackComStep')
        self._insertFunctionStep('generateFiducialSeedStep')
        self._insertFunctionStep('generateFiducialModelStep')
        self._insertFunctionStep('translateFiducialPointModelStep')
        self._insertFunctionStep('computeFiducialAlignmentStep')
        self._insertFunctionStep('mergeAlignmentsStep')
        if self.computeAlignment.get() == 0:
            self._insertFunctionStep('computeInterpolatedStackStep')
        self._insertFunctionStep('_createOutputStep')

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)
            prefix = os.path.join(workingFolder, tsId)
            pw.utils.makePath(workingFolder)

            tiList = [ti.clone() for ti in ts]
            tiList.sort(key=lambda ti: ti.getTiltAngle())
            tiList.reverse()

            writeTiStack(tiList,
                         outputStackFn=prefix + '.st',
                         outputTltFn=prefix + '.rawtlt')

    def generateTrackComStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            paramsDict = {
                          'tsId': ts.getTsId(),
                          'rotationAngle': self.rotationAngle.get()
                          }
            self.translateTrackCom(ts.getTsId(), paramsDict)

    def generateFiducialSeedStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)
            paramsAutofidseed = {
                'trackCommandFile': '%s_track.com' % tsId,
                'minSpacing': 0.85,
                'peakStorageFraction': 1.0,
                'RotationAngle': self.rotationAngle.get(),
                'targetNumberOfBeads': self.numberFiducial.get()
            }

            argsAutofidseed = "-TrackCommandFile %(trackCommandFile)s " \
                              "-MinSpacing %(minSpacing)f " \
                              "-PeakStorageFraction %(peakStorageFraction)f " \
                              "-TargetNumberOfBeads %(targetNumberOfBeads)d"

            if self.twoSurfaces.get() == 0:
                argsAutofidseed += " -TwoSurfaces "

            self.runJob('autofidseed', argsAutofidseed % paramsAutofidseed, cwd=workingFolder)

    def generateFiducialModelStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)
            paramsBeadtrack = {
                'inputSeedModel': '%s.seed' % tsId,
                'outputModel': '%s.fid' % tsId,
                'imageFile': '%s.st' % tsId,
                'imagesAreBinned': 1,
                'tiltFile': '%s.rawtlt' %tsId,
                'tiltDefaultGrouping': 7,
                'magDefaultGrouping': 5,
                'rotDefaultGrouping': 1,
                'minViewsForTiltalign': 4,
                'beadDiameter': 4.95,
                'fillGaps': 1,
                'maxGapSize': 5,
                'minTiltRangeToFindAxis': 10.0,
                'minTiltRangeToFindAngles': 20.0,
                'boxSizeXandY': '32,32',
                'roundsOfTracking': 2,
                'localAreaTracking': 1,
                'localAreaTargetSize': 1000,
                'minBeadsInArea': 8,
                'minOverlapBeads': 5,
                'maxBeadsToAverage': 4,
                'sobelFilterCentering': 1,
                'pointsToFitMaxAndMin': '7,3',
                'densityRescueFractionAndSD': '0.6,1.0',
                'distanceRescueCriterion': 10.0,
                'rescueRelaxationDensityAndDistance': '0.7,0.9',
                'postFitRescueResidual': 2.5,
                'densityRelaxationPostFit': 0.9,
                'maxRescueDistance': 2.5,
                'residualsToAnalyzeMaxAndMin': '9,5',
                'deletionCriterionMinAndSD': '0.04,2.0'
            }

            argsBeadtrack = "-InputSeedModel %(inputSeedModel)s " \
                            "-OutputModel %(outputModel)s " \
                            "-ImageFile %(imageFile)s " \
                            "-ImagesAreBinned %(imagesAreBinned)d " \
                            "-TiltFile %(tiltFile)s " \
                            "-TiltDefaultGrouping %(tiltDefaultGrouping)d " \
                            "-MagDefaultGrouping %(magDefaultGrouping)d " \
                            "-RotDefaultGrouping %(rotDefaultGrouping)d " \
                            "-MinViewsForTiltalign %(minViewsForTiltalign)d " \
                            "-BeadDiameter %(beadDiameter)f " \
                            "-FillGaps %(fillGaps)d " \
                            "-MaxGapSize %(maxGapSize)d " \
                            "-MinTiltRangeToFindAxis %(minTiltRangeToFindAxis)f " \
                            "-MinTiltRangeToFindAngles %(minTiltRangeToFindAngles)f " \
                            "-BoxSizeXandY %(boxSizeXandY)s " \
                            "-RoundsOfTracking %(roundsOfTracking)d " \
                            "-LocalAreaTracking %(localAreaTracking)d " \
                            "-LocalAreaTargetSize %(localAreaTargetSize)d " \
                            "-MinBeadsInArea %(minBeadsInArea)d " \
                            "-MinOverlapBeads %(minOverlapBeads)d " \
                            "-MaxBeadsToAverage %(maxBeadsToAverage)d " \
                            "-SobelFilterCentering %(sobelFilterCentering)d " \
                            "-PointsToFitMaxAndMin %(pointsToFitMaxAndMin)s " \
                            "-DensityRescueFractionAndSD %(densityRescueFractionAndSD)s " \
                            "-DistanceRescueCriterion %(distanceRescueCriterion)f " \
                            "-RescueRelaxationDensityAndDistance %(rescueRelaxationDensityAndDistance)s " \
                            "-PostFitRescueResidual %(postFitRescueResidual)f " \
                            "-DensityRelaxationPostFit %(densityRelaxationPostFit)f " \
                            "-MaxRescueDistance %(maxRescueDistance)f " \
                            "-ResidualsToAnalyzeMaxAndMin %(residualsToAnalyzeMaxAndMin)s " \
                            "DeletionCriterionMinAndSD %(deletionCriterionMinAndSD)s"

            self.runJob('beadtrack', argsBeadtrack % paramsBeadtrack, cwd=workingFolder)

    def translateFiducialPointModelStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)

            paramsPoint2Model = {
                                'inputFile': '%s.fid' % tsId,
                                'outputFile': '%s_fid.txt' % tsId
                                }
            argsPoint2Model = "-InputFile %(inputFile)s " \
                              "-OutputFile %(outputFile)s"

            self.runJob('model2point', argsPoint2Model % paramsPoint2Model, cwd=workingFolder)

    def computeFiducialAlignmentStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)
            paramsTiltAlign = {
                'modelFile': '%s.fid' % tsId,
                'imageFile': '%s.st' %tsId,
                'imagesAreBinned': 1,
                'outputModelFile': '%s.3dmod' %tsId,
                'outputResidualFile': '%s.resid' %tsId,
                'outputFidXYZFile': '%sfid.xyz' %tsId,
                'outputTiltFile': '%s.tlt' %tsId,
                'outputXAxisTiltFile': '%s.xtilt' %tsId,
                'outputTransformFile': '%s.tltxf' %tsId,
                'outputFilledInModel': '%s_nogaps.fid' %tsId,
                'rotationAngle': self.rotationAngle.get(),
                'tiltFile': '%s.rawtlt' %tsId,
                'angleOffset': 0.0,
                'rotOption': 1,
                'rotDefaultGrouping': 5,
                'tiltOption': 5,
                'tiltDefaultGrouping': 5,
                'magReferenceView': 1,
                'magOption': 1,
                'magDefaultGrouping': 4,
                'xStretchOption': 0,
                'skewOption': 0,
                'xStretchDefaultGrouping': 7,
                'skewDefaultGrouping': 11,
                'beamTiltOption': 0,
                'xTiltOption': 0,
                'xTiltDefaultGrouping': 2000,
                'residualReportCriterion': 3.0,
                'surfacesToAnalyze': 2,
                'metroFactor': 0.25,
                'maximumCycles': 1000,
                'kFactorScaling': 1.0,
                'noSeparateTiltGroups': 1,
                'axisZShift': 0.0,
                'shiftZFromOriginal': 1,
                'localAlignments': 0,
                'outputLocalFile': '%slocal.xf' %tsId,
                'targetPatchSizeXandY': '700,700',
                'minSizeOrOverlapXandY':	'0.5,0.5',
                'minFidsTotalAndEachSurface': '8,3',
                'fixXYZCoordinates': 0,
                'localOutputOptions': '1,0,1',
                'localRotOption': 3,
                'localRotDefaultGrouping': 6,
                'localTiltOption': 5,
                'localTiltDefaultGrouping': 6,
                'localMagReferenceView': 1,
                'localMagOption': 3,
                'localMagDefaultGrouping': 7,
                'localXStretchOption': 0,
                'localXStretchDefaultGrouping': 7,
                'localSkewOption': 0,
                'localSkewDefaultGrouping': 11,
            }

            argsTiltAlign = "-ModelFile %(modelFile)s " \
                            "-ImageFile %(imageFile)s " \
                            "-ImagesAreBinned %(imagesAreBinned)d " \
                            "-OutputModelFile %(outputModelFile)s " \
                            "-OutputResidualFile %(outputResidualFile)s " \
                            "-OutputFidXYZFile %(outputFidXYZFile)s " \
                            "-OutputTiltFile %(outputTiltFile)s " \
                            "-OutputXAxisTiltFile %(outputXAxisTiltFile)s " \
                            "-OutputTransformFile %(outputTransformFile)s " \
                            "-OutputFilledInModel %(outputFilledInModel)s " \
                            "-RotationAngle %(rotationAngle)f " \
                            "-TiltFile %(tiltFile)s " \
                            "-AngleOffset %(angleOffset)f " \
                            "-RotOption %(rotOption)d " \
                            "-RotDefaultGrouping %(rotDefaultGrouping)d " \
                            "-TiltOption %(tiltOption)d " \
                            "-TiltDefaultGrouping %(tiltDefaultGrouping)d " \
                            "-MagReferenceView %(magReferenceView)d " \
                            "-MagOption %(magOption)d " \
                            "-MagDefaultGrouping %(magDefaultGrouping)d " \
                            "-XStretchOption %(xStretchOption)d " \
                            "-SkewOption %(skewOption)d " \
                            "-XStretchDefaultGrouping %(xStretchDefaultGrouping)d " \
                            "-BeamTiltOption %(beamTiltOption)d " \
                            "-XTiltOption %(xTiltOption)d " \
                            "-XTiltDefaultGrouping %(xTiltDefaultGrouping)d " \
                            "-ResidualReportCriterion %(residualReportCriterion)f " \
                            "-SurfacesToAnalyze %(surfacesToAnalyze)d " \
                            "-MetroFactor %(metroFactor)f " \
                            "-MaximumCycles %(maximumCycles)d " \
                            "-KFactorScaling %(kFactorScaling)f " \
                            "-NoSeparateTiltGroups %(noSeparateTiltGroups)d " \
                            "-AxisZShift %(axisZShift)f " \
                            "-ShiftZFromOriginal %(shiftZFromOriginal)d " \
                            "-LocalAlignments %(localAlignments)d " \
                            "-OutputLocalFile %(outputLocalFile)s " \
                            "-TargetPatchSizeXandY %(targetPatchSizeXandY)s " \
                            "-MinSizeOrOverlapXandY %(minSizeOrOverlapXandY)s " \
                            "-MinFidsTotalAndEachSurface %(minFidsTotalAndEachSurface)s " \
                            "-FixXYZCoordinates %(fixXYZCoordinates)d " \
                            "-LocalOutputOptions %(localOutputOptions)s " \
                            "-LocalRotOption %(localRotOption)d " \
                            "-LocalRotDefaultGrouping %(localRotDefaultGrouping)d " \
                            "-LocalTiltOption %(localTiltOption)d " \
                            "-LocalTiltDefaultGrouping %(localTiltDefaultGrouping)d " \
                            "-LocalMagReferenceView %(localMagReferenceView)d " \
                            "-LocalMagOption %(localMagOption)d " \
                            "-LocalMagDefaultGrouping %(localMagDefaultGrouping)d " \
                            "-LocalXStretchOption %(localXStretchOption)d " \
                            "-LocalXStretchDefaultGrouping %(localXStretchDefaultGrouping)s " \
                            "-LocalSkewOption %(localSkewOption)d " \
                            "-LocalSkewDefaultGrouping %(localSkewDefaultGrouping)d"

            self.runJob('tiltalign', argsTiltAlign % paramsTiltAlign, cwd=workingFolder)

    def mergeAlignmentsStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            matrix = []
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)
            for ti in ts:
                applyTransform = True
                if ti.hasTransform():
                    transform = ti.getTransform()
                    matrix.append(transform.getMatrixAsList())
                else:
                    applyTransform = False
                    break
            if applyTransform == True:
                self.createTransformFile(tsId, matrix)
                paramsXfproduct = {
                    'inputFile1': '%s.prexg' % tsId,
                    'inputFile2': '%s.tltxf' % tsId,
                    'outputFile': '%s_fid.xf' % tsId
                }

                argsXfproduct = "-InputFile1 %(inputFile1)s " \
                                "-InputFile2 %(inputFile2)s " \
                                "-OutputFile %(outputFile)s"

                self.runJob('xfproduct', argsXfproduct % paramsXfproduct, cwd=workingFolder)

            else:
                fileIn = os.path.join(workingFolder, '%s.tltxf' % tsId)
                fileOut = os.path.join(workingFolder, '%s_fid.xf' % tsId)
                print(fileIn)
                print(fileOut)
                os.rename(fileIn, fileOut)

    def computeInterpolatedStackStep(self):
        outputInterpolatedSetOfTiltSeries = self.getOutputInterpolatedSetOfTiltSeries()

        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            outputInterpolatedSetOfTiltSeries.append(newTs)
            workingFolder = self._getExtraPath(tsId)
            paramsAlginment = {
                'input': "%s.st" % tsId,
                'output': '%s_fidali.st' % tsId,
                'xform': "%s_fid.xf" % tsId,
                'bin': int(self.binning.get()),
                'mode': 0,
                'float': 2,
                'imagebinned': 1.0}
            argsAlignment = "-input %(input)s " \
                            "-output %(output)s " \
                            "-xform %(xform)s " \
                            "-bin %(bin)d " \
                            "-mode %(mode)s " \
                            "-float %(float)s " \
                            "-imagebinned %(imagebinned)s"

            self.runJob('newstack', argsAlignment % paramsAlginment, cwd=workingFolder)
            for index, tiltImage in enumerate(ts):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)
                newTi.setLocation(index + 1, (os.path.join(workingFolder, '%s_fidali.st' % tsId)))
                if self.binning > 1:
                    newTi.setSamplingRate(tiltImage.getSamplingRate() * int(self.binning.get()))
                newTs.append(newTi)
            if self.binning > 1:
                newTs.setSamplingRate(ts.getSamplingRate() * int(self.binning.get()))
            newTs.write()
            outputInterpolatedSetOfTiltSeries.update(newTs)  # update items and size info
            outputInterpolatedSetOfTiltSeries.write()
        self._store()

    def _createOutputStep(self):
        self.newSetOfLandmarks = self._createSetOfLandmarks()
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            fiducialList = self.parseFiducialModelFile(tsId)
            prevTiltIm = 0
            chainId = 0

            for index, fiducial in enumerate(fiducialList):
                if int(fiducial[2]) <= prevTiltIm:
                    chainId += 1
                prevTiltIm = int(fiducial[2])
                landmark = Landmark(fiducial[0], fiducial[1], fiducial[2], chainId, tsId)
                self.newSetOfLandmarks.append(landmark)
        self._defineOutputs(outputFiducialModel=self.newSetOfLandmarks)
        self._defineSourceRelation(self.inputSetOfTiltSeries, self.newSetOfLandmarks)

    # --------------------------- UTILS functions ----------------------------
    def createTransformFile(self, tsId, matrix):
        matrixFile = self._getExtraPath('%s/%s.prexg' % (tsId, tsId))
        with open(matrixFile, "w") as fOut:
            for line in matrix:
                fOut.write('%12.7f' % line[0])
                fOut.write('%12.7f' % line[3])
                fOut.write('%12.7f' % line[1])
                fOut.write('%12.7f' % line[4])
                fOut.write('%12.3f' % line[2])
                fOut.write('%12.3f' % line[5])
                fOut.write('\n')

    def translateTrackCom(self, tsId, paramsDict):
        trackFileName = tsId + "/" + tsId + "_track.com"
        trackFilePath = os.path.join(self._getExtraPath(), trackFileName)

        if self.importTrackFile == 0:
            with open(self.trackFilePath.get(), 'r') as f:
                userTrackFile = f.read()
            with open(trackFilePath, 'w') as f:
                f.write(userTrackFile)
        else:
            template = """# Command file for running BEADTRACK
#
####CreatedVersion####4.9.12
#
# For beads lighter than background, add a line with "LightBeads"
#
# To restrict tilt alignment to a subset of views, add a line with:
# "MaxViewsInAlign #_of_views"
#
# To exclude views, add a line "SkipViews view_list"; with the list of views
#
# To specify sets of views to be grouped separately in automapping, add a line
# "SeparateGroup view_list" with the list of views, one line per group
#
$beadtrack -StandardInput
ImageFile	%(tsId)s.st
ImagesAreBinned	1
InputSeedModel	%(tsId)s.seed
OutputModel	%(tsId)s.fid
RotationAngle	%(rotationAngle)f
TiltFile	%(tsId)s.rawtlt
TiltDefaultGrouping	7
MagDefaultGrouping	5
RotDefaultGrouping	1
BeadDiameter	4.95
FillGaps
MaxGapSize	5
RoundsOfTracking	2
#
# Set this to 1 to track in local areas
LocalAreaTracking	1
LocalAreaTargetSize	1000
MinBeadsInArea	8
MinOverlapBeads	5
#
# CONTROL PARAMETERS FOR EXPERTS, EXPERIMENTATION, OR SPECIAL CASES
#
# minimum range of tilt angles for finding axis and for finding tilts
MinViewsForTiltalign	4
MinTiltRangeToFindAxis	10.0
MinTiltRangeToFindAngles	20.0
BoxSizeXandY	32,32
MaxBeadsToAverage	4
# points and minimum for extrapolation
PointsToFitMaxAndMin	7,3
# fraction of mean, and # of SD below mean: density criterion for rescue
DensityRescueFractionAndSD	0.6,1.0
# distance criterion for rescue
DistanceRescueCriterion	10.0
# relaxation of criterion for density and distance rescues
RescueRelaxationDensityAndDistance	0.7,0.9
# distance for rescue after fit
PostFitRescueResidual	2.5
# relaxation of density criterion, maximum radius to search
DensityRelaxationPostFit	0.9
MaxRescueDistance	2.5
# Max and min residual changes to use to get mean and SD change
ResidualsToAnalyzeMaxAndMin	9,5
# minimum residual difference, criterion # of sd's
DeletionCriterionMinAndSD	0.04,2.0
SobelFilterCentering
$if (-e ./savework) ./savework
"""
            with open(trackFilePath, 'w') as f:
                f.write(template % paramsDict)

    def parseFiducialModelFile(self, tsId):
        fiducialFileName = tsId + "/" + tsId + "_fid.txt"
        fiducialFilePath = os.path.join(self._getExtraPath(), fiducialFileName)
        fiducialList = []
        with open(fiducialFilePath) as f:
            fiducialText = f.read().splitlines()
            for line in fiducialText:
                vector = line.split()
                fiducialList.append(vector)
        return fiducialList

    def getOutputInterpolatedSetOfTiltSeries(self):
        if not hasattr(self, "outputInterpolatedSetOfTiltSeries"):
            outputInterpolatedSetOfTiltSeries = self._createSetOfTiltSeries()
            outputInterpolatedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputInterpolatedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            if self.binning > 1:
                samplingRate = self.inputSetOfTiltSeries.get().getSamplingRate()
                samplingRate *= self.binning.get()
                outputInterpolatedSetOfTiltSeries.setSamplingRate(samplingRate)
            self._defineOutputs(outputInterpolatedSetOfTiltSeries=outputInterpolatedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputInterpolatedSetOfTiltSeries)
        return self.outputInterpolatedSetOfTiltSeries

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'newSetOfLandmarks'):
            summary.append("Input Tilt-Series: %d.\nFiducial models generated: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.newSetOfLandmarks.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'newSetOfLandmarks'):
            methods.append("The fiducial model has been computed for %d "
                           "Tilt-series using the IMOD procedure.\n"
                           % (self.outputSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods