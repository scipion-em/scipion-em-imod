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
import numpy as np
import pyworkflow as pw
import pyworkflow.em as pyem
import pyworkflow.protocol.params as params
import pyworkflow.em.data as data
import tomo.objects as tomoObj
from tomo.objects import LandmarkModel
from tomo.protocols import ProtTomoBase
from tomo.convert import writeTiStack


class ProtFiducialModel(pyem.EMProtocol, ProtTomoBase):
    """
    Construction of a fiducial model based on the IMOD procedure.
    More info:
        https://bio3D.colorado.edu/imod/doc/etomoTutorial.html
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

        form.addParam('fiducialDiameter',
                      params.FloatParam,
                      label='Fiducial diameter',
                      default='4.95',
                      help="Fiducials diameter to be tracked for alignment.")

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
        self._insertFunctionStep('computeFiducialAlignmentStep')
        self._insertFunctionStep('translateFiducialPointModelStep')
        self._insertFunctionStep('computeOutputStackStep')
        if self.computeAlignment.get() == 0:
            self._insertFunctionStep('computeInterpolatedStackStep')
        self._insertFunctionStep('createOutputStep')
        self._insertFunctionStep('cleanDirectory')

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

            """Apply the transformation form the input tilt-series"""
            inputStack = prefix + '.st'
            transformedStack = prefix + '_transformed.st'
            newStack = True
            for index, ti in enumerate(ts):
                if ti.hasTransform():
                    ih = pyem.ImageHandler()
                    if newStack:
                        ih.createEmptyImage(fnOut=transformedStack,
                                            xDim=ti.getXDim(),
                                            yDim=ti.getYDim(),
                                            nDim=ts.getSize())
                        newStack = False
                    transform = ti.getTransform().getMatrix()
                    transformArray = np.array(transform)
                    ih.applyTransform(inputFile=str(index+1)+'@'+inputStack,
                                      outputFile=str(ts.getSize()-index)+'@'+transformedStack,
                                      transformMatrix=transformArray,
                                      shape=(ti.getXDim(), ti.getYDim()),
                                      borderAverage=True)
                else:
                    os.rename(inputStack, transformedStack)
                    break

    def generateTrackComStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            paramsDict = {
                          'tsId': ts.getTsId(),
                          'rotationAngle': self.rotationAngle.get(),
                          'fiducialDiameter': self.fiducialDiameter.get()
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
                argsAutofidseed += " -TwoSurfaces"

            self.runJob('autofidseed', argsAutofidseed % paramsAutofidseed, cwd=workingFolder)

    def generateFiducialModelStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)
            paramsBeadtrack = {
                'inputSeedModel': '%s.seed' % tsId,
                'outputModel': '%s.fid' % tsId,
                'imageFile': '%s_transformed.st' % tsId,
                'imagesAreBinned': 1,
                'tiltFile': '%s.rawtlt' %tsId,
                'tiltDefaultGrouping': 7,
                'magDefaultGrouping': 5,
                'rotDefaultGrouping': 1,
                'minViewsForTiltalign': 4,
                'beadDiameter': self.fiducialDiameter.get(),
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
                            "-DeletionCriterionMinAndSD %(deletionCriterionMinAndSD)s"

            self.runJob('beadtrack', argsBeadtrack % paramsBeadtrack, cwd=workingFolder)

    def translateFiducialPointModelStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)

            paramsGapPoint2Model = {
                                'inputFile': '%s.fid' % tsId,
                                'outputFile': '%s_fid.txt' % tsId
                                }
            argsGapPoint2Model = "-InputFile %(inputFile)s " \
                                 "-OutputFile %(outputFile)s"
            self.runJob('model2point', argsGapPoint2Model % paramsGapPoint2Model, cwd=workingFolder)

            paramsNoGapPoint2Model = {
                                'inputFile': '%s_noGaps.fid' % tsId,
                                'outputFile': '%s_noGaps_fid.txt' % tsId
                                }
            argsNoGapPoint2Model = "-InputFile %(inputFile)s " \
                                   "-OutputFile %(outputFile)s"
            self.runJob('model2point', argsNoGapPoint2Model % paramsNoGapPoint2Model, cwd=workingFolder)

    def computeFiducialAlignmentStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)
            paramsTiltAlign = {
                'modelFile': '%s.fid' % tsId,
                'imageFile': '%s_transformed.st' %tsId,
                'imagesAreBinned': 1,
                'outputModelFile': '%s_fidxyz.mod' %tsId,
                'outputResidualFile': '%s_resid.txt' %tsId,
                'outputFidXYZFile': '%s_fid.xyz' %tsId,
                'outputTiltFile': '%s_interpolated.tlt' %tsId,
                'outputTransformFile': '%s.fidxf' %tsId,
                'outputFilledInModel': '%s_noGaps.fid' %tsId,
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

    def computeOutputStackStep(self):
        outputSetOfTiltSeries = self.getOutputSetOfTiltSeries()
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            outputSetOfTiltSeries.append(newTs)

            tltFileName = tsId + "/" + tsId + "_interpolated.tlt"
            tltFilePath = os.path.join(self._getExtraPath(), tltFileName)
            tltList = self.parseAngleTltFile(tltFilePath)

            transformationMatricesFile = tsId + "/" + tsId + ".fidxf"
            transformationMatricesFilePath = os.path.join(self._getExtraPath(), transformationMatricesFile)
            newTransformationMatricesList = self.formatTransformationMatrix(transformationMatricesFilePath)

            for index, ti in enumerate(ts):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(ti, copyId=True)
                newTi.setLocation(index + 1, ti.getLocation()[1])
                newTi.setTiltAngle(float(tltList[index]))

                if ti.hasTransform():
                    transform = data.Transform()
                    previousTransform = ti.getTransform().getMatrix()
                    newTransform = newTransformationMatricesList[:, :, index]
                    previousTransformArray = np.array(previousTransform)
                    newTransformArray = np.array(newTransform)
                    outputTransformMatrix = np.matmul(previousTransformArray, newTransformArray)
                    transform.setMatrix(outputTransformMatrix)
                    newTi.setTransform(transform)

                newTs.append(newTi)

            newTs.write()
            outputSetOfTiltSeries.update(newTs)
            outputSetOfTiltSeries.write()
        self._store()

        for ts in outputSetOfTiltSeries:
            for ti in ts:
                print(ti.getTransform())
                print(ti.getTiltAngle())


    def computeInterpolatedStackStep(self):
        outputInterpolatedSetOfTiltSeries = self.getOutputInterpolatedSetOfTiltSeries()

        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            outputInterpolatedSetOfTiltSeries.append(newTs)

            workingFolder = self._getExtraPath(tsId)
            paramsAlginment = {
                'input': "%s_transformed.st" % tsId,
                'output': '%s_fidali.st' % tsId,
                'xform': "%s.fidxf" % tsId,
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

            tltFileName = tsId + "/" + tsId + "_interpolated.tlt"
            tltFilePath = os.path.join(self._getExtraPath(), tltFileName)
            tltList = self.parseAngleTltFile(tltFilePath)
            for index, ti in enumerate(ts):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(ti, copyId=True)
                newTi.setLocation(index + 1, (os.path.join(workingFolder, '%s_fidali.st' % tsId)))
                newTi.setTiltAngle(float(tltList[index]))
                if self.binning > 1:
                    newTi.setSamplingRate(ti.getSamplingRate() * int(self.binning.get()))
                newTs.append(newTi)
            if self.binning > 1:
                newTs.setSamplingRate(ts.getSamplingRate() * int(self.binning.get()))
            newTs.write()
            outputInterpolatedSetOfTiltSeries.update(newTs)
            outputInterpolatedSetOfTiltSeries.write()
        self._store()

    def createOutputStep(self):
        """Create the output set of landmark models with gaps"""
        self.newSetOfLandmarkModelsGaps = self._createSetOfLandmarkModels(suffix='Gaps')
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()

            fiducialGapFileName = tsId + "/" + tsId + "_fid.txt"
            fiducialGapFilePath = os.path.join(self._getExtraPath(), fiducialGapFileName)
            fiducialGapList = self.parseFiducialModelFile(fiducialGapFilePath)

            fiducialModelGapFileName = tsId + "/" + tsId + ".fid"
            fiducialModelGapPath = os.path.join(self._getExtraPath(), fiducialModelGapFileName)

            landmarkModelGapsFileName = tsId + "/" + tsId + "_gaps.sfid"
            landmarkModelGapsFilePath = os.path.join(self._getExtraPath(), landmarkModelGapsFileName)

            landmarkModel = LandmarkModel(tsId, landmarkModelGapsFilePath, fiducialModelGapPath)

            prevTiltIm = 0
            chainId = 0
            for index, fiducial in enumerate(fiducialGapList):
                if int(fiducial[2]) <= prevTiltIm:
                    chainId += 1
                prevTiltIm = int(fiducial[2])
                landmarkModel.addLandmark(fiducial[0], fiducial[1], fiducial[2], chainId)
            self.newSetOfLandmarkModelsGaps.append(landmarkModel)

        self._defineOutputs(outputFiducialModelGaps=self.newSetOfLandmarkModelsGaps)
        self._defineSourceRelation(self.inputSetOfTiltSeries, self.newSetOfLandmarkModelsGaps)

        """Create the output set of landmark models with no gaps"""
        self.newSetOfLandmarkModelsNoGaps = self._createSetOfLandmarkModels(suffix='NoGaps')
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()

            fiducialNoGapFileName = tsId + "/" + tsId + "_noGaps_fid.txt"
            fiducialNoGapFilePath = os.path.join(self._getExtraPath(), fiducialNoGapFileName)
            fiducialNoGapList = self.parseFiducialModelFile(fiducialNoGapFilePath)

            fiducialModelNoGapFileName = tsId + "/" + tsId + "_noGaps.fid"
            fiducialModelNoGapPath = os.path.join(self._getExtraPath(), fiducialModelNoGapFileName)

            landmarkModelNoGapsFileName = tsId + "/" + tsId + "_noGaps.sfid"
            landmarkModelNoGapsFilePath = os.path.join(self._getExtraPath(), landmarkModelNoGapsFileName)

            landmarkModel = LandmarkModel(tsId, landmarkModelNoGapsFilePath, fiducialModelNoGapPath)

            prevTiltIm = 0
            chainId = 0
            for index, fiducial in enumerate(fiducialNoGapList):
                if int(fiducial[2]) <= prevTiltIm:
                    chainId += 1
                prevTiltIm = int(fiducial[2])
                landmarkModel.addLandmark(fiducial[0], fiducial[1], fiducial[2], chainId)
            self.newSetOfLandmarkModelsNoGaps.append(landmarkModel)

        self._defineOutputs(outputFiducialModelNoGaps=self.newSetOfLandmarkModelsNoGaps)
        self._defineSourceRelation(self.inputSetOfTiltSeries, self.newSetOfLandmarkModelsNoGaps)

        """Create the output set of coordinates 3D from the fiducials in the tilt series"""
        self.newSetOfCoordinates3D = self._createSetOfCoordinates3D(volSet=self.getOutputSetOfTiltSeries(),
                                                                    suffix='FiducialModels')
        for index, ts in enumerate(self.getOutputSetOfTiltSeries()):
            tsId = ts.getTsId()
            coorFileName = tsId + "/" + tsId + "_fid.xyz"
            xDim = ts.getFirstItem().getXDim()
            yDim = ts.getFirstItem().getYDim()
            coorList = self.parse3DCoordinatesFile(coorFileName, xDim, yDim)
            print (coorList)
            for element in coorList:
                newCoor3D = tomoObj.Coordinate3D(x=element[0],
                                                 y=element[1],
                                                 z=element[2])
                newCoor3D.setVolId(index+1)
                newCoor3D.setVolName(tsId)
                self.newSetOfCoordinates3D.append(newCoor3D)

        self._defineOutputs(outputSetOfCoordinates3D=self.newSetOfCoordinates3D)
        self._defineSourceRelation(self.inputSetOfTiltSeries, self.newSetOfCoordinates3D)

    def cleanDirectory(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)
            if os.path.exists(os.path.join(workingFolder, "%s.st" % tsId)):
                os.remove(os.path.join(workingFolder, "%s.st" % tsId))
            os.remove(os.path.join(workingFolder, "%s.rawtlt" % tsId))
            os.remove(os.path.join(workingFolder, "%s_transformed.st" % tsId))

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
ImageFile	%(tsId)s_transformed.st
ImagesAreBinned	1
InputSeedModel	%(tsId)s.seed
OutputModel	%(tsId)s.fid
RotationAngle	%(rotationAngle)f
TiltFile	%(tsId)s.rawtlt
TiltDefaultGrouping	7
MagDefaultGrouping	5
RotDefaultGrouping	1
BeadDiameter	%(fiducialDiameter)f
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

    def parseFiducialModelFile(self, fiducialFilePath):
        fiducialList = []
        with open(fiducialFilePath) as f:
            fiducialText = f.read().splitlines()
            for line in fiducialText:
                vector = line.split()
                fiducialList.append(vector)
        return fiducialList

    def parseAngleTltFile(self, tltFilePath):
        angleList = []
        with open(tltFilePath) as f:
            tltText = f.read().splitlines()
            for line in tltText:
                angleList.append(float(line))
        angleList.reverse()
        return angleList

    def formatTransformationMatrix(self, matrixFile):
        with open(matrixFile, "r") as matrix:
            lines = matrix.readlines()
        numberLines = len(lines)
        frameMatrix = np.empty([3, 3, numberLines])
        i = 0
        for line in lines:
            values = line.split()
            frameMatrix[0, 0, i] = float(values[0])
            frameMatrix[1, 0, i] = float(values[1])
            frameMatrix[0, 1, i] = float(values[2])
            frameMatrix[1, 1, i] = float(values[3])
            frameMatrix[0, 2, i] = float(values[4])
            frameMatrix[1, 2, i] = float(values[5])
            frameMatrix[2, 0, i] = 0.0
            frameMatrix[2, 1, i] = 0.0
            frameMatrix[2, 2, i] = 1.0
            i += 1
        return frameMatrix

    def parse3DCoordinatesFile(self, coorFileName, xDim, yDim):
        coorList = []
        coorFilePath = os.path.join(self._getExtraPath(),coorFileName)
        with open(coorFilePath) as f:
            coorText = f.read().splitlines()
            for line in coorText:
                vector = line.split()
                coorList.append([float(vector[1]) - xDim/2, float(vector[2]) - yDim/2, float(vector[3])])
        return coorList

    def getOutputInterpolatedSetOfTiltSeries(self):
        if not hasattr(self, "outputInterpolatedSetOfTiltSeries"):
            outputInterpolatedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Interpolated')
            outputInterpolatedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputInterpolatedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            if self.binning > 1:
                samplingRate = self.inputSetOfTiltSeries.get().getSamplingRate()
                samplingRate *= self.binning.get()
                outputInterpolatedSetOfTiltSeries.setSamplingRate(samplingRate)
            self._defineOutputs(outputInterpolatedSetOfTiltSeries=outputInterpolatedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputInterpolatedSetOfTiltSeries)
        return self.outputInterpolatedSetOfTiltSeries

    def getOutputSetOfTiltSeries(self):
        if not hasattr(self, "outputSetOfTiltSeries"):
            outputSetOfTiltSeries = self._createSetOfTiltSeries()
            outputSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            self._defineOutputs(outputSetOfTiltSeries=outputSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfTiltSeries)
        return self.outputSetOfTiltSeries

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputFiducialModelGaps'):
            summary.append("Input Tilt-Series: %d.\nFiducial models generated presenting gaps: %d."
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputFiducialModelGaps.getSize()))

        if hasattr(self, 'outputFiducialModelNoGaps'):
            summary.append("Fiducial models generated with no gaps: %d."
                           % (self.outputFiducialModelNoGaps.getSize()))

        if hasattr(self, 'outputSetOfTiltSeries'):
            summary.append("Transformation matrices updated from the input Tilt-Series: %d."
                           % (self.outputSetOfTiltSeries.getSize()))

        if hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            summary.append("Interpolated Tilt-Series calculated: %d."
                           % (self.outputInterpolatedSetOfTiltSeries.getSize()))

        if hasattr(self, 'outputSetOfCoordinates3D'):
            summary.append("Fiducial 3D coordinates calculated for %d Tilt-series: %d."
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputSetOfCoordinates3D.getSize()))

        if not summary:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputFiducialModelGaps'):
            methods.append("The fiducial model (presenting gaps) has been computed for %d "
                           "Tilt-series using the IMOD procedure."
                           % (self.outputFiducialModelGaps.getSize()))

        if hasattr(self, 'outputFiducialModelNoGaps'):
            methods.append("The fiducial model (with no gaps) has been computed for %d "
                           "Tilt-series using the IMOD procedure."
                           % (self.outputFiducialModelNoGaps.getSize()))

        if hasattr(self, 'outputSetOfTiltSeries'):
            methods.append("The transformation matrices has been computed for %d "
                           "Tilt-series using the IMOD procedure."
                           % (self.outputSetOfTiltSeries.getSize()))

        if hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            methods.append("%d Tilt-Series have been interpolated using the IMOD procedure."
                           % (self.outputInterpolatedSetOfTiltSeries.getSize()))

        if hasattr(self, 'outputSetOfCoordinates3D'):
            methods.append("%d fiducial 3D coordinates have been calculated for %d Tilt-series."
                           % (self.outputSetOfCoordinates3D.getSize(),
                              self.inputSetOfTiltSeries.get().getSize()))

        if not methods:
            methods.append("Output classes not ready yet.")
        return methods
