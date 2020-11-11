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
import imod.utils as utils
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.objects import Transform
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from pyworkflow.object import Set
from tomo.objects import LandmarkModel
from tomo.protocols import ProtTomoBase
from imod import Plugin
from pwem.emlib.image import ImageHandler


class ProtImodFiducialAlignment(EMProtocol, ProtTomoBase):
    """
    Construction of a fiducial model and alignment of tilt-series based on the IMOD procedure.
    More info:
        https://bio3D.colorado.edu/imod/doc/etomoTutorial.html
    """

    _label = 'fiducial alignment'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series.')

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
                      label='Fiducial diameter (nm)',
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

        groupInterpolation.addParam('binning',
                                    params.FloatParam,
                                    default=1.0,
                                    label='Binning',
                                    help='Binning to be applied to the interpolated tilt-series in IMOD convention. '
                                         'Images will be binned by the given factor. Must be an integer bigger than 1')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('generateTrackComStep', ts.getObjId())
            self._insertFunctionStep('generateFiducialSeedStep', ts.getObjId())
            self._insertFunctionStep('generateFiducialModelStep', ts.getObjId())
            self._insertFunctionStep('computeFiducialAlignmentStep', ts.getObjId())
            self._insertFunctionStep('translateFiducialPointModelStep', ts.getObjId())
            self._insertFunctionStep('computeOutputStackStep', ts.getObjId())
            if self.computeAlignment.get() == 0:
                self._insertFunctionStep('computeOutputInterpolatedStackStep', ts.getObjId())
            self._insertFunctionStep('computeOutputModelsStep', ts.getObjId())
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        """Apply the transformation form the input tilt-series"""
        outputTsFileName = os.path.join(tmpPrefix, ts.getFirstItem().parseFileName())
        ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, ts.getFirstItem().parseFileName(extension=".tlt"))
        ts.generateTltFile(angleFilePath)

        """"Link to input tilt-series (needed for fiducial model viewer)"""
        # TODO: there is no need to come from a prealigned stack
        inputTS = os.path.join(extraPrefix, ts.getFirstItem().parseFileName())
        if ts.getFirstItem().hasTransform():
            path.copyFile(outputTsFileName, inputTS)

        else:
            path.createLink(ts.getFirstItem().getLocation()[1], inputTS)

    def generateTrackComStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsDict = {
            'imageFile': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName()),
            'inputSeedModel': os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".seed")),
            'outputModel': os.path.join(extraPrefix, ts.getFirstItem().parseFileName(suffix="_gaps", extension=".fid")),
            'tiltFile': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName(extension=".tlt")),
            'rotationAngle': self.rotationAngle.get(),
            'fiducialDiameter': self.fiducialDiameter.get()
        }

        self.translateTrackCom(ts, paramsDict)

    def generateFiducialSeedStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        paramsAutofidseed = {
            'trackCommandFile': os.path.join(extraPrefix,
                                             ts.getFirstItem().parseFileName(suffix="_track", extension=".com")),
            'minSpacing': 0.85,
            'peakStorageFraction': 1.0,
            'RotationAngle': self.rotationAngle.get(),
            'targetNumberOfBeads': self.numberFiducial.get()
        }
        argsAutofidseed = "-TrackCommandFile %(trackCommandFile)s " \
                          "-MinSpacing %(minSpacing)f " \
                          "-PeakStorageFraction %(peakStorageFraction)f " \
                          "-TargetNumberOfBeads %(targetNumberOfBeads)d "
        if self.twoSurfaces.get() == 0:
            argsAutofidseed += " -TwoSurfaces"

        Plugin.runImod(self, 'autofidseed', argsAutofidseed % paramsAutofidseed)

        autofidseedDirPath = os.path.join(self._getExtraPath(tsId), "autofidseed.dir")
        path.makePath(autofidseedDirPath)
        path.moveTree("autofidseed.dir", autofidseedDirPath)
        path.moveFile("autofidseed.info", self._getExtraPath(tsId))

    def generateFiducialModelStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        paramsBeadtrack = {
            'inputSeedModel': os.path.join(extraPrefix, ts.getFirstItem().parseFileName(extension=".seed")),
            'outputModel': os.path.join(extraPrefix, ts.getFirstItem().parseFileName(suffix="_gaps", extension=".fid")),
            'imageFile': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName()),
            'imagesAreBinned': 1,
            'tiltFile': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName(extension=".tlt")),
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

        Plugin.runImod(self, 'beadtrack', argsBeadtrack % paramsBeadtrack)

    def computeFiducialAlignmentStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        paramsTiltAlign = {
            'modelFile': os.path.join(extraPrefix, ts.getFirstItem().parseFileName(suffix="_gaps", extension=".fid")),
            'imageFile': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName()),
            'imagesAreBinned': 1,
            'outputModelFile': os.path.join(extraPrefix,
                                            ts.getFirstItem().parseFileName(suffix="_fidxyz", extension=".mod")),
            'outputResidualFile': os.path.join(extraPrefix,
                                               ts.getFirstItem().parseFileName(suffix="_resid", extension=".txt")),
            'outputFidXYZFile': os.path.join(extraPrefix,
                                             ts.getFirstItem().parseFileName(suffix="_fid", extension=".xyz")),
            'outputTiltFile': os.path.join(extraPrefix,
                                           ts.getFirstItem().parseFileName(suffix="_interpolated", extension=".tlt")),
            'outputTransformFile': os.path.join(extraPrefix,
                                                ts.getFirstItem().parseFileName(suffix="_fid", extension=".xf")),
            'outputFilledInModel': os.path.join(extraPrefix,
                                                ts.getFirstItem().parseFileName(suffix="_noGaps", extension=".fid")),
            'rotationAngle': self.rotationAngle.get(),
            'tiltFile': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName(extension=".tlt")),
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
            'outputLocalFile': os.path.join(extraPrefix,
                                            ts.getFirstItem().parseFileName(suffix="_local", extension=".xf")),
            'targetPatchSizeXandY': '700,700',
            'minSizeOrOverlapXandY': '0.5,0.5',
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

        Plugin.runImod(self, 'tiltalign', argsTiltAlign % paramsTiltAlign)

    def translateFiducialPointModelStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        paramsGapPoint2Model = {
            'inputFile': os.path.join(extraPrefix,
                                      ts.getFirstItem().parseFileName(suffix="_gaps", extension=".fid")),
            'outputFile': os.path.join(extraPrefix,
                                       ts.getFirstItem().parseFileName(suffix="_gaps_fid", extension=".txt"))
        }
        argsGapPoint2Model = "-InputFile %(inputFile)s " \
                             "-OutputFile %(outputFile)s"

        Plugin.runImod(self, 'model2point', argsGapPoint2Model % paramsGapPoint2Model)

        paramsNoGapPoint2Model = {
            'inputFile': os.path.join(extraPrefix,
                                      ts.getFirstItem().parseFileName(suffix="_noGaps", extension=".fid")),
            'outputFile': os.path.join(extraPrefix,
                                       ts.getFirstItem().parseFileName(suffix="_noGaps_fid", extension=".txt"))
        }
        argsNoGapPoint2Model = "-InputFile %(inputFile)s " \
                               "-OutputFile %(outputFile)s"

        Plugin.runImod(self, 'model2point', argsNoGapPoint2Model % paramsNoGapPoint2Model)

    def computeOutputStackStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        tltFilePath = os.path.join(extraPrefix,
                                   ts.getFirstItem().parseFileName(suffix="_interpolated", extension=".tlt"))
        tltList = utils.formatAngleList(tltFilePath)

        transformationMatricesFilePath = os.path.join(
                                            extraPrefix,
                                            ts.getFirstItem().parseFileName(suffix="_fid", extension=".xf"))
        newTransformationMatricesList = utils.formatTransformationMatrix(transformationMatricesFilePath)

        outputSetOfTiltSeries = self.getOutputSetOfTiltSeries()
        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputSetOfTiltSeries.append(newTs)

        for index, ti in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(ti, copyId=True)
            newTi.setLocation(ti.getLocation())
            newTi.setTiltAngle(float(tltList[index]))

            if ti.hasTransform():
                transform = Transform()
                previousTransform = ti.getTransform().getMatrix()
                newTransform = newTransformationMatricesList[:, :, index]
                previousTransformArray = np.array(previousTransform)
                newTransformArray = np.array(newTransform)
                outputTransformMatrix = np.matmul(previousTransformArray, newTransformArray)
                transform.setMatrix(outputTransformMatrix)
                newTi.setTransform(transform)
            else:
                transform = Transform()
                newTransform = newTransformationMatricesList[:, :, index]
                newTransformArray = np.array(newTransform)
                transform.setMatrix(newTransformArray)
                newTi.setTransform(transform)

            newTs.append(newTi)

        newTs.write(properties=False)

        outputSetOfTiltSeries.update(newTs)
        outputSetOfTiltSeries.write()

        self._store()

    def computeOutputInterpolatedStackStep(self, tsObjId):
        outputInterpolatedSetOfTiltSeries = self.getOutputInterpolatedSetOfTiltSeries()

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsAlignment = {
            'input': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName()),
            'output': os.path.join(extraPrefix, ts.getFirstItem().parseFileName()),
            'xform': os.path.join(extraPrefix, ts.getFirstItem().parseFileName(suffix="_fid", extension=".xf")),
            'bin': int(self.binning.get()),
            'imagebinned': 1.0}
        argsAlignment = "-input %(input)s " \
                        "-output %(output)s " \
                        "-xform %(xform)s " \
                        "-bin %(bin)d " \
                        "-imagebinned %(imagebinned)s"

        Plugin.runImod(self, 'newstack', argsAlignment % paramsAlignment)

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputInterpolatedSetOfTiltSeries.append(newTs)

        tltFilePath = os.path.join(extraPrefix,
                                   ts.getFirstItem().parseFileName(suffix="_interpolated", extension=".tlt"))

        tltList = utils.formatAngleList(tltFilePath)

        if self.binning > 1:
            newTs.setSamplingRate(ts.getSamplingRate() * int(self.binning.get()))

        for index, ti in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(ti, copyId=True)
            newTi.setLocation(index + 1, os.path.join(extraPrefix, ti.parseFileName()))
            newTi.setTiltAngle(float(tltList[index]))
            if self.binning > 1:
                newTi.setSamplingRate(ti.getSamplingRate() * int(self.binning.get()))
            newTs.append(newTi)

        ih = ImageHandler()
        x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
        newTs.setDim((x, y, z))
        newTs.write(properties=False)

        outputInterpolatedSetOfTiltSeries.update(newTs)
        outputInterpolatedSetOfTiltSeries.updateDim()
        outputInterpolatedSetOfTiltSeries.write()
        self._store()

    def computeOutputModelsStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        """Create the output set of landmark models with gaps"""
        outputSetOfLandmarkModelsGaps = self.getOutputFiducialModelGaps()

        fiducialModelGapPath = os.path.join(extraPrefix,
                                            ts.getFirstItem().parseFileName(suffix="_gaps", extension=".fid"))

        landmarkModelGapsResidPath = os.path.join(extraPrefix,
                                                  ts.getFirstItem().parseFileName(suffix="_resid", extension=".txt"))
        fiducialGapResidList = utils.formatFiducialResidList(landmarkModelGapsResidPath)

        landmarkModelGaps = LandmarkModel(tsId, landmarkModelGapsResidPath, fiducialModelGapPath)

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
        outputSetOfLandmarkModelsNoGaps = self.getOutputFiducialModelNoGaps()

        fiducialNoGapFilePath = os.path.join(extraPrefix,
                                             ts.getFirstItem().parseFileName(suffix="_noGaps_fid", extension=".txt"))
        fiducialNoGapList = utils.formatFiducialList(fiducialNoGapFilePath)

        fiducialModelNoGapPath = os.path.join(extraPrefix,
                                              ts.getFirstItem().parseFileName(suffix="_noGaps", extension=".fid"))

        landmarkModelNoGapsFilePath = os.path.join(extraPrefix,
                                                   ts.getFirstItem().parseFileName(suffix="_noGaps", extension=".sfid"))

        landmarkModelNoGapsResidPath = os.path.join(extraPrefix,
                                                    ts.getFirstItem().parseFileName(suffix="_resid", extension=".txt"))
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

        """Create the output set of coordinates 3D from the fiducials in the tilt series"""
        outputSetOfCoordinates3D = self.getOutputSetOfCoordinates3Ds()

        coordFilePath = os.path.join(extraPrefix,
                                     ts.getFirstItem().parseFileName(suffix="_fid", extension=".xyz"))

        xDim = ts.getFirstItem().getXDim()
        yDim = ts.getFirstItem().getYDim()
        coordList = utils.format3DCoordinatesList(coordFilePath, xDim, yDim)
        for element in coordList:
            newCoord3D = tomoObj.Coordinate3D(x=element[0],
                                              y=element[1],
                                              z=element[2])
            newCoord3D.setVolume(ts)
            newCoord3D.setVolId(tsObjId)
            outputSetOfCoordinates3D.append(newCoord3D)
            outputSetOfCoordinates3D.update(newCoord3D)
        outputSetOfCoordinates3D.write()
        self._store()

    def createOutputStep(self):
        self.getOutputSetOfTiltSeries().setStreamState(Set.STREAM_CLOSED)
        if self.computeAlignment.get() == 0:
            self.getOutputInterpolatedSetOfTiltSeries().setStreamState(Set.STREAM_CLOSED)
        self.getOutputFiducialModelGaps().setStreamState(Set.STREAM_CLOSED)
        self.getOutputFiducialModelNoGaps().setStreamState(Set.STREAM_CLOSED)
        self.getOutputSetOfCoordinates3Ds().setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def translateTrackCom(self, ts, paramsDict):
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        trackFilePath = os.path.join(extraPrefix,
                                     ts.getFirstItem().parseFileName(suffix="_track", extension=".com"))

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
ImageFile	%(imageFile)s
ImagesAreBinned	1
InputSeedModel	%(inputSeedModel)s
OutputModel	%(outputModel)s
RotationAngle	%(rotationAngle)f
TiltFile	%(tiltFile)s
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

    def getOutputSetOfTiltSeries(self):
        if hasattr(self, "outputSetOfTiltSeries"):
            self.outputSetOfTiltSeries.enableAppend()
        else:
            outputSetOfTiltSeries = self._createSetOfTiltSeries()
            outputSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfTiltSeries=outputSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfTiltSeries)
        return self.outputSetOfTiltSeries

    def getOutputInterpolatedSetOfTiltSeries(self):
        if hasattr(self, "outputInterpolatedSetOfTiltSeries"):
            self.outputInterpolatedSetOfTiltSeries.enableAppend()
        else:
            outputInterpolatedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Interpolated')
            outputInterpolatedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputInterpolatedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            if self.binning > 1:
                samplingRate = self.inputSetOfTiltSeries.get().getSamplingRate()
                samplingRate *= self.binning.get()
                outputInterpolatedSetOfTiltSeries.setSamplingRate(samplingRate)
            outputInterpolatedSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputInterpolatedSetOfTiltSeries=outputInterpolatedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputInterpolatedSetOfTiltSeries)
        return self.outputInterpolatedSetOfTiltSeries

    def getOutputFiducialModelGaps(self):
        if hasattr(self, "outputFiducialModelGaps"):
            self.outputFiducialModelGaps.enableAppend()
        else:
            outputFiducialModelGaps = self._createSetOfLandmarkModels(suffix='Gaps')
            outputFiducialModelGaps.copyInfo(self.inputSetOfTiltSeries.get())
            outputFiducialModelGaps.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputFiducialModelGaps=outputFiducialModelGaps)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputFiducialModelGaps)
        return self.outputFiducialModelGaps

    def getOutputFiducialModelNoGaps(self):
        if hasattr(self, "outputFiducialModelNoGaps"):
            self.outputFiducialModelNoGaps.enableAppend()
        else:
            outputFiducialModelNoGaps = self._createSetOfLandmarkModels(suffix='NoGaps')
            outputFiducialModelNoGaps.copyInfo(self.inputSetOfTiltSeries.get())
            outputFiducialModelNoGaps.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputFiducialModelNoGaps=outputFiducialModelNoGaps)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputFiducialModelNoGaps)
        return self.outputFiducialModelNoGaps

    def getOutputSetOfCoordinates3Ds(self):
        if hasattr(self, "outputSetOfCoordinates3D"):
            self.outputSetOfCoordinates3D.enableAppend()
        else:
            outputSetOfCoordinates3D = self._createSetOfCoordinates3D(volSet=self.getOutputSetOfTiltSeries(),
                                                                      suffix='LandmarkModel')
            outputSetOfCoordinates3D.setSamplingRate(self.inputSetOfTiltSeries.get().getSamplingRate())
            outputSetOfCoordinates3D.setPrecedents(self.inputSetOfTiltSeries)
            outputSetOfCoordinates3D.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfCoordinates3D=outputSetOfCoordinates3D)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfCoordinates3D)
        return self.outputSetOfCoordinates3D

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
