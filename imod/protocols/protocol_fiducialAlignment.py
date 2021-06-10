# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
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
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.objects import Transform
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from pyworkflow.object import Set
from tomo.objects import LandmarkModel
from tomo.protocols import ProtTomoBase
import tomo.constants as constants
from imod import Plugin
from pwem.emlib.image import ImageHandler


class ProtImodFiducialAlignment(EMProtocol, ProtTomoBase):
    """
    Construction of a fiducial model and alignment of tilt-series based on the IMOD procedure.
    More info:
        https://bio3D.colorado.edu/imod/doc/etomoTutorial.html
    """

    _label = 'fiducial alignment'
    _devStatus = BETA

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
                      display=params.EnumParam.DISPLAY_HLIST,
                      help="Track fiducials differentiating in which side of the sample are located.")

        form.addParam('fiducialRadius',
                      params.FloatParam,
                      label='Fiducial radius (nm)',
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

        form.addSection('Global variables')

        form.addParam('refineSobelFilter',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Refine center with Sobel filter',
                      expertLevel=params.LEVEL_ADVANCED,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Use edge-detecting Sobel filter to refine the bead positions.')

        form.addParam('scalableSigmaForSobelFilter',
                      params.FloatParam,
                      default=0.5,
                      condition='refineSobelFilter==0',
                      label='Sobel sigma relative to bead size',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Sigma for gaussian kernel filtering of single beads before Sobel filtering, as fraction of '
                           'bead diameter. The default sigma is 0.5 pixels regardless of bead size. '
                           'A value of around 0.12 diameters is needed for higher noise (eg. cryo) data.')

        form.addParam('rotationSolutionType',
                      params.EnumParam,
                      choices=['No rotation', 'One rotation', 'Group rotations', 'Solve for all rotations'],
                      default=3,
                      label='Rotation solution type',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Type of rotation solution.')

        form.addParam('groupRotationSize',
                      params.IntParam,
                      default=5,
                      condition='rotationSolutionType==2',
                      label='Group size',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Size of the rotation group')

        form.addParam('magnificationSolutionType',
                      params.EnumParam,
                      choices=['Fixed magnification at 1.0', 'Group magnifications', 'Solve for all magnifications'],
                      default=1,
                      label='Magnification solution type',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Type of magnification solution.')

        form.addParam('groupMagnificationSize',
                      params.IntParam,
                      default=4,
                      condition='magnificationSolutionType==1',
                      label='Group size',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Size of the magnification group')

        form.addParam('tiltAngleSolutionType',
                      params.EnumParam,
                      choices=['Fixed tilt angles', 'Group tilt angles', 'Solve for all except minimum tilt'],
                      default=1,
                      label='Tilt angle solution type',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Type of tilt angle solution.')

        form.addParam('groupTiltAngleSize',
                      params.IntParam,
                      default=5,
                      condition='tiltAngleSolutionType==1',
                      label='Group size',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Size of the tilt angle group')

        form.addParam('distortionSolutionType',
                      params.EnumParam,
                      choices=['Disabled', 'Full solution', 'Skew only'],
                      default=0,
                      label='Distortion solution type',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Type of distortion solution.')

        form.addParam('xStretchGroupSize',
                      params.IntParam,
                      default=7,
                      condition='distortionSolutionType==1',
                      label='X stretch group size',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Basic grouping size for X stretch')

        form.addParam('skewGroupSize',
                      params.IntParam,
                      default=11,
                      condition='tiltAngleSolutionType==1 or tiltAngleSolutionType==2',
                      label='Skew group size',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Size of the skew group')

        form.addSection('Erase gold beads')

        form.addParam('eraseGoldBeads',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Erase gold beads',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Remove the gold beads detected during fiducial alignment with ccderaser program. This '
                           'option will generate an interpolated tilt series with the gold beads erased and '
                           'interpolated with the calculated transformation matrices form the alignment. ')

        groupEraseGoldBeads = form.addGroup('Gold bead eraser',
                                            condition='eraseGoldBeads==0')

        groupEraseGoldBeads.addParam('betterRadius',
                                     params.IntParam,
                                     default=10,
                                     label='Bead diameter (pixels)',
                                     help="For circle objects, this entry specifies a radius to use for points without "
                                          "an individual point size instead of the object's default sphere radius. "
                                          "This entry is floating point and can be used to overcome the limitations of "
                                          "having an integer default sphere radius. If there are multiple circle "
                                          "objects, enter one value to apply to all objects or a value for each "
                                          "object.")

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):

        self._failedTs = []

        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('generateTrackComStep', ts.getObjId())
            self._insertFunctionStep('generateFiducialSeedStep', ts.getObjId())
            self._insertFunctionStep('generateFiducialModelStep', ts.getObjId())
            self._insertFunctionStep('computeFiducialAlignmentStep', ts.getObjId())
            self._insertFunctionStep('translateFiducialPointModelStep', ts.getObjId())
            self._insertFunctionStep('computeOutputStackStep', ts.getObjId())

            if self.computeAlignment.get() == 0 or self.eraseGoldBeads.get() == 0:
                self._insertFunctionStep('computeOutputInterpolatedStackStep', ts.getObjId())

            if self.eraseGoldBeads.get() == 0:
                self._insertFunctionStep('eraseGoldBeadsStep', ts.getObjId())

            self._insertFunctionStep('computeOutputModelsStep', ts.getObjId())
            self._insertFunctionStep('createOutputFailedSet', ts.getObjId())

        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ----------------------------
    def tryExceptDecorator(func):
        """ This decorator wraps the step in a try/except module which adds the tilt series ID to the failed TS array
        in case the step fails"""
        def wrapper(self, tsId):
            try:
                func(self, tsId)
            except:
                self._failedTs.append(tsId)
        return wrapper

    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        firstItem = ts.getFirstItem()

        """Apply the transformation form the input tilt-series"""
        outputTsFileName = os.path.join(tmpPrefix, firstItem.parseFileName())
        ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt"))
        ts.generateTltFile(angleFilePath)

        """"Link to input tilt-series (needed for fiducial model viewer)"""
        # TODO: there is no need to come from a prealigned stack
        inputTS = os.path.join(extraPrefix, firstItem.parseFileName())
        if firstItem.hasTransform():
            path.copyFile(outputTsFileName, inputTS)

        else:
            path.createLink(firstItem.getLocation()[1], inputTS)

    def generateTrackComStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        firstItem = ts.getFirstItem()

        fiducialRadiusPixel = self.fiducialRadius.get() / (self.inputSetOfTiltSeries.get().getSamplingRate() / 10)

        boxSizeXandY = int(3.3 * self.fiducialRadius.get() / (self.inputSetOfTiltSeries.get().getSamplingRate() / 10))

        # Make boxSizeXandY parameter even due to computational efficiency
        if boxSizeXandY % 2 == 1:
            boxSizeXandY += 1

        paramsDict = {
            'imageFile': os.path.join(tmpPrefix, firstItem.parseFileName()),
            'inputSeedModel': os.path.join(extraPrefix, firstItem.parseFileName(extension=".seed")),
            'outputModel': os.path.join(extraPrefix, firstItem.parseFileName(suffix="_gaps", extension=".fid")),
            'tiltFile': os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt")),
            'rotationAngle': self.rotationAngle.get(),
            'fiducialRadius': fiducialRadiusPixel,
            'samplingRate': self.inputSetOfTiltSeries.get().getSamplingRate() / 10,
            'scalableSigmaForSobelFilter': self.scalableSigmaForSobelFilter.get(),
            'boxSizeXandY': boxSizeXandY,
            'distanceRescueCriterion': 0.75 * fiducialRadiusPixel,
            'postFitRescueResidual': 0.2 * fiducialRadiusPixel,
            'maxRescueDistance': 0.2 * fiducialRadiusPixel,
            'minDiamForParamScaling': 12.5,
            'deletionCriterionMinAndSD': '0.3,2.0'
        }

        self.translateTrackCom(ts, paramsDict)

    @tryExceptDecorator
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

    @tryExceptDecorator
    def generateFiducialModelStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        firstItem = ts.getFirstItem()

        fiducialRadiusPixel = self.fiducialRadius.get() / (self.inputSetOfTiltSeries.get().getSamplingRate() / 10)

        boxSizeXandY = int(3.3 * self.fiducialRadius.get() / (self.inputSetOfTiltSeries.get().getSamplingRate() / 10))

        # Make boxSizeXandY parameter even due to computational efficiency
        if boxSizeXandY % 2 == 1:
            boxSizeXandY += 1

        paramsBeadtrack = {
            'inputSeedModel': os.path.join(extraPrefix, firstItem.parseFileName(extension=".seed")),
            'outputModel': os.path.join(extraPrefix, firstItem.parseFileName(suffix="_gaps", extension=".fid")),
            'imageFile': os.path.join(tmpPrefix, firstItem.parseFileName()),
            'imagesAreBinned': 1,
            'tiltFile': os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt")),
            'tiltDefaultGrouping': 7,
            'magDefaultGrouping': 5,
            'rotDefaultGrouping': 1,
            'minViewsForTiltalign': 4,
            'beadDiameter': fiducialRadiusPixel,
            'fillGaps': 1,
            'maxGapSize': 5,
            'minTiltRangeToFindAxis': 10.0,
            'minTiltRangeToFindAngles': 20.0,
            'boxSizeXandY': "%d,%d" % (boxSizeXandY, boxSizeXandY),
            'roundsOfTracking': 2,
            'localAreaTracking': 1,
            'localAreaTargetSize': 1000,
            'minBeadsInArea': 8,
            'minOverlapBeads': 5,
            'maxBeadsToAverage': 4,
            'sobelFilterCentering': 1,
            'pointsToFitMaxAndMin': '7,3',
            'densityRescueFractionAndSD': '0.6,1.0',
            'distanceRescueCriterion': 0.75 * fiducialRadiusPixel,
            'rescueRelaxationDensityAndDistance': '0.7,0.9',
            'postFitRescueResidual': 0.2 * fiducialRadiusPixel,
            'densityRelaxationPostFit': 0.9,
            'maxRescueDistance': 0.2 * fiducialRadiusPixel,
            'residualsToAnalyzeMaxAndMin': '9,5',
            'deletionCriterionMinAndSD': '0.3,2.0',
            'minDiamForParamScaling': 12.5
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
                        "-BeadDiameter %(beadDiameter).2f " \
                        "-FillGaps %(fillGaps)d " \
                        "-MaxGapSize %(maxGapSize)d " \
                        "-MinTiltRangeToFindAxis %(minTiltRangeToFindAxis).2f " \
                        "-MinTiltRangeToFindAngles %(minTiltRangeToFindAngles).2f " \
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
                        "-DistanceRescueCriterion %(distanceRescueCriterion).2f " \
                        "-RescueRelaxationDensityAndDistance %(rescueRelaxationDensityAndDistance)s " \
                        "-PostFitRescueResidual %(postFitRescueResidual).2f " \
                        "-DensityRelaxationPostFit %(densityRelaxationPostFit).2f " \
                        "-MaxRescueDistance %(maxRescueDistance).2f " \
                        "-ResidualsToAnalyzeMaxAndMin %(residualsToAnalyzeMaxAndMin)s " \
                        "-DeletionCriterionMinAndSD %(deletionCriterionMinAndSD)s " \
                        "-MinDiamForParamScaling %(minDiamForParamScaling).1f"

        Plugin.runImod(self, 'beadtrack', argsBeadtrack % paramsBeadtrack)

    @tryExceptDecorator
    def computeFiducialAlignmentStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        firstItem = ts.getFirstItem()

        paramsTiltAlign = {
            'modelFile': os.path.join(extraPrefix, firstItem.parseFileName(suffix="_gaps", extension=".fid")),
            'imageFile': os.path.join(tmpPrefix, firstItem.parseFileName()),
            'imagesAreBinned': 1,
            'outputModelFile': os.path.join(extraPrefix,
                                            firstItem.parseFileName(suffix="_fidxyz", extension=".mod")),
            'outputResidualFile': os.path.join(extraPrefix,
                                               firstItem.parseFileName(suffix="_resid", extension=".txt")),
            'outputFidXYZFile': os.path.join(extraPrefix,
                                             firstItem.parseFileName(suffix="_fid", extension=".xyz")),
            'outputTiltFile': os.path.join(extraPrefix,
                                           firstItem.parseFileName(suffix="_interpolated", extension=".tlt")),
            'outputTransformFile': os.path.join(extraPrefix,
                                                firstItem.parseFileName(suffix="_fid", extension=".xf")),
            'outputFilledInModel': os.path.join(extraPrefix,
                                                firstItem.parseFileName(suffix="_noGaps", extension=".fid")),
            'rotationAngle': self.rotationAngle.get(),
            'tiltFile': os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt")),
            'angleOffset': 0.0,
            'rotOption': self.getRotationType(),
            'rotDefaultGrouping': self.groupRotationSize.get(),
            'tiltOption': self.getTiltAngleType(),
            'tiltDefaultGrouping': self.groupTiltAngleSize.get(),
            'magReferenceView': 1,
            'magOption': self.getMagnificationType(),
            'magDefaultGrouping': self.groupMagnificationSize.get(),
            'xStretchOption': self.getStretchType(),
            'skewOption': self.getSkewType(),
            'xStretchDefaultGrouping': self.xStretchGroupSize.get(),
            'skewDefaultGrouping': self.skewGroupSize.get(),
            'beamTiltOption': 0,
            'xTiltOption': 0,
            'xTiltDefaultGrouping': 2000,
            'residualReportCriterion': 3.0,
            'surfacesToAnalyze': self.getSurfaceToAnalyze(),
            'metroFactor': 0.25,
            'maximumCycles': 1000,
            'kFactorScaling': 1.0,
            'noSeparateTiltGroups': 1,
            'axisZShift': 0.0,
            'shiftZFromOriginal': 1,
            'localAlignments': 0,
            'outputLocalFile': os.path.join(extraPrefix,
                                            firstItem.parseFileName(suffix="_local", extension=".xf")),
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
            'outputTiltAlignFileText': os.path.join(extraPrefix, "outputTiltAlign.txt"),
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
                        "-SkewDefaultGrouping %(skewDefaultGrouping)d " \
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
                        "-LocalSkewDefaultGrouping %(localSkewDefaultGrouping)d " \
                        "> %(outputTiltAlignFileText)s "

        Plugin.runImod(self, 'tiltalign', argsTiltAlign % paramsTiltAlign)

        self.generateTaSolutionText(os.path.join(extraPrefix, "outputTiltAlign.txt"),
                                    os.path.join(extraPrefix, "taSolution.log"),
                                    ts.getSize(),
                                    ts.getSamplingRate())

    def translateFiducialPointModelStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        firstItem = ts.getFirstItem()

        # Check that previous steps have been completed satisfactorily
        # if os.path.exists(os.path.join(extraPrefix,
        #                                firstItem.parseFileName(suffix="_gaps", extension=".fid"))):
        #     paramsGapModel2Point = {
        #         'inputFile': os.path.join(extraPrefix,
        #                                   firstItem.parseFileName(suffix="_gaps", extension=".fid")),
        #         'outputFile': os.path.join(extraPrefix,
        #                                    firstItem.parseFileName(suffix="_gaps_fid", extension=".txt"))
        #     }
        #     argsGapModel2Point = "-InputFile %(inputFile)s " \
        #                          "-OutputFile %(outputFile)s"
        #
        #     Plugin.runImod(self, 'model2point', argsGapModel2Point % paramsGapModel2Point)

        # Check that previous steps have been completed satisfactorily
        if os.path.exists(os.path.join(extraPrefix,
                                       firstItem.parseFileName(suffix="_noGaps", extension=".fid"))):
            paramsNoGapModel2Point = {
                'inputFile': os.path.join(extraPrefix,
                                          firstItem.parseFileName(suffix="_noGaps", extension=".fid")),
                'outputFile': os.path.join(extraPrefix,
                                           firstItem.parseFileName(suffix="_noGaps_fid", extension=".txt"))
            }
            argsNoGapModel2Point = "-InputFile %(inputFile)s " \
                                   "-OutputFile %(outputFile)s"

            Plugin.runImod(self, 'model2point', argsNoGapModel2Point % paramsNoGapModel2Point)

    def computeOutputStackStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        firstItem = ts.getFirstItem()

        # Check that previous steps have been completed satisfactorily
        if os.path.exists(os.path.join(extraPrefix, firstItem.parseFileName(suffix="_fid", extension=".xf"))):
            tltFilePath = os.path.join(
                extraPrefix,
                firstItem.parseFileName(suffix="_interpolated", extension=".tlt")
            )

            tltList = utils.formatAngleList(tltFilePath)

            transformationMatricesFilePath = os.path.join(
                extraPrefix,
                firstItem.parseFileName(suffix="_fid", extension=".xf")
            )

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
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        firstItem = ts.getFirstItem()

        # Check that previous steps have been completed satisfactorily
        if os.path.exists(os.path.join(extraPrefix, firstItem.parseFileName(suffix="_fid", extension=".xf"))):
            outputInterpolatedSetOfTiltSeries = self.getOutputInterpolatedSetOfTiltSeries()

            paramsAlignment = {
                'input': os.path.join(tmpPrefix, firstItem.parseFileName()),
                'output': os.path.join(extraPrefix, firstItem.parseFileName(suffix="_interpolated")),
                'xform': os.path.join(extraPrefix, firstItem.parseFileName(suffix="_fid", extension=".xf")),
                'bin': int(self.binning.get()),
                'imagebinned': 1.0}

            argsAlignment = "-input %(input)s " \
                            "-output %(output)s " \
                            "-xform %(xform)s " \
                            "-bin %(bin)d " \
                            "-imagebinned %(imagebinned)s "

            rotationAngleAvg = utils.calculateRotationAngleFromTM(self.getOutputSetOfTiltSeries()[tsObjId])

            # Check if rotation angle is greater than 45º. If so, swap x and y dimensions to adapt output image sizes to
            # the final sample disposition.
            if rotationAngleAvg > 45 or rotationAngleAvg < -45:
                paramsAlignment.update({
                    'size': "%d,%d" % (firstItem.getYDim(), firstItem.getXDim())
                })

                argsAlignment += " -size %(size)s "

            Plugin.runImod(self, 'newstack', argsAlignment % paramsAlignment)

            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            outputInterpolatedSetOfTiltSeries.append(newTs)

            tltFilePath = os.path.join(
                extraPrefix,
                firstItem.parseFileName(suffix="_interpolated", extension=".tlt")
            )

            tltList = utils.formatAngleList(tltFilePath)

            if self.binning > 1:
                newTs.setSamplingRate(ts.getSamplingRate() * int(self.binning.get()))

            for index, ti in enumerate(ts):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(ti, copyId=True)
                newTi.setLocation(index + 1, os.path.join(extraPrefix, ti.parseFileName(suffix="_interpolated")))
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

    @tryExceptDecorator
    def eraseGoldBeadsStep(self, tsObjId):
            ts = self.inputSetOfTiltSeries.get()[tsObjId]
            tsId = ts.getTsId()

            extraPrefix = self._getExtraPath(tsId)
            tmpPrefix = self._getTmpPath(tsId)

            firstItem = ts.getFirstItem()

            # Move interpolated tilt-series to tmp folder and generate a new one with the gold beads erased back in the
            # extra folder
            path.moveFile(os.path.join(extraPrefix, ts.getFirstItem().parseFileName(suffix='_interpolated')),
                          os.path.join(tmpPrefix, ts.getFirstItem().parseFileName(suffix='_interpolated')))

            # Generate interpolated model
            paramsImodtrans = {
                'inputFile': os.path.join(extraPrefix,
                                          firstItem.parseFileName(suffix="_noGaps", extension=".fid")),
                'outputFile': os.path.join(extraPrefix,
                                           firstItem.parseFileName(suffix="_noGaps_ali", extension=".fid")),
                'transformFile': os.path.join(extraPrefix,
                                              firstItem.parseFileName(suffix="_fid", extension=".xf"))
            }

            argsImodtrans = "-2 %(transformFile)s " \
                            "%(inputFile)s " \
                            "%(outputFile)s "

            Plugin.runImod(self, 'imodtrans', argsImodtrans % paramsImodtrans)

            # Erase gold beads
            paramsCcderaser = {
                'inputFile': os.path.join(tmpPrefix, firstItem.parseFileName(suffix='_interpolated')),
                'outputFile': os.path.join(extraPrefix, firstItem.parseFileName(suffix='_interpolated')),
                'modelFile': os.path.join(extraPrefix,
                                          firstItem.parseFileName(suffix="_noGaps_ali", extension=".fid")),
                'betterRadius': self.betterRadius.get(),
                'polynomialOrder': 0,
                'circleObjects': "/"
            }

            argsCcderaser = "-InputFile %(inputFile)s " \
                            "-OutputFile %(outputFile)s " \
                            "-ModelFile %(modelFile)s " \
                            "-BetterRadius %(betterRadius)d " \
                            "-PolynomialOrder %(polynomialOrder)d " \
                            "-CircleObjects %(circleObjects)s " \
                            "-MergePatches " \
                            "-ExcludeAdjacent"

            Plugin.runImod(self, 'ccderaser', argsCcderaser % paramsCcderaser)

    def computeOutputModelsStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        firstItem = ts.getFirstItem()

        # Create the output set of landmark models with gaps
        # Check that previous steps have been completed satisfactorily
        # if os.path.exists(
        #         os.path.join(extraPrefix, firstItem.parseFileName(suffix="_gaps", extension=".fid"))):
        #
        #     outputSetOfLandmarkModelsGaps = self.getOutputFiducialModelGaps()
        #
        #     landmarkModelNoGapsFilePath = os.path.join(
        #         extraPrefix,
        #         firstItem.parseFileName(suffix="_gaps", extension=".sfid")
        #     )
        #
        #     fiducialModelGapPath = os.path.join(
        #         extraPrefix,
        #         firstItem.parseFileName(suffix="_gaps", extension=".fid")
        #     )
        #
        #     landmarkModelGapsResidPath = os.path.join(
        #         extraPrefix,
        #         firstItem.parseFileName(suffix="_resid", extension=".txt")
        #     )
        #
        #     fiducialGapResidList = utils.formatFiducialResidList(landmarkModelGapsResidPath)
        #
        #     landmarkModelGaps = LandmarkModel(tsId=tsId,
        #                                       fileName=landmarkModelNoGapsFilePath,
        #                                       modelName=fiducialModelGapPath)
        #
        #     prevTiltIm = 0
        #     chainId = 0
        #     for index, fiducial in enumerate(fiducialGapResidList):
        #         if int(fiducial[2]) <= prevTiltIm:
        #             chainId += 1
        #         prevTiltIm = int(fiducial[2])
        #         landmarkModelGaps.addLandmark(xCoor=fiducial[0],
        #                                       yCoor=fiducial[1],
        #                                       tiltIm=fiducial[2],
        #                                       chainId=chainId,
        #                                       xResid=fiducial[3],
        #                                       yResid=fiducial[4])
        #
        #     outputSetOfLandmarkModelsGaps.append(landmarkModelGaps)
        #     outputSetOfLandmarkModelsGaps.update(landmarkModelGaps)
        #     outputSetOfLandmarkModelsGaps.write()

        # Create the output set of landmark models with no gaps
        if os.path.exists(
                os.path.join(extraPrefix, ts.getFirstItem().parseFileName(suffix="_noGaps_fid", extension=".txt"))):

            outputSetOfLandmarkModelsNoGaps = self.getOutputFiducialModelNoGaps()

            fiducialNoGapFilePath = os.path.join(
                extraPrefix,
                firstItem.parseFileName(suffix="_noGaps_fid", extension=".txt")
            )

            fiducialNoGapList = utils.formatFiducialList(fiducialNoGapFilePath)

            fiducialModelNoGapPath = os.path.join(
                extraPrefix,
                firstItem.parseFileName(suffix="_noGaps", extension=".fid")
            )

            landmarkModelNoGapsFilePath = os.path.join(
                extraPrefix,
                firstItem.parseFileName(suffix="_noGaps", extension=".sfid")
            )

            landmarkModelNoGapsResidPath = os.path.join(
                extraPrefix,
                firstItem.parseFileName(suffix="_resid", extension=".txt")
            )

            fiducialNoGapsResidList = utils.formatFiducialResidList(landmarkModelNoGapsResidPath)

            landmarkModelNoGaps = LandmarkModel(tsId=tsId,
                                                fileName=landmarkModelNoGapsFilePath,
                                                modelName=fiducialModelNoGapPath)

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

        # Create the output set of landmark models with no gaps
        if os.path.exists(
                os.path.join(extraPrefix, ts.getFirstItem().parseFileName(suffix="_fid", extension=".xyz"))):

            outputSetOfCoordinates3D = self.getOutputSetOfCoordinates3Ds()

            coordFilePath = os.path.join(
                extraPrefix,
                firstItem.parseFileName(suffix="_fid", extension=".xyz")
            )

            xDim = firstItem.getXDim()
            yDim = firstItem.getYDim()
            coordList = utils.format3DCoordinatesList(coordFilePath)

            for element in coordList:
                newCoord3D = tomoObj.Coordinate3D()
                newCoord3D.setVolume(ts)
                newCoord3D.setX(element[0], constants.BOTTOM_LEFT_CORNER)
                newCoord3D.setY(element[1], constants.BOTTOM_LEFT_CORNER)
                newCoord3D.setZ(element[2], constants.BOTTOM_LEFT_CORNER)

                newCoord3D.setVolId(tsObjId)
                outputSetOfCoordinates3D.append(newCoord3D)
                outputSetOfCoordinates3D.update(newCoord3D)
            outputSetOfCoordinates3D.write()
            self._store()

    def createOutputFailedSet(self, tsObjId):
        # Check if the tilt-series ID is in the failed tilt-series list to add it to the set
        if tsObjId in self._failedTs:
            outputFailedSetOfTiltSeries = self.getOutputFailedSetOfTiltSeries()

            ts = self.inputSetOfTiltSeries.get()[tsObjId]
            tsId = ts.getTsId()

            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            outputFailedSetOfTiltSeries.append(newTs)

            for index, ti in enumerate(ts):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(ti, copyId=True)
                newTi.setLocation(ti.getLocation())
                if self.binning > 1:
                    newTi.setSamplingRate(ti.getSamplingRate() * int(self.binning.get()))
                newTs.append(newTi)

            ih = ImageHandler()
            x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
            newTs.setDim((x, y, z))
            newTs.write(properties=False)

            outputFailedSetOfTiltSeries.update(newTs)
            outputFailedSetOfTiltSeries.updateDim()
            outputFailedSetOfTiltSeries.write()
            self._store()

    def createOutputStep(self):
        if hasattr(self, "outputSetOfTiltSeries"):
            self.getOutputSetOfTiltSeries().setStreamState(Set.STREAM_CLOSED)
        if hasattr(self, "outputInterpolatedSetOfTiltSeries"):
            self.getOutputInterpolatedSetOfTiltSeries().setStreamState(Set.STREAM_CLOSED)
        # if hasattr(self, "outputFiducialModelGaps"):
        #     self.getOutputFiducialModelGaps().setStreamState(Set.STREAM_CLOSED)
        if hasattr(self, "outputFiducialModelNoGaps"):
            self.getOutputFiducialModelNoGaps().setStreamState(Set.STREAM_CLOSED)
        if hasattr(self, "outputSetOfCoordinates3D"):
            self.getOutputSetOfCoordinates3Ds().setStreamState(Set.STREAM_CLOSED)
        if hasattr(self, "outputFailedSetOfTiltSeries"):
            self.getOutputFailedSetOfTiltSeries().setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getRotationType(self):
        if self.rotationSolutionType.get() == 0:
            return 0
        elif self.rotationSolutionType.get() == 1:
            return -1
        elif self.rotationSolutionType.get() == 2:
            return 3
        elif self.rotationSolutionType.get() == 3:
            return 1

    def getMagnificationType(self):
        if self.magnificationSolutionType.get() == 0:
            return 0
        elif self.magnificationSolutionType.get() == 1:
            return 3
        elif self.magnificationSolutionType.get() == 2:
            return 1

    def getTiltAngleType(self):
        if self.tiltAngleSolutionType.get() == 0:
            return 0
        elif self.tiltAngleSolutionType.get() == 1:
            return 5
        elif self.tiltAngleSolutionType.get() == 2:
            return 2

    def getSkewType(self):
        if self.distortionSolutionType.get() == 0:
            return 0
        elif self.distortionSolutionType.get() == 1 or self.distortionSolutionType.get() == 2:
            return 3

    def getStretchType(self):
        if self.distortionSolutionType.get() == 0 or self.distortionSolutionType.get() == 2:
            return 0
        elif self.distortionSolutionType.get() == 1:
            return 3

    def getSurfaceToAnalyze(self):
        if self.twoSurfaces.get() == 0:
            return 2
        elif self.twoSurfaces.get() == 1:
            return 1

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
RotationAngle	%(rotationAngle).1f 
TiltFile	%(tiltFile)s
TiltDefaultGrouping	7
MagDefaultGrouping	5
RotDefaultGrouping	1
PixelSize   %(samplingRate)f
BeadDiameter	%(fiducialRadius).2f
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
BoxSizeXandY	%(boxSizeXandY)d,%(boxSizeXandY)d
MaxBeadsToAverage	4
# points and minimum for extrapolation
PointsToFitMaxAndMin	7,3
# fraction of mean, and # of SD below mean: density criterion for rescue
DensityRescueFractionAndSD	0.6,1.0
# distance criterion for rescue
DistanceRescueCriterion	%(distanceRescueCriterion).2f
# relaxation of criterion for density and distance rescues
RescueRelaxationDensityAndDistance	0.7,0.9
# distance for rescue after fit
PostFitRescueResidual	%(postFitRescueResidual).2f
# relaxation of density criterion, maximum radius to search
DensityRelaxationPostFit	0.9
MaxRescueDistance	%(maxRescueDistance).2f
# Max and min residual changes to use to get mean and SD change
ResidualsToAnalyzeMaxAndMin	9,5
# minimum residual difference, criterion # of sd's
DeletionCriterionMinAndSD	%(deletionCriterionMinAndSD)s
MinDiamForParamScaling %(minDiamForParamScaling).1f
"""

        if self.refineSobelFilter.get() == 0:
            template += """SobelFilterCentering
ScalableSigmaForSobel   %(scalableSigmaForSobelFilter)f
$if (-e ./savework) ./savework
"""
        elif self.refineSobelFilter.get() == 1:
            template += """$if (-e ./savework) ./savework"""

        with open(trackFilePath, 'w') as f:
            f.write(template % paramsDict)

    def generateTaSolutionText(self, tiltAlignOutputLog, taSolutionLog, numberOfTiltImages, pixelSize):
        """ This method generates a text file containing the TA solution from the tiltalign output log. """

        searchingPassword = "deltilt"

        with open(tiltAlignOutputLog, 'r') as fRead:
            lines = fRead.readlines()

            counts = []

            for index, line in enumerate(lines):
                if searchingPassword in line:
                    counts.append([index])

        lastApparition = max(counts)[0]

        outputLinesAsMatrix = []

        # Take only the lines that compose the table containing the ta solution info (until blank line)
        # Convert lines into numpy array for posterior operation

        index = lastApparition + 1
        while True:
            vector = lines[index].split()
            vector = [float(i) for i in vector]
            outputLinesAsMatrix.append(vector)
            if int(vector[0]) == numberOfTiltImages:
                break
            index += 1

        matrixTaSolution = np.array(outputLinesAsMatrix)

        # Find the position in table of the minimum tilt angle image
        _, indexAng = min((abs(val), idx) for (idx, val) in enumerate(matrixTaSolution[:, 2]))

        # Multiply last column by the sampling rate in nanometer
        matrixTaSolution[:, -1] = matrixTaSolution[:, -1] * pixelSize / 10

        # Get minimum rotation to write in file
        minimumRotation = matrixTaSolution[indexAng][1]

        # Save new matrixTaSolution info into file
        np.savetxt(fname=taSolutionLog,
                   X=matrixTaSolution,
                   fmt=" %i\t%.1f\t%.1f\t%.2f\t%.4f\t%.4f\t%.2f\t%.2f",
                   header=" At minimum tilt, rotation angle is %.2f\n\n"
                          " view   rotation    tilt    deltilt     mag      dmag      skew    resid-nm"
                          % minimumRotation,
                   comments='')

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

    def getOutputFailedSetOfTiltSeries(self):
        if hasattr(self, "outputFailedSetOfTiltSeries"):
            self.outputFailedSetOfTiltSeries.enableAppend()
        else:
            outputFailedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Failed')
            outputFailedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputFailedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            outputFailedSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputFailedSetOfTiltSeries=outputFailedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputFailedSetOfTiltSeries)
        return self.outputFailedSetOfTiltSeries

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

    # def getOutputFiducialModelGaps(self):
    #     if hasattr(self, "outputFiducialModelGaps"):
    #         self.outputFiducialModelGaps.enableAppend()
    #     else:
    #         outputFiducialModelGaps = self._createSetOfLandmarkModels(suffix='Gaps')
    #         outputFiducialModelGaps.copyInfo(self.inputSetOfTiltSeries.get())
    #         outputFiducialModelGaps.setStreamState(Set.STREAM_OPEN)
    #         self._defineOutputs(outputFiducialModelGaps=outputFiducialModelGaps)
    #         self._defineSourceRelation(self.inputSetOfTiltSeries, outputFiducialModelGaps)
    #     return self.outputFiducialModelGaps

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
                                                                      suffix='Fiducials3D')
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

        if hasattr(self, 'outputFailedSetOfTiltSeries'):
            summary.append("Failed tilt-series: %d."
                           % (self.outputFailedSetOfTiltSeries.getSize()))

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

        if hasattr(self, 'outputFailedSetOfTiltSeries'):
            methods.append("%d tilt-series have failed during the fiducial alignment protocol execution."
                           % (self.outputFailedSetOfTiltSeries.getSize()))

        if not methods:
            methods.append("Output classes not ready yet.")
        return methods
