# *****************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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

from pyworkflow import BETA
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
import tomo.objects as tomoObj

from .. import Plugin, utils
from .protocol_base import ProtImodBase, TLT_EXT, XF_EXT


class ProtImodFiducialModel(ProtImodBase):
    """
    Construction of a fiducial model and alignment of tilt-series based
    on the IMOD procedure.

    More info:
        https://bio3d.colorado.edu/imod/doc/man/autofidseed.html
        https://bio3d.colorado.edu/imod/doc/man/beadtrack.html
        https://bio3d.colorado.edu/imod/doc/man/model2point.html
    """

    _label = 'Generate fiducial model'
    _devStatus = BETA

    FIDUCIAL_MODEL = 0
    PATCH_TRACKING = 1

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('typeOfModel',
                      params.EnumParam,
                      choices=["Make seed and Track", "Patch Tracking"],
                      default=0,
                      important=True,
                      label='Model generation')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt Series')

        self._patchTrackingForm(form, 'typeOfModel==%d' % self.PATCH_TRACKING)
        self._fiducialSeedForm(form, 'typeOfModel==%d' % self.FIDUCIAL_MODEL)

    def _patchTrackingForm(self, form, condition, levelType=params.LEVEL_NORMAL):
        patchtrack = form.addGroup('Patch Tracking', expertLevel=levelType, condition=condition)

        patchtrack.addParam('sizeOfPatches', params.NumericListParam,
                            label='Size of the patches (X,Y)',
                            default='100 100',
                            expertLevel=levelType,
                            help="Size of the  patches to track by correlation. In imod documentation "
                                 "(tiltxcorr: SizeOfPatchesXandY)")

        patchtrack.addParam('patchLayout',
                            params.EnumParam,
                            choices=['Fractional overlap of patches',
                                     'Number of patches'],
                            default=0,
                            label='Patch layout',
                            display=params.EnumParam.DISPLAY_HLIST,
                            help='To be added')

        patchtrack.addParam('overlapPatches',
                            params.NumericListParam,
                            default='0.33 0.33',
                            condition='patchLayout==0',
                            label='Fractional overlap of the patches (X,Y)',
                            help="Fractional overlap in X and Y to track by correlation. In imod documentation"
                                 "(tiltxcorr: NumberOfPatchesXandY)")

        patchtrack.addParam('numberOfPatches',
                            params.NumericListParam,
                            condition='patchLayout==1',
                            label='Number of patches (X,Y)',
                            help="Number of patches in X and Y of the patches. In imod documentation"
                                 "(tiltxcorr: OverlapOfPatchesXandY)")

        patchtrack.addParam('iterationsSubpixel',
                            params.IntParam,
                            default=1,
                            label='Iterations to increase subpixel accuracy',
                            help="Number of iteration of each correlation to reduce interpolation of the peak position"
                                 "In imod documentation: (tiltxcorr: IterateCorrelations)")

        self.trimimgForm(patchtrack, pxTrimCondition='True', correlationCondition='True',
                         levelType=params.LEVEL_ADVANCED)
        self.filteringParametersForm(form, condition=condition, levelType=params.LEVEL_ADVANCED)

    def _fiducialSeedForm(self, form, condition, levelType=params.LEVEL_NORMAL):
        seedModel = form.addGroup('"Make seed and Track', expertLevel=levelType, condition=condition)
        seedModel.addParam('fiducialDiameter',
                           params.FloatParam,
                           label='Fiducial diameter (nm)',
                           default='10',
                           important=True,
                           help="Fiducials diameter to be tracked for alignment.")

        seedModel.addParam('twoSurfaces',
                           params.EnumParam,
                           choices=['Yes', 'No'],
                           default=1,
                           label='Find beads on two surfaces?',
                           display=params.EnumParam.DISPLAY_HLIST,
                           help="Track fiducials differentiating in which side "
                                "of the sample are located.")

        seedModel.addParam('numberFiducial',
                           params.IntParam,
                           label='Number of fiducials',
                           default='25',
                           expertLevel=params.LEVEL_ADVANCED,
                           help="Number of fiducials to be tracked for alignment.")

        seedModel.addParam('doTrackWithModel', params.BooleanParam,
                           default=True,
                           label="Track with fiducial model as seed",
                           help="Turn the tracked model into new seed and "
                                "repeat tracking.")

        seedModel.addParam('shiftsNearZeroFraction',
                           params.FloatParam,
                           label='Shifts near zero fraction',
                           default='0.2',
                           expertLevel=params.LEVEL_ADVANCED,
                           help="Fraction of the tracking box size above which to "
                                "supply shifts near zero tilt to Beadtrack. The "
                                "dominant net shifts in the bead positions between "
                                "views are found as described above, and if one of "
                                "the shifts is larger than this fraction of the "
                                "-BoxSizeXandY entry to Beadtrack, then the shifts "
                                "are provided when running Beadtrack on the initial "
                                "seed models. Also, a command file will be written "
                                "with modified parameters, named as the root name "
                                "of the input command file followed by '_adjusted' "
                                "and its extension. Enter 0 or a large value to "
                                "disable this analysis.")

        groupGlobalVariables = form.addGroup('Filter variables',
                                             expertLevel=params.LEVEL_ADVANCED, condition=condition)

        groupGlobalVariables.addParam('refineSobelFilter',
                                      params.EnumParam,
                                      choices=['Yes', 'No'],
                                      default=0,
                                      label='Refine center with Sobel filter?',
                                      expertLevel=params.LEVEL_ADVANCED,
                                      display=params.EnumParam.DISPLAY_HLIST,
                                      help='Use edge-detecting Sobel filter '
                                           'to refine the bead positions.')

        groupGlobalVariables.addParam('scalableSigmaForSobelFilter',
                                      params.FloatParam,
                                      default=0.12,
                                      condition='refineSobelFilter==0',
                                      label='Sobel sigma relative to bead size',
                                      expertLevel=params.LEVEL_ADVANCED,
                                      help='Sigma for gaussian kernel filtering '
                                           'of single beads before Sobel '
                                           'filtering, as fraction of bead '
                                           'diameter. The default sigma is 0.5 '
                                           'pixels regardless of bead size. '
                                           'A value of around 0.12 diameters is '
                                           'needed for higher noise (eg. cryo) '
                                           'data.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        for tsId in self.tsDict.keys():
            self._insertFunctionStep(self.convertInputStep, tsId)
            if self.typeOfModel == self.FIDUCIAL_MODEL:
                self._insertFunctionStep(self.generateTrackComStep, tsId)
                self._insertFunctionStep(self.generateFiducialSeedStep, tsId)
                self._insertFunctionStep(self.generateFiducialModelStep, tsId)
            else:
                self._insertFunctionStep(self.xcorrStep, tsId)
                self._insertFunctionStep(self.chopcontsStep, tsId)
            self._insertFunctionStep(self.translateFiducialPointModelStep, tsId)
            self._insertFunctionStep(self.computeOutputModelsStep, tsId)
            self._insertFunctionStep(self.createOutputFailedSetStep, tsId)

        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self._failedTs = []
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self.inputSetOfTiltSeries.get()}

    def generateTrackComStep(self, tsId):
        ts = self.tsDict[tsId]
        fiducialDiameterPixel = self.fiducialDiameter.get() / (self.inputSetOfTiltSeries.get().getSamplingRate() / 10)
        boxSizeXandY = max(3.3 * fiducialDiameterPixel + 2, 2 * fiducialDiameterPixel + 20, 32)
        boxSizeXandY = min(512, 2 * int(boxSizeXandY / 2))
        scaling = fiducialDiameterPixel / 12.5 if fiducialDiameterPixel > 12.5 else 1

        paramsDict = {
            'imageFile': self.getTmpOutFile(tsId),
            'inputSeedModel': self.getExtraOutFile(tsId, ext="seed"),
            'outputModel': self.getExtraOutFile(tsId, suffix="gaps", ext="fid"),
            'tiltFile': self.getExtraOutFile(tsId, ext=TLT_EXT),
            'rotationAngle': ts.getAcquisition().getTiltAxisAngle(),
            'fiducialDiameter': fiducialDiameterPixel,
            'samplingRate': self.inputSetOfTiltSeries.get().getSamplingRate() / 10,
            'scalableSigmaForSobelFilter': self.scalableSigmaForSobelFilter.get(),
            'boxSizeXandY': boxSizeXandY,
            'distanceRescueCriterion': 10 * scaling,
            'postFitRescueResidual': 2.5 * scaling,
            'maxRescueDistance': 2.5 * scaling,
            'minDiamForParamScaling': 12.5,
            'deletionCriterionMinAndSD': f"{0.04 * scaling:0.3f},2.0",
        }

        self.translateTrackCom(ts, paramsDict)

    @ProtImodBase.tryExceptDecorator
    def generateFiducialSeedStep(self, tsId):
        paramsAutofidseed = {
            'trackCommandFile': self.getExtraOutFile(tsId, suffix="track", ext="com"),
            'minSpacing': 0.85,
            'peakStorageFraction': 1.0,
            'targetNumberOfBeads': self.numberFiducial.get(),
            'shiftsNearZeroFraction': self.shiftsNearZeroFraction.get()
        }

        argsAutofidseed = "-TrackCommandFile %(trackCommandFile)s " \
                          "-MinSpacing %(minSpacing)f " \
                          "-PeakStorageFraction %(peakStorageFraction)f " \
                          "-TargetNumberOfBeads %(targetNumberOfBeads)d " \
                          "-ShiftsNearZeroFraction %(shiftsNearZeroFraction)f " \
                          "-AdjustSizes "

        if self.twoSurfaces.get() == 0:
            argsAutofidseed += " -TwoSurfaces "

        Plugin.runImod(self, 'autofidseed', (argsAutofidseed % paramsAutofidseed))

        autofidseedDirPath = os.path.join(self._getExtraPath(tsId), "autofidseed.dir")
        path.makePath(autofidseedDirPath)
        path.moveTree("autofidseed.dir", autofidseedDirPath)
        path.moveFile("autofidseed.info", self._getExtraPath(tsId))

    @ProtImodBase.tryExceptDecorator
    def generateFiducialModelStep(self, tsId):
        ts = self.tsDict[tsId]
        firstItem = ts.getFirstItem()

        fiducialDiameterPixel = self.fiducialDiameter.get() / (self.inputSetOfTiltSeries.get().getSamplingRate() / 10)
        boxSizeXandY = max(3.3 * fiducialDiameterPixel + 2, 2 * fiducialDiameterPixel + 20, 32)
        boxSizeXandY = min(512, 2 * int(boxSizeXandY / 2))
        scaling = fiducialDiameterPixel / 12.5 if fiducialDiameterPixel > 12.5 else 1

        paramsBeadtrack = {
            'inputSeedModel': self.getExtraOutFile(tsId, ext="seed"),
            'outputModel': self.getExtraOutFile(tsId, suffix="gaps", ext="fid"),
            'imageFile': self.getTmpOutFile(tsId),
            'imagesAreBinned': 1,
            'tiltFile': self.getExtraOutFile(tsId, ext=TLT_EXT),
            'tiltDefaultGrouping': 7,
            'magDefaultGrouping': 5,
            'rotDefaultGrouping': 1,
            'minViewsForTiltalign': 4,
            'beadDiameter': fiducialDiameterPixel,
            'fillGaps': 1,
            'maxGapSize': 5,
            'minTiltRangeToFindAxis': 10.0,
            'minTiltRangeToFindAngles': 20.0,
            'boxSizeXandY': f"{boxSizeXandY},{boxSizeXandY}",
            'roundsOfTracking': 2,
            'localAreaTracking': 1,
            'localAreaTargetSize': 1000,
            'minBeadsInArea': 8,
            'minOverlapBeads': 5,
            'maxBeadsToAverage': 4,
            'sobelFilterCentering': 1,
            'pointsToFitMaxAndMin': '7,3',
            'densityRescueFractionAndSD': '0.6,1.0',
            'distanceRescueCriterion': 10 * scaling,
            'rescueRelaxationDensityAndDistance': '0.7,0.9',
            'postFitRescueResidual': 2.5 * scaling,
            'densityRelaxationPostFit': 0.9,
            'maxRescueDistance': 2.5 * scaling,
            'residualsToAnalyzeMaxAndMin': '9,5',
            'deletionCriterionMinAndSD': f"{0.04 * scaling:0.3f},2.0",
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
                        "-MinDiamForParamScaling %(minDiamForParamScaling).1f "

        # Excluded views
        excludedViews = ts.getExcludedViewsIndex(caster=str)
        if len(excludedViews):
            argsBeadtrack += f"-SkipViews {','.join(excludedViews)} "

        if firstItem.hasTransform():
            XfFileName = self.getExtraOutFile(tsId, ext=XF_EXT)
            argsBeadtrack += f"-prexf {XfFileName} "

        Plugin.runImod(self, 'beadtrack', argsBeadtrack % paramsBeadtrack)

        if self.doTrackWithModel:
            # repeat tracking with the current model as seed
            path.copyFile(paramsBeadtrack['inputSeedModel'],
                          self.getExtraOutFile(tsId, suffix="orig", ext="seed"))
            path.moveFile(paramsBeadtrack['outputModel'],
                          paramsBeadtrack['inputSeedModel'])

            Plugin.runImod(self, 'beadtrack', argsBeadtrack % paramsBeadtrack)

    def translateFiducialPointModelStep(self, tsId):
        ts = self.tsDict[tsId]

        # Check that previous steps have been completed satisfactorily
        gapsFidFile = self.getExtraOutFile(tsId, suffix='gaps', ext='fid')
        if os.path.exists(gapsFidFile):
            paramsGapModel2Point = {
                'inputFile': gapsFidFile,
                'outputFile': self.getExtraOutFile(tsId, suffix="gaps_fid", ext="txt")
            }
            argsGapModel2Point = "-InputFile %(inputFile)s " \
                                 "-OutputFile %(outputFile)s"

            Plugin.runImod(self, 'model2point', argsGapModel2Point % paramsGapModel2Point)

    def computeOutputModelsStep(self, tsId):
        ts = self.tsDict[tsId]

        # Create the output set of landmark models with gaps
        # Check that previous steps have been completed satisfactorily
        fiducialModelGapPath = self.getExtraOutFile(tsId, suffix='gaps', ext='fid')
        if os.path.exists(fiducialModelGapPath):
            output = self.getOutputFiducialModelGaps()
            landmarkModelGapsFilePath = self.getExtraOutFile(tsId, suffix='gaps', ext='sfid')
            fiducialModelGapTxtPath = self.getExtraOutFile(tsId, suffix="gaps_fid", ext="txt")

            fiducialGapList = utils.formatFiducialList(fiducialModelGapTxtPath)
            fiducialDiameterPixel = self.fiducialDiameter.get() / (
                    self.inputSetOfTiltSeries.get().getSamplingRate() / 10)

            landmarkModelGaps = tomoObj.LandmarkModel(tsId=tsId,
                                                      tiltSeriesPointer=ts,
                                                      fileName=landmarkModelGapsFilePath,
                                                      modelName=fiducialModelGapPath,
                                                      size=fiducialDiameterPixel,
                                                      hasResidualInfo=False)

            landmarkModelGaps.setTiltSeries(ts)

            prevTiltIm = 0
            chainId = 0

            for index, fiducial in enumerate(fiducialGapList):
                if int(fiducial[2]) <= prevTiltIm:
                    chainId += 1

                prevTiltIm = int(fiducial[2])

                landmarkModelGaps.addLandmark(xCoor=fiducial[0],
                                              yCoor=fiducial[1],
                                              tiltIm=fiducial[2] + 1,
                                              chainId=chainId,
                                              xResid=0,
                                              yResid=0)

            output.append(landmarkModelGaps)
            output.update(landmarkModelGaps)
            output.write()

    def xcorrStep(self, tsId):
        """
        Imod uses the next command line for the xcorr alignment
        $tiltxcorr                  -StandardInput
        InputFile	                cryo_preali.mrc
        OutputFile	                cryo_pt.fid
        RotationAngle	            -12.6
        TiltFile	                cryo.rawtlt
        FilterRadius2	            0.125
        FilterSigma1	            0.03
        FilterSigma2	            0.03
        BordersInXandY	            102,102
        IterateCorrelations	        1
        SizeOfPatchesXandY	        680,680
        OverlapOfPatchesXandY	    0.33,0.33
        PrealignmentTransformFile	cryo.prexg
        ImagesAreBinned	            1
        """
        ts = self.tsDict[tsId]
        angleFilePath = self.getExtraOutFile(tsId, ext=TLT_EXT)
        xfFile = self.getExtraOutFile(tsId, ext=XF_EXT)
        ts.writeXfFile(xfFile)

        borders = self.pxTrim.getListFromValues()
        sizePatches = self.sizeOfPatches.getListFromValues()

        BordersInXandY = '%d,%d' % (borders[0], borders[1])
        SizeOfPatchesXandY = '%d,%d' % (sizePatches[0], sizePatches[1])

        paramsTiltXCorr = {
            'inputFile': self.getTmpOutFile(tsId),
            'outputFile': self.getExtraOutFile(tsId, suffix="pt", ext=".fid"),
            'RotationAngle': ts.getAcquisition().getTiltAxisAngle(),
            'TiltFile': angleFilePath,
            'FilterRadius2': self.filterRadius2.get(),
            'FilterSigma1': self.filterSigma1.get(),
            'FilterSigma2': self.filterSigma2.get(),
            'BordersInXandY': BordersInXandY,
            'IterateCorrelations': self.iterationsSubpixel.get(),
            'SizeOfPatchesXandY': SizeOfPatchesXandY,
            'PrealignmentTransformFile': xfFile,

            'ImagesAreBinned': 1,
        }
        argsTiltXCorr = " " \
                        "-InputFile %(inputFile)s " \
                        "-OutputFile %(outputFile)s " \
                        "-RotationAngle %(RotationAngle)s " \
                        "-TiltFile %(TiltFile)s " \
                        "-FilterRadius2 %(FilterRadius2)s " \
                        "-FilterSigma1 %(FilterSigma1)s " \
                        "-FilterSigma2 %(FilterSigma2)s " \
                        "-BordersInXandY %(BordersInXandY)s " \
                        "-IterateCorrelations %(IterateCorrelations)s " \
                        "-SizeOfPatchesXandY %(SizeOfPatchesXandY)s " \
                        "-PrealignmentTransformFile %(PrealignmentTransformFile)s " \
                        "-ImagesAreBinned %(ImagesAreBinned)s "

        if self.patchLayout.get() == 0:
            patchesXY = self.overlapPatches.getListFromValues(caster=float)
            OverlapOfPatchesXandY = '%f,%f' % (patchesXY[0], patchesXY[1])
            argsTiltXCorr += ' -OverlapOfPatchesXandY %s ' % OverlapOfPatchesXandY
        else:
            numberPatchesXY = self.numberOfPatches.getListFromValues()
            argsTiltXCorr += ' -NumberOfPatchesXandY %d,%d ' % (numberPatchesXY[0], numberPatchesXY[1])

        Plugin.runImod(self, 'tiltxcorr', argsTiltXCorr % paramsTiltXCorr)

    def chopcontsStep(self, tsId):
        """
        $imodchopconts -StandardInput
        InputModel cryo_pt.fid
        OutputModel cryo.fid
        MinimumOverlap	4
        AssignSurfaces 1
        """
        MinimumOverlap = 4
        AssignSurfaces = 1
        LengthOfPieces = -1

        paramschopconts = {
            'inputFile': self.getExtraOutFile(tsId, suffix="pt", ext=".fid"),
            'outputFile': self.getExtraOutFile(tsId, suffix="gaps", ext=".fid"),
            'MinimumOverlap': MinimumOverlap,
            'AssignSurfaces': AssignSurfaces,
            'LengthOfPieces': LengthOfPieces
        }
        argschopconts = " " \
                        "-InputModel %(inputFile)s " \
                        "-OutputModel %(outputFile)s " \
                        "-MinimumOverlap %(MinimumOverlap)s " \
                        "-AssignSurfaces %(AssignSurfaces)s " \
                        "-LengthOfPieces %(LengthOfPieces)s "

        Plugin.runImod(self, 'imodchopconts', argschopconts % paramschopconts)

    def createOutputStep(self):
        if self.FiducialModelGaps:
            self.FiducialModelGaps.setStreamState(Set.STREAM_CLOSED)
        if self.FailedTiltSeries:
            self.FailedTiltSeries.setStreamState(Set.STREAM_CLOSED)

        self._store()

    def createOutputFailedSetStep(self, tsId):
        ts = self.tsDict[tsId]
        super().createOutputFailedSet(ts)

    # --------------------------- UTILS functions -----------------------------
    def translateTrackCom(self, ts, paramsDict):
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        firstItem = ts.getFirstItem()
        trackFilePath = os.path.join(extraPrefix,
                                     firstItem.parseFileName(suffix="_track",
                                                             extension=".com"))

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
LowPassCutoffInverseNm  0.3
ImageFile	%(imageFile)s
ImagesAreBinned	1
InputSeedModel	%(inputSeedModel)s
OutputModel	%(outputModel)s
RotationAngle	%(rotationAngle).2f 
TiltFile	%(tiltFile)s
TiltDefaultGrouping	7
MagDefaultGrouping	5
RotDefaultGrouping	1
PixelSize   %(samplingRate)f
BeadDiameter	%(fiducialDiameter).2f
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
"""
        # Excluded views
        excludedViews = ts.getExcludedViewsIndex(caster=str)
        if len(excludedViews):
            template += f"SkipViews {','.join(excludedViews)}"

        if firstItem.hasTransform():
            XfFileName = self.getExtraOutFile(tsId, ext=XF_EXT)
            template += f"PrealignTransformFile {XfFileName}"

        with open(trackFilePath, 'w') as f:
            f.write(template % paramsDict)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.FiducialModelGaps:
            summary.append("Input tilt-series: %d\nFiducial models "
                           "(with gaps) generated: %d"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.FiducialModelGaps.getSize()))

        if self.FailedTiltSeries:
            summary.append("Failed tilt-series: %d"
                           % (self.FailedTiltSeries.getSize()))

        if not summary:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.FiducialModelGaps:
            methods.append("The fiducial model (with gaps) has "
                           "been computed for %d "
                           "tilt-series using the IMOD *beadtrack* command."
                           % (self.FiducialModelGaps.getSize()))

        if self.FailedTiltSeries:
            methods.append("%d tilt-series have failed during the "
                           "protocol execution."
                           % (self.FailedTiltSeries.getSize()))

        return methods
