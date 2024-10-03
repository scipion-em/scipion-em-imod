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

import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import STEPS_SERIAL
import pyworkflow.utils.path as path
import tomo.objects as tomoObj

from imod import utils
from imod.protocols import ProtImodBase
from imod.constants import (TLT_EXT, XF_EXT, FID_EXT, TXT_EXT, SEED_EXT,
                            SFID_EXT, OUTPUT_FIDUCIAL_GAPS_NAME,
                            FIDUCIAL_MODEL, PATCH_TRACKING, PT_FRACTIONAL_OVERLAP, PT_NUM_PATCHES)
from imod.protocols.protocol_base import IN_TS_SET
from pyworkflow.object import Set
from pyworkflow.utils import Message


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
    _possibleOutputs = {OUTPUT_FIDUCIAL_GAPS_NAME: tomoObj.SetOfLandmarkModels}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.stepsExecutionMode = STEPS_SERIAL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)

        form.addParam('typeOfModel',
                      params.EnumParam,
                      choices=["Make seed and Track", "Patch Tracking"],
                      default=FIDUCIAL_MODEL,
                      important=True,
                      label='Model generation')

        form.addParam(IN_TS_SET,
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt Series')

        self._patchTrackingForm(form, 'typeOfModel == %i' % PATCH_TRACKING)
        self._fiducialSeedForm(form, 'typeOfModel == %i' % FIDUCIAL_MODEL)

    def _patchTrackingForm(self, form, condition, levelType=params.LEVEL_NORMAL):
        patchtrack = form.addGroup('Patch Tracking',
                                   expertLevel=levelType,
                                   condition=condition)

        patchtrack.addParam('sizeOfPatches', params.NumericListParam,
                            label='Size of the patches (X,Y)',
                            default='680 680',
                            expertLevel=levelType,
                            help="Size of the  patches to track by correlation. "
                                 "In imod documentation "
                                 "(tiltxcorr: SizeOfPatchesXandY)")

        patchtrack.addParam('patchLayout',
                            params.EnumParam,
                            choices=['Fractional overlap of patches',
                                     'Number of patches'],
                            default=PT_FRACTIONAL_OVERLAP,
                            label='Patch layout',
                            display=params.EnumParam.DISPLAY_HLIST)
        # TODO: # help='To be added')

        patchtrack.addParam('overlapPatches',
                            params.NumericListParam,
                            default='0.33 0.33',
                            condition='patchLayout == %i' % PT_FRACTIONAL_OVERLAP,
                            label='Fractional overlap of the patches (X,Y)',
                            help="Fractional overlap in X and Y to track by correlation. "
                                 "In imod documentation"
                                 "(tiltxcorr: OverlapOfPatchesXandY)")

        patchtrack.addParam('numberOfPatches',
                            params.NumericListParam,
                            condition='patchLayout == %i' % PT_NUM_PATCHES,
                            label='Number of patches (X,Y)',
                            help="Number of patches in X and Y of the patches. "
                                 "In imod documentation"
                                 "(tiltxcorr: NumberOfPatchesXandY)")

        patchtrack.addParam('iterationsSubpixel',
                            params.IntParam,
                            default=1,
                            label='Iterations to increase subpixel accuracy',
                            help="Number of iteration of each correlation to reduce "
                                 "interpolation of the peak position"
                                 "In imod documentation: (tiltxcorr: IterateCorrelations)")

        self.addTrimingParams(patchtrack, pxTrimCondition=True,
                              correlationCondition=True,
                              levelType=params.LEVEL_ADVANCED)
        self.filteringParametersForm(form, condition=condition,
                                     levelType=params.LEVEL_ADVANCED)
        form.getParam("filterRadius2").setDefault(0.125)
        form.getParam("filterSigma2").setDefault(0.03)

    @staticmethod
    def _fiducialSeedForm(form, condition, levelType=params.LEVEL_NORMAL):
        seedModel = form.addGroup('"Make seed and Track', expertLevel=levelType, condition=condition)
        seedModel.addParam('fiducialDiameter',
                           params.FloatParam,
                           label='Fiducial diameter (nm)',
                           default=10.,
                           important=True,
                           help="Fiducials diameter to be tracked for alignment.")

        seedModel.addParam('twoSurfaces',
                           params.BooleanParam,
                           default=False,
                           label='Find beads on two surfaces?',
                           help="Track fiducials differentiating in which side "
                                "of the sample are located.")

        seedModel.addParam('numberFiducial',
                           params.IntParam,
                           label='Number of fiducials',
                           default=10,
                           help="Number of fiducials to be tracked for alignment.")

        seedModel.addParam('doTrackWithModel', params.BooleanParam,
                           default=True,
                           label="Track with fiducial model as seed",
                           help="Turn the tracked model into new seed and "
                                "repeat tracking.")

        groupGlobalVariables = form.addGroup('Filter variables',
                                             expertLevel=params.LEVEL_ADVANCED, condition=condition)

        groupGlobalVariables.addParam('refineSobelFilter',
                                      params.BooleanParam,
                                      default=True,
                                      label='Refine center with Sobel filter?',
                                      expertLevel=params.LEVEL_ADVANCED,
                                      help='Use edge-detecting Sobel filter '
                                           'to refine the bead positions.')

        groupGlobalVariables.addParam('scalableSigmaForSobelFilter',
                                      params.FloatParam,
                                      default=0.12,
                                      condition='refineSobelFilter',
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

            if self.typeOfModel.get() == FIDUCIAL_MODEL:
                self._insertFunctionStep(self.generateTrackComStep, tsId)
                self._insertFunctionStep(self.generateFiducialSeedStep, tsId)
                self._insertFunctionStep(self.generateFiducialModelStep, tsId)
            else:
                self._insertFunctionStep(self.xcorrStep, tsId)
                self._insertFunctionStep(self.chopcontsStep, tsId)

            self._insertFunctionStep(self.translateFiducialPointModelStep, tsId)
            self._insertFunctionStep(self.computeOutputModelsStep, tsId)

        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        tsSet = self.getInputSet()
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in tsSet}
        self.sRate = tsSet.getSamplingRate()
        self.acq = tsSet.getAcquisition()

    def generateTrackComStep(self, tsId):
        ts = self.tsDict[tsId]
        fiducialDiameterPixel, boxSizeXandY, scaling = self.getFiducialParams()
        paramsDict = {
            'imageFile': self.getTmpOutFile(tsId),
            'inputSeedModel': self.getExtraOutFile(tsId, ext=SEED_EXT),
            'outputModel': self.getExtraOutFile(tsId, suffix="gaps", ext=FID_EXT),
            'tiltFile': self.getExtraOutFile(tsId, ext=TLT_EXT),
            'rotationAngle': self.acq.getTiltAxisAngle(),
            'fiducialDiameter': fiducialDiameterPixel,
            'samplingRate': self.sRate / 10,
            'scalableSigmaForSobelFilter': self.scalableSigmaForSobelFilter.get(),
            'boxSizeXandY': boxSizeXandY,
            'distanceRescueCriterion': 10 * scaling,
            'postFitRescueResidual': 2.5 * scaling,
            'maxRescueDistance': 2.5 * scaling,
            'minDiamForParamScaling': 12.5,
            'deletionCriterionMinAndSD': f"{0.04 * scaling:0.3f},2.0",
        }

        excludedViews = ts.getExcludedViewsIndex(caster=str)
        hasAlign = ts.hasAlignment()
        self.translateTrackCom(tsId, paramsDict, excludedViews, hasAlign)

    def generateFiducialSeedStep(self, tsId):
        try:
            trackFile = self.getExtraOutFile(tsId, suffix="track", ext="com")
            if os.path.exists(trackFile):
                paramsAutofidseed = {
                    "-TrackCommandFile": trackFile,
                    "-MinSpacing": 0.85,
                    "-AdjustSizes": "",
                    "-PeakStorageFraction": 1.0,
                    "-TargetNumberOfBeads": self.numberFiducial.get(),
                }

                if self.twoSurfaces:
                    paramsAutofidseed["-TwoSurfaces"] = ""

                self.runProgram('autofidseed', paramsAutofidseed)

                autofidseedDirPath = self._getExtraPath(tsId, "autofidseed.dir")
                path.makePath(autofidseedDirPath)
                path.moveTree("autofidseed.dir", autofidseedDirPath)
                path.moveFile("autofidseed.info", self._getExtraPath(tsId))
        except Exception as e:
            self._failedItems.append(tsId)
            self.error(f'autofidseed execution failed for tsId {tsId} -> {e}')

    def generateFiducialModelStep(self, tsId):
        ts = self.tsDict[tsId]
        if tsId not in self._failedItems:
            try:
                fiducialDiameterPixel, boxSizeXandY, scaling = self.getFiducialParams()

                paramsBeadtrack = {
                    "-ImageFile": self.getTmpOutFile(tsId),
                    "-InputSeedModel": self.getExtraOutFile(tsId, ext=SEED_EXT),
                    "-OutputModel": self.getExtraOutFile(tsId, suffix="gaps", ext=FID_EXT),
                    "-TiltFile": self.getExtraOutFile(tsId, ext=TLT_EXT),
                    "-RotationAngle": ts.getAcquisition().getTiltAxisAngle(),
                    "-PixelSize": self.sRate / 10,
                    "-ImagesAreBinned": 1,
                    "-TiltDefaultGrouping": 7,
                    "-MagDefaultGrouping": 5,
                    "-RotDefaultGrouping": 1,
                    "-MinViewsForTiltalign": 4,
                    "-BeadDiameter": fiducialDiameterPixel,
                    "-FillGaps": 1,
                    "-MaxGapSize": 5,
                    "-MinTiltRangeToFindAxis": 10.0,
                    "-MinTiltRangeToFindAngles": 20.0,
                    "-BoxSizeXandY": f"{boxSizeXandY},{boxSizeXandY}",
                    "-RoundsOfTracking": 2,
                    "-LocalAreaTracking": 1,
                    "-LocalAreaTargetSize": 1000,
                    "-MinBeadsInArea": 8,
                    "-MinOverlapBeads": 5,
                    "-MaxBeadsToAverage": 4,
                    "-PointsToFitMaxAndMin": '7,3',
                    "-DensityRescueFractionAndSD": '0.6,1.0',
                    "-DistanceRescueCriterion": 10 * scaling,
                    "-RescueRelaxationDensityAndDistance": '0.7,0.9',
                    "-PostFitRescueResidual": 2.5 * scaling,
                    "-DensityRelaxationPostFit": 0.9,
                    "-MaxRescueDistance": 2.5 * scaling,
                    "-ResidualsToAnalyzeMaxAndMin": '9,5',
                    "-DeletionCriterionMinAndSD": f"{0.04 * scaling:0.3f},2.0",
                    "-MinDiamForParamScaling": 12.5,
                    "-LowPassCutoffInverseNm": 0.3
                }

                if self.refineSobelFilter:
                    paramsBeadtrack["-SobelFilterCentering"] = ""
                    paramsBeadtrack["-ScalableSigmaForSobel"] = self.scalableSigmaForSobelFilter.get()

                # Excluded views
                excludedViews = ts.getExcludedViewsIndex(caster=str)
                if len(excludedViews):
                    paramsBeadtrack["-SkipViews"] = ",".join(excludedViews)

                if ts.hasAlignment():
                    paramsBeadtrack["-prexf"] = self.getExtraOutFile(tsId, ext=XF_EXT)

                self.runProgram('beadtrack', paramsBeadtrack)

                if self.doTrackWithModel:
                    # repeat tracking with the current model as seed
                    path.copyFile(paramsBeadtrack['-InputSeedModel'],
                                  self.getExtraOutFile(tsId, suffix="orig", ext=SEED_EXT))
                    path.moveFile(paramsBeadtrack['-OutputModel'],
                                  paramsBeadtrack['-InputSeedModel'])

                    self.runProgram('beadtrack', paramsBeadtrack)

            except Exception as e:
                self._failedItems.append(tsId)
                self.error(f'beadtrack execution failed for tsId {tsId} -> {e}')

    def xcorrStep(self, tsId):
        try:
            ts = self.tsDict[tsId]
            angleFilePath = self.getExtraOutFile(tsId, ext=TLT_EXT)
            xfFile = self.getExtraOutFile(tsId, ext=XF_EXT)
            ts.writeXfFile(xfFile)

            borders = self.pxTrim.getListFromValues(caster=str)
            sizePatches = self.sizeOfPatches.getListFromValues(caster=str)

            paramsTiltXCorr = {
                "-InputFile": self.getTmpOutFile(tsId),
                "-OutputFile": self.getExtraOutFile(tsId, suffix="pt", ext=FID_EXT),
                "-RotationAngle": self.acq.getTiltAxisAngle(),
                "-TiltFile": angleFilePath,
                "-FilterRadius2": self.filterRadius2.get(),
                "-FilterSigma1": self.filterSigma1.get(),
                "-FilterSigma2": self.filterSigma2.get(),
                "-BordersInXandY": ",".join(borders),
                "-IterateCorrelations": self.iterationsSubpixel.get(),
                "-SizeOfPatchesXandY": ",".join(sizePatches),
                "-PrealignmentTransformFile": xfFile,
                "-ImagesAreBinned": 1,
            }

            if self.patchLayout.get() == PT_FRACTIONAL_OVERLAP:
                patchesXY = self.overlapPatches.getListFromValues(caster=str)
                paramsTiltXCorr["-OverlapOfPatchesXandY"] = ",".join(patchesXY)
            else:
                numberPatchesXY = self.numberOfPatches.getListFromValues(caster=str)
                paramsTiltXCorr["-NumberOfPatchesXandY"] = ",".join(numberPatchesXY)

            self.runProgram('tiltxcorr', paramsTiltXCorr)

        except Exception as e:
            self._failedItems.append(tsId)
            self.error(f'tiltxcorr execution failed for tsId {tsId} -> {e}')

    def chopcontsStep(self, tsId):
        if tsId not in self._failedItems:
            try:
                paramschopconts = {
                    "-InputModel": self.getExtraOutFile(tsId, suffix="pt", ext=FID_EXT),
                    "-OutputModel": self.getExtraOutFile(tsId, suffix="gaps", ext=FID_EXT),
                    "-MinimumOverlap": 4,
                    "-AssignSurfaces": 1,
                    "-LengthOfPieces": -1
                }
                self.runProgram('imodchopconts', paramschopconts)
            except Exception as e:
                self._failedItems.append(tsId)
                self.error(f'imodchopconts execution failed for tsId {tsId} -> {e}')

    def translateFiducialPointModelStep(self, tsId):
        if tsId not in self._failedItems:
            gapsFidFile = self.getExtraOutFile(tsId, suffix='gaps', ext=FID_EXT)

            if os.path.exists(gapsFidFile):
                paramsGapModel2Point = {
                    "-InputFile": gapsFidFile,
                    "-OutputFile": self.getExtraOutFile(tsId, suffix="gaps_fid", ext=TXT_EXT)
                }
                self.runProgram('model2point', paramsGapModel2Point)

    def computeOutputModelsStep(self, tsId):
        """ Create the output set of landmark models with gaps. """
        ts = self.tsDict[tsId]
        if tsId in self._failedItems:
            self.createOutputFailedSet(ts)
        else:
            fiducialModelGapPath = self.getExtraOutFile(tsId, suffix='gaps', ext=FID_EXT)

            if os.path.exists(fiducialModelGapPath):
                output = self.getOutputFiducialModel(self.getInputSet(pointer=True),
                                                     attrName=OUTPUT_FIDUCIAL_GAPS_NAME, suffix="Gaps")
                landmarkModelGapsFilePath = self.getExtraOutFile(tsId, suffix='gaps', ext=SFID_EXT)
                fiducialModelGapTxtPath = self.getExtraOutFile(tsId, suffix="gaps_fid", ext=TXT_EXT)

                fiducialGapList = utils.formatFiducialList(fiducialModelGapTxtPath)
                fiducialDiameter = self.fiducialDiameter.get() * 10  # Angstroms

                landmarkModelGaps = tomoObj.LandmarkModel(tsId=tsId,
                                                          tiltSeriesPointer=ts,
                                                          fileName=landmarkModelGapsFilePath,
                                                          modelName=fiducialModelGapPath,
                                                          size=fiducialDiameter,
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
                output.write(output)
                self._store(output)
            else:
                self.createOutputFailedSet(ts)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        fidModelGaps = getattr(self, OUTPUT_FIDUCIAL_GAPS_NAME, None)
        if fidModelGaps is not None:
            summary.append(f"Input tilt-series: {self.getInputSet().getSize()}\n"
                           "Fiducial models (with gaps) generated: "
                           f"{fidModelGaps.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []

        fidModelGaps = getattr(self, OUTPUT_FIDUCIAL_GAPS_NAME, None)
        if fidModelGaps is not None:
            methods.append("The fiducial model (with gaps) has been computed for "
                           f"{fidModelGaps.getSize()} tilt-series using "
                           "the IMOD *beadtrack* command.")

        return methods

    # --------------------------- UTILS functions -----------------------------
    def translateTrackCom(self, tsId, paramsDict, excludedViews,
                          hasAlignment=False):
        trackFilePath = self.getExtraOutFile(tsId, suffix="track", ext="com")

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

        if self.refineSobelFilter:
            template += "\nSobelFilterCentering"
            template += "\nScalableSigmaForSobel   %(scalableSigmaForSobelFilter)f"

        if len(excludedViews):
            template += f"\nSkipViews {','.join(excludedViews)}"

        if hasAlignment:
            XfFileName = self.getExtraOutFile(tsId, ext=XF_EXT)
            template += f"\nPrealignTransformFile {XfFileName}"

        with open(trackFilePath, 'w') as f:
            f.write(template % paramsDict)

    def getFiducialParams(self):
        """ Precalculate a few params. """
        fiducialDiameterPixel = self.fiducialDiameter.get() / (self.sRate/10)
        # formulas from bin/copytomos
        boxSizeXandY = max(3.3 * fiducialDiameterPixel + 2, 2 * fiducialDiameterPixel + 20, 32)
        boxSizeXandY = min(512, 2 * int(boxSizeXandY / 2))
        scaling = fiducialDiameterPixel / 12.5 if fiducialDiameterPixel > 12.5 else 1

        return fiducialDiameterPixel, boxSizeXandY, scaling

    def getOutputFiducialModelGaps(self, inputSet):
        if self.FiducialModelGaps:
            self.FiducialModelGaps.enableAppend()
        else:
            fidModelGaps = self._createSetOfLandmarkModels(suffix='Gaps')

            fidModelGaps.copyInfo(inputSet)
            fidModelGaps.setSetOfTiltSeries(inputSet)
            fidModelGaps.setHasResidualInfo(False)
            fidModelGaps.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_FIDUCIAL_GAPS_NAME: fidModelGaps})
            self._defineSourceRelation(inputSet, fidModelGaps)

        return self.FiducialModelGaps
