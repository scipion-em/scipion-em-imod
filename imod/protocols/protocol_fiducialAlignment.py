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
import numpy as np

import pyworkflow.protocol.params as params
from imod.protocols.protocol_base import BINNING_FACTOR
from pwem import ALIGN_NONE
from pyworkflow.protocol.constants import STEPS_SERIAL
from pwem.objects import Transform
from pyworkflow.utils import Message
from tomo.objects import (LandmarkModel, SetOfLandmarkModels, SetOfTiltSeries,
                          TiltSeries, TiltSeriesCoordinate)

from imod import utils
from imod.protocols import ProtImodBase
from imod.constants import (TLT_EXT, XF_EXT, FID_EXT, TXT_EXT, XYZ_EXT,
                            MOD_EXT, SFID_EXT, OUTPUT_TILTSERIES_NAME,
                            OUTPUT_FIDUCIAL_NO_GAPS_NAME,
                            OUTPUT_TS_INTERPOLATED_NAME,
                            OUTPUT_TS_COORDINATES_NAME)

# Rotation solution types
NO_ROTATION = 0
ONE_ROTATION = 1
GROUP_ROTATIONS = 2
ALL_ROTATIONS = 3
ROT_SOLUTION_CHOICES = ['No rotation',
                        'One rotation',
                        'Group rotations',
                        'Solve for all rotations']

# Magnification solution types
FIXED_MAG = 0
GROUP_MAGS = 1
ALL_MAGS = 2
MAG_SOLUTION_CHOICES = ['Fixed magnification at 1.0',
                        'Group magnifications',
                        'Solve for all magnifications']

# Tilt angle solution types
FIXED_TILT = 0
GROUP_TILTS = 1
ALL_EXCEPT_MIN = 2
TILT_SOLUTION_CHOICES = ['Fixed tilt angles',
                         'Group tilt angles',
                         'Solve for all except minimum tilt']

# Distortion solution types
DIST_DISABLED = 0
DIST_FULL_SOLUTION = 1
DIST_SKEW_ONLY = 2
DISTORTION_SOLUTION_CHOICES = ['Disabled',
                               'Full solution',
                               'Skew only']


class ProtImodFiducialAlignment(ProtImodBase):
    """
    Construction of a fiducial model and alignment of tilt-series based
    on the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/tiltalign.html

    This program will solve for the displacements, rotations, tilts, and
    magnification differences relating a set of tilted views of an object.
    It uses a set of fiducial points that have been identified in a series
    of views. These input data are read from a model in which each fiducial
    point is a separate contour.

    This program has several notable features:

    1) Any given fiducial point need not be present in every view. Thus,
    one can track each fiducial point only through the set of views in
    which it can be reliably identified, and one can even skip views in the
    middle of that set.

    2) The program can solve for distortion (stretching) in the plane of
    the section.

    3) It is possible to constrain several views to have the same unknown
    value of rotation, tilt angle, magnification, compression, or distor-
    tion.  This can reduce the number of unknowns and can give more accu-
    rate overall solutions.

    4) If the fiducial points are supposed to lie in one or two planes,
    then after the minimization procedure is complete, the program can ana-
    lyze the solved point positions and determine the slope of this plane.
    It uses this slope to estimate how to adjust tilt angles so as to make
    the planes be horizontal in a reconstruction.

    5) The program can use a robust fitting method to give different
    weights to different modeled points based on their individual fitting
    errors.  Points with the most extreme errors are eliminated from the
    fit, and ones with high but less extreme errors are down-weighted.
    This fitting provides a substitute for fixing many modeled point posi-
    tions based on their errors.

    _On the alignment model_:\n
    The program implements the following model for the imaging of the spec-
    imen in each individual view:
       1) The specimen itself changes by
         a) an isotropic size change (magnification variable);
         b) additional thinning in the Z dimension (compression variable); and
         c) linear stretch along one axis in the specimen plane, implemented by
            variables representing stretch along the X axis and skew between
            the X and Y axes;
       2) The specimen is tilted slightly around the X axis (X tilt variable)
       3) The specimen is tilted around the X axis by the negative of the beam
          tilt, if any (one variable for all views)
       4) The specimen is tilted around the Y axis (tilt variable)
       5) The specimen is tilted back around the X axis by the beam tilt, if any
       6) The projected image rotates in the plane of the camera (rotation
          variable)
       7) The projected image may stretch along an axis midway between the
          original X and Y axes (one variable for all views)
       8) The image shifts on the camera

    The complete model is summarized in:
       Mastronarde, D. N. 2008.  Correction for non-perpendicularity of beam
       and tilt axis in tomographic reconstructions with the IMOD package. J.
       Microsc.  230: 212-217.
    The version of the model prior to the addition of beam tilt is described
    in more detail in:
       Mastronarde, D. N.  2007.  Fiducial marker and hybrid alignment methods
       for single- and double-axis tomography.  In: Electron Tomography, Ed.
       J. Frank, 2nd edition, pp 163-185. Springer, New York.
    """

    _label = 'Fiducial alignment'
    _possibleOutputs = {
        OUTPUT_TILTSERIES_NAME: SetOfTiltSeries,
        OUTPUT_FIDUCIAL_NO_GAPS_NAME: SetOfLandmarkModels
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.stepsExecutionMode = STEPS_SERIAL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)

        form.addParam('inputSetOfLandmarkModels',
                      params.PointerParam,
                      pointerClass='SetOfLandmarkModels',
                      important=True,
                      label='Fiducial model')

        form.addParam('twoSurfaces',
                      params.BooleanParam,
                      default=False,
                      label='Assume beads on two surfaces?',
                      help="Track fiducials differentiating in which side of the sample are located.\n"
                           "IMPORTANT: It is highly recommended to match the option selected in the "
                           "generation of the fiducial models. In case they  do not match, it is not "
                           "intended to fail but could be missing the whole potential of the algorithm. "
                           "In case the algorithm used for the calculation of the fiducial models does "
                           "not consider this option it is algo recommended to set this "
                           "option to 'No'.")

        form.addParam('computeAlignment',
                      params.BooleanParam,
                      default=False,
                      label='Generate interpolated tilt-series?',
                      important=True,
                      help='Generate and save the interpolated tilt-series applying the obtained transformation '
                           'matrices.\n'
                           'By default, the output of this protocol will be a tilseries that will have associated'
                           'the alignment information as a transformation matrix. When this option is set as Yes, '
                           'then a second output, called interpolated tilt series, is generated. The interpolated tilt '
                           'series should be used for visualization purpose but not for image processing')

        groupInterpolation = form.addGroup('Interpolated tilt-series',
                                           condition='computeAlignment')

        groupInterpolation.addParam(BINNING_FACTOR,
                                    params.IntParam,
                                    default=1,
                                    label='Binning for the interpolated',
                                    help='Binning to be applied to the interpolated  tilt-series in IMOD '
                                         'convention. \n'
                                         'Binning is an scaling factor given by an integer greater than 1. '
                                         'IMOD uses ordinary binning to reduce images in size by the given factor. '
                                         'The value of a binned pixel is the average of pixel values in each block '
                                         'of pixels being binned. Binning is applied before all other image '
                                         'transformations.')

        form.addSection('Global variables')

        form.addParam('rotationSolutionType',
                      params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=ROT_SOLUTION_CHOICES,
                      default=GROUP_ROTATIONS,
                      label='Rotation solution type',
                      help='Type of rotation solution: See rotOption in tiltalign IMOD command \n'
                           '* No rotation: The in-plane rotation will not be estimated\n'
                           '* One rotation: To solve for a single rotation variable \n'
                           '* Group rotations: Group views to solve for fewer rotations variables. Automapping of '
                           'rotation variables linearly changing values\n'
                           '* Solve for all rotations: for each view having an independent rotation\n')

        form.addParam('groupRotationSize',
                      params.IntParam,
                      default=5,
                      condition='rotationSolutionType == %i' % GROUP_ROTATIONS,
                      label='Group size',
                      help='Default group size when automapping rotation variables')

        form.addParam('magnificationSolutionType',
                      params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=MAG_SOLUTION_CHOICES,
                      default=FIXED_MAG,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Magnification solution type',
                      help='Type of magnification solution: See MagOption in tiltaling IMOD command\n'
                           '* Fixed magnification: Do not solve magnification. This fixes all magnifications at 1.0.\n'
                           '* Group magnifications: Group views to solve for fewer magnifications variables. '
                           'Automapping of variables (linearly changing values)  \n'
                           '* Solve for all magnifications: to vary all magnifications  of each view independently\n')

        form.addParam('groupMagnificationSize',
                      params.IntParam,
                      default=4,
                      condition='magnificationSolutionType == %i' % GROUP_MAGS,
                      label='Group size',
                      help='Group size when automapping magnification variables')

        form.addParam('tiltAngleSolutionType',
                      params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=TILT_SOLUTION_CHOICES,
                      default=GROUP_TILTS,
                      label='Tilt angle solution type',
                      help='Type of tilt angle solution: See TiltOption in tiltalign IMOD command\n'
                           ' * Fixed tilt angles: To fix all tilt angles at their initial (input) values \n'
                           ' * Group tilt angles: To automap groups of tilt angles (linearly changing values) \n'
                           ' * Solve for all except minimum tilt:to solve for all tilt angles except for the view '
                           'at minimum tilt \n')

        form.addParam('groupTiltAngleSize',
                      params.IntParam,
                      default=5,
                      condition='tiltAngleSolutionType == %i' % GROUP_TILTS,
                      label='Group size',
                      help='Average default group size when automapping tilt variables')

        form.addParam('distortionSolutionType',
                      params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=DISTORTION_SOLUTION_CHOICES,
                      default=DIST_DISABLED,
                      label='Distortion solution type',
                      help='Type of skew solution:'
                           '* 0 to fix all skew angles at 0.0 \n'
                           '* 1 to vary all skew angles independently\n '
                           '* 2 to specify a mapping of skew variables, or\n '
                           '* 3 or 4 for automapping of variables (3 for linearly changing values or 4 for values all '
                           'the same within a group)..')

        form.addParam('xStretchGroupSize',
                      params.IntParam,
                      default=7,
                      condition='distortionSolutionType == %i' % DIST_FULL_SOLUTION,
                      label='X stretch group size',
                      help='Basic grouping size for X stretch')

        form.addParam('skewGroupSize',
                      params.IntParam,
                      default=11,
                      condition='distortionSolutionType in [%i, %i]' % (DIST_FULL_SOLUTION, DIST_SKEW_ONLY),
                      label='Skew group size',
                      help='Size of the skew group')

    # NOTE:
    # The gold bead eraser only remove the fiducial markers in the interpolated TS and only the ones used for the
    #  TS alignment, so for now this functionality will be removed. Jorge (07/08/2024)
    #
    # form.addSection('Erase gold beads')
    #
    # form.addParam('eraseGoldBeads',
    #               params.BooleanParam,
    #               default=False,
    #               label='Erase gold beads',
    #               help='Remove the gold beads detected during fiducial '
    #                    'alignment with *ccderaser* program. This option '
    #                    'will generate an interpolated tilt series with '
    #                    'the gold beads erased and interpolated with '
    #                    'the calculated transformation matrices form '
    #                    'the alignment.')
    #
    # form.addParam('betterRadius',  # actually diameter
    #               params.IntParam,
    #               default=18,
    #               label='Bead diameter (px)',
    #               help="For circle objects, this entry "
    #                    "specifies a radius to use for points "
    #                    "without an individual point size "
    #                    "instead of the object's default sphere "
    #                    "radius. This entry is floating point "
    #                    "and can be used to overcome the "
    #                    "limitations of having an integer "
    #                    "default sphere radius. If there are "
    #                    "multiple circle objects, enter one "
    #                    "value to apply to all objects or a "
    #                    "value for each object.")

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        binning = self.binning.get()

        for tsId in self.tsDict.keys():
            self.fiducialDiameterPixel = self.lmDict[tsId][0]
            self._insertFunctionStep(self.convertInputStep, tsId)
            self._insertFunctionStep(self.computeFiducialAlignmentStep, tsId)
            self._insertFunctionStep(self.translateFiducialPointModelStep, tsId)
            self._insertFunctionStep(self.computeOutputStackStep, tsId)

            if self.computeAlignment:  # or self.eraseGoldBeads:
                self._insertFunctionStep(self.computeOutputInterpolatedStackStep,
                                         tsId, binning)
            # if self.eraseGoldBeads:
            #     self._insertFunctionStep(self.eraseGoldBeadsStep, tsId)
            self._insertFunctionStep(self.computeOutputModelsStep, tsId)

        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        inFiduSet = self.getInputFiduSet()
        self.inTsSetPointer = inFiduSet.getSetOfTiltSeries(pointer=True)
        # There can be failed fiducial models, so the TsIds used as reference must be the ones present in the input set
        # of landmark models
        fidTsIds = inFiduSet.getUniqueValues(TiltSeries.TS_ID_FIELD)
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self.inTsSetPointer.get()
                       if ts.getTsId() in fidTsIds}

        lms = inFiduSet.aggregate(["COUNT"], TiltSeries.TS_ID_FIELD,
                                  [TiltSeries.TS_ID_FIELD, "_size", "_modelName"])
        self.lmDict = {lm[TiltSeries.TS_ID_FIELD]: (lm["_size"], lm["_modelName"]) for lm in lms}

    def computeFiducialAlignmentStep(self, tsId):
        try:
            ts = self.tsDict[tsId]

            paramsTiltAlign = {
                "-ModelFile": self.lmDict[tsId][1],
                "-ImageFile": self.getTmpOutFile(tsId),
                "-ImagesAreBinned": 1,
                "-UnbinnedPixelSize": ts.getSamplingRate() / 10,
                "-OutputModelFile": self.getExtraOutFile(tsId, suffix="fidxyz", ext=MOD_EXT),
                "-OutputResidualFile": self.getExtraOutFile(tsId, suffix="resid", ext=TXT_EXT),
                "-OutputFidXYZFile": self.getExtraOutFile(tsId, suffix="fid", ext=XYZ_EXT),
                "-OutputTiltFile": self.getExtraOutFile(tsId, suffix="interpolated", ext=TLT_EXT),
                "-OutputXAxisTiltFile": self.getExtraOutFile(tsId, ext="xtilt"),
                "-OutputTransformFile": self.getExtraOutFile(tsId, suffix="fid", ext=XF_EXT),
                "-OutputFilledInModel": self.getExtraOutFile(tsId, suffix="noGaps", ext=FID_EXT),
                "-RotationAngle": ts.getAcquisition().getTiltAxisAngle(),
                "-TiltFile": self.getExtraOutFile(tsId, ext=TLT_EXT),
                "-AngleOffset": 0.0,
                "-RotOption": self.getRotationType(),
                "-RotDefaultGrouping": self.groupRotationSize.get(),
                "-TiltOption": self.getTiltAngleType(),
                "-TiltDefaultGrouping": self.groupTiltAngleSize.get(),
                "-MagReferenceView": 1,
                "-MagOption": self.getMagnificationType(),
                "-MagDefaultGrouping": self.groupMagnificationSize.get(),
                "-XStretchOption": self.getStretchType(),
                "-SkewOption": self.getSkewType(),
                "-XStretchDefaultGrouping": self.xStretchGroupSize.get(),
                "-SkewDefaultGrouping": self.skewGroupSize.get(),
                "-BeamTiltOption": 0,
                "-XTiltOption": 0,
                "-XTiltDefaultGrouping": 2000,
                "-ResidualReportCriterion": 3.0,
                "-SurfacesToAnalyze": self.getSurfaceToAnalyze(),
                "-MetroFactor": 0.25,
                "-MaximumCycles": 1000,
                "-KFactorScaling": 1.0,
                "-NoSeparateTiltGroups": 1,
                "-AxisZShift": 0.0,
                "-ShiftZFromOriginal": 1,
                "-TargetPatchSizeXandY": '700,700',
                "-MinSizeOrOverlapXandY": '0.5,0.5',
                "-MinFidsTotalAndEachSurface": '8,3',
                "-FixXYZCoordinates": 0,
                "-RobustFitting": ""
            }

            # Excluded views
            excludedViews = ts.getExcludedViewsIndex(caster=str)
            if len(excludedViews):
                paramsTiltAlign["-ExcludeList"] = ",".join(excludedViews)

            paramsTiltAlign["2>&1 | tee "] = self._getExtraPath("align.log")

            self.runProgram('tiltalign', paramsTiltAlign)
            self.runProgram('alignlog', {'-s': "> taSolution.log"},
                            cwd=self._getExtraPath())

        except Exception as e:
            self._failedItems.append(tsId)
            self.error(f'tiltalign execution failed for tsId {tsId} -> {e}')

    def translateFiducialPointModelStep(self, tsId):
        # Check that previous steps have been completed satisfactorily
        if tsId not in self._failedItems:
            noGapsFid = self.getExtraOutFile(tsId, suffix="noGaps", ext=FID_EXT)
            if os.path.exists(noGapsFid):
                paramsNoGapModel2Point = {
                    "-InputFile": noGapsFid,
                    "-OutputFile": self.getExtraOutFile(tsId, suffix="noGaps_fid", ext=TXT_EXT)
                }
                self.runProgram('model2point', paramsNoGapModel2Point)

    def computeOutputStackStep(self, tsId):
        ts = self.tsDict[tsId]
        if tsId not in self._failedItems:
            tmFilePath = self.getExtraOutFile(tsId, suffix="fid", ext=XF_EXT)
            if os.path.exists(tmFilePath) and os.stat(tmFilePath).st_size != 0:
                tltFilePath = self.getExtraOutFile(tsId, suffix="interpolated", ext=TLT_EXT)
                tltList = utils.formatAngleList(tltFilePath)
                alignmentMatrix = utils.formatTransformationMatrix(tmFilePath)
                output = self.getOutputSetOfTS(self.inTsSetPointer)
                self.copyTsItems(output, ts, tsId,
                                 updateTsCallback=self.updateTsNonInterp,
                                 updateTiCallback=self.updateTiNonInterp,
                                 copyDisabledViews=True,
                                 copyId=True,
                                 copyTM=False,
                                 alignmentMatrix=alignmentMatrix,
                                 tltList=tltList)
            else:
                self.createOutputFailedSet(ts)

    def computeOutputInterpolatedStackStep(self, tsId, binning):
        """ Generate interpolated stack. """
        if tsId not in self._failedItems:
            tmpFileName = self.getExtraOutFile(tsId, suffix="fid", ext=XF_EXT)
            if os.path.exists(tmpFileName) and os.stat(tmpFileName).st_size != 0:
                ts = self.tsDict[tsId]
                output = self.getOutputSetOfTS(self.inTsSetPointer, binning,
                                               attrName=OUTPUT_TS_INTERPOLATED_NAME,
                                               suffix="Interpolated")

                firstItem = ts.getFirstItem()
                tsExcludedIndices = ts.getExcludedViewsIndex()
                paramsDict = self.getBasicNewstackParams(
                    ts,
                    self.getExtraOutFile(tsId),
                    inputTsFileName=self.getTmpOutFile(tsId),
                    xfFile=tmpFileName,
                    firstItem=firstItem,
                    binning=binning,
                    doSwap=True,
                    doTaper=True,
                    tsExcludedIndices=tsExcludedIndices,
                )
                self.runProgram('newstack', paramsDict)

                tltFilePath = self.getExtraOutFile(tsId, suffix="interpolated", ext=TLT_EXT)
                tltList = utils.formatAngleList(tltFilePath)
                self.copyTsItems(output, ts, tsId,
                                 updateTsCallback=self.updateTsInterp,
                                 updateTiCallback=self.updateTiInterp,
                                 copyId=True,
                                 copyTM=False,
                                 excludedViews=len(tsExcludedIndices) > 0,
                                 tltList=tltList)

    # def eraseGoldBeadsStep(self, tsId):
    #     """ Erase gold beads on aligned stack. """
    #     if tsId not in self._failedItems:
    #         try:
    #             paramsCcderaser = {
    #                 "-InputFile": self.getTmpOutFile(tsId),
    #                 "-OutputFile": self.getExtraOutFile(tsId),
    #                 "-ModelFile": self.getExtraOutFile(tsId, suffix="noGaps", ext=FID_EXT),
    #                 "-BetterRadius": self.betterRadius.get() / 2,
    #                 "-PolynomialOrder": 0,
    #                 "-CircleObjects": "/",
    #                 "-MergePatches": 1,
    #                 "-ExcludeAdjacent": "",
    #                 "-SkipTurnedOffPoints": 1,
    #                 "-ExpandCircleIterations": 3
    #             }
    #             self.runProgram('ccderaser', paramsCcderaser)
    #         except Exception as e:
    #             self._failedItems.append(tsId)
    #             self.error(f'ccderaser execution failed for tsId {tsId} -> {e}')

    def computeOutputModelsStep(self, tsId):
        """ Create output sets of landmarks and 3D coordinates. """
        ts = self.tsDict[tsId]
        if tsId in self._failedItems:
            self.createOutputFailedSet(ts)
        else:
            # Create the output set of landmark models with no gaps
            fiducialNoGapFilePath = self.getExtraOutFile(tsId, suffix="noGaps_fid", ext=TXT_EXT)
            if os.path.exists(fiducialNoGapFilePath):
                output = self.getOutputFiducialModel(self.inTsSetPointer)
                fiducialNoGapList = utils.formatFiducialList(fiducialNoGapFilePath)
                fiducialModelNoGapPath = self.getExtraOutFile(tsId, suffix="noGaps", ext=FID_EXT)
                landmarkModelNoGapsFilePath = self.getExtraOutFile(tsId, suffix="noGaps", ext=SFID_EXT)
                landmarkModelNoGapsResidPath = self.getExtraOutFile(tsId, suffix="resid", ext=TXT_EXT)

                fiducialNoGapsResidList = utils.formatFiducialResidList(landmarkModelNoGapsResidPath)

                landmarkModelNoGaps = LandmarkModel(tsId=tsId,
                                                    fileName=landmarkModelNoGapsFilePath,
                                                    modelName=fiducialModelNoGapPath,
                                                    size=self.fiducialDiameterPixel,
                                                    hasResidualInfo=True)

                prevTiltIm = 0
                chainId = 0
                indexFake = 0
                firstExec = True

                for fiducial in fiducialNoGapList:
                    if (int(float(fiducial[2])) <= prevTiltIm) or firstExec:
                        chainId += 1
                        firstExec = False
                    prevTiltIm = int(float(fiducial[2]))

                    if indexFake < len(fiducialNoGapsResidList) and fiducial[2] == fiducialNoGapsResidList[indexFake][2]:
                        landmarkModelNoGaps.addLandmark(xCoor=fiducial[0],
                                                        yCoor=fiducial[1],
                                                        tiltIm=fiducial[2] + 1,
                                                        chainId=chainId,
                                                        xResid=fiducialNoGapsResidList[indexFake][3],
                                                        yResid=fiducialNoGapsResidList[indexFake][4])
                        indexFake += 1

                    else:
                        landmarkModelNoGaps.addLandmark(xCoor=fiducial[0],
                                                        yCoor=fiducial[1],
                                                        tiltIm=fiducial[2] + 1,
                                                        chainId=chainId,
                                                        xResid=float('nan'),
                                                        yResid=float('nan'))

                output.append(landmarkModelNoGaps)
                output.update(landmarkModelNoGaps)
                self._store(output)

        # Create the output set of 3D coordinates
        coordFilePath = self.getExtraOutFile(tsId, suffix="fid", ext=XYZ_EXT)
        if os.path.exists(coordFilePath):
            output = self.getOutputSetOfTiltSeriesCoordinates(self.inTsSetPointer)
            coordList, xDim, yDim = utils.format3DCoordinatesList(coordFilePath)

            for element in coordList:
                newCoord3D = TiltSeriesCoordinate(tsId=tsId)
                newCoord3D.setPosition(element[0] - (xDim / 2),
                                       element[1] - (yDim / 2),
                                       element[2],
                                       sampling_rate=ts.getSamplingRate())
                output.append(newCoord3D)

            self._store(output)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        interpTS = getattr(self, OUTPUT_TS_INTERPOLATED_NAME, None)
        fidModelNoGaps = getattr(self, OUTPUT_FIDUCIAL_NO_GAPS_NAME, None)
        tsCoords = getattr(self, OUTPUT_TS_COORDINATES_NAME, None)

        if fidModelNoGaps is not None:
            summary.append("Fiducial models generated with no gaps: "
                           f"{fidModelNoGaps.getSize()}")

        if self.TiltSeries:
            summary.append("Transformation matrices updated from the "
                           f"input tilt-series: {self.TiltSeries.getSize()}")

        if interpTS is not None:
            summary.append("Interpolated tilt-series calculated: "
                           f"{interpTS.getSize()}")

        if tsCoords is not None:
            summary.append("Fiducial 3D coordinates calculated: "
                           f"{tsCoords.getSize()}")

        if not summary:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []

        fidModelNoGaps = getattr(self, OUTPUT_FIDUCIAL_NO_GAPS_NAME, None)
        if fidModelNoGaps is not None:
            methods.append("Solved fiducials alignment for "
                           f"{fidModelNoGaps.getSize()} "
                           "tilt-series using IMOD *tiltalign* command.")

        return methods

    # --------------------------- UTILS functions -----------------------------
    def getInputFiduSet(self, pointer=False):
        return self.inputSetOfLandmarkModels.get() if not pointer else self.inputSetOfLandmarkModels
    
    def getInputSet(self, pointer=False):
        inFiduSet = self.getInputFiduSet()
        return inFiduSet.getSetOfTiltSeries(pointer=pointer)

    def getRotationType(self):
        return {
            0: 0,
            1: -1,
            2: 3,
            3: 1
        }[self.rotationSolutionType.get()]

    def getMagnificationType(self):
        return {
            0: 0,
            1: 3,
            2: 1
        }[self.magnificationSolutionType.get()]

    def getTiltAngleType(self):
        return {
            0: 0,
            1: 5,
            2: 2
        }[self.tiltAngleSolutionType.get()]

    def getSkewType(self):
        return {
            0: 0,
            1: 3,
            2: 3
        }[self.distortionSolutionType.get()]

    def getStretchType(self):
        return {
            0: 0,
            2: 0,
            1: 3
        }[self.distortionSolutionType.get()]

    def getSurfaceToAnalyze(self):
        return 2 if self.twoSurfaces else 1

    @staticmethod
    def updateTsNonInterp(tsId, ts, tsOut, **kwargs):
        tsOut.setAlignment2D()

    @staticmethod
    def updateTiNonInterp(origIndex, index, tsId, ts, ti, tsOut, tiOut, alignmentMatrix=None, tltList=None, **kwargs):
        transform = Transform()
        newTransform = alignmentMatrix[:, :, index]
        newTransformArray = np.array(newTransform)

        if ti.hasTransform():
            previousTransform = ti.getTransform().getMatrix()
            previousTransformArray = np.array(previousTransform)
            outputTransformMatrix = np.matmul(newTransformArray, previousTransformArray)
            transform.setMatrix(outputTransformMatrix)
        else:
            transform.setMatrix(newTransformArray)

        tiOut.setTransform(transform)
        tiOut.setTiltAngle(float(tltList[index]))

    @staticmethod
    def updateTsInterp(tsId, ts, tsOut, **kwargs):
        tsOut.getAcquisition().setTiltAxisAngle(0.)
        tsOut.setAlignment(ALIGN_NONE)
        tsOut.setInterpolated(True)

    def updateTiInterp(self, origIndex, index, tsId, ts, ti, tsOut, tiOut, tltList=None, **kwargs):
        super().updateTi(origIndex, index, tsId, ts, ti, tsOut, tiOut, **kwargs)
        tiOut.setTiltAngle(float(tltList[index]))