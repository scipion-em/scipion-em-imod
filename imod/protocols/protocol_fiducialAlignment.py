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
from pwem.objects import Transform
from pwem.emlib.image import ImageHandler
from pyworkflow.object import Set
from tomo.objects import (LandmarkModel, SetOfTiltSeries, TiltImage,
                          TiltSeries, TiltSeriesCoordinate)

from .. import Plugin, utils
from .protocol_base import ProtImodBase, TLT_EXT, XF_EXT, FID_EXT, TXT_EXT, XYZ_EXT, MOD_EXT, SFID_EXT


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
    _possibleOutputs = {"outputSetOfTiltSeries": SetOfTiltSeries}

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfLandmarkModels',
                      params.PointerParam,
                      pointerClass='SetOfLandmarkModels',
                      important=True,
                      label='Fiducial model')

        # TODO: Allow for a different set of tilt-series input source than the one from the landmark model. This is not
        # TODO: possible due to a change of data type when applying the transformation with scipion applyTransform
        # TODO: method due to in a change in the output datatype (always float) which triggers the following error in
        # TODO: the imod tiltalign program:
        # TODO: ERROR: TILTALIGN - TWO POINTS (#    2 AND    3) ON VIEW    2 IN CONTOUR    1 OF OBJECT   1

        # form.addParam('setOfTiltSeriesSource',
        #               params.EnumParam,
        #               choices=['Yes', 'No'],
        #               default=0,
        #               label='Use same set of tilt-series form model',
        #               display=params.EnumParam.DISPLAY_HLIST,
        #               help="By default the set of tilt-series to be algined is the same from which the fiducial models"
        #                    "have been obtained. If the user wants to sepecify the set of tilt series to be aligned "
        #                    "then select No.")
        #
        # form.addParam('inputSetOfTiltSeries',
        #               params.PointerParam,
        #               pointerClass='SetOfTiltSeries',
        #               condition='setOfTiltSeriesSource==1',
        #               label='Input set of tilt-series.')

        form.addParam('twoSurfaces',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Assume beads on two surfaces?',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help="Track fiducials differentiating in which side of the sample are located.\n"
                           "IMPORTANT: It is highly recommended to match the option selected in the "
                           "generation of the fiducial models. In case they  do not match, it is not "
                           "intended to fail but could be missing the whole potential of the algorithm. "
                           "In case the algorithm used for the calculation of the fiducial models does "
                           "not consider this option it is algo recommended to set this "
                           "option to 'No'.")

        form.addParam('computeAlignment',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Generate interpolated tilt-series?',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Generate and save the interpolated tilt-series applying the obtained transformation '
                           'matrices.\n'
                           'By default, the output of this protocol will be a tilseries that will have associated'
                           'the alignment information as a transformation matrix. When this option is set as Yes, '
                           'then a second output, called interpolated tilt series, is generated. The interpolated tilt '
                           'series should be used for visualization purpose but not for image processing')

        groupInterpolation = form.addGroup('Interpolated tilt-series',
                                           condition='computeAlignment==0')

        groupInterpolation.addParam('binning',
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
                      choices=['No rotation', 'One rotation',
                               'Group rotations', 'Solve for all rotations'],
                      default=3,
                      label='Rotation solution type',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Type of rotation solution: See rotOption in tiltalign IMOD command \n'
                           '* No rotation: The in-plane rotation will not be estimated\n'
                           '* One rotation: To solve for a single rotation variable \n'
                           '* Group rotations: Group views to solve for fewer rotations variables. Automapping of '
                           'rotation variables linearly changing values\n'
                           '* Solve for all rotations: for each view having an independent rotation\n')

        form.addParam('groupRotationSize',
                      params.IntParam,
                      default=5,
                      condition='rotationSolutionType==2',
                      label='Group size',
                      help='Default group size when automapping rotation variables')

        form.addParam('magnificationSolutionType',
                      params.EnumParam,
                      choices=['Fixed magnification at 1.0',
                               'Group magnifications',
                               'Solve for all magnifications'],
                      default=1,
                      label='Magnification solution type',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Type of magnification solution: See MagOption in tiltaling IMOD command\n'
                           '* Fixed magnification: Do not solve magnification. This fixes all magnifications at 1.0.\n'
                           '* Group magnifications: Group views to solve for fewer magnifications variables. '
                           'Automapping of variables (linearly changing values)  \n'
                           '* Solve for all magnifications: to vary all magnifications  of each view independently\n')

        form.addParam('groupMagnificationSize',
                      params.IntParam,
                      default=4,
                      condition='magnificationSolutionType==1',
                      label='Group size',
                      help='Group size when automapping magnification variables')

        form.addParam('tiltAngleSolutionType',
                      params.EnumParam,
                      choices=['Fixed tilt angles', 'Group tilt angles',
                               'Solve for all except minimum tilt'],
                      default=1,
                      label='Tilt angle solution type',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Type of tilt angle solution: See TiltOption in tiltalign IMOD command\n'
                           ' * Fixed tilt angles: To fix all tilt angles at their initial (input) values \n'
                           ' * Group tilt angles: To automap groups of tilt angles (linearly changing values) \n'
                           ' * Solve for all except minimum tilt:to solve for all tilt angles except for the view '
                           'at minimum tilt \n')

        form.addParam('groupTiltAngleSize',
                      params.IntParam,
                      default=5,
                      condition='tiltAngleSolutionType==1',
                      label='Group size',
                      help='Average default group size when automapping tilt variables')

        form.addParam('distortionSolutionType',
                      params.EnumParam,
                      choices=['Disabled', 'Full solution', 'Skew only'],
                      default=0,
                      label='Distortion solution type',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Type of skew solution:'
                           '* 0 to fix all skew angles at 0.0 \n'
                           '* 1 to vary all skew angles independently\n '
                           '* 2 to specify a mapping of skew variables, or\n '
                           '* 3 or 4 for automapping of variables (3 for linearly changing values or 4 for values all '
                           'the same within a group)..')

        form.addParam('xStretchGroupSize',
                      params.IntParam,
                      default=7,
                      condition='distortionSolutionType==1',
                      label='X stretch group size',
                      help='Basic grouping size for X stretch')

        form.addParam('skewGroupSize',
                      params.IntParam,
                      default=11,
                      condition='tiltAngleSolutionType==1 or tiltAngleSolutionType==2',
                      label='Skew group size',
                      help='Size of the skew group')

        form.addSection('Erase gold beads')

        form.addParam('eraseGoldBeads',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Erase gold beads',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Remove the gold beads detected during fiducial '
                           'alignment with *ccderaser* program. This option '
                           'will generate an interpolated tilt series with '
                           'the gold beads erased and interpolated with '
                           'the calculated transformation matrices form '
                           'the alignment.')

        groupEraseGoldBeads = form.addGroup('Gold bead eraser',
                                            condition='eraseGoldBeads==0')

        groupEraseGoldBeads.addParam('betterRadius',  # actually diameter
                                     params.IntParam,
                                     default=18,
                                     label='Bead diameter (px)',
                                     help="For circle objects, this entry "
                                          "specifies a radius to use for points "
                                          "without an individual point size "
                                          "instead of the object's default sphere "
                                          "radius. This entry is floating point "
                                          "and can be used to overcome the "
                                          "limitations of having an integer "
                                          "default sphere radius. If there are "
                                          "multiple circle objects, enter one "
                                          "value to apply to all objects or a "
                                          "value for each object.")

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        for lm in self.inputSetOfLandmarkModels.get():
            lmTsId = lm.getTsId()
            self.fiducialDiameterPixel = lm.getSize()
            self._insertFunctionStep(self.convertInputStep, lmTsId)
            self._insertFunctionStep(self.computeFiducialAlignmentStep, lmTsId)
            self._insertFunctionStep(self.translateFiducialPointModelStep, lmTsId)
            self._insertFunctionStep(self.computeOutputStackStep, lmTsId)

            if self.computeAlignment.get() == 0 or self.eraseGoldBeads.get() == 0:
                self._insertFunctionStep(self.computeOutputInterpolatedStackStep,
                                         lmTsId)

            if self.eraseGoldBeads.get() == 0:
                self._insertFunctionStep(self.eraseGoldBeadsStep, lmTsId)

            self._insertFunctionStep(self.computeOutputModelsStep, lmTsId)
            self._insertFunctionStep(self.createOutputFailedStep, lmTsId)

        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.inputSetOfTiltSeries = self.inputSetOfLandmarkModels.get().getSetOfTiltSeries(pointer=True)

        tsIds = self.inputSetOfLandmarkModels.get().aggregate(["COUNT"], "_tsId", ["_tsId"])
        tsIds = set([d['_tsId'] for d in tsIds])

        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in
                       self._getSetOfInputTS() if
                       ts.getTsId() in tsIds}

        self._failedTs = []

    @ProtImodBase.tryExceptDecorator
    def computeFiducialAlignmentStep(self, tsId):
        ts = self.tsDict[tsId]
        lm = self.inputSetOfLandmarkModels.get().getLandmarkModelFromTsId(tsId=tsId)

        paramsTiltAlign = {
            'modelFile': lm.getModelName(),
            'imageFile': self.getTmpOutFile(tsId),
            'imagesAreBinned': 1,
            'unbinnedPixelSize': ts.getSamplingRate() / 10,
            'outputModelFile': self.getExtraOutFile(tsId, suffix="fidxyz", ext=MOD_EXT),
            'outputResidualFile': self.getExtraOutFile(tsId, suffix="resid", ext=TXT_EXT),
            'outputFidXYZFile': self.getExtraOutFile(tsId, suffix="fid", ext=XYZ_EXT),
            'outputTiltFile': self.getExtraOutFile(tsId, suffix="interpolated", ext=TLT_EXT),
            'outputXAxisTiltFile': self.getExtraOutFile(tsId, ext="xtilt"),
            'outputTransformFile': self.getExtraOutFile(tsId, suffix="fid", ext=XF_EXT),
            'outputFilledInModel': self.getExtraOutFile(tsId, suffix="noGaps", ext=FID_EXT),
            'rotationAngle': ts.getAcquisition().getTiltAxisAngle(),
            'tiltFile': self.getExtraOutFile(tsId, ext=TLT_EXT),
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
            'targetPatchSizeXandY': '700,700',
            'minSizeOrOverlapXandY': '0.5,0.5',
            'minFidsTotalAndEachSurface': '8,3',
            'fixXYZCoordinates': 0,
            'outputTiltAlignFileText': self._getExtraPath("align.log"),
        }

        argsTiltAlign = "-ModelFile %(modelFile)s " \
                        "-ImageFile %(imageFile)s " \
                        "-ImagesAreBinned %(imagesAreBinned)d " \
                        "-UnbinnedPixelSize %(unbinnedPixelSize)f " \
                        "-OutputModelFile %(outputModelFile)s " \
                        "-OutputResidualFile %(outputResidualFile)s " \
                        "-OutputFidXYZFile %(outputFidXYZFile)s " \
                        "-OutputTiltFile %(outputTiltFile)s " \
                        "-OutputXAxisTiltFile %(outputXAxisTiltFile)s " \
                        "-OutputTransformFile %(outputTransformFile)s " \
                        "-OutputFilledInModel %(outputFilledInModel)s " \
                        "-RotationAngle %(rotationAngle).2f " \
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
                        "-TargetPatchSizeXandY %(targetPatchSizeXandY)s " \
                        "-MinSizeOrOverlapXandY %(minSizeOrOverlapXandY)s " \
                        "-MinFidsTotalAndEachSurface %(minFidsTotalAndEachSurface)s " \
                        "-FixXYZCoordinates %(fixXYZCoordinates)d " \
                        "-RobustFitting "

        # Excluded views
        excludedViews = ts.getExcludedViewsIndex(caster=str)
        if len(excludedViews):
            argsTiltAlign += f"-ExcludeList {','.join(excludedViews)} "

        argsTiltAlign += "2>&1 | tee %(outputTiltAlignFileText)s "

        Plugin.runImod(self, 'tiltalign', argsTiltAlign % paramsTiltAlign)
        Plugin.runImod(self, 'alignlog', '-s > taSolution.log', cwd=self._getExtraPath())

    @ProtImodBase.tryExceptDecorator
    def translateFiducialPointModelStep(self, tsId):
        # Check that previous steps have been completed satisfactorily
        noGapsFid = self.getExtraOutFile(tsId, suffix="noGaps", ext=FID_EXT)
        if os.path.exists(noGapsFid):
            paramsNoGapModel2Point = {
                'inputFile': noGapsFid,
                'outputFile': self.getExtraOutFile(tsId, suffix="noGaps_fid", ext=TXT_EXT)
            }
            argsNoGapModel2Point = "-InputFile %(inputFile)s " \
                                   "-OutputFile %(outputFile)s"

            Plugin.runImod(self, 'model2point', argsNoGapModel2Point % paramsNoGapModel2Point)

    @ProtImodBase.tryExceptDecorator
    def computeOutputStackStep(self, tsId):
        ts = self.tsDict[tsId]

        # Check that previous steps have been completed satisfactorily
        transformationMatricesFilePath = self.getExtraOutFile(tsId, suffix="fid", ext=XF_EXT)
        if os.path.exists(transformationMatricesFilePath) and os.stat(transformationMatricesFilePath).st_size != 0:
            tltFilePath = self.getExtraOutFile(tsId, suffix="interpolated", ext=TLT_EXT)
            tltList = utils.formatAngleList(tltFilePath)
            newTransformationMatricesList = utils.formatTransformationMatrix(transformationMatricesFilePath)
            output = self.getOutputSetOfTiltSeries(self._getSetOfInputTS())
            newTs = TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            output.append(newTs)

            for index, tiltImage in enumerate(ts):
                newTi = TiltImage()
                newTi.copyInfo(tiltImage, copyId=True, copyTM=False)
                newTi.setLocation(tiltImage.getLocation())
                newTi.setTiltAngle(float(tltList[index]))
                newTi.setAcquisition(tiltImage.getAcquisition())

                if tiltImage.hasTransform():
                    transform = Transform()
                    previousTransform = tiltImage.getTransform().getMatrix()
                    newTransform = newTransformationMatricesList[:, :, index]
                    previousTransformArray = np.array(previousTransform)
                    newTransformArray = np.array(newTransform)
                    outputTransformMatrix = np.matmul(newTransformArray, previousTransformArray)
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

            output.update(newTs)
            output.write()

            self._store()
        else:
            raise FileNotFoundError(
                "Error (computeOutputStackStep): \n Imod output file "
                "%s does not exist or it is empty" % transformationMatricesFilePath)

    @ProtImodBase.tryExceptDecorator
    def computeOutputInterpolatedStackStep(self, tsId):
        tsIn = self.tsDict[tsId]
        tsId = tsIn.getTsId()
        firstItem = tsIn.getFirstItem()

        # Check that previous steps have been completed satisfactorily
        tmpFileName = self.getExtraOutFile(tsId, suffix="fid", ext=XF_EXT)
        if os.path.exists(tmpFileName) and os.stat(tmpFileName).st_size != 0:
            output = self.getOutputInterpolatedTS(self._getSetOfInputTS())

            paramsAlignment = {
                'input': self.getTmpOutFile(tsId),
                'output': self.getExtraOutFile(tsId),
                'xform': tmpFileName,
                'bin': self.binning.get(),
                'imagebinned': 1.0}

            argsAlignment = "-input %(input)s " \
                            "-output %(output)s " \
                            "-xform %(xform)s " \
                            "-bin %(bin)d " \
                            "-antialias -1 " \
                            "-imagebinned %(imagebinned)s " \
                            "-taper 1,1 "

            rotationAngle = tsIn.getAcquisition().getTiltAxisAngle()

            # Check if rotation angle is greater than 45º. If so, swap x
            # and y dimensions to adapt output image sizes to
            # the final sample disposition.
            if 45 < abs(rotationAngle) < 135:
                paramsAlignment.update({
                    'size': "%d,%d" %
                            (firstItem.getYDim() // self.binning.get(),
                             firstItem.getXDim() // self.binning.get())
                })

                argsAlignment += " -size %(size)s "

            Plugin.runImod(self, 'newstack', argsAlignment % paramsAlignment)

            newTs = TiltSeries(tsId=tsId)
            newTs.copyInfo(tsIn)
            newTs.setInterpolated(True)
            output.append(newTs)

            tltFilePath = self.getExtraOutFile(tsId, suffix="interpolated", ext=TLT_EXT)
            tltList = utils.formatAngleList(tltFilePath)

            if self.binning > 1:
                newTs.setSamplingRate(tsIn.getSamplingRate() * self.binning.get())

            for index, tiltImage in enumerate(tsIn):
                newTi = TiltImage()
                newTi.copyInfo(tiltImage, copyId=True, copyTM=False)
                newTi.setAcquisition(tiltImage.getAcquisition())
                newTi.setLocation(index + 1, self.getExtraOutFile(tsId))
                newTi.setTiltAngle(float(tltList[index]))
                if self.binning > 1:
                    newTi.setSamplingRate(tiltImage.getSamplingRate() * self.binning.get())
                newTs.append(newTi)

            dims = self._getOutputDim(newTi.getFileName())
            newTs.setDim(dims)
            newTs.write(properties=False)

            output.update(newTs)
            output.write()
            self._store()
        else:
            raise FileNotFoundError(
                "Error (computeOutputInterpolatedStackStep): \n "
                "Imod output file %s does not exist or it is empty" % tmpFileName)

    @ProtImodBase.tryExceptDecorator
    def eraseGoldBeadsStep(self, tsId):
        ts = self.tsDict[tsId]
        tsId = ts.getTsId()

        # Erase gold beads on aligned stack
        paramsCcderaser = {
            'inputFile': self.getTmpOutFile(tsId),
            'outputFile': self.getExtraOutFile(tsId),
            'modelFile': self.getExtraOutFile(tsId, suffix="noGaps", ext=FID_EXT),
            'betterRadius': self.betterRadius.get() / 2,
            'polynomialOrder': 0,
            'circleObjects': "/"
        }

        argsCcderaser = "-InputFile %(inputFile)s " \
                        "-OutputFile %(outputFile)s " \
                        "-ModelFile %(modelFile)s " \
                        "-BetterRadius %(betterRadius)f " \
                        "-PolynomialOrder %(polynomialOrder)d " \
                        "-CircleObjects %(circleObjects)s " \
                        "-MergePatches 1 " \
                        "-ExcludeAdjacent " \
                        "-SkipTurnedOffPoints 1 " \
                        "-ExpandCircleIterations 3 "

        Plugin.runImod(self, 'ccderaser', argsCcderaser % paramsCcderaser)

    @ProtImodBase.tryExceptDecorator
    def computeOutputModelsStep(self, tsId):
        ts = self.tsDict[tsId]
        tsId = ts.getTsId()

        # Create the output set of landmark models with no gaps
        fiducialNoGapFilePath = self.getExtraOutFile(tsId, suffix="noGaps_fid", ext=TXT_EXT)
        if os.path.exists(fiducialNoGapFilePath):
            output = self.getOutputFiducialModelNoGaps()
            output.setSetOfTiltSeries(self._getSetOfInputTS())
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
            output.write()

        # Create the output set of 3D coordinates
        coordFilePath = self.getExtraOutFile(tsId, suffix="fid", ext=XYZ_EXT)

        if os.path.exists(coordFilePath):

            output = self.getOutputSetOfTiltSeriesCoordinates(self._getSetOfInputTS())

            coordList, xDim, yDim = utils.format3DCoordinatesList(coordFilePath)

            for element in coordList:
                newCoord3D = TiltSeriesCoordinate()
                newCoord3D.setTsId(ts.getTsId())
                newCoord3D.setPosition(element[0] - (xDim / 2),
                                       element[1] - (yDim / 2),
                                       element[2],
                                       sampling_rate=ts.getSamplingRate())

                output.append(newCoord3D)
            output.write()
            self._store()

    def createOutputStep(self):
        if self.TiltSeries:
            self.TiltSeries.setStreamState(Set.STREAM_CLOSED)
        if self.InterpolatedTiltSeries:
            self.InterpolatedTiltSeries.setStreamState(Set.STREAM_CLOSED)
        if self.FiducialModelNoGaps:
            self.FiducialModelNoGaps.setStreamState(Set.STREAM_CLOSED)
        if self.TiltSeriesCoordinates:
            self.TiltSeriesCoordinates.setStreamState(Set.STREAM_CLOSED)
        if self.FailedTiltSeries:
            self.FailedTiltSeries.setStreamState(Set.STREAM_CLOSED)

        self._store()

    def createOutputFailedStep(self, tsId):
        ts = self.tsDict[tsId]
        super().createOutputFailedSet(ts)

    # --------------------------- UTILS functions -----------------------------
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

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if self.FiducialModelNoGaps:
            summary.append("Fiducial models generated with no gaps: %d"
                           % (self.FiducialModelNoGaps.getSize()))

        if self.TiltSeries:
            summary.append("Transformation matrices updated from the "
                           "input tilt-series: %d"
                           % (self.TiltSeries.getSize()))

        if self.InterpolatedTiltSeries:
            summary.append("Interpolated tilt-series calculated: %d"
                           % (self.InterpolatedTiltSeries.getSize()))

        if self.TiltSeriesCoordinates:
            summary.append("Fiducial 3D coordinates calculated: %d"
                           % (self.TiltSeriesCoordinates.getSize()))

        if self.FailedTiltSeries:
            summary.append("Failed tilt-series: %d"
                           % (self.FailedTiltSeries.getSize()))

        if not summary:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []

        if self.TiltSeriesCoordinates:
            methods.append("Solved fiducials alignment for %d "
                           "tilt-series using IMOD *tiltalign* command."
                           % (self.FiducialModelNoGaps.getSize()))

        return methods
