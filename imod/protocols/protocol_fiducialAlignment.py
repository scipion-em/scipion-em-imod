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
import tomo.objects as tomoObj
from pyworkflow.object import Set
from tomo.objects import LandmarkModel
import tomo.constants as constants
from pwem.emlib.image import ImageHandler
from imod import Plugin
from imod.protocols.protocol_base import ProtImodBase

class ProtImodFiducialAlignment(ProtImodBase):
    """
    Construction of a fiducial model and alignment of tilt-series based on the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/tiltalign.html
        https://bio3d.colorado.edu/imod/doc/man/model2point.html
        https://bio3d.colorado.edu/imod/doc/man/imodtrans.html
        https://bio3d.colorado.edu/imod/doc/man/newstack.html
        https://bio3d.colorado.edu/imod/doc/man/ccderaser.html
    """

    _label = 'Fiducial alignment'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfLandmarkModels',
                      params.PointerParam,
                      pointerClass='SetOfLandmarkModels',
                      important=True,
                      label='Input set of fiducial models.')

        # TODO: Allow for a different set of tilt-series input source than the one from the landmark model. This is not
        # TODO: possible due to a change of data type when applying the transformation with scipion applyTransform
        # TODO: mthod due to in a change in the output datatype (always float) which triggers the following error in
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
                      default=0,
                      label='Find on two surfaces',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help="Track fiducials differentiating in which side of the sample are located.\n"
                           "IMPORTANT: It is highly recmended to match the option selected in the generation of the "
                           "fiducial models. In case they do not match, it is not intended to fail but could be "
                           "missing the whole potential of the algorithm. In case the algorithm used fot he calculation"
                           "of the fiducial models does not consider this option it is algo recomended to set this "
                           "option to 'No'.")

        form.addParam('rotationAngle',
                      params.FloatParam,
                      label='Tilt rotation angle (deg)',
                      default='0.0',
                      important=True,
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

        self.inputSetOfTiltSeries = self.inputSetOfLandmarkModels.get().getSetOfTiltSeries(pointer=True)

        for ts in self.inputSetOfTiltSeries.get():
            tsObjId = ts.getObjId()
            self._insertFunctionStep(self.convertInputStep, tsObjId, True, True)
            self._insertFunctionStep(self.computeFiducialAlignmentStep, tsObjId)
            self._insertFunctionStep(self.translateFiducialPointModelStep, tsObjId)
            self._insertFunctionStep(self.computeOutputStackStep, tsObjId)

            if self.computeAlignment.get() == 0 or self.eraseGoldBeads.get() == 0:
                self._insertFunctionStep(self.computeOutputInterpolatedStackStep, tsObjId)

            if self.eraseGoldBeads.get() == 0:
                self._insertFunctionStep(self.eraseGoldBeadsStep, tsObjId)

            self._insertFunctionStep(self.computeOutputModelsStep, tsObjId)
            self._insertFunctionStep(self.createOutputFailedSet, tsObjId)

        self._insertFunctionStep(self.createOutputStep)

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

    @tryExceptDecorator
    def computeFiducialAlignmentStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        lm = self.inputSetOfLandmarkModels.get().getLandmarkModelFromTsId(tsId=tsId)

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        firstItem = ts.getFirstItem()

        paramsTiltAlign = {
            'modelFile': lm.getModelName(),
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
                        "2>&1 | tee %(outputTiltAlignFileText)s "

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

            self.getOutputSetOfTiltSeries(self.inputSetOfTiltSeries.get())
            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            self.outputSetOfTiltSeries.append(newTs)

            for index, tiltImage in enumerate(ts):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)
                newTi.setLocation(tiltImage.getLocation())
                newTi.setTiltAngle(float(tltList[index]))

                if tiltImage.hasTransform():
                    transform = Transform()
                    previousTransform = tiltImage.getTransform().getMatrix()
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

            self.outputSetOfTiltSeries.update(newTs)
            self.outputSetOfTiltSeries.write()

            self._store()

    def computeOutputInterpolatedStackStep(self, tsObjId):
        tsIn = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = tsIn.getTsId()
        outputTsIdList = [ts.getTsId() for ts in self.inputSetOfTiltSeries.get()]

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        firstItem = tsIn.getFirstItem()

        # Check that previous steps have been completed satisfactorily
        if os.path.exists(os.path.join(extraPrefix, firstItem.parseFileName(suffix="_fid", extension=".xf"))):
            self.getOutputInterpolatedSetOfTiltSeries(self.inputSetOfTiltSeries.get())

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

            rotationAngleAvg = utils.calculateRotationAngleFromTM(
                self.outputSetOfTiltSeries[outputTsIdList.index(tsId) + 1])

            # Check if rotation angle is greater than 45º. If so, swap x and y dimensions to adapt output image sizes to
            # the final sample disposition.
            if rotationAngleAvg > 45 or rotationAngleAvg < -45:
                paramsAlignment.update({
                    'size': "%d,%d" % (firstItem.getYDim(), firstItem.getXDim())
                })

                argsAlignment += " -size %(size)s "

            Plugin.runImod(self, 'newstack', argsAlignment % paramsAlignment)

            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(tsIn)
            self.outputInterpolatedSetOfTiltSeries.append(newTs)

            tltFilePath = os.path.join(
                extraPrefix,
                firstItem.parseFileName(suffix="_interpolated", extension=".tlt")
            )

            tltList = utils.formatAngleList(tltFilePath)

            if self.binning > 1:
                newTs.setSamplingRate(tsIn.getSamplingRate() * int(self.binning.get()))

            for index, tiltImage in enumerate(tsIn):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)
                newTi.setAcquisition(tiltImage.getAcquisition())
                newTi.setLocation(index + 1, os.path.join(extraPrefix, tiltImage.parseFileName(suffix="_interpolated")))
                newTi.setTiltAngle(float(tltList[index]))
                if self.binning > 1:
                    newTi.setSamplingRate(tiltImage.getSamplingRate() * int(self.binning.get()))
                newTs.append(newTi)

            ih = ImageHandler()
            x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
            newTs.setDim((x, y, z))
            newTs.write(properties=False)

            self.outputInterpolatedSetOfTiltSeries.update(newTs)
            self.outputInterpolatedSetOfTiltSeries.updateDim()
            self.outputInterpolatedSetOfTiltSeries.write()
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

            outputSetOfCoordinates3D = \
                self.getOutputSetOfCoordinates3Ds(self.inputSetOfTiltSeries.get(), self.outputSetOfTiltSeries)

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
            self.getOutputFailedSetOfTiltSeries(self.inputSetOfTiltSeries.get())

            ts = self.inputSetOfTiltSeries.get()[tsObjId]
            tsId = ts.getTsId()

            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            self.outputFailedSetOfTiltSeries.append(newTs)

            for index, tiltImage in enumerate(ts):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)
                newTi.setAcquisition(tiltImage.getAcquisition())
                newTi.setLocation(tiltImage.getLocation())
                if self.binning > 1:
                    newTi.setSamplingRate(tiltImage.getSamplingRate() * int(self.binning.get()))
                newTs.append(newTi)

            ih = ImageHandler()
            x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
            newTs.setDim((x, y, z))
            newTs.write(properties=False)

            self.outputFailedSetOfTiltSeries.update(newTs)
            self.outputFailedSetOfTiltSeries.updateDim()
            self.outputFailedSetOfTiltSeries.write()
            self._store()

    def createOutputStep(self):
        if hasattr(self, "outputSetOfTiltSeries"):
            self.outputSetOfTiltSeries.setStreamState(Set.STREAM_CLOSED)
        if hasattr(self, "outputInterpolatedSetOfTiltSeries"):
            self.outputInterpolatedSetOfTiltSeries.setStreamState(Set.STREAM_CLOSED)
        if hasattr(self, "outputFiducialModelNoGaps"):
            self.getOutputFiducialModelNoGaps().setStreamState(Set.STREAM_CLOSED)
        if hasattr(self, "outputSetOfCoordinates3D"):
            self.getOutputSetOfCoordinates3Ds().setStreamState(Set.STREAM_CLOSED)
        if hasattr(self, "outputFailedSetOfTiltSeries"):
            self.outputFailedSetOfTiltSeries.setStreamState(Set.STREAM_CLOSED)

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

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []

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
