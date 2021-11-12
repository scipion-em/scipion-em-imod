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


class ProtImodFiducialModel(ProtImodBase):
    """
    Construction of a fiducial model and alignment of tilt-series based on the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/autofidseed.html
        https://bio3d.colorado.edu/imod/doc/man/beadtrack.html
        https://bio3d.colorado.edu/imod/doc/man/model2point.html
    """

    _label = 'Generate fiducial model'
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
                      important=True,
                      help="Fiducials diameter to be tracked for alignment.")

        form.addParam('numberFiducial',
                      params.IntParam,
                      label='Number of fiducials',
                      default='25',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Number of fiducials to be tracked for alignment.")

        form.addParam('shiftsNearZeroFraction',
                      params.FloatParam,
                      label='Shifts near zero fraction',
                      default='0.2',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Fraction of the tracking box size above which to supply shifts near zero tilt to "
                           "Beadtrack. The dominant net shifts in the bead positions between views are found as "
                           "described above, and if one of the shifts is larger than this fraction of the "
                           "-BoxSizeXandY entry to Beadtrack, then the shifts are provided when running Beadtrack on "
                           "the initial seed models. Also, a command file will be written with modified parameters, "
                           "named as the root name of the input command file followed by '_adjusted' and its "
                           "extension. Enter 0 or a large value to disable this analysis.")

        groupGlobalVariables = form.addGroup('Filter variables',
                                             expertLevel=params.LEVEL_ADVANCED)

        groupGlobalVariables.addParam('refineSobelFilter',
                                      params.EnumParam,
                                      choices=['Yes', 'No'],
                                      default=1,
                                      label='Refine center with Sobel filter',
                                      expertLevel=params.LEVEL_ADVANCED,
                                      display=params.EnumParam.DISPLAY_HLIST,
                                      help='Use edge-detecting Sobel filter to refine the bead positions.')

        groupGlobalVariables.addParam('scalableSigmaForSobelFilter',
                                      params.FloatParam,
                                      default=0.5,
                                      condition='refineSobelFilter==0',
                                      label='Sobel sigma relative to bead size',
                                      expertLevel=params.LEVEL_ADVANCED,
                                      help='Sigma for gaussian kernel filtering of single beads before Sobel '
                                           'filtering, as fraction of bead diameter. The default sigma is 0.5 pixels '
                                           'regardless of bead size. A value of around 0.12 diameters is needed for '
                                           'higher noise (eg. cryo) data.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._failedTs = []

        for ts in self.inputSetOfTiltSeries.get():
            tsObjId = ts.getObjId()
            self._insertFunctionStep(self.convertInputStep, tsObjId, True, True)
            self._insertFunctionStep(self.generateTrackComStep, tsObjId)
            self._insertFunctionStep(self.generateFiducialSeedStep, tsObjId)
            self._insertFunctionStep(self.generateFiducialModelStep, tsObjId)
            self._insertFunctionStep(self.translateFiducialPointModelStep, tsObjId)
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
            except Exception as e:
                self._failedTs.append(tsId)

        return wrapper

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
            'rotationAngle': ts.getAcquisition().getTiltAxisAngle(),
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
            'targetNumberOfBeads': self.numberFiducial.get(),
            'shiftsNearZeroFraction': self.shiftsNearZeroFraction.get()
        }

        argsAutofidseed = "-TrackCommandFile %(trackCommandFile)s " \
                          "-MinSpacing %(minSpacing)f " \
                          "-PeakStorageFraction %(peakStorageFraction)f " \
                          "-TargetNumberOfBeads %(targetNumberOfBeads)d " \
                          "-ShiftsNearZeroFraction %(shiftsNearZeroFraction)f"

        if self.twoSurfaces.get() == 0:
            argsAutofidseed += " -TwoSurfaces"

        Plugin.runImod(self, 'autofidseed', (argsAutofidseed % paramsAutofidseed))

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

    def translateFiducialPointModelStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        firstItem = ts.getFirstItem()

        # Check that previous steps have been completed satisfactorily
        if os.path.exists(os.path.join(extraPrefix,
                                       firstItem.parseFileName(suffix="_gaps", extension=".fid"))):
            paramsGapModel2Point = {
                'inputFile': os.path.join(extraPrefix,
                                          firstItem.parseFileName(suffix="_gaps", extension=".fid")),
                'outputFile': os.path.join(extraPrefix,
                                           firstItem.parseFileName(suffix="_gaps_fid", extension=".txt"))
            }
            argsGapModel2Point = "-InputFile %(inputFile)s " \
                                 "-OutputFile %(outputFile)s"

            Plugin.runImod(self, 'model2point', argsGapModel2Point % paramsGapModel2Point)

    def computeOutputModelsStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        firstItem = ts.getFirstItem()

        # Create the output set of landmark models with gaps
        # Check that previous steps have been completed satisfactorily
        if os.path.exists(
                os.path.join(extraPrefix, firstItem.parseFileName(suffix="_gaps", extension=".fid"))):

            self.getOutputFiducialModelGaps()

            self.outputFiducialModelGaps.setSetOfTiltSeries(self.inputSetOfTiltSeries.get())

            landmarkModelGapsFilePath = os.path.join(
                extraPrefix,
                firstItem.parseFileName(suffix="_gaps", extension=".sfid")
            )

            fiducialModelGapPath = os.path.join(
                extraPrefix,
                firstItem.parseFileName(suffix="_gaps", extension=".fid")
            )

            fiducialModelGapTxtPath = os.path.join(
                extraPrefix,
                firstItem.parseFileName(suffix="_gaps_fid", extension=".txt")
            )

            fiducialGapList = utils.formatFiducialList(fiducialModelGapTxtPath)

            landmarkModelGaps = LandmarkModel(tsId=tsId,
                                              tiltSeriesPointer=ts,
                                              fileName=landmarkModelGapsFilePath,
                                              modelName=fiducialModelGapPath)

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

            self.outputFiducialModelGaps.append(landmarkModelGaps)
            self.outputFiducialModelGaps.update(landmarkModelGaps)
            self.outputFiducialModelGaps.write()

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
        if hasattr(self, "outputFiducialModelGaps"):
            self.outputFiducialModelGaps.setStreamState(Set.STREAM_CLOSED)
        if hasattr(self, "outputFailedSetOfTiltSeries"):
            self.outputFailedSetOfTiltSeries.setStreamState(Set.STREAM_CLOSED)

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
RotationAngle	%(rotationAngle).2f 
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

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputFiducialModelGaps'):
            summary.append("Input Tilt-Series: %d.\nFiducial models generated presenting gaps: %d."
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputFiducialModelGaps.getSize()))

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

        if hasattr(self, 'outputFailedSetOfTiltSeries'):
            methods.append("%d tilt-series have failed during the fiducial alignment protocol execution."
                           % (self.outputFailedSetOfTiltSeries.getSize()))

        if not methods:
            methods.append("Output classes not ready yet.")
        return methods
