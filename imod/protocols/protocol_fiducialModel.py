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
import numpy as np

from tomo.objects import TiltSeriesDict, TiltSeries, Tomogram
from tomo.protocols import ProtTomoBase
from tomo.convert import writeTiStack


class ProtFiducialModel(pyem.EMProtocol, ProtTomoBase):
    """
    Construction of a fiducial model based on the IMOD procedure.

    More info:
        https://bio3d.colorado.edu/imod/doc/etomoTutorial.html
    """

    ######################################################################CLEAN
    commandline = "autofidseed -TrackCommandFile track.com -MinSpacing 0.85 -PeakStorageFraction 1.0 -TwoSurfaces " \
                  "-TargetNumberOfBeads 25"
    commandline = "autofidseed -StandardInput BB.preali -TrackCommandFile track.com -MinSpacing 0.85 " \
                  "-PeakStorageFraction 1.0 -TwoSurfaces -TargetNumberOfBeads 25"

    missingLines = "PrealignTransformFile	BBa.prexg --------------------------"

    ##########################################################################################################################

    _label = 'fiducial model'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries', params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-Series')

        form.addParam('twoSurfaces', params.EnumParam,
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
                      help="Number of fiducials to be tracked for alignment.")

        form.addParam('rotationAngle',
                      params.FloatParam,
                      label='Tilt rotation angle (deg)',
                      default='0.0',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Angle from the vertical to the tilt axis in raw images.")

        form.addParam('importTrackFile', params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Import track file',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Import a customized track file for the execution of the program.")

        group = form.addGroup('Track file',
                              condition='importTrackFile==0')

        group.addParam('trackFilePath',
                       params.PathParam,
                       label='File path',
                       expertLevel=params.LEVEL_ADVANCED,
                       help="Customized file path location.")

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('generateTrackComStep')
        self._insertFunctionStep('generateFiducialSeedStep')
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
                              "-TargetNumberOfBeads %(targetNumberOfBeads)d "

            if self.twoSurfaces.get() == 0:
                argsAutofidseed += "-TwoSurfaces "

            self.runJob('autofidseed', argsAutofidseed % paramsAutofidseed, cwd=workingFolder)

    def _createOutputStep(self):
        pass

    # --------------------------- UTILS functions ----------------------------
    def translateTrackCom(self, tsId, paramsDict):
        trackFileName = tsId + "/" + tsId + "_track.com"
        trackFilePath = os.path.join(self._getExtraPath(), trackFileName)

        if self.importTrackFile == 0:
            with open(self.trackFilePath.get(), 'r') as f:
                template = f.read()
            with open(trackFilePath, 'w') as f:
                f.write(template)
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
# To exclude views, add a line "SkipViews view_list" with the list of views
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
