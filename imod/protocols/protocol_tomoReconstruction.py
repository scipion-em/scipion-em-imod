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
from pwem.protocols import EMProtocol
import pyworkflow as pw
import pyworkflow.protocol.params as params
from tomo.protocols import ProtTomoBase
from tomo.convert import writeTiStack


class ProtTomoReconstruction(EMProtocol, ProtTomoBase):
    """
    Tomogram reconstruction procedure based on the IMOD procedure.

    More info:
        https://bio3d.colorado.edu/imod/doc/etomoTutorial.html
    """

    _label = 'tomogram reconstruction'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries', params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-Series')

        form.addParam('tomoThickness', params.EnumParam,
                      default=100,
                      label='Tomogram thickness', important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Size in pixels of the tomogram in the z axis (beam direction).')

        form.addParam('binning', params.FloatParam,
                       default=1.0,
                       label='Binning',
                       help='Binning to be applied to the interpolated tilt-series. '
                            'Must be a integer bigger than 1')

        form.addParam('rotationAngle',
                      params.FloatParam,
                      label='Tilt rotation angle (deg)',
                      default='0.0',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Angle from the vertical to the tilt axis in raw images.")

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('computeReconstructionStep')

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

    def computeReconstructionStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)
            # self.runJob('cp', '/home/fede/Downloads/Etomo_tutorialData/BBa.ali /home/fede/ScipionUserData/projects/Tomo_IMOD/Runs/019116_ProtTomoReconstruction/extra/BBa/BBa.ali')
            self.runJob('tilt', '-InputProjections %s.st -OutputFile %s.mrc -TILTFILE %s.rawtlt -THICKNESS 100' % (
            tsId, tsId, tsId), cwd=workingFolder)
