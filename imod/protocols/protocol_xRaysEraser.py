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
import pyworkflow.utils.path as path
import pyworkflow.protocol.params as params
from tomo.protocols import ProtTomoBase
import imod.utils as utils


class ProtImodXraysEraser(EMProtocol, ProtTomoBase):
    """
    Erase X-rays, defects, and fiducial markers from aligned tilt-series based on the IMOD procedure.
    More info:
            https://bio3d.colorado.edu/imod/doc/man/ccderaser.html
    """

    _label = 'x-rays eraser'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

        # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfLandmarkModels',
                      params.PointerParam,
                      pointerClass='SetOfLandmarkModels',
                      important=True,
                      label='Input set of tilt-series',
                      help='Input set of landmark models containing the location of the gold beads through the series')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for lm in self.inputSetOfLandmarkModels.get():
            self._insertFunctionStep('convertInputStep', lm.getObjId())

    def convertInputStep(self, tsObjId):
        lm = self.inputSetOfLandmarkModels.get()[tsObjId]
        tsId = lm.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        path.makePath(extraPrefix)
        tmpPrefix = self._getTmpPath(tsId)
        utils.generateIMODFiducialTextFile(lm,
                                           os.path.join(extraPrefix, "%s_fid.txt" % lm.getTsId()))
