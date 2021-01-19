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

import pyworkflow.protocol.params as params
from pyworkflow.utils import path
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from imod import Plugin

import os


class ImodProtGoldBeadPicker3d(EMProtocol, ProtTomoBase):
    """
    3-dimensional gold bead picker using the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/findbeads3d.html
    """

    _label = 'Gold bead picker 3D'

# -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTomograms',
                      params.PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Input set of tilt-series',
                      help='Input set of tomograms from which gold beads will be picked.')

        form.addParam('beadDiameter',
                      params.IntParam,
                      label='Fiducial diameter (pixels)',
                      default='10',
                      help="Diameter of beads in pixels.")

        form.addParam('beadsColor',
                      params.EnumParam,
                      choices=['Dark', 'Light'],
                      label='Bead contrast',
                      default='0',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help="Contrast of the gold beads:\n"
                           "-Dark: beads are dark on light background.\n"
                           "-Light: beads are light on dark background.")

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTomograms.get():
            self._insertFunctionStep('pickGoldBeadsStep', ts.getObjId())
            #self._insertFunctionStep('createOutputStep', ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def pickGoldBeadsStep(self, tsObjId):
        tomo = self.inputSetOfTomograms.get()[tsObjId]
        location = tomo.getLocation()[1]
        fileName, fileExtension = os.path.splitext(location)


        # FEDE DEBUG ***
        path.copyFile(location, )

        extraPrefix = self._getExtraPath(os.path.basename(fileName))
        print(extraPrefix)
        #tmpPrefix = self._getTmpPath(os.path.basename(fileName))
        path.makePath(extraPrefix)
        #path.makePath(tmpPrefix)

        """ Run findbeads3d IMOD program """
        paramsFindbeads3d = {
            'inputFile': location,
            'outputFile': os.path.join(extraPrefix, "%s.mod" % os.path.basename(fileName)),
            'beadSize': self.beadDiameter.get(),
        }

        argsFindbeads3d = "-InputFile %(inputFile)s " \
                          "-OutputFile %(outputFile)s " \
                          "-BeadSize %(beadSize)d "

        Plugin.runImod(self, 'findbeads3d', argsFindbeads3d % paramsFindbeads3d)
