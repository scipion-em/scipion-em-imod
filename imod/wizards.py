# *****************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Centro Nacional de Biotecnologia, CSIC, Spain
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

from pyworkflow.wizard import Wizard

from imod.protocols import ProtImodXcorrPrealignment


class ProtImodXcorrPrealignmentWizard(Wizard):
    _targets = [(ProtImodXcorrPrealignment, ['tiltAxisAngle'])]

    @staticmethod
    def _getTiltAxisAngle(protocol):
        tiltAxisAngle = None
        if protocol.inputSetOfTiltSeries.hasValue():
            tiltAxisAngle = protocol.inputSetOfTiltSeries.get().getAcquisition().getTiltAxisAngle()

        return tiltAxisAngle

    def show(self, form, *args):
        angle = self._getTiltAxisAngle(form.protocol)
        if angle is not None:
            form.setVar('tiltAxisAngle', angle)
