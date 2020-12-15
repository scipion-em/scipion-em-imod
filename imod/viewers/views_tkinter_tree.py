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

import pyworkflow.viewer as pwviewer
from pyworkflow.gui.dialog import ListDialog
from imod import Plugin


class ImodDialog(ListDialog):
    # pasar un parametro modo y poner varios inits con ifs que lacen diferentes "doubleClickOn"
    def __init__(self, parent, **kwargs):
        self.provider = kwargs.get("provider", None)
        ListDialog.__init__(self, parent,
                            "Tilt-series List",
                            allowsEmptySelection=False,
                            itemDoubleClick=self.doubleClickOnTs,
                            **kwargs)

    def doubleClickOnTs(self, e=None):
        # Remove :mrc if present
        fn = " " + e.getFirstItem().getFileName().split(':')[0]
        pwviewer.CommandView("%s %s" % (Plugin.getImodCmd('3dmod'), fn))

