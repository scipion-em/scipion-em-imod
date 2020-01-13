# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import pwem


from .constants import IMOD_HOME, ETOMO_CMD
from distutils.spawn import find_executable

_logo = ""
_references = ['Kremer1996', 'Mastronarde2017']



class Plugin(pwem.Plugin):

    _validationMsg = None

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(IMOD_HOME, 'imod-4.9')

    @classmethod
    def getEnviron(cls):
        return None

    @classmethod
    def validateInstallation(cls):
        """
        Check if imod is in the path
        """

        if not cls._validationMsg:
            cls._validationMsg = ["imod's %s command not found in path, please install it." % ETOMO_CMD] if not find_executable(ETOMO_CMD) else []

        return cls._validationMsg


