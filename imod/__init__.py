# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Centro Nacional de Biotecnologia, CSIC, Spain
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
import pwem
import pyworkflow as pw

from .constants import IMOD_HOME, ETOMO_CMD
from distutils.spawn import find_executable

_logo = ""
_references = ['Kremer1996', 'Mastronarde2017']


class Plugin(pwem.Plugin):
    _homeVar = IMOD_HOME
    _validationMsg = None

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(IMOD_HOME, 'imod-4.10.42/IMOD')

    @classmethod
    def getEnviron(cls):
        return None

    @classmethod
    def validateInstallation(cls):
        """
        Check if imod is in the path
        """

        if not cls._validationMsg:
            cls._validationMsg = [
                "imod's %s command not found in path, please install it." % ETOMO_CMD] if not find_executable(
                ETOMO_CMD) else []

        return cls._validationMsg

    @classmethod
    def getDependencies(cls):
        neededPrograms = ['libjpeg62', 'java', 'python']

        return neededPrograms

    @classmethod
    def defineBinaries(cls, env):
        version = '4.10.42'
        IMOD_INSTALLED = 'imod_%s_installed' % version

        if 'ubuntu' in os.getenv('DESKTOP_SESSION', 'unknown'):
            # Download .sh
            installationCmd = 'wget https://bio3d.colorado.edu/ftp/latestIMOD/RHEL6-64_CUDA8.0/' \
                              'imod_4.10.42_RHEL6-64_CUDA8.0.sh && '

            # Run .sh skipping copying startup scripts (avoid sudo permissions to write to /etc/profile.d)
            installationCmd += 'sh imod_4.10.42_RHEL6-64_CUDA8.0.sh -dir . -skip && '

            # Create installation finished flag file
            installationCmd += 'touch %s' % IMOD_INSTALLED

            env.addPackage('imod',
                           version=version,
                           tar='void.tgz',
                           neededProgs=cls.getDependencies(),
                           commands=[(installationCmd, IMOD_INSTALLED)],
                           default=True)

    @classmethod
    def runImod(cls, protocol, program, args, cwd=None):
        """ Run IMOD command from a given protocol. """
        fullProgram = '%s/bin/%s' % (cls.getVar(IMOD_HOME), program)
        protocol.runJob('sh', os.path.join(cls.getVar(IMOD_HOME), "IMOD-linux.sh"))
        protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
