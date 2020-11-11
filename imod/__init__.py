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

from .constants import IMOD_HOME, ETOMO_CMD, DEFAULT_VERSION
from distutils.spawn import find_executable
from pyworkflow.gui.project.utils import OS

__version__ = '3.0.5'
_logo = ""
_references = ['Kremer1996', 'Mastronarde2017']


class Plugin(pwem.Plugin):
    _homeVar = IMOD_HOME
    _validationMsg = None

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(IMOD_HOME, cls._getIMODFolder(DEFAULT_VERSION))

    @classmethod
    def _getEMFolder(cls, version, *paths):
        return os.path.join("imod-%s" % version, *paths)

    @classmethod
    def _getIMODFolder(cls, version, *paths):
        return  os.path.join(cls._getEMFolder(version, "IMOD"), *paths)

    @classmethod
    def _getProgram(cls, program):
        """ Returns the same program  if config missing
        or the path to the program based on the config file."""
        # Compose path based on config
        progFromConfig = cls.getHome("bin", program)

        # Check if IMOD from config exists
        if os.path.exists(progFromConfig):
            return progFromConfig
        else:
            return program

    @classmethod
    def getEnviron(cls):
        env=pwem.pwutils.Environ(os.environ)
        if 'IMOD_DIR' in env:
            del env['IMOD_DIR']
        if 'IMOD_PATH' in env:
            del env['IMOD_PATH']
        return env

    @classmethod
    def validateInstallation(cls):
        """
        Check if imod is in the path
        """

        if not cls._validationMsg:

            etomo = cls._getProgram(ETOMO_CMD)

            cls._validationMsg = [
                "imod's %s command not found in path, please install it." % etomo] if not find_executable(
                ETOMO_CMD) and not os.path.exists(etomo) else []

        return cls._validationMsg

    @classmethod
    def getDependencies(cls):
        neededPrograms = ['java', 'python']

        return neededPrograms

    @classmethod
    def defineBinaries(cls, env):
        IMOD_INSTALLED = 'imod_%s_installed' % DEFAULT_VERSION

        if 'linux' in OS.getPlatform().lower():

            # Add jpg lib
            jpeg = env.addLibrary(
                'jpeg',
                tar='libjpeg-turbo-1.3.1.tgz',
                flags=['--without-simd'],
                default=False)

            # Download .sh
            installationCmd = 'wget --continue https://bio3d.colorado.edu/ftp/latestIMOD/RHEL6-64_CUDA8.0/' \
                              'imod_4.10.42_RHEL6-64_CUDA8.0.sh --no-check-certificate && '

            # Run .sh skipping copying startup scripts (avoid sudo permissions to write to /etc/profile.d)
            installationCmd += 'sh imod_4.10.42_RHEL6-64_CUDA8.0.sh -dir . -yes -skip && '

            # Create installation finished flag file
            installationCmd += 'touch %s' % IMOD_INSTALLED

            env.addPackage('imod',
                           deps=[jpeg],
                           version=DEFAULT_VERSION,
                           tar='void.tgz',
                           createBuildDir=True,
                           buildDir=cls._getEMFolder(DEFAULT_VERSION),
                           neededProgs=cls.getDependencies(),
                           libChecks="libjpeg62",
                           commands=[(installationCmd, IMOD_INSTALLED)],
                           default=True)

    @classmethod
    def runImod(cls, protocol, program, args, cwd=None):
        """ Run IMOD command from a given protocol. """

        # Get the command
        cmd = cls.getImodCmd(program)

        # Run the protocol with that command
        protocol.runJob(cmd, args, env=cls.getEnviron(), cwd=cwd)

    @classmethod
    def getImodCmd(cls, program):
        """ Composes an IMOD command for a given program. """

        # Program to run
        program = cls._getProgram(program)

        # Command to run
        cmd = ""

        # If absolute ... (then it is based on the config)
        if os.path.isabs(program):
            cmd += ". " + cls.getHome("IMOD-linux.sh") + " && "

        cmd += program

        return cmd
