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

import os
from shutil import which

from pyworkflow.gui import FileTreeProvider
from pyworkflow.gui.project.utils import OS
import pwem

from .constants import IMOD_HOME, ETOMO_CMD, DEFAULT_VERSION, VERSIONS


__version__ = '3.1.8'
_logo = "icon.png"
_references = ['Kremer1996', 'Mastronarde2017']


def getImodEnv():
    """ This function allows to call imod outside this plugin. """

    return Plugin.getHome("IMOD-linux.sh && ")


class Plugin(pwem.Plugin):
    _homeVar = IMOD_HOME
    _validationMsg = None
    _url = "https://github.com/scipion-em/scipion-em-imod"

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(IMOD_HOME, cls._getIMODFolder(DEFAULT_VERSION))

    @classmethod
    def _getEMFolder(cls, version, *paths):
        return os.path.join("imod-%s" % version, *paths)

    @classmethod
    def _getIMODFolder(cls, version, *paths):
        return os.path.join(cls._getEMFolder(version, "IMOD"), *paths)

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
        env = pwem.pwutils.Environ(os.environ)
        if 'IMOD_DIR' in env:
            del env['IMOD_DIR']
        if 'IMOD_PATH' in env:
            del env['IMOD_PATH']
        return env

    @classmethod
    def validateInstallation(cls):
        """ Check if imod is in the path """

        if not cls._validationMsg:
            etomo = cls._getProgram(ETOMO_CMD)

            cls._validationMsg = [
                "IMOD's %s command not found in path, please "
                "install it." % etomo] if not which(
                ETOMO_CMD) and not os.path.exists(etomo) else []

        return cls._validationMsg

    @classmethod
    def getDependencies(cls):
        neededPrograms = ['java', 'python']

        return neededPrograms

    @classmethod
    def defineBinaries(cls, env):
        version = VERSIONS[-1]
        cls.installImod(env, version, version == DEFAULT_VERSION)

    @classmethod
    def installImod(cls, env, version, default):
        IMOD_INSTALLED = 'imod_%s_installed' % version
        if 'linux' in OS.getPlatform().lower():
            # Add jpg lib, once
            JPEG_NAME = 'jpeg'

            if not env.hasTarget(JPEG_NAME):
                jpeg = env.addLibrary(
                    JPEG_NAME,
                    tar='libjpeg-turbo-1.3.1.tgz',
                    flags=['--without-simd'],
                    default=False)
            else:
                jpeg = env.getTarget(JPEG_NAME)

            # Download .sh
            installationCmd = 'wget --continue http://bio3d.colorado.edu/imod/AMD64-RHEL5/' \
                              'imod_%s_RHEL7-64_CUDA10.1.sh --no-check-certificate && ' % version

            # Run .sh skipping copying startup scripts (avoid sudo permissions to write to /etc/profile.d)
            installationCmd += 'sh imod_%s_RHEL7-64_CUDA10.1.sh -dir . -yes -skip && ' % version

            # Create installation finished flag file
            installationCmd += 'touch %s' % IMOD_INSTALLED

            env.addPackage('imod',
                           deps=[jpeg],
                           version=version,
                           tar='void.tgz',
                           createBuildDir=True,
                           buildDir=cls._getEMFolder(version),
                           neededProgs=cls.getDependencies(),
                           commands=[(installationCmd, IMOD_INSTALLED)],
                           default=default)

    @classmethod
    def runImod(cls, protocol, program, args, cwd=None):
        """ Run IMOD command from a given protocol. """

        ncpus = protocol.numberOfThreads.get()

        # Get the command
        cmd = cls.getImodCmd(program, ncpus)

        # Run the protocol with that command
        protocol.runJob(cmd, args, env=cls.getEnviron(), cwd=cwd,
                        numberOfMpi=1, numberOfThreads=1)

    @classmethod
    def getImodCmd(cls, program, ncpus=1):
        """ Composes an IMOD command for a given program. """

        # Program to run
        program = cls._getProgram(program)

        # Command to run
        cmd = ""

        if ncpus > 1:
            cmd += f"export IMOD_PROCESSORS={ncpus} && "

        # If absolute ... (then it is based on the config)
        if os.path.isabs(program):
            cmd += ". " + cls.getHome("IMOD-linux.sh") + " && "

        cmd += program

        return cmd


from .file_handlers import *

register = FileTreeProvider.registerFileHandler
register(ImodHandler(), '.ali', '.st', '.rec', '.mrc', '.mrcs')
