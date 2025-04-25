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

from pyworkflow.protocol.constants import STEPS_SERIAL
from pyworkflow.gui import FileTreeProvider
from pyworkflow.gui.project.utils import OS
import pwem

from imod.constants import (IMOD_HOME, ETOMO_CMD, DEFAULT_VERSION,
                            VERSIONS, IMOD_VIEWER_BINNING, BRT_ENV_ACTIVATION, BRT_DEFAULT_ACTIVATION_CMD, BRT_CUDA_LIB,
                            BRT, BRT_PROGRAM_DEFAULT_VERSION, BRT_ENV_NAME, BRT_PROGRAM)

__version__ = '3.7.0'

from pyworkflow.utils import Environ

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
        cls._defineVar(IMOD_VIEWER_BINNING, 1)
        cls._defineVar(BRT_ENV_ACTIVATION, BRT_DEFAULT_ACTIVATION_CMD)
        cls._defineVar(BRT_CUDA_LIB, pwem.Config.CUDA_LIB)

    @classmethod
    def getBRTEnvActivation(cls):
        return cls.getVar(BRT_ENV_ACTIVATION)

    @classmethod
    def getViewerBinning(cls):
        return cls.getVar(IMOD_VIEWER_BINNING)

    @classmethod
    def _getEMFolder(cls, version, *paths):
        return os.path.join("imod-%s" % version, *paths)

    @classmethod
    def _getIMODFolder(cls, version, *paths):
        return cls._getEMFolder(version, "IMOD", *paths)

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
        for version in VERSIONS:
            cls.installImod(env, version, version == DEFAULT_VERSION)
        # Install yet-another-imod-wrapper
        cls.installBatchRunTomo(env)

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
            installationCmd = 'wget --continue https://bio3d.colorado.edu/imod/AMD64-RHEL5/' \
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
    def installBatchRunTomo(cls, env):
        BRT_INSTALLED = '%s_%s_installed' % (BRT, BRT_PROGRAM_DEFAULT_VERSION)
        installationCmd = cls.getCondaActivationCmd()
        # Create the environment
        installationCmd += ' conda create -y -n %s -c conda-forge python=3.8 && ' % BRT_ENV_NAME

        # Activate new the environment
        installationCmd += 'conda activate %s && ' % BRT_ENV_NAME

        # Install BRT
        installationCmd += f'pip install {BRT_PROGRAM}=={BRT_PROGRAM_DEFAULT_VERSION} && '

        # Flag installation finished
        installationCmd += 'touch %s' % BRT_INSTALLED

        BRT_commands = [(installationCmd, BRT_INSTALLED)]
        envPath = os.environ.get('PATH', "")  # keep path since conda likely in there
        installEnvVars = {'PATH': envPath} if envPath else None

        env.addPackage(BRT,
                       version=BRT_PROGRAM_DEFAULT_VERSION,
                       tar='void.tgz',
                       commands=BRT_commands,
                       neededProgs=cls.getDependenciesBRT(),
                       vars=installEnvVars,
                       default=True)

    @classmethod
    def getDependenciesBRT(cls):
        # try to get CONDA activation command
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = []
        if not condaActivationCmd:
            neededProgs.append('conda')
        return neededProgs
    
    @classmethod
    def runBRT(cls, protocol, args, cwd=None, numberOfMpi=1):
        """ Run yet-another-imod-wrapper (batchruntomo) command from a given protocol. """
        cmd = cls.getCondaActivationCmd() + " "
        cmd += cls.getBRTEnvActivation()
        cmd += f"&& export PATH={cls.getHome('bin')}:$PATH "
        cmd += f"&& export IMOD_DIR={cls.getHome()} "
        cmd += f"&& {BRT_PROGRAM} "
        protocol.runJob(cmd, args, env=cls.getEnviron(), cwd=cwd, numberOfMpi=numberOfMpi)

    @classmethod
    def runImod(cls, protocol, program, args, cwd=None):
        """ Run IMOD command from a given protocol. """

        if protocol.stepsExecutionMode == STEPS_SERIAL:
            ncpus = protocol.numberOfThreads.get()
        else:
            ncpus = 1

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
