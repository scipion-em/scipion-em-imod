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
import logging
logger = logging.getLogger(__name__)

import pyworkflow.viewer as pwviewer
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pwem.viewers import DataViewer
import tomo.objects as tomoObj

import imod.protocols
from imod.constants import (OUTPUT_TILTSERIES_NAME, OUTPUT_FIDUCIAL_NO_GAPS_NAME,
                            OUTPUT_TS_COORDINATES_NAME)
from imod.viewers.views_tkinter_tree import ImodGenericView
from imod import Plugin
from imod.utils import generateIMODFidFile


class ImodViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Imod program 3dmod
    """
    _environments = [pwviewer.DESKTOP_TKINTER, Plugin.getEnviron()]
    _targets = [
        tomoObj.TiltSeries,
        tomoObj.Tomogram,
        tomoObj.SetOfTomograms,
        tomoObj.SetOfTiltSeries,
        tomoObj.SetOfLandmarkModels,
        tomoObj.LandmarkModel
    ]
    _name = 'Imod'

    def _visualize(self, obj, **kwargs):
        env = Plugin.getEnviron()
        cls = type(obj)

        if issubclass(cls, (tomoObj.TiltSeries, tomoObj.Tomogram, tomoObj.LandmarkModel)):
            binning = kwargs.get("binning", 1)
            view = ImodObjectView(obj, protocol=self.protocol, binning=binning)
        else:  # Set object
            view = ImodGenericView(self.getTkRoot(), self.protocol, obj)

        view._env = env
        return [view]


class ImodObjectView(pwviewer.CommandView):
    """ Wrapper to visualize different type of objects with the 3dmod """

    def __init__(self, obj, protocol=None, binning=1, **kwargs):
        """
        :param obj: Object to deal with, a single item of a set
        :param protocol: protocol owner of obj
        :param kwargs: extra kwargs
        """
        # Get default binning level has been defined for 3dmod
        binningstr = str(binning)
        cmd = f"{Plugin.getImodCmd('3dmod')} "

        if isinstance(obj, tomoObj.TiltSeries):
            prj = protocol.getProject()
            inputFn = obj.getFirstItem().getFileName()
            angleFilePath = prj.getTmpPath(pwutils.replaceBaseExt(inputFn, "tlt"))
            obj.generateTltFile(angleFilePath)

            cmd += f"-b {binningstr},1 "
            cmd += f"-a {angleFilePath} {inputFn.split(':')[0]}"

        elif isinstance(obj, tomoObj.LandmarkModel):
            prj = protocol.getProject()
            ts = obj.getTiltSeries()
            tsFn = ts.getFirstItem().getFileName()
            if ts.hasAlignment() and obj.applyTSTransformation():
                # Input and output extensions must match if we want to apply the transform with Xmipp
                extension = pwutils.getExt(tsFn)
                outputTSPath = prj.getTmpPath("ts_interpolated_%s_%s_%s%s" % (
                                                prj.getShortName(),
                                                protocol.getObjId(),
                                                obj.getObjId(),
                                                extension))

                if not os.path.exists(outputTSPath):
                    ts.applyTransform(outputTSPath)

            else:
                outputTSPath = tsFn

            fidFileName = obj.getModelName()

            if fidFileName is None:
                fidFileName = generateIMODFidFile(protocol, obj)

            angleFilePath = prj.getTmpPath(pwutils.replaceBaseExt(tsFn, "tlt"))
            ts.generateTltFile(angleFilePath)

            cmd += f"-a {angleFilePath} -m {outputTSPath} {fidFileName}"

        # A path called from the object browser
        elif isinstance(obj, str):
            cmd += f"{obj}"

        else:  # Tomogram
            cmd += f"-b {binningstr},{binningstr} "
            cmd += f"{obj.getFileName()}"

        logger.info(f"Executing command: {cmd}")
        pwviewer.CommandView.__init__(self,  cmd)


class ImodEtomoViewer(pwviewer.ProtocolViewer):
    """ Viewer form for Etomo interactive results. """

    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [imod.protocols.ProtImodEtomo]
    _label = 'viewer etomo'

    def _defineParams(self, form):
        form.addSection(label='Visualization')

        group = form.addGroup('Tilt Series Alignment')
        group.addParam('savedTsPreAli', params.LabelParam,
                       label="Pre-aligned tilt-series")
        group.addParam('savedTsAli', params.LabelParam,
                       label="Aligned tilt-series")
        group.addParam('saved3DCoord', params.LabelParam,
                       label="3D Coordinates")
        group.addParam('savedFiducials', params.LabelParam,
                       label="Landmark models with no gaps")

        group = form.addGroup('Tomogram')
        group.addParam('savedReconsTomo', params.LabelParam,
                       label="Reconstructed raw tomogram")
        group.addParam('savedPostProcessTomo', params.LabelParam,
                       label="Post-processed tomogram")

        self.defineOutputsSetNames()

    def defineOutputsSetNames(self, **kwargs):
        self.outputSetName = {'savedTsPreAli': 'PrealignedTiltSeries',
                              'savedTsAli': OUTPUT_TILTSERIES_NAME,
                              'saved3DCoord': OUTPUT_TS_COORDINATES_NAME,
                              'savedFiducials': OUTPUT_FIDUCIAL_NO_GAPS_NAME,
                              'savedReconsTomo': 'FullTomograms',
                              'savedPostProcessTomo': 'PostProcessTomograms'}

    def _getVisualizeDict(self):
        return {'savedTsPreAli': self._showOutputSet,
                'savedTsAli': self._showOutputSet,
                'saved3DCoord': self._showOutputSet,
                'savedFiducials': self._showOutputSet,
                'savedReconsTomo': self._showOutputSet,
                'savedPostProcessTomo': self._showOutputSet,
                }

    def _showOutputSet(self, param=None):
        try:
            outputName = self.outputSetName.get(param)
            if hasattr(self.protocol, outputName):
                outputSet = getattr(self.protocol, outputName)
                if param == 'saved3DCoord':
                    dataviewer = DataViewer(protocol=self.protocol,
                                            project=self.protocol.getProject())
                    dataviewer._visualize(outputSet)[0].show()
                else:
                    ImodGenericView(self.getTkRoot(), self.protocol, outputSet).show()
            else:
                self._notGenerated()

        except Exception as e:
            return [self.errorMessage(str(e), "Error displaying the output")]

    def _notGenerated(self, param=None):
        return [self.infoMessage('Outputs are not generated yet.', 'Info').show()]
