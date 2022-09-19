# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [2]
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


import tempfile
import os

import pyworkflow.viewer as pwviewer
from imod.protocols.protocol_base import OUTPUT_TILTSERIES_NAME, OUTPUT_COORDINATES_3D_NAME, \
    OUTPUT_FIDUCIAL_NO_GAPS_NAME
from imod.viewers.views_tkinter_tree import ImodGenericViewer, ImodSetView, \
    ImodSetOfLandmarkModelsView, ImodSetOfTomogramsView
import pyworkflow.protocol.params as params

import tomo.objects
import imod.protocols
from imod import Plugin
from pwem.viewers import DataViewer


class ImodViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Imod program 3dmod
    """
    _environments = [pwviewer.DESKTOP_TKINTER, Plugin.getEnviron()]
    _targets = [
        tomo.objects.TiltSeries,
        tomo.objects.Tomogram,
        tomo.objects.SetOfTomograms,
        tomo.objects.SetOfTiltSeries,
        tomo.objects.SetOfLandmarkModels,
        tomo.objects.LandmarkModel
    ]

    def _visualize(self, obj, **kwargs):
        env = Plugin.getEnviron()
        cls = type(obj)

        if issubclass(cls, tomo.objects.TiltSeries):
            view = ImodObjectView(obj.getFirstItem())
        elif issubclass(cls, tomo.objects.Tomogram):
            view = ImodObjectView(obj)
        elif issubclass(cls, tomo.objects.LandmarkModel):
            view = ImodObjectView(obj)
        else:
            view = ImodGenericViewer(self.getTkRoot(), self.protocol, obj)

        view._env = env
        return [view]


class ImodObjectView(pwviewer.CommandView):
    """ Wrapper to visualize different type of objects with the 3dmod """

    def __init__(self, obj, **kwargs):
        # Accept file paths
        if isinstance(obj, str):
            fn = Plugin.getImodCmd('3dmod') + ' ' + obj

        elif isinstance(obj, tomo.objects.LandmarkModel):
            if obj.getTiltSeries().getFirstItem().hasTransform():
                # Input and output extensions must match if we want to apply the transform with Xmipp
                _, extension = os.path.splitext(obj.getTiltSeries().getFirstItem().getFileName())

                outputTSInterpolatedPath = os.path.join(tempfile.gettempdir(), "ts_interpolated." + extension)
                obj.getTiltSeries().applyTransform(outputTSInterpolatedPath)

                fn = Plugin.getImodCmd('3dmod') + " -m " + outputTSInterpolatedPath + " " + \
                      obj.getModelName() + " ; "

            else:
                fn = Plugin.getImodCmd('3dmod') + " -m " + obj.getTiltSeries().getFirstItem().getFileName() + \
                      " " + obj.getModelName() + " ; "

        else:
            fn = Plugin.getImodCmd('3dmod') + ' ' + obj.getFileName().split(':')[0]

        pwviewer.CommandView.__init__(self,  fn)


class ImodEtomoViewer(pwviewer.ProtocolViewer):
    """ Viewer form for Etomo interactive results. """

    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [imod.protocols.ProtImodEtomo]
    _label = 'viewer etomo'

    def _defineParams(self, form):
        form.addSection(label='Visualization')

        group = form.addGroup('Tilt Series Alignment')
        group.addParam('savedTsPreAli', params.LabelParam,
                       label="Pre-aligned tilt-series",
                       help="Through this option the intermediate pre-aligned "
                            "tilt-series can be shown.")
        group.addParam('savedTsAli', params.LabelParam,
                       label="Aligned tilt-series",
                       help="Through this option the intermediate aligned "
                            "tilt-series can be shown.")
        group.addParam('saved3DCoord', params.LabelParam,
                       label="3D Coordinates",
                       help="Through this option the 3D coordinates can "
                            "be shown.")
        group.addParam('savedFiducials', params.LabelParam,
                       label="Landmark models no gaps",
                       help="Through this option the obtained fiducial model "
                            "can be shown.")

        group = form.addGroup('Tomogram')
        group.addParam('savedReconsTomo', params.LabelParam,
                       label="Reconstructed full tomogram",
                       help="Through this option the final reconstructed "
                            "tomogram can be shown.")
        group.addParam('savedPostProcessTomo', params.LabelParam,
                       label="Postprocess tomogram",
                       help="Through this option the postprocess "
                            "tomogram can be shown.")

        self.defineOutputsSetNames()

    def defineOutputsSetNames(self, **kwargs):
        self.outputSetName = {'savedTsPreAli': 'PrealignedTiltSeries',
                              'savedTsAli': OUTPUT_TILTSERIES_NAME,
                              'saved3DCoord': OUTPUT_TILTSERIES_NAME,
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
                if param == 'savedTsPreAli' or param == 'savedTsAli':
                    ImodSetView(outputSet).show()
                elif param == 'savedFiducials':
                    ImodSetOfLandmarkModelsView(outputSet).show()
                elif param == 'savedReconsTomo' or param == 'savedPostProcessTomo':
                    ImodSetOfTomogramsView(outputSet).show()
                elif param == 'saved3DCoord':
                    dataviewer = DataViewer(protocol=self.protocol,
                                            project=self.protocol.getProject())
                    dataviewer._visualize(outputSet)[0].show()
            else:
                self._notGenerated()

        except Exception as e:
            return [self.errorMessage(str(e), "Error displaying the output")]

    def _notGenerated(self, param=None):
       return [self.infoMessage('Output not generated yet. ', 'Info').show()]
