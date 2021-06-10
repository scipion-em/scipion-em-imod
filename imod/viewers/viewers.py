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
import os

import pyworkflow.viewer as pwviewer
from imod.viewers.views_tkinter_tree import (CtfEstimationTreeProvider,
                                             CtfEstimationListDialog,
                                             ImodGenericViewer)
import pyworkflow.protocol.params as params

import tomo.objects
import imod.protocols
from imod import Plugin


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
        # Remove :mrc if present
        if isinstance(obj, tomo.objects.LandmarkModel):
            tsId = os.path.basename(obj.getFileName()).split('_')[0]
            if os.path.exists(os.path.join(os.path.split(obj.getModelName())[0],
                                           "%s_preali.st" % tsId)):
                prealiTSPath = os.path.join(os.path.split(obj.getModelName())[0],
                                            "%s_preali.st" % tsId)
            elif os.path.exists(os.path.join(os.path.split(obj.getModelName())[0],
                                "%s.preali" % tsId)):
                prealiTSPath = os.path.join(os.path.split(obj.getModelName())[0],
                                            "%s.preali" % tsId)
            else:
                prealiTSPath = ""

            fn = Plugin.getImodCmd('3dmod') + " -m " + prealiTSPath + " " + obj.getModelName() + " ; "

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
        group.addParam('saveTsPreAli', params.LabelParam,
                       label="Register pre-aligned tilt-series (.preali)",
                       help="Through this option the intermediate pre-aligned "
                            "tilt-series can be registered as an output of "
                            "this protocol.")
        group.addParam('saveTsAli', params.LabelParam,
                       label="Register aligned tilt-series (.ali)",
                       help="Through this option the intermediate aligned "
                            "tilt-series can be registered as an output of "
                            "this protocol.")
        group.addParam('saveTsOriginal', params.LabelParam,
                       label="Register original tilt-series (+ alignment info)",
                       help="Through this option the original tilt-series can "
                            "be registered as an output of this protocol. The "
                            "obtained alignment parameters will be stored as "
                            "metadata. ")
        group.addParam('saveFiducials', params.LabelParam,
                       label="Register fiducials model",
                       help="Through this option the obtained fiducial model "
                            "can be registered as an output of this protocol. "
                            "The obtained alignment parameters will be stored "
                            "as metadata. ")

        group = form.addGroup('Tomogram')
        group.addParam('saveReconsTomo', params.LabelParam,
                       label="Register reconstructed tomogram (_full.rec)",
                       help="Through this option the final reconstructed "
                            "tomogram can be registered as an output of "
                            "this protocol.")

        form.addParam('doShowResHistogram', params.LabelParam,
                      label="Show resolution histogram")

    def _getVisualizeDict(self):
        return {'saveTsPreAli': self._registerOutput,
                'saveTsAli': self._registerOutput,
                'saveTsOriginal': self._notImplemented,
                'saveFiducials': self._notImplemented,
                'saveReconsTomo': self._registerOutput,
                }

    def _registerOutput(self, param=None):
        try:
            self.protocol.registerOutput(param)
            return [self.infoMessage('Output registered successfully. ')]
        except Exception as e:
            return [self.errorMessage(str(e))]

    def _notImplemented(self, param=None):
        return [self.errorMessage('Output not implemented yet. ')]


class CtfEstimationTomoViewer(pwviewer.Viewer):
    """ This class implements a view using Tkinter ListDialog
    and the CtfEstimationTreeProvider.
    """
    _label = 'ctf estimation viewer'
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [tomo.objects.SetOfCTFTomoSeries]

    def __init__(self, parent, protocol, **kwargs):
        self._tkParent = parent.root
        self._protocol = protocol
        self._title = 'ctf estimation viewer'

    def visualize(self, obj, windows=None, protocol=None):
        objName = obj.getObjName().split('.')[1]
        for output in self._protocol._iterOutputsNew():
            if output[0] == objName:
                self._outputSetOfCTFTomoSeries = output[1]
                break
        self._inputSetOfTiltSeries = self._outputSetOfCTFTomoSeries.getSetOfTiltSeries()
        self._provider = CtfEstimationTreeProvider(self._tkParent,
                                                   self._protocol,
                                                   self._outputSetOfCTFTomoSeries)
        CtfEstimationListDialog(self._tkParent, self._title, self._provider,
                                self._protocol, self._inputSetOfTiltSeries)