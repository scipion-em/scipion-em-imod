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

import pyworkflow.viewer as pwviewer
import pyworkflow.protocol.params as params

import tomo.objects
import imod.protocols


class ImodViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        tomo.objects.TiltSeries,
        tomo.objects.Tomogram,
        tomo.objects.SetOfTomograms,
    ]

    def _visualize(self, obj, **kwargs):
        views = []
        cls = type(obj)

        if issubclass(cls, tomo.objects.TiltSeries):
            views.append(ImodObjectView(obj.getFirstItem()))

        elif issubclass(cls, tomo.objects.Tomogram):
            views.append(ImodObjectView(obj))
        elif issubclass(cls, tomo.objects.SetOfTomograms):
            for t in obj:
                views.append(ImodObjectView(t))

        return views
    

class ImodObjectView(pwviewer.CommandView):
    """ Wrapper to visualize different type of objects with the 3dmod.
    """
    def __init__(self, obj, **kwargs):
        # Remove :mrc if present
        fn = obj.getFileName().split(':')[0]
        pwviewer.CommandView.__init__(self, '3dmod "%s"' % fn)


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

