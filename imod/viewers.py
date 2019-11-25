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

import pyworkflow.viewer as pwviewer

import tomo.objects


class ImodViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        tomo.objects.TiltSeries,
        tomo.objects.Tomogram,
        tomo.objects.SetOfTomograms,
        tomo.objects.SetOfTiltSeries,
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
        elif issubclass(cls, tomo.objects.SetOfTiltSeries):
            views.append(ImodSetView(obj))

        return views
    

class ImodObjectView(pwviewer.CommandView):
    """ Wrapper to visualize different type of objects with the 3dmod.
    """
    def __init__(self, obj, **kwargs):
        # Remove :mrc if present
        fn = obj.getFileName().split(':')[0]
        pwviewer.CommandView.__init__(self, '3dmod %s' % fn)


class ImodSetView(pwviewer.CommandView):
    """ Wrapper to visualize different type of objects with the 3dmod.
    """
    def __init__(self, set, **kwargs):
        fn = ""
        for item in set:
            # Remove :mrc if present
            fn += " " + item.getFirstItem().getFileName().split(':')[0]
        print(fn)
        pwviewer.CommandView.__init__(self, '3dmod%s' % fn)



