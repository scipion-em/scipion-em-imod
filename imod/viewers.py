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
from tkinter import *
import matplotlib.pyplot as plt

import pyworkflow.viewer as pwviewer
from pyworkflow.gui import *
from pyworkflow.gui.tree import TreeProvider
from pyworkflow.gui.dialog import ListDialog
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
    ]

    def _visualize(self, obj, **kwargs):
        env = Plugin.getEnviron()
        view = []
        cls = type(obj)

        if issubclass(cls, tomo.objects.TiltSeries):
            view = ImodObjectView(obj.getFirstItem())
        elif issubclass(cls, tomo.objects.Tomogram):
            view = ImodObjectView(obj)
        elif issubclass(cls, tomo.objects.SetOfTomograms):
            view = ImodSetOfTomogramsView(obj)
        elif issubclass(cls, tomo.objects.SetOfTiltSeries):
            view = ImodSetView(obj)
        elif issubclass(cls, tomo.objects.SetOfLandmarkModels):
            view = ImodSetOfLandmarkModelsView(obj)

        view._env = env
        return [view]


class ImodObjectView(pwviewer.CommandView):
    """ Wrapper to visualize different type of objects with the 3dmod """

    def __init__(self, obj, **kwargs):
        # Remove :mrc if present
        fn = obj.getFileName().split(':')[0]
        pwviewer.CommandView.__init__(self, Plugin.getImodCmd('3dmod') + ' ' + fn)


class ImodSetView(pwviewer.CommandView):
    """ Wrapper to visualize different type of objects with the 3dmod """

    def __init__(self, set, **kwargs):
        fn = ""
        for item in set:
            # Remove :mrc if present
            fn += " " + item.getFirstItem().getFileName().split(':')[0]
        pwviewer.CommandView.__init__(self, "%s %s" % (Plugin.getImodCmd('3dmod'), fn))


class ImodSetOfLandmarkModelsView(pwviewer.CommandView):
    """ Wrapper to visualize landmark models with 3dmod """

    def __init__(self, set, **kwargs):
        fn = ""
        for item in set:
            tsId = os.path.basename(item.getFileName()).split('_')[0]
            if os.path.exists(os.path.join(os.path.split(item.getModelName())[0], "%s_preali.st" % tsId)):
                prealiTSPath = os.path.join(os.path.split(item.getModelName())[0], "%s_preali.st" % tsId)
            elif os.path.exists(os.path.join(os.path.split(item.getModelName())[0], "%s.preali" % tsId)):
                prealiTSPath = os.path.join(os.path.split(item.getModelName())[0], "%s.preali" % tsId)
            else:
                prealiTSPath = ""
            fn += Plugin.getImodCmd('3dmod') + " -m " + prealiTSPath + " " + item.getModelName() + " ; "
        pwviewer.CommandView.__init__(self, fn)


class ImodSetOfTomogramsView(pwviewer.CommandView):
    """ Wrapper to visualize set of tomograms with 3dmod """

    def __init__(self, set, **kwargs):
        fn = " -s 0,0 "
        for item in set:
            fn += " " + item.getLocation()[1]
        pwviewer.CommandView.__init__(self, Plugin.getImodCmd('3dmod') + fn)


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


class CTFSerieStates:
    UNCHECKED = 'unchecked'
    CHECKED = 'checked'


class CtfEstimationTreeProvider(TreeProvider, ttk.Treeview):
    """ Model class that will retrieve the information from TiltSeries and
    prepare the columns/rows models required by the TreeDialog GUI.
    """
    COL_CTF_SERIE = 'CTF Series'
    COL_CTF_SERIE_ENABLE = 'enabled'
    CRITERIA_1 = 'criteria1'
    CRITERIA_2 = 'criteria2'
    COL_CTF_EST_IX = 'index'
    COL_CTF_EST_DEFOCUS_U = 'defocusU'
    COL_CTF_EST_DEFOCUS_V = 'defocusV'
    COL_CTF_EST_DEFOCUS_RATIO = 'defocusRatio'
    COL_CTF_EST_DEFOCUS_ANGLE = 'defocusAngle'
    COL_CTF_EST_AST = 'ast'
    COL_CTF_EST_RES = 'resolution'
    COL_CTF_EST_FIT = 'fitQuality'

    ORDER_DICT = {COL_CTF_EST_DEFOCUS_U: '_defocusU',
                  COL_CTF_EST_DEFOCUS_V: '_defocusV',
                  COL_CTF_EST_DEFOCUS_RATIO: '_defocusRatio',
                  COL_CTF_EST_DEFOCUS_ANGLE: '_defocusAngle',
                  COL_CTF_EST_RES: '_resolution',
                  COL_CTF_EST_FIT: '_fitQuality',
    }

    def __init__(self, master, protocol, outputSetOfCTFTomoSeries, **kw):

        ttk.Treeview.__init__(self, master, **kw)

        self.im_checked = gui.getImage(Icon.CHECKED)
        self.im_unchecked = gui.getImage(Icon.UNCHECKED)

        self.tag_configure(CTFSerieStates.UNCHECKED, image=self.im_unchecked)
        self.tag_configure(CTFSerieStates.CHECKED, image=self.im_checked)

        self.protocol = protocol
        self.ctfSeries = outputSetOfCTFTomoSeries
        TreeProvider.__init__(self, sortingColumnName=self.COL_CTF_SERIE)
        self.selectedDict = {}
        self.mapper = protocol.mapper
        self.maxNum = 200


    def insert(self, parent, index, iid=None, **kw):
        """ same method as for standard treeview but add the tag 'unchecked'
            automatically if no tag among ('checked', 'unchecked')
            is given """
        if "tags" not in kw:
            kw["tags"] = (CTFSerieStates.UNCHECKED,)
        ttk.Treeview.insert(self, parent, index, iid, **kw)

    def getObjects(self):
        # Retrieve all objects of type className
        project = self.protocol.getProject()
        objects = []

        orderBy = self.ORDER_DICT.get(self.getSortingColumnName(), 'id')
        direction = 'ASC' if self.isSortingAscending() else 'DESC'

        for ctfSerie in self.ctfSeries:
            ctfEstObj = ctfSerie.clone()
            ctfEstObj._allowsSelection = True
            ctfEstObj._parentObject = None
            objects.append(ctfEstObj)
            for item in ctfSerie.iterItems(orderBy=orderBy, direction=direction):
                ctfEstItem = item.clone()
                ctfEstItem._allowsSelection = False
                ctfEstItem._parentObject = ctfEstObj
                objects.append(ctfEstItem)

        return objects

    def _sortObjects(self, objects):
        pass

    def objectKey(self, pobj):
        pass

    def getColumns(self):
        cols = [
            (self.COL_CTF_SERIE, 300),
            (self.COL_CTF_SERIE_ENABLE, 70),
            (self.CRITERIA_1, 150),
            (self.CRITERIA_2, 150),
            (self.COL_CTF_EST_IX, 150),
            (self.COL_CTF_EST_DEFOCUS_U, 150),
            (self.COL_CTF_EST_DEFOCUS_V, 150),
            (self.COL_CTF_EST_DEFOCUS_RATIO, 150),
            (self.COL_CTF_EST_DEFOCUS_ANGLE, 150)
        ]
        return cols

    def isSelected(self, obj):
        """ Check if an object is selected or not. """
        return False

    @staticmethod
    def _getParentObject(pobj, default=None):
        return getattr(pobj, '_parentObject', default)

    def getObjectInfo(self, obj):
        objId = obj.getObjId()

        if isinstance(obj, tomo.objects.CTFTomoSeries):
            key =  obj.getTsId()
            text = obj.getTsId()
            values = ['', '', '', '', '', '']
            opened = False
        else:  # CTFTomo
            key = str(obj.getObjId())
            text = ''
            values = ['', '', '', obj.getIndex(), str(obj.getDefocusU()),
                      str(obj.getDefocusV()), str(obj.getDefocusRatio()),
                      str(obj.getDefocusAngle())]
            opened = False

        item = {
            'key': key, 'text': text,
            'values': tuple(values),
            'open': opened,
            'selected': False,
            'parent': obj._parentObject
        }
        if isinstance(obj, tomo.objects.CTFTomoSeries):
            item['tags'] = (CTFSerieStates.UNCHECKED,)

        return item

    def getObjectActions(self, obj):
        actions = []
        defocusUList = []
        defocusVList = []
        defocusRatioList = []
        if isinstance(obj, tomo.objects.CTFTomoSeries):
            for ctfSerie in self.ctfSeries:
                if ctfSerie.getObjId() == obj.getObjId():
                    for item in ctfSerie.iterItems(orderBy='id'):
                        defocusUList.append(item.getDefocusU())
                        defocusVList.append(item.getDefocusV())

                    plt.figure('CTF Estimation Plotter')
                    plt.plot(defocusUList, marker='o', label='DefocusU')
                    plt.plot(defocusVList, marker='o', label='DefocusV')
                    plt.legend()
                    plt.title("CTF Estimation Plotter")
                    plt.show()
                    break
        return actions

class CtfEstimationDialogView:
    """ This class implements a view using Tkinter ListDialog
    and the CtfEstimationTreeProvider.
    """

    def __init__(self, parent, protocol, outputSetOfCTFTomoSeries, **kwargs):
        self._tkParent = parent
        self._protocol = protocol
        self._title = 'ctf estimation viewer'
        self._outputSetOfCTFTomoSeries = outputSetOfCTFTomoSeries
        self._provider = CtfEstimationTreeProvider(self._tkParent, self._protocol, self._outputSetOfCTFTomoSeries)

    def show(self):
        dlg = ListDialog(self._tkParent, self._title, self._provider)


class CtfEstimationViewer(pwviewer.ProtocolViewer):
    """ Wrapper to visualize outputs of tilt series motion correction protocols
    """

    _label = 'ctf estimation viewer'
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [imod.protocols.ProtImodCtfEstimation]

    def _defineParams(self, form):
            form.addSection(label='Visualization of ctf estimation series')
            form.addParam('displayFullCTFEstSeries', params.LabelParam,
                          label='Display full ctf estimation series',
                          help='Shows full ctf estimation series'
                          )

    def getOutputSetOfCTFTomoSeries(self):
        return getattr(self.protocol, 'outputSetOfCTFTomoSeries')

    def _displayFullTiltSeries(self, param=None):
        return self._visualize(self.getOutputSetOfCTFTomoSeries())

    def _getVisualizeDict(self):
        return {
            'displayFullCTFEstSeries': self._displayFullTiltSeries
        }

    def _visualize(self, outputSetOfCTFTomoSeries):

        setCTFEstView = CtfEstimationDialogView(self.getTkRoot(), self.protocol,
                                                 outputSetOfCTFTomoSeries)


        return [setCTFEstView]

