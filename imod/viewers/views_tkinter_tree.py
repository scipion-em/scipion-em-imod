# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import threading
import tkinter
from tkinter import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from imod.protocols import ProtImodEtomo
from pyworkflow.gui import *
from pyworkflow.gui.tree import TreeProvider
from pyworkflow.gui.dialog import ListDialog, showInfo
import pyworkflow.viewer as pwviewer
from pyworkflow.plugin import Domain

import tomo.objects
from imod import Plugin


class protClass:
    protImodEtomoClass = 1
    protImodCTFEstimation = 2

class ImodGenericTreeProvider(TreeProvider):
    """ Model class that will retrieve the information from TiltSeries,
        Tomogram, SetOfTomograms, SetOfTiltSeries and  prepare the
        columns/rows models required by the TreeDialog GUI.
    """
    COL_TS = 'Tilt series'
    COL_INFO = 'Info'
    COL_STATUS = 'Status'
    COL_PREALIGNED = 'Prealigned'
    COL_ALIGNED = 'Aligned'
    COL_COOR3D = 'Coordinates 3D'
    COL_LANDMODEL_NO_GAPS = 'Landmark models no gaps'
    COL_RECONST_TOMOGRAM = 'Full tomograms'
    COL_PREPROCESS_RECONST_TOMOGRAM = 'Postprocess tomograms'

    ORDER_DICT = {COL_TS: 'id'}

    def __init__(self, protocol, objs, isInteractive=False):
        self.title = 'TiltSeries display'
        if isinstance(objs, tomo.objects.SetOfTomograms):
            self.COL_TS = 'Tomograms'
            self.title = 'Tomograms display'
        elif isinstance(objs, tomo.objects.SetOfLandmarkModels):
            self.COL_TS = 'LandmarkModels'
            self.title = 'LandmarkModels display'
        self.protocol = protocol
        self.objs = objs
        self.isInteractive = isInteractive
        TreeProvider.__init__(self, sortingColumnName=self.COL_TS)
        self.selectedDict = {}
        self.mapper = protocol.mapper
        self.maxNum = 200

    def getObjects(self):
        # Retrieve all objects of type className
        objects = []

        orderBy = self.ORDER_DICT.get(self.getSortingColumnName(), 'id')
        direction = 'ASC' if self.isSortingAscending() else 'DESC'

        for obj in self.objs.iterItems(orderBy=orderBy, direction=direction):
            if isinstance(obj, tomo.objects.TiltSeries):
                item = obj.clone(ignoreAttrs=('_mapperPath',))
            else:
                item = obj.clone()
            item._allowsSelection = True
            item._parentObject = None
            objects.append(item)

        return objects

    def _sortObjects(self, objects):
        pass

    def objectKey(self, pobj):
        pass

    def getColumns(self):
        cols = [
            (self.COL_TS, 100),
            (self.COL_INFO, 350)]
        if self.isInteractive:
            protocolClass = self.getProtocolClass()
            if protocolClass == protClass.protImodEtomoClass:
                cols.append((self.COL_PREALIGNED, 80))
                cols.append((self.COL_ALIGNED, 70))
                cols.append((self.COL_COOR3D, 110))
                cols.append((self.COL_LANDMODEL_NO_GAPS, 190))
                cols.append((self.COL_RECONST_TOMOGRAM, 120))
                cols.append((self.COL_PREPROCESS_RECONST_TOMOGRAM, 180))
            else:
                cols.append((self.COL_STATUS, 80))

        return cols

    def isSelected(self, obj):
        """ Check if an object is selected or not. """
        return False

    @staticmethod
    def _getParentObject(pobj, default=None):
        return getattr(pobj, '_parentObject', default)

    def getObjectInfo(self, obj):
        itemId = obj.getTsId()
        if itemId is None:
            itemId = str(obj.getObjId())

        key = obj.getObjId()
        text = itemId
        values = [str(obj)]
        tags = ''
        if self.isInteractive:
            status = self.getObjStatus(obj, values)
            tags = (status,)

        opened = True

        item = {
            'key': key, 'text': text,
            'values': tuple(values),
            'open': opened,
            'selected': False,
            'parent': obj._parentObject,
            'tags': tags
        }
        return item

    def getProtocolClass(self):
        if issubclass(self.protocol.__class__, ProtImodEtomo):
            protocolClass = protClass.protImodEtomoClass
        else:
            protocolClass = protClass.protImodCTFEstimation
        return protocolClass

    def getImodEtomoColumnValues(self, obj, values):
        status = 'pending'
        for item in  self.protocol.inputSetOfTiltSeries.get():
            if item.getTsId() == obj.getTsId():
                """Prealigned tilt-series"""
                prealiFilePath = self.protocol.getFilePath(item,
                                                           extension=".preali")
                if os.path.exists(prealiFilePath):
                    values.append('Yes')
                    status = 'done'
                else:
                    values.append('No')

                """Aligned tilt-series"""
                aligFilePath = self.protocol.getFilePath(item, extension=".ali")
                if os.path.exists(aligFilePath):
                    values.append('Yes')
                    status = 'done'
                else:
                    values.append('No')

                coordFilePath = self.protocol.getFilePath(item, suffix='fid',
                                                          extension=".xyz")
                if os.path.exists(coordFilePath):
                    values.append('Yes')
                    status = 'done'
                else:
                    values.append('No')

                """Landmark models with no gaps"""
                if (os.path.exists(self.protocol.getFilePath(item, suffix="_nogaps",
                                                             extension=".fid")) and
                        os.path.exists(
                            self.protocol.getFilePath(item, extension=".resid"))):
                    values.append('Yes')
                    status = 'done'
                else:
                    values.append('No')

                """Full reconstructed tomogram"""
                reconstructTomoFilePath = self.protocol.getFilePath(item,
                                                                    suffix="_full",
                                                                    extension=".rec")
                if os.path.exists(reconstructTomoFilePath):
                    values.append('Yes')
                    status = 'done'
                else:
                    values.append('No')

                """Post-processed reconstructed tomogram"""
                posprocessedRecTomoFilePath = self.protocol.getFilePath(item,
                                                                        extension=".rec")
                if os.path.exists(posprocessedRecTomoFilePath):
                    values.append('Yes')
                    status = 'done'
                else:
                    values.append('No')
                break
        return status

    def getImodCTFEstimationColumnValues(self, obj, values):
        status = 'pending'
        for item in self.protocol.inputSetOfTiltSeries:
            if item.getTsId() == obj.getTsId():
                extraPrefix = self.protocol._getExtraPath(item.getTsId())
                defocusFilePath = os.path.join(extraPrefix, item.getFirstItem().parseFileName(extension=".defocus"))
                if os.path.exists(defocusFilePath):
                    values.append('DONE')
                    status = 'done'
                else:
                    values.append('TODO')
                break
        return status

    def getObjStatus(self, obj, values):
        protocolClass = self.getProtocolClass()

        if protocolClass == protClass.protImodEtomoClass:
            status = self.getImodEtomoColumnValues(obj, values)
        else:
            status = self.getImodCTFEstimationColumnValues(obj, values)

        return status

    def getObjectActions(self, obj):
        actions = []
        if not self.isInteractive:
            viewers = Domain.findViewers(obj.getClassName(),
                                         pwviewer.DESKTOP_TKINTER)
            for viewerClass in viewers:
                def createViewer(viewerClass, obj):
                    proj = self.protocol.getProject()
                    item = self.objs[obj.getObjId()]  # to load mapper
                    return lambda : viewerClass(project=proj).visualize(item)
                actions.append(('Open with %s' % viewerClass.__name__,
                                createViewer(viewerClass, obj)))
        return actions

    def configureTags(self, tree):
        tree.tag_configure("pending", foreground="red")
        tree.tag_configure("done", foreground="green")


class ImodListDialog(ListDialog):
    def __init__(self, parent, title, provider, displayAllButton=True,
                 createSetButton=False,
                 itemDoubleClick=False, **kwargs):
        self.displayAllButton = displayAllButton
        self.createSetButton = createSetButton
        self._itemDoubleClick = itemDoubleClick
        self.provider = provider
        ListDialog.__init__(self, parent, title, provider, message=None,
                            allowSelect=False,  cancelButton=True, **kwargs)

    def body(self, bodyFrame):
        bodyFrame.config()
        gui.configureWeigths(bodyFrame)
        dialogFrame = tk.Frame(bodyFrame)
        dialogFrame.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        dialogFrame.config()
        gui.configureWeigths(dialogFrame, row=1)
        self._createFilterBox(dialogFrame)
        self._col = 0
        if self.displayAllButton:
            self.displayAll = self._addButton(dialogFrame,
                                              'Display all at once',
                                              pwutils.Icon.ACTION_VISUALIZE,
                                              self._displayAll,
                                              sticky='ne',
                                              state=tk.NORMAL)
        if self.createSetButton:
            self.createSet = self._addButton(dialogFrame,
                                             'CTFTomo',
                                             pwutils.Icon.PLUS_CIRCLE,
                                             self._createOutput,
                                             sticky='ne',
                                             state=tk.NORMAL)
        self._createTree(dialogFrame)
        self.initial_focus = self.tree
        if self._itemDoubleClick:
            self.tree.itemDoubleClick = self.doubleClickOnItem

    def _addButton(self, frame, text, image, command, sticky='news',
                   state=tk.NORMAL):

        defaults = {'activebackground': gui.cfgButtonActiveBgColor,
                    'bg': gui.cfgButtonBgColor,
                    'fg': gui.cfgButtonFgColor,
                    'activeforeground': gui.cfgButtonActiveFgColor,
                    'compound': tk.LEFT}

        btn = tk.Button(frame, text=text, image=self.getImage(image),
                        cursor='hand2', state=state,
                        **defaults)
        btn.bind('<Button-1>', command)
        btn.grid(row=0, column=self._col, sticky=sticky,
                 padx=(0, 5), pady=5)
        self._col += 1
        return btn

    def _displayAll(self, e=None):
        set = self.provider.objs
        if isinstance(set, tomo.objects.SetOfTiltSeries):
            ImodSetView(set)
        elif isinstance(set, tomo.objects.SetOfLandmarkModels):
            ImodSetOfLandmarkModelsView(set)
        elif isinstance(set, tomo.objects.SetOfTomograms):
            ImodSetOfTomogramsView(set)

    def _createOutput(self, e=None):
        self.provider.protocol.createOutput()

    def doubleClickOnItem(self, e=None):
        ts = e
        protocol = self.provider.protocol
        self.proc = threading.Thread(target=protocol.runAllSteps,
                                     args=(ts,))
        self.proc.start()
        self.after(1000, self.refresh_gui)

    def refresh_gui(self):
        self.tree.update()
        if self.proc.isAlive():
            self.after(1000, self.refresh_gui)


class ImodSetView(pwviewer.CommandView):
    """ Wrapper to visualize different type of objects with the 3dmod """

    def __init__(self, set, **kwargs):
        fn = ""
        for item in set:
            # Remove :mrc if present
            fn += " " + item.getFirstItem().getFileName().split(':')[0]
        pwviewer.CommandView.__init__(self, "%s %s" % (Plugin.getImodCmd('3dmod'), fn))
        self.show()


class ImodSetOfTomogramsView(pwviewer.CommandView):
    """ Wrapper to visualize set of tomograms with 3dmod """

    def __init__(self, set, **kwargs):
        fn = " -s 0,0 "
        for item in set:
            fn += " " + item.getLocation()[1]
        pwviewer.CommandView.__init__(self, Plugin.getImodCmd('3dmod') + fn)
        self.show()


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
        self.show()


class ImodGenericViewer(pwviewer.View):
    """ This class implements a view using Tkinter ListDialog
    and the ImodTreeProvider.
    """
    def __init__(self, parent, protocol, objs, displayAllButton=True,
                 createSetButton=False, isInteractive=False,
                 itemDoubleClick=False, **kwargs):
        """
         Params:
            parent: Tkinter parent widget


        From kwargs:
                message: message tooltip to show when browsing.
                selected: the item that should be selected.
                validateSelectionCallback:
                    a callback function to validate selected items.
                allowSelect: if set to False, the 'Select' button will not
                    be shown.
                allowsEmptySelection: if set to True, it will not validate
                    that at least one element was selected.
        """
        self._tkParent = parent
        self._protocol = protocol
        self._provider = ImodGenericTreeProvider(protocol, objs, isInteractive)
        self.title = self._provider.title
        self.displayAllButton = displayAllButton
        self.createSetButton = createSetButton
        self.itemDoubleClick = itemDoubleClick

    def show(self):
        ImodListDialog(self._tkParent, self.title, self._provider,
                       displayAllButton=self.displayAllButton,
                       createSetButton=self.createSetButton,
                       itemDoubleClick=self.itemDoubleClick)


