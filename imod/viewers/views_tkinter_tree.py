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


class CTFSerieStates:
    UNCHECKED = 'unchecked'
    CHECKED = 'checked'
    ODD = 'odd'
    EVEN = 'even'
    FAILED = 'Failed'
    OK = 'Ok'


class CtfEstimationTreeProvider(TreeProvider, ttk.Treeview):
    """ Model class that will retrieve the information from TiltSeries and
    prepare the columns/rows models required by the TreeDialog GUI.
    """
    COL_CTF_SERIE = 'CTF Series'
    CRITERIA_1 = 'defocusUDeviation'
    CRITERIA_2 = 'defocusVDeviation'
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
                  COL_CTF_EST_FIT: '_fitQuality'}

    def __init__(self, master, protocol, outputSetOfCTFTomoSeries, **kw):
        ttk.Treeview.__init__(self, master, **kw)
        self.protocol = protocol
        self.ctfSeries = outputSetOfCTFTomoSeries
        TreeProvider.__init__(self, sortingColumnName=self.COL_CTF_SERIE)
        self.selectedDict = {}
        self.mapper = protocol.mapper
        self.maxNum = 200
        self._checkedItems = 0

    def getObjects(self):
        # Retrieve all objects of type className
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

    def getCTFSeries(self):
        return self.ctfSeries

    def _sortObjects(self, objects):
        pass

    def objectKey(self, pobj):
        pass

    def getColumns(self):
        cols = [
            (self.COL_CTF_SERIE, 200),
            (self.CRITERIA_1, 150),
            (self.CRITERIA_2, 150),
            (self.COL_CTF_EST_IX, 100),
            (self.COL_CTF_EST_DEFOCUS_U, 100),
            (self.COL_CTF_EST_DEFOCUS_V, 100),
            (self.COL_CTF_EST_DEFOCUS_RATIO, 100),
            (self.COL_CTF_EST_DEFOCUS_ANGLE, 100)
        ]
        return cols

    def isSelected(self, obj):
        """ Check if an object is selected or not. """
        return False

    @staticmethod
    def _getParentObject(pobj, default=None):
        return getattr(pobj, '_parentObject', default)

    def getObjectInfo(self, obj):
        if isinstance(obj, tomo.objects.CTFTomoSeries):
            key = obj.getTsId()
            text = obj.getTsId()
            values = [CTFSerieStates.OK if obj.getIsDefocusUDeviationInRange()
                      else CTFSerieStates.FAILED,
                      CTFSerieStates.OK if obj.getIsDefocusVDeviationInRange()
                      else CTFSerieStates.FAILED,
                      '', '', '']
            opened = False
            selected = obj.isEnabled()
        else:  # CTFTomo
            key = "%s.%s" % (obj._parentObject.getTsId(), str(obj.getObjId()))
            text = ''
            values = [CTFSerieStates.OK if obj.getIsDefocusUDeviationInRange()
                      else CTFSerieStates.FAILED,
                      CTFSerieStates.OK if obj.getIsDefocusVDeviationInRange()
                      else CTFSerieStates.FAILED,
                      obj.getIndex(), str(obj.getDefocusU()),
                      str(obj.getDefocusV()), str(obj.getDefocusRatio()),
                      str(obj.getDefocusAngle())]
            opened = False
            selected = False

        item = {
            'key': key, 'text': text,
            'values': tuple(values),
            'open': opened,
            'selected': selected,
            'parent': obj._parentObject
        }
        if isinstance(obj, tomo.objects.CTFTomoSeries):
            tags = CTFSerieStates.UNCHECKED
            if not (obj.getIsDefocusUDeviationInRange() and obj.getIsDefocusVDeviationInRange()):
                obj.setEnabled(True)
                tags = CTFSerieStates.CHECKED
                self._checkedItems += 1

            if obj.getObjId() % 2 == 0:
                item['tags'] = (tags, CTFSerieStates.ODD,)
            else:
                item['tags'] = (tags,  CTFSerieStates.EVEN)
        else:
            if obj.getObjId() % 2 == 0:
                item['tags'] = (CTFSerieStates.ODD,)
            else:
                item['tags'] = (CTFSerieStates.EVEN,)
        return item


class CTFEstimationTree(BoundTree):
    def __init__(self, master, provider,  **opts):
        BoundTree.__init__(self, master, provider, frame=True, **opts)
        self.selectedItem = None
        self._checkedItems = provider._checkedItems

    def check_item(self, item):
        """ check the box of item and change the state of the boxes of item's
            ancestors accordingly """
        tags = CTFSerieStates.EVEN
        if CTFSerieStates.ODD in self.item(item, 'tags'):
            tags = CTFSerieStates.ODD

        if CTFSerieStates.UNCHECKED in self.item(item, 'tags'):
            self.item(item, tags=(CTFSerieStates.CHECKED, tags,))
            self._checkedItems += 1
            self.getSelectedObj().setEnabled(False)
            self.item(item)['selected'] = True
        else:
            self.item(item, tags=(CTFSerieStates.UNCHECKED, tags,))
            self.getSelectedObj().setEnabled(True)
            self._checkedItems -= 1
            self.item(item)['selected'] = False

    def _onClick(self, event=None):
        self._unpostMenu()
        x, y, widget = event.x, event.y, event.widget
        elem = widget.identify("element", x, y)
        self.selectedItem = self.identify_row(y)
        self.focus(self.selectedItem)
        if "image" in elem:
            self.check_item(self.selectedItem)

    def getSelectedItem(self):
        return self.selectedItem

    def getSelectedObj(self):
        obj = None
        if self.selectedItem:
            selected = self.getFirst()
            if selected is None:
                obj = self._objDict[self.selectedItem]
            else:
                obj = self._objDict[selected]
        return obj


class CtfEstimationListDialog(ListDialog):
    def __init__(self, parent, title, provider, protocol, inputTS, **kwargs):
        self._project = protocol.getProject()
        self._protocol = protocol
        self._inputSetOfTiltSeries = inputTS
        self._checkedItems = provider._checkedItems
        self.suffix = 0
        ListDialog.__init__(self, parent, title, provider, message=None,
                            allowSelect=False, **kwargs)

    def body(self, bodyFrame):
        bodyFrame.config()
        self._col = 1
        self._fillCTFEstimationGUI(bodyFrame)

    def _addButton(self, frame, text, image, command, sticky='news', state=tk.NORMAL):
        btn = tk.Label(frame, text=text, image=self.getImage(image),
                        compound=tk.LEFT, cursor='hand2', state=state)
        btn.bind('<Button-1>', command)
        btn.grid(row=0, column=self._col, sticky=sticky,
                 padx=(0, 5), pady=5)
        self._col += 1
        return btn

    def _fillCTFEstimationGUI(self, bodyFrame):
        # Create a top panel to put the filter box and bottoms
        topPanel = tk.Frame(bodyFrame)
        topPanel.grid(row=0, column=0, padx=0, pady=0, sticky='news')
        self._createTopPanel(topPanel)

        # Create a bottom panel to put the tree and the plotter
        bottomPanel = tk.Frame(bodyFrame)
        bottomPanel.grid(row=1, column=0, padx=0, pady=0, sticky='news')
        self._createBottomPanel(bottomPanel)

    def _createTopPanel(self, topPanel):
        self._createFilterBox(topPanel)

        topRigthPanel = tk.Frame(topPanel)
        topRigthPanel.grid(row=0, column=1, padx=0, pady=0, sticky='news')
        self._createRecalculateBottom(topRigthPanel)
        self._createViewerHelp(topRigthPanel)

    def _createRecalculateBottom(self, topRigthPanel):

        state = tk.NORMAL
        if self._checkedItems or self._checkedItems == len(self.provider.getCTFSeries()):
            state = tk.DISABLED
        self.generateSubsetButton = self._addButton(topRigthPanel,
                                                    'Generate subsets',
                                                    pwutils.Icon.PROCESSING,
                                                    self._actionCreateSets,
                                                    sticky='ne',
                                                    state=state)

    def _createViewerHelp(self, topRigthPanel):
        self._addButton(topRigthPanel, pwutils.Message.LABEL_HELP,
                        pwutils.Icon.ACTION_HELP, self._showHelp, sticky='ne')

    def _actionCreateSets(self, event=None):

        if self.generateSubsetButton['state'] == tk.NORMAL:
            protocol = self.provider.protocol
            ctfSeries = self.provider.getCTFSeries()
            suffix = self._getSuffix(protocol)
            goodCTFName = 'goodCtf%s' % suffix
            badCTFName = 'badCtf%s' % suffix

            outputSetOfgoodCTFTomoSeries = ctfSeries.createCopy(protocol._getPath(),
                                                                prefix=goodCTFName,
                                                                copyInfo=True)
            outputSetOfbadCTFTomoSeries = ctfSeries.createCopy(protocol._getPath(),
                                                               prefix=badCTFName,
                                                               copyInfo=True)
            for ctfSerie in ctfSeries:
                ctfSerieClon = ctfSerie.clone()
                if CTFSerieStates.UNCHECKED in self.tree.item(ctfSerie.getTsId(),
                                                              'tags'):
                    # Adding the ctfSerie to the good set of ctfTomoSeries
                    outputSetOfgoodCTFTomoSeries.append(ctfSerieClon)
                    outputSetOfgoodCTFTomoSeries.setSetOfTiltSeries(self._inputSetOfTiltSeries)

                else:
                    # Adding the ctfSerie to the bad set of ctfTomoSeries
                    outputSetOfbadCTFTomoSeries.append(ctfSerieClon)
                    outputSetOfbadCTFTomoSeries.setSetOfTiltSeries(self._inputSetOfTiltSeries)

                for item in ctfSerie.iterItems():
                    ctfEstItem = item.clone()
                    ctfSerieClon.append(ctfEstItem)

            outputgoodCTFSetName = 'goodSetOfCTFTomoSeries%s' % suffix
            outputbadCTFSetName = 'badSetOfCTFTomoSeries%s' % suffix

            if len(outputSetOfgoodCTFTomoSeries) > 0:
                protocol._defineOutputs(**{outputgoodCTFSetName: outputSetOfgoodCTFTomoSeries})

            if len(outputSetOfbadCTFTomoSeries) > 0:
                protocol._defineOutputs(**{outputbadCTFSetName: outputSetOfbadCTFTomoSeries})

            protocol._store()
            self.cancel()

    def _showHelp(self, event=None):
        showInfo('CTFTomoSeries viewer Help',
                 'This viewer calculates the standard deviation with respect '
                 'to the mean of the defocusU and defocusV values. If the '
                 'values of the images are not in the range they are marked '
                 'as Failed and therefore the CTFTomoSerie is marked as '
                 'Failed as well. '
                 'On the other hand, the viewer allows you to create two '
                 'subsets of CTFTomoSeries which are classified as good '
                 'and bad respectively. '
                 'Note: The ctfseries that are checked are the ones that '
                 'represent the bad ctfseries', self.parent)

    def _getSuffix(self, protocol):
        """
        Return the number of the last output in order to complete the new
        output with a suffix
        """
        if self.suffix == 0:
            self.suffix = sum(1 for _ in protocol.iterOutputAttributes()) + 1
        else:
            self.suffix = 1
        return self.suffix

    def _createBottomPanel(self, bottomPanel):
        self._createCTFEstimationGUI(bottomPanel)
        self.initial_focus = self.tree

    def _createCTFEstimationGUI(self, bottomPanel):
        # Create a division Paned
        pw = tk.PanedWindow(bottomPanel, orient=tk.HORIZONTAL)
        # Create a left panel to put the tree
        bottomleftPanel = tk.Frame(pw)
        bottomleftPanel.grid(row=0, column=0, padx=0, pady=0, sticky='news')
        self._createTree(bottomleftPanel)
        pw.add(bottomleftPanel)
        # Panel to put the plotter
        self.bottomRightPanel = ttk.Frame(pw)
        self.bottomRightPanel.grid(row=0, column=1, padx=0, pady=0, sticky='news')
        self._createPloter(self.bottomRightPanel)
        pw.add(self.bottomRightPanel)
        pw.pack(fill=BOTH, expand=True)
        # This method is used to show sash
        pw.configure(sashrelief=RAISED)

    def _createTree(self, parent):

        gui.configureWeigths(parent)

        self.tree = CTFEstimationTree(parent, self.provider,
                                      selectmode=self._selectmode)
        item = self.tree.identify_row(0)
        self.tree.selection_set(item)
        self.tree.focus(item)
        self.tree.selectedItem = item
        self.im_checked = gui.getImage(Icon.CHECKED)
        self.im_unchecked = gui.getImage(Icon.UNCHECKED)
        self.tree.tag_configure(CTFSerieStates.UNCHECKED,
                                image=self.im_unchecked)
        self.tree.tag_configure(CTFSerieStates.CHECKED,
                                image=self.im_checked)
        self.tree.tag_configure(CTFSerieStates.EVEN, background='#F2F2F2',
                                foreground='black')
        self.tree.tag_configure(CTFSerieStates.ODD, background='#E6E6E6',
                                foreground='black')
        self.tree.bind("<Button-1>", self._createPloter, True)

    def _createPloter(self, event):
        obj = self.tree.getSelectedObj()
        self._checkedItems = self.tree._checkedItems
        if self._checkedItems and self._checkedItems != len(self.provider.getCTFSeries()):
            self.generateSubsetButton['state'] = tk.NORMAL
        else:
            self.generateSubsetButton['state'] = tk.DISABLED
        if obj is not None:
            plotterPanel = tk.Frame(self.bottomRightPanel)
            defocusUList = []
            defocusVList = []
            itemSelected = self.tree.selectedItem

            if self.tree.parent(self.tree.selectedItem):
                itemSelected = self.tree.parent(self.tree.selectedItem)

            for ctfSerie in self.provider.getCTFSeries():
                if ctfSerie.getTsId() == itemSelected:
                    for item in ctfSerie.iterItems(orderBy='id'):
                        defocusUList.append(item.getDefocusU())
                        defocusVList.append(item.getDefocusV())

                    fig = Figure(figsize=(7, 7), dpi=100)
                    defocusPlot = fig.add_subplot(111)
                    defocusPlot.plot(defocusUList, marker='o', label='DefocusU')
                    defocusPlot.plot(defocusVList, marker='o', label='DefocusV')
                    fig.legend()
                    canvas = FigureCanvasTkAgg(fig, master=plotterPanel)
                    canvas.draw()
                    canvas.get_tk_widget().pack(fill=tkinter.BOTH, expand=0)
                    plotterPanel.grid(row=0, column=1, sticky='news')
                    break


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
            cols.append((self.COL_PREALIGNED, 80))
            cols.append((self.COL_ALIGNED, 70))
            cols.append((self.COL_COOR3D, 110))
            cols.append((self.COL_LANDMODEL_NO_GAPS, 190))
            cols.append((self.COL_RECONST_TOMOGRAM, 120))
            cols.append((self.COL_PREPROCESS_RECONST_TOMOGRAM, 180))

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

    def getObjStatus(self, obj, values):
        status = 'pending'
        for item in self.protocol.inputSetOfTiltSeries.get():
            if item.getTsId() == obj.getTsId():
                """Prealigned tilt-series"""
                prealiFilePath = self.protocol.getFilePath(item, extension=".preali")
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
                if (os.path.exists(self.protocol.getFilePath(item, suffix="_nogaps", extension=".fid")) and
                        os.path.exists(self.protocol.getFilePath(item, extension=".resid"))):
                    values.append('Yes')
                    status = 'done'
                else:
                    values.append('No')

                """Full reconstructed tomogram"""
                reconstructTomoFilePath = self.protocol.getFilePath(item, suffix="_full", extension=".rec")
                if os.path.exists(reconstructTomoFilePath):
                    values.append('Yes')
                    status = 'done'
                else:
                    values.append('No')

                """Post-processed reconstructed tomogram"""
                posprocessedRecTomoFilePath = self.protocol.getFilePath(item, extension=".rec")
                if os.path.exists(posprocessedRecTomoFilePath):
                    values.append('Yes')
                    status = 'done'
                else:
                    values.append('No')
                break
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
                 itemDoubleClick=False, **kwargs):
        self.displayAllButton = displayAllButton
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
        self._createTree(dialogFrame)
        self.initial_focus = self.tree
        if self._itemDoubleClick:
            self.tree.itemDoubleClick = self.doubleClickOnItem

    def _addButton(self, frame, text, image, command, sticky='news', state=tk.NORMAL):
        btn = tk.Button(frame, text=text, image=self.getImage(image),
                        compound=tk.LEFT, cursor='hand2', state=state)
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

    def doubleClickOnItem(self, e=None):
        ts = e
        protocol = self.provider.protocol
        if issubclass(protocol.__class__, ProtImodEtomo):
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
                 isInteractive=False, itemDoubleClick=False, **kwargs):
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
        self.itemDoubleClick = itemDoubleClick

    def show(self):
        ImodListDialog(self._tkParent, self.title, self._provider,
                       displayAllButton=self.displayAllButton,
                       itemDoubleClick=self.itemDoubleClick)


