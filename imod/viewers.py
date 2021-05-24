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
import tkinter
from tkinter import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import pyworkflow.viewer as pwviewer
from pyworkflow.gui import *
from pyworkflow.gui.tree import TreeProvider
from pyworkflow.gui.dialog import ListDialog, showInfo
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