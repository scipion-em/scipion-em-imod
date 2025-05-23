# *****************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import os.path
import threading

from pyworkflow.gui import *
from pyworkflow.gui.tree import TreeProvider
from pyworkflow.gui.dialog import ListDialog
import pyworkflow.viewer as pwviewer
import tomo.objects

from imod import Plugin
from imod.constants import MRC_EXT, XYZ_EXT, FID_EXT, RESID_EXT, DEFOCUS_EXT
from imod.protocols import ProtImodEtomo


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
    COL_LANDMODEL_NO_GAPS = 'Fiducial models w/o gaps'
    COL_RECONST_TOMOGRAM = 'Raw tomograms'
    COL_PREPROCESS_RECONST_TOMOGRAM = 'Post-processed tomograms'

    ORDER_DICT = {COL_TS: 'id'}

    def __init__(self, protocol, objs, isInteractive=False):
        self.title = 'Imod set viewer'
        if isinstance(objs, tomo.objects.SetOfTomograms):
            self.COL_TS = 'Tomograms'
            self.title = 'Tomograms display'
        elif isinstance(objs, tomo.objects.SetOfLandmarkModels):
            self.COL_TS = 'Landmark models'
            self.title = 'Landmark models display'
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
                item = obj.clone(ignoreAttrs=['_mapperPath'])
            elif isinstance(obj, tomo.objects.LandmarkModel):
                self.objs.completeLandmarkModel(obj)
                item = obj.clone()
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
            (self.COL_TS, 200),
            (self.COL_INFO, 400)]
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
        for item in self.protocol.inputSetOfTiltSeries.get():
            if item.getTsId() == obj.getTsId():
                tsId = item.getTsId()
                """Prealigned tilt-series"""
                prealiFilePath = self.protocol.getExtraOutFile(tsId, suffix="preali",
                                                               ext=MRC_EXT)
                if os.path.exists(prealiFilePath):
                    values.append('Yes')
                    status = 'done'
                else:
                    values.append('No')

                """Aligned tilt-series"""
                aligFilePath = self.protocol.getExtraOutFile(tsId, suffix="ali",
                                                             ext=MRC_EXT)
                if os.path.exists(aligFilePath):
                    values.append('Yes')
                    status = 'done'
                else:
                    values.append('No')

                coordFilePath = self.protocol.getExtraOutFile(tsId, suffix='fid',
                                                              ext=XYZ_EXT)
                coordFilePath = coordFilePath.replace("_fid", "fid")  # due to etomo bug
                if os.path.exists(coordFilePath):
                    values.append('Yes')
                    status = 'done'
                else:
                    values.append('No')

                """Landmark models with no gaps"""
                if (os.path.exists(self.protocol.getExtraOutFile(tsId, suffix="nogaps", ext=FID_EXT)) and
                        os.path.exists(self.protocol.getExtraOutFile(tsId, ext=RESID_EXT))):
                    values.append('Yes')
                    status = 'done'
                else:
                    values.append('No')

                """Full reconstructed tomogram"""
                reconstructTomoFilePath = self.protocol.getExtraOutFile(tsId, suffix="full_rec",
                                                                        ext=MRC_EXT)
                if os.path.exists(reconstructTomoFilePath):
                    values.append('Yes')
                    status = 'done'
                else:
                    values.append('No')

                """Post-processed reconstructed tomogram"""
                posprocessedRecTomoFilePath = self.protocol.getExtraOutFile(tsId, suffix="rec",
                                                                            ext=MRC_EXT)
                if os.path.exists(posprocessedRecTomoFilePath):
                    values.append('Yes')
                    status = 'done'
                else:
                    values.append('No')
                break
        return status

    def getImodCTFEstimationColumnValues(self, obj, values):
        status = 'pending'
        for item in self.protocol.inputSetOfTiltSeries.get():
            if item.getTsId() == obj.getTsId():
                defocusFilePath = self.protocol.getExtraOutFile(item.getTsId(), ext=DEFOCUS_EXT)
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

    def configureTags(self, tree):
        tree.tag_configure("pending", foreground="red")
        tree.tag_configure("done", foreground="green")


class ImodListDialog(ListDialog):
    def __init__(self, parent, title, provider, createSetButton=False,
                 itemDoubleClick=True, lockGui=False, **kwargs):
        self.createSetButton = createSetButton
        self._itemDoubleClick = itemDoubleClick
        self.provider = provider
        self.binningVar = None
        ListDialog.__init__(self, parent, title, provider,
                            message=None,
                            allowSelect=False,
                            lockGui=lockGui,
                            cancelButton=True,
                            **kwargs)

    def body(self, bodyFrame):
        bodyFrame.config()
        gui.configureWeigths(bodyFrame)
        dialogFrame = tk.Frame(bodyFrame)
        dialogFrame.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        dialogFrame.config()
        gui.configureWeigths(dialogFrame, row=1)
        self._createFilterBox(dialogFrame)
        self._col = 0
        if self.createSetButton:
            self.createSet = self._addButton(dialogFrame,
                                             'CTFTomo',
                                             pwutils.Icon.PROCESSING,
                                             self._createOutput,
                                             sticky='ne',
                                             state=tk.NORMAL)
        else:
            self._addBinningBox()
        self._createTree(dialogFrame)
        self.initial_focus = self.tree
        if self._itemDoubleClick:
            if self.provider.isInteractive:  # etomo, ctf estimation
                self.tree.itemDoubleClick = self.runProtocolSteps
            else:
                self.tree.itemDoubleClick = self.openImodViewer

    def _addBinningBox(self):
        self.binningVar = tk.StringVar(value=str(Plugin.getViewerBinning()))
        frame = self.searchBoxframe
        label = tk.Label(frame, text="Display binning")
        label.grid(row=0, column=2, sticky='nw')

        entry = tk.Entry(frame, bg=Config.SCIPION_BG_COLOR,
                         width=3, textvariable=self.binningVar,
                         font=gui.getDefaultFont())
        entry.grid(row=0, column=3, sticky='news')

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

    def _createOutput(self, e=None):
        self.provider.protocol.createOutput()

    def openImodViewer(self, item=None):
        from imod.viewers import ImodObjectView
        prot = self.provider.protocol
        item = self.provider.objs[item.getObjId()]  # to load mapper

        ImodObjectView(item, protocol=prot, binning=self.binningVar.get()).show()

    def runProtocolSteps(self, e=None):
        ts = e
        protocol = self.provider.protocol
        self.proc = threading.Thread(target=protocol.runAllSteps,
                                     args=(ts,))
        self.proc.start()
        self.after(1000, self.refresh_gui)

    def refresh_gui(self):
        self.tree.update()
        if self.proc.is_alive():
            self.after(1000, self.refresh_gui)


class ImodGenericView(pwviewer.View):
    """ This class implements a view using Tkinter ListDialog
    and the ImodTreeProvider.
    """
    def __init__(self, parent, protocol, objs, createSetButton=False,
                 isInteractive=False, itemDoubleClick=True, **kwargs):
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
        self.createSetButton = createSetButton
        self.itemDoubleClick = itemDoubleClick

    def show(self):
        ImodListDialog(self._tkParent, self.title, self._provider,
                       createSetButton=self.createSetButton,
                       itemDoubleClick=self.itemDoubleClick)
