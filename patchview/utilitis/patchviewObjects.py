# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 11:41:59 2021

@author: MHu
"""
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import *
# from PyQt5.QtWidgets import QSpinBox, QDoubleSpinBox, QGridLayout, QGroupBox, QFormLayout, QTabWidget, QTableWidget, QHeaderView
# from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout, QSizePolicy, QScrollArea, QFrame, QSplitter, QToolButton, QMenu, QAction
import pyqtgraph as pg
from pyqtgraph import QtGui
from pyqtgraph import GraphicsLayoutWidget
from pyqtgraph.parametertree import Parameter, ParameterTree
import os, glob
import numpy as np
import pandas as pd
from pyqtgraph.widgets import MatplotlibWidget
import pdb
translate = QtCore.QCoreApplication.translate


class FileModel(pg.QtWidgets.QFileSystemModel):
    """
    Class for file system model
    """

    def __init__(self, root_path):
        super(FileModel, self).__init__()

        # hide system files
        self.setFilter(
            QtCore.QDir.AllDirs | QtCore.QDir.NoDotAndDotDot | QtCore.QDir.AllEntries
        )

        # filter out non dats and disable showing
        self.setNameFilters(["*.dat", "*.abf","*.asc","*.ASC"])
        self.setNameFilterDisables(False)

        # set root
        self.setRootPath(root_path)


class sliceSelectionDialog(pg.QtWidgets.QWidget):
    def __init__(self, parent=None, items=None):
        super(sliceSelectionDialog, self).__init__(parent)
        self.items = items
        layout = QtGui.QFormLayout()
        self.btn = QtGui.QPushButton("Choose from list")
        self.btn.clicked.connect(self.getItem)  ## list all items

        self.sliceName = QtGui.QLineEdit()
        layout.addRow(self.btn, self.sliceName)
        self.btn1 = QtGui.QPushButton("Edit slice name:")
        self.btn1.clicked.connect(self.gettext)

        self.setLayout(layout)
        self.setWindowTitle("Input Dialog demo")

    def getItem(self):
        item, ok = QtGui.QInputDialog.getItem(
            self, "select input dialog", "list of .dat files", self.items, 0, False
        )

        if ok and item:
            self.le.setText(item)
            self.sliceName.setText(item)

    def gettext(self):
        text, ok = QtGui.QInputDialog.getText(
            self, "Slice to import", "Edit slice name:"
        )
        if ok:
            self.sliceName.setText(str(text))


class FileView(pg.QtWidgets.QTreeView):
    """
    Class for view of file system model
    """

    def __init__(self, parent, root_path):
        super(FileView, self).__init__()
        self.frame = parent

        # set model
        self.model = FileModel(root_path)
        self.setModel(self.model)

        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.clicked.connect(self.clickedItem)
        self.myMenu = QtWidgets.QMenu("Menu", self)
        self.addFileAction = pg.QtGui.QAction("Import .dat/.abf files for this slice")
        self.addFileAction.triggered.connect(self.importAllDats_clicked)
        self.myMenu.addAction(self.addFileAction)
        self.currentPath = root_path
        # set root
        self.setRootIndex(self.model.index(root_path))

        # hide unnecessary columns
        self.hideColumn(1)
        self.hideColumn(2)
        self.hideColumn(3)

        # bind event
        self.doubleClicked.connect(self.on_double_click)
        self.customContextMenuRequested.connect(self.openContextMenu)

    def importAllDats_clicked(self):
        """We load all .dat files with current folder.
        And we sort these dat according to their protocols
        And sorted files are presented into each tab!
        """
        self.frame.clearAllTrees()
        datFiles = glob.glob(os.path.join(self.currentPath, "*.dat"))
        print(f"{len(datFiles)} .dat files found")
        if len(datFiles) >= 1:
            Items = []
            for d in datFiles:
                self.currentPath, dat0 = os.path.split(d)
                Items.append(dat0[:-4])
            datFiles = Items
            fileNameLen = [len(t) for t in datFiles]
            longestName = np.argmax(fileNameLen)
            shortestName = np.argmin(fileNameLen)
            ## make sure we really got two different files
            if shortestName == longestName:
                if shortestName != 0:
                    shortestName = shortestName - 1
                else:
                    if len(datFiles) == 1:
                        shortestName = 0
                    else:
                        shortestName = shortestName + 1
            longestName = datFiles[longestName]
            shortestName = datFiles[shortestName]
            fileNamePrefix = []  ## prefix for slice name
            for idx, c in enumerate(longestName):  ## common characters between files
                if idx < len(shortestName) and c == shortestName[idx]:
                    fileNamePrefix.append(c)
                else:
                    break
            sliceNameStart = -1
            for jdex, c in enumerate(fileNamePrefix):
                if c == "s" or c == "S":  ## reach the begining of slice name
                    sliceNameStart = jdex
                if sliceNameStart > 0 and c.isdigit():
                    sliceNameStop = jdex
                    break
            else:
                sliceNameStop = sliceNameStart + 1  ## no digit attached

            sliceName = "".join(fileNamePrefix[:sliceNameStop])
            sliceNumber = []
            for f in datFiles:
                sliceNumber_ = []
                for c in f[sliceNameStop:]:
                    if c.isdigit():
                        sliceNumber_.append(c)
                    else:
                        break
                sliceNumber_ = "".join(sliceNumber_)
                if sliceNumber_ != "" and (sliceNumber_ not in sliceNumber):
                    sliceNumber.append(sliceNumber_)

            items = ["All slices"]
            sliceNumber.sort()
            for sliceN in sliceNumber:
                items.append(sliceName + sliceN)
            print("Unique slices found:")
            print(items)
            sliceName, ok = QtWidgets.QInputDialog.getItem(
                self, "Input Dialog", "choose slice:", items, 0, False
            )
            if ok and sliceName:
                print(sliceName, " selected")
                self.loadSliceDats(self.currentPath, datFiles, sliceName)

            else:
                print("no slice file loaded")

    def loadSliceDats(self, filePath, datFiles, sliceName):
        if sliceName == "All slices":
            sliceFiles = datFiles  ## load all!
        else:
            n = len(sliceName)
            sliceFiles = []

            for f in datFiles:
                if f[:n] == sliceName:
                    sliceFiles.append(f)  ## gather files that belong to this slice

        for f in sliceFiles:
            fullpath = os.path.join(filePath, f)
            self.addPulseTree(fullpath)

    def addPulseTree(self, filePath):
        self.frame.update_sortedTrees(filePath)

    def clickedItem(self, index):
        model_index = self.model.index(index.row(), 0, index.parent())

        # get file from index
        file_path = os.path.abspath(self.model.filePath(model_index))
        self.frame.root, f = os.path.split(file_path)
        #        os.chdir(self.frame.root)
        # check extension
        _, ext = os.path.splitext(file_path)
        if ext == "":
            # print("this is a folder")
            self.currentPath = file_path
        elif ext == ".dat":
            self.currentFile = file_path

    def openContextMenu(self, index):
        self.myMenu.exec_(QtGui.QCursor.pos())

    @QtCore.pyqtSlot(QtCore.QModelIndex)
    def on_double_click(self, index):
        """
        Event for changing .dat file
        :param index:
        :return:
        """
        model_index = self.model.index(index.row(), 0, index.parent())

        # get file from index
        file_path = os.path.abspath(self.model.filePath(model_index))
        self.frame.root, f = os.path.split(file_path)
        #        os.chdir(self.frame.root)
        # check extension
        _, ext = os.path.splitext(file_path)
        if ext in ['.asc','.ASC']:
            self.frame.prepareTree(file_path)
        else:
            self.frame.clearAllTrees()
            self.frame.update_pul(file_path, dat_index=None, ext=ext)


class SelectionView(pg.QtWidgets.QTreeWidget):
    """Class for viewing selected .dat series"""

    def __init__(self, parent):
        super(SelectionView, self).__init__()
        self.frame = parent
        self.seriesTypes = {}
        self.seriesGroups = []
        self.current_dat = None
        self.current_datFileName = []
        self.bundleGroup = {}
        self.currentPath = []
        self.current_index = []
        self.seriesNode = []
        self.dat_files = (
            {}
        )  ## this dictionary stores selected seiries, using dat file name as key.
        self.stimChannelID = (
            {}
        )  ## this dictionary stores the stimChan ID in  selected seiries, using dat file name as key.
        self.seleted_index = []
        self.bundle = []
        self.indices = None
        self.sweepCount = 0
        self.traceCount = 0

        # self.setHeaderLabels(['Node', 'Label','dat', 'dat_Index'])
        self.setHeaderLabels(["Cell", "File", "Level", "Index"])
        self.setColumnWidth(0, 200)
        self.itemSelectionChanged.connect(self.on_selection_changed2)


    def openContextMenu(self, index):
        self.myMenu.exec_(QtGui.QCursor.pos())

    def update_treeSeries(self, root_item, index):
        """
            New recursive function. Skip the first two levels

        Parameters
        ----------
        root_item : TYPE
            DESCRIPTION.
        index : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        root = self.bundle.pul[self.current_index[0]][self.current_index[1]]
        node = root
        for i in index:
            node = node[i]
        node_type = node.__class__.__name__
        if node_type.endswith("Record"):
            node_type = node_type[:-6]
        node_label = node.Label
        node_datIndex = str(self.current_index[0]) + " " + str(self.current_index[1])
        if index == []:
            node_dat = self.current_datFileName
        else:
            node_dat = ""
            extraIndex = [str(j) for j in index]
            node_datIndex = node_datIndex + " " + " ".join(extraIndex)
        try:
            if node_type[:2] == "V9":
                node_type += str(getattr(node, node_type[3:] + "Count"))
            else:
                node_type += str(getattr(node, node_type + "Count"))
        except AttributeError:
            pass

        # self.setHeaderLabels(['Node', 'Label','dat', 'dat_Index'])
        # self.setHeaderLabels(['Cell label', 'File','Series', 'Index'])
        if len(index) == 1:
            item = pg.QtGui.QTreeWidgetItem(
                [node_type, node_dat, node_label, node_datIndex]
            )
        else:
            item = pg.QtGui.QTreeWidgetItem(
                [node_label, node_dat, node_type, node_datIndex]
            )

        item.setCheckState(0, QtCore.Qt.Checked)
        item.node = node
        item.index = index
        root_item.addChild(item)

        for i in range(len(node.children)):
            self.update_treeSeries(item, index + [i])

    def update_tree_recursive(self, root_item, index):
        """
        Parameters
        ----------
        root_item : TYPE
            DESCRIPTION.
        index : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        """Recursively read tree information from the bundle's embedded .pul file
        and add items into the GUI tree to allow browsing.
        """
        root = self.bundle.pul
        node = root
        for i in index:
            node = node[i]
        node_type = node.__class__.__name__
        if node_type.endswith("Record"):
            node_type = node_type[:-6]
        try:
            if node_type[:2] == "V9":
                node_type += str(getattr(node, node_type[3:] + "Count"))
            else:
                node_type += str(getattr(node, node_type + "Count"))
        except AttributeError:
            pass

        try:
            node_label = node.Label
            if (
                len(index) == 3
                and self.frame.parameters["Protocols"]["This serie"]["Type"]
                == "Spontaneous"
            ):
                node_label = "Channel name"
                node_type = "Concatenated sweeps"
        except AttributeError:
            if node_type[:5] == "Pulse":
                node_label = self.current_dat
            else:
                node_label = ""
        if node_type[:2] == "V9":
            item = pg.QtGui.QTreeWidgetItem([node_type[3:], node_label])
        else:
            item = pg.QtGui.QTreeWidgetItem([node_type, node_label])
        item.setCheckState(0, QtCore.Qt.Checked)
        if len(index) < 2:
            item.setExpanded(True)
        if len(index) == 2:  ## if this series in the selection list!
            self.seriesNode.append(node)
            if index in self.current_index or (
                index == self.current_index
            ):  ## only add this series

                root_item.addChild(item)
                item.node = node
                item.index = index
                for i in range(len(node.children)):
                    self.update_tree_recursive(item, index + [i])
        else:
            root_item.addChild(item)
            item.node = node
            item.index = index
            for i in range(len(node.children)):
                self.update_tree_recursive(item, index + [i])

    #    @QtCore.pyqtSlot(QtCore.QModelIndex)

    def getAllTrees(self):
        root = self.invisibleRootItem  ## get the root of this tree!
        child_count = root.childCount()
        for i in range(child_count):
            item = root.child(i)
            url = print(item.text(1))  # text at first (0) column

    def on_selection_changed2(self):

        selected = self.selectedItems()[0]  # self.seleted_index
        if selected == []:
            return
        if selected.text(1) != "":
            currentDatFile = selected.text(1) + ".dat"
            dat_index = [int(j) for j in selected.text(3).split(" ")]
            currentDatFilePath = os.path.join(self.currentPath, currentDatFile)
            print(currentDatFilePath)
            self.frame.update_pul(currentDatFilePath, dat_index=dat_index, ext=".dat")

    def on_selection_changed(self):
        """
        Event for browsing through wave traces
        :return:
        """
        selected = self.selectedItems()[0]  # self.seleted_index
        if selected == []:
            return
        currentTreeRoot = selected.parent()
        currentDatFile = []
        print(selected.index)
        # pdb.set_trace()
        for j in range(len(selected.index)):

            if currentTreeRoot.text(0) == "Pulsed":
                print(currentTreeRoot.text(1))
                currentDatFile = currentTreeRoot.text(1)
                break
            else:
                currentTreeRoot = currentTreeRoot.parent()
        # self.frame.update_selectView(currentDatFile, self.indices)
        if currentDatFile != []:
            self.frame.update_pul(currentDatFile)

            indices = []
            ### for muliple selection
            index = selected.index
            self.frame.checkSeriesType(selected)
            if len(index) == 4:  ## grab series level of data
                indices.append(index)
            self.indices = indices
            self.frame.update_trace_plot_selectionTree()


class datTree(object):
    def __init__(self, node_type=None, node_label=""):
        self.item = []
        self.node = []
        self.node_type = node_type
        self.node_label = node_label

    def addChild(self, item):
        self.item.append(item)


class cellROI(pg.EllipseROI):
    def __init__(self, parent, mainWindow, pos, size, label, **args):
        super(cellROI, self).__init__(pos, size, **args)
        self.parent = parent
        self.mainWindow = mainWindow  ## main window
        self.label = label
        self.pos = pos
        self.addTranslateHandle(
            pos=(0.5, 0.5), name=label
        )  ## so we could have name for this ROI
        self.sigRegionChangeFinished.connect(self.update)
        self.sigRemoveRequested.connect(self.remove)
        self.TextLabel = None

    def update(self):
        if hasattr(self, "handles"):
            pass  # print(self.getSceneHandlePositions())

        if hasattr(self, "TextLabel"):
            if self.TextLabel:
                pos = self.getSceneHandlePositions()[-1][1]
                # idx = self.getSceneHandlePositions()[-1][0]
                self.pos = self.parent.view.mapToView(pos)
                self.mainWindow.updateROIlabel(self.TextLabel, self.pos, self.label)
                # self.TextLabel.setPos(self.pos[0], self.pos[1]-10)
                # print('update text position')

    def remove(self):
        self.parent.removeItem(self)
        self.sigRegionChanged.disconnect()
        self.sigRegionChangeFinished.disconnect()
        self.mainWindow.roidata.pop(
            self.getSceneHandlePositions()[-1][0]
        )  ## key for this ROI??
        self.parent.removeItem(self.TextLabel)  ## remove associated text label
        self.mainWindow.updateROI_table()


class PulView(pg.QtWidgets.QTreeWidget):
    """
    Class for viewing tree of pul file
    """

    def __init__(self, parent):
        super(PulView, self).__init__()
        self.frame = parent
        self.dat_file = None
        self.bundle = None
        self.abf = None  ## for .abf file
        self.abf_stimInfo = None
        self.pul = None  ## fir .dat file
        self.indices = None
        self.sweepCount = 0
        self.traceCount = 0
        self.seriesNode = []
        self.setHeaderLabels(["Node", "Label"])
        self.setColumnWidth(0, 200)
        # allow multi selection
        self.setSelectionMode(QAbstractItemView.ExtendedSelection)

        # bind event
        # self.itemSelectionChanged.connect(self.on_selection_changed)

        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.openContextMenu)
        # bind event
        self.itemSelectionChanged.connect(self.on_selection_changed)
        self.myMenu = QtWidgets.QMenu("Menu", self)

        self.fpAnAction = pg.QtGui.QAction("Firing pattern analysis")
        self.fpAnAction.triggered.connect(
            self.frame.extraEphyAction_clicked
        )  # self.add_clicked)
        self.myMenu.addAction(self.fpAnAction)

        self.eventAnAction = pg.QtGui.QAction("Spontaneous event analysis")
        self.eventAnAction.triggered.connect(
            self.frame.eventDetectionAction_clicked
        )  # self.add_clicked)
        self.myMenu.addAction(self.eventAnAction)

    def concatenatingSweeps_clicked(self):
        print(self.indices)
        self.add_clicked()

    def openContextMenu(self, index):
        if self.indices:
            if len(self.indices) != 0:
                self.myMenu.exec_(QtGui.QCursor.pos())

    def add_clicked(self):
        """
        Add selected tree to tree list
        """
        # get file from index
        file_path = self.dat_file
        self.frame.root, f = os.path.split(file_path)
        #        os.chdir(self.frame.root)
        # check extension
        _, ext = os.path.splitext(file_path)
        if ext == ".dat":
            # print(self.indices)
            self.frame.update_selectView(file_path, self.indices)

    def update_tree_recursive(self, root_item, index, dat_index=None):
        """Recursively read tree information from the HEKA bundle's embedded .pul file
        and add items into the GUI tree to allow browsing.
        """
        root = self.bundle.pul
        self.filetype = ".dat"
        node = root
        if node == None:
            return
        for i in index:
            node = node[i]
        node_type = node.__class__.__name__
        if node_type.endswith("Record"):
            node_type = node_type[:-6]
        try:
            if node_type[:2] == "V9":
                node_type += str(getattr(node, node_type[3:] + "Count"))
            else:
                node_type += str(getattr(node, node_type + "Count"))
        except AttributeError:
            pass

        try:
            node_label = node.Label
        except AttributeError:
            if node_type[:5] == "Pulse":
                node_label = self.bundle.header.Version
            else:
                node_label = ""
        if node_type[:2] == "V9":
            item = pg.QtGui.QTreeWidgetItem([node_type[3:], node_label])
        else:
            item = pg.QtGui.QTreeWidgetItem([node_type, node_label])
        if len(index) == 2:
            self.seriesNode.append(node)

        if dat_index == None:  ## update all
            root_item.addChild(item)
            item.node = node
            item.index = index

            if len(index) < 2:
                item.setExpanded(True)
            for i in range(len(node.children)):
                self.update_tree_recursive(item, index + [i])
        else:
            if len(index) == 2:  ## if this series in the selection list!
                self.seriesNode.append(node)
                if index == dat_index:  ## only add this series

                    root_item.addChild(item)
                    item.node = node
                    item.index = index
                    for i in range(len(node.children)):
                        self.update_tree_recursive(item, index + [i])
            else:
                root_item.addChild(item)
                root_item.setExpanded(True)
                item.node = node
                item.index = index
                if len(index) < 2:
                    item.setExpanded(True)
                for i in range(len(node.children)):
                    self.update_tree_recursive(item, index + [i], dat_index)

    def update_ABFtree(self, root_item):
        """read information from the Axon ABF file using neo.AxonIO interface
        and add items into the GUI tree to allow browsing.
        """
        r = self.abf  ## root
        self.filetype = ".abf"
        self.dat_file = self.frame.currentPulseTree.dat_file

        node = r
        if node == None:
            return
        self.frame.parameters["fs"] = r.get_signal_sampling_rate()
        metaInfo = r._axon_info  ## more information in lit
        self.abf.yUnits = r.header["signal_channels"]["units"][0]  ## 'mV' for signal
        self.abf.xUnits = "s"
        nBlock = r.header["nb_block"]  ## could have multiple blocks
        nSegments = r.header[
            "nb_segment"
        ]  ## this is a list; number of segements for each block
        ## root level
        node_type = "ABF"
        node_label = (
            metaInfo["fFileSignature"].decode()
            + "_v"
            + str(metaInfo["fFileVersionNumber"])
        )  # abfVersion
        abfNode = pg.QtGui.QTreeWidgetItem([node_type, node_label])
        root_item.addChild(abfNode)
        root_item.setExpanded(True)
        abfNode.node = r
        abfNode.index = 0

        self.abfBlocks = r.read()  # read the entire file > a list of Blocks
        ## block level
        for b in range(nBlock):
            block = pg.QtGui.QTreeWidgetItem(["Block", str(b + 1)])
            block.setExpanded(True)
            block.index = [0, b]
            abfNode.addChild(block)
            for seg in range(nSegments[b]):
                segment = pg.QtGui.QTreeWidgetItem(["Sweep" + str(seg + 1), ""])
                segment.setExpanded(True)
                segment.index = [0, b, seg]
                block.addChild(segment)
        rawP = r.read_raw_protocol()
        stimChanIdx = 0  ## channel index where stimuli is applied to
        if len(rawP[0]) > 0:
            self.abf_stimOn = True
            nSamples = len(rawP[0][0][0])
            stimData = np.zeros((nSamples, nSegments[0]))
            for j, stim in enumerate(rawP[0]):
                stimData[:, j] = stim[stimChanIdx]

            self.abf_stimData = stimData
            self.abf_stimTime = np.arange(nSamples) / r.get_signal_sampling_rate()
            self.abf_stimUnit = rawP[2][stimChanIdx]  ## 'pA'
        else:
            self.abf_stimOn = False  ## spontaneous events!
            print("spontaneous events!")
            ## evoke event analysis panel automatically here!
            # self.frame.eventDetectionAction_clicked()

    def get_plot_params(self):
        """
        Gets information for plotting.
        :return:
        """
        if self.indices is not None:

            ret = []

            for trace in self.indices:
                data = self.bundle.data[trace]

                pul = self.pul[trace[0]][trace[1]][trace[2]][trace[3]]
                y_label = pul.Label
                y_units = pul.YUnit
                x_interval = pul.XInterval

                ret.append((data, x_interval, y_label, y_units))

            return ret

        else:
            return [(None, None, None, None)]

    #    @QtCore.pyqtSlot(QtCore.QModelIndex)
    def on_selection_changed(self):
        """
        Event for browsing through wave traces
        :return:
        """
        selected = self.selectedItems()
        indices = []

        for j, item in enumerate(selected):  ### for muliple selection
            index = item.index
            # print(index)
            if len(index) == 2:  ## grab series level of data
                indices.append(index)

        self.indices = indices
        self.frame.currentPulseTree = self
        self.frame.update_trace_plot()


class OptionsView(pg.QtWidgets.QWidget):
    """
    Class for view containing plotting options.
    """

    def __init__(self, parent):
        super(OptionsView, self).__init__()
        self.frame = parent

        # init inputs
        self.extract_spikes_toggle = QCheckBox("Extract spikes")
        self.extract_spikes_toggle.setChecked(False)

        self.group_toggle = QCheckBox("Group spikes")
        self.group_toggle.setChecked(False)

        self.arg_type =QComboBox()
        self.arg_type.addItems(["max", "min"])
        self.arg_type_label = QLabel("Center on window")

        self.spike_edge = QComboBox()
        self.spike_edge.addItems(["rising", "falling"])
        self.spike_edge_label = QLabel("Threshold edge")

        self.spike_thresh = QSpinBox()
        self.spike_thresh.setKeyboardTracking(True)
        self.spike_thresh.setRange(1, 10)  # ms
        self.spike_thresh.setValue(2)
        self.spike_thresh.setSingleStep(1)
        self.spike_thresh_label = QLabel("Spike thresh")

        self.group_window = QDoubleSpinBox()
        self.group_window.setKeyboardTracking(False)
        self.group_window.setRange(0.01, 1000)  # ms
        self.group_window.setValue(2.5)
        self.group_window.setSingleStep(2.5)
        self.group_window_label =QLabel("Group window (ms)")

        # layout and add
        layout = QVBoxLayout()
        grid_layout = QGridLayout()

        grid_layout.addWidget(self.extract_spikes_toggle, 0, 0, 1, 2)
        grid_layout.addWidget(self.group_toggle, 1, 0, 1, 2)
        grid_layout.addWidget(self.arg_type_label, 2, 0, 1, 1)
        grid_layout.addWidget(self.arg_type, 2, 1, 1, 1)
        grid_layout.addWidget(self.spike_edge_label, 3, 0, 1, 1)
        grid_layout.addWidget(self.spike_edge, 3, 1, 1, 1)
        grid_layout.addWidget(self.spike_thresh_label, 4, 0, 1, 1)
        grid_layout.addWidget(self.spike_thresh, 4, 1, 1, 1)
        grid_layout.addWidget(self.group_window_label, 5, 0, 1, 1)
        grid_layout.addWidget(self.group_window, 5, 1, 1, 1)

        layout.addLayout(grid_layout)
        layout.addStretch(1)

        self.setLayout(layout)

        # event binders
        self.extract_spikes_toggle.toggled.connect(self.on_extract_spikes_toggle)
        self.group_toggle.toggled.connect(self.on_group_button)
        self.arg_type.currentIndexChanged.connect(self.arg_type_changed)
        self.spike_edge.currentIndexChanged.connect(self.spike_edge_changed)
        self.spike_thresh.valueChanged.connect(self.on_spike_thresh_changed)
        self.group_window.valueChanged.connect(self.on_group_window_changed)

    def on_extract_spikes_toggle(self, state):
        """
        Draw extracted spikes.
        :param state:
        :return:
        """
        if state:
            if not self.frame.spike_view.isVisible():
                self.frame.spike_view.resize(
                    self.frame.trace_view.width(), self.frame.trace_view.height() // 2
                )
                self.frame.spike_view.show()
            self.frame.update_spike_plot()

        elif not state:
            if self.frame.spike_view.isVisible():
                self.frame.spike_view.hide()

    def on_group_button(self, state):
        """
        Toggle grouping of spikes.
        :param state:
        :return:
        """
        self.frame.toggle_grouping(state)

    def on_group_window_changed(self, value):
        """
        Toggle grouping of spikes.
        :param value:
        :return:
        """
        self.frame.update_group_window(value)

    def on_spike_thresh_changed(self, value):
        """
        Toggle grouping of spikes.
        :param value:
        :return:
        """
        self.frame.update_spike_threshold(value)

    def arg_type_changed(self, index):
        """
        Toggle grouping of spikes.
        :param value:
        :return:
        """
        options = ["max", "min"]

        self.frame.update_arg_type(options[index])

    def spike_edge_changed(self, index):
        """
        Toggle grouping of spikes.
        :param value:
        :return:
        """
        options = ["rising", "falling"]

        self.frame.update_spike_edge(options[index])


class TableView(pg.TableWidget):
    """
    Class for a table
    """

    def __init__(self, parent, title="", editable=False, sortable=True):
        super(TableView, self).__init__(editable=editable, sortable=sortable)
        self.frame = parent
        self.name = title
        self.setWindowTitle(title)
        self.contextMenu.addAction(
            translate("TableWidget", "Load data")
        ).triggered.connect(self.LoadData)
        self.contextMenu.addAction(
            translate("TableWidget", "Clear all")
        ).triggered.connect(self.clear)

    def LoadData(self):
        file_name = pg.QtGui.QFileDialog.getOpenFileName(
            self, "Load saved file", "", "TSV Files (*.tsv)"
        )
        # print(file_name)
        if file_name[0] == "":
            return
        else:
            x = pd.read_table(file_name[0])
            x = x.set_index("ChanID")
            data = x.to_records()
            self.setData(
                np.array(data, dtype=[("ChanID", object), ("X", object), ("Y", object)])
            )
            self.frame.roidata = list(data)
            self.frame.redrawAllROIs()

    def keyPressEvent(self, event):
        key = event.key()
        if key == QtCore.Qt.Key_Return or key == QtCore.Qt.Key_Enter:
            # Process current item here
            self.frame.eventTable_toggleState(self.currentRow(), self.currentColumn())
        else:
            super(TableView, self).keyPressEvent(event)


class TabView(pg.QtWidgets.QTabWidget):
    """
    Class for tab view.
    """

    def __init__(self, parent, enableContextMenu=False):
        super().__init__()
        self.frame = parent
        self.enableContextMenu = enableContextMenu
        if enableContextMenu:
            self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
            self.customContextMenuRequested.connect(self.openContextMenu)
            self.actionItems = []

    def openContextMenu(self):
        idx = self.currentIndex()
        try:
            myMenu = QtWidgets.QMenu("Menu", self)
            myMenu.addAction(self.actionItems[idx])
            myMenu.exec_(QtGui.QCursor.pos())
        except IndexError:
            print("action item not set for this tab")

    def insertTabWithAction(self, tabObj, actions, label=""):
        ## inherit parent addTab function
        if isinstance(actions, list):
            for a in actions:
                assert isinstance(
                    a, pg.QtGui.QAction
                ), "second argument should be list of instance of QAction"
        else:
            assert isinstance(
                actions, pg.QtGui.QAction
            ), "second argument should be instance of QAction"
        if self.enableContextMenu:
            super().addTab(tabObj, label)
            ## store action item so we can add it on the fly
            self.actionItems.append(actions)
        else:
            raise AttributeError(
                "Set enableContextMenu as True before using this method"
            )


class MatplotView(MatplotlibWidget.MatplotlibWidget):
    """
    Matplotlib figure widget
    """

    def __init__(self, parent=None, size=(5.0, 4.0), title="", dpi=100, plotGrid=None):
        """plotGrid is a tuple with (nrow, ncol)"""
        super(MatplotView, self).__init__(size=size, dpi=dpi)
        self.setParent(parent)
        self.figure = self.getFigure()
        self.clf()

    def subplots(self, nrow, ncol):
        ## set up subplots
        self.axes = [
            self.figure.add_subplot(nrow, ncol, j)
            for j in np.arange(1, nrow * ncol + 1)
        ]
        return self.axes

    def setPlotGrid(self, nRow, nCol):
        self.grid = self.figure.add_gridspec(nRow, nCol)

    def plotline(self, i, j, x, y, color="r", title="", linewidth=0.5, ax=None):
        if i >= self.grid.nrows or j >= self.grid.ncols:
            raise ValueError("row index or col index out of bound")
        if ax is None:
            ax = self.figure.add_subplot(self.grid[i, j])
        ax.plot(x, y, color=color, linewidth=linewidth)
        if len(title) > 1:
            ax.set_title(title, fontsize=12, color="k")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
        return ax

    def adjustSubplot(
        self, left=0.08, bottom=0.1, top=0.95, right=0.95, wspace=0.35, hspace=0.35
    ):
        self.figure.subplots_adjust(
            left=left, bottom=bottom, top=top, right=right, wspace=wspace, hspace=hspace
        )

    def clf(self):
        self.figure.clear()


class PlotView(GraphicsLayoutWidget):
    """
    Class for plot view.
    """

    def __init__(self, parent, title="", background="k"):
        super(PlotView, self).__init__()
        self.frame = parent
        self.setWindowTitle(title)
        self.setBackground(background)
        #        self.layout = pg.GraphicsLayoutWidget(title = title)
        self.ci.setBorder()
        self.ci.setBorder(pg.mkPen(None))
        self.ci.setContentsMargins(10, 10, 10, 20)
        self.setCentralItem(self.ci)
        self.show()

    def refresh(self):
        self.ci.setContentsMargins(10, 10, 10, 20)


class ParameterView:
    """
    Class for enhance pyqtgraph.parametertreeParameter class
    """

    def __init__(
        self, parent, name="params", type="group", children=None, readonly=False
    ):
        self.parent = parent
        self.params = children  ## actual parameters
        self.state = []  ## save and restore
        self.p = Parameter.create(
            name=name, type=type, children=children, readonly=readonly
        )

    def setChildValue(self, childPath, childValue):
        """
        Set value of a child
        ----------
        childPath : list of string
            the full path for this child
        childValue :
        Returns
        -------
        None.
        """
        vals = self.p.getValues()
        v = None
        depths = 0
        expression = "vals"
        while v == None:
            child = vals[childPath[depths]]
            v = child[0]
            if v == None:
                index = "[1]"
            else:
                index = "[0]"
            expression = "".join("[" + childPath[depths] + "]" + index)
            depths += 1
        eval(expression)
        self.p.setValue(vals)

    def change(self, param, changes):
        for param, change, data in changes:
            path = self.p.childPath(param)
            if path is not None:
                childName = ".".join(path)
            else:
                childName = param.name()
            print("  Parameter: %s" % childName)
            print("  Value:      %s" % str(data))
            print("  ----------")
        self.parent.rePlotSpikePanels()

    def valueChanging(param, value):
        print("Value changing (not finalized): %s %s" % (param, value))

    def save(self):
        self.state = self.p.saveState()

    def restore(self):
        add = self.p["Save/Restore parameters", "Restore State", "Add missing items"]
        rem = self.p["Save/Restore parameters", "Restore State", "Remove extra items"]
        self.p.restoreState(self.state, addChildren=add, removeChildren=rem)


class miniEvent:
    def __init__(self):
        self.init()

    def init(self):
        self.node = []  ## store series name
        self.isConcat = False
        self.eventType = []
        if not hasattr(
            self, "datFile"
        ):  ## we do not want to lose this information when change sweep or tracw within the same file!
            self.datFile = []
        self.seriesName = []
        self.traceName = []
        self.nSweep = []
        self.currentSweep = []
        self.time = []  ## current whole traces
        self.data = []  ## current whole traces. filtered.
        self.data_raw = []  ## raw data
        self.dataMin = -1
        self.dataMax = 1
        self.m = 0
        self.std = 0
        self.D = []  ## deconvoled
        self.sweepType = "single"  ## or 'all' for all sweeps; 'single' for single sweep
        self.eventMarker = None  ## scatterplotItem for marking detected events
        self.traceXUnit = []
        self.traceYUnit = []
        self.eventTime = {}  ## event only
        self.eventData = {}  ## event only
        self.eventTimeStamps = {}
        self.peakStart = []
        self.peakIndex = []
        self.globalBaseline = []
        self.downstroke_samples = 0
        self.upstroke_samples = 0
        self.eventTableDF = []
        self.template = []
        self.templateTwinLine = []  # list of template twin line
        self.eventLineReloaded = False  ## a flag to avoid repeatly loading event
        self.decovPlot = ""

        if hasattr(self, "tracePlot"):
            self.tracePlot.clear()
        else:
            self.tracePlot = []

    def setSeriesSampleRate(self, fs):
        self.seriesFS = fs

    def updateCurrentSweep(self):
        """
        save events data in dictionary using sweep number as key
        """
        eventTime = []
        eventData = []
        eventList = []
        # time = self.time.get(self.currentSweep)
        # data = self.data.get(self.currentSweep)
        for l in self.tracePlot.vb.allChildItems():
            if hasattr(l, "name"):
                if l.name() == "template1":
                    startPos = l.pos().x()
                    endPos = l.twin.pos().x()
                    index1 = int(startPos * self.seriesFS)
                    index2 = int(endPos * self.seriesFS)
                    time_ = self.time[index1:index2]
                    data_ = self.data[index1:index2]
                    eventTime.append(time_)
                    eventData.append(data_)
                    eventList.append([startPos, endPos])
        self.eventTimeStamps.update({self.currentSweep: eventList})
        self.eventTime.update({self.currentSweep: eventTime})
        self.eventData.update({self.currentSweep: eventData})
        # print('Event cache updated')

    def get_events(self):
        eventTime = []
        eventData = []
        for sweep in self.eventTime.keys():
            eventTime.extend(self.eventTime.get(sweep))
            eventData.extend(self.eventData.get(sweep))
        return eventTime, eventData

    def clear(self):
        self.init()


class SplitView(pg.QtWidgets.QSplitter):
    """A boiler plate for constructing multipanel tab.
    Top row may contain multiple figures with matplotlib backend
    Bottom row split horinzontally into two.
    Left half can be used to add smaller tabs, such as a parameter view.
    Right Half can be used to add tabs needing more space, such add tables
    """

    def __init__(self, title, FigureTitle="Figures"):
        super().__init__(QtCore.Qt.Vertical)
        self.title = title
        self.matplotViews = {}  ## dictionry to add matplot views
        self.tables = {}  ## dictionry to add table
        self.parTreeViews = {}  ## dictionry to add parameter trees
        self.parsTree_values_init = {}
        self.func = ""
        self.topSplitter = QSplitter(QtCore.Qt.Horizontal)  ## bottom row
        self.bottomSplitter = QSplitter(QtCore.Qt.Horizontal)  ## bottom row
        self.addWidget(self.topSplitter)
        self.addWidget(self.bottomSplitter)

        self.setStretchFactor(0, 3)
        self.setStretchFactor(1, 1)

        self.top_tabs = TabView(self)  ## top row
        self.bottom_tabs = TabView(self)  ## bottom row

        # top row
        self.topSplitter.addWidget(self.top_tabs)
        self.plotSplitter = QSplitter(QtCore.Qt.Horizontal)  ## bottom row
        self.top_tabs.addTab(self.plotSplitter, FigureTitle)

        # bottom row
        ## add a splitter. Left for parameters, Right for tables
        self.bottomSplitter = QSplitter(QtCore.Qt.Horizontal)  ## bottom row
        self.addWidget(self.bottomSplitter)

        ## bottom left tab areas
        self.bottomLeft_tabview = TabView(self)
        self.bottomRight_tabview = TabView(self)
        self.bottomSplitter.addWidget(self.bottomLeft_tabview)
        ## bottom right tab areas
        self.bottomSplitter.addWidget(self.bottomRight_tabview)
        self.bottomSplitter.setStretchFactor(0, 1)
        self.bottomSplitter.setStretchFactor(1, 4)

    def addImageviewatTop(self):
        self.slice_view = pg.ImageView(
            view=pg.PlotItem()
        )  ## add an image window for slice image
        self.plotSplitter.addWidget(self.slice_view)

    def addMatPlotviewAtTop(self, plotNames, size=(5.0, 4.0), title="", dpi=100):
        """
        Use this function to add new matplotlib widget which will be saved
        into 'matplotViews' dictionary using 'plotName' as key
        """
        if isinstance(plotNames, list):
            for name in plotNames:
                matplotView = MatplotView()
                self.plotSplitter.addWidget(matplotView)
                self.matplotViews.update({name: matplotView})
        else:
            matplotView = MatplotView(size=size, title=title, dpi=dpi)
            self.plotSplitter.addWidget(matplotView)
            self.matplotViews.update({plotNames: matplotView})

    def addTablesAtBottomRight(self, tableNames, editable=False, sortable=True):
        """
         Use this function to add table widget which will be saved
        into 'tables' dictionary using 'tableName' as key
        """
        if isinstance(tableNames, list):
            for name in tableNames:
                table = TableView(self, "", editable, sortable)
                self.tables.update({name: table})
                self.bottomRight_tabview.addTab(table, name)
        else:
            table = TableView(self)
            self.tables.update({tableNames: table})
            self.bottomRight_tabview.addTab(table, tableNames)

    def addParameterToLeftTab(self, title, pars, func=None, tooltips=None):
        """use this to add a parameter tree.
        title: tab name for this table
        pars: parameters to be loaded into this tree widget
        """
        parsTree = ParameterTree()
        parsTree.setHeaderLabels(["Parameter                 ", "Value"])
        parsTree_view = ParameterView(self, name="params", type="group", children=pars)
        if not func == None:
            parsTree_view.p.sigTreeStateChanged.connect(func)
        # self.fpParTree_data_view.p.sigTreeStateChanged.connect(self.fpParTree_stateChange)
        # self.thresholdPar.sigValueChanged.connect(self.thresholdParChange)
        parsTree.setParameters(parsTree_view.p, showTop=False)
        if tooltips is not None:
            parsTree = self.addTooltips(parsTree, tooltips)
        self.bottomLeft_tabview.addTab(parsTree, title)  ##
        self.parTreeViews[title] = parsTree_view
        self.parsTree_values_init[title] = False
    
    def addTooltips(self, parTree, tooltips):
        ''' Add tooltips for parameter tree'''
        for idx, tt in enumerate(tooltips):
            parTree.children()[0].children()[idx].setToolTip(tt['value'])
        return parTree
    
    def setParTreePars(self, pars, index=0, title="Data selection", func=None):
        self.bottomLeft_tabview.setCurrentIndex(index)
        parsTree = ParameterTree()
        parsTree.setHeaderLabels(["Parameter", "Value"])
        parsTree_view = ParameterView(self, name="params", type="group", children=pars)
        if not func == None:
            parsTree_view.p.sigTreeStateChanged.connect(func)
        # self.fpParTree_data_view.p.sigTreeStateChanged.connect(self.fpParTree_stateChange)
        # self.thresholdPar.sigValueChanged.connect(self.thresholdParChange)
        parsTree.setParameters(parsTree_view.p, showTop=False)
        self.bottomLeft_tabview.removeTab(index)
        self.bottomLeft_tabview.insertTab(
            index, parsTree, QtGui.QIcon("Data/connection.png"), title
        )
        self.parTreeViews.update({title: parsTree_view})
        self.parsTree_values_init.update({title: True})
        # self.parsTree_view.p.setValue(pars)
        # self.parsTree.setParameters(self.parsTree_view.p, showTop=False)

    def getParTreePars(self, title):
        return self.parTreeViews[title].p.getValues()

    def setParTreeParsStateChangeFunc(self, title, func):
        self.parTreeViews[title].p.sigTreeStateChanged.connect(func)
        ## TODO: do we need this?
        # self.parsTree.setParameters(self.parsTree_view.p, showTop=False)

    def reset(self):
        """reset this widget"""

        for table in self.tables:
            self.tables[table].clear()

        for plot in self.matplotViews:
            self.matplotViews[plot].figure.clear()
            self.matplotViews[plot].clf()
        # self.parTreeViews = {}
        # self.parsTree_values_init = {}


class infiniteLinePair(pg.InfiniteLine):
    """
    modified ininite line class to for paried lines and support Alt-click removable ability
    Plus auto label! How cool is that!
    """

    def __init__(
        self, parent=None, offset=0, pos=None, index=None, Frame=None, **kwargs
    ):
        super(infiniteLinePair, self).__init__(pos=pos, **kwargs)
        self.parent = parent
        self.Frame = Frame  ## main frame object
        self.index = index  ## index of the data
        self.offset = offset
        self.twin = pg.InfiniteLine(
            pos=pos + offset,
            span=(0, 0.8),
            angle=90,
            movable=True,
            pen="g",
            hoverPen="w",
            name="templateLine",
        )

        self.lineLable = pg.InfLineLabel(
            self, text="%0.3f" % (pos), movable=True
        )  ## add a label
        self.active = True

    def mouseDragEventTwin(self):
        self.Frame.updateEvents()

    def mouseDragEvent(self, ev):
        if self.movable and ev.button() == QtCore.Qt.LeftButton:
            if ev.isStart():
                self.moving = True
                self.cursorOffset = self.pos() - self.mapToParent(ev.buttonDownPos())
                self.startPosition = self.pos()
            ev.accept()
            if not self.moving:
                return
            self.setPos(self.cursorOffset + self.mapToParent(ev.pos()))
            self.lineLable.setFormat("%0.3f" % (self.pos().x()))  ## update label text
            self.sigDragged.emit(self)
            if ev.isFinish():
                self.moving = False
                self.sigPositionChangeFinished.emit(self)
                self.Frame.updateEvents()

    def mouseClickEvent(self, ev):
        if ev == "":
            return
        if self.moving and ev.button() == QtCore.Qt.RightButton:
            ev.accept()
            self.setPos(self.startPosition)
            self.moving = False
            self.sigDragged.emit(self)
            self.sigPositionChangeFinished.emit(self)
        if ev.modifiers() == QtCore.Qt.AltModifier:  ## remove this item
            self.parent.removeItem(self)
            self.parent.removeItem(self.twin)
            self.active = False
            self.Frame.updateEvents()

    def updateWindowLength(self, newOffset):
        self.twin.setValue(self.pos() + newOffset)


class EventMark(pg.InfiniteLine):
    """
    modified ininite line class to for event line
    """

    def __init__(
        self,
        pos=None,
        timeIdx=None,
        penc="g",
        userGen=False,
        userMarker=True,
        spanc=(0.88, 0.92),
        hoverPenc="w",
        markerSize=7,
        **kwargs,
    ):
        super(EventMark, self).__init__(
            pos=pos,
            angle=90,
            pen=penc,
            hoverPen=hoverPenc,
            movable=True,
            name="eventLine",
            span=spanc,
            **kwargs,
        )
        self.active = True
        self.timeIdx = timeIdx
        self.userGen = userGen  ## generated by user
        if userMarker:
            self.addMarker("o", position=0.9, size=markerSize)

    def mouseClickEvent(self, ev):
        if ev == "":
            return
        # if ev.modifiers() == QtCore.Qt.AltModifier:  ## switch states
        if self.active:
            self.setPen((155, 155, 155, 255))
            self.setHoverPen("g")
            print("Event invalidated")
        else:
            # self.clearMarkers()
            # self.addMarker('^')
            self.setPen("g")
            self.setHoverPen((155, 155, 155, 255))
            print("Event recovered")
        self.active = not (self.active)

    def mouseDragEvent(self, ev):
        pass