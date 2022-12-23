"""
Patchviewer by M.H.
"""
from PyQt5 import QtGui, QtWidgets, QtCore
from PyQt5.QtWidgets import QMessageBox, QTableWidgetItem
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from pyqtgraph.parametertree import Parameter, ParameterTree
from patchview.utilitis import AllMyParsHere as AllMyPars
import sys
import os
import numpy as np
import pandas as pd
from patchview.HekaIO import HEKA_Reader_MAIN as HEKA
from patchview.HekaIO.HekaHelpers import HekaBundleInfo
from scipy import signal
import scipy.stats as stats
from scipy.optimize import curve_fit
import seaborn as sns
sns.despine()
from patchview.ephys import ephys_features as ephys_ft
from patchview.ephys import ephys_extractor as efex
from patchview.ephys import extraEhpys_PV
from patchview.ephys.extraEhpys_PV import extract_sweep_feature
from scipy.stats import pearsonr
import itertools
import time as sysTime
from copy import deepcopy
import glob
from patchview.utilitis.debugHelpers import debugInfo
from patchview.utilitis.AnalysisMethods import (
    loadYAML,
    filterDatSeries,
    calculateConnectionTraces,
    cleanASCfile,
    smooth, smooth2D, padding_
)
from neurom.geom.transform import (translate, rotate)
from patchview.utilitis import fitFuncs
from patchview.utilitis.morphorFeatureExtrator import (
    getSomaStats, extractMorphFeatures, sholl_analysis,
     sholl_single_axis,sholl_2D_density, sholl_polar)
from patchview.utilitis.patchviewObjects import *
import networkx as nx

## reading neuroluscida file
from morphio import SomaType
import neurom as morphor_nm
from neurom import viewer as morphor_viewer
from neurom.io.multiSomas import MultiSoma
# from neurom.features.utilities import (getSomaStats, extractMorhporFeatures, sholl_analysis)
from neurom.core.morphology import Morphology
from neurom.core.types import NeuriteType
import neo
from neo.io import AxonIO
import warnings
import matplotlib.pyplot as MATPLT
from matplotlib.ticker import FormatStrFormatter
from matplotlib.collections import LineCollection
from sklearn import linear_model
from matplotlib import image as ReadImage
from patchview.utilitis.PVdat2NWB import dat2NWB
import pynwb
from hdmf.spec import NamespaceCatalog
from hdmf.utils import docval, getargs, popargs, call_docval_func, get_docval
from hdmf.backends.io import HDMFIO
from hdmf.backends.hdf5 import HDF5IO as _HDF5IO
from hdmf.validate import ValidatorMap
from hdmf.build import BuildManager, TypeMap
import hdmf.common
from hdmf.spec import NamespaceCatalog  # noqa: E402
from hdmf.utils import docval, getargs, call_docval_func, get_docval, fmt_docval_args  # noqa: E402
from hdmf.build import BuildManager, TypeMap  # noqa: E402
warnings.filterwarnings("ignore")
from appdirs import *
patchview_dir, this_filename = os.path.split(__file__)
appname = "Patchview"
__version__ = "0.2.6.0"
NeutriteColors = {NeuriteType.apical_dendrite:'m',
NeuriteType.basal_dendrite: 'b',
NeuriteType.axon:'r'}
class MainWindow(QtWidgets.QMainWindow):
    """
    Main frame.
    """

    def __init__(self, app, parent=None):
        super(MainWindow, self).__init__(parent)
        self.app = app
        self.setUserProfile()
        self.EphyFeaturesObj = []
        self.spikeTableSavePath = ""
        self.create_mainWindow()
        self.setWindowTitle("PatchView")
        self.setWindowIcon(
            pg.QtGui.QIcon(os.path.join(patchview_dir, "Data", "icons", "PatchViewer.ico"))
             
        )
        #        self.resize(1200, 800)
        self.showMaximized()
        


    def create_mainWindow(self):
        """Configure Qt GUI:
        Main window + splitters to let user resize panes
        """
        
        self.make_layout()
        ## menu bar
        self.add_menubar()
        self.add_toolbar()

        if not self.showSpikePanelsFlag:
            self.hideSpikePanels()
        self.statusBar = QtWidgets.QStatusBar()
        self.setStatusBar(self.statusBar)
        self.MouseMode = pg.ViewBox.RectMode
        self.setCentralWidget(self.mainFrame)
        self.plotItemList = []  ## collect all plot handles
        self.currentAnMorphView= []
        self.mainFrame.show()

    def add_selectedDatFile(self, filePath):
        """adding file to the file list if not exist already"""
        i1 = QtGui.QTreeWidgetItem([filePath])
        self.fileList.addTopLevelItem(i1)

    def make_layout(self):
        self.mainFrame = pg.LayoutWidget() 
        self.MainLayout =  self.mainFrame.layout
        self.frame_splitter = QSplitter(QtCore.Qt.Horizontal) ## main splitter
        self.MainLayout.addWidget(self.frame_splitter, 1, 0, 1, 1)
        self.file_view = FileView(self, self.root)
        self.options_view = OptionsView(self)
        self.cellFeatures_view = TableView(self)
        self.sweepFeatures_view = TableView(self)
        self.spikeTable_view = TableView(self)  # a talbe of detected spikes
        # self.sliceView.tables['ROI list'] = TableView(self)  # a talbe of detected spikes
        self.roidata = []  ## roi coordinates
        self.neuronMorph =  None 
        # list .dat files
        self.files_view = TabView(self)
        self.files_view.addTab(self.file_view, "Dat file")
        # self.tab_view.addTab(self.options_view, 'Options')
        ## Create two ParameterTree widgets, both accessing the same data
        self.globalSettings = ParameterView(
            self, name="params", type="group", children=AllMyPars.params
        )
        self.parTree = ParameterTree()
        self.parTree.setParameters(self.globalSettings.p, showTop=False)
        self.parTree.setWindowTitle("Global parameter")
        self.files_view.addTab(self.parTree, "global parameters")
        ## list of pulse tree
        self.trees_view = TabView(self)
        self.PulseTrees = []
        for k in range(1):  ## Eight channels!
            self.pul_view = PulView(self)
            self.PulseTrees.append(
                self.pul_view
            )  ## PulseTrees object has the actual pulview
            self.trees_view.addTab(self.pul_view, "Current tree")
        self.trees_view.setCurrentIndex(0)
        self.currentPulseTree = self.PulseTrees[
            0
        ]  ## this is to update current pulse tree to be manipulated
        self.trees_view.currentChanged.connect(self.updateCurrentPulseTree)
        # splitter for file browser (top) and  pul view (bottom)
        self.tree_splitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.tree_splitter.setStretchFactor(0, 1)
        self.tree_splitter.setStretchFactor(1, 9)
        ### top panel: file browser (left), list of selected files (right)
        self.tree_splitter.addWidget(self.files_view)

        ### middle panel: current dat file tree (left), selected series/sweeps (right)
        self.selectedTrees_splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.tree_splitter.addWidget(self.selectedTrees_splitter)

        self.selectedFilesTabView = TabView(self, enableContextMenu=True)

        self.selTrees = SelectionView(self)  ## this is to select connection
        # self.selTrees.setColumnCount(2)
        self.selectedTrees_splitter.addWidget(self.selectedFilesTabView)
        self.selectedTrees_splitter.addWidget(self.trees_view)

        ## Add Firing pattern protocol to the selection region
        self.selTreesFP = SelectionView(self)  ## this is to select FPs
        ctActionFP = pg.QtGui.QAction("&Firing patterns (Batch)")
        ctActionFP.triggered.connect(self.checkFP_action_clicked)
        self.selectedFilesTabView.insertTabWithAction(
            self.selTreesFP, ctActionFP, "Firing pattern"
        )

        ctAction = pg.QtGui.QAction("&Synaptic connections (Batch)")
        ctAction.triggered.connect(self.checkConnectionAction_clicked)
        self.selectedFilesTabView.insertTabWithAction(
            self.selTrees, ctAction, "Connection"
        )

        self.selTreesSP = SelectionView(
            self
        )  ## this is to select .dat for Spontaniouss
        ctActionSP = pg.QtGui.QAction("&Event detection")
        ctActionSP.triggered.connect(self.checkEvent_action_clicked)
        self.selectedFilesTabView.insertTabWithAction(
            self.selTreesSP, ctActionSP, "Spontaneous"
        )

        self.selTreesDispatcher = {
            "Firing pattern": self.selTreesFP,
            "Connection": self.selTrees,
            "Spontaneous": self.selTreesSP,
        }

        ### Bottom panel: recording parameters
        self.parameters_view = TabView(self)
        self.tree_splitter.addWidget(self.parameters_view)
        self.parameter_Tab = TableView(self)
        self.parameters_view.addTab(self.parameter_Tab, "Recording parameters")

        # splitter for plots
        self.plot_splitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)

        self.topPlot_splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)  ## first row
        self.plot_splitter.addWidget(self.topPlot_splitter)

        self.trace_view = PlotView(self.topPlot_splitter)  ## first row, first column
        self.topPlot_splitter.addWidget(self.trace_view)
        self.trace_view2 = PlotView(
            self.topPlot_splitter
        )  ## first row, second column. Show spikes of selected sweep
        self.topPlot_splitter.addWidget(self.trace_view2)
        self.topPlot_splitter.setStretchFactor(0, 2)
        self.topPlot_splitter.setStretchFactor(1, 1)

        self.matplot_splitter = QtWidgets.QSplitter(
            QtCore.Qt.Horizontal
        )  ##  second row. i want to split the bottom into tow colomns
        self.plot_splitter.addWidget(self.matplot_splitter)
        self.matplot_splitter.setStretchFactor(0, 2)
        self.matplot_splitter.setStretchFactor(1, 2)

        self.spikes_view = MatplotView()  ## group of all spikes
        self.matplot_splitter.addWidget(self.spikes_view)
        self.spikeSummary_view = MatplotView()
        self.matplot_splitter.addWidget(self.spikeSummary_view)
        self.showSpikePanelsFlag = 0

        ### main plotting and data view panels
        self.visulization_view = TabView(self)
        self.visulization_view.addTab(self.plot_splitter, "Browser")

        ## Firing pattern tabls
        splitViewTab_FP = SplitView("Firing pattern")

        splitViewTab_FP.addMatPlotviewAtTop(
            "Firing pattern"
        )  ## show FP representive curves
        splitViewTab_FP.addParameterToLeftTab("Data selection", AllMyPars.fp_analysis)
        spdet_toolTips = AllMyPars.fp_analysis_spikeDetection_help[0]["children"]
        splitViewTab_FP.addParameterToLeftTab(
            "Spike detection",
            AllMyPars.fp_analysis_spikeDetection,
            self.spike_detection_event,
            spdet_toolTips
        )
        # splitViewTab_FP.addParameterToLeftTab(
        #     "Spike detection help", AllMyPars.fp_analysis_spikeDetection_help
        # )
        splitViewTab_FP.bottomLeft_tabview.setCurrentIndex(
            1
        )  # default showing tab for "spike detection"

        splitViewTab_FP.addTablesAtBottomRight(
            ["Sweep features", "Spike features", "Cell features"]
        )
        splitViewTab_FP.bottomRight_tabview.setCurrentIndex(
            2
        )  # default showing tab for "Cell features"
        self.splitViewTab_FP = splitViewTab_FP
        self.visulization_view.addTab(self.splitViewTab_FP, self.splitViewTab_FP.title)

        splitViewTab_connection = SplitView("Connections")
        splitViewTab_connection.addMatPlotviewAtTop(["Average traces", "Graph"])
        splitViewTab_connection.addParameterToLeftTab(
            "TBD", AllMyPars.connection_analysis
        )
        splitViewTab_connection.addTablesAtBottomRight("Detected connections")
        self.splitViewTab_connection = splitViewTab_connection
        self.visulization_view.addTab(
            self.splitViewTab_connection, self.splitViewTab_connection.title
        )

        self.sliceView = SplitView("Slice view", "Slice")  # self.sliceView.slice_view
        self.sliceView.addImageviewatTop()
        # self.selectedFilesTabView.addTab(self.sliceView.tables['ROI list'], 'ROI list')
        self.sliceView.addTablesAtBottomRight(["ROI list"])

        self.visulization_view.addTab(self.sliceView, self.sliceView.title)
        ## Event detection GUI
        self.event_Vsplitter = QSplitter(QtCore.Qt.Vertical)
        self.visulization_view.addTab(self.event_Vsplitter, "Event detection")

        ## first row, first column
        self.event_TopSplitter = QSplitter(QtCore.Qt.Horizontal)  ## top row
        self.event_Vsplitter.addWidget(self.event_TopSplitter)
        self.event_traceView = PlotView(self.event_TopSplitter)  ## plot original trace

        self.event_TopSplitter.addWidget(self.event_traceView)
        self.event_traceView2 = PlotView(self.event_TopSplitter)  ## plot extract event
        self.event_traceView2.addPlot(row=0, col=0)
        self.event_TopSplitter.addWidget(self.event_traceView2)
        self.event_TopSplitter.setStretchFactor(
            0, 6
        )  ## we want long and narrwow window
        self.event_TopSplitter.setStretchFactor(1, 4)
        ## second row, horinzontal s;itter
        self.event_BottomSplitter = QSplitter(QtCore.Qt.Horizontal)  ## bottom row

        ## add tabs for parameters
        self.eventPar_tabview = TabView(self)
        self.event_BottomSplitter.addWidget(self.eventPar_tabview)

        ## preprocesing tab
        self.eventParTree_data = ParameterTree()
        self.eventParTree_data_view = ParameterView(
            self,
            name="params",
            type="group",
            children=AllMyPars.event_detection_preprocessing,
        )
        self.eventParTree_data_view.p.sigTreeStateChanged.connect(
            self.eventParTree_stateChange
        )
        self.eventParTree_data.setParameters(
            self.eventParTree_data_view.p, showTop=False
        )
        self.eventPar_tabview.addTab(self.eventParTree_data, "Preprocessing")  ##

        ## template tab
        par_template = ParameterTree()  ## the widget
        self.event_template_view = ParameterView(
            self,
            name="params",
            type="group",
            children=AllMyPars.event_detection_template,
        )  # the parameter
        self.event_template_view.p.sigTreeStateChanged.connect(
            self.event_template_stateChange
        )  ## call back
        par_template.setParameters(
            self.event_template_view.p, showTop=False
        )  ## wiget paramter
        self.eventPar_tabview.addTab(par_template, "Template")  ##

        ## peak detection tab
        par_pd = ParameterTree()  ## the widget
        self.event_pd_view = ParameterView(
            self,
            name="params",
            type="group",
            children=AllMyPars.event_deconvPeak_parameters,
        )  # the parameter
        self.event_pd_view.p.sigTreeStateChanged.connect(
            self.event_pd_stateChange
        )  ## call back
        par_pd.setParameters(self.event_pd_view.p, showTop=False)  ## wiget paramter
        self.eventPar_tabview.addTab(par_pd, "Peak detection")  ##

        ## peak detection tab
        par_eventpp = ParameterTree()  ## the widget
        self.event_eventpp_view = ParameterView(
            self, name="params", type="group", children=AllMyPars.event_post_processing
        )  # the parameter
        self.event_eventpp_view.p.sigTreeStateChanged.connect(
            self.event_pd_stateChange
        )  ## call back
        par_eventpp.setParameters(
            self.event_eventpp_view.p, showTop=False
        )  ## wiget paramter
        self.eventPar_tabview.addTab(par_eventpp, "Event postprocessing")  ##

        ## help tab
        self.eventParTree_help = ParameterTree()
        self.eventParTree_help_view = ParameterView(
            self,
            name="params",
            type="group",
            children=AllMyPars.event_detection_helps,
            readonly=True,
        )
        self.eventParTree_help.setParameters(
            self.eventParTree_help_view.p, showTop=False
        )
        self.eventPar_tabview.addTab(self.eventParTree_help, "Helps")  ##

        ## add event data table
        self.eventData_tabview = TabView(self)
        self.event_BottomSplitter.addWidget(self.eventData_tabview)

        ## list of current template
        self.templateList = TableView(self)
        self.eventData_tabview.addTab(self.templateList, "Template List")

        self.templateFit = TableView(self)
        self.eventData_tabview.addTab(self.templateFit, "Template Fit")

        self.eventTable = TableView(self, "miniEvents")
        # self.eventTable.cellClicked.connect(self.eventTable_enterAction)
        self.eventTable.currentCellChanged.connect(self.eventTable_enterAction)
        self.eventTable.cellDoubleClicked.connect(self.eventTable_doubleClickAction)
        self.eventData_tabview.addTab(self.eventTable, "Event List")

        ## figure tab for event detection
        event_plotVsplitter = QSplitter(
            QtCore.Qt.Vertical
        )  ## split this if we want more plots
        self.event_matplotView1 = MatplotView()
        event_plotVsplitter.addWidget(self.event_matplotView1)
        self.eventData_tabview.addTab(event_plotVsplitter, "Figures")

        self.event_BottomSplitter.setStretchFactor(0, 1)
        self.event_BottomSplitter.setStretchFactor(1, 1)
        ## Done with bottom half. Add it to the frame
        self.event_Vsplitter.addWidget(self.event_BottomSplitter)
        self.event_Vsplitter.setStretchFactor(0, 1)
        self.event_Vsplitter.setStretchFactor(1, 3)

        splitViewTab_morph = SplitView("Morphology", "Trees")
        splitViewTab_morph.addMatPlotviewAtTop(["2D"], size=(100, 100), dpi=1000)

        splitViewTab_morph.addTablesAtBottomRight(
            ["Summary", "Distance (um)"],
            editable=True,
            sortable=False,
        )

        self.morphAnaFigs_matplotView = MatplotView() ## host for morph analysis figures
        splitViewTab_morph.bottomRight_tabview.addTab(self.morphAnaFigs_matplotView,"Figures")

        splitViewTab_morph.addParameterToLeftTab(
            "Analysis", AllMyPars.Morphor_analysis, self.morph_analysis_event
        )
        splitViewTab_morph.addParameterToLeftTab("Legend", AllMyPars.Morphor_legend)

        self.splitViewTab_morph = splitViewTab_morph
        self.visulization_view.addTab(splitViewTab_morph, splitViewTab_morph.title)


        ## add a right tool bar!
        self.frame_splitter.addWidget(self.tree_splitter)
        self.frame_splitter.addWidget(self.visulization_view)
        self.frame_splitter.setStretchFactor(0, 4)
        self.frame_splitter.setStretchFactor(1, 6)

        self.tree_splitter.setStretchFactor(0, 3)
        self.tree_splitter.setStretchFactor(1, 5)
        self.tree_splitter.setStretchFactor(2, 1)


    def event_exportViewedEvents(self):
        pv = self.event_pd_view.p.getValues()
        prePeakW = int(
            pv["Visualization choices"][1]["Pre event time (before peak)"][0]
            * self.parameters["fs"]
        )
        afterPeakW = int(
            pv["Visualization choices"][1]["Post event time (after peak)"][0]
            * self.parameters["fs"]
        )
        wlen_pre = int(0.01 * self.parameters["fs"])
        event_peakIdxs = self.events.peakIndex
        if len(event_peakIdxs) > 0:
            peakIdx = event_peakIdxs[0]
            wTime = (
                self.events.time[peakIdx - prePeakW : peakIdx + afterPeakW]
                - self.events.time[peakIdx]
            )  ## align at peak
            wdata = np.zeros((len(wTime), len(event_peakIdxs)))
            rmColsIdx = []
            for idx, peakIdx in enumerate(event_peakIdxs):
                vv = np.squeeze(
                    self.events.data_raw[peakIdx - prePeakW : peakIdx + afterPeakW]
                )
                if len(vv) == wdata.shape[0]:
                    wdata[:, idx] = vv
                else:
                    rmColsIdx.append(idx)  ## remove non-valid columns
            if len(rmColsIdx) > 0:
                wdata = np.delete(wdata, np.s_[rmColsIdx], 1)
            return wdata, wTime
        else:
            return [], []

    def event_viewSingleEvents(
        self, upStrokeIdx, onsetIdx, peakIdx, donwstrokeIdx, noMove=False
    ):
        pv = self.event_pd_view.p.getValues()
        prePeakW = int(
            pv["Visualization choices"][1]["Pre event time (before peak)"][0]
            * self.parameters["fs"]
        )
        afterPeakW = int(
            pv["Visualization choices"][1]["Post event time (after peak)"][0]
            * self.parameters["fs"]
        )

        wlen_pre = int(0.01 * self.parameters["fs"])
        # wlen_post = int(0.10 * self.parameters['fs']) ## show larger window to give context
        eventTime = self.events.time[peakIdx]
        peakHeight = self.events.data[peakIdx]
        eventOnsetTime = self.events.time[onsetIdx]
        wTime = self.events.time[onsetIdx - prePeakW : peakIdx + afterPeakW]
        wdata = self.events.data[onsetIdx - prePeakW : peakIdx + afterPeakW]
        wTime2 = self.events.time[onsetIdx - wlen_pre : donwstrokeIdx]
        wdata2 = self.events.data[onsetIdx - wlen_pre : donwstrokeIdx]
        self.eventMarkerLine_v.setPos(eventTime)
        self.eventMarkerLine_h.setPos(peakHeight)
        plotHandle = self.event_traceView2.getItem(0, 0)
        plotHandle.clear()
        plotHandle.plot(wTime - eventOnsetTime, wdata)
        plotHandle.plot(wTime2 - eventOnsetTime, wdata2, pen=pg.mkPen("r", width=2))
        # plotHandle.setLabels(bottom=('Time', traceXUnit), left=('Voltage', traceYUnit))
        plotHandle.autoRange()

    def eventTable_enterAction(self, row1, col1, row0, col0):
        """
        update event trace plot and template plot!

        """
        item = self.eventTable.item(row1, 4)  ## onset index
        if item == None:
            return
        onsetIdx = item.value
        peakIdx = self.eventTable.item(row1, 5).value  ## peak index
        downStrokeIdx = self.eventTable.item(row1, 6).value
        # peakVal = self.eventTable.item(row1, 8).value ## peak height
        donwstrokeIdx = peakIdx + self.events.downstroke_samples
        upstrokeIdx = peakIdx - self.events.upstroke_samples
        self.event_viewSingleEvents(upstrokeIdx, onsetIdx, peakIdx, donwstrokeIdx)

    def switchStateLine(self, l, manul="False"):
        if l.active:
            # self.clearMarkers()
            # self.addMarker('v')
            l.setPen((155, 155, 155, 255))
            l.active = False
            print("Event invalidated")
        else:
            l.setPen("g")
            l.active = True
            print("Event recovered")

    def eventTable_toggleState(self, row, column):
        ## toogle inUse state!
        InUseColumn = 11
        if column == InUseColumn:  ## for the InUse column
            item = self.eventTable.takeItem(row, column)
            i = item.index
            v = not (item.value)
            newItem = pg.widgets.TableWidget.TableWidgetItem(v, i)
            self.eventTable.setItem(row, column, newItem)
            if newItem.value:
                print(f"Row {row} activated!")
            else:
                print(f"Row {row} deactivated!")
            lineIdx = self.events.eventTableDF.index[row]
            linetype = self.eventTable.item(row, 3).value
            self.updateEventMarker(linetype, lineIdx, newItem.value)

    def eventTable_doubleClickAction(self, row, column):
        self.eventTable_toggleState(row, column)

    def updateEventMarker(self, ltype, lineIdx, lval):
        if self.events.sweepType == "single":
            lauto, lman = self.getEventLine()
            if ltype == "auto":
                self.switchStateLine(lauto[lineIdx], manul="False")
            else:
                self.switchStateLine(lman[lineIdx], manul="False")
        else:
            if self.events.eventMarker != None:
                scatterPlot = self.events.eventMarker
                penColor = (
                    scatterPlot.getSpotOpts(scatterPlot.data[lineIdx])[2]
                    .color()
                    .getRgb()
                )
                if penColor[0] > 0:
                    pc = "g"  ## change r to g
                else:
                    pc = "r"  ## change g to r
                spotsPen = [scatterPlot.getSpotOpts(yy)[2] for yy in scatterPlot.data]
                spotsBrush = [scatterPlot.getSpotOpts(yy)[3] for yy in scatterPlot.data]
                spotsPen[lineIdx] = pg.mkPen(pc)
                spotsBrush[lineIdx] = pg.mkBrush(pc)
                scatterPlot.setPen(spotsPen)
                scatterPlot.setBrush(spotsBrush)

    def updateROI_table(self):
        data = np.array(
            self.roidata, dtype=[("ChanID", object), ("X", object), ("Y", object)]
        )
        self.sliceView.tables["ROI list"].clear()
        self.sliceView.tables["ROI list"].setData(data)  ## update QtTable widget

    def updateROIlabel(self, tlabel, pos, idx):
        tlabel.setPos(pos.x(), pos.y() - 10)
        for j, roi_idx in enumerate(self.roidata):
            if idx == roi_idx[0]:
                self.roidata.pop(j)
                self.roidata.append((idx, pos.x(), pos.y()))
                self.updateROI_table()
                return
        self.roidata.append((idx, pos.x(), pos.y()))
        self.updateROI_table()

    def showDialog(self):
        text, ok = QtWidgets.QInputDialog.getText(
            self, "Input Dialog", "Enter channel Name:"
        )
        if ok:
            self.currentROI_label = str(text)
            for roi_idx in self.roidata:
                if text == roi_idx[0]:
                    print("Duplicate channel label!")
                    self.currentROI_label = str("")  ## prevent duplicate chan label
                    return
        else:
            self.currentROI_label = str("")

    def redrawAllROIs(self):
        for roiID, i, j in self.roidata:
            roiCenter = [i - 5, j - 5]
            newroi = cellROI(
                self.sliceView.slice_view,
                self,
                roiCenter,
                [10, 10],
                label=roiID,
                pen=(1, 20),
                removable=True,
            )
            Tlabel = self.addTextlabel(roiID, roiCenter)

            newroi.TextLabel = Tlabel
            self.sliceView.slice_view.addItem(newroi.TextLabel)  ## add text label
            self.sliceView.slice_view.addItem(newroi)  ## add ROI handles

            self.slice_viewROIs.append(newroi)
            self.slice_viewROIsText.append(Tlabel)

    def addTextlabel(self, label, pos):
        TextLabel = pg.TextItem(text=label, color=(255, 255, 255), anchor=(0, 0))
        TextLabel.setFont(QtGui.QFont("serif", 14))
        TextLabel.setPos(pos[0], pos[1] - 10)
        return TextLabel

    def updateTreeMorphView(self, fname):
        # try:
        pv = self.splitViewTab_morph.getParTreePars("Analysis")
        if pv["Options"][1]["Draw contour"][0]:
            drawContour = True
        else:
            drawContour = False
        if pv["Options"][1]["Ignore diameters"][0]:
            realDia = True
        else:
            realDia= False
        self.neuronMorph =  None 
        neurons = morphor_nm.load_morphology(fname, somaType = SomaType.SOMA_CYLINDERS)
        fig2D = self.splitViewTab_morph.matplotViews["2D"].getFigure()
        fig2D.set_size_inches(
            5, 5, forward=False)
        fig2D.clf()
        ax2D = fig2D.add_subplot(111)
        ax2D.cla()
        if isinstance(neurons, MultiSoma):
            print('Multi soma detected!')
            ## population of neurons with soma only
            centers = {
                "Name": [],
                "X": [],
                "Y": [],
                "Z": [],
                "average radius": [],
                "maximal radius": [],
                "minimal radius": [],
                "min_max_ratio": [],
                "maximal diameter": [],
                "intercell distance": [],
            }

            for idx, n in enumerate(neurons.pop_neurons):
                maxDia, soma_center, soma_radius, soma_avgRadius = getSomaStats(n)
                centers["X"].append(soma_center[0])
                centers["Y"].append(soma_center[1])
                centers["Z"].append(soma_center[2])
                centers["average radius"].append(soma_avgRadius)
                centers["maximal radius"].append(np.max(soma_radius))
                centers["minimal radius"].append(np.min(soma_radius))
                centers["min_max_ratio"].append(np.min(soma_radius) / np.max(soma_radius))
                centers["maximal diameter"].append(maxDia)
                centers["intercell distance"].append(neurons.getDistMatrix(centers))
                centers["Name"].append(str(idx + 1))

            morphor_viewer.draw(
                neurons, mode="2d", fig=fig2D, ax=ax2D,
                contour_on=drawContour, contour_color='g', contour_linewidth=0.2,
                 labels=None
            )
            self.splitViewTab_morph.matplotViews["2D"].draw()
            xlim1 = ax2D.xaxis.get_data_interval().copy()
            ylim1 = ax2D.yaxis.get_data_interval().copy()

            if len(neurons.custom_data) > 0:
                for datablock in neurons.custom_data:
                    d = datablock.data_block
                    if d[0][4] == 8:  ## PIA, just draw lines
                        ax2D.plot(d[:, 0], d[:, 1], "k", label="Pia")
                        self.pia = d
            else:
                self.pia = None

            xlim2 = ax2D.xaxis.get_data_interval()
            ylim2 = ax2D.yaxis.get_data_interval()
            xlim0 = np.concatenate([xlim1, xlim2])
            ylim0 = np.concatenate([ylim1, ylim2])
            ax2D.set_title("{0}".format(fname), fontsize=12, color="k")

            df = pd.DataFrame.from_dict(centers, orient="index").transpose()
            self.neuronsData = df  ## store neurons data
            self.morphorAxes2D = ax2D

            df = np.array(
                df.to_records(index=False)
            )  ## format dataframe for using in QtTable wiget
            self.splitViewTab_morph.tables["Summary"].clear()
            self.splitViewTab_morph.tables["Summary"].appendData(df)
            self.splitViewTab_morph.tables["Summary"].show()
    
            self.updateInterCellDistance()
            self.neuronMorph =  None ## no dendrites for further analysis
        else:  ## single neuron with dendtrite and/or axons
            fig2D = self.splitViewTab_morph.matplotViews["2D"].getFigure()
            fig2D.clf()
            widths = [3,1]
            # heights = [5,1]
            gs0 = fig2D.add_gridspec(1, 2, width_ratios = widths)
            ax2D = fig2D.add_subplot(gs0[0,0])
            ax2D.cla()
            pv = self.splitViewTab_morph.getParTreePars("Analysis")
            rotateAngle = pv["Options"][1]["Rotate tree (degree)"][0]
            step_size=pv['Options'][1]['Bin size (um)'][0]
            smoothBins=pv['Options'][1]['Gaussian window size (num. bins)'][0]
            smoothStandardDeviation=pv['Options'][1]['Std of Gaussian kernel (num. bins)'][0]
            if rotateAngle!=0:
                neurons = rotate(neurons, [0,0,1], rotateAngle*np.pi/180)
            self.neuronMorph = neurons
            morphor_viewer.draw(
                neurons, mode="2d", realistic_diameters=realDia, fig=fig2D, ax=ax2D,
                contour_on=drawContour, contour_color='g', contour_linewidth=0.2,rotationContour=rotateAngle
            )
            ax2D.axis("equal")  ## only works for 2D
            # ax2D.view_init(azim=0, elev=90)
            ax2D.set_title("{0}".format(fname), fontsize=12, color="k")

            # ax3D.set_zlim(zmin=max_scaling[0], zmax=max_scaling[1])
            print("Neuron morphology loaded!")
            df_summary = {}
            df_summary["ASC file"] = [fname]
            df_summary = extractMorphFeatures(neurons, df_summary)
            df_summary = {k:[df_summary[k]] for k in df_summary} ## weird requr. of pandas
            df_summary  = pd.DataFrame.from_dict(df_summary, orient="index").transpose()
            df_summary  = np.array(
                df_summary.to_records(index=False)
            )  ## format dataframe for using in QtTable widget
            self.splitViewTab_morph.tables["Summary"].clear()
            self.splitViewTab_morph.tables["Summary"].appendData(df_summary)
            self.splitViewTab_morph.tables["Summary"].show()
            ax_densityX = fig2D.add_subplot(gs0[0,1], sharey = ax2D)
            ax_densityX.cla()
            self.update_density(axis='y', step_size=step_size, ax = ax_densityX,
             flipAxes=True, addTitle=False,smoothBins=smoothBins,smoothStandardDeviation=smoothStandardDeviation)
            ax_densityX.xaxis.set_ticks_position('top')
            ax_densityX.xaxis.set_label_position('top') 
            ax_densityX.tick_params(labelbottom=False,labeltop=True)
            ax_densityX.spines["right"].set_visible(False)
            ax_densityX.spines["bottom"].set_visible(False)
            fig2D.subplots_adjust(top=0.904, bottom=0.04, left=0.0, right=0.90, hspace=0.0, wspace=0.254)

        ax2D.axis("off")
        self.splitViewTab_morph.matplotViews['2D'].draw()
        self.splitViewTab_morph.matplotViews["2D"].canvas.draw()
        fig2D.tight_layout()
        self.splitViewTab_morph.bottomRight_tabview.setCurrentIndex(0)
        self.splitViewTab_morph.bottomRight_tabview.setCurrentIndex(2)
        # self.splitViewTab_morph.matplotViews['3D'].draw()

        self.splitViewTab_morph.show()

    def minmaxDist(self, ps):
        """calculate minmial and maximal diameter
        ps: 2-D numpy array
        """
        nPoints = len(ps)
        maxDia = 0
        for j in range(nPoints - 1):
            for k in range(j + 1, nPoints):
                diameter = np.sqrt(np.sum((ps[j] - ps[k]) ** 2))
                if diameter > maxDia:
                    maxDia = diameter
        return maxDia

    def update_2D_polar_density(self, ax=None, addTitle =True, angle_step=np.pi/15, neurite_type = NeuriteType.all, step_size=5):
            if self.neuronMorph is not None:
                if ax is None:
                    fig = self.morphAnaFigs_matplotView.getFigure()
                    fig.clf()
                    ax = fig.add_subplot(111,projection="polar")
                    ax.cla()                       
                hist, A, R = sholl_polar(self.neuronMorph, step_size=step_size, angle_step=angle_step)
                pc = ax.pcolormesh(A, R, hist.T, cmap="magma_r")
                fig.colorbar(pc, orientation='vertical')
                if addTitle:
                    ax.set(title= "XY plane density")
                self.morphAnaFigs_matplotView.figure.subplots_adjust(top=0.861,bottom=0.02,left=0.0,right=0.7,hspace=0.2,wspace=0.0)
                self.morphAnaFigs_matplotView.draw()
                self.morphAnaFigs_matplotView.canvas.draw()
            else:
                print('Not a morphological object for sholl analysis')
            self.splitViewTab_morph.bottomRight_tabview.setCurrentIndex(2)
            
    def update_2D_density(self, ax=None, addTitle =True, neurite_type = NeuriteType.all,
     step_size=5,  useFullRange=True, showColorbar=False, showAxisValues=False,smoothBins=11,
                smoothStandardDeviation=2):
            if self.neuronMorph is not None:
                if ax is None:
                    fig = self.morphAnaFigs_matplotView.getFigure()
                    fig.clf()
                axes = []
                images = []
                for ntype in [NeuriteType.all, NeuriteType.axon, NeuriteType.basal_dendrite, NeuriteType.apical_dendrite, 'All dendrite']:
                    if ntype == NeuriteType.all:
                        ax = fig.add_subplot(151, aspect='equal')
                        ntypeName = "All neurites"
                    elif ntype == NeuriteType.axon:
                        ax = fig.add_subplot(152, aspect='equal')
                        ntypeName = "Axon"                
                    elif ntype == NeuriteType.basal_dendrite:
                        ax = fig.add_subplot(153, aspect='equal')
                        ntypeName = "Basla dendrite"
                    elif ntype == NeuriteType.apical_dendrite:
                        ax = fig.add_subplot(154, aspect='equal')
                        ntypeName = "Apical dendrite"
                    else:
                        ax = fig.add_subplot(155, aspect='equal')
                        ntypeName = "All dendrite"                       
                    axes.append(ax)
                    if ntype != 'All dendrite':                     
                        d2d, xedges, yedges, centerH, centerV = sholl_2D_density(self.neuronMorph, step_size=step_size,
                         neurite_type=ntype,useFullRange=useFullRange)
                    else:
                        d2d, xedges, yedges, centerH, centerV = sholl_2D_density(self.neuronMorph, step_size=step_size,
                         neurite_type=[NeuriteType.basal_dendrite, NeuriteType.apical_dendrite], useFullRange=useFullRange)
                    ax.cla() 

                    if len(d2d) > 0 :
                        d2d = smooth2D(d2d, smoothBins, smoothStandardDeviation)
                        xedges -=centerH
                        yedges -=centerV
                        im = ax.imshow(np.transpose(d2d), extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],interpolation='bilinear', 
                        origin='lower')
                        images.append(im)
                        soma = MATPLT.Circle((0, 0), 5, color='r')
                        ax.add_patch(soma)
                        if addTitle:
                            ax.set(title= f"{ntypeName} density")
                        if showColorbar:
                            fig.colorbar(im, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
                    else:
                        ax.axis('off')
                    if not showAxisValues:
                        ax.axis('off')
                if showAxisValues:
                    axes[0].set_xlabel('X (um)')
                    axes[0].set_ylabel('Y (um)')

                fig.tight_layout()
                self.morphAnaFigs_matplotView.draw()
                self.morphAnaFigs_matplotView.canvas.draw()

            else:
                print('Not a morphological object for sholl analysis')
            self.splitViewTab_morph.bottomRight_tabview.setCurrentIndex(2)

    def update_density(self, axis='y', ax=None,step_size=1,
     flipAxes=False, addTitle =True, neurite_type = NeuriteType.all,
     smoothBins=11,smoothStandardDeviation=2):
        if self.neuronMorph is not None:
            if ax is None:
                fig = self.morphAnaFigs_matplotView.getFigure()
                fig.clf()
                ax = fig.add_subplot(111)
                ax.cla()
            else:
                fig = ''

            if neurite_type == NeuriteType.all:
                neurite_type  = [NeuriteType.apical_dendrite, NeuriteType.basal_dendrite, NeuriteType.axon]
            else:
                neurite_type  = [neurite_type]
            for neurite_type in neurite_type :
                c = NeutriteColors[neurite_type]
                density, bin_edges, centerVal, _ = sholl_single_axis(self.neuronMorph, step_size=step_size, axis=axis,
                 neurite_type=neurite_type)
                if density==[]:
                    continue
                density = smooth(density, window_len=smoothBins, std=smoothStandardDeviation,
                window='gaussian')
                ndiff = len(bin_edges) - len(density)
                if ndiff > 0:
                    bin_edges = bin_edges[ndiff:]
                else:
                    density = density[-ndiff:]
                if flipAxes:
                    ax.plot(density, bin_edges-centerVal, color=c)
                    xlabel = "um per bin"
                    ylabel = 'Distance from soma (um)'
                else:
                    ax.plot(bin_edges, density-centerVal, color=c)
                    ylabel = "um per bin"
                    xlabel = 'Distance from soma (um)'

            ax.set(xlabel=xlabel, ylabel=ylabel)
            if addTitle:
                ax.set(title= axis+" density")
            self.morphAnaFigs_matplotView.draw()
            self.morphAnaFigs_matplotView.canvas.draw()  
            if fig != '':
                fig.tight_layout()
        else:
            print('Not a morphological object for sholl analysis')
        self.splitViewTab_morph.bottomRight_tabview.setCurrentIndex(2)


    def update_sholl(self, step_size=1, smoothBins=11,smoothStandardDeviation=2):
        if self.neuronMorph is not None:        
            fig = self.morphAnaFigs_matplotView.getFigure()
            fig.clf()
            ax = fig.add_subplot(111)
            ax.cla()
            sns.set_style("whitegrid")
            for neurite in self.neuronMorph.neurites:
                sholl_dist, sholl_bins = sholl_analysis(self.neuronMorph, step_size=step_size, neurite_type=neurite.type)
                c = NeutriteColors[neurite.type]
                sholl_bins = smooth(sholl_bins, window_len=smoothBins, std=smoothStandardDeviation,
                window='gaussian')
                ndiff = len(sholl_bins) - len(sholl_dist)
                if ndiff > 0:
                    sholl_bins = sholl_bins[ndiff:]
                else:
                    sholl_dist = sholl_dist[-ndiff:]
                ax = sns.lineplot(x=sholl_bins, y=sholl_dist, ax=ax, color=c)
            ax.set(xlabel='Distance from soma (um)', ylabel="Num. points",\
            title="Sholl frequency")
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)
            fig.tight_layout()
            self.morphAnaFigs_matplotView.draw()
            self.morphAnaFigs_matplotView.canvas.draw()

        else:
            print('Not a morphological object for sholl analysis')
        self.splitViewTab_morph.bottomRight_tabview.setCurrentIndex(2)

    def updateSliceView(self):
        if hasattr(self, "slice_viewROIs"):
            for roi in self.slice_viewROIs:
                roi.remove()

        self.currentROI_label = ""
        self.selectedFilesTabView.setCurrentIndex(3)
        self.visulization_view.setCurrentIndex(3)  ## index for slice_view is 6

        self.roidata = []
        self.sliceView.tables["ROI list"].clear()

        dialog = QtWidgets.QFileDialog(self)
        dialog.setWindowTitle("Import Slice image")
        dialog.setNameFilter("Image files (*.png *.jpg)")
        if hasattr(self, "cwd"):
            dialog.setDirectory(self.cwd)
        else:
            dialog.setDirectory(QtCore.QDir.currentPath())
        dialog.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        if dialog.exec_() == QtWidgets.QDialog.Accepted:
            fname = str(dialog.selectedFiles()[0])
        else:
            return

        img = ReadImage.imread(fname)
        self.sliceView.slice_view.clear()
        label = pg.TextItem(text="", color=(200, 0, 0), anchor=(0, 0))
        label.setFont(QtGui.QFont("serif", 12))

        self.sliceView.slice_view.addItem(label)
        ## Display the data and assign each frame a time value from 1.0 to 3.0
        self.sliceView.slice_view.setImage(
            img,
            pos=[0, 0],
            axes={"x": 1, "y": 0, "c": 2},
            xvals=np.linspace(1, 3, img.shape[0]),
        )
        # self.sliceView.slice_view.setCurrentIndex(2)
        vLine = pg.InfiniteLine(angle=90, movable=False)
        hLine = pg.InfiniteLine(angle=0, movable=False)
        self.sliceView.slice_view.addItem(vLine, ignoreBounds=True)
        self.sliceView.slice_view.addItem(hLine, ignoreBounds=True)
        vb = self.sliceView.slice_view.view
        vb.layout.setContentsMargins(0, 0, 0, 0)

        def mouseMoved(evt):
            if self.sliceView.slice_view.view.sceneBoundingRect().contains(
                evt.x(), evt.y()
            ):
                mousePoint = vb.mapToView(evt)
                label.setText("x=%0.0f, y=%0.0f" % (mousePoint.x(), mousePoint.y()))
                vLine.setPos(mousePoint.x())
                hLine.setPos(mousePoint.y())
                label.setPos(mousePoint.x(), mousePoint.y())

        def mouseClicked(evt):
            if (
                evt.modifiers() == QtCore.Qt.ControlModifier
            ):  ## press control key to add ROI. to avoid accidently clicking
                pos1 = evt.scenePos()
                if self.sliceView.slice_view.view.sceneBoundingRect().contains(
                    pos1.x(), pos1.y()
                ):
                    mousePoint = vb.mapToView(pos1)
                    self.showDialog()
                    roiID = self.currentROI_label  # str(len(self.roidata)+1)
                    if roiID:
                        roiCenter = [mousePoint.x() - 5, mousePoint.y() - 5]
                        # self.roidata.append((roiID, roiCenter[0], roiCenter[1]))  ## addd item to roidata dictionary
                        newroi = cellROI(
                            self.sliceView.slice_view,
                            self,
                            roiCenter,
                            [10, 10],
                            label=roiID,
                            pen=(1, 20),
                            removable=True,
                        )
                        Tlabel = self.addTextlabel(roiID, roiCenter)

                        newroi.TextLabel = Tlabel
                        self.sliceView.slice_view.addItem(
                            newroi.TextLabel
                        )  ## add text label
                        self.sliceView.slice_view.addItem(newroi)  ## add ROI handles
                        self.updateROIlabel(Tlabel, mousePoint, roiID)
                    # self.updateROI_table()

        self.sliceView.slice_view.scene.sigMouseMoved.connect(mouseMoved)
        self.sliceView.slice_view.scene.sigMouseClicked.connect(mouseClicked)
        ## Set a custom color map
        colors = [(0, 0, 0), (0, 0, 255), (255, 255, 255)]
        cmap = pg.ColorMap(pos=[0.3, 0.5, 0.9], color=colors)
        self.sliceView.slice_view.setColorMap(cmap)
        self.sliceView.slice_view.autoRange()

    def updateCurrentPulseTree(self):
        self.currentPulseTree = self.PulseTrees[self.trees_view.currentIndex()]
        self.updateStatusBar_fileName()
        self.hideSpikePanels()

    def add_toolbar(self):
        global_icon = os.path.join(patchview_dir, "Data", "icons", "settings.png")
        rect_icon = os.path.join(patchview_dir, "Data", "icons", "rectangle.png")
        spk_icon = os.path.join(patchview_dir, "Data", "icons", "spikes.png")
        connection2_icon = os.path.join(patchview_dir, "Data", "icons", "settings.png")
        neuron_icon = os.path.join(patchview_dir, "Data", "icons", "neuron.png")

        self.toolbar = QtWidgets.QToolBar("Main toolbar")
        self.addToolBar(
            2, self.toolbar
        )  # https://doc.qt.io/qt-5/qt.html#ToolBarArea-enum
 
        self.toogleMouseModePanAction = pg.QtGui.QAction(
            pg.QtGui.QIcon(rect_icon), "Mouse Mode"
        )
        self.toogleMouseModePanAction.setShortcut("Ctrl+m")
        self.toogleMouseModePanAction.setStatusTip("switch between pan and rect mode")
        self.toogleMouseModePanAction.triggered.connect(self.MouseMode_clicked)
        self.toolbar.addAction(self.toogleMouseModePanAction)

        # self.toogleDetectSpikeAction = pg.QtGui.QAction(
        #     pg.QtGui.QIcon(spk_icon), "Firing pattern"
        # )
        # self.toogleDetectSpikeAction.setShortcut("Alt+a")
        # self.toogleDetectSpikeAction.setStatusTip("Firing pattern")
        # self.toogleDetectSpikeAction.triggered.connect(self.detectSpikes_clicked)
        # self.toolbar.addAction(self.toogleDetectSpikeAction)

        # self.exportAction = pg.QtGui.QAction(
        #     pg.QtGui.QIcon(connection2_icon), "Connectivtiy"
        # )
        # self.exportAction.setShortcut("Alt+e")
        # self.exportAction.setStatusTip("Connectivtiy")
        # self.exportAction.triggered.connect(self.exportConnectionFig_clicked)
        # self.toolbar.addAction(self.exportAction)

        self.importSlicetAction = pg.QtGui.QAction(
            pg.QtGui.QIcon(neuron_icon), "Import slice image"
        )
        # self.exportAction.setShortcut("Alt+e")
        self.importSlicetAction.setStatusTip("import slice image")
        self.importSlicetAction.triggered.connect(self.importSlice_clicked)
        self.toolbar.addAction(self.importSlicetAction)

        gp1_icon = os.path.join(patchview_dir, "Data", "icons", "GP1.png")
        self.gpAction = pg.QtGui.QAction(pg.QtGui.QIcon(gp1_icon), "Coupling ratio")
        # self.exportAction.setShortcut("Alt+e")
        self.gpAction.setStatusTip("Analysis gap junction")
        self.gpAction.triggered.connect(self.gp_clicked)
        self.toolbar.addAction(self.gpAction)

        GPcorr_icon = os.path.join(patchview_dir, "Data", "icons", "GPcorr.png")
        self.gpAction2 = pg.QtGui.QAction(pg.QtGui.QIcon(GPcorr_icon), "correlations")
        # self.exportAction.setShortcut("Alt+e")
        self.gpAction2.setStatusTip("Analysis gap junction")
        self.gpAction2.triggered.connect(self.gpc_clicked)
        self.toolbar.addAction(self.gpAction2)

    def importTree_clicked(self):        
        dialog = QtWidgets.QFileDialog(self)
        dialog.setWindowTitle("Import Slice image")
        dialog.setNameFilter("Image files (*.asc *.ASC)")
        if hasattr(self, "cwd"):
            dialog.setDirectory(self.cwd)
        else:
            dialog.setDirectory(QtCore.QDir.currentPath())
        dialog.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        if dialog.exec_() == QtWidgets.QDialog.Accepted:
            fname = str(dialog.selectedFiles()[0])
            self.cwd = os.path.dirname(os.path.abspath(fname))
        else:
            return
        self.prepareTree(fname)
        self.currentAnMorphView = ''

    def prepareTree(self, fname):
        self.visulization_view.setCurrentIndex(5)
        self.currentMorphTreeFile = fname
        if glob.glob(fname[:-4] + "_mod.ASC") == []:
            fname = cleanASCfile(fname)  ## to clean not-want sections
        self.splitViewTab_morph.matplotViews["2D"].getFigure().clf()
        self.morphAnaFigs_matplotView.getFigure().clf()
        self.updateTreeMorphView(fname)
        print('morphfile', fname, fname[-8:])
        if fname[-8:] == "_mod.ASC":
            os.remove(fname)


    def importSlice_clicked(self):
        self.updateSliceView()

    def MouseMode_clicked(self):
        rect_icon = os.path.join(patchview_dir, "Data", "icons", "rectangle.png")
        nav_icon = os.path.join(patchview_dir, "Data", "icons", "navigation.png")
        if self.MouseMode == pg.ViewBox.RectMode:
            self.MouseMode = pg.ViewBox.PanMode
            self.toogleMouseModePanAction.setIcon(pg.QtGui.QIcon(nav_icon))
            self.toogleMouseModePanAction.setStatusTip("Your Mouse is in pan mode")
        else:
            self.MouseMode = pg.ViewBox.RectMode
            self.toogleMouseModePanAction.setIcon(pg.QtGui.QIcon(rect_icon))
            self.toogleMouseModePanAction.setStatusTip(
                "Your Mouse is in selection mode"
            )

        self.setViewboxMouseMode()

    def setViewboxMouseMode(self):
        for plt in self.plotItemList:
            vb = plt.getViewBox()
            vb.setMouseMode(self.MouseMode)  ## set as rect mode

    def add_menubar(self):
        """
        set up manu bar for the main window
        """
        self.mbar = QMenuBar()
        self.MainLayout.addWidget(self.mbar, 0, 0)

        self.fileMenu = self.mbar.addMenu("&File")
        self.saveFileAction = QAction("&Export .plk")
        #        self.saveFileAction.setShortcut("Ctrl+s")
        self.saveFileAction.setStatusTip("Export data to pickle object")
        self.saveFileAction.triggered.connect(self.save_clicked)
        self.fileMenu.addAction(self.saveFileAction)

        self.saveFileAction2 = QAction("&Export .NWB")
        self.saveFileAction2.setShortcut("Ctrl+s")
        self.saveFileAction2.setStatusTip("Export data NWB format")
        self.saveFileAction2.triggered.connect(self.saveNWB_clicked)
        self.fileMenu.addAction(self.saveFileAction2)

        self.fileMenu.addSeparator()

        self.resetAction = QAction("&Reset")
        #        self.saveFileAction.setShortcut("Ctrl+s")
        self.resetAction.setStatusTip("Reset all wigets")
        self.resetAction.triggered.connect(self.resetAll_clicked)
        self.fileMenu.addAction(self.resetAction)

        self.exitAction = QAction("&Exit")
        self.exitAction.triggered.connect(self.exit_clicked)
        self.fileMenu.addAction(self.exitAction)

        ## TODO
        # self.AnalysisMenu = self.mbar.addMenu("&Analysis")
        # self.AnaAction4 = pg.QtGui.QAction("&Gap junction (Batch)")
        # # self.AnaAction4.triggered.connect(self.checkGapJunctionActionBatch_clicked)
        # self.AnaAction4.setStatusTip("Check gap connections in selected series")
        # self.AnalysisMenu.addAction(self.AnaAction4)

        self.VisualMenu = self.mbar.addMenu("&Visualization")
        self.VisualAction2 = QAction("&Show averaged traces")
        self.VisualAction2.setShortcut("Ctrl+m")
        self.VisualAction2.setStatusTip("Plot averaged stimulated response")
        self.VisualAction2.triggered.connect(self.VisualAction2_clicked)
        self.VisualMenu.addAction(self.VisualAction2)

        self.VisualAction3 = QAction("&Show aligned spikes")
        self.VisualAction3.setShortcut("Alt+s")
        self.VisualAction3.setStatusTip("Plot aligned spikes during each sweep")
        self.VisualAction3.triggered.connect(
            self.plotSpikes_clicked
        )  ## fine tune firing pattern plotting with user define parameters
        self.VisualMenu.addAction(self.VisualAction3)

        self.VisualAction4 = QAction("&Show concatenated trace")
        self.VisualAction4.setStatusTip("Plot concatenated trace")
        self.VisualAction4.triggered.connect(
            self.choose2PlotSP
        )  ## fine tune firing pattern plotting with user define parameters
        self.VisualMenu.addAction(self.VisualAction4)


        self.OptionMenu = self.mbar.addMenu("&Option")
        self.OptionAction2 = QAction("&Switch background color")
        # OptionAction1.setShortcut("Ctrl+a")
        self.OptionAction2.triggered.connect(self.switchBackground_clicked)
        self.OptionMenu.addAction(self.OptionAction2)

        self.OptionAction3 = QAction("&Spike panels", self, checkable=True)
        # OptionAction1.setShortcut("Ctrl+a")
        self.OptionAction3.triggered.connect(self.hideSpikePanels_clicked)
        self.OptionMenu.addAction(self.OptionAction3)
        self.OptionAction3.setChecked(False)  # .isChecked()

        # self.OptionAction4 = QAction(
        #     "&Show 3D firing pattern", self, checkable=True
        # )
        # self.OptionAction4.triggered.connect(self.FiringPattern3D_clicked)
        # self.OptionMenu.addAction(self.OptionAction4)
        # self.OptionAction4.setChecked(False)

        self.OptionAction5 = QAction(
            "&Remove stimuli artifacts", self, checkable=True
        )
        self.OptionAction5.triggered.connect(self.removeStimArtifacts_clicked)
        self.OptionAction5.setChecked(False)
        self.removeArtifact = 1
        self.OptionMenu.addAction(self.OptionAction5)

        self.dv2dt2 = False  ## by default, do not calucate second order
        self.OptionAction6 = QAction(
            "&dv^2/dt^2 vs dv/dt", self, checkable=True
        )
        # OptionAction1.setShortcut("Ctrl+a")
        self.OptionAction6.triggered.connect(self.dvdt2_clicked)
        self.OptionMenu.addAction(self.OptionAction6)
        self.OptionAction6.setChecked(False)  # .isChecked()

        self.batchFPana = False  ## by default, do not plot when batch processing
        self.OptionAction7 = QAction(
            "Show plots for firing pattern batch analysis", self, checkable=True
        )
        # OptionAction1.setShortcut("Ctrl+a")
        self.OptionAction7.triggered.connect(
            self.batchFPana_clicked
        )  # self.dvdt2_clicked
        self.OptionMenu.addAction(self.OptionAction7)
        self.OptionAction7.setChecked(self.batchFPana)

        # self.OptionAction8= pg.QtGui.QAction(
        #     "Configuration", self, checkable=False
        # )
        # # OptionAction1.setShortcut("Ctrl+a")
        # self.OptionAction8.triggered.connect(
        #     self.pvSetting_clicked
        # )  # self.dvdt2_clicked
        # self.OptionMenu.addAction(self.OptionAction8)

        self.HelpMenu = self.mbar.addMenu("&About")
        self.HelpAction1 = QAction("&Documentation")
        self.HelpAction1.triggered.connect(self.link2Doc_clicked)
        self.HelpMenu.addAction(self.HelpAction1)

        self.HelpAction3 = QAction("&LICENSE")
        self.HelpAction3.setStatusTip("BSD-3")
        self.HelpAction3.triggered.connect(self.License_clicked)
        self.HelpMenu.addAction(self.HelpAction3)

    def pvSetting_clicked(self):
        DATA_PATH = os.path.join(patchview_dir, "Data", "uis","batchDownloadOptionsDialog.ui")
        WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(DATA_PATH)
        self.bd_dialog  = TemplateBaseClass()
        self.bd_dialog_form = WindowTemplate()
        self.bd_dialog_form.setupUi(self.bd_dialog)
        self.bd_dialog.exec_()

    def updateEvents(self):
        self.events.updateCurrentSweep()
        self.event_showTemplate(False)

    def addEventLinePair(self, tracePlot, evt):
        """TODO
        Add back event line pairs from saved events if there is any

        Returns
        -------
        None.

        """
        if self.events.eventTimeStamps:
            if self.events.currentSweep in self.events.eventTimeStamps.keys():
                eventList = self.events.eventTimeStamps.get(self.events.currentSweep)
                # self.events.eventTimeStamps.update({self.currentSweep: [startPos, endPos]})

                # offsetTime = self.eventParTree_data_view.p.getValues()['Template'][1]['Window length'][0]
                for e in eventList:
                    vlineL = infiniteLinePair(
                        parent=tracePlot,
                        offset=e[1] - e[0],
                        span=(0, 0.8),
                        pos=e[0],
                        index=None,
                        Frame=self,
                        angle=90,
                        movable=True,
                        pen="r",
                        hoverPen="w",
                        name="template1",
                    )
                    vlineL.mouseClickEvent(evt)
                    vlineL.twin.sigPositionChangeFinished.connect(
                        vlineL.mouseDragEventTwin
                    )
                    # if self.events.currentSweep in self.events.templateTwinLine.keys():
                    #     self.events.templateTwinLine[self.events.currentSweep].append(vlineL)
                    # else:
                    #     self.events.templateTwinLine[self.events.currentSweep]=[vlineL]
                    tracePlot.addItem(vlineL, ignoreBounds=True)
                    tracePlot.addItem(vlineL.twin, ignoreBounds=True)
                    self.events.templateTwinLine.append(vlineL)
                    self.events.eventLineReloaded = True  ## turn off reloading flag
                    self.events.tracePlot = tracePlot

    def checkNodeIndex(self):
        selected = self.selTreesSP.selectedItems()
        if selected:
            sel = selected[0]
            sel_index = [int(s) for s in sel.text(3).split(" ")]
            if len(sel_index) != 2:
                self.showdialog("choose the series that you want to analysis")
                return ["", ""]
            else:
                return [sel.text(1), sel_index]
        else:  ## it may come from selection area!
            if self.currentPulseTree.dat_file:
                index = self.currentPulseTree.selectedItems()[0].index
                if len(index) != 2:
                    self.showdialog("choose the series that you want to analysis")
                    return ["", ""]
                else:
                    return [self.currentPulseTree.dat_file, index]  # self.seleted_index
            else:
                self.showdialog(
                    "Load a dat file & choose the series that you want to analyze"
                )
            return ["", ""]

    def eventDetection_prepareWindow(self):
        """
        setup the event detection plotitem
        """
        # self.event_traceView.clear()
        def mouseMoved(pos):
            tracePlot = self.event_traceView.getItem(row=0, col=0)
            # if tracePlot.sceneBoundingRect().contains(pos):
            mousePoint = tracePlot.vb.mapSceneToView(pos)
            # index = int(mousePoint.x())
            # t'>x=%0.1f,   <span style='color: red'>y1=%0.1f</span>,   <span style='color: green'>y2=%0.1f</span>" % (mousePoint.x(), data1[index], data2[index]))
            vLine.setPos(mousePoint.x())
            hLine.setPos(mousePoint.y())
            vLine.setZValue(1000)
            hLine.setZValue(1000)

        def mouseClicked(evt):
            tracePlot = self.event_traceView.getItem(row=0, col=0)
            if not self.events.eventLineReloaded:
                self.addEventLinePair(tracePlot, evt)
            if evt == "":
                return
            pos = evt.scenePos()  ## this is for MouseClicked event
            # button = evt.button()
            # pdb.set_trace()
            if tracePlot.sceneBoundingRect().contains(pos):
                mousePoint = tracePlot.vb.mapSceneToView(pos)
                index = int(mousePoint.x() * self.parameters["fs"])
                print(f"x: {mousePoint.x()}, y: {self.events.data[index]}")
                time = self.events.time
                if index >= 0 and index < len(time):
                    """ADDing template selection twin-lines"""
                    if (
                        evt.modifiers() == QtCore.Qt.ControlModifier
                        and evt.button() == QtCore.Qt.LeftButton
                    ):
                        offsetTime = self.event_template_view.p.getValues()["Template"][
                            1
                        ]["Window length"][0]
                        vlineL = infiniteLinePair(
                            parent=tracePlot,
                            offset=offsetTime,
                            span=(0, 0.8),
                            pos=time[index],
                            index=index,
                            Frame=self,
                            angle=90,
                            movable=True,
                            pen="r",
                            hoverPen="w",
                            name="template1",
                        )
                        vlineL.mouseClickEvent(evt)
                        vlineL.twin.sigPositionChangeFinished.connect(
                            vlineL.mouseDragEventTwin
                        )
                        tracePlot.addItem(vlineL, ignoreBounds=True)
                        tracePlot.addItem(vlineL.twin, ignoreBounds=True)
                        self.events.templateTwinLine.append(vlineL)
                        self.updateEvents()

                vLine.setPos(mousePoint.x())
                hLine.setPos(mousePoint.y())

        vLine = pg.InfiniteLine(angle=90, movable=False, name="cursorV")
        hLine = pg.InfiniteLine(angle=0, movable=False, name="cursorH")
        self.eventMarkerLine_v = vLine
        self.eventMarkerLine_h = hLine
        self.event_traceView.clear()
        tracePlot = self.event_traceView.addPlot(row=0, col=0)
        # tracePlot.setMenuEnabled(False)
        tracePlot.addItem(vLine, ignoreBounds=True)
        tracePlot.addItem(hLine, ignoreBounds=True)

        tracePlot.scene().sigMouseMoved.connect(mouseMoved)
        tracePlot.scene().sigMouseClicked.connect(mouseClicked)
        tracePlot.scene().sigMouseClicked.emit("")

        self.events.tracePlot = tracePlot

    def makeNewEventWindow(self, eventObjReset, isConcat=None):
        # self.events = miniEvent()
        if not hasattr(self, "events"):
            self.events = miniEvent()
            self.events.data_raw = []

        else:
            if eventObjReset:
                self.events.init()
        plotHandle = self.event_traceView.getItem(0, 0)
        self.event_matplotView1.clf()
        if plotHandle != None:
            ## clean up
            self.events.tracePlot = []
            try:
                plotHandle.scene().sigMouseClicked.disconnect()
                plotHandle.scene().sigMouseMoved.disconnect()
            except:
                print("No connection")
            # plotHandle.unregister()
            # pdb.set_trace()
            plotHandle.clear()
            self.event_traceView.removeItem(plotHandle)

        self.eventDetection_prepareWindow()
        self.plotConcatenatedSP(isConcat)

    def eventDetectionAction_clicked(self):
        self.makeNewEventWindow(True)
        self.visulization_view.setCurrentIndex(4)

    def event_clearTemplate(self):
        print("to be implemented!")

    def event_template_stateChange(self, param, changes):
        pv = self.event_template_view.p.getValues()
        for param, change, data in changes:
            childName = param.name()
            if childName == "Window length":
                newWindow = pv["Template"][1]["Window length"][0]
                for l in self.events.templateTwinLine:
                    l.updateWindowLength(newWindow)
                self.event_showTemplate(False)
            elif childName in [
                "Show template",
                "Normalizing",
                "Average",
                "Window length",
            ]:
                self.event_showTemplate(toFit=False)
            elif childName == "Fit":
                self.event_showTemplate(toFit=True)
            elif childName == "Clear all templates":
                self.event_clearTemplate()

    def event_pd_stateChange(self, param, changes):
        # pv = self.event_pd_view.p.getValues()
        for param, change, data in changes:
            childName = param.name()
            if childName in ["Threshold (stdev)", "Peak Threshold (stdev)"]:
                if self.events.sweepType == "single":
                    self.event_showDetected()
                    self.event_save2Table()
                else:
                    self.event_caluclateAllSweeps()
                    self.event_save2TableAllSweep()
                self.event_plotHIST()
            elif childName == "Detect current sweep":
                self.event_showDetected()
                self.event_save2Table()
            elif childName == "Detect events for all sweeps":
                self.event_caluclateAllSweeps()
                self.event_save2TableAllSweep()
            elif childName == "fit each event":
                self.event_fitEachEvent()
            elif childName == "show event amplitude statistics":
                self.event_plotHIST()
            elif childName == "show averaged event waveform":
                self.event_averagePlot()
            elif childName == "export event waveforms":
                currentFile = self.currentPulseTree.dat_file
                pv = self.eventParTree_data_view.p.getValues()  ##
                currentTrace = str(pv["Data selection"][1]["Trace"][0])
                fileName = self.exportFile(
                    title="Save event waveforms",
                    defaultName=currentFile[:-4]
                    + "_chan"
                    + currentTrace
                    + "_waveforms.csv",
                    extension="CSV files (*.csv)",
                )
                if fileName != []:
                    # events_waveforms, T = self.event_extractWaveforms()
                    events_waveforms, T = self.event_exportViewedEvents()
                    nSweep = events_waveforms.shape[1]
                    header = ["Event" + str(j) for j in range(nSweep)]

                    df = pd.DataFrame(events_waveforms)
                    df.columns = header
                    df.insert(0, "Time", T)
                    with open(fileName, "w") as f:
                        df.to_csv(f, header=True, line_terminator="\n")

    def exportFile(self, title=None, defaultName=None, extension="all files (*.*)"):
        fileName = QtGui.QFileDialog.getSaveFileName(
            None, title, defaultName, extension
        )
        if isinstance(fileName, tuple):
            return fileName[0]
        else:
            return []

    def eventParTree_stateChange(self, param, changes):
        pv = self.eventParTree_data_view.p.getValues()
        for param, change, data in changes:
            childName = param.name()
            # print('  Parameter: %s    Value:%s'% (childName, str(data)))
        if childName in [
            "Trace",
            "Concatenate all sweeps for current channel",
            "High frequency cutoff",
            "Low frequency cutoff",
        ]:
            self.makeNewEventWindow(True)  ## make new window and reset miniEvent object
            self.event_showTemplate(False)
        elif childName == "Sweep":
            self.events.eventLineReloaded = (
                False  ## turn off flag so we can reload event lines
            )
            self.events.currentSweep = pv["Data selection"][1]["Sweep"][0]
            self.makeNewEventWindow(False)  ## retain miniEvent data
            # self.events.eventLineReloaded = True  ## turn on to supress further reloading

    def showEventRewriteWarning(self):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        ans = msg.question(
            self,
            "",
            "Continue will rewrite events in the table. Click No to cancel?",
            msg.Yes | msg.No,
        )
        msg.setWindowTitle("Event rewrite warning")
        return ans == msg.Yes

    def event_caluclateAllSweeps(self):
        ## It would be too slow to plot all these events. Just show the final table
        ## and plot in a sweep-bases
        # self.events = miniEvent()
        if not hasattr(self, "events"):
            self.events = miniEvent()
        else:
            self.events.init()
        plotHandle = self.event_traceView.getItem(0, 0)
        if plotHandle != None:
            ## clean up
            self.events.tracePlot = []
            try:
                plotHandle.scene().sigMouseClicked.disconnect()
                plotHandle.scene().sigMouseMoved.disconnect()
            except:
                print("No connection")
            # plotHandle.unregister()
            # pdb.set_trace()
            plotHandle.clear()
            self.event_traceView.removeItem(plotHandle)

        self.eventDetection_prepareWindow()
        if self.currentPulseTree.filetype == ".dat":
            self.detectEventForAllSweeps()
        elif self.currentPulseTree.filetype == ".abf":
            self.detectEventForAllSweeps_abf()

    def getEventLine(self):
        p0 = self.events.tracePlot.vb.allChildItems()
        namedObjs = [l for l in p0 if hasattr(l, "name")]
        autoLines = [
            l for l in namedObjs if l.name() == "eventLine" and l.userGen == False
        ]
        manLines = [
            l for l in namedObjs if l.name() == "eventLine" and l.userGen == True
        ]
        return autoLines, manLines

    def event_risingSlope(self, x, y, w=None):

        ransac = linear_model.RANSACRegressor()
        ransac.fit(x, y, sample_weight=w)
        prv = ransac.predict(x)
        slope = ransac.estimator_.coef_[0][0]
        return slope, x, prv

    def event_extractWaveforms(self):
        event_peakIdxs = self.events.peakIndex
        if len(event_peakIdxs) > 0:
            onset_idx = self.events.peakStart
            wl = (
                self.events.upstroke_samples + self.events.downstroke_samples
            )  ## window length
            events_waveforms = np.zeros((len(event_peakIdxs), wl))
            for j, peakIdx in enumerate(event_peakIdxs):
                idx1 = peakIdx - self.events.upstroke_samples
                idx2 = peakIdx + self.events.downstroke_samples
                events_waveforms[j, :] = self.events.data[idx1:idx2]
            T = (
                self.events.time[idx1:idx2] - self.events.time[onset_idx[j]]
            )  ## aligned at onset
            return events_waveforms, T
        else:
            return [], []

    def event_fitEachEvent(self):
        events_waveforms, T = self.event_extractWaveforms()
        if not isinstance(events_waveforms, list):
            fitSlopes = []
            with pg.ProgressDialog(
                "linear fit for event onset",
                0,
                events_waveforms.shape[0] - 1,
                busyCursor=True,
                nested=False,
            ) as dlg:
                for j in range(events_waveforms.shape[0]):
                    dlg += 1
                    slope, *_ = self.event_getLinearSlope(T, events_waveforms[j, :])
                    fitSlopes.append(slope)
            return fitSlopes
        else:
            return []

    def event_getLinearSlope(self, time, data):
        """
        Linear fit of the first phase of event

        Parameters
        ----------
        time : 1-D array
            time of event. zero at event onset
        data : 1-D array
            event waveform

        Returns
        -------
        slope : float
            slope of liinear fit
        lineX : 1-D array
            range of time for the fit
        prv : 1-D array
            fitted value.
        t0_idx : float
            onset time index.
        peak_idx : float
            peak time index.

        """
        peak_idx = np.argmax(np.abs(data - data[0]))
        t0_idx = np.argmin(np.abs(time))
        if t0_idx >= peak_idx:
            return np.nan, np.nan, np.nan, t0_idx, peak_idx, np.nan
        else:
            # if len(data.shape)>1:
            data1 = data.reshape(-1, 1)
            time1 = time.reshape(-1, 1)

            x = (
                time1[t0_idx:peak_idx] * 1000
            )  ## only the intial section: between onset and peak
            y = data1[t0_idx:peak_idx] * 1000
            ##  centered at steepest slope
            dy = np.diff(data)[t0_idx:peak_idx] * 1000  ## changes
            maxDvIdx = np.argmax(np.abs(dy))  ## steep slop time index
            steepIdx = maxDvIdx + t0_idx
            lb0 = np.max([3, maxDvIdx - 15])  ## lower bound
            lb_range = list(np.arange(lb0, maxDvIdx))
            ub0 = np.min([maxDvIdx + 15, len(x) - 3])  ## upper bound
            ub_range = list(np.arange(maxDvIdx, ub0))
            whole_range = lb_range + ub_range

            x = np.array([x[j] for j in whole_range])
            y = np.array([y[j] for j in whole_range])
            if len(x) > 1:
                slope, lineX, prv = self.event_risingSlope(x, y)
                return slope, lineX, prv, t0_idx, peak_idx, steepIdx
            else:
                return np.nan, np.nan, np.nan, t0_idx, peak_idx, np.nan

    def event_averagePlot(self):
        """average event waveform with linear fit of inital slope"""
        events_waveforms, T = self.event_exportViewedEvents()
        events_waveforms = np.transpose(events_waveforms)
        data_ = self.events.data_raw.copy()
        data_[np.abs(data_) > np.mean(data_) + 3 * np.std(data_)] = np.mean(data_)
        raw_baseline = np.mean(data_) * 1000.0
        if not isinstance(events_waveforms, list):
            mean_event = np.mean(events_waveforms, axis=0)
            (
                slope,
                lineX,
                prv,
                extremIdxTime,
                extremIdxData,
                steepIdx,
            ) = self.event_getLinearSlope(T * 1000, mean_event)
            v = list(
                events_waveforms.reshape(
                    (np.product(events_waveforms.shape),), order="C"
                )
                * 1000.0
            )
            t = list(T * 1000) * events_waveforms.shape[0]
            df = pd.DataFrame(list(zip(v, t)), columns=["Voltage", "Time"])
            self.eventData_tabview.setCurrentIndex(3)

            from matplotlib.collections import PolyCollection

            mplWidget = self.event_matplotView1
            mplWidget.clf()
            gs = mplWidget.figure.add_gridspec(1, 1)
            axes = mplWidget.figure.add_subplot(gs[0, 0])
            ax = sns.lineplot(data=df, x="Time", y="Voltage", ax=axes, ci="sd")
            for child in ax.findobj(PolyCollection):
                child.set_color([1, 0, 0, 0.7])
            ax.set(xlabel="Time (ms)", ylabel="Voltage (mV)")
            ## slope
            axes.plot(lineX, prv, "--", color="r", linewidth=1.5)
            ## baseline
            axes.axhline(y=raw_baseline, color="gray", linestyle="--")

            datName = self.events.datFile[:-4].split("\\")[-1]
            axes.set_title(f"{datName} wavform: averaged \u00B1 standard deviation")
            axes.set_xlim([T[0] * 1000 - 3, T[-1] * 1000 + 3])
            mplWidget.draw()

    def event_plotHIST(self):
        df = self.events.eventTableDF
        mplWidget = self.event_matplotView1
        mplWidget.clf()
        self.eventData_tabview.setCurrentIndex(3)
        datName = self.events.datFile[:-4].split("\\")[-1]
        gs = mplWidget.figure.add_gridspec(2, 2)
        axes = mplWidget.figure.add_subplot(gs[0, 0])
        if self.events.traceYUnit == "mV":
            scalingFactor = 1000
        else:
            scalingFactor = 1
        abs_peakhist, abs_binVal, _ = axes.hist(
            df["peak abs-height"] * scalingFactor, bins=25
        )
        axes.spines["right"].set_visible(False)
        axes.spines["top"].set_visible(False)
        axes.set_xlabel(f"Peak height ({self.events.traceYUnit})")
        axes.set_ylabel("# Peaks")
        axes.set_title(f"{datName}: absolute peak height")

        axes = mplWidget.figure.add_subplot(gs[0, 1])
        abs_peakcum = np.cumsum(abs_peakhist) / np.sum(abs_peakhist) * 100
        axes.plot(abs_binVal[1:], abs_peakcum)
        axes.spines["right"].set_visible(False)
        axes.spines["top"].set_visible(False)
        axes.set_xlabel(f"Peak height ({self.events.traceYUnit})")
        axes.set_ylabel("Cumulative percentage (%)")
        axes.set_title(f"{datName}: culmulative absolute peak height")

        axes = mplWidget.figure.add_subplot(gs[1, 0])
        loc_peakhist, loc_binVal, _ = axes.hist(
            df["peak loc-height"] * scalingFactor, bins=25
        )
        axes.spines["right"].set_visible(False)
        axes.spines["top"].set_visible(False)
        axes.set_xlabel(f"Peak height ({self.events.traceYUnit})")
        axes.set_ylabel("# Peaks")
        axes.set_title(f"{datName}: local peak height")

        axes = mplWidget.figure.add_subplot(gs[1, 1])
        loc_peakcum = np.cumsum(loc_peakhist) / np.sum(loc_peakhist) * 100
        axes.plot(loc_binVal[1:], loc_peakcum)
        axes.spines["right"].set_visible(False)
        axes.spines["top"].set_visible(False)
        axes.set_xlabel(f"Peak height ({self.events.traceYUnit})")
        axes.set_ylabel("Cumulative percentage (%)")
        axes.set_title(f"{datName}: culmulative local peak height")

        mplWidget.figure.tight_layout()
        mplWidget.draw()

    def event_save2TableAllSweep(self):
        ## grab events from marked event

        ss = self.events.datFile
        if self.events.eventMarker == None:
            print("No events are marked!")
            return
        eventMarkers = self.events.eventMarker
        eventPenRvalue = [
            eventMarkers.getSpotOpts(yy)[2].color().getRgb()[0]
            for yy in eventMarkers.data
        ]  # red is for disabled one
        autoType = [
            True if eventPenRvalue[idx] == 0 else False
            for idx, index in enumerate(self.events.peakIndex)
        ]
        # autoIdx = [index for idx, index in enumerate(self.events.peakIndex) if eventStatus[idx]==0] ## only grab green ones
        autoIdx = self.events.peakIndex  ## only grab green ones

        autoTime = [self.events.time[idx] for idx in autoIdx]
        peak_onset = (
            self.events.peakStart
        )  ## peak starting point; corresponding to peak of deconvoled peak

        # pdb.set_trace()
        peak_onset = (
            self.events.peakStart
        )  ## peak starting point; corresponding to peak of deconvoled peak
        peak_downStroke50 = [idx + self.events.downstroke_samples for idx in autoIdx]
        ## abosolute height
        peak_heights = [self.events.data[p] for p in autoIdx]

        baseline_ = np.mean(self.events.data)
        median_peakHeight = np.median(peak_heights)
        data_ = self.events.data.copy()
        rcduration = [self.events.time[-1]] * len(data_)
        data_[data_ > median_peakHeight] = baseline_
        baseline_ = np.mean(data_)
        peak_heights_abs = [
            peak - baseline_ for peak in peak_heights
        ]  ## peak hight relative to baseline
        ## local height
        peak_heights_local = [
            peak - self.events.data[peak_onset[idx]]
            for idx, peak in enumerate(peak_heights)
        ]  ## peak hight relative to baseline
        pv = self.event_pd_view.p.getValues()
        isFitSlope = pv["Peak parameters for raw trace"][1]["Linear fit during onset"][
            0
        ]
        if isFitSlope:
            slopFit = self.event_fitEachEvent()  ## this fit the inital segemnt of event
        else:
            slopFit = [np.nan] * len(peak_heights_abs)

        df1 = pd.DataFrame(
            list(
                zip(
                    peak_onset,
                    autoIdx,
                    slopFit,
                    autoTime,
                    peak_heights_local,
                    peak_heights_abs,
                    rcduration,
                    autoType,
                )
            ),
            columns=[
                "onset",
                "peak",
                "slope",
                "peak time",
                "peak loc-height",
                "peak abs-height",
                "recordingTime",
                "InUse",
            ],
        )
        df1.insert(loc=0, column="type", value="auto")
        df1.insert(loc=0, column="trace", value=self.events.traceName)
        df1.insert(loc=0, column="series", value=self.events.seriesName)
        df1.insert(loc=0, column="dat", value=ss.split("\\")[-1])

        df = df1
        df.sort_values(
            by="peak time", axis=0, ascending=True, inplace=True, ignore_index=False
        )
        self.events.eventTableDF = df  ## this is to cmmunicate the data back to trace!

        dft = np.array(
            df.to_records(index=False)
        )  ## format dataframe for using in QtTable wiget
        self.eventTable.clear()
        self.eventTable.appendData(dft.copy())  ## update QtTable widget
        self.eventTable.show()
        self.eventData_tabview.setCurrentIndex(2)

    def event_save2Table(self):
        ## grab events from marked event
        # eventsIdx = self.retrieveEventIdx()
        l_auto, l_man = self.getEventLine()
        ss = self.events.datFile
        autoTime = [x.value() for x in l_auto]
        autoType = [x.active for x in l_auto]
        autoIdx = [x.timeIdx for x in l_auto]  ## peak
        peak_onset = (
            self.events.peakStart
        )  ## peak starting point; corresponding to peak of deconvoled peak
        peak_downStroke50 = [idx + self.events.downstroke_samples for idx in autoIdx]
        ## abosolute height
        peak_heights = [self.events.data[p] for p in autoIdx]
        baseline_ = np.mean(self.events.data)

        median_peakHeight = np.median(peak_heights)
        data_ = self.events.data.copy()
        rcduration = [self.events.time[-1]] * len(data_)
        data_[data_ > median_peakHeight] = baseline_
        baseline_ = np.mean(data_)
        peak_heights_abs = [
            peak - baseline_ for peak in peak_heights
        ]  ## peak hight relative to baseline
        ## local height
        peak_heights_local = [
            peak - self.events.data[peak_onset[idx]]
            for idx, peak in enumerate(peak_heights)
        ]  ## peak hight relative to baseline
        # peak_idx, peak_heights = fitFuncs.getRawPeaks(self.events.data, autoTime)
        # time = self.events.time[peak_idx]
        df1 = pd.DataFrame(
            list(
                zip(
                    peak_onset,
                    autoIdx,
                    peak_downStroke50,
                    autoTime,
                    peak_heights_local,
                    peak_heights_abs,
                    rcduration,
                    autoType,
                )
            ),
            columns=[
                "onset",
                "peak",
                "downStroke",
                "peak time",
                "peak loc-height",
                "peak abs-height",
                "recordingTime",
                "InUse",
            ],
        )
        df1.insert(loc=0, column="type", value="auto")
        df1.insert(loc=0, column="trace", value=self.events.traceName)
        df1.insert(loc=0, column="series", value=self.events.seriesName)
        df1.insert(loc=0, column="dat", value=ss.split("\\")[-1])

        if l_man != []:
            manTimeIdx = [x.timeIdx for x in l_man]
            manType = [x.active for x in l_auto]
            time2 = self.events.time[manTimeIdx]
            peak2 = self.events.data[manTimeIdx] - baseline_
            peak2_local = [
                peak - self.events.data[manTimeIdx[idx - 50]]
                for idx, peak in enumerate(peak2)
            ]  ## peak hight relative to baseline
            peak_onset = [idx - self.event.upstroke_samples for idx in manTimeIdx]
            peak_downStroke50 = [
                idx - self.event.downstroke_samples for idx in manTimeIdx
            ]
            df2 = pd.DataFrame(
                list(
                    zip(
                        peak_onset,
                        manTimeIdx,
                        peak_downStroke50,
                        time2,
                        peak2_local,
                        peak2,
                        rcduration,
                        manType,
                    )
                ),
                columns=[
                    "onset",
                    "peak",
                    "downStroke",
                    "peak time",
                    "peak loc-height",
                    "peak abs-height",
                    "recordingTime",
                    "InUse",
                ],
            )
            df2.drop_duplicates(keep="first", inplace=True, ignore_index=False)
            df2.insert(loc=0, column="type", value="manual")
            df2.insert(loc=0, column="trace", value=self.events.traceName)
            df2.insert(loc=0, column="series", value=self.events.seriesName)
            df2.insert(loc=0, column="dat", value=ss.split("\\")[-1])

            df = df2.append(df1, ignore_index=False)
        else:
            df = df1
        df.sort_values(
            by="peak time", axis=0, ascending=True, inplace=True, ignore_index=False
        )
        self.events.eventTableDF = df  ## this is to cmmunicate the data back to trace!

        dft = np.array(
            df.to_records(index=False)
        )  ## format dataframe for using in QtTable wiget
        # self.eventTable.clear()
        self.eventTable.appendData(dft.copy())  ## update QtTable widget
        self.eventTable.show()
        self.eventData_tabview.setCurrentIndex(2)

    def event_saveToDisk(self):
        pass

    def setTresholdLine(self, t):
        self.setMarkerLineValue("threshold", t)

    def mouseDragEventthreholdLine(self, evt):
        tracePlot = self.events.decovPlot
        pos = evt.scenePos()  ## this is for MouseClicked event
        mousePoint = tracePlot.vb.mapSceneToView(pos)
        # button = evt.button()
        # pdb.set_trace()
        if tracePlot.sceneBoundingRect().contains(pos):
            mousePoint = tracePlot.vb.mapSceneToView(pos)
            # t = self.events.threholdLine.pos().y()
            t = mousePoint.y()
            m, std0 = self.events.m, self.events.std
            nSTD = (t - m) / std0
            # self.eventParTree_data_view.blockTreeChangeEmit = 1
            # self.eventParTree_data_view.p.blockTreeChangeEmit = 1
            self.events.threholdPar.setValue(nSTD)

    def event_getPeaks(self, D, nSTD, miniDist):
        m, th, std0 = self.calBaselineValue(D, nSTD)
        peak_idx, peak_heights, contour_heights = fitFuncs.getPeaks(
            D, height=th, width=(10, 80), wlen=150, distance=miniDist
        )
        return peak_idx, peak_heights, contour_heights

    def event_deconvHist(self, D):
        plotHandle = self.event_traceView2.getItem(0, 0)
        plotHandle.clear()
        self.event_traceView2.removeItem(plotHandle)
        plotHandle = self.event_traceView2.addPlot(0, 0)

        threholdLine = pg.InfiniteLine(
            angle=90,
            pen=pg.mkPen(color="g"),
            hoverPen="r",
            movable=True,
            name="threshold",
        )
        threholdLine.sigPositionChangeFinished.connect(self.mouseDragEventthreholdLine)
        eventHist, bin_edges = np.histogram(D, bins=250)
        wd = bin_edges[2] - bin_edges[1]
        bargraph = pg.BarGraphItem(x=bin_edges[:-1], height=eventHist, width=wd)
        plotHandle.addItem(bargraph)
        plotHandle.addItem(threholdLine)
        self.events.threholdLine = threholdLine
        self.events.threholdPar = self.eventParTree_data_view.p.param(
            "Detection"
        ).children()[0]

    def event_showDetected(self, rescale=True):
        if self.events.time == []:
            return
        if self.templateFit.items == []:
            self.showdialog("Fit a template first!")
            return
        # TT0 = sysTime.clock()
        time = self.events.time.copy()
        D, upstroke_samples, downstroke_samples, fitPars = self.event_detection()
        p0 = self.events.tracePlot
        p0.setLabel("bottom", "", units="")
        p0.showAxis("bottom", False)
        if self.events.eventMarker != None:
            self.events.eventMarker.clear()
            p0.removeItem(self.events.eventMarker)
            self.events.eventMarker = None

        if self.events.decovPlot != "":
            plotHandle = self.events.decovPlot
            # self.events.threholdLine.disconnect()
            plotHandle.clear()
            try:
                plotHandle.removeItem(self.events.threholdLine)
            except:
                print("line 2737: nothing to be removed")
            self.events.decovPlot = ""
            try:
                self.event_traceView.removeItem(plotHandle)
            except:
                print("line 2742: Nothing to be removed")
            try:
                l1, l2 = self.getEventLine()
                if l2 != []:
                    l = l1.extend(l2)
                else:
                    l = l1
                for x in l:
                    p0.removeItem(x)
            except:
                print("line 2752: Nothing to be removed")
        plotHandle = self.event_traceView.addPlot(row=1, col=0)
        plotHandle.showGrid(x=True, y=True)

        # plotHandle = self.event_traceView.addPlot(row=1, col=0)
        def mouseMoved(pos):
            tracePlot = self.events.decovPlot
            mousePoint = tracePlot.vb.mapSceneToView(pos)
            vLine.setPos(mousePoint.x())

        def mouseClicked(evt):
            tracePlot = self.events.tracePlot
            pos = evt.scenePos()  ## this is for MouseClicked event
            # button = evt.button()
            # pdb.set_trace()
            if tracePlot.sceneBoundingRect().contains(pos):
                mousePoint = tracePlot.vb.mapSceneToView(pos)
                index = int(mousePoint.x() * self.parameters["fs"])
                # print(mousePoint.x(), self.events.data[index])
                # print('X index', index)
                time = self.events.time
                # pdb.set_trace()
                if index >= 0 and index < len(time):
                    """ADDing template selection twin-lines"""
                    if (
                        evt.modifiers() == QtCore.Qt.ShiftModifier
                        and evt.button() == QtCore.Qt.LeftButton
                    ):
                        vpeakline = EventMark(
                            mousePoint.x(), index, penc="r", userGen=True
                        )
                        p0.addItem(vpeakline)

        vLine = pg.InfiniteLine(angle=90, movable=False, name="cursorH_D")
        # baselineLine = pg.InfiniteLine(angle=0, pen = pg.mkPen(color = 'w', style=QtCore.Qt.DashLine),
        #                                hoverPen = 'm', movable=False, name='baseline')
        pv = self.event_pd_view.p.getValues()
        nSTD = pv["Peak parameters for deconvoled trace"][1]["Threshold (stdev)"][0]
        m, th, std0 = self.calBaselineValue(D, nSTD)

        g = plotHandle.plot(time, D, pen=pg.mkPen("w"), name="Deconvolution")

        ## parameters for detecting peaks from the deconvoluted trace
        miniDist = int(pv["Peak parameters for deconvoled trace"][1]["Distance"][0])
        wlen = int(pv["Peak parameters for deconvoled trace"][1]["Wlen (samples)"][0])
        if self.events.traceYUnit == "A":
            prominence = int(
                pv["Peak parameters for deconvoled trace"][1]["Prominence (pA)"][0]
            )
        else:
            prominence = int(
                pv["Peak parameters for deconvoled trace"][1]["Prominence (mV)"][0]
            )
        width_min = int(
            pv["Peak parameters for deconvoled trace"][1]["Width"][1]["minimal"][0]
        )
        width_max = int(
            pv["Peak parameters for deconvoled trace"][1]["Width"][1]["maximal"][0]
        )
        # pdb.set_trace()
        peak_onset_idx, peak_heights_D, peak_props = fitFuncs.getPeaks(
            D,
            height=th,
            width=(width_min, width_max),
            wlen=wlen,
            distance=miniDist,
            prominence=prominence,
        )

        # peakOnset_raw = self.events.data[peak_idx0]
        # onset_nSD = pv['Peak parameters for raw trace'][1]['Onset Threshold (stdev)'][0]
        data_ = self.events.data.copy()
        data_[np.abs(data_) > np.mean(data_) + 3 * np.std(data_)] = np.mean(data_)
        raw_baseline = np.mean(data_)
        x_th = [self.events.time[0], self.events.time[-1]]
        deltD = np.median(np.diff(self.events.data))
        peak_nSD = pv["Peak parameters for raw trace"][1]["Peak Threshold (stdev)"][0]
        peak_threh = raw_baseline + peak_nSD * np.sign(fitPars[2]) * np.std(data_)
        try:
            p0.removeItem(self.events.threholdLine_raw)  ## removing the old line
        except:
            print("")
        threholdLine_raw = p0.plot(
            x_th,
            [peak_threh, peak_threh + deltD],
            pen=pg.mkPen("r", style=QtCore.Qt.DashLine),
        )
        p0.plot(
            x_th,
            [raw_baseline, raw_baseline + deltD],
            pen=pg.mkPen("w", style=QtCore.Qt.DashLine),
        )
        peak_idx, selected_index = fitFuncs.getRawPeaks3(
            self.events.data, D, peak_onset_idx, peak_threh, np.sign(fitPars[2])
        )
        self.events.globalBaseline = raw_baseline
        self.events.peakIndex = peak_idx  ## the actual peak index
        self.events.peakStart = [peak_onset_idx[j] for j in selected_index]
        self.events.downstroke_samples = downstroke_samples
        self.events.upstroke_samples = upstroke_samples
        # peak_idx = [int(peakIdx_) for peakIdx_ in peak_props['right_ips']]

        xunit = self.events.traceXUnit
        yunit = self.events.traceYUnit
        plotHandle.setLabels(bottom=("Time", xunit), left=(yunit))

        plotHandle.setXLink(p0)
        p0.setZValue(0)
        g.setZValue(0)
        plotHandle.addItem(vLine, ignoreBounds=True)
        for p in peak_idx:  # contour_heights[j], peak_heights[j]
            vpeakline = EventMark(time[p], p)
            p0.addItem(vpeakline)
        p0.autoRange()
        self.events.sweepType = "single"  ## or 'all'

        y_th = [th, th + deltD]
        threholdLine = plotHandle.plot(x_th, y_th, pen=pg.mkPen("r"))
        threholdLine.setZValue(1000)
        self.events.threholdLine = threholdLine
        self.events.threholdLine_raw = threholdLine_raw
        self.events.threholdPar = self.event_pd_view.p.param(
            "Peak parameters for deconvoled trace"
        ).children()[0]
        plotHandle.scene().sigMouseMoved.connect(mouseMoved)
        plotHandle.scene().sigMouseClicked.connect(mouseClicked)
        plotHandle.setYRange(D.min(), D.max())
        self.events.decovPlot = plotHandle

    def event_detection(self, useData=None, useTime=None):

        # get selected template
        pv = self.event_template_view.p.getValues()
        funcName = pv["Template"][1]["Fit options"][1]["function"][0]
        templateId = pv["Template"][1]["Template to use"][0]
        # get template data from the template table
        data0 = self.templateFit.serialize(useSelection=False)
        data = data0.splitlines()[templateId].split("\t")
        popts = data[3:-1]
        popts = [float(f) for f in popts]
        func = fitFuncs.curveFitFuncs.get(funcName)

        fs = self.parameters["fs"]
        if isinstance(
            useData, type(None)
        ):  ## use sweep data. Need to reorganize this part!
            time = self.events.time.copy()
            template = func(time, *popts)
            deconv_data = fitFuncs.signal_deconvolution(self.events.data, template, fs)
        else:
            time = useTime
            template = func(time, *popts)
            deconv_data = fitFuncs.signal_deconvolution(useData, template, fs)

        self.events.D = deconv_data
        upstroke_samples = int(popts[0] * self.parameters["fs"])
        downstroke_samples = int(popts[1] * self.parameters["fs"])
        return deconv_data, upstroke_samples, downstroke_samples, popts

    def normEvent(self, data, imax, peakSign):
        data0 = data[:imax]
        if peakSign > 0:
            lowQuntile = np.quantile(data0, 0.1)
            baseline = np.mean(data0[data0 <= lowQuntile])
        else:
            highQuntile = np.quantile(data0, 0.9)
            baseline = np.mean(data0[data0 >= highQuntile])
        data0 = data - baseline
        peak = np.abs(data0[imax])
        data0 = data0 / peak
        return data0

    def alignWithPeak(self, time, data, isNorm):
        baseline = np.mean(data[:20])
        d = data - baseline
        imax = np.argmax(np.abs(d))
        peakSign = np.sign(d[imax])
        if isNorm:
            data = self.normEvent(data, imax, peakSign)
        time = time - time[imax]
        return time, data, peakSign

    def getEventData(self, l):
        startPos = l.pos().x()
        endPos = l.twin.pos().x()
        index1 = int(startPos * self.events.seriesFS)
        index2 = int(endPos * self.events.seriesFS)
        time_ = self.events.time[index1:index2]
        data_ = self.events.data[index1:index2]
        return time_, data_

    def updateTemplateTable(self, dt, popt):
        ## convert dataframe into a record that can be populated into a QtTable
        dtKeys = list(dt.keys())
        # pdb.set_trace()
        for idx, k in enumerate(dtKeys):
            dt[k] = popt[idx]
        dtKeys.insert(0, "trace")
        dtKeys.insert(0, "series")
        dtKeys.insert(0, "session")
        ss = self.events.datFile
        dt["trace"] = self.events.traceName
        dt["session"] = ss.split("\\")[-1]
        dt["series"] = self.events.seriesName
        funcPars = {}
        for k in dtKeys:
            funcPars[k] = dt[k]
        df = pd.DataFrame.from_dict(funcPars, orient="index").transpose()
        df = np.array(
            df.to_records(index=False)
        )  ## format dataframe for using in QtTable wiget
        self.templateFit.clear()  ## clear table for new template
        self.templateFit.appendData(df.copy())  ## update QtTable widget
        self.templateFit.show()

    def event_showTemplate(self, toFit):
        plotHandle = self.event_traceView2.getItem(0, 0)
        plotHandle.clear()
        if self.events.eventTime:
            time, data = self.events.get_events()
            if len(time) == 0:
                return  ## no events

            index = []
            t = []
            d = []
            isNorm = self.event_template_view.p.getValues()["Template"][1][
                "Normalizing"
            ][0]
            if not toFit:
                isAvg = self.event_template_view.p.getValues()["Template"][1][
                    "Average"
                ][0]
            else:
                isAvg = True  ## get avg is we want to fit a template

            if isAvg:
                isNorm = True

            nEvent = len(time)
            for idx in range(nEvent):
                time_, data_, peakSign = self.alignWithPeak(
                    time[idx], data[idx], isNorm
                )
                index.extend([idx] * len(time_))
                t.extend(time_)
                d.extend(data_)
                if not toFit:
                    plotHandle.plot(time_, data_, pen=(idx, nEvent))
            df = {"index": index, "time": t, "data": d}
            df = pd.DataFrame(df)
            minT = df[["index", "time"]].groupby("index").min().max().values[0]
            # minT = -0.01
            maxT = df[["index", "time"]].groupby("index").max().min().values[0]

            df = df.loc[
                (df["time"] > minT) & (df["time"] < maxT)
            ]  ## limit the overlap range
            nbins = int((maxT - minT) * self.parameters["fs"] / 4)
            # print(f'Number of bins: {nbins}')
            qCut = pd.cut(df["time"], bins=np.linspace(minT, maxT, nbins))  ## binning
            groups = df.groupby(qCut)
            if isAvg:
                tx, ty = groups.mean().time.values, groups.mean().data.values
                if peakSign > 0:
                    offset1 = (65, 25)
                else:
                    offset1 = (-40, -20)
                if toFit:
                    pv = self.event_template_view.p.getValues()
                    funcName = pv["Template"][1]["Fit options"][1]["function"][0]
                    riseTime = pv["Template"][1]["Fit options"][1][
                        "fitting parameters"
                    ][1]["Tau fast"][0]
                    decayTime = pv["Template"][1]["Fit options"][1][
                        "fitting parameters"
                    ][1]["Tau slow"][0]
                    bounds = (
                        [0.00001, 0.001, -10.0, -10.0, -10.0],
                        [0.01, 0.03, 10, 15.0, 10.0],
                    )
                    p0 = np.array([riseTime / 1000, decayTime / 1000, 1.0, 1.0, 0.01])
                    func = fitFuncs.curveFitFuncs.get(funcName)
                    funcPars0 = fitFuncs.curveFitPars.get(funcName).copy()
                    # #fitFunc =
                    if peakSign > 0:
                        j = np.where(ty >= 0.15)[0][0]
                    else:
                        j = np.where(ty < -0.12)[0][0]
                    tx2 = tx - tx[j]
                    plotHandle.addLegend(offset=offset1)
                    g = plotHandle.plot(
                        tx2[j:], ty[j:], pen=pg.mkPen("w", width=5), name="Average"
                    )  ## average rsponses

                    # pdb.set_trace()
                    try:
                        popt, pcov = self.fitTemplate(func, tx2[j:], ty[j:], p0, bounds)
                    except Exception as e:
                        print(e)
                    else:
                        print(f"fitting parameters:{popt}")
                        self.updateTemplateTable(funcPars0, popt)
                        template = func(tx2[j:], *popt)
                        self.events.template.append(template)
                        g = plotHandle.plot(
                            tx2[j:],
                            template,
                            pen=pg.mkPen("r", width=2),
                            name="Fit: Tau fast = %3.2f ms, Tau slow = %3.2f ms"
                            % tuple(popt[:2] * 1000),
                        )
                        g.setAlpha(0.6, False)
                        self.eventData_tabview.setCurrentIndex(1)

                else:
                    plotHandle.addLegend(offset=offset1)
                    g = plotHandle.plot(tx, ty, pen=pg.mkPen("w", width=5), name="Avg.")
                    g.setAlpha(0.6, False)

                plotHandle.setLabels(bottom=("Time", "S"), left=(""))
            plotHandle.autoRange()

    def fitTemplate(self, func, tx, ty, p0, bounds):
        popt, pcov = curve_fit(func, tx, ty, p0, bounds=bounds)
        return popt, pcov

    def removeStimulationArtifacts(self, data, th):
        dvdt = (
            np.diff(np.abs(data)) * self.parameters["fs"]
        )  ## dv/dt. ABS operations make sure "spiking" are positive
        peaks_idx, peak_heights_D, _ = fitFuncs.getPeaks(dvdt, height=th)
        return peaks_idx

    def detectEventForAllSweeps_abf(self):
        datName, series_index = self.checkNodeIndex()
        if series_index == "":
            print("Nothing to process!")
            return
        plotHandle = self.events.tracePlot  # self.event_traceView.getItem(0, 0)
        if plotHandle == []:
            return
        pv = self.eventParTree_data_view.p.getValues()  ##
        lfcut = pv["Data selection"][1]["Low frequency cutoff"][0]
        hfcut = pv["Data selection"][1]["High frequency cutoff"][0]

        self.events.isConcat = True
        self.events.traceYUnit = self.currentPulseTree.abf.yUnits
        if self.events.traceYUnit == "V" or self.events.traceYUnit == "mV":
            outlierCutoff_LV = pv["PSP Outliers"][1][
                "Outlier voltage (mV) - lower bound"
            ][0]
            outlierCutoff_UV = pv["PSP Outliers"][1][
                "Outlier voltage (mV) - upper bound"
            ][0]
            outlierCutoff_rv = pv["PSP Outliers"][1]["replacement value"][0]
        else:
            outlierCutoff_LV = pv["PSC Outliers"][1][
                "Outlier voltage (pA) - lower bound"
            ][0]
            outlierCutoff_UV = pv["PSC Outliers"][1][
                "Outlier voltage (pA) - upper bound"
            ][0]
            outlierCutoff_rv = pv["PSC Outliers"][1]["replacement value"][0]

        ## read data
        time, data = self.extractSingleSeries_ABF(series_index)
        nSweep = data.shape[1]
        data = data.reshape(
            (np.product(data.shape),), order="F"
        )  ## flatten from 2D to 1D concatenated
        self.events.data_raw = data  ## save a copy of raw data
        data = self.bandPass_signal(
            data, highCutOff=hfcut, lowCutOff=lfcut, useButter=True
        )
        time = np.arange(len(data)) / self.parameters["fs"]

        ## remove spiking and/or stimuli artifacts
        data_ = data.copy()
        baseline_ = np.mean(data_)
        std_ = np.std(data_)
        data_[data > baseline_ + 3 * std_] = baseline_
        data_[data < baseline_ - 3 * std_] = baseline_
        baseline_ = np.mean(data_)  ## iterative estamtion of baseline
        if pv["Spikes"][1]["Removing spikes"][0]:
            dvdt_th = pv["Spikes"][1]["dv/dt (V/s) - threhold"][0]
            peaks_lv = self.removeStimulationArtifacts(data.copy(), dvdt_th)
            for p in peaks_lv:
                data[p - 20 : p + 20] = baseline_

        currentTrace = pv["Data selection"][1]["Trace"][0] - 1
        self.events.node = []
        self.events.nSweep = nSweep
        self.events.seriesName = series_index[1] + 1
        self.events.traceName = currentTrace + 1
        self.events.datFile = datName  # self.currentPulseTree.dat_file
        self.events.setSeriesSampleRate(self.parameters["fs"])  ## series sampling rate
        self.events.globalBaseline = baseline_
        self.events.time = time
        self.events.data = data
        mean_data = np.mean(data)
        median_data = np.median(data)

        if outlierCutoff_rv == "bound":
            data[data <= outlierCutoff_LV] = outlierCutoff_LV
            data[data >= outlierCutoff_UV] = outlierCutoff_UV
        elif outlierCutoff_rv == "median":
            data[data <= outlierCutoff_LV] = median_data
            data[data >= outlierCutoff_UV] = median_data
        elif outlierCutoff_rv == "mean":
            data[data <= outlierCutoff_LV] = mean_data
            data[data >= outlierCutoff_UV] = mean_data
        currentTrace = pv["Data selection"][1]["Trace"][0] - 1
        myPen = pg.mkPen(
            color=pg.intColor(currentTrace, hues=8)
        )  ## the pen to draw this
        self.events.tracePlot.plot(time, data, pen=myPen)
        if self.events.traceYUnit == "mV":
            self.events.tracePlot.setLabels(bottom=("Time", "S"), left=("V"))
        self.verifyAllSweepEvents(data, time)

    def detectEventForAllSweeps(self):
        datName, series_index = self.checkNodeIndex()
        if series_index == "":
            print("Nothing to process!")
            return
        (
            bundleClass,
            stimChanLabels,
            serieIndex,
            treeChildrenIdx,
        ) = self.getDatFilesFromTreeWiget("Spontaneous", datName, series_index)
        plotHandle = self.events.tracePlot  # self.event_traceView.getItem(0, 0)
        if plotHandle == []:
            return
        if bundleClass == []:
            print("data not found! line 3241")
            # pdb.set_trace()
            bundleClass = HekaBundleInfo(
                self.currentPulseTree.bundle.file_name
            )  ## a bundle with some extra functions
            nSweep = bundleClass.countSweeps(series_index)
            self.parameters["fs"] = bundleClass.getSeriesSamplingRate(series_index)  ##
            self.events.setSeriesSampleRate(
                self.parameters["fs"]
            )  ## series sampling rate
            serieIndex = series_index
            # pdb.set_trace()
            self.events.node = self.currentPulseTree
        else:
            self.events.node = self.selTreesSP.selectedItems()[0]
        nSweep = bundleClass.countSweeps(series_index)

        pv = self.eventParTree_data_view.p.getValues()  ##
        currentSweep = pv["Data selection"][1]["Sweep"][0]
        currentTrace = pv["Data selection"][1]["Trace"][0] - 1
        lfcut = pv["Data selection"][1]["Low frequency cutoff"][0]
        highCutOff = pv["Data selection"][1]["High frequency cutoff"][0]

        self.events.nSweep = nSweep
        self.events.datFile = self.currentPulseTree.dat_file
        self.events.seriesName = series_index[1] + 1
        self.events.traceName = currentTrace + 1
        self.events.setSeriesSampleRate(
            bundleClass.getSeriesSamplingRate(serieIndex)
        )  ## series sampling rate
        self.parameters["fs"] = bundleClass.getSeriesSamplingRate(
            serieIndex
        )  ## this need to be change!
        # self.currentPulseTree= self.selTreesSP.selectedItems()[0]
        self.currentPulseTree.bundle = bundleClass.bundle

        seriesIdx = list(series_index.copy())  ## get series level index
        seriesIdx.append(0)  ## sweep 0
        seriesIdx.append(currentTrace)
        ## loops through all the sweeps
        SweepIdx = range(nSweep)
        self.events.isConcat = True
        data_ = []

        ## concatenating
        with pg.ProgressDialog(
            "Concatenating traces", 0, nSweep - 1, busyCursor=True, nested=False
        ) as dlg:
            for sweep in SweepIdx:
                dlg += 1
                seriesIdx[2] = sweep  ## change sweeep level index
                data = self.currentPulseTree.bundle.data[seriesIdx]
                data_.extend(data)
            if dlg.wasCanceled():
                raise Exception("Processing canceled by user")

        ## preprocessing
        lfcut = pv["Data selection"][1]["Low frequency cutoff"][0]
        hfcut = pv["Data selection"][1]["High frequency cutoff"][0]

        self.events.data_raw = np.array(data_)
        data = self.bandPass_signal(
            data_, highCutOff=hfcut, lowCutOff=lfcut, useButter=True
        )
        mean_data = np.mean(data)
        if pv["Spikes"][1]["Removing spikes"][0]:
            dvdt_th = pv["Spikes"][1]["dv/dt (V/s) - threhold"][0]
            peaks_lv = self.removeStimulationArtifacts(data.copy(), dvdt_th)
            for p in peaks_lv:
                data[p - 20 : p + 20] = mean_data
        time = np.arange(len(data)) / self.parameters["fs"]
        self.events.time = time
        self.events.data = data
        mean_data = np.mean(data)
        median_data = np.median(data)
        trace = bundleClass.bundle.pul[seriesIdx[0]][seriesIdx[1]][seriesIdx[2]][
            seriesIdx[3]
        ]  ## get trace meta information
        self.events.traceYUnit = trace.YUnit
        if trace.YUnit == "V":
            outlierCutoff_LV = pv["PSP Outliers"][1][
                "Outlier voltage (mV) - lower bound"
            ][0]
            outlierCutoff_UV = pv["PSP Outliers"][1][
                "Outlier voltage (mV) - upper bound"
            ][0]
            outlierCutoff_rv = pv["PSP Outliers"][1]["replacement value"][0]
        else:
            outlierCutoff_LV = pv["PSC Outliers"][1][
                "Outlier voltage (pA) - lower bound"
            ][0]
            outlierCutoff_UV = pv["PSC Outliers"][1][
                "Outlier voltage (pA) - upper bound"
            ][0]
            outlierCutoff_rv = pv["PSC Outliers"][1]["replacement value"][0]

        if outlierCutoff_rv == "bound":
            data[data <= outlierCutoff_LV] = outlierCutoff_LV
            data[data >= outlierCutoff_UV] = outlierCutoff_UV
        elif outlierCutoff_rv == "median":
            data[data <= outlierCutoff_LV] = median_data
            data[data >= outlierCutoff_UV] = median_data
        elif outlierCutoff_rv == "mean":
            data[data <= outlierCutoff_LV] = mean_data
            data[data >= outlierCutoff_UV] = mean_data
        currentTrace = pv["Data selection"][1]["Trace"][0] - 1
        myPen = pg.mkPen(
            color=pg.intColor(currentTrace, hues=8)
        )  ## the pen to draw this
        self.events.tracePlot.plot(time, data, pen=myPen, name=trace.Label)
        self.verifyAllSweepEvents(data, time)

    def verifyAllSweepEvents(self, data, time):
        ## Deconvolution and peak detection
        with pg.ProgressDialog(
            "Detecting events. Please wait...", 0, 120, busyCursor=True, nested=False
        ) as dlg:
            D, upstroke_samples, downstroke_samples, fitPars = self.event_detection(
                useData=data, useTime=time
            )
            dlg += 10
            p0 = self.events.tracePlot
            p0.setLabel("bottom", "", units="")
            p0.showAxis("bottom", False)
            if self.events.decovPlot != "":
                self.events.peakIndex = []  ## the actual peak index
                plotHandle = self.events.decovPlot
                # self.events.threholdLine.disconnect()
                plotHandle.clear()
                try:
                    plotHandle.removeItem(self.events.threholdLine)
                except:
                    print("line 2737: nothing to be removed")
                try:
                    plotHandle.removeItem(self.events.threholdLine_raw)
                except:
                    print("")
                self.events.decovPlot = ""
                try:
                    self.event_traceView.removeItem(plotHandle)
                except:
                    print("line 2742: Nothing to be removed")
                try:
                    l1, l2 = self.getEventLine()
                    if l2 != []:
                        l = l1.extend(l2)
                    else:
                        l = l1
                    for x in l:
                        p0.removeItem(x)
                except:
                    print("line 2752: Nothing to be removed")
            plotHandle = self.event_traceView.addPlot(row=1, col=0)
            plotHandle.showGrid(x=True, y=True)

            # plotHandle = self.event_traceView.addPlot(row=1, col=0)
            def mouseMoved(pos):
                tracePlot = self.events.decovPlot
                mousePoint = tracePlot.vb.mapSceneToView(pos)
                vLine.setPos(mousePoint.x())

            def mouseClicked(evt):
                tracePlot = self.events.tracePlot
                pos = evt.scenePos()  ## this is for MouseClicked event
                # button = evt.button()
                # pdb.set_trace()
                if tracePlot.sceneBoundingRect().contains(pos):
                    mousePoint = tracePlot.vb.mapSceneToView(pos)
                    index = int(mousePoint.x() * self.parameters["fs"])
                    # print(mousePoint.x(), self.events.data[index])
                    # print('X index', index)
                    time = self.events.time
                    # pdb.set_trace()
                    if index >= 0 and index < len(time):
                        """ADDing template selection twin-lines"""
                        if (
                            evt.modifiers() == QtCore.Qt.ShiftModifier
                            and evt.button() == QtCore.Qt.LeftButton
                        ):
                            vpeakline = EventMark(
                                mousePoint.x(), index, penc="r", userGen=True
                            )
                            p0.addItem(vpeakline)

            vLine = pg.InfiniteLine(angle=90, movable=False, name="cursorH_D")
            # baselineLine = pg.InfiniteLine(angle=0, pen = pg.mkPen(color = 'w', style=QtCore.Qt.DashLine),
            #                                hoverPen = 'm', movable=False, name='baseline')
            pv = self.event_pd_view.p.getValues()
            # isFitSlope = pv['Peak parameters for raw trace'][1]['Linear fit during onset'][0]
            nSTD = pv["Peak parameters for deconvoled trace"][1]["Threshold (stdev)"][0]
            m, th, std0 = self.calBaselineValue(D, nSTD)

            g = plotHandle.plot(
                self.events.time, D, pen=pg.mkPen("w"), name="Deconvolution"
            )
            dlg += 10
            ## parameters for detecting peaks from the deconvoluted trace
            miniDist = int(pv["Peak parameters for deconvoled trace"][1]["Distance"][0])
            wlen = int(
                pv["Peak parameters for deconvoled trace"][1]["Wlen (samples)"][0]
            )  ## can use template window length
            if self.events.traceYUnit == "A":
                prominence = int(
                    pv["Peak parameters for deconvoled trace"][1]["Prominence (pA)"][0]
                )
            else:
                prominence = int(
                    pv["Peak parameters for deconvoled trace"][1]["Prominence (mV)"][0]
                )
            width_min = int(
                pv["Peak parameters for deconvoled trace"][1]["Width"][1]["minimal"][0]
            )
            width_max = int(
                pv["Peak parameters for deconvoled trace"][1]["Width"][1]["maximal"][0]
            )
            # pdb.set_trace()
            peak_onset_idx, peak_heights_D, peak_props = fitFuncs.getPeaks(
                D,
                height=th,
                width=(width_min, width_max),
                wlen=wlen,
                distance=miniDist,
                prominence=prominence,
            )

            data_ = self.events.data.copy()
            data_[np.abs(data_) > np.mean(data_) + 3 * np.std(data_)] = np.mean(data_)
            raw_baseline = np.mean(data_)
            self.events.globalBaseline = raw_baseline
            # get raw peaks from original trace
            x_th = [self.events.time[0], self.events.time[-1]]
            deltD = np.median(np.diff(self.events.data))
            peak_nSD = pv["Peak parameters for raw trace"][1]["Peak Threshold (stdev)"][
                0
            ]
            waveformDuration = pv["Peak parameters for raw trace"][1][
                "Waveform post peak duration"
            ][0]

            peak_threh = raw_baseline + peak_nSD * np.sign(fitPars[2]) * np.std(data_)
            p0.plot(
                x_th,
                [peak_threh, peak_threh + deltD],
                pen=pg.mkPen("r", style=QtCore.Qt.DashLine),
            )
            p0.plot(
                x_th,
                [raw_baseline, raw_baseline + deltD],
                pen=pg.mkPen("w", style=QtCore.Qt.DashLine),
            )
            peak_idx, selected_index = fitFuncs.getRawPeaks3(
                self.events.data, D, peak_onset_idx, peak_threh, np.sign(fitPars[2])
            )
            self.events.peakIndex = peak_idx  ## the actual peak index auto detected
            self.events.peakStart = [peak_onset_idx[j] for j in selected_index]
            self.events.downstroke_samples = int(
                waveformDuration * self.parameters["fs"]
            )
            self.events.upstroke_samples = upstroke_samples * 4

            # peak_idx = [int(peakIdx_) for peakIdx_ in peak_props['right_ips']]
            xunit = self.events.traceXUnit
            yunit = self.events.traceYUnit
            plotHandle.setLabels(bottom=("Time", xunit), left=(yunit))

            plotHandle.setXLink(p0)
            p0.setZValue(0)
            g.setZValue(0)
            plotHandle.addItem(vLine, ignoreBounds=True)
            # nPeaks = len(peak_idx0)
            # y_peaks = np.ones_like(time[peak_idx])*np.median(self.events.data[peak_idx])
            y_peaks = self.events.data[peak_idx]
            ## much quicker!
            invalidColor = pg.mkColor(0.5)

            def marker_mouseClicked(scatterPlot, spots):
                penColor = spots[0].pen().color().getRgb()
                if penColor[0] > 0:
                    pc = "g"  # change r to g
                else:
                    pc = invalidColor  ## change g to r
                spotsPen = [scatterPlot.getSpotOpts(yy)[2] for yy in scatterPlot.data]
                spotsBrush = [scatterPlot.getSpotOpts(yy)[3] for yy in scatterPlot.data]
                spotsPen[spots[0].index()] = pg.mkPen(pc)
                spotsBrush[spots[0].index()] = pg.mkBrush(pc)
                scatterPlot.setPen(spotsPen)
                scatterPlot.setBrush(spotsBrush)

            def marker_mouseHovered(scatterPlot, spot):  # update wave plot
                if len(spot) > 0:
                    spotIdx = spot[0].index()
                    onsetIdx = self.events.peakStart[spotIdx]
                    peakIdx = self.events.peakIndex[spotIdx]
                    donwstrokeIdx = peakIdx + self.events.downstroke_samples
                    upstrokeIdx = peakIdx - self.events.upstroke_samples
                    self.event_viewSingleEvents(
                        upstrokeIdx, onsetIdx, peakIdx, donwstrokeIdx, noMove=True
                    )

            if self.events.eventMarker != None:
                self.events.eventMarker.clear()
                p0.removeItem(self.events.eventMarker)
                self.events.eventMarker = None

            eventMarker = pg.ScatterPlotItem(
                x=time[peak_idx],
                y=y_peaks,
                pen="g",
                hoverPen="r",
                brush="g",
                hoverBrush="r",
                symbol="o",
                size=6,
                pxMode=True,
                hoverable=True,
            )
            eventMarker.sigClicked.connect(marker_mouseClicked)
            eventMarker.sigHovered.connect(marker_mouseHovered)

            p0.addItem(eventMarker)
            self.events.eventMarker = eventMarker
            x_th = [self.events.time[0], self.events.time[-1]]
            deltD = np.median(np.diff(self.events.data))
            y_th = [th, th + deltD]
            threholdLine = plotHandle.plot(
                x_th, y_th, pen=pg.mkPen("r", width=1.5, style=QtCore.Qt.DashLine)
            )
            threholdLine.setZValue(1)
            self.events.threholdLine = threholdLine
            self.events.threholdPar = self.event_pd_view.p.param(
                "Peak parameters for deconvoled trace"
            ).children()[0]
            plotHandle.scene().sigMouseMoved.connect(mouseMoved)
            plotHandle.scene().sigMouseClicked.connect(mouseClicked)
            plotHandle.setYRange(D.min(), D.max())
            self.events.decovPlot = plotHandle
            dlg.setValue(120)  # done
            self.events.sweepType = "all"  ## or 'single'
            if dlg.wasCanceled():
                raise Exception("Processing canceled by user")

    def plotConcatenatedSP(self, toConcat=None):
        """special functions for plotting concatenated spontaus traces"""
        # pdb.set_trace()
        datName, series_index = self.checkNodeIndex()
        if series_index == "":
            print("Nothing to process!")
            return
        plotHandle = self.events.tracePlot  # self.event_traceView.getItem(0, 0)
        if plotHandle == []:
            return

        pv = self.eventParTree_data_view.p.getValues()
        currentSweep = pv["Data selection"][1]["Sweep"][0]
        currentTrace = pv["Data selection"][1]["Trace"][0] - 1
        if toConcat == None:
            isConcat = pv["Data selection"][1][
                "Concatenate all sweeps for current channel"
            ][0]
        else:
            isConcat = toConcat
        lowCutOff = pv["Data selection"][1]["Low frequency cutoff"][0]
        highCutOff = pv["Data selection"][1]["High frequency cutoff"][0]
        isRemoveSpike = pv["Spikes"][1]["Removing spikes"][0]
        if self.currentPulseTree.filetype == ".dat":
            (
                bundleClass,
                stimChanLabels,
                serieIndex,
                treeChildrenIdx,
            ) = self.getDatFilesFromTreeWiget("Spontaneous", datName, series_index)
            if bundleClass == []:
                print(f"Current dat file: {self.currentPulseTree.bundle.file_name}!")
                # pdb.set_trace()
                bundleClass = HekaBundleInfo(
                    self.currentPulseTree.bundle.file_name
                )  ## a bundle with some extra functions
                nSweep = bundleClass.countSweeps(series_index)
                self.parameters["fs"] = bundleClass.getSeriesSamplingRate(
                    series_index
                )  ##
                self.events.setSeriesSampleRate(
                    self.parameters["fs"]
                )  ## series sampling rate

            else:
                nSweep = bundleClass.countSweeps(series_index)
                self.events.node = self.selTreesSP.selectedItems()[0]
                self.events.setSeriesSampleRate(
                    bundleClass.getSeriesSamplingRate(serieIndex)
                )  ## series sampling rate
                self.parameters["fs"] = bundleClass.getSeriesSamplingRate(
                    serieIndex
                )  ## this need to be change!
                # self.currentPulseTree= self.selTreesSP.selectedItems()[0]
                self.currentPulseTree.bundle = bundleClass.bundle
        else:
            block = self.currentPulseTree.abfBlocks[
                series_index[1]
            ]  ## choose the block selected
            ## loops through all the sweeps.
            nSweep = len(block.segments)
            self.events.setSeriesSampleRate(
                self.parameters["fs"]
            )  ## series sampling rate

        self.events.nSweep = nSweep
        self.events.datFile = datName  # self.currentPulseTree.dat_file
        self.events.seriesName = series_index[1] + 1
        self.events.traceName = currentTrace + 1

        if isConcat:
            SweepIdx = range(nSweep)
            self.events.isConcat = True
        else:
            if currentSweep > nSweep:
                print("Exceeding available number of sweeps")
                return
            else:
                SweepIdx = [currentSweep - 1]
                self.events.isConcat = False

        data_ = []
        time_ = []
        if self.currentPulseTree.filetype == ".dat":  ## reading '.dat' files
            with pg.ProgressDialog(
                "Concatenating traces",
                maximum=len(SweepIdx),
                busyCursor=True,
                nested=False,
            ) as dlg:

                seriesIdx = list(series_index.copy())  ## get series level index
                seriesIdx.append(0)  ## sweep 0
                seriesIdx.append(currentTrace)
                # plotHandle = self.event_traceView.addPlot(col = 0, colspan=4)  ## add a subplot on the graphic layout

                ## loops through all the sweeps
                timeoffset = 0
                setLabel = True
                for sweep in SweepIdx:
                    seriesIdx[2] = sweep  ## change sweeep level index
                    myPen = pg.mkPen(
                        color=pg.intColor(currentTrace, hues=8)
                    )  ## the pen to draw this
                    try:
                        trace = bundleClass.bundle.pul[seriesIdx[0]][seriesIdx[1]][
                            seriesIdx[2]
                        ][
                            seriesIdx[3]
                        ]  ## get trace meta information
                    except:
                        print("Invalid sweep index")
                        return
                    if not isConcat:
                        title = ""
                    else:
                        title = " all sweeps"
                    # print(f'High cutoff: {highCutOff}')
                    p, deltaT, time, data = self.plotSingleTrace(
                        plotHandle,
                        trace,
                        seriesIdx,
                        myPen,
                        setLabel,
                        False,
                        False,
                        timeoffset,
                        highCutOff,
                        title,
                        isRemoveSpike,
                        lowCutOff,
                    )
                    timeoffset += deltaT
                    time_.extend(time)

                    mean_data = np.mean(data)
                    median_data = np.median(data)
                    if trace.YUnit == "V":
                        outlierCutoff_LV = pv["PSP Outliers"][1][
                            "Outlier voltage (mV) - lower bound"
                        ][
                            0
                        ]  ## lower bound:data point beyond this would be replace with global mean
                        outlierCutoff_UV = pv["PSP Outliers"][1][
                            "Outlier voltage (mV) - upper bound"
                        ][
                            0
                        ]  ## upper bound:data point beyond this would be replace with global mean
                        outlierCutoff_rv = pv["PSP Outliers"][1]["replacement value"][
                            0
                        ]  # replacement type
                    else:
                        outlierCutoff_LV = pv["PSC Outliers"][1][
                            "Outlier voltage (pA) - lower bound"
                        ][
                            0
                        ]  ## lower bound:data point beyond this would be replace with global mean
                        outlierCutoff_UV = pv["PSC Outliers"][1][
                            "Outlier voltage (pA) - upper bound"
                        ][
                            0
                        ]  ## upper bound:data point beyond this would be replace with global mean
                        outlierCutoff_rv = pv["PSC Outliers"][1]["replacement value"][
                            0
                        ]  # replacementt type

                    if outlierCutoff_rv == "bound":
                        data[data <= outlierCutoff_LV] = outlierCutoff_LV
                        data[data >= outlierCutoff_UV] = outlierCutoff_UV
                    elif outlierCutoff_rv == "median":
                        data[data <= outlierCutoff_LV] = median_data
                        data[data >= outlierCutoff_UV] = median_data
                    elif outlierCutoff_rv == "mean":
                        data[data <= outlierCutoff_LV] = mean_data
                        data[data >= outlierCutoff_UV] = mean_data
                    data_.extend(data)
                    dlg += 1
                    if dlg.wasCanceled():
                        print("Canceled stage %s" % sweep)
                        break
            traceXUnit = trace.XUnit
            traceYUnit = trace.YUnit
        elif self.currentPulseTree.filetype == ".abf":
            print("reading abf files")
            for idx, sweepIdx_ in enumerate(SweepIdx):  # enumerate(block.segments):
                myPen = pg.mkPen(
                    color=pg.intColor(idx, hues=3)
                )  ## the pen to draw this
                trace = (
                    block.segments[sweepIdx_].analogsignals[0].transpose()[0]
                )  ## numpy array for current sweep  ## get trace meta information
                segmentIdx = [series_index[1], sweepIdx_]
                data, time = self.plotSingleTrace_ABF(
                    plotHandle,
                    segmentIdx,
                    trace,
                    myPen,
                    title=False,
                    scaleBar=False,
                    analaysisSpike=False,
                    plotStim=False,
                )
                data_.extend(data)  ## filtered version
                time_.extend(time)
            # time_, data_ = self.extractSingleSeries_ABF(series_index)
            traceXUnit = self.currentPulseTree.abf.xUnits
            traceYUnit = self.currentPulseTree.abf.yUnits

        if isConcat:
            self.events.currentSweep = 0
        else:
            self.events.currentSweep = currentSweep
        self.events.traceXUnit = traceXUnit
        self.events.traceYUnit = traceYUnit
        self.events.time = np.array(time_)
        out_ = np.array(data_)
        self.events.data = out_
        self.events.dataMin = out_.min()
        self.events.dataMax = out_.max()

        ## update paramters for events object
        ## make it a little bit pretty
        plotHandle2 = self.event_traceView2.getItem(0, 0)
        plotHandle2.showGrid(x=True, y=True)
        if traceYUnit == "mV":
            traceYUnit = "V"
        plotHandle2.setLabels(bottom=("Time", traceXUnit), left=(traceYUnit))

    def calBaselineValue(self, data, nSTD):
        if self.events.std == 0:
            hist, bin_edges = np.histogram(data, bins=25, density=False, normed=False)
            bindex = np.argmax(hist)
            m = bin_edges[bindex]
            std0 = np.std(data)
            self.events.m = m
            self.events.std = std0
        else:
            m = self.events.m
            std0 = self.events.std
        # return , (bin_edges[bindex]+ M)*0.4
        # x = data[(data > data.min()*0.6) & (data < data.max()*0.6)]
        # m = np.mean(x)
        print(m, nSTD, std0)
        t = m + nSTD * std0
        return m, t, std0

    def setMarkerLineValue(self, lineName, p):
        l = self.getMarkerline(lineName)
        l.setPos(p)  ## set baseline position

    def getMarkerLineValue(self, lineName):
        l = self.getMarkerline(lineName)
        return l.pos().y()  ## just for horizontal line!

    def getMarkerline(self, lineName):
        for l in self.events.decovPlot.vb.allChildItems():
            if hasattr(l, "name"):
                if l.name() == lineName:
                    return l

    def choose2PlotSP(self):
        selected = self.currentPulseTree.selectedItems()
        sel = selected[0]
        # self.plotConcatenatedSP(sel, [0, 2, 5]) ## SP trace will be automatically concatenated

    def showdialog(self, msgText, details=""):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText(msgText)
        msg.setWindowTitle("Help information")
        msg.setDetailedText(details)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

    def License_clicked(self):
        DATA_PATH = os.path.join(patchview_dir, "Data", "LICENSE.txt")
        with open(DATA_PATH) as f:
            BSD_3 = f.readlines()
        BSD_3b = ""
        for l in BSD_3:
            BSD_3b = BSD_3b + l
        self.showdialog(
            "This program is under BSD-3 license.\nCopyright (c) 2020-2022, Ming Hu. All rights reserved.",
            BSD_3b,
        )
    def link2Doc_clicked(self):
         self.showdialog(
            "For source code & documentation, see details below: ",
            "Source code & issue reporting:\n https://github.com/ZeitgeberH/patchview \n\nDocumentation:\n https://patchview.readthedocs.io/en/latest/index.html",
        )       
    def makeNWBobject(self):
        selected = self.currentPulseTree.selectedItems()
        sel = selected[0]
        selIdx = sel.index
        if len(selIdx) < 2:
            return
        nwbObj = dat2NWB(self.currentPulseTree.dat_file, [selIdx[0], selIdx[1]])
        return nwbObj

    def getSaveFileNameNWB(self):
        """save file dialog. Return the path and filename for saving user selected data"""
        fileName = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save File",
            self.root + self.QtFileNameLabel.text()[:-4],
            "NWB (*.nwb)",
        )
        return fileName[0]

    def saveNWB_clicked(self):
        fileName = self.getSaveFileNameNWB()
        if fileName !='':
            self.nwbObj = self.makeNWBobject()
            if self.nwbObj !=None:
                self.statusBar.showMessage(" " + fileName, 3000)
                self.nwbObj.saveNBF(fileName)

    def exit_clicked(self):
        self.reset()
        self.close()

    #        self.app.closeAllWindows()

    def reset(self):
        """this wipe out every thing. Be careful"""
        self.currentSelectedIndex = []
        self.trace_view.clear()  ## refresh the layout
        self.trees_view.setCurrentIndex(0)
        self.cellFeatures_view.clear()
        self.sweepFeatures_view.clear()
        self.spikeTable_view.clear()
        self.sliceView.tables["ROI list"].clear()
        self.splitViewTab_connection.reset()
        self.splitViewTab_FP.reset()
        self.spikes_view.clf()
        self.spikes_view.draw()
        self.spikeSummary_view.clf()
        self.spikeSummary_view.draw()
        for t in self.PulseTrees:
            t.clear()  ## clear all trees views
        self.currentPulseTree = self.PulseTrees[0]
        self.selTrees.clear()
        self.plotItemList = []  ## collect all plot handles
        self.fpParTree_data_initilized = False
        self.EphyFeaturesObj = []  ## clear ephy object
        if self.showSpikePanelsFlag == 1:
            self.hideSpikePanels()
        self.pia = None
        self.neuronsData = None
        self.neuronMorph =  None 
        self.spikeTableSavePath = ""

    def resetAll_clicked(self):
        self.reset()

    def gatherParameters(self, ext=".dat"):
        self.parameters = {}
        self.setUserProfile([])
        if ext == ".dat":
            trace = self.currentPulseTree.bundle.pul[0][0][0][0]
            self.parameters["XUnit"] = trace.XUnit
            if trace.XUnit == "s":
                self.parameters["fs"] = 1.0 / trace.XInterval  ## sampling rate
        elif ext == ".abf":
            self.parameters["XUnit"] = "s"
            self.parameters["fs"] = self.currentPulseTree.abf.get_signal_sampling_rate()

    def setUserProfile(self, user=None):
        """read the user profile and set the active file as the last user"""
        if user == None:
            self.user = "MH"
        else:
            self.user = user
        self.useConfigDir = user_config_dir(appname)
        self.loadUserParamters(user)

    def loadUserParamters(self, user):
        confilePath  = os.path.join(self.useConfigDir,"patchview.yaml")
        if not os.path.exists(self.useConfigDir):
            import shutil 
            os.makedirs(self.useConfigDir)
            DATA_PATH = os.path.join(patchview_dir, "Data", "patchview.yaml")    
            shutil.copy(DATA_PATH, confilePath)

        pars = loadYAML(confilePath)
        parameters = {}
        # parameters['HF'] = pars['Filters']['High cutoff']
        parameters["filter_option"] = pars["Filters"]["Option"]
        parameters["HF"] = pars["Filters"]["High cutoff"]
        parameters["CleanUp"] = pars["CleanUp"]
        parameters["Protocols"] = pars["Protocols"]
        parameters["Protocols"]["This serie"] = {"Type": [], "StimChanID": []}
        parameters["plotStim"] = pars["Visulization"][
            "Plot stim"
        ]  ## plot stimulation current
        parameters["Downsampling"] = pars["Visulization"]["Downsampling"]
        parameters["RootDir"] = pars["RootDir"]
        self.root = pars["RootDir"]
        self.parameters = parameters
        # self.getEphySpikePars()
        return parameters

    def updateUserParamters(self):
        ## this is the parameters changed by the user in the global paramter tree box

        pv = self.globalSettings.p.getValues()
        self.parameters["Notch"] = pv["data preprocessing"][1][
            "Notch filter frequency"
        ][0]

        self.parameters["HF"] = pv["data preprocessing"][1]["High frequency cutoff"][0]
        self.parameters["CleanUp"]["minimalSweeps"] = pv["data preprocessing"][1][
            "minimal number of sweeps"
        ][0]

    def switchBackground_clicked(self):
        if self.trace_view._background == "k":
            self.trace_view.setBackground("w")
        else:
            self.trace_view.setBackground("k")

    def showSpikePanels(self):
        self.showSpikePanelsFlag = 1
        self.spikes_view.show()
        self.spikeSummary_view.show()
        self.trace_view2.show()
        self.plot_splitter.setStretchFactor(0, 1)
        self.plot_splitter.setStretchFactor(1, 1)
        self.topPlot_splitter.setStretchFactor(0, 2)
        self.topPlot_splitter.setStretchFactor(1, 1)
        self.OptionAction3.setChecked(True)

    def hideSpikePanels(self):
        self.showSpikePanelsFlag = 0
        self.spikes_view.hide()
        self.spikeSummary_view.hide()
        self.trace_view2.hide()
        self.OptionAction3.setChecked(False)

    def removeStimArtifacts_clicked(self):
        if self.removeArtifact == 1:
            self.OptionAction5.setChecked(True)
            self.removeArtifact = 0
            # print('38',self.OptionAction5.isChecked())
            return
        else:
            self.OptionAction5.setChecked(False)
            self.removeArtifact = 1
            # print('43',self.OptionAction5.isChecked())
            return

    # def FiringPattern3D_clicked(self):
    #     self.OptionAction4.setChecked(not self.OptionAction4.isChecked())

    def calcuatedvdt2(self, checkstatus):
        self.OptionAction6.setChecked(checkstatus)
        print(checkstatus)

    def setbatchFPana(self, checkstatus):
        self.OptionAction7.setChecked(checkstatus)
        print(checkstatus)

    def dvdt2_clicked(self):

        if self.dv2dt2:
            self.calcuatedvdt2(False)
            self.dv2dt2 = False
        else:
            self.calcuatedvdt2(True)
            self.dv2dt2 = True

    def batchFPana_clicked(self):

        if self.batchFPana:
            self.setbatchFPana(False)
            self.batchFPana = False
        else:
            self.setbatchFPana(True)
            self.batchFPana = True

    def hideSpikePanels_clicked(self):
        if self.showSpikePanelsFlag == 1:
            self.hideSpikePanels()
        else:
            self.showSpikePanels()

    def getSaveFileName(self):
        """save file dialog. Return the path and filename for saving user selected data"""
        fileName = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save File",
            self.root + self.QtFileNameLabel.text()[:-4] + "_PV",
            "pickle (*.dat)",
        )
        return fileName[0]

    def save_clicked(self):
        selected = self.currentPulseTree.selectedItems()
        if len(selected) < 1:
            return
        sel = selected[0]
        # plot selected dependent on what level is selected
        for sel in selected:
            treeLevel = len(sel.index)  ## group 1; series, 2; sweep: 3; trace: 4
            if treeLevel < 1:
                return  ## do nothing at group level
            if treeLevel == 2:  ##  series level: all sweeps, all channels
                self.checkSeriesType(sel)  ## check protocol type for this selection
                self.saveSingleSeries(sel)
            if treeLevel == 3:  ##  sweep level: all channels at this sweep
                pass
                self.saveSingleSweep(sel)
            if treeLevel == 4:  ##  series level: all sweeps, all channels
                self.saveSingleTrace(
                    sel.node, sel.index
                )  ## third arugument is for pen. if empty use default pen

    def saveSingleSeries(self, sel):
        fileName = self.getSaveFileName()
        import pickle

        if fileName != "":
            self.statusBar.showMessage("Saving " + fileName[:-4], 3000)
            (
                time,
                data,
                stimTime,
                stimData,
                stimInfo,
                serieInfo,
            ) = self.extractSingleSeries(sel)
            pdata = {
                "time": time,
                "data": data,
                "stimTime": stimTime,
                "stimData": stimData,
                "stimInfo": stimInfo,
            }
            with open(fileName[:-4], "wb") as f:
                pickle.dump(pdata, f)
                print("Data saved")

    def updateRecordingParameterTable(self, stimInfo, traceInfo):
        data = np.array(
            [
                ("Recording mode", traceInfo["Recording mode"], ""),
                ("Holding Current", stimInfo[0]["Vholding"], "pA"),
                ("Injected Current", stimInfo[1]["amplitude"], "pA"),
                ("Stimulus Duration", stimInfo[1]["duration"] * 1e3, "ms"),
                (
                    "Seal Resistance",
                    stimInfo[1]["sealresistance"] / 1e9,
                    "G" + "\u03A9",
                ),
                ("Stimulus->Dac", stimInfo[0]["StimToDac"], ""),  ##
                ("X-inteval", stimInfo[0]["sampleInteval"], "second"),
            ],
            dtype=[("Key", object), ("Value", object), ("Unit", object)],
        )
        self.parameter_Tab.setData(data)

    #        self.parameter_Tab.show()

    def extractStimData_ABF(self):
        rawP = self.currentPulseTree.abf.read_raw_protocol()
        nProtocol = len(rawP[0])  ## should be the same as nSegments for firing pattern.
        stimChanIdx = 0  ## channel index where stimuli is applied to
        stimUnit = rawP[2][stimChanIdx]  ## 'pA'
        stimData = rawP[0]
        self.currentPulseTree.stimUnit = stimUnit
        stimInfo = []
        for stim in stimData:
            stim_diff = np.diff(stim[stimChanIdx])
            idx = np.argwhere(stim_diff != 0).flatten().tolist()
            if len(idx) == 0:
                stim_start = int(0.1 * self.parameters["fs"])
                stim_end = int(0.7 * self.parameters["fs"])
                idx = [stim_start + 1, stim_end + 1]
            else:
                stim_start = idx[0] + 1
                stim_end = idx[1] + 1
            stim_ = [
                [],
                {
                    "start": stim_start,
                    "end": stim_end,
                    "amplitude": stim[stimChanIdx][idx[0] + 1],
                    "Vholding": stim[stimChanIdx][0],
                    "sampleInteval": 1 / self.parameters["fs"],
                },
                [],
            ]

            stimInfo.append(stim_)
        return stimInfo

    def extractSingleSeries_ABF(self, selIdx):
        block = self.currentPulseTree.abfBlocks[selIdx[1]]  ## choose the block selected
        ## loops through all the sweeps.
        nSweep = len(block.segments)
        nSamples = self.currentPulseTree.abf.get_signal_size(0, 0)
        data = np.zeros((nSamples, nSweep))
        for idx, seg in enumerate(block.segments):  # enumerate(block.segments):
            data[:, idx] = seg.analogsignals[0].transpose()[0]
        time = np.arange(nSamples) / self.parameters["fs"]
        return time, data / 1000.0  ## convert it back to Volt for downstream analaysis

    def extractSingleSeries(self, sel):
        """Extract single series (multiple sweep) data"""
        print(self.parameters["Protocols"]["This serie"]["Type"])
        if self.parameters["Protocols"]["This serie"]["Type"] == "Firing pattern":
            stimChanIndex = self.parameters["Protocols"]["This serie"]["StimChanID"]
            traceID = self.getCellTraceID(sel, stimChanIndex)

            ## get number of sweeps
            #            nSweep = sel.node.NumberSweeps
            nSweep = len(sel.node.children)
            seriesIdx = list(sel.index.copy())  ## get series level index
            seriesIdx.append(0)  ## sweep 0
            seriesIdx.append(traceID)
            data0 = self.currentPulseTree.bundle.data[
                [seriesIdx[0], seriesIdx[1], 0, seriesIdx[3]]
            ]
            data = np.zeros((data0.shape[0], nSweep))
            stimTime0, stimData0, stimInfo0 = self.currentPulseTree.bundle.stim(
                seriesIdx
            )  ## assume only 1 channel is stimulating
            #            self.updateRecordingParameterTable(stimInfo0)
            stimData = np.zeros(
                (len(stimData0), nSweep)
            )  ## ndrarry, time X traces X sweep
            ## loops through all the sweeps
            stimInfo = []
            for sweep in range(nSweep):
                seriesIdx[2] = sweep  ## change sweeep level index
                trace = self.currentPulseTree.bundle.pul[seriesIdx[0]][seriesIdx[1]][
                    seriesIdx[2]
                ][
                    seriesIdx[3]
                ]  ## get trace meta information
                time, data[:, sweep] = self.extractSingleTrace(trace, seriesIdx)
                (
                    stimTime,
                    stimData[:, sweep],
                    stimLevel,
                ) = self.currentPulseTree.bundle.stim(
                    seriesIdx
                )  ## assume only 1 channel is stimulating
                stimInfo.append(stimLevel)

            return time, data, stimTime, stimData, stimInfo, sel.node
        else:
            # if self.parameters['Protocols']['This serie']['Type'] == 'Connection' or self.parameters['Protocols']['This serie']['Type'] == 'Spontaneous':
            #            nSweep = sel.node.NumberSweeps
            nSweep = len(sel.node.children)
            nTraces = len(sel.node.children[0].children)
            seriesIdx = list(sel.index.copy())  ## get series level index
            seriesIdx.append(0)  ## sweep 0
            seriesIdx.append(0)  ## trace 0
            data0 = self.currentPulseTree.bundle.data[
                [seriesIdx[0], seriesIdx[1], 0, 0]
            ]
            data = np.zeros(
                (data0.shape[0], nTraces, nSweep)
            )  ## ndrarry, time X traces X sweep
            stimTime0, stimData0, stimInfo0 = self.currentPulseTree.bundle.stim(
                seriesIdx
            )  ## assume only 1 channel is stimulating
            #            self.updateRecordingParameterTable(stimInfo0)

            stimData = np.zeros(
                (len(stimData0), nSweep)
            )  ## ndrarry, time X traces X sweep
            ## loops through all the sweeps
            stimInfo = []
            for sweep in range(nSweep):
                seriesIdx[2] = sweep  ## change sweeep level index
                for t in range(nTraces):
                    seriesIdx[3] = t  ## change sweeep level index
                    trace = self.currentPulseTree.bundle.pul[seriesIdx[0]][
                        seriesIdx[1]
                    ][seriesIdx[2]][
                        seriesIdx[3]
                    ]  ## get trace meta information
                    time, data[:, t, sweep] = self.extractSingleTrace(trace, seriesIdx)
                    (
                        stimTime,
                        stimData[:, sweep],
                        stimLevel,
                    ) = self.currentPulseTree.bundle.stim(
                        seriesIdx
                    )  ## assume only 1 channel is stimulating
                    stimInfo.append(stimLevel)
            return time, data, stimTime, stimData, stimInfo, sel.node

    def extractSingleTrace(self, trace, index):
        ### to extract single traces for exporting
        data = self.currentPulseTree.bundle.data[index]
        data = self.bandPass_signal(data)
        samplePeroid = 1 / self.parameters["fs"]
        time = np.linspace(
            trace.XStart, trace.XStart + samplePeroid * (len(data) - 1), len(data)
        )  ## trace.Interval varies for some reason
        return time, data

    def plotSpikes_clicked(self):
        if not self.EphyFeaturesObj:
            returnVal = self.getEphyFeatures()
        else:
            returnVal = 0
        if returnVal:
            self.EphyFeaturesObj.plot_alignedSpikes(mplWidget=self.spikes_view)

    def getDatAndePhyFeatures(self, seriesIdx0):
        ''' Extract voltage data  and stimuli data
        '''
        # bundle =  HEKA.Bundle(dat_file)
        bundle = self.currentPulseTree.bundle
        ## number of sweeps for this series
        nSweep = len(
            bundle.pul[seriesIdx0[0]].children[seriesIdx0[1]].children
        ) 
        if nSweep  <2: ## single channel recording
            self.showdialog("Not enough number of sweeps!")
            return None, None, None, None
        ## number of channels for this series
        nTraces = len(bundle.pul[seriesIdx0[0]].children[seriesIdx0[1]].children[0])
        if nTraces  <2: ## single channel recording
            self.showdialog("Only applicable for multichannel recordings!")
            return None, None, None, None
        if self.parameters["Protocols"]["This serie"]["StimChanID"] != []:
            stimChanIndex = (
                self.parameters["Protocols"]["This serie"]["StimChanID"] - 1
            )  ## this is the command channel
        else:
            selected = self.currentSelectedIndex
            sel = selected[0]
            stimChanIndex = self.searchForStimchan(sel.node.Label)
        seriesIdx = list(seriesIdx0.copy())  ## get series level index
        seriesIdx.append(0)  ## sweep 0
        seriesIdx.append(stimChanIndex)  ## trace 0

        _, _, stimInfo = self.currentPulseTree.bundle.stim(
            seriesIdx
        )  ## assume only 1 channel is stimulating

        data0 = bundle.data[[seriesIdx[0], seriesIdx[1], 0, 0]]
        data = np.zeros(
            (data0.shape[0], nTraces, nSweep)
        )  ## ndrarry, time X traces X sweep
        ## loops through all the sweeps
        for sweep in range(nSweep):
            seriesIdx[2] = sweep  ## change sweeep level index
            for t in range(nTraces):
                seriesIdx[3] = t  ## change trace level index
                da_ = bundle.data[seriesIdx]
                da_ = self.bandPass_signal(da_, 1000)
                data[:, t, sweep] = da_
                if t == 0:
                    trace = bundle.pul[seriesIdx[0]][seriesIdx[1]][seriesIdx[2]][
                        seriesIdx[3]
                    ]  ## get trace meta information
                    time = np.linspace(
                        trace.XStart,
                        trace.XStart + trace.XInterval * (len(da_) - 1),
                        len(da_),
                    )
        return time, data, stimInfo, stimChanIndex

    def gapJunctionAnalysisCorrelations(
        self, time, data, stimChanIndex, stimInfo, delay=0.001, windowLength=0.2
    ):
        """
        check correlations between cell voltage in stimulated channel and other simutaneouly recorded channel
        delay: staring time of window relative to stimulus onset
        windowLength: duration of window for calcuating CC
        """
        nChan = data.shape[1]
        nSweep = data.shape[2]
        stimStart = stimInfo[1]["start"] * stimInfo[1]["sampleInteval"]
        # print(f"Stimuli onset at {stimStart} second")
        t1 = stimStart + delay
        t2 = t1 + windowLength
        t0Idx = np.flatnonzero(time >= stimStart)[0]
        t1Idx = np.flatnonzero(time >= t1)[0]
        t2Idx = np.flatnonzero(time >= t2)[0]
        data_ = data[t1Idx:t2Idx, :, :]
        corrsR = np.zeros((nSweep, nChan))
        corrsP = np.zeros((nSweep, nChan))
        relativeDeflections = np.zeros((nSweep, nChan))

        for ichan in range(nChan):
            for iSweep in range(nSweep):
                baseline_stim = np.nanmean(data[:t0Idx, stimChanIndex, iSweep])
                baseline = np.nanmean(
                    data[:t0Idx, ichan, iSweep]
                )  ## get average baseline
                x1 = data_[:, ichan, iSweep] - baseline
                x2 = data_[:, stimChanIndex, iSweep] - baseline_stim
                corrsR[iSweep, ichan], corrsP[iSweep, ichan] = pearsonr(x1, x2)

                relativeDeflections[iSweep, ichan] = (
                    np.nanmean(data_[:, ichan, iSweep]) - baseline
                )

        return relativeDeflections, corrsR, corrsP

    def gapJunctionPlot(
        self, relativeDeflections, stimChanID, corrsR, option=0, mplWidget=None
    ):
        # time2, stim, stimInfo= self.currentPulseTree.bundle.stim(index)    ## assume only 1 channel is stimulating
        # stimStartSample,delaySamples, responseSamples, avg2, peak, peak_idx,significantTest, baseline, meanCorrelation =  self.calculateConnectionTraces(time, data, stimInfo)
        relativeDeflections = relativeDeflections * 1000
        nRow = 2  # nSignificant //nCol
        nCol = relativeDeflections.shape[1] // 2
        gs = mplWidget.figure.add_gridspec(nRow, nCol)
        colorMap = [
            pg.glColor(pg.intColor(j, hues=relativeDeflections.shape[1]))
            for j in range(relativeDeflections.shape[1])
        ]

        # print(data.shape[1])
        mplWidget.figure.clear()
        mplWidget.show()
        mplWidget.clf()
        if option == 0:
            data = relativeDeflections
        else:
            data = corrsR

        for j in range(data.shape[1]):
            iCol = np.mod(j, nCol)
            iRow = j // nCol
            if j != stimChanID:
                plt0 = mplWidget.figure.add_subplot(gs[iRow, iCol])

                x, y = data[:, stimChanID], data[:, j]
                r, p = pearsonr(x, y)
                # pdb.set_trace()
                plt0.plot(
                    relativeDeflections[:, stimChanID],
                    y,
                    color=colorMap[j],
                    marker=".",
                    markersize=6,
                    markerfacecolor=colorMap[j],
                    linestyle="None",
                )

                try:
                    m, b = np.polyfit(x, y, 1)
                    plt0.plot(x, m * x + b, "k")
                except Exception as e:
                    print(e)

                plt0.axhline(0, 0, 1.0, color="gray", linestyle="--")
                if p < 0.01:
                    fontColor = "r"
                else:
                    fontColor = "k"
                plt0.set_title(
                    "{0:d} vs {1:d}, r:{2:0.2f}, p:{3:.3f}".format(
                        stimChanID + 1, j + 1, r, p
                    ),
                    fontsize=12,
                    color=fontColor,
                )
                if iRow > 0:
                    plt0.set_xlabel(
                        "Avg."
                        + " \u0394"
                        + "voltage (mV) in Chan"
                        + str(stimChanID + 1),
                        fontsize=8,
                    )
                if option == 1:
                    plt0.set_ylabel("Correlation", fontsize=8)
                else:
                    plt0.set_ylabel("Avg." + " \u0394" + "mV", fontsize=8)
        mplWidget.figure.subplots_adjust(
            left=0.08, bottom=0.1, top=0.95, right=0.95, wspace=0.35, hspace=0.35
        )
        self.showSpikePanelsFlag = 1
        self.spikes_view.draw()

    def getEphyFeatures(self):

        selected = self.currentPulseTree.selectedItems()
        self.visulization_view.setCurrentIndex(1)
        if not len(selected):
            self.showdialog("Load a pulse tree first!")
            return 0
        sel = selected[0]

        if len(sel.index) != 2:
            self.showdialog("To extract spike fetures, please select a series!")
            return 0
        if self.currentPulseTree.filetype == ".dat":
            title = (
                self.sliceName[:-4]
                + " "
                + "Series"
                + str(sel.index[1] + 1)
                + " "
                + sel.node.Label
            )
            (
                time,
                data,
                stimTime,
                stimData,
                stimInfo,
                serieInfo,
            ) = self.extractSingleSeries(sel)
        elif self.currentPulseTree.filetype == ".abf":
            title = (
                self.currentPulseTree.dat_file[:-4].split("\\")[-1]
                + " "
                + "Block"
                + str(sel.index[0] + 1)
                + " Sweep"
                + str(sel.index[1] + 1)
            )
            time, data = self.extractSingleSeries_ABF(sel.index)
            stimInfo = self.currentPulseTree.abf_stimInfo

        (
            dv_cutoff1,
            min_height1,
            min_peak1,
            thresh_frac1,
            baseline_interval1,
            baseline_detect_thresh1,
            max_interval1,
            filterHighCutFreq_spike,
        ) = self.getEphySpikePars()
        ## Here we plug in all the available parameters
        self.EphyFeaturesObj = extraEhpys_PV.extraEphys(
            time,
            data,
            datIndex=sel.index,
            stimInfo=stimInfo,
            title=title,
            filterHighCutFreq=filterHighCutFreq_spike / 1000,
            dv_cutoff=dv_cutoff1,
            max_interval=max_interval1,
            min_height=min_height1,
            min_peak=min_peak1,
            thresh_frac=thresh_frac1,
            baseline_interval=baseline_interval1,
            baseline_detect_thresh=baseline_detect_thresh1,
        )

        # fileName = self.currentPulseTree.dat_file[:-4]+'_Series' + str(sel.index[1]+1) + '_'+ sel.node.Label+'_PvSpikeFeatures.pdf'
        self.plot_splitter.setStretchFactor(1, 3)
        self.topPlot_splitter.setStretchFactor(1, 1)
        self.trace_view2.clear()
        self.showSpikePanels()
        self.EphyFeaturesObj.plot_w_style(
            mplWidget=self.spikeSummary_view
        )  # filename = fileName)    ## plot some surmmary features
        self.EphyFeaturesObj.plot_alignedSpikes(
            mplWidget=self.splitViewTab_FP.matplotViews["Firing pattern"]
        )
        self.splitViewTab_FP.matplotViews["Firing pattern"].draw()

        sw = self.EphyFeaturesObj.df_related_features.copy()
        sw = sw[sw.spike_count>0]
        cellFeature = self.EphyFeaturesObj.Cell_Features
        for n in ['current','spike_count']:
            cellFeature[n]  = 0
            cellFeature[n]  = cellFeature[n].astype('object')
            cellFeature.at[0,n] = list(sw[n].values)

        self.updateEphyTable(
            cellFeature,
            self.splitViewTab_FP.tables["Cell features"],
        )
        spike_df = self.moveSweepCountFront(self.EphyFeaturesObj.df, "sweepCurrent")
        spike_df = self.moveSweepCountFront(spike_df, "sweepCount")
        self.updateEphyTable(spike_df, self.splitViewTab_FP.tables["Spike features"])
        sw = self.EphyFeaturesObj.df_related_features.copy()
        sw = self.moveSweepCountFront(sw, "current")
        sw.insert(loc=0, column="sweepCount", value=np.arange(1, len(sw) + 1))
        self.updateEphyTable(sw, self.splitViewTab_FP.tables["Sweep features"])
        self.splitViewTab_FP.bottomLeft_tabview.setCurrentIndex(1)
        return 1

    @staticmethod
    def moveSweepCountFront(df, feature):
        """feature: such as 'sweepCount'"""
        if feature in df.columns:
            sw = df[feature]
            df.drop(labels=[feature], axis=1, inplace=True)
            df.insert(0, feature, sw)
        return df

    def updateConnectionTable(self, df):
        df = np.array(
            df.to_records(index=False)
        )  ## format dataframe for using in QtTable wiget
        self.splitViewTab_connection.tables["Detected connections"].appendData(df)
        self.splitViewTab_connection.tables["Detected connections"].show()

    def updateEphyTable(self, df, table, title=None):
        ## convert dataframe into a record that can be populated into a QtTable
        if title is None:
            title = self.EphyFeaturesObj.title
        df.insert(0, "cellName", title)
        df = np.array(
            df.to_records(index=False)
        )  ## format dataframe for using in QtTable wiget
        table.appendData(df.copy())  ## update QtTable widget
        table.show()

    def extraEphyAction_clicked(self):
        self.updateUserParamters()
        pv = self.splitViewTab_FP.getParTreePars("Spike detection")
        resetTable = pv["Tables"][1]["Reset table automatically"][
            0
        ]  # this may different than general filter options
        if resetTable:
            self.splitViewTab_FP.tables["Sweep features"].clear()
            self.splitViewTab_FP.tables["Spike features"].clear()
            self.splitViewTab_FP.tables["Cell features"].clear()
        self.getEphyFeatures()

    def plotAveragedTraces(self, sel):
        '''
         plot sweep averaged traces
        '''
        selected = self.currentPulseTree.selectedItems()
        if not len(selected):
            self.showdialog("Load a pulse tree first!")
            return 0
        if len(sel.index) != 2:
            return 0
        self.checkSeriesType(sel)
        #        self.parameters['Protocols']['This serie']={'Type':'Firing pattern','StimChanID':stimChanIndex}
        if self.parameters["Protocols"]["This serie"]["Type"] in [
            "Connection",
            "Retina",
        ]:
            ## data should be ndarray, time X trace X sweep
            time, data, seriesIdx, traces = self.extractConnectionTraces(sel)
            #            df, seriesIdx = self.extractConnectionTraces(selected[0])
            #            avg = np.squeeze(np.mean(data, axis = 2))
            #            self.plotSeriesAverage(time, avg, seriesIdx)
            self.plotSeriesVariations(time, data, seriesIdx, traces)
        return 1

    def VisualAction2_clicked(self):
        """plot averaged traces at series level. when the protocol is for connection pattern"""
        if self.currentPulseTree.filetype == ".dat":
            datName, series_index = self.checkNodeIndex()
            if series_index == "":
                print("Nothing to process!")
                return
            (
                bundleClass,
                stimChanLabels,
                serieIndex,
                treeChildrenIdx,
            ) = self.getDatFilesFromTreeWiget("Spontaneous", datName, series_index)
            plotHandle = self.events.tracePlot  # self.event_traceView.getItem(0, 0)
            if plotHandle == []:
                return
            if bundleClass == []:
                print("data not found! line 3241")
                return
            nSweep = bundleClass.countSweeps(series_index)
        elif self.currentPulseTree.filetype == ".abf":
            datName, series_index = self.checkNodeIndex()
            time, data = self.extractSingleSeries_ABF(series_index)
            avg = np.mean(data, 1)
            pdb.set_trace()

    def gpAna(self, option):
        ''' Custom module for cross correlation between voltage of pair of channels in a multi-patch experiment
        '''
        selected = self.currentSelectedIndex
        sel = selected[0]
        if len(sel.index) != 2:
            return
        self.checkSeriesType(sel)
        if self.parameters["Protocols"]["This serie"]["Type"] in ["Firing pattern"]:
            time, data, stimInfo, stimChanIndex = self.getDatAndePhyFeatures(sel.index)
            if time is None: ## single channel recording
                return
            relativeDeflections, corrsR, corrsP = self.gapJunctionAnalysisCorrelations(
                time, data, stimChanIndex, stimInfo, delay=0.001, windowLength=0.2
            )
            self.gapJunctionPlot(
                relativeDeflections,
                stimChanIndex,
                corrsR,
                option=option,
                mplWidget=self.spikes_view,
            )

    def gp_clicked(self):
        self.gpAna(option=0)

    def gpc_clicked(self):
        self.gpAna(option=1)

    def checkGapJunctionAction_clicked(self):
        self.gpAna(option=0)

    def getDatFilesFromTreeWiget(self, sessionType: str, datName=None, seriesIdx=None):
        """
        Extract data from tree widget. unchecked items are ignored.
        """
        ## check if there's uncheck items
        ## also check duplicate items. Strategy to deal tih duplicate?
        tree = self.selTreesDispatcher[sessionType]
        # dat_files = tree.dat_files  ## this is a ditionary contain the full path and .dat file name
        path = tree.currentPath
        # keys = dat_files.keys() ## key is dat path+file name, value:list is series index
        # path, _ = os.path.split(list(keys)[0]) ## these dat lives in a same folder thus share path
        root = tree.invisibleRootItem()  ## root of current tree
        nSeries = root.childCount()
        usedBundles = []
        stimChanLabel = []
        serieIndex = []
        treeChildIdx = []  #             if cellDatKey == datName and cellLa
        for idx in range(nSeries):
            child = root.child(idx)
            cellLabel = child.text(0)  ## value in the columns of the widget
            cellDatKey = child.text(1)  ## value in the columns of the widget
            seriesidx = [int(s) for s in child.text(3).split(" ")]
            checkState = child.checkState(0) == QtCore.Qt.Checked
            if not checkState:
                print(f"**{cellLabel}** is not selected!")
                continue
            if datName == None:
                treeChildIdx.append(idx)
                usedBundles.append(
                    tree.bundleGroup[os.path.join(path, cellDatKey + ".dat")]
                )
                chanID = self.searchForStimchan(cellLabel)
                stimChanLabel.append(chanID)
                serieIndex.append([int(j) for j in seriesidx])
            else:  ## only spontaneous protocol need datName variable?
                if datName == cellDatKey and seriesidx == seriesIdx:
                    treeChildIdx = idx
                    usedBundles = tree.bundleGroup[
                        os.path.join(path, cellDatKey + ".dat")
                    ]
                    chanID = self.searchForStimchan(cellLabel)
                    stimChanLabel = chanID
                    serieIndex = [int(j) for j in seriesidx]
                    break

        return usedBundles, stimChanLabel, serieIndex, treeChildIdx

    def setFP_Data_selectionPars(
        self, tabName, bundleFiles, stimChanLabels, serieIndex, treeChildrenIdx
    ):
        if self.splitViewTab_FP.parsTree_values_init[tabName] == False:
            fp_par = deepcopy(AllMyPars.fp_analysis)
            fp_pars = []
            for idx, chanName in enumerate(stimChanLabels):
                bundle = bundleFiles[idx]
                sel = serieIndex[idx]
                nSweeps = bundle.countSweeps(sel)
                fp_par_ = deepcopy(fp_par[0])
                fp_par_[
                    "name"
                ] = f"{bundle.fileName[:-4]} Cell{str(chanName)} Ser.{str(sel[1]+1)}"
                fp_par_["children"][0]["limits"] = (1, nSweeps)
                fp_par_["children"][0]["value"] = nSweeps  ##  deploaroze sweep number
                fp_par_["children"][1]["value"] = int(
                    nSweeps / 2
                )  ## threhold sweep number
                fp_par_["children"][2]["value"] = 1  ## hyperpolize sweep number
                fp_pars.append(fp_par_)
            self.splitViewTab_FP.setParTreePars(
                fp_pars, 0, "Data selection", self.fpParTree_stateChange
            )

        else:
            fp_pars = []
            pars = self.splitViewTab_FP.getParTreePars(
                "Data selection"
            )  # self.fpParTree_data_view.p
            ## need to add a couple of sentence here to preserve changes
            for idx, chanName in enumerate(stimChanLabels):
                bundle = bundleFiles[idx]
                sel = serieIndex[idx]
                nSweeps = bundle.countSweeps(sel)
                cellName = (
                    f"{bundle.fileName[:-4]} Cell{str(chanName)} Ser.{str(sel[1]+1)}"
                )
                fp_par_ = deepcopy(AllMyPars.fp_analysis[0])
                fp_par_["name"] = cellName
                fp_par_["children"][0]["limits"] = (1, nSweeps)
                if cellName in pars.keys():
                    fp_par_["children"][0]["value"] = pars[cellName][1][
                        "Depolarize sweep"
                    ][
                        0
                    ]  ##  deploaroze sweep number
                    fp_par_["children"][1]["value"] = pars[cellName][1][
                        "Threshold sweep"
                    ][
                        0
                    ]  ## threhold sweep number
                    fp_par_["children"][2]["value"] = pars[cellName][1][
                        "Hyperpolar sweep"
                    ][
                        0
                    ]  ## hyperpolize sweep number
                else:
                    fp_par_["children"][0][
                        "value"
                    ] = nSweeps  ##  deploaroze sweep number
                    fp_par_["children"][1]["value"] = int(
                        nSweeps / 2
                    )  ## threhold sweep number
                    fp_par_["children"][2]["value"] = 1  ## hyperpolize sweep number
                fp_pars.append(fp_par_)
            self.splitViewTab_FP.setParTreePars(
                fp_pars, 0, "Data selection", self.fpParTree_stateChange
            )

    def fpParTree_stateChange(self, x, y):
        for idx, child in enumerate(x.children()):
            for gc in child.children():
                if gc == y[0][0]:
                    self.updateFP_batchPlot(idx, child)
                    return

    def updateFP_batchPlot(self, file_idx, child):
        """idx is the index where the plot should change"""
        pars = [p.value() for p in child.children()]
        (
            bundleFiles,
            stimChanLabels,
            serieIndex,
            treeChildrenIdx,
        ) = self.getDatFilesFromTreeWiget("Firing pattern")
        sel = serieIndex[file_idx]
        stimChanIndex = stimChanLabels[file_idx]
        nChan = len(stimChanLabels)
        nCol = 4
        if nChan % nCol != 0:
            nRow = nChan // nCol + 1
        else:
            nRow = nChan // nCol
        if nRow < 2:
            nRow = 2
        bundle = bundleFiles[file_idx]
        irow, icol = file_idx // nCol, file_idx % nCol
        axes = (
            self.splitViewTab_FP.matplotViews["Firing pattern"].getFigure().get_axes()
        )
        axes[file_idx].cla()
        nSweeps = bundle.countSweeps(sel)
        time = bundle.getSweepTimeStamps(sel)

        if pars[0] - 1 >= 0:
            line2Idx = sel.copy()
            line2Idx.extend([pars[0] - 1, stimChanIndex - 1])
            data2 = bundle.getSingleTraceData(line2Idx)
            ax = self.splitViewTab_FP.matplotViews["Firing pattern"].plotline(
                irow, icol, time * 1000.0, data2 * 1000.0, "r", "", 0.5, axes[file_idx]
            )

        if pars[1] - 1 >= 0:
            line0Idx = sel.copy()
            line0Idx.extend([pars[1] - 1, stimChanIndex - 1])
            data0 = bundle.getSingleTraceData(line0Idx)
            # axes = self.splitViewTab_FP.matplotViews['Firing pattern'].getFigure().get_axes()

            self.splitViewTab_FP.matplotViews["Firing pattern"].plotline(
                irow, icol, time * 1000.0, data0 * 1000.0, "g", "", 0.5, ax
            )

        if pars[2] - 1 >= 0:
            line1Idx = sel.copy()
            line1Idx.extend([pars[2] - 1, stimChanIndex - 1])
            data1 = bundle.getSingleTraceData(line1Idx)
            # title = 'cell '+str(stimChanIndex) + ' Ser.'+str(sel[1]+1)
            title = (
                f"{bundle.fileName[:-4]} Cell{str(stimChanIndex)} Ser.{str(sel[1]+1)}"
            )
            self.splitViewTab_FP.matplotViews["Firing pattern"].plotline(
                irow, icol, time * 1000.0, data1 * 1000.0, "b", title, 0.5, ax
            )
        # ax.set_xlabel('Time (mS)')
        # ax.set_ylabel('Voltage (mV)')
        self.splitViewTab_FP.matplotViews["Firing pattern"].draw()

    def checkEvent_action_clicked(self):
        self.makeNewEventWindow(True)
        self.visulization_view.setCurrentIndex(4)

    def checkFP_action_clicked(self):
        # debugInfo('', True)
        self.visulization_view.setCurrentIndex(1)
        (
            bundleFiles,
            stimChanLabels,
            serieIndex,
            treeChildrenIdx,
        ) = self.getDatFilesFromTreeWiget("Firing pattern")

        ## this needs to be updated!!
        self.setFP_Data_selectionPars(
            "Data selection", bundleFiles, stimChanLabels, serieIndex, treeChildrenIdx
        )
        nChan = len(stimChanLabels)
        nCol = 4
        if nChan % nCol != 0:
            nRow = nChan // nCol + 1
        else:
            nRow = nChan // nCol
        if nRow < 2:
            nRow = 2
        self.splitViewTab_FP.matplotViews["Firing pattern"].clf()
        axes = self.splitViewTab_FP.matplotViews["Firing pattern"].setPlotGrid(
            nRow, nCol
        )
        pars = self.splitViewTab_FP.getParTreePars(
            "Data selection"
        )  # self.fpParTree_data_view.p

        for file_idx, bundle in enumerate(bundleFiles):
            sel = serieIndex[file_idx]
            stimChanIndex = stimChanLabels[file_idx]
            irow, icol = file_idx // nCol, file_idx % nCol
            fs = bundle.getSeriesSamplingRate(sel)
            nSweeps = bundle.countSweeps(sel)
            time = bundle.getSweepTimeStamps(sel)
            # cellName = 'cell '+str(stimChanIndex) + ' Ser.'+str(sel[1]+1)
            cellName = (
                f"{bundle.fileName[:-4]} Cell{str(stimChanIndex)} Ser.{str(sel[1]+1)}"
            )
            ## plot three examplar traces: hypolarization, baseline-like, depolarization
            ## depolarziation
            if self.batchFPana:
                line2Idx = sel.copy()
                sweepIdx = pars[cellName][1]["Depolarize sweep"][0] - 1
                if stimChanIndex == None:
                    stimChanIndex = 1
                line2Idx.extend([sweepIdx, stimChanIndex - 1])
                data2 = bundle.getSingleTraceData(line2Idx)
                ax = self.splitViewTab_FP.matplotViews["Firing pattern"].plotline(
                    irow, icol, time * 1000.0, data2 * 1000.0, "r"
                )
                #  baseline-like
                line0Idx = sel.copy()
                sweepIdx = pars[cellName][1]["Threshold sweep"][0] - 1
                line0Idx.extend([sweepIdx, stimChanIndex - 1])
                data0 = bundle.getSingleTraceData(line0Idx)
                self.splitViewTab_FP.matplotViews["Firing pattern"].plotline(
                    irow, icol, time * 1000.0, data0 * 1000.0, "g", "", 0.5, ax
                )
                # hyperpolzation
                line1Idx = sel.copy()
                sweepIdx = pars[cellName][1]["Hyperpolar sweep"][0] - 1
                line1Idx.extend([0, stimChanIndex - 1])
                data1 = bundle.getSingleTraceData(line1Idx)
                self.splitViewTab_FP.matplotViews["Firing pattern"].plotline(
                    irow, icol, time * 1000.0, data1 * 1000.0, "b", cellName, 0.5, ax
                )

            stimSel = sel.copy()
            stimSel.extend([0, 0])
            stimInfo = []
            for j in range(nSweeps):
                stimSel[2] = j
                _, _, stimInfo_ = bundle.getStim(stimSel)
                stimInfo.append(stimInfo_)
            data = bundle.getSeriesData(sel)
            ## extract the channel that get stimuli
            if data.shape[1] >1:
                data = np.squeeze(data[:,stimChanIndex - 1,:])
            else:
                data = np.squeeze(data[:,0,:])
            if data.shape[0] < time.shape[0]:
                time = time[:data.shape[0]]
            else:
                ndiff = data.shape[0]-time.shape[0]
                time = np.hstack((time, time[1:ndiff+1]+time[-1]))
            assert time.shape[0]==data.shape[0], "number of samples not match!"
            # try:
            ephObj = self.doFPanalysis(time, data, sel, stimInfo, cellName)
            sw = ephObj.df_related_features.copy()
            sw = sw[sw.spike_count>0]
            cellFeature = ephObj.Cell_Features
            for n in ['current','spike_count']:
                cellFeature[n]  = 0
                cellFeature[n]  = cellFeature[n].astype('object')
                cellFeature.at[0,n] = list(sw[n].values)
            # pdb.set_trace()
            self.updateEphyTable(
                cellFeature,
                self.splitViewTab_FP.tables["Cell features"],
                ephObj.title,
            )


        
            self.updateEphyTable(sw, self.splitViewTab_FP.tables["Sweep features"], ephObj.title)

            print(cellName, "done!")
                # except:
                #     print(cellName, 'failed!')
        if self.batchFPana:
            self.splitViewTab_FP.matplotViews["Firing pattern"].adjustSubplot(
                left=0.06, right=0.98, top=0.95, bottom=0.06, wspace=0.2, hspace=0.6
            )
            ax.set_xlabel("Time (mS)")
            ax.set_ylabel("Voltage (mV)")
            self.splitViewTab_FP.matplotViews["Firing pattern"].draw()
        self.splitViewTab_FP.bottomLeft_tabview.setCurrentIndex(0)
        self.eventData_tabview.setCurrentIndex(2)
        # self.splitViewTab_FP.matplotViews['Firing pattern'].figure.tight_layout()

    def doFPanalysis(self, time, data, datIndex, stimInfo, title):
        (
            dv_cutoff1,
            min_height1,
            min_peak1,
            thresh_frac1,
            baseline_interval1,
            baseline_detect_thresh1,
            max_interval1,
            filterHighCutFreq_spike,
        ) = self.getEphySpikePars()
        ## Here we plug in all the available parameters

        ephyFeaturesObj = extraEhpys_PV.extraEphys(
            time,
            data,
            datIndex=datIndex,
            stimInfo=stimInfo,
            title=title,
            filterHighCutFreq=filterHighCutFreq_spike / 1000,
            dv_cutoff=dv_cutoff1,
            max_interval=max_interval1,
            min_height=min_height1,
            min_peak=min_peak1,
            thresh_frac=thresh_frac1,
            baseline_interval=baseline_interval1,
            baseline_detect_thresh=baseline_detect_thresh1,
        )
        return ephyFeaturesObj

    def checkConnectionAction_clicked(self):
        """
         Analysis connections in selected series in self.selTrees.
        Returns
        -------
        None.
        """
        self.visulization_view.setCurrentIndex(2)
        (
            bundleFiles,
            stimChanLabels,
            serieIndex,
            treeIdx,
        ) = self.getDatFilesFromTreeWiget(
            "Connection"
        )  ## A dictionary which stores name of dat files
        df = {}
        df["Dat file"] = []
        df["SerieIdx"] = []
        df["Source Chan"] = []
        df["Target Chan"] = []
        df["Baseline(mV)"] = []
        df["Peak value(mV)"] = []
        df["Delta(mV)"] = []
        df["Delay(mS)"] = []
        df["Trial consistency"] = []
        nSig = self.splitViewTab_connection.tables[
            "Detected connections"
        ].rowCount()  ## number of significant traces (estimated from connection table)
        print("nsig", nSig)

        print("nsig " + str(nSig))
        # print(data.shape[1])

        self.splitViewTab_connection.matplotViews["Average traces"].figure.clear()
        self.splitViewTab_connection.matplotViews["Average traces"].show()
        self.splitViewTab_connection.matplotViews["Average traces"].clf()

        nSig = 0
        allDats = []
        for dat_file_idx, bundle in enumerate(bundleFiles):
            sel = serieIndex[dat_file_idx]
            print(f"Processing {bundle.fileName}, series {sel}")
            fs = bundle.getSeriesSamplingRate(sel)
            time = bundle.getSweepTimeStamps(sel)
            data = bundle.getSeriesData(sel)
            data = filterDatSeries(data, fs, hcutFreq=1000)
            stimChanID = stimChanLabels[
                dat_file_idx
            ]  # self.searchForStimchan(bundle.getSeriesLabel(sel)) #, seriesIdx0)
            time_stim, stimData, stimInfo = bundle.getStim(
                sel + [0]
            )  ## assume only 1 channel is stimulating
            stimStartSample = stimInfo[1]["start"]
            # delaySamples, responseSamples, avg2, peak, peak_idx, significantTest, baselineMean,\
            #     meanCorrelation =  calculateConnectionTraces(time, fs, data, stimStartSample)
            d_ = calculateConnectionTraces(time, fs, data, stimStartSample)
            allDats.append(d_)
            for idx, j in enumerate(d_[5]):
                if (
                    j == 1 and idx != stimChanID - 1
                ):  ## do not iinclude stimulation channel itself
                    nSig += 1

        # pdb.set_trace()
        nCol = 1  ## five cols per row
        if np.mod(nSig, nCol) == 0:
            nRow = nSig // nCol + 1
        else:
            nRow = nSig // nCol + 1
        gs = self.splitViewTab_connection.matplotViews[
            "Average traces"
        ].figure.add_gridspec(nRow + 1, nCol)

        axes = []
        nodeColor = []
        arrowStyl = []
        xlimRange = [0.15, 0.8]
        colorMap = [pg.glColor(pg.intColor(j, hues=8)) for j in range(8)]
        G = nx.DiGraph()
        edges = {}
        nSig_count = 0
        pltCount = 0
        colIdx = 0
        nFiles = len(bundleFiles)
        for dat_file_idx, bundle in enumerate(bundleFiles):
            sel = serieIndex[dat_file_idx]
            # print(f'Processing {bundle.fileName}, series {sel}')
            fs = bundle.getSeriesSamplingRate(sel)
            time = bundle.getSweepTimeStamps(sel)
            data = bundle.getSeriesData(sel)
            # debugInfo('', True)
            data = filterDatSeries(data, fs, hcutFreq=1000)
            stimChanID = stimChanLabels[
                dat_file_idx
            ]  # self.searchForStimchan(bundle.getSeriesLabel(sel)) #, seriesIdx0)
            _, _, stimInfo = bundle.getStim(
                sel + [0]
            )  ## assume only 1 channel is stimulating
            stimStartSample = stimInfo[1]["start"]
            # delaySamples, responseSamples, avg2, peak, peak_idx, significantTest, baselineMean,\
            #     meanCorrelation =  calculateConnectionTraces(time, fs, data, stimStartSample)
            (
                delaySamples,
                responseSamples,
                avg2,
                peak,
                peak_idx,
                significantTest,
                baselineMean,
                meanCorrelation,
            ) = allDats[dat_file_idx]
            for idx, j in enumerate(significantTest):
                nSig_count += 1
                if idx == stimChanID - 1:
                    stimData = avg2[:, idx] * 1000.0
                if (
                    j == 1 and idx != stimChanID - 1
                ):  ## do not iinclude stimulation channel itself
                    df["Dat file"].append(bundle.fileName)
                    df["SerieIdx"].append(sel)
                    df["Target Chan"].append(idx + 1)
                    df["Source Chan"].append(stimChanID)
                    df["Delay(mS)"].append(
                        (
                            time[stimStartSample + delaySamples + peak_idx[idx]]
                            - time[stimStartSample]
                        )
                        * 1000
                    )
                    df["Peak value(mV)"].append(peak[idx] * 1000)
                    df["Delta(mV)"].append((peak[idx] - baselineMean) * 1000)
                    df["Baseline(mV)"].append(baselineMean * 1000)
                    df["Trial consistency"].append(meanCorrelation)
                    G.add_node(stimChanID)  # , color =  colorMap[stimChanID])
                    G.add_node(idx + 1)
                    edges[(stimChanID, idx + 1)] = peak[idx] - baselineMean
                    G.add_edge(stimChanID, idx + 1, weight=peak[idx] - baselineMean)
                    rowIdx = pltCount // nCol
                    colIdx = np.mod(pltCount, nCol)
                    if len(axes) > 0:
                        plt0 = self.splitViewTab_connection.matplotViews[
                            "Average traces"
                        ].figure.add_subplot(gs[rowIdx, colIdx], sharex=axes[0])
                    else:
                        plt0 = self.splitViewTab_connection.matplotViews[
                            "Average traces"
                        ].figure.add_subplot(gs[rowIdx, colIdx])
                    axes.append(plt0)
                    axes[pltCount].plot(
                        time, avg2[:, idx] * 1000.0, color=colorMap[idx]
                    )
                    axes[pltCount].plot(
                        time[stimStartSample + delaySamples + peak_idx[idx]],
                        peak[idx] * 1000,
                        marker=(8, 2, 0),
                        markersize=8,
                        markerfacecolor="k",
                    )
                    axes[pltCount].axvline(
                        time[stimStartSample + delaySamples],
                        0,
                        0.6,
                        color="gray",
                        linestyle="--",
                    )
                    axes[pltCount].axvline(
                        time[stimStartSample + delaySamples + responseSamples],
                        0,
                        0.6,
                        color="gray",
                        linestyle="--",
                    )
                    # axes[pltCount].axis('off')

                    # axes[pltCount].set_xlim(xlimRange)
                    #                    axes[j].legend("Ch {:d}  p={:.2f}".format(j+1, p), frameon=False, fancybox = False, numpoints = 6)
                    # axes[j].set_title("Ch {:d}  p={:.2f}".format(j+1, p), fontsize = 8)
                    axes[pltCount].set_title(
                        "Ch {:d} -> Ch {:d}".format(stimChanID, idx + 1), fontsize=10
                    )
                    axes[pltCount].spines["right"].set_visible(False)
                    axes[pltCount].spines["top"].set_visible(False)
                    axes[pltCount].spines["bottom"].set_visible(False)
                    axes[pltCount].yaxis.set_ticks_position("left")
                    axes[pltCount].tick_params(bottom=False)
                    axes[pltCount].tick_params(labelbottom=False)
                    pltCount = pltCount + 1
        rowIdx = pltCount // nCol
        plt0 = self.splitViewTab_connection.matplotViews[
            "Average traces"
        ].figure.add_subplot(gs[rowIdx, colIdx], sharex=axes[0])
        axes.append(plt0)
        axes[-1].plot(time, stimData, color="k")
        axes[-1].spines["right"].set_visible(False)
        axes[-1].spines["top"].set_visible(False)
        axes[-1].set_xlabel("Time (S)")
        edgeStyles = []
        for e in list(G.edges):
            if edges[e] < 0:
                edgeStyles.append("dashed")
            else:
                edgeStyles.append("solid")

        self.splitViewTab_connection.matplotViews["Average traces"].draw()
        self.splitViewTab_connection.matplotViews[
            "Average traces"
        ].figure.tight_layout()
        self.splitViewTab_connection.tables["Detected connections"].clear()
        df = pd.DataFrame(df)
        self.updateConnectionTable(
            df
        )  ## update connection table for significant connections
        self.visualizingConnectionsGraph(G, colorMap, edgeStyles)
        return G, colorMap, edgeStyles

    def getNodePositions(self):
        ## get node position from user defined ROIs
        ## otherwise return ''
        if self.roidata == []:
            return []
        else:
            chanID = {}
            for idx, x, y in self.roidata:  ## extract chan numeric index
                z0 = 0
                z = [int(t) for t in idx if np.char.isdigit(t)]  ## digit only
                for j, n in enumerate(z[::-1]):
                    z0 += n * 10**j  ## back to acutal decimal number
                chanID.update({z0: [x, y]})
            return chanID

    def visualizingConnectionsGraph(self, G, colorMap, s):

        self.splitViewTab_connection.matplotViews["Graph"].clf()
        self.splitViewTab_connection.topSplitter.setStretchFactor(0, 2)
        self.splitViewTab_connection.topSplitter.setStretchFactor(1, 1)
        gs = self.splitViewTab_connection.matplotViews["Graph"].figure.add_gridspec(
            1, 1
        )
        plt0 = self.splitViewTab_connection.matplotViews["Graph"].figure.add_subplot(
            gs[0, 0]
        )
        pos = self.getNodePositions()
        pos_default = nx.circular_layout(G)

        if pos == []:
            pos = pos_default  ## use default if no usable cordinates
        else:
            print(pos.keys())
            print(pos_default.keys())
            for k in pos.keys():
                if k not in list(pos_default.keys()):
                    G.add_node(k)
                pos_default.update(
                    {k: np.array(pos[k])}
                )  ## update postion from user defined cordinates
            pos = pos_default
        c = [colorMap[j - 1] for j in list(G.nodes)]
        nx.draw_networkx_nodes(G, pos, ax=plt0, node_size=160, node_color=c)
        # Ns.set_color(c)
        nx.draw_networkx_labels(G, pos, font_size=13, font_color="w", ax=plt0)
        arcs = nx.draw_networkx_edges(
            G,
            pos,
            ax=plt0,
            connectionstyle="arc3,rad=0.15",
            min_target_margin=15,
            arrowsize=25,
            width=1.5,
        )
        for i, arc in enumerate(arcs):  # change alpha values of arcs
            arc.set_linestyle(s[i])

        self.splitViewTab_connection.matplotViews["Graph"].draw()
        plt0.invert_yaxis()
        plt0.spines["right"].set_visible(False)
        plt0.spines["top"].set_visible(False)
        plt0.spines["bottom"].set_visible(False)
        plt0.spines["left"].set_visible(False)
        plt0.tick_params(
            axis="both",
            bottom=False,
            top=False,
            left=False,
            right=False,
            labelbottom=False,
            labeltop=False,
            labelleft=False,
            labelright=False,
        )
        plt0.axis("off")

    def getDatSeries(self, dat_file, seriesIdx0):
        bundle = HEKA.Bundle(dat_file)
        nSweep = len(
            bundle.pul[seriesIdx0[0]].children[seriesIdx0[1]].children
        )  ## number of sweeps for this series
        nTraces = len(bundle.pul[seriesIdx0[0]].children[seriesIdx0[1]].children[0])
        serieLabel = bundle.pul[seriesIdx0[0]].children[seriesIdx0[1]].Label
        # stimChanID = self.label2StimChanID(serieLabel)-1 ## alway pass around 0 based value!!
        stimChanID = self.searchForStimchan(serieLabel)  # , seriesIdx0)
        seriesIdx = list(seriesIdx0.copy())  ## get series level index
        seriesIdx.append(0)  ## sweep 0
        seriesIdx.append(0)  ## trace 0
        data0 = bundle.data[[seriesIdx[0], seriesIdx[1], 0, 0]]
        #        d = {'trace': [], 'time': [], 'sweep': [], 'data': []}
        #        DF0 = pd.DataFrame(d)
        #        DF0 = []
        data = np.zeros(
            (data0.shape[0], nTraces, nSweep)
        )  ## ndrarry, time X traces X sweep
        ## loops through all the sweeps
        for sweep in range(nSweep):
            seriesIdx[2] = sweep  ## change sweeep level index
            for t in range(nTraces):
                seriesIdx[3] = t  ## change sweeep level index
                da_ = bundle.data[seriesIdx]
                da_ = self.bandPass_signal(da_)
                data[:, t, sweep] = da_
                if t == 0:
                    trace = bundle.pul[seriesIdx[0]][seriesIdx[1]][seriesIdx[2]][
                        seriesIdx[3]
                    ]  ## get trace meta information
                    time = np.linspace(
                        trace.XStart,
                        trace.XStart + trace.XInterval * (len(da_) - 1),
                        len(da_),
                    )

        return time, data, bundle, stimChanID

    def replotForExport(self, time, data, index, mplWidget=None, SignificantOnly=True):
        time2, stim, stimInfo = self.currentPulseTree.bundle.stim(
            index
        )  ## assume only 1 channel is stimulating
        fs = self.parameters["fs"]
        stimStartSample = stimInfo[1]["start"]
        (
            delaySamples,
            responseSamples,
            avg2,
            peak,
            peak_idx,
            significantTest,
            baseline,
            meanCorrelation,
        ) = calculateConnectionTraces(time, fs, data, stimStartSample)
        nSignificant = np.sum(significantTest)
        nCol = nSignificant
        nRow = 2  # nSignificant //nCol
        if nRow == 0:
            nRow = 1
        gs = mplWidget.figure.add_gridspec(nRow, nCol)
        axes = []
        axes2 = []
        colorMap = [
            pg.glColor(pg.intColor(j, hues=data.shape[1])) for j in range(data.shape[1])
        ]
        # font = {'family': 'serif',
        # 'color':  'darkred',
        # 'weight': 'normal',
        # 'size': 16,
        # }
        xlimRange = [0.15, 0.35]
        stimChanID = self.parameters["Protocols"]["This serie"]["StimChanID"] - 1
        c = 0
        # print(data.shape[1])
        mplWidget.figure.clear()
        mplWidget.show()
        mplWidget.clf()
        for j in range(data.shape[1]):
            if significantTest[j] == 1 and stimChanID != j:
                c = c + 1
                plt0 = mplWidget.figure.add_subplot(gs[0, c - 1])
                axes.append(plt0)
                axes[c - 1].plot(time, avg2[:, j] * 1000, color=colorMap[j])
                axes[c - 1].plot(
                    time[stimStartSample + delaySamples + peak_idx[j]],
                    peak[j] * 1000,
                    marker=(8, 2, 0),
                    markersize=8,
                    markerfacecolor="k",
                )
                axes[c - 1].axvline(
                    time[stimStartSample + delaySamples],
                    0,
                    0.6,
                    color="gray",
                    linestyle="--",
                )
                axes[c - 1].axvline(
                    time[stimStartSample + delaySamples + responseSamples],
                    0,
                    0.6,
                    color="gray",
                    linestyle="--",
                )
                # axes[c-1].axis('off')
                axes[c - 1].set_xlim(xlimRange)
                axes[c - 1].spines["right"].set_visible(False)
                axes[c - 1].spines["top"].set_visible(False)
                axes[c - 1].spines["bottom"].set_visible(False)
                axes[c - 1].yaxis.set_ticks_position("left")
                axes[c - 1].yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
                axes[c - 1].tick_params(bottom=False)
                axes[c - 1].tick_params(labelbottom=False)
                #                    axes[j].legend("Ch {:d}  p={:.2f}".format(j+1, p), frameon=False, fancybox = False, numpoints = 6)
                # axes[j].set_title("Ch {:d}  p={:.2f}".format(j+1, p), fontsize = 8)
                axes[c - 1].set_title(
                    "Ch {:d} -> Ch {:d}".format(stimChanID + 1, j + 1), fontsize=10
                )

                plt1 = mplWidget.figure.add_subplot(gs[1, c - 1])
                axes2.append(plt1)
                axes2[c - 1].plot(
                    time, avg2[:, stimChanID] * 1000, color=colorMap[stimChanID]
                )
                axes2[c - 1].axvline(
                    time[stimStartSample + delaySamples],
                    0,
                    0.6,
                    color="gray",
                    linestyle="--",
                )
                axes2[c - 1].axvline(
                    time[stimStartSample + delaySamples + responseSamples],
                    0,
                    0.6,
                    color="gray",
                    linestyle="--",
                )
                # axes2[c-1].axis('off')
                axes2[c - 1].set_xlim(xlimRange)
                axes2[c - 1].spines["right"].set_visible(False)
                axes2[c - 1].spines["top"].set_visible(False)
                axes2[c - 1].spines["bottom"].set_visible(False)
                axes2[c - 1].yaxis.set_ticks_position("left")
                axes2[c - 1].yaxis.set_major_formatter(FormatStrFormatter("%.1f"))

        self.showSpikePanelsFlag = 1

    def replotAll(
        self,
        time,
        stimStartSample,
        responseSamples,
        avg2,
        significantTest,
        mplWidget=None,
        SignificantOnly=False,
    ):
        if mplWidget != None:  ## pure matplotlib
            print("")

    def replotEverything0(
        self, time, data, index, mplWidget=None, SignificantOnly=False
    ):
        if mplWidget != None:  ## pure matplotlib
            self.spikes_view.show()
            self.plot_splitter.setStretchFactor(0, 1)
            self.plot_splitter.setStretchFactor(1, 1)
            self.topPlot_splitter.setStretchFactor(1, 1)
            time2, stim, stimInfo = self.currentPulseTree.bundle.stim(
                index
            )  ## assume only 1 channel is stimulating
            # print(stimInfo)
            stimTiming = stimInfo[1]
            stimStartSample = stimTiming["start"]
            # print(stimStartSample)
            avg = np.squeeze(np.mean(data, axis=2))
            mplWidget.figure.clear()
            mplWidget.show()
            mplWidget.clf()
            self.plot_splitter.setStretchFactor(0, 1)
            self.plot_splitter.setStretchFactor(1, 1)
            nTrace = avg.shape[1]
            colorMap = [pg.glColor(pg.intColor(j, hues=nTrace)) for j in range(nTrace)]
            nCol = 1
            nRow = nTrace // nCol
            gs = mplWidget.figure.add_gridspec(nRow, nCol)
            axes = []
            hcutFreq = 1000
            responseTime = (
                0.040  ## 100 ms for calcuating baseline And to compare after baseline
            )
            responseOnsetDelay = 0.006
            responseSamples = np.int(self.parameters["fs"] * responseTime)
            delaySamples = np.int(self.parameters["fs"] * responseOnsetDelay)
            xlimRange = [0.15, 0.8]
            nStd = 4
            font = {
                "family": "serif",
                "color": "darkred",
                "weight": "normal",
                "size": 16,
            }
            df = {}
            df["Dat file"] = []
            df["Source Chan"] = []
            df["Target Chan"] = []
            df["Baseline(mV)"] = []
            df["Peak value(mV)"] = []
            df["Delta(mV)"] = []
            df["Delay(mS)"] = []
            df["Trial consistency"] = []
            for j in range(nTrace):
                # data_j = np.squeeze(data[:,j,:])  ## data from channel j: time X sweeps
                # baselineMean =
                data_ = np.squeeze(data[:, j, :])
                data_ = self.bandPass_signal(data_.transpose(), hcutFreq)
                # print(data_.shape)
                avg2 = np.mean(data_, axis=0)
                # avg2 = np.pad(avg[:,j], 1000, 'median')
                # avg2 = self.bandPass_signal(avg2, hcutFreq)
                # avg2 = avg2[1000:-1000]
                t0 = stimStartSample - delaySamples
                baselineData = np.mean(data_[:, :t0], axis=0)
                baselineMean = np.mean(baselineData)
                baselineStd = np.std(baselineData)
                t1 = stimStartSample + delaySamples
                t2 = stimStartSample + delaySamples + responseSamples
                responseData = np.squeeze(np.mean(data_[:, t1:t2], axis=0))
                peak_idx = np.argmax(
                    np.abs(responseData - baselineMean)
                )  ## of course it should relative to baseline!
                # print(t1+peak_idx)
                peak = responseData[peak_idx]  ## get peak during this window
                # print('Chan'+str(j+1))
                # print(baselineMean, baselineStd, peak)
                trialConsistency = np.corrcoef(data_[:, t1:t2], rowvar=True)
                selIdx = np.triu_indices(data_.shape[0], 1)
                meanCorrelation = np.mean(trialConsistency[selIdx])
                if (
                    np.abs(peak - baselineMean) >= nStd * baselineStd
                    and meanCorrelation >= 0.1
                ):
                    s = (8, 2, 0)
                    if (
                        j
                        != self.parameters["Protocols"]["This serie"]["StimChanID"] - 1
                    ):
                        existDat = [
                            self.splitViewTab_connection.tables[
                                "Detected connections"
                            ].item(0, kkk)
                            for kkk in range(
                                self.splitViewTab_connection.tables[
                                    "Detected connections"
                                ].rowCount()
                            )
                        ]  ## list of existing dat file names
                        if self.currentPulseTree.dat_file not in existDat:
                            df["Dat file"].append(self.currentPulseTree.dat_file)
                            df["Target Chan"].append(j + 1)
                            df["Source Chan"].append(
                                self.parameters["Protocols"]["This serie"]["StimChanID"]
                            )
                            df["Delay(mS)"].append(
                                (time[t1 + peak_idx] - time[stimStartSample]) * 1000
                            )
                            df["Peak value(mV)"].append(peak * 1000)
                            df["Delta(mV)"].append((peak - baselineMean) * 1000)
                            df["Baseline(mV)"].append(baselineMean * 1000)
                            df["Trial consistency"].append(meanCorrelation)
                else:
                    s = ""
                ## temp, p = stats.ranksums(responseData, baselineData)  ## Not very reliable!!

                # print('chan ' + str(j+1)+': p=', str(p))
                # print("Chan {:d}  p-value {:.2f}".format(j+1, p))
                plt0 = mplWidget.figure.add_subplot(gs[j // nCol, np.mod(j, nCol)])
                axes.append(plt0)
                axes[j].plot(time, avg2 * 1000.0, color=colorMap[j])
                axes[j].plot(
                    time[t1 + peak_idx],
                    peak * 1000.0,
                    marker=s,
                    markersize=8,
                    markerfacecolor="k",
                )
                axes[j].axvline(time[t1], 0, 0.6, color="gray", linestyle="--")
                axes[j].axvline(time[t2], 0, 0.6, color="gray", linestyle="--")
                # axes[j].axis('off')
                axes[j].set_xlim(xlimRange)
                axes[j].spines["right"].set_visible(False)
                axes[j].spines["top"].set_visible(False)
                axes[j].spines["bottom"].set_visible(False)
                axes[j].yaxis.set_ticks_position("left")
                axes[j].yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
                # axes[pltCount].tick_params(bottom=False)
                # axes[pltCount].tick_params(labelbottom=False)
                #                    axes[j].legend("Ch {:d}  p={:.2f}".format(j+1, p), frameon=False, fancybox = False, numpoints = 6)
                # axes[j].set_title("Ch {:d}  p={:.2f}".format(j+1, p), fontsize = 8)
                axes[j].set_title("Ch {:d}".format(j + 1), fontsize=8)
            mplWidget.figure.subplots_adjust(
                left=0.08, bottom=0.1, top=0.95, right=0.95, wspace=0.35, hspace=0.35
            )
            self.updateConnectionTable(
                pd.DataFrame(df)
            )  ## update connection table for significant connections
            self.showSpikePanelsFlag = 1

    def plot_averagedTraces_matplotlib(self, time, data, index, mplWidget=None):
        def onclick(event):
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Exporting significant traces only")
            msg.setWindowTitle("Exporting")
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            msg.buttonClicked.connect(msgbtn)
            msg.exec_()

        def msgbtn():
            print("test")

        self.replotEverything0(time, data, index, mplWidget)
        mplWidget.draw()

        # cid = mplWidget.figure.canvas.mpl_connect('button_press_event', onclick)

    def exportConnectionFig_clicked(self):
        selected = self.currentSelectedIndex
        print(selected)
        sel = selected[0]
        if len(sel.index) != 2:
            return
        self.checkSeriesType(sel)
        #        self.parameters['Protocols']['This serie']={'Type':'Firing pattern','StimChanID':stimChanIndex}
        if self.parameters["Protocols"]["This serie"]["Type"] in [
            "Connection",
            "Retina",
        ]:
            ## data should be ndarray, time X trace X sweep
            time, data, seriesIdx, traces = self.extractConnectionTraces(selected[0])
            #            df, seriesIdx = self.extractConnectionTraces(selected[0])
            #            avg = np.squeeze(np.mean(data, axis = 2))
            #            self.plotSeriesAverage(time, avg, seriesIdx)
            self.replotForExport(
                time, data, seriesIdx, self.spikes_view, SignificantOnly=True
            )
            #            self.spikes_view.refresh()
            # self.spikes_view.show()
            self.spikes_view.draw()

    def plotSeriesVariations(self, time, data, sel, traces):
        plt = []
        self.trace_view.clear()  ## refresh the layout
        nTraces = data.shape[1]
        avg = np.squeeze(np.mean(data, axis=2))

        self.plot_splitter.setStretchFactor(0, 1)
        self.plot_splitter.setStretchFactor(1, 1)
        self.topPlot_splitter.setStretchFactor(1, 1)
        self.spikes_view.show()
        ## This is the connection analysis!!!
        self.plot_averagedTraces_matplotlib(
            time, data, sel, self.spikes_view
        )  ## seperate matplotlib window
        #        data_std = np.squeeze(np.std(data, axis = 2))
        #        print(avg.shape)
        #        print(data_std.shape)
        for j in range(nTraces):

            myPen = pg.mkPen(
                color=pg.intColor(j, hues=nTraces)
            )  ## the pen to draw this
            if j == self.parameters["Protocols"]["This serie"]["StimChanID"] - 1:
                stimPen = myPen
            plt.append(
                self.trace_view.addPlot(col=0, colspan=4)
            )  ## add a subplot on the graphic layout

            self.trace_view.nextRow()  ## get ready for next subplot
            #            if j < nTraces-1:  ## disable xlable to avoid cluttering
            #                xlabelOn = False
            #            else:
            #                xlabelOn = True
            plt[j].showGrid(x=True, y=True)
            if j == nTraces - 1:
                plt[j].setLabels(
                    bottom=("Time", traces[0].XUnit), left=(traces[0].YUnit)
                )
            #            err = pg.ErrorBarItem(x = time, y = avg[:,j] , top = data_std[:,j], bottom = data_std[:,j], beam = 0.5)
            g = plt[j].plot(time, avg[:, j], pen=myPen, name=traces[j].Label)
            legend = pg.LegendItem(offset=(70, 20))
            legend.setParentItem(plt[j])
            legend.addItem(g, traces[j].Label)
            self.plotItemList.append(plt[j])
            plt[j].hideAxis("bottom")
            # plt[j].hideAxis('left')
        #            plt[j].addItem(err)

        if self.parameters["plotStim"]:
            self.trace_view.nextRow()  ## for stimulation
            plt_stim = self.trace_view.addPlot(col=0, colspan=4)
            plt_stim.setXRange(0.1, 0.45, padding=0)
            self.plotSingleStimTrace(
                plt_stim, sel, stimPen
            )  ## third arugument is for pen. if empty use default pen
            for p in plt:
                p.setXLink(plt_stim)  ## link x axis for all subplots
            #            plt[0].showAxis('bottom', False)
            title = (
                self.sliceName[:-4] + " " + "Series" + str(sel[1] + 1)
            )  # + ' sweep'+ str(sel[2]+1)
            plt_stim.setTitle(title)
            # plt[-1].setLabel('bottom', "", units='')
            self.trace_view.ci.layout.setRowStretchFactor(0, 3)
            self.trace_view.ci.layout.setRowStretchFactor(1, 1)
        else:
            for p in plt[:-1]:
                p.setXLink(plt[-1])  ## link x axis for all subplots
        if self.plotItemList:  ## update mouse mode if not empty
            self.setViewboxMouseMode()

    def plotSingleTrace3(
        self, plotHandle, time, data, myPen, traceLabel, xlabelOn=True, scaleBarOn=True
    ):
        ## with data alread
        ## plot a single trace
        plotHandle.showGrid(x=True, y=True)
        g = plotHandle.plot(time, data, pen=myPen, name=traceLabel)
        plotHandle.autoRange()
        if not xlabelOn:  ## no x labels
            plotHandle.setLabel("bottom", "", units="")
            #            ax=plotHandle.getAxis('bottom')    #This is the trick
            #            ax.setTicks("")
            plotHandle.showAxis("bottom", False)
        if scaleBarOn:
            scale = pg.ScaleBar(size=0.02, suffix="s")
            scale.setParentItem(plotHandle.getViewBox())
            scale.anchor((1, 1), (1, 1), offset=(-20, -20))
        legend = pg.LegendItem(offset=(70, 20))
        legend.setParentItem(plotHandle)
        legend.addItem(g, traceLabel)

    def extractConnectionTraces(self, sel):
        """Main function to get the series traces out"""
        ## get number of sweeps
        nSweep = sel.node.NumberSweeps
        nTraces = len(sel.node.children[0].children)
        seriesIdx = list(sel.index.copy())  ## get series level index
        seriesIdx.append(0)  ## sweep 0
        seriesIdx.append(0)  ## trace 0
        data0 = self.currentPulseTree.bundle.data[[seriesIdx[0], seriesIdx[1], 0, 0]]
        #        d = {'trace': [], 'time': [], 'sweep': [], 'data': []}
        #        DF0 = pd.DataFrame(d)
        #        DF0 = []
        data = np.zeros(
            (data0.shape[0], nTraces, nSweep)
        )  ## ndrarry, time X traces X sweep
        ## loops through all the sweeps
        traces = []
        for sweep in range(nSweep):
            seriesIdx[2] = sweep  ## change sweeep level index
            for t in range(nTraces):

                seriesIdx[3] = t  ## change trace level index
                trace = self.currentPulseTree.bundle.pul[seriesIdx[0]][seriesIdx[1]][
                    seriesIdx[2]
                ][
                    seriesIdx[3]
                ]  ## get trace meta information
                time, data_ = self.extractSingleTrace(trace, seriesIdx)
                if len(data_) != data0.shape[0]:
                    print(
                        "Sweep:",
                        sweep,
                        "Trace:",
                        t + 1,
                        "len of time and data: ",
                        len(time),
                        len(data_),
                    )
                    data_ = 0  # np.nan
                data[:, t, sweep] = data_
                traces.append(trace)
        return time, data, seriesIdx, traces

    def plotSeriesAverage(self, time, data, seriesIdx):
        plt = []
        nTraces = data.shape[1]
        for j in range(nTraces):

            myPen = pg.mkPen(
                color=pg.intColor(j, hues=nTraces)
            )  ## the pen to draw this
            plt.append(
                self.trace_view.addPlot(col=0, colspan=4)
            )  ## add a subplot on the graphic layout

            self.trace_view.nextRow()  ## get ready for next subplot
            #            if j < nTraces-1:  ## disable xlable to avoid cluttering
            #                xlabelOn = False
            #            else:
            #                xlabelOn = True
            self.plotSingleTrace2(
                plt[j], time, data[:, j], myPen, "trace" + str(j), True, False
            )

        if self.parameters["plotStim"]:
            self.trace_view.nextRow()  ## for stimulation
            plt_stim = self.trace_view.addPlot(col=0, colspan=4)
            self.plotSingleStimTrace(
                plt_stim, seriesIdx, []
            )  ## third arugument is for pen. if empty use default pen
            for p in plt:
                p.setXLink(plt_stim)  ## link x axis for all subplots
            #            plt[0].showAxis('bottom', False)

            plt[-1].setLabel("bottom", "", units="")
            title = (
                self.sliceName[:-4]
                + " "
                + "Series"
                + str(seriesIdx[1] + 1)
                + " sweep"
                + str(seriesIdx[2] + 1)
            )
            plt[-1].setTitle(title)
            self.trace_view.ci.layout.setRowStretchFactor(0, 3)
            self.trace_view.ci.layout.setRowStretchFactor(1, 1)
        else:
            for p in plt[:-1]:
                p.setXLink(plt[-1])  ## link x axis for all subplots

    def plotSingleTrace2(
        self, plotHandle, time, data, myPen, traceLabel, xlabelOn=True, scaleBarOn=True
    ):
        ## with data alread
        ## plot a single trace
        plotHandle.showGrid(x=True, y=True)
        g = plotHandle.plot(time, data, pen=myPen, name=traceLabel)
        plotHandle.autoRange()
        if not xlabelOn:  ## no x labels
            plotHandle.setLabel("bottom", "", units="")
            #            ax=plotHandle.getAxis('bottom')    #This is the trick
            #            ax.setTicks("")
            plotHandle.showAxis("bottom", False)
        if scaleBarOn:
            scale = pg.ScaleBar(size=0.02, suffix="s")
            scale.setParentItem(plotHandle.getViewBox())
            scale.anchor((1, 1), (1, 1), offset=(-20, -20))
        legend = pg.LegendItem(offset=(70, 20))
        legend.setParentItem(plotHandle)
        legend.addItem(g, traceLabel)

    def plotFiringPattern_clicked(self):
        """fine tune firing pattern plotting with user define parameters. To be finished"""
        selected = self.currentPulseTree.selectedItems()

    def update_newTree(self, bundle, root_item, index):
        """Recursively read tree information from the bundle's embedded .pul file
        and add items into the GUI tree to allow browsing.
        """
        root = bundle.pul
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
            item = datTree(node_type[3:], node_label)
        else:
            item = datTree(node_type, node_label)
        root_item.addChild(item)
        item.node = node
        item.index = index
        if len(index) == 2:
            self.seriesNode.append(node)
        if len(index) < 2:
            item.setExpanded(True)
        for i in range(len(node.children)):
            self.update_newTree(item, index + [i])

    def update_sortedTrees(self, dat_file):
        bundleClass = HekaBundleInfo(
            dat_file + ".dat"
        )  ## a bundle with some extra functions
        # bundle = HEKA.Bundle(dat_file +'.dat') #current bundle. raw class
        print(f"\n.dat file {dat_file}")
        print("-" * 20)
        if bundleClass.bundle.data == None:
            print("skip")
            return

        selectedSereis = []
        seriesTypes = []  ## could be different series type
        for g in range(bundleClass.countGroups()):  ## group level
            for seriesIdx in range(bundleClass.countSeries([g])):

                seriesLabel = bundleClass.getSeriesLabel(
                    [g, seriesIdx]
                )  # bundle.pul[g][seriesIdx]
                seriesType = self.querySessionProtocolByLabel(seriesLabel)
                if (
                    seriesType != "Spontaneous"
                ):  ## non-SP series should have minimal number of sweeps
                    if (
                        bundleClass.countSweeps([g, seriesIdx])
                        < self.parameters["CleanUp"]["minimalSweeps"]
                    ):
                        print(
                            f"   Series {seriesIdx:>2}: Label:",
                            seriesLabel,
                            f" has {bundleClass.countSweeps([g, seriesIdx])} sweep < minimal sweep specified. Excluding...",
                        )
                        continue
                selectedSereis.append([g, seriesIdx])
                seriesTypes.append(seriesType)
                print(
                    f'***Series {seriesIdx:>2}: Label: "{seriesLabel}", Type: {seriesType} added!'
                )
        self.loadQualifiedTree(dat_file, bundleClass, selectedSereis, seriesTypes)

    def clearAllTrees(self):
        self.reset()
        self.fpParTree_data_initilized = False
        for tree in self.selTreesDispatcher.keys():
            self.selTreesDispatcher[tree].clear()

    def loadQualifiedTree(self, dat_file0, bundle, seriesIdxList, seriesTypes):
        seriesTypeCounts = {"Firing pattern": 0, "Connection": 0, "Spontaneous": 0}
        for idx, (seriesIdx, seriesType) in enumerate(zip(seriesIdxList, seriesTypes)):
            self.selTreesDispatcher[
                seriesType
            ].bundle = bundle.bundle  ## current bundle
            self.selTreesDispatcher[
                seriesType
            ].current_dat = dat_file0  ## but we need this
            (
                self.selTreesDispatcher[seriesType].currentPath,
                self.selTreesDispatcher[seriesType].current_datFileName,
            ) = os.path.split(dat_file0)
            self.selTreesDispatcher[
                seriesType
            ].current_index = seriesIdx  ## and this to properly load our tree
            self.selTreesDispatcher[seriesType].update_treeSeries(
                self.selTreesDispatcher[seriesType].invisibleRootItem(), []
            )

            dat_file = dat_file0 + ".dat"
            self.selTreesDispatcher[seriesType].bundleGroup.update({dat_file: bundle})
            seriesTypeCounts[seriesType] += 1
            if dat_file in self.selTreesDispatcher[seriesType].dat_files.keys():
                self.selTreesDispatcher[seriesType].dat_files[dat_file].append(
                    seriesIdx
                )
            else:
                self.selTreesDispatcher[seriesType].dat_files.update(
                    {dat_file: [seriesIdx]}
                )

            if dat_file in self.selTreesDispatcher[seriesType].seriesTypes.keys():
                self.selTreesDispatcher[seriesType].seriesTypes[dat_file].append(
                    seriesType
                )
            else:
                self.selTreesDispatcher[seriesType].seriesTypes.update(
                    {dat_file: [seriesType]}
                )

        # max_key = max(seriesTypeCounts, key=seriesTypeCounts.get)
        seriesCounts = list(seriesTypeCounts.values())
        whichMax = np.argmax(seriesCounts)
        if isinstance(
            whichMax, np.int64
        ):  ## set the type with maximal series as current
            self.selectedFilesTabView.setCurrentIndex(whichMax)
        else:
            self.selectedFilesTabView.setCurrentIndex(0)

    def update_selectView(self, dat_file, Index):

        """
        Updating tree-selecting widgets. Book keeping added trees together with
        indexes. These can be used for batch-processing

        Parameters
        ----------
        dat_file : TYPE
            DESCRIPTION.
        Index : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.selTrees.bundle = HEKA.Bundle(dat_file)  ## current bundle
        self.selTrees.current_dat = dat_file
        self.selTrees.current_index = Index
        if dat_file not in self.selTrees.dat_files.keys():
            self.selTrees.dat_files[dat_file] = []
        if dat_file not in self.selTrees.stimChannelID.keys():
            self.selTrees.stimChannelID[dat_file] = []
        self.selTrees.dat_files[dat_file].extend(Index)
        if self.parameters["Protocols"]["This serie"]["Type"] in [
            "Firing pattern",
            "Connection",
        ]:
            self.selTrees.stimChannelID[dat_file].append(
                self.parameters["Protocols"]["This serie"]["StimChanID"] - 1
            )  ## ID for stimulation channel
        self.selTrees.update_tree_recursive(self.selTrees.invisibleRootItem(), [])

    def update_pul(self, dat_file, dat_index=None, ext=".dat"):
        """
        Makes call to update pul view.
        :param dat_file:
        :return:
        """
        if dat_file == "":
            return
        self.hideSpikePanels()
        self.currentPulseTree.clear()
        self.currentPulseTree.sweepCount = 0
        self.currentPulseTree.traceCount = 0
        self.currentPulseTree.dat_file = dat_file  ## data file
        if ext == ".dat":
            self.currentPulseTree.bundle = HEKA.Bundle(dat_file)
            self.currentPulseTree.update_tree_recursive(
                self.currentPulseTree.invisibleRootItem(), [], dat_index
            )
        elif ext == ".abf":
            self.currentPulseTree.abf = AxonIO(filename=dat_file)

            self.currentPulseTree.update_ABFtree(
                self.currentPulseTree.invisibleRootItem()
            )
            self.parameters["fs"] = self.currentPulseTree.abf.get_signal_sampling_rate()
            self.currentPulseTree.abf_stimInfo = self.extractStimData_ABF()
        # self.currentPulseTree.expandAll()
        self.gatherParameters(ext)
        self.updateStatusBar_fileName()

    def updateStatusBar_fileName(self):
        if self.currentPulseTree.dat_file != None:
            if hasattr(self, "QtFileNameLabel"):
                self.statusBar.removeWidget(self.QtFileNameLabel)
            self.QtFileNameLabel = QtWidgets.QLabel()
            self.QtFileNameLabel.setText(self.currentPulseTree.dat_file)
            self.statusBar.addWidget(self.QtFileNameLabel)
            self.root, self.sliceName = os.path.split(self.currentPulseTree.dat_file)

    def resettrace3DView(self):
        """For resetting 3D view of firing pattern traces"""
        #        self.trace3DView.reset()
        for x in self.trace3DView.items:
            x._setView(None)
        self.trace3DView.items = []
        self.trace3DView.opts["fov"] = 90
        self.trace3DView.opts["distance"] = 1.2
        self.trace3DView.opts["elevation"] = 30
        self.trace3DView.opts["azimuth"] = -100
        self.trace3DView.update()

    #        self.trace3DView.pan(0.2,-0.,-0.2)

    def update_trace_plot_selectionTree(self):
        """Show data associated with the selected tree node.
        For all nodes, the meta-data is updated in the bottom tree.
        For trace nodes, the data is plotted.
        """
        self.plt = []  ## plot handles
        self.plotItemList = []  ## collect all plot handles
        self.trace_view.clear()  ## refresh the layout
        #        self.trace3DView.clear()
        # self.trace3DView.hide()
        selected = self.selTrees.selectedItems()
        if not selected:
            return
        self.currentSelectedIndex = selected
        # update data tree
        sel = selected[0]
        serieIndex = sel.index
        if (
            self.EphyFeaturesObj
        ):  ## hide spike panels if a new series or new group is selected
            if (
                serieIndex[0] != self.EphyFeaturesObj.datIndex[0]
                or serieIndex[1] != self.EphyFeaturesObj.datIndex[1]
            ):
                #                self.EphyFeaturesObj = []  ## maybe we still want to get it back?
                self.hideSpikePanels()
        if len(sel.index) < 1:
            return

        self.vds = (
            1  ## to fine control whether to downsampling or not when visulaization
        )
        # plot selected dependent on what level is selected
        for kkk, sel in enumerate(selected):
            treeLevel = len(sel.index)  ## group 1; series, 2; sweep: 3; trace: 4

            if treeLevel < 1:
                return  ## do nothing at group level
            if treeLevel == 2:  ##  series level: all sweeps, all channels
                self.checkSeriesType(sel)  ## check protocol type for this selection
                self.plotSingleSeries(sel)
                self.plotSingleSeries3D(sel)
                self.parameter_Tab.clear()
                self.plotAveragedTraces(sel)
            if treeLevel == 3:  ##  sweep level: all channels at this sweep
                self.plotSingleSweep(sel)
                self.parameter_Tab.clear()
            if treeLevel == 4:  ##  trace level: all sweeps, all channels
                self.vds = 0  ## no downsampling
                if kkk == 0:
                    plt = self.trace_view.addPlot()
                self.plotSingleTrace(
                    plt, sel.node, sel.index, [], True, True
                )  ## third arugument is for pen. if empty use default pen
        if self.plotItemList:  ## update mouse mode if not empty
            self.setViewboxMouseMode()

    def update_trace_plot(self, selected=None):
        """Show data associated with the selected tree node.
        For all nodes, the meta-data is updated in the bottom tree.
        For trace nodes, the data is plotted.
        """
        self.plt = []  ## plot handles
        self.plotItemList = []  ## collect all plot handles
        ### BUILD A plot!
        # pdb.set_trace()
        self.trace_view.clear()  ## refresh the layout
        #        self.trace3DView.clear()
        # self.trace3DView.hide()
        self.updateUserParamters()
        selected = self.currentPulseTree.selectedItems()
        if not selected:
            return
        self.currentSelectedIndex = selected
        self.visulization_view.setCurrentIndex(0)
        # update data tree
        sel = selected[0]
        serieIndex = sel.index
        if (
            self.EphyFeaturesObj
        ):  ## hide spike panels if a new series or new group is selected
            if (
                serieIndex[0] != self.EphyFeaturesObj.datIndex[0]
                or serieIndex[1] != self.EphyFeaturesObj.datIndex[1]
            ):
                #                self.EphyFeaturesObj = []  ## maybe we still want to get it back?
                self.hideSpikePanels()
        if len(sel.index) < 1:
            return

        self.vds = (
            1  ## to fine control whether to downsampling or not when visulaization
        )
        # plot selected dependent on what level is selected
        Nsel = len(selected)
        if Nsel > 1:
            getSpikeFeature = False
            selected.sort(key=lambda x: x.index[2])  ## sorted
            multiSweeps = True
        else:
            getSpikeFeature = True  ## Default not to show spike features
            multiSweeps = False
        for kkk, sel in enumerate(selected):
            treeLevel = len(sel.index)  ## group 1; series, 2; sweep: 3; trace: 4
            self.checkSeriesType(sel)
            if treeLevel < 1:
                return  ## do nothing at group level
            if treeLevel == 2:  ##  series level: all sweeps, all channels
                # self.checkSeriesType(sel)  ## check protocol type for this selection
                if self.currentPulseTree.filetype == ".dat":
                    if (
                        self.parameters["Protocols"]["This serie"]["Type"]
                        != "Spontaneous"
                    ):
                        if self.showSpikePanelsFlag != 1:
                            self.showSpikePanels()
                        self.plotSingleSeries(sel)
                        # self.plotSingleSeries3D(sel)
                        self.parameter_Tab.clear()
                        self.plotAveragedTraces(sel)
                elif self.currentPulseTree.filetype == ".abf":
                    self.plotSingleSeries_ABF(sel)

            if treeLevel == 3:  ##  sweep level: all channels at this sweep
                if self.currentPulseTree.filetype == ".dat":
                    self.plotSingleSweep(sel, kkk, Nsel, getSpikeFeature, multiSweeps)
                    self.parameter_Tab.clear()
                else:
                    if kkk == 0:
                        plotHandle = self.trace_view.addPlot(row=0, col=0)
                        self.plotItemList.append(plotHandle)
                    segmentIdx = [sel.index[1], sel.index[2]]
                    trace = (
                        self.currentPulseTree.abfBlocks[segmentIdx[0]]
                        .segments[segmentIdx[1]]
                        .analogsignals[0]
                        .transpose()[0]
                    )
                    self.plotSingleTrace_ABF(
                        plotHandle,
                        segmentIdx,
                        trace,
                        None,
                        highCutOff=None,
                        analaysisSpike=True,
                    )

            if treeLevel == 4:  ##  trace level: all sweeps, all channels
                if self.currentPulseTree.filetype == ".dat":
                    self.vds = 0  ## no downsampling
                    if kkk == 0:
                        # self.checkSeriesType(self, sel)
                        plt = self.trace_view.addPlot(row=0, col=0)
                        # self.plotSingleTrace(plt, sel.node, sel.index, [], True, True, False) ## third arugument is for pen. if empty use default pen
                        if self.parameters["Protocols"]["This serie"]["Type"] in [
                            "Firing pattern",
                            "Connection",
                        ]:
                            self.trace_view.nextRow()  ## for stimulation
                            plt_stim = self.trace_view.addPlot(row=1, col=0)
                            # plt_stim.setXRange(0.1, 0.45, padding = 0)
                            stimPen = pg.mkPen(color="w")
                            self.trace_view.ci.layout.setRowStretchFactor(0, 3)
                            self.trace_view.ci.layout.setRowStretchFactor(1, 1)

                    if self.parameters["Protocols"]["This serie"]["Type"] in [
                        "Firing pattern",
                        "Connection",
                    ]:
                        traceHandle, *p = self.plotSingleTrace(
                            plt, sel.node, sel.index, [], False, True, False
                        )
                        self.plotSingleStimTrace(
                            plt_stim, sel.index, stimPen
                        )  ## third arugument is for pen. if empty use default pen
                        plt_stim.setXLink(traceHandle)

                    else:
                        self.plotSingleTrace(
                            plt, sel.node, sel.index, [], True, True, False
                        )
                elif self.currentPulseTree.filetype == ".abf":
                    print("Need to figure out if there are multichannels")

        if self.plotItemList:  ## update mouse mode if not empty
            self.setViewboxMouseMode()

    def plotSingleStimTrace(self, plt, index, myPen):
        """only at series level!"""
        time, stim, stimInfo = self.currentPulseTree.bundle.stim(
            index
        )  ## assume only 1 channel is stimulating
        #        tr = self.currentPulseTree.bundle.pgf

        if myPen == []:
            myPen = pg.mkPen(
                "w"
            )  ## default 8 channels. May set a global variable to set channle number
        plt.showGrid(x=True, y=True)
        plt.plot(time, stim, pen=myPen)
        plt.setLabels(bottom=("Time (s)"))
        plt.setLabels(left=("Current (pA)"))
        return stimInfo[1]["amplitude"]

    def label2StimChanID(self, label):
        label = label.replace("-", " ")
        for x in label.split():
            if (
                x in self.parameters["Protocols"]["Firing pattern"]
            ):  ## check for firing pattern protocl
                for y in label.split():
                    for yy in y:
                        if yy.isdigit():
                            stimChanIndex = int(yy)  ## extract channel index
                            return stimChanIndex

            if (
                x in self.parameters["Protocols"]["Connection"]
            ):  ## check for firing pattern protocl
                for y in label.split():
                    for yy in y:
                        if yy.isdigit():
                            stimChanIndex = int(yy)  ## extract channel index
                            # print(self.parameters['Protocols'])
                            return stimChanIndex

    def searchForStimchan(self, label):
        ## series label = sel.node.Label
        label = label.replace("-", " ")
        for x in label.split():
            if (
                x in self.parameters["Protocols"]["Firing pattern"]
            ):  ## check for firing pattern protocl
                self.parameters["Protocols"]["This serie"]["Type"] = "Firing pattern"
                for y in label.split():
                    if y.isdigit():
                        stimChanIndex = int(y)  ## extract channel index
                        self.parameters["Protocols"]["This serie"][
                            "StimChanID"
                        ] = stimChanIndex
                        # print(f'Stimuli channel found at {stimChanIndex}')
                        return stimChanIndex
                    else:
                        for yy in y:
                            if yy.isdigit():
                                stimChanIndex = int(yy)  ## extract channel index
                                self.parameters["Protocols"]["This serie"][
                                    "StimChanID"
                                ] = stimChanIndex
                                # print(f'Stimuli channel found at {stimChanIndex}')
                                return stimChanIndex

            if (
                x in self.parameters["Protocols"]["Connection"]
            ):  ## check for firing pattern protocl
                for y in label.split():
                    if y.isdigit():
                        stimChanIndex = int(y)  ## extract channel index
                        self.parameters["Protocols"]["This serie"] = {
                            "Type": "Connection",
                            "StimChanID": stimChanIndex,
                        }
                        # print(self.parameters['Protocols'])
                        # print(f'Stimuli channel found at {stimChanIndex}')
                        return stimChanIndex
                    else:
                        for yy in y:
                            if yy.isdigit():
                                stimChanIndex = int(yy)  ## extract channel index
                                self.parameters["Protocols"]["This serie"] = {
                                    "Type": "Connection",
                                    "StimChanID": stimChanIndex,
                                }
                                # print(self.parameters['Protocols'])
                                # print(f'Stimuli channel found at {stimChanIndex}')
                                return stimChanIndex
            if (
                x in self.parameters["Protocols"]["Spontaneous"]
            ):  ## check for firing pattern protocl
                self.parameters["Protocols"]["This serie"] = {
                    "Type": "Spontaneous",
                    "StimChanID": "",
                }
                # print('No timuli channel!')
                return ""

    def querySessionProtocolByLabel(self, label):
        ## series lable: label = sel.node.Label
        label = label.replace("-", " ")
        for x in label.split():
            if (
                x in self.parameters["Protocols"]["Firing pattern"]
            ):  ## check for firing pattern protocl
                return "Firing pattern"
            if (
                x in self.parameters["Protocols"]["Connection"]
            ):  ## check for firing pattern protocl
                return "Connection"
            if (
                x in self.parameters["Protocols"]["Spontaneous"]
            ):  ## check for firing pattern protocl
                return "Spontaneous"

    def checkSeriesType(self, sel):
        """Check the protocol type of the selected serie.
        and extract useful experiment related (stimulation channel etc.) informations from the label
        """
        #        print('checking series')
        if len(sel.index) < 2:  ## not a series
            return
        elif len(sel.index) == 3:
            sel = sel.parent()
        elif len(sel.index) == 4:
            sel = sel.parent()
            sel = sel.parent()  ## to reach series level!
        if self.currentPulseTree.filetype == ".dat":
            protocolType = self.querySessionProtocolByLabel(sel.node.Label)
            stimChanID = self.searchForStimchan(sel.node.Label)
            # print('stimchanID 5913', stimChanID)
            self.parameters["Protocols"]["This serie"] = {
                "Type": protocolType,
                "StimChanID": stimChanID,
            }

    def getCellTraceID(self, sel, chanID):
        """trace ID is not always consistent with channel ID.
        For example, one stim at Chan2, and record only chan2. Then traceID will be 1.
        The best way is to match the trace label with chan ID, for example for cell 2, trace label should be "Vmon-2"
        """
        sweep = sel.node.children[0]
        ## get all traces header information from this sweep
        traces = sweep.children
        traceLabel = "Vmon-" + str(
            chanID
        )  ## this is need to be changed if there's I-mon records
        target = None
        for j, t in enumerate(
            traces
        ):  ## loop through all traces in the first sweep to find the matching trace
            #            print(j, t.Label)
            if t.Label == traceLabel:
                target = j
                break
        if target == None:
            target = 0
        return target

    def getCellTraceID_sweep(self, sel, chanID):
        """trace ID is not always consistent with channel ID.
        For example, one stim at Chan2, and record only chan2. Then traceID will be 1.
        The best way is to match the trace label with chan ID, for example for cell 2, trace label should be "Vmon-2"
        """
        traces = sel.node.children
        ## get all traces header information from this sweep
        traceLabel = "Vmon-" + str(
            chanID
        )  ## this is need to be changed if there's I-mon records
        target = None
        for j, t in enumerate(
            traces
        ):  ## loop through all traces in the first sweep to find the matching trace
            #            print(j, t.Label)
            if t.Label == traceLabel:
                target = j
                break
        if target == None:
            target = 0
        return target

    def plotSingleSeries3D(self, sel):
        """quick overview of firing patterns of selected channels"""
        #        print('Plooting single series')
        if self.parameters["Protocols"]["This serie"]["Type"] == "Firing pattern":
            stimChanIndex = self.parameters["Protocols"]["This serie"]["StimChanID"]
            #            print(stimChanIndex)
            traceID = self.getCellTraceID(sel, stimChanIndex)
            ## get number of sweeps
            nSweep = len(sel.node.children)
            if sel.node.NumberSweeps != nSweep:
                print(
                    "Header-nSweep: %g, tree-nSweep: %g"
                    % (sel.node.NumberSweeps, nSweep)
                )
            seriesIdx = list(sel.index.copy())  ## get series level index
            seriesIdx.append(0)  ## sweep 0
            seriesIdx.append(traceID)
            ## loops through all the sweeps
            for sweep in range(0, nSweep):
                seriesIdx[2] = sweep  ## change sweeep level index
                # print(seriesIdx)
                #                myPen =  pg.mkPen(color=pg.intColor(sweep, hues = nSweep)) ## the pen to draw this
                try:
                    trace = self.currentPulseTree.bundle.pul[seriesIdx[0]][
                        seriesIdx[1]
                    ][seriesIdx[2]][
                        seriesIdx[3]
                    ]  ## get trace meta information
                except:
                    print("extra sweep index")
                # self.plotSingleTrace(plotHandle, trace, seriesIdx, myPen, True, True)
                data = self.currentPulseTree.bundle.data[seriesIdx]
                time = np.linspace(
                    trace.XStart,
                    trace.XStart + trace.XInterval * (len(data) - 1),
                    len(data),
                )

                pts = np.vstack(
                    [time, np.ones(len(time)) * sweep * 4e-2, data]
                ).transpose()
                plt = gl.GLLinePlotItem(
                    pos=pts, color=pg.glColor((sweep, nSweep * 1.1)), antialias=True
                )
                self.trace3DView.addItem(plt)

    def plotSingleSeries_ABF(self, sel):
        block = self.currentPulseTree.abfBlocks[
            sel.index[1]
        ]  ## choose the block selected
        plotHandle = self.trace_view.addPlot(
            col=0, colspan=4
        )  ## add a subplot on the graphic layout
        self.plotItemList.append(plotHandle)
        ## loops through all the sweeps.
        nSweep = len(block.segments)
        if nSweep > 3:  ## no need to plot all of them:
            sweepIdx = [0, int(nSweep / 2), nSweep - 1]
        else:
            sweepIdx = range(nSweep)
        for idx, sweepIdx_ in enumerate(sweepIdx):  # enumerate(block.segments):
            myPen = pg.mkPen(color=pg.intColor(idx, hues=3))  ## the pen to draw this
            trace = (
                block.segments[sweepIdx_].analogsignals[0].transpose()[0]
            )  ## numpy array for current sweep  ## get trace meta information
            segmentIdx = [sel.index[1], sweepIdx_]
            self.plotSingleTrace_ABF(
                plotHandle,
                segmentIdx,
                trace,
                myPen,
                title=False,
                scaleBar=False,
                analaysisSpike=False,
                plotStim=False,
            )
        scale = pg.ScaleBar(size=0.02, suffix="s")
        scale.setParentItem(plotHandle.getViewBox())
        scale.anchor((1, 1), (1, 1), offset=(-20, -20))
        plotHandle.setTitle(
            self.currentPulseTree.dat_file[:-4].split("\\")[-1]
            + " block"
            + str(sel.index[1] + 1)
        )

        if self.currentPulseTree.abf_stimOn:
            self.trace_view.nextRow()  ## for stimulation
            plt_stim = self.trace_view.addPlot(col=0, colspan=4)
            for idx, sweep in enumerate(sweepIdx):
                myPen = pg.mkPen(
                    color=pg.intColor(idx, hues=3)
                )  ## the pen to draw this
                stim = self.currentPulseTree.abf_stimData[:, sweep]
                plt_stim.plot(self.currentPulseTree.abf_stimTime, stim, pen=myPen)
            plt_stim.showGrid(x=True, y=True)

            plt_stim.setLabels(bottom=("Time (s)"))
            plt_stim.setLabels(
                left=(f"Current ({self.currentPulseTree.abf_stimUnit })")
            )
            plt_stim.setXLink(plotHandle)  ## link x axis for all subplots
            self.trace_view.ci.layout.setRowStretchFactor(0, 3)
            self.trace_view.ci.layout.setRowStretchFactor(1, 1)
            plotHandle.setLabel("bottom", "", units="")

    def plotSingleTrace_ABF(
        self,
        plotHandle,
        segmentIdx,
        trace,
        myPen,
        highCutOff=None,
        title=True,
        scaleBar=True,
        analaysisSpike=True,
        plotStim=True,
    ):
        plotHandle.showGrid(x=True, y=True)
        plotHandle.setLabels(
            bottom=("Time", self.currentPulseTree.abf.xUnits),
            left=(self.currentPulseTree.abf.yUnits),
        )

        data = self.bandPass_signal(trace, highCutOff)
        time = np.arange(len(data)) / self.parameters["fs"]
        if myPen == None:
            myPen = pg.mkPen("r")
        if self.OptionAction5.isChecked():
            print("remove stimuli artifact", segmentIdx)
            (
                t,
                v,
                dvdt,
                peaks,
                waveTime,
                waves,
                waves_dvdt,
                spike_df,
                dur,
                dv_cutoff,
            ) = self.getSingleTraceFeature_ABF(trace, time, segmentIdx)
            for p in peaks:
                data[p - 20 : p + 20] = np.mean(data[:100])
        g = plotHandle.plot(time, data, pen=myPen)
        plotHandle.autoRange()
        if scaleBar:
            scale = pg.ScaleBar(size=0.02, suffix="s")
            scale.setParentItem(plotHandle.getViewBox())
            scale.anchor((1, 1), (1, 1), offset=(-20, -20))
        if title:
            title_ = (
                self.currentPulseTree.dat_file[:-4].split("\\")[-1]
                + " "
                + "Block"
                + str(segmentIdx[0] + 1)
                + " Sweep"
                + str(segmentIdx[1] + 1)
            )
            plotHandle.setTitle(title_)

        if self.currentPulseTree.abf_stimOn:
            if analaysisSpike and (not self.OptionAction5.isChecked()):
                self.trace_view2.show()
                self.plot_alignedSpikes_SingleTrace(
                    segmentIdx,
                    mplWidget=self.spikes_view,
                    plotWidget2=self.trace_view2,
                    isABF=True,
                    trace=trace,
                    time=time,
                )
            if plotStim:
                self.trace_view.nextRow()  ## for stimulation
                plt_stim = self.trace_view.addPlot(col=0, colspan=4)
                # myPen = pg.mkPen('w')  ## the pen to draw this
                stim = self.currentPulseTree.abf_stimData[:, segmentIdx[1]]
                plt_stim.plot(self.currentPulseTree.abf_stimTime, stim, pen="w")
                plt_stim.showGrid(x=True, y=True)

                plt_stim.setLabels(bottom=("Time (s)"))
                plt_stim.setLabels(
                    left=(f"Current ({self.currentPulseTree.abf_stimUnit })")
                )
                plt_stim.setXLink(plotHandle)  ## link x axis for all subplots
                self.trace_view.ci.layout.setRowStretchFactor(0, 3)
                self.trace_view.ci.layout.setRowStretchFactor(1, 1)
                plotHandle.setLabel("bottom", "", units="")
        return data, time

    def plotSingleSeries(self, sel):
        """quick overview of firing patterns of selected channels"""
        #        print('Plooting single series')
        if self.parameters["Protocols"]["This serie"]["Type"] == "Firing pattern":
            stimChanIndex = self.parameters["Protocols"]["This serie"]["StimChanID"]
            traceID = self.getCellTraceID(sel, stimChanIndex)
            ## get number of sweeps
            nSweep = len(sel.node.children)
            if sel.node.NumberSweeps != nSweep:
                print(
                    "Header-nSweep: %g, tree-nSweep: %g"
                    % (sel.node.NumberSweeps, nSweep)
                )
            seriesIdx = list(sel.index.copy())  ## get series level index
            seriesIdx.append(0)  ## sweep 0
            seriesIdx.append(traceID)
            plotHandle = self.trace_view.addPlot(
                col=0, colspan=4
            )  ## add a subplot on the graphic layout
            ## loops through all the sweeps.
            if nSweep > 6:  ## no need to plot all of them:
                sweepIdx = [0, int(nSweep / 2), nSweep - 1]
            else:
                sweepIdx = range(nSweep)

            for idx, sweep in enumerate(sweepIdx):
                seriesIdx[2] = sweep  ## change sweeep level index
                # print(seriesIdx)
                myPen = pg.mkPen(
                    color=pg.intColor(idx, hues=3)
                )  ## the pen to draw this
                try:
                    trace = self.currentPulseTree.bundle.pul[seriesIdx[0]][
                        seriesIdx[1]
                    ][seriesIdx[2]][
                        seriesIdx[3]
                    ]  ## get trace meta information
                except:
                    print("extra sweep index")
                self.plotSingleTrace(plotHandle, trace, seriesIdx, myPen, True, True)
            if self.parameters["plotStim"] and self.parameters["Protocols"][
                "This serie"
            ]["Type"] in ["Firing pattern", "Connection"]:
                self.trace_view.nextRow()  ## for stimulation
                plt_stim = self.trace_view.addPlot(col=0, colspan=4)
                for idx, sweep in enumerate(sweepIdx):
                    seriesIdx[2] = sweep  ## change sweeep level index
                    myPen = pg.mkPen(
                        color=pg.intColor(idx, hues=3)
                    )  ## the pen to draw this
                    self.plotSingleStimTrace(
                        plt_stim, seriesIdx, myPen
                    )  ## third arugument is for pen. if empty use default pen
                plt_stim.setXLink(plotHandle)  ## link x axis for all subplots
                self.trace_view.ci.layout.setRowStretchFactor(0, 3)
                self.trace_view.ci.layout.setRowStretchFactor(1, 1)
                plotHandle.setLabel("bottom", "", units="")
            title_sufix = self.sliceName[:-4] + " " + "Series" + str(seriesIdx[1] + 1)
            sweepNames = ", ".join([str(s + 1) for s in sweepIdx])
            plotHandle.setTitle(title_sufix + " sweep:" + sweepNames)

    #            scale = pg.ScaleBar(size=0.1, suffix='s')
    #            scale.setParentItem(plt.getViewBox())
    #            scale.anchor((1, 1), (1, 1), offset=(-20, -20))

    def plotSingleSweep(
        self, sel, kkk=None, Nsel=None, getSpikeFeature=False, multiSweep=False
    ):
        ## series level object
        traces = sel.node.children  ## get all traces header information from this sweep
        plt = []
        self.checkSeriesType(sel)
        # print(self.parameters['Protocols'])
        #        print(self.parameters['Protocols']['This serie']['Type'])
        if self.parameters["Protocols"]["This serie"]["Type"] == "Firing pattern":
            stimChanIndex = self.parameters["Protocols"]["This serie"]["StimChanID"]
            traceID = self.getCellTraceID_sweep(sel, stimChanIndex)
            seriesIdx = list(sel.index.copy())  ## get series level index
            seriesIdx.append(traceID)
            #            print(seriesIdx)
            plt0 = self.trace_view.getItem(row=0, col=0)
            if plt0 == None:
                plt0 = self.trace_view.addPlot(row=0, col=0)
            # plt.append(self.trace_view.addPlot(row= 0, col = 0))  ## add a subplot on the graphic layout
            # print(seriesIdx)
            if Nsel == None:
                myPen = []
            else:
                myPen = pg.mkPen(
                    color=pg.intColor(Nsel - kkk - 1, hues=Nsel)
                )  ## the pen to draw this
                if Nsel < 2:
                    title = ""
                else:
                    title = " multiple sweeps"
            try:
                trace = self.currentPulseTree.bundle.pul[seriesIdx[0]][seriesIdx[1]][
                    seriesIdx[2]
                ][
                    seriesIdx[3]
                ]  ## get trace meta information
            except:
                print("extra sweep index")
            self.plotSingleTrace(
                plt0, trace, seriesIdx, myPen, True, True, False, 0, None, title
            )
            plt.append(plt0)
            if getSpikeFeature:
                self.plot_splitter.setStretchFactor(0, 1)
                self.plot_splitter.setStretchFactor(1, 1)
                self.trace_view2.show()
                self.plot_alignedSpikes_SingleTrace(
                    seriesIdx, mplWidget=self.spikes_view, plotWidget2=self.trace_view2
                )
        else:
            for j, t in enumerate(traces):
                seriesIdx = list(sel.index.copy())  ## get series level index
                seriesIdx.append(
                    j
                )  ## add trace level index. This is the final index to get trace data!
                myPen = pg.mkPen(
                    color=pg.intColor(j, hues=len(traces))
                )  ## the pen to draw this
                plt.append(
                    self.trace_view.addPlot(col=0, colspan=4)
                )  ## add a subplot on the graphic layout
                if len(traces) > 1:
                    self.trace_view.nextRow()  ## get ready for next subplot
                if j < len(traces) - 1:  ## disable xlable to avoid cluttering
                    xlabelOn = False
                else:
                    if self.parameters["plotStim"] and self.parameters["Protocols"][
                        "This serie"
                    ]["Type"] in ["Firing pattern", "Connection"]:
                        xlabelOn = (
                            False  ## no need to plot this axis if plotStim is on!
                        )
                    else:
                        xlabelOn = True
                self.plotSingleTrace(plt[j], t, seriesIdx, myPen, xlabelOn, xlabelOn)

        if self.parameters["plotStim"] and self.parameters["Protocols"]["This serie"][
            "Type"
        ] in ["Firing pattern", "Connection"]:
            plt_stim = self.trace_view.getItem(row=1, col=0)
            if plt_stim == None:
                singleStim = True
                self.trace_view.nextRow()  ## for stimulation
                plt_stim = self.trace_view.addPlot(row=1, col=0)
            else:
                singleStim = False
            if myPen == []:
                myPen = pg.mkPen(
                    "w"
                )  ## default 8 channels. May set a global variable to set channle number
            else:
                myPen = pg.mkPen(
                    color=pg.intColor(Nsel - kkk - 1, hues=Nsel)
                )  ## the pen to draw this
            stimAmp = self.plotSingleStimTrace(
                plt_stim, sel.index, myPen
            )  ## third arugument is for pen. if empty use default pen
            if singleStim:
                title = (
                    self.sliceName[:-4]
                    + " "
                    + "Series"
                    + str(seriesIdx[1] + 1)
                    + " Sweep"
                    + str(seriesIdx[2] + 1)
                    + " Trace"
                    + str(seriesIdx[3] + 1)
                    + f" {stimAmp:.0f} pA"
                )

            else:
                title = (
                    self.sliceName[:-4]
                    + " "
                    + "Series"
                    + str(seriesIdx[1] + 1)
                    + " Trace"
                    + str(seriesIdx[3] + 1)
                    + " multiple sweeps"
                )

            for p in plt:
                p.setXLink(plt_stim)  ## link x axis for all subplots
            plt_stim.setTitle(title)
            self.trace_view.ci.layout.setRowStretchFactor(0, 3)
            self.trace_view.ci.layout.setRowStretchFactor(1, 1)
        else:
            for p in plt[:-1]:
                p.setXLink(plt[-1])  ## link x axis for all subplots

    def getEphySpikePars(self):
        pv = self.splitViewTab_FP.getParTreePars("Spike detection")
        # This is just to get spikes peaks. Need to use higher cutoff for dvdt-v phase plot
        HF = pv["Spike detection parameters"][1]["High frequency cutoff"][
            0
        ]  # this may different than general filter options
        dv_cutoff = pv["Spike detection parameters"][1][
            "dv/dt cut off (min=1;max=100)"
        ][0]
        max_interval = pv["Spike detection parameters"][1][
            "max_interval (min=0.001;max=0.02)"
        ][0]
        min_height = pv["Spike detection parameters"][1]["peak height (min=2;max=10)"][
            0
        ]
        min_peak = pv["Spike detection parameters"][1][
            "peak voltage (min=-30;max=150)"
        ][0]
        thresh_frac = pv["Spike detection parameters"][1][
            "thresh_frac (min=0.02;max=0.08)"
        ][0]
        baseline_interval = pv["Spike detection parameters"][1][
            "baseline_interval (min=0.05;max=0.15)"
        ][0]
        baseline_detect_thresh = pv["Spike detection parameters"][1][
            "baseline_detect_thresh (min=-30;max=150)"
        ][0]

        return (
            dv_cutoff,
            min_height,
            min_peak,
            thresh_frac,
            baseline_interval,
            baseline_detect_thresh,
            max_interval,
            HF,
        )

    def getSingleTraceFeature_ABF(self, v, t, seriesIdx):
        self.updateUserParamters()
        v = np.array(v)
        stimInfo = self.currentPulseTree.abf_stimInfo
        sweepIdx = seriesIdx[1]
        start = (
            stimInfo[sweepIdx][1]["start"] * stimInfo[sweepIdx][1]["sampleInteval"]
        )  ## stimuli start
        end = (
            stimInfo[sweepIdx][1]["end"] * stimInfo[sweepIdx][1]["sampleInteval"]
        )  ## stimuli end
        print("start, end", start, end)
        vhold = stimInfo[sweepIdx][1]["Vholding"]
        curr = stimInfo[sweepIdx][1]["amplitude"]
        sampleRate = 1 / stimInfo[sweepIdx][1]["sampleInteval"]
        ## this is to get smoothy signals, but not over-filtered
        if self.parameters["HF"] < sampleRate / 2:
            filterHighCutFreq_signal = self.parameters["HF"]
        else:
            filterHighCutFreq_signal = np.floor((sampleRate - 100.0)) / 2  ## KHz
        (
            dv_cutoff1,
            min_height1,
            min_peak1,
            thresh_frac1,
            baseline_interval1,
            baseline_detect_thresh1,
            max_interval1,
            filterHighCutFreq_spike,
        ) = self.getEphySpikePars()

        ## this is to get spike peak.
        if filterHighCutFreq_spike >= sampleRate / 2:
            filterHighCutFreq_spike = np.floor((sampleRate - 100.0)) / 2  ## KHz

        stimSeries = self.currentPulseTree.abf_stimData[:, seriesIdx[1]] - vhold
        EphysObject = efex.EphysSweepFeatureExtractor(
            t=t,
            v=v,
            i=stimSeries,
            start=start,
            end=end,
            filter=filterHighCutFreq_spike / 1000,
            dv_cutoff=dv_cutoff1,
            max_interval=max_interval1,
            min_height=min_height1,
            min_peak=min_peak1,
            thresh_frac=thresh_frac1,
            baseline_interval=baseline_interval1,
            baseline_detect_thresh=baseline_detect_thresh1,
        )
        EphysObject.process_spikes()
        spike_df, sweep_df = extract_sweep_feature(t, v, curr, start, end, EphysObject)
        v_filtered = ephys_ft.calculate_v_filter(v, t, filterHighCutFreq_signal / 1000)
        dvdt = ephys_ft.calculate_dvdt(v, t, filterHighCutFreq_signal / 1000)
        dvdt1 = np.insert(dvdt, 0, 0)
        if EphysObject._spikes_df.size:
            peaks = spike_df["peak_index"].to_list()  ## list of peak index
            try:
                maxSPwidth = int(spike_df["width"].max(skipna=True) * sampleRate * 2.5)
            except:
                maxSPwidth = int(
                    (spike_df["peak_index"] - spike_df["threshold_index"]).max() * 5
                )
                if not self.OptionAction5.isChecked():  ## that's OK for non-spike type
                    print("exception happens for spike width detection")
            try:
                before_ = int(
                    (spike_df["peak_index"] - spike_df["threshold_index"]).max()
                ) + int(1 * sampleRate / 1000.0)
            except:
                maxSPwidth = int(maxSPwidth // 2)
            try:
                after_ = int((spike_df["trough_index"] - spike_df["peak_index"]).max())
                if after_ > 3 * before_:
                    after_ = 3 * before_
            except:
                after_ = int(maxSPwidth - maxSPwidth // 20)

            waves = np.zeros((after_ + before_, len(peaks)))
            waves_dvdt = np.zeros((after_ + before_, len(peaks)))
            waveTime = np.arange(-before_, after_) / sampleRate * 1e3
            for j, p in enumerate(peaks):
                start_index = p - before_
                end_index = p + after_
                # print(start_index, end_index, sweepCount)
                waves[:, j] = v_filtered[start_index:end_index]
                waves_dvdt[:, j] = dvdt1[start_index:end_index]
        else:
            peaks = []
            waves = []
            waves_dvdt = []
            waveTime = []
        return (
            t,
            v_filtered,
            dvdt1,
            peaks,
            waveTime,
            waves,
            waves_dvdt,
            spike_df,
            end - start,
            dv_cutoff1,
        )

    def getSingleTraceFeature(self, seriesIdx):
        self.updateUserParamters()
        trace = self.currentPulseTree.bundle.pul[seriesIdx[0]][seriesIdx[1]][
            seriesIdx[2]
        ][
            seriesIdx[3]
        ]  ## get trace meta information
        v = self.currentPulseTree.bundle.data[seriesIdx]
        v = v * 1000.0
        t = np.linspace(
            trace.XStart, trace.XStart + trace.XInterval * (len(v) - 1), len(v)
        )

        stimTime, stim, stimInfo = self.currentPulseTree.bundle.stim(
            seriesIdx
        )  ## assume only 1 channel is stimulating
        start = stimInfo[1]["start"] * stimInfo[1]["sampleInteval"]  ## stimuli start
        end = stimInfo[1]["end"] * stimInfo[1]["sampleInteval"]  ## stimuli end
        vhold = stimInfo[1]["Vholding"]
        curr = stimInfo[1]["amplitude"]  ## absolute
        stim = np.array(stim) - vhold  ## convert it to relative to holding current
        sampleRate = 1 / stimInfo[1]["sampleInteval"]

        ## this is to get smoothy signals, but not over-filtered
        if self.parameters["HF"] < sampleRate / 2:
            filterHighCutFreq_signal = self.parameters["HF"]  ## 10 KHz
        else:
            filterHighCutFreq_signal = np.floor((sampleRate - 100.0)) / 2  ## KHz
        (
            dv_cutoff1,
            min_height1,
            min_peak1,
            thresh_frac1,
            baseline_interval1,
            baseline_detect_thresh1,
            max_interval1,
            filterHighCutFreq_spike,
        ) = self.getEphySpikePars()

        ## this is to get spike peak.
        if filterHighCutFreq_spike >= sampleRate / 2:
            filterHighCutFreq_signal = np.floor((sampleRate - 100.0)) / 2  ## KHz

        EphysObject = efex.EphysSweepFeatureExtractor(
            t=t,
            v=v,
            i=stim,
            start=start,
            end=end,
            filter=filterHighCutFreq_spike / 1000,
            dv_cutoff=dv_cutoff1,
            max_interval=max_interval1,
            min_height=min_height1,
            min_peak=min_peak1,
            thresh_frac=thresh_frac1,
            baseline_interval=baseline_interval1,
            baseline_detect_thresh=baseline_detect_thresh1,
        )
        EphysObject.process_spikes()
        spike_df, sweep_df = extract_sweep_feature(t, v, curr, start, end, EphysObject)
        v_filtered = ephys_ft.calculate_v_filter(v, t, filterHighCutFreq_signal / 1000)
        dvdt = ephys_ft.calculate_dvdt(v, t, filterHighCutFreq_signal / 1000)
        dvdt1 = np.insert(dvdt, 0, 0)
        if EphysObject._spikes_df.size:
            peaks = spike_df["peak_index"].to_list()  ## lis tof peak index
            try:
                maxSPwidth = int(spike_df["width"].max(skipna=True) * sampleRate * 2.5)
            except:
                maxSPwidth = int(
                    (spike_df["peak_index"] - spike_df["threshold_index"]).max() * 5
                )
                print("exception happens")
            try:
                before_ = (
                    int((spike_df["peak_index"] - spike_df["threshold_index"]).max())
                    + 100
                )
            except:
                maxSPwidth = int(maxSPwidth // 2)
            try:
                after_ = int((spike_df["trough_index"] - spike_df["peak_index"]).max())
                if after_ > 3 * before_:
                    after_ = 3 * before_
            except:
                after_ = int(maxSPwidth - maxSPwidth // 20)

            waves = np.zeros((after_ + before_, len(peaks)))
            waves_dvdt = np.zeros((after_ + before_, len(peaks)))
            waveTime = np.arange(-before_, after_) / sampleRate * 1e3
            for j, p in enumerate(peaks):
                start_index = p - before_
                end_index = p + after_
                # print(start_index, end_index, sweepCount)
                waves[:, j] = v_filtered[start_index:end_index]
                waves_dvdt[:, j] = dvdt1[start_index:end_index]
        else:
            peaks = []
            waves = []
            waves_dvdt = []
            waveTime = []
        return (
            t,
            v_filtered,
            dvdt1,
            peaks,
            waveTime,
            waves,
            waves_dvdt,
            spike_df,
            end - start,
            dv_cutoff1,
        )

    def plot_dvdt_v_phasePlot(self, v_spike, dvdt_spike, sweepNumber, plotWidget):
        plotWidget.clear()  ## refresh the layout
        plt = plotWidget.addPlot()
        plt.showGrid(x=True, y=True)
        if not isinstance(v_spike, list):
            for j in range(v_spike.shape[1]):
                myPen = pg.mkPen(
                    color=pg.intColor(j, hues=v_spike.shape[1] + 1)
                )  ## the pen to draw this
                # dvdt = ft.calculate_dvdt(data[:,j], time, self.filter)
                # dvdt = np.diff(data[:,j])*self.sampleRate/1000
                if self.dv2dt2:
                    dv2dt2 = np.diff(dvdt_spike[:, j])
                    plt.plot(
                        dvdt_spike[1:, j],
                        dv2dt2,
                        pen=myPen,
                        name="Sweep" + str(sweepNumber),
                    )
                    plt.setLabel("left", "dv^2/dt^2", units="mV/mS^2")
                    plt.setLabel("bottom", "dv/dt", units="mV/mS")
                else:
                    plt.plot(
                        v_spike[:, j],
                        dvdt_spike[:, j],
                        pen=myPen,
                        name="Sweep" + str(sweepNumber),
                    )
                    plt.setLabel("left", "dv/dt", units="mV/mS")
                    plt.setLabel("bottom", "Voltage", units="mV")
            plt.setTitle(
                "Sweep "
                + str(sweepNumber + 1)
                + " Spike count:"
                + str(v_spike.shape[1])
            )
        return plt

    def plot_alignedSpikes_SingleTrace(
        self,
        traceIdx,
        axis=None,
        mplWidget=None,
        plotWidget2=None,
        isABF=False,
        trace=None,
        time=None,
    ):
        if isABF:
            (
                t,
                v,
                dvdt,
                peaks,
                waveTime,
                waves,
                waves_dvdt,
                spike_df,
                dur,
                dv_cutoff,
            ) = self.getSingleTraceFeature_ABF(trace, time, traceIdx)
            title = (
                self.currentPulseTree.dat_file[:-4].split("\\")[-1]
                + " "
                + "Block"
                + str(traceIdx[0] + 1)
                + " Sweep"
                + str(traceIdx[1] + 1)
            )
            sweepNumber = traceIdx[1]
        else:
            (
                t,
                v,
                dvdt,
                peaks,
                waveTime,
                waves,
                waves_dvdt,
                spike_df,
                dur,
                dv_cutoff,
            ) = self.getSingleTraceFeature(traceIdx)
            title = (
                self.sliceName[:-4]
                + " "
                + "Series"
                + str(traceIdx[1] + 1)
                + " sweep"
                + str(traceIdx[2] + 1)
            )
            sweepNumber = traceIdx[2]

        print(title)
        if mplWidget == None:  ## pure matplotlib
            fig = MATPLT.figure(figsize=(10, 10))
            fig.canvas.set_window_title(title)
            MATPLT.subplot(2, 1, 1)
            MATPLT.plot(t, v)
            MATPLT.plot(t, dvdt)
            for p in peaks:
                MATPLT.plot(t[p], v[p], "xr")
            MATPLT.xlabel("Time (s)")
            MATPLT.ylabel("Voltage (mV)")
            MATPLT.subplot(2, 1, 2)
            for j in range(len(peaks)):
                MATPLT.plot(waveTime, waves[:, j])
            MATPLT.xlabel("Time (ms)")
            MATPLT.ylabel("Voltage (mV)")
        else:  ## pyqtgraph MatplotlibWidget.MatplotlibWidget objects defined in PatchViewer 3
            mplWidget.figure.clear()
            mplWidget.show()
            mplWidget.clf()
            gs = mplWidget.figure.add_gridspec(2, 2)
            axes = []
            plt0 = mplWidget.figure.add_subplot(gs[0, 0])
            axes.append(plt0)
            # a secondear axis for dvdt
            ax_y2 = axes[0].twinx()
            ax_y2.plot(t, dvdt, linewidth=1.0, alpha=0.7)

            ax_y2.axhline(y=dv_cutoff, color="blue", linestyle="--")
            axes[0].plot(t, v, "r", linewidth=2.5, alpha=0.6, zorder=10)
            peakHeight_threhold = self.splitViewTab_FP.getParTreePars(
                "Spike detection"
            )["Spike detection parameters"][1]["peak voltage (min=-30;max=150)"][0]
            axes[0].axhline(y=peakHeight_threhold, color="r", linestyle="--")

            axes[0].title.set_text(title)
            axes[0].set_xlabel("Time (S)")
            axes[0].set_ylabel("Voltage (mV)")
            if len(peaks) > 0:
                for p in peaks:
                    axes[0].plot(t[p], v[p], "og", zorder=1)
                # axes[0].legend(['dv/dt','v', 'peak'])

                plt1 = mplWidget.figure.add_subplot(gs[0, 1])
                axes.append(plt1)
                for j in range(len(peaks)):
                    axes[1].plot(waveTime, waves[:, j])
                axes[1].title.set_text(
                    f"N={len(peaks)}, avg.rate: {len(peaks)/dur:.1f} Hz"
                )
                axes[1].set_xlabel("Time (mS)")
                axes[1].set_ylabel("Voltage (mV)")
                axes[1].grid("on")

                plt2 = mplWidget.figure.add_subplot(gs[1, 0])
                axes.append(plt2)
                peakTime = [t[p] for p in peaks]
                axes[2].plot(peakTime, spike_df["width"] * 1e3, "k", linewidth=1.0)
                axes[2].plot(peakTime, spike_df["width"] * 1e3, "g.")
                axes[2].set_xlabel("Time (S)")
                axes[2].set_ylabel("Spike width (mS)")
                xlim_range = axes[0].get_xlim()
                axes[2].grid("on")
                axes[2].set_xlim([xlim_range[0], xlim_range[-1]])

                plt3 = mplWidget.figure.add_subplot(gs[1, 1])
                axes.append(plt3)
                peakTime = [t[p] for p in peaks]
                peakV = [v[p] for p in peaks]
                axes[3].plot(peakTime, peakV, "k", linewidth=1.0)
                axes[3].plot(peakTime, peakV, "g.")
                axes[3].grid("on")
                axes[3].set_xlabel("Time (S)")
                axes[3].set_ylabel("Peak voltage (mv)")
                axes[3].set_xlim([xlim_range[0], xlim_range[-1]])
            else:
                axes[0].legend(["dv/dt", "v"])

            mplWidget.draw()
            #                    print(c)
            mplWidget.figure.tight_layout()

        if plotWidget2 != None:
            self.plot_dvdt_v_phasePlot(waves, waves_dvdt, sweepNumber, plotWidget2)

    def getCurrentIdx(self):
        #        if hasattr(self, 'currentPulseTree'):
        try:
            selected = self.currentPulseTree.selectedItems()
            # update data tree
            sel = selected[0]
            seriesIdx = sel.index
            seriesIdx = list(seriesIdx.copy())
            return seriesIdx
        except:
            return []

    def getSimulatedTraceInASweep(self):
        selected = self.currentPulseTree.selectedItems()
        # update data tree
        sel = selected[0]
        seriesIdx = list(sel.index.copy())  ## get series level index
        if self.parameters["Protocols"]["This serie"]["Type"] == "Firing pattern":
            stimChanIndex = self.parameters["Protocols"]["This serie"]["StimChanID"]
            traceID = self.getCellTraceID_sweep(sel, stimChanIndex)
            seriesIdx.append(traceID)
        else:
            seriesIdx.append(0)

    def rePlotSpikePanels(self):
        seriesIdx = self.getCurrentIdx()
        if len(seriesIdx) > 0:
            if len(seriesIdx) < 3:
                self.showdialog("Must select a trace or a sweep first!")
                return
            elif len(seriesIdx) == 3:
                seriesIdx = self.getSimulatedTraceInASweep()
        else:
            return
        self.plot_splitter.setStretchFactor(0, 1)
        self.plot_splitter.setStretchFactor(1, 1)
        self.plot_alignedSpikes_SingleTrace(seriesIdx, mplWidget=self.spikes_view)

    def detectSpikes_clicked(self):
        self.rePlotSpikePanels()

    def plotSingleTrace(
        self,
        plotHandle,
        trace,
        index,
        myPen,
        xlabelOn=True,
        scaleBarOn=True,
        plotStim=False,
        timeoffset=0,
        highCutOff=None,
        title="",
        removeSpike=False,
        lowCutOff=None,
    ):
        ## plot a single trace
        plotHandle.showGrid(x=True, y=True)
        plotHandle.setLabels(bottom=("Time", trace.XUnit), left=(trace.YUnit))
        title_ = self.sliceName[:-4] + " " + "Series" + str(index[1] + 1)
        data = self.currentPulseTree.bundle.data[index]
        time = np.linspace(
            trace.XStart, trace.XStart + trace.XInterval * (len(data) - 1), len(data)
        )
        time = time + timeoffset
        deltaT = trace.XInterval * (len(data) - 1)
        if self.parameters["Downsampling"] < 1.0 and self.vds == 1:  ## no downsampling
            ## do downsampling for visulization purpose
            time = time[:: int(1 / self.parameters["Downsampling"])]
            data = data[:: int(1 / self.parameters["Downsampling"])]
        if lowCutOff != None:
            useButter = True
        else:
            useButter = False
        if hasattr(self, "events"):
            self.events.data_raw = data
        data = self.bandPass_signal(data, highCutOff, lowCutOff, useButter)
        if removeSpike:
            pv = self.eventParTree_data_view.p.getValues()  ##
            baseline_ = np.mean(data)
            dvdt_th = pv["Spikes"][1]["dv/dt (V/s) - threhold"][0]
            peaks_lv = self.removeStimulationArtifacts(data.copy(), dvdt_th)
            for p in peaks_lv:
                data[p - 20 : p + 20] = baseline_

        pens = []
        if not myPen:
            NChan = len(
                self.currentPulseTree.bundle.pul[index[0]][index[1]][index[2]].children
            )

            myPen = pg.mkPen(
                color=pg.intColor(trace.TraceCount - 1, hues=NChan)
            )  ## default 8 channels. May set a global variable to set channle number
            pens.append(myPen)
        g = plotHandle.plot(time, data, pen=myPen, name=trace.Label)
        plotHandle.autoRange()
        if not xlabelOn:  ## no x labels
            plotHandle.setLabel("bottom", "", units="")
            plotHandle.showAxis("bottom", False)
        if scaleBarOn:
            scale = pg.ScaleBar(size=0.02, suffix="s")
            scale.setParentItem(plotHandle.getViewBox())
            scale.anchor((1, 1), (1, 1), offset=(-20, -20))
        if title == "":
            plotHandle.setTitle(title_ + " sweep" + str(index[2] + 1))
        else:
            plotHandle.setTitle(title_ + title)

        legend = pg.LegendItem(offset=(70, 20))
        legend.setParentItem(plotHandle)
        legend.addItem(g, trace.Label)
        self.plotItemList.append(plotHandle)

        if plotStim and self.parameters["Protocols"]["This serie"]["Type"] in [
            "Firing pattern",
            "Connection",
        ]:
            if self.parameters["plotStim"]:
                self.trace_view.nextRow()  ## for stimulation
                plt_stim = self.trace_view.addPlot(col=0, colspan=4)
                plt_stim.setXRange(0.1, 0.45, padding=0)
                myPen = pg.mkPen(color="w")
                self.plotSingleStimTrace(
                    plt_stim, index, myPen
                )  ## third arugument is for pen. if empty use default pen
                self.trace_view.ci.layout.setRowStretchFactor(0, 3)
                self.trace_view.ci.layout.setRowStretchFactor(1, 1)
                plt_stim.setXLink(plotHandle)

        time_stim, stim, stimInfo = self.currentPulseTree.bundle.stim(
            index
        )  ## assume only 1 channel is stimulating
        traceInfo = self.gatherUsefulTraceInfo(trace)
        if self.parameters["Protocols"]["This serie"]["Type"] in [
            "Firing pattern",
            "Connection",
        ]:
            self.updateRecordingParameterTable(stimInfo, traceInfo)
        return plotHandle, deltaT, time, data

    def gatherUsefulTraceInfo(self, trace):
        RECORDMODE = [
            "In out",
            "On Cell",
            "Out Out",
            "Whole cell",
            "C-Clamp",
            "------",
            "------",
        ]
        modeIdx = int.from_bytes(
            trace.RecordingMode, byteorder=self.currentPulseTree.bundle.endian
        )
        #        print('Mode Idx: %g' % modeIdx)
        traceInfo = {"Recording mode": RECORDMODE[modeIdx]}
        return traceInfo

    def bandPass_signal(self, data, highCutOff=None, lowCutOff=None, useButter=False):
        if self.globalSettings.p.getValues()["data preprocessing"][1][
            "Apply Notch filter"
        ][0]:
            w0 = self.parameters["Notch"]
            # print('notch filter', w0, self.parameters['fs'])
            b, a = signal.iirnotch(w0, 1000.0, self.parameters["fs"])
            f = signal.filtfilt(b, a, data)
        else:
            f = data
        if (
            self.parameters["filter_option"] == 0 and useButter == False
        ):  ## default fourth order Bessel-Thomson filter
            self.Bessel_lowpass(highCutOff, None, lowCutOff)
            f = signal.sosfiltfilt(self.parameters["bessel_filter_sos"], f)
        else:  ## Butter
            self.butter_bandpass(highCutOff, lowCutOff, order=2)
            # print("second order Butter bandpass")
            f = signal.sosfiltfilt(self.parameters["butter_sos"], f)

        return f

    def Bessel_lowpass(self, highCutOff=None, fs=None, lowCutOff=None):
        """
        Lowpass Bessel/Thomson filter
        """
        # pdb.set_trace()
        if fs == None:
            sample_freq = self.parameters["fs"]
        else:
            sample_freq = fs
        if highCutOff == None:
            highCutOff = self.parameters["HF"]
        if highCutOff > 1500:
            filt_coeff = (highCutOff - 500) / (
                sample_freq / 2.0
            )  # filter kHz -> Hz, then get fraction of Nyquist frequency
        else:
            filt_coeff = highCutOff / (sample_freq / 2.0)
        #        print(filt_coeff, highCutOff, sample_freq)
        if filt_coeff < 0 or filt_coeff >= 1:
            print(
                "bessel coeff ({:f}) is outside of valid range [0,1); \
                            cannot filter sampling frequency {:.1f} kHz with \
                            cutoff frequency {:.1f} kHz.".format(
                    filt_coeff, sample_freq / 1e3, highCutOff / 1e3
                )
            )
            print("Using Nyqst frequency for high cutoff")
            filt_coeff = 0.95
        if lowCutOff != None:
            lowCutoff_coeff = lowCutOff / (sample_freq / 2.0)
            filt_coeff = [filt_coeff]
            filt_coeff.insert(0, lowCutoff_coeff)
            self.parameters["bessel_filter_sos"] = signal.bessel(
                2, filt_coeff, btype="bandpass", output="sos"
            )
        else:
            self.parameters["bessel_filter_sos"] = signal.bessel(
                4, filt_coeff, btype="low", output="sos"
            )

    def butter_bandpass(self, highCutOff, lowCutOff, order=2):
        """
        Low pass Butter filter.

        Args:
            - highCutOff (float) : the high cutoff frequency of the filter.
            - fs       (float) : the sampling rate.
            - order      (int) : order of the filter, by default defined to 5.
        """
        # calculate the Nyquist frequency
        nyq = 0.5 * self.parameters["fs"]
        # design filter
        high = highCutOff / nyq
        low = lowCutOff / nyq
        if high >= 1:
            raise ValueError(
                "butter filter high cutoff frequency is higher than Nyquist frequency"
            )
        # returns the filter coefficients: numerator and denominator
        self.parameters["butter_sos"] = signal.butter(
            order, [low, high], btype="band", output="sos"
        )

    def butter_lowpass(self):
        """
        Low pass Butter filter.

        Args:
            - high_cut (float) : the high cutoff frequency of the filter.
            - fs       (float) : the sampling rate.
            - order      (int) : order of the filter, by default defined to 5.
        """
        # calculate the Nyquist frequency
        order = 2
        nyq = 0.5 * self.parameters["fs"]
        # design filter
        high = self.parameters["HF"] / nyq
        if high >= 1:
            raise ValueError(
                "butter filter high cutoff frequency is higher than Nyquist frequency"
            )
        # returns the filter coefficients: numerator and denominator
        self.parameters["butter_sos"] = signal.butter(
            order, high, btype="low", output="sos"
        )

    def updateInterCellDistance(self):
        summary = self.splitViewTab_morph.tables["Summary"]
        nCells = summary.rowCount()
        cellNames = []
        for j in range(nCells):
            cellNames.append(summary.item(j, 0).value)
        distanceList = {"Pair": [], "Distance": []}
        for j in range(nCells):
            for p in range(j):
                distanceList["Pair"].append(cellNames[j] + "-" + cellNames[p])
                distanceList["Distance"].append(
                    summary.item(j, 9).value[p]
                )  ## hard coded column index!
        df_D = pd.DataFrame.from_dict(distanceList, orient="index").transpose()
        df_D = np.array(
            df_D.to_records(index=False)
        )  ## format dataframe for using in QtTable wiget
        self.splitViewTab_morph.tables["Distance (um)"].clear()
        self.splitViewTab_morph.tables["Distance (um)"].appendData(df_D)
        self.splitViewTab_morph.tables["Distance (um)"].show()

    def interceptPerpendecularLine(self, a, b, x0, y0):
        x = (x0 + a * y0 - a * b) / (1 + a**2)
        y = a * x + b
        D = np.sqrt((x - x0) ** 2 + (y - y0) ** 2)
        return x, y, D

    def measureDist2Pia(self):
        if self.pia is not None:
            df = pd.DataFrame(self.pia[:, :2], columns=["x", "y"]).sort_values(
                by="x", inplace=False
            )
            ransac = linear_model.RANSACRegressor()
            ransac.fit(df["x"].values.reshape(-1, 1), df["y"].values.reshape(-1, 1))
            slope = ransac.estimator_.coef_[0][0]
            intercept = ransac.estimator_.intercept_[0]
            xl = self.morphorAxes2D.xaxis.get_data_interval().copy()
            x = np.linspace(xl[0], xl[1], 10).reshape(-1, 1)
            y = ransac.predict(x)
            self.morphorAxes2D.plot(x, y, "r--", label="Pia-line")

            summary = self.splitViewTab_morph.tables["Summary"]
            nCells = summary.rowCount()
            df = self.neuronsData.copy()
            df["Distance to Pia"] = 0
            for j in range(nCells):
                x0, y0 = summary.item(j, 1).value, summary.item(j, 2).value
                x1, y1, D1 = self.interceptPerpendecularLine(slope, intercept, x0, y0)
                df.loc[j, "Distance to Pia"] = D1
                self.morphorAxes2D.plot(
                    [x0, x1], [y0, y1], "g--", label="neuron-pia-line"
                )
            self.splitViewTab_morph.matplotViews["2D"].draw()
            self.splitViewTab_morph.tables["Summary"].clear()
            df = np.array(
                df.to_records(index=False)
            )  ## format dataframe for using in QtTable wiget
            self.splitViewTab_morph.tables["Summary"].appendData(df)
            self.splitViewTab_morph.tables["Summary"].show()

    def morph_analysis_event(self, param, changes):
        pv = self.splitViewTab_morph.getParTreePars("Analysis")
        step_size=pv['Options'][1]['Bin size (um)'][0]
        anlge_step_size=pv['Options'][1]['Angle bin (degree)'][0]*np.pi/180
        useFullRange=pv['Options'][1]['Use full range for density plot'][0]
        showColorBar=pv['Options'][1]['Show color bar for density plot'][0]
        showAxisVal=pv['Options'][1]['Show axis for density plot'][0]
        smoothBins= pv['Options'][1]['Gaussian window size (num. bins)'][0]
        smoothStandardDeviation = pv['Options'][1]['Std of Gaussian kernel (num. bins)'][0]
        for param, change, data in changes:
            childName = param.name()
            if childName in [ "Rotate tree (degree)", "Draw contour", "Ignore diameters", "Scale bar length"]:
                fname = cleanASCfile(
                    self.currentMorphTreeFile
                )  ## to clean not-wanted sections
                self.updateTreeMorphView(fname)
                if fname[-8:] == "_mod.ASC":
                    os.remove(fname)
            if childName in ["Rotate tree (degree)", "Bin size (um)", "Gaussian window size (num. bins)",
             "Std of Gaussian kernel (num. bins)", "Angle bin (degree)", "Use full range for density plot",
             "Show color bar for density plot", "Show axis for density plot"]:
                if self.currentAnMorphView!="":
                    childName = self.currentAnMorphView
            if childName == "Sholl analysis":
                self.update_sholl(step_size=step_size,smoothBins=smoothBins,
                smoothStandardDeviation=smoothStandardDeviation)
                self.currentAnMorphView = childName
            if childName in ["X axis density"]:
                self.update_density('x',step_size=step_size,smoothBins=smoothBins,
                smoothStandardDeviation=smoothStandardDeviation)
                self.currentAnMorphView = childName
            if childName in ["Y axis density"]:
                self.update_density('y',step_size=step_size,smoothBins=smoothBins,
                smoothStandardDeviation=smoothStandardDeviation)
                self.currentAnMorphView = childName
            if childName in ["XY plane density"]:
                self.update_2D_density(step_size=step_size, useFullRange=useFullRange,
                showColorbar=showColorBar, showAxisValues=showAxisVal,smoothBins=smoothBins,
                smoothStandardDeviation=smoothStandardDeviation)
                self.currentAnMorphView = childName
            if childName in ["XY polar density"]:
                self.update_2D_polar_density(step_size=step_size, angle_step = anlge_step_size)
                self.currentAnMorphView = childName
            if childName == "Update cell names":
                self.updateInterCellDistance()
            elif childName == "Distance to Pia":
                self.measureDist2Pia()
            elif childName == "Export High resolution figure":
                dpi = pv["Save options"][1]["DPI"][0]
                format = pv["Save options"][1]["Format"][0]
                figsize = pv["Save options"][1]["Size"][0]
                fig = self.splitViewTab_morph.matplotViews["2D"].getFigure()
                self.saveHighResFigure(fig, dpi, figsize, format)

    def saveHighResFigure(self, fig, dpi, figsize, format):
        fig.set_size_inches(figsize, figsize, forward=True)
        fig.savefig(
            self.currentMorphTreeFile[:-4] + "." + format,
            dpi=dpi,
            format=format,
            metadata=None,
            bbox_inches=None,
            pad_inches=0.1,
            facecolor="auto",
            edgecolor="auto",
            transparent=True,
            backend=None,
        )
        self.showdialog(f"{self.currentMorphTreeFile[:-4]} saved!")

    def spike_detection_event(self, param, changes):
        from pathlib import Path
        for param, change, data in changes:
            childName = param.name()
            if childName == "Clear all tables":
                self.splitViewTab_FP.tables["Sweep features"].clear()
                self.splitViewTab_FP.tables["Spike features"].clear()
                self.splitViewTab_FP.tables["Cell features"].clear()

            if childName == "Save all tables":
                seriesName = (
                    self.splitViewTab_FP.tables["Cell features"].item(0, 0).value
                )
                if self.spikeTableSavePath == "":
                    fileName = self.exportFile(
                        title="Save all tables with prefix",
                        defaultName=seriesName,
                        extension="TSV files (*.tsv)",
                    )
                    if fileName == "":
                        return
                    fparts = Path(fileName).parts
                    self.spikeTableSavePath = "".join(fparts[:-1])
                else:
                    fileName = os.path.join(
                        self.spikeTableSavePath, seriesName + ".tsv"
                    )

                if fileName != "":
                    for f in ["Sweep features", "Spike features", "Cell features"]:
                        fileName1 = fileName[:-4] + "_" + f + ".tsv"
                        with open(fileName1, "w") as fd:
                            fd.write(
                                self.splitViewTab_FP.tables[f].serialize(
                                    useSelection=False
                                )
                            )
                    self.showdialog(f"{seriesName} saved!")

def main(app):
    main = MainWindow(app)
    main.mainFrame.show() 
    sys.exit(app.exec_())

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    main(app)