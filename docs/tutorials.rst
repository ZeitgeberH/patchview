================
GUI Tutorials
================
Here, we walk through the main steps to use PatchView.

**Loading data**
-----------------
Patchview provides multiple ways to import data.

*import single file*
^^^^^^^^^^^^^^^^^^^^^^
Use the file browser to expand (click the arrow) directories. Once located the file (.dat or .abf), double clicking
the file. The contents of this file will be represented as a hierarchical tree in a tab named "Current tree" in the
left panel. Click the arrow at start of the name to expand. 

*import multiple files simultaneously*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Note**: this method assume that you follow a naming convention: date-slice#-cell#. For example, 20220202s1c1.dat, means
this file was recorded on date 2022/02/02 from slice 1 and cell 1. 

Sometimes, you may record multiple cells from a single slice. Each cell's recording may be stored in a separate file. 
Or, they were stored in a single file. And You may also record multiple slices from a single day.
The following method would enable you importing cells from a single slice (possibly in multiple files). 

Right click the folder that contains data files, a context menu would pop out. Clicking "import files", a dialog box with
a dropdown list of each individual slice name and "All slice". Selecting the former will load files from same slice.
Selecting 'All slices' will load all files in the current folder. 

An advantage of loading data in this way is that PatchView will try to sort series by their protocol. Some common protocols
in a multi-patch experiment test for firing pattern, connection and spontaneous event. In Heka patchmaster, you can specify
a unique label for each protocol: for example, 'FP' for firing pattern; 'EPSP' for excitatory post synaptic potential.
PatchView would read those labels and sort each recorded series into corresponding tabs shown in the first middle panel.
Currently three tabs are available: "Firing pattern", "Connection", "Spontaneous".

**Interactive with figure planel**
-------------------------------------

**Left button**:  Left mouse button have two modes of interaction: Rect mode |mouseMode1| and Pan mode |mouseMode2|. 

* Rect mode: Left click and drag a rectangle around the region of interest to zoom in.
* Pan mode:  Interacts with items in the scene (select/move objects, etc). If there are no movable objects under the mouse cursor, then dragging with the left button will pan the scene instead. 

Rect mode is default. Click the corresponding icon in the toolbar (far right column in GUI) to switch between
those two modes.

**Right Button**:

* Right button drag: Scales the scene. Dragging left/right scales horizontally; dragging up/down scales vertically (although some scenes will have their x/y scales locked together). If there are x/y axes visible in the scene, then right-dragging over the axis will _only_ affect that axis. 
* Right button clicks: Clicking the right button in most cases will show a context menu with a variety of options depending on the object(s) under the mouse cursor. 

.. |mouseMode1| image:: resources/images/rectangle.png
    :height: 25px
.. |mouseMode2| image:: resources/images/navigation.png
    :height: 25px

**Firing Pattern analysis**
----------------------------
*single series FP analysis*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Right clicking the Series in the "Current tree" tab, choose "Firing pattern analysis". This will trigger firing pattern
analysis for the whole series. The results are shown in two formats in the "Firing pattern" tab: 

* the top panel shows aligned spike waveform in the first 16 sweeps that have at least one spikes.
* the bottom left panel have three tabs. Each tab holds a table for sweep-wise, spike-wise and cell-wise features.

*multiple series FP analysis*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Load all series of interests using method described in **import multiple files simultaneously** section. Invoke firing 
pattern analysis from the contex menu of the "Firing pattern" tab in the first middle left panel.

**Monosynaptic connection detection**
-------------------------------------------
Same loading method as the multiple series FP analysis. Invoke **Synaptic connections (Batch)** from the context menu of 
the "Connection" tab in the first middle left panel.

**Postsynaptic event (PSP/PSC) analysis**
-------------------------------------------
Load file that has spontaneous event series. Right click the series of interest to invoke context menu. Left click
"Spontaneous event analysis" will load the series and activate the "Event detection" tab. 

The top two panels show the original trace (left) and template or detected events (right). The bottom left panel have
several tabs, which hold parameters widget and actionable buttons. Major operations are done through this panel. The
bottom right panel show tables of various intermediate results and figures.

To start, activate the "Preprocessing" tab at the bottom left panel. This tab enables select sweep and traces within
currently loaded series. It also has various parameters for preprocessing.

*Specify templates*
^^^^^^^^^^^^^^^^^^^^^^
Next, activate "Template" tab. Move cursor into the left upper panel to specify several typical events as template. See the following
picture for how to do this. Once templates are collected, click "Fit" in the tab to get a bi-exponential fit of the averaged event.
The top right panel will now show the averaged the template and its fit together with time constant of rising and decay phase time constants.
The "Template Fit" tab in the bottom right area will represent extra information about the fit.

.. image:: resources/images/event_template.png
    :width: 800
    :alt: Alternative text

*detect events*
^^^^^^^^^^^^^^^^^^^^^^
Once the template fit is finished, activate "Peak detection" tab. Then click "Detect current sweep" button to analyze events
for current sweep, or "Detect events for all sweeps" for all sweeps.  See the following graphic guide. 

.. image:: resources/images/event_sweep.png
    :width: 800
    :alt: Alternative text

*visualize and manually curate events*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. image:: resources/images/event_curate.png
    :width: 800
    :alt: Alternative text

*Postprocessing and exporting*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The "Event postprocessing" tab summarize the results with event histograms. The wave forms and tables can also be exported
from here.

.. image:: resources/images/event_histExport.png
    :width: 800
    :alt: Alternative text

**Morphological data**
---------------------------
Click the |tree| icon in the toolbar section to load ASC file.  

The bottom left panels have options for draw contours or update cell names (for multiple neurons).
To export a high resolution image, use the export option in this panel,
in stead of the plot widget's build-in save button.

.. |tree| image:: resources/images/tree.png
    :height: 25px

.. image:: resources/images/morphor_tree.png
    :width: 800
    :alt: Alternative text
