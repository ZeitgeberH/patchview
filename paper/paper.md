---
title: 'PatchView: A Python Package for Patch-clamp Data Analysis and Visualization'
tags:
  - Python
  - neuroscience
  - patch-clamp
  - mini-event
  - Heka
  - Axon instruments
authors:
  - name: Ming Hu
    orcid: 0000-0002-5254-5315
    affiliation: "1, 2"
  - name: Xiaolong Jiang
    affiliation: "1, 2, 3"
affiliations:
  - name: Department of Neuroscience, Baylor College of Medicine, Houston, TX
    index: 1
  - name: Jan and Dan Duncan Neurological Research Institute at Texas Children’s Hospital, Houston,TX, USA
    index: 2
  - name: Department of Ophthalmology, Baylor College of Medicine, Houston, TX
    index: 3
date: 19 July 2022
bibliography: paper.bib
---

# Summary

``PatchView`` is a Python package for visualizing and analyzing multi channel patch-clamp (multi-patch) data [@Neher:1992], including intrinsic membrane properties and firing pattern analysis [@Komendantov:2019; @McKinney:2022], mini-event analysis[@Clements:1997], synaptic connection detection [@Jiang:2015], morphological analysis [@NeuroM:2022] and more. 

``PatchView`` integrates multiple open-source tools and provides an intuitive graphic user interface (GUI)  for navigating through the data files and analyses pipeline. It is aimed to enable users to perform most analysis quickly for the data collected in a typical patch-clamp experiment without managing a Python environment or writing any Python scripts.

Please refer to the [full documentation](https://patchview-doc.readthedocs.io/en/latest/) and the [source code](https://github.com/ZeitgeberH/patchview) for detailed information.

# Statement of Need
**Research purpose**

 ``PatchView`` was designed to be used by neuroscientists for handling electrophysiology data recorded from cells in alive tissues (such as brain slice or cultured cells) in patch-clamp experiments. The target audience is anyone working with patch-clamp data, with or without programming skills. 

**Problems solved**

 Main functionalities of ``PatchView``:

* Importing both [Heka](https://www.heka.com/about/about_main.html#smart-ephys) data and [Axon Instruments](https://www.moleculardevices.com/products/axon-patch-clamp-system#gref) data. Exporting to Python pickle file or NWB (Neurodata Without Borders) file format [@NWB:2022].
* Visualizing single and multiple traces with zoom, pan operations.
* Automatically sorting experiments data according to predefined labels.
* Performing analysis on intrinsic membrane properties, action potential detection, firing pattern analysis  (\autoref{fig:gui}).
* Synaptic connection analysis (\autoref{fig:connection}).
* Mini-event analysis (\autoref{fig:miniEvent-gui}).
* Visualizing and quantification of neuron's morphological reconstruction from [Neurolucida](https://www.mbfbioscience.com/neurolucida360)

**Compares to other commonly-used packages** 

Commercial software such as [Patchmaster](https://www.heka.com/downloads/downloads_main.html#down_patchmaster_next), [Clampfit](https://www.moleculardevices.com/products/axon-patch-clamp-system/acquisition-and-analysis-software/pclamp-software-suite#Overview) provide rich functions for handling this type of data. However, the former only supports Heka data, while the latter only support Axon Instruments data. PatchView supports both .dat (from Heka) and .abf format (from Axon Instruments). To facilitate data sharing, PatchView could export imported data as NWB file. 

``Stimfit``[@Guzman:2014] is a well-known python package for dealing with pre/post synaptic events in single channel. Compared to ``Stimfit``, ``PatchView`` also provides more intuitive user interface (\autoref{fig:miniEvent-gui}) and more native support for Heka ``dat`` file. For instance, most of the time, a Heka ``dat`` file may host data recorded in multiple experiments from a single cell. These experiments may need to be analyzed with different pipelines. ``PatchView`` leverages the labels (those are usually predefined by experimenters for each protocols) associated with each experiment to automatically sort data into its corresponding category. Sorted data can be directly submitted to downstream pipeline, such as mini-event analysis. 

In addition to that, other software mentioned above does not support analysis for data recorded from multiple channels simultaneously. PatchView supports multi-channel analysis, such as synaptic connection analysis (\autoref{fig:connection}). 

# Figures

![`PatchView` graphic user interface. Top left: file browsers. Middle left: experiments data inner trial selection.
Right: multiple plots during a typical analysis). The toolbar is seen on the far right edge of the interface.
\label{fig:gui}](gui.png){ width=100% }

![`PatchView` mini-event analysis GUI. Top left: large-scale view of the trace. The top yellow line is the original trace. The bottom white trace is a deconvolved trace with separate threshold (red horizontal line). Top right: currently selected event. Bottom left: parameters panel. Bottom right: table for detected events' detailed quantification. 
\label{fig:miniEvent-gui}](event.png){ width=100% }

![`PatchView` synaptic connection analysis GUI. Middle left: data series browser. It lists series recorded in a set of synaptic connection experiments from a group of neurons simultaneously recorded in the same slice. Top right: the left plot shows averaged traces of stimuli invoked responses. titles show the channel names of pre and post neurons. The right plot is a graph representation of detected connections. Bottom right panel: quantification of detected connections: including peak voltage, peak height, peak delay, trial consistency. 
\label{fig:connection}](connection.png){ width=100% }

# Acknowledgements

We acknowledge open-source packages (please see the [documentation]((https://patchview-doc.readthedocs.io/en/latest/)) for the detailed list) that are used in PatchView. We would like to thank lab members in Dr. Xiaolong Jiang's lab for their feedback in using PatchView for their data analysis. The work is supported by the grants R01 MH122169 and R01 MH120404 (to XJ) from the National Institutes of Health.

# References
