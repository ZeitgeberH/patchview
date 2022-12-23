# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 15:42:29 2020

@author: mhu
"""

## global setting in the file tab
params = [
    {
        "name": "data preprocessing",
        "type": "group",
        "children": [
            {
                "name": "Notch filter frequency",
                "type": "float",
                "value": 60,
                "limits": (20, 100),
                "default": 60,
                "step": 10,
                "siPrefix": True,
                "suffix": "hz",
            },
            {"name": "Apply Notch filter", "type": "bool", "value": False},
            {
                "name": "High frequency cutoff",
                "type": "float",
                "value": 5e3,
                "limits": (10, 1.2e4),
                "default": 5e3,
                "step": 100,
                "siPrefix": True,
                "suffix": "hz",
            },
            {
                "name": "minimal number of sweeps",
                "type": "int",
                "value": 3,
                "limits": (1, 20),
                "default": 3,
                "step": 1,
            },
        ],
    },
    {
        "name": "Basic parameters",
        "type": "group",
        "children": [
            {"name": "Root Directory", "type": "str", "value": "D:\Mhu\Projects"},
            {"name": "Color map", "type": "colormap"},
        ],
    },
    {
        "name": "Save/Restore parameters",
        "type": "group",
        "children": [
            {"name": "Save State", "type": "action"},
            {
                "name": "Restore State",
                "type": "action",
                "children": [
                    {"name": "Add missing items", "type": "bool", "value": True},
                    {"name": "Remove extra items", "type": "bool", "value": True},
                ],
            },
        ],
    },
]


## Data selection tab for firing pattern
fp_analysis = [
    {
        "name": "Cell 1 Ser.1",
        "type": "group",
        "children": [
            {
                "name": "Depolarize sweep",
                "type": "int",
                "value": 1,
                "limits": (0, 128),
                "default": 1,
                "step": 1,
                "siPrefix": False,
            },
            {
                "name": "Threshold sweep",
                "type": "int",
                "value": 1,
                "limits": (0, 128),
                "default": 1,
                "step": 1,
                "siPrefix": False,
            },
            {
                "name": "Hyperpolar sweep",
                "type": "int",
                "value": 1,
                "limits": (0, 128),
                "default": 1,
                "step": 1,
                "siPrefix": False,
            },
        ],
    },
]

fp_analysis_spikeDetection = [
    {
        "name": "Spike detection parameters",
        "type": "group",
        "children": [
            {
                "name": "High frequency cutoff",
                "type": "float",
                "value": 5e3,
                "limits": (1e3, 1.2e4),
                "default": 5e3,
                "step": 500,
                "siPrefix": True,
                "suffix": "hz",
            },
            {
                "name": "dv/dt cut off (min=1;max=100)",
                "type": "float",
                "value": 20.0,
                "limits": (1, 200),
                "default": 20,
                "step": 1.0,
                "siPrefix": True,
                "suffix": "V/S",
            },
            {
                "name": "peak height (min=2;max=10)",
                "type": "float",
                "value": 2.0,
                "default": 2.0,
                "limits": (2.0, 50.0),
                "step": 1,
                "siPrefix": True,
                "suffix": "mV",
            },
            {
                "name": "peak voltage (min=-30;max=150)",
                "type": "float",
                "value": -20.0,
                "limits": (-30.0, 150.0),
                "default": -20.0,
                "siPrefix": True,
                "suffix": "mV",
                "tip": "Minimal absolute voltage at the peak of detected spikes",
            },
            {
                "name": "max_interval (min=0.001;max=0.02)",
                "type": "float",
                "value": 0.01,
                "limits": (0.001, 0.02),
                "step": 0.001,
                "default": 0.01,
                "siPrefix": True,
                "suffix": "S",
                "tip": "maximum acceptable time between start of spike and time of peak",
            },
            {
                "name": "thresh_frac (min=0.02;max=0.08)",
                "type": "float",
                "value": 0.05,
                "limits": (0.02, 0.08),
                "step": 0.01,
                "default": 0.05,
                "tip": "fraction of average upstroke for threshold calculation",
            },
            {
                "name": "baseline_interval (min=0.05;max=0.15)",
                "type": "float",
                "value": 0.1,
                "limits": (0.05, 0.15),
                "step": 0.01,
                "default": 0.1,
                "siPrefix": True,
                "suffix": "S",
                "tip": "interval length for baseline voltage calculation (before start if start is defined)",
            },
            {
                "name": "baseline_detect_thresh (min=-30;max=150)",
                "type": "float",
                "value": 0.3,
                "limits": (0.1, 0.5),
                "step": 0.1,
                "default": 0.3,
                "siPrefix": True,
                "suffix": "V/S",
                "tip": "Maximal dV/dT threshold for evaluating flatness of baseline region",
            },
            {
                "name": "window length before peak ",
                "type": "float",
                "value": 0.002,
                "limits": (0.001, 0.01),
                "step": 0.001,
                "default": 0.002,
                "siPrefix": True,
                "suffix": "S",
            },
            {
                "name": "window length after peak ",
                "type": "float",
                "value": 0.006,
                "limits": (0.001, 0.01),
                "step": 0.001,
                "default": 0.005,
                "siPrefix": True,
                "suffix": "S",
            },
        ],
    },
    {
        "name": "Visual aesthetics",
        "type": "group",
        "children": [
            {
                "name": "line width",
                "type": "float",
                "value": 1.5,
                "limits": (0.5, 2.5),
                "default": 1.5,
                "step": 0.1,
            },
        ],
    },
    {
        "name": "Tables",
        "type": "group",
        "children": [
            {"name": "Reset table automatically", "type": "bool", "value": True},
            {"name": "Clear all tables", "type": "action"},
            {"name": "Save all tables", "type": "action"},
        ],
    },
]

fp_analysis_spikeDetection_help = [
    {
        "name": "Spike detection parameters",
        "type": "group",
        "children": [
            {
                "name": "High frequency cutoff",
                "type": "str",
                "value": "for spike detection only. Spike waveform high cutoff is set in global parameters",
            },
            {
                "name": "dv/dt cut off (min=1;max=100)",
                "type": "str",
                "value": "minimum dV/dt to qualify as a spike in V/s",
            },
            {
                "name": "Peak height (min=2;max=10)",
                "type": "str",
                "value": "minimum acceptable height from threshold to peak in mV",
            },
            {
                "name": "Peak voltage (min=-30;max=150)",
                "type": "str",
                "value": "minimum acceptable absolute peak level in mV ",
            },
            {
                "name": "Max_interval (min=0.001;max=0.02)",
                "type": "str",
                "value": "maximum acceptable time between start of spike and time of peak in sec",
            },
            {
                "name": "Thresh_frac (min=0.02;max=0.08)",
                "type": "str",
                "value": " fraction of average upstroke for threshold calculation ",
            },
            {
                "name": "Baseline_interval (min=0.05;max=0.15)",
                "type": "str",
                "value": "interval length for baseline voltage calculation",
            },
            {
                "name": "Baseline_detect_thresh (min=-30;max=150)",
                "type": "str",
                "value": "dV/dt threshold for evaluating flatness of baseline region",
            },
            {
                "name": "Window length before peak ",
                "type": "str",
                "value": "duration for upstroke region",
            },
            {
                "name": "Window length after peak ",
                "type": "str",
                "value": "duration for downstroke region",
            },
        ],
    },
]

## This is just a place holder for now
connection_analysis = [
    {
        "name": "Data",
        "type": "group",
        "children": [
            {
                "name": "Depolarize sweep",
                "type": "int",
                "value": 1,
                "limits": (0, 128),
                "default": 1,
                "step": 1,
                "siPrefix": False,
            },
            {
                "name": "Threshold sweep",
                "type": "int",
                "value": 1,
                "limits": (0, 128),
                "default": 1,
                "step": 1,
                "siPrefix": False,
            },
            {
                "name": "Hyperpolar sweep",
                "type": "int",
                "value": 1,
                "limits": (0, 128),
                "default": 1,
                "step": 1,
                "siPrefix": False,
            },
        ],
    },
]


## Data selection tab for event detection
event_detection_preprocessing = [
    {
        "name": "Data selection",
        "type": "group",
        "children": [
            {
                "name": "Sweep",
                "type": "int",
                "value": 1,
                "limits": (1, 1000),
                "default": 1,
                "step": 1,
                "siPrefix": False,
            },
            {
                "name": "Trace",
                "type": "int",
                "value": 1,
                "limits": (1, 1000),
                "default": 1,
                "step": 1,
                "siPrefix": False,
            },
            {
                "name": "Concatenate all sweeps for current channel",
                "type": "bool",
                "value": False,
            },
            {
                "name": "Low frequency cutoff",
                "type": "float",
                "value": 1,
                "limits": (0.1, 99),
                "default": 1,
                "step": 1,
                "siPrefix": True,
                "suffix": "HZ",
            },
            {
                "name": "High frequency cutoff",
                "type": "float",
                "value": 1000,
                "limits": (100, 5000),
                "default": 1000,
                "step": 10,
                "siPrefix": True,
                "suffix": "HZ",
            },
        ],
    },
    {
        "name": "Detrend",
        "type": "group",
        "children": [
            {"name": "To be done", "type": "action"},
        ],
    },
    {
        "name": "PSP Outliers",
        "type": "group",
        "children": [
            {
                "name": "Outlier voltage (mV) - lower bound",
                "type": "float",
                "value": -0.15,
                "limits": (-0.12, -0.001),
                "default": -0.15,
                "step": 0.005,
                "siPrefix": True,
                "suffix": "V",
            },
            {
                "name": "Outlier voltage (mV) - upper bound",
                "type": "float",
                "value": 0.1,
                "limits": (0, 0.5),
                "default": 0.1,
                "step": 0.005,
                "siPrefix": True,
                "suffix": "V",
            },
            {
                "name": "replacement value",
                "type": "list",
                "values": ["bound", "mean", "median", -0.065],
                "value": "bound",
            },
        ],
    },
    {
        "name": "PSC Outliers",
        "type": "group",
        "children": [
            {
                "name": "Outlier voltage (pA) - lower bound",
                "type": "float",
                "value": -500e-12,
                "limits": (-500e-12, 500e-12),
                "default": -500e-12,
                "step": 10e-12,
                "siPrefix": True,
                "suffix": "A",
            },
            {
                "name": "Outlier voltage (pA) - upper bound",
                "type": "float",
                "value": 500e-12,
                "limits": (-500e-12, 500e-12),
                "default": 500e-12,
                "step": 10e-12,
                "siPrefix": True,
                "suffix": "A",
            },
            {
                "name": "replacement value",
                "type": "list",
                "values": ["bound", "mean", "median", 0],
                "value": "bound",
            },
        ],
    },
    {
        "name": "Spikes",
        "type": "group",
        "children": [
            {"name": "Removing spikes", "type": "bool", "value": False},
            {
                "name": "peak voltage(mV) - threshold",
                "type": "float",
                "value": -0.02,
                "limits": (-0.05, 0),
                "default": -0.02,
                "step": 0.002,
                "siPrefix": True,
                "suffix": "V",
            },
            {
                "name": "dv/dt (V/s) - threhold",
                "type": "float",
                "value": 3,
                "limits": (1, 50),
                "default": 3,
                "step": 1,
                "siPrefix": True,
                "suffix": "V/S",
            },
        ],
    },
    {
        "name": "Save/Restore parameters",
        "type": "group",
        "children": [
            {"name": "Save State", "type": "action"},
            {
                "name": "Restore State",
                "type": "action",
                "children": [
                    {"name": "Add missing items", "type": "bool", "value": True},
                    {"name": "Remove extra items", "type": "bool", "value": True},
                ],
            },
        ],
    },
]

## Template tab for event detection
event_detection_template = [
    {
        "name": "Template",
        "type": "group",
        "children": [
            {"name": "Fit", "type": "action"},
            {"name": "Show template", "type": "action"},
            {
                "name": "Window length",
                "type": "float",
                "value": 0.10,
                "limits": (0.001, 0.5),
                "default": 0.10,
                "step": 0.01,
                "siPrefix": True,
                "suffix": "mS",
            },
            # {'name': 'IPSP', 'type': 'bool', 'value': False},
            {
                "name": "Fit options",
                "type": "group",
                "children": [
                    {
                        "name": "function",
                        "type": "list",
                        "values": ["BiExponential", "Hodgkin-Huxley g_Na function"],
                        "value": 0,
                    },
                    {
                        "name": "fitting parameters",
                        "type": "group",
                        "children": [
                            {
                                "name": "Tau fast",
                                "type": "float",
                                "value": 2.5,
                                "limits": (0.01, 100),
                                "default": 2.5,
                                "step": 0.1,
                                "siPrefix": True,
                                "suffix": "mS",
                            },
                            {
                                "name": "Tau slow",
                                "type": "float",
                                "value": 15.0,
                                "limits": (0.01, 100),
                                "default": 15.0,
                                "step": 0.5,
                                "siPrefix": True,
                                "suffix": "mS",
                            },
                            {
                                "name": "Baseline",
                                "type": "float",
                                "value": 0,
                                "limits": (-100, 100),
                                "default": 0,
                                "step": 0.1,
                                "siPrefix": True,
                                "suffix": "mV",
                            },
                            {
                                "name": "Max # of passes",
                                "type": "int",
                                "value": 16,
                                "limits": (10, 100),
                                "default": 16,
                                "step": 2,
                                "siPrefix": False,
                            },
                            {
                                "name": "Max # of iterations per pass",
                                "type": "int",
                                "value": 16,
                                "limits": (10, 100),
                                "default": 16,
                                "step": 2,
                                "siPrefix": False,
                            },
                        ],
                    },
                ],
            },
            {"name": "Average", "type": "bool", "value": True},
            {"name": "Normalizing", "type": "bool", "value": True},
            {"name": "Template to use", "type": "int", "value": 1, "limits": (1, 100)},
            {"name": "Clear all templates", "type": "action"},
        ],
    },
]

## parameters for deconvolution algrithoms and find_peak function
"""
height : number or ndarray or sequence, optional
    Required height of peaks. Either a number, ``None``, an array matching
    `x` or a 2-element sequence of the former. The first element is
    always interpreted as the  minimal and the second, if supplied, as the
    maximal required height.
threshold : number or ndarray or sequence, optional
    Required threshold of peaks, the vertical distance to its neighboring
    samples. Either a number, ``None``, an array matching `x` or a
    2-element sequence of the former. The first element is always
    interpreted as the  minimal and the second, if supplied, as the maximal
    required threshold.
distance : number, optional
    Required minimal horizontal distance (>= 1) in samples between
    neighbouring peaks. Smaller peaks are removed first until the condition
    is fulfilled for all remaining peaks.
prominence : number or ndarray or sequence, optional
    Required prominence of peaks. Either a number, ``None``, an array
    matching `x` or a 2-element sequence of the former. The first
    element is always interpreted as the  minimal and the second, if
    supplied, as the maximal required prominence.
width : number or ndarray or sequence, optional
    Required width of peaks in samples. Either a number, ``None``, an array
    matching `x` or a 2-element sequence of the former. The first
    element is always interpreted as the  minimal and the second, if
    supplied, as the maximal required width.
wlen : int, optional
    Used for calculation of the peaks prominences, thus it is only used if
    one of the arguments `prominence` or `width` is given. See argument
    `wlen` in `peak_prominences` for a full description of its effects.
rel_height : float, optional
    Used for calculation of the peaks width, thus it is only used if `width`
    is given. See argument  `rel_height` in `peak_widths` for a full
    description of its effects.
plateau_size : number or ndarray or sequence, optional
    Required size of the flat top of peaks in samples. Either a number,
    ``None``, an array matching `x` or a 2-element sequence of the former.
    The first element is always interpreted as the minimal and the second,
    if supplied as the maximal required plateau size.
"""
event_deconvPeak_parameters = [
    {
        "name": "Detection method",
        "type": "group",
        "children": [
            {
                "name": "Event Criterion",
                "type": "list",
                "values": ["Template match", "Advanced"],
                "value": 0,
            },
            # {'name': 'Show histogram of deconvolved', 'type': 'action'},
            {"name": "Detect current sweep", "type": "action"},
            {"name": "Detect events for all sweeps", "type": "action"},
            ## this is used when user made some manul changes to auto detected events
            {"name": "Update events", "type": "action"},
        ],
    },
    {
        "name": "Peak parameters for deconvoled trace",
        "type": "group",
        "children": [
            # corresponding to the height in std of voltage. NOT the threhold in the find_peak function
            {
                "name": "Threshold (stdev)",
                "type": "float",
                "value": 3,
                "limits": (1, 100),
                "default": 3,
                "step": 0.2,
                "siPrefix": False,
            },
            # distance: minimal interval between ajacent peaks
            {
                "name": "Distance",
                "type": "int",
                "limits": (1, 200),
                "value": 5,
                "step": 1,
                "default": 5,
                "siPrefix": False,
            },
            # wlen: range to look for prominences
            {
                "name": "Wlen (samples)",
                "type": "int",
                "limits": (5, 1000),
                "value": 50,
                "step": 1,
                "default": 50,
                "siPrefix": False,
            },
            # Prominence: minimal for PSP
            {
                "name": "Prominence (mV)",
                "type": "float",
                "limits": (0.1, 10),
                "value": 0.5,
                "step": 0.1,
                "default": 1,
                "siPrefix": True,
                "suffix": "mV",
            },
            ## for PSC
            {
                "name": "Prominence (pA)",
                "type": "float",
                "limits": (0.1, 10),
                "value": 0.5,
                "step": 0.1,
                "default": 1,
                "siPrefix": True,
                "suffix": "pA",
            },
            {
                "name": "Width",
                "type": "group",
                "children": [
                    {
                        "name": "minimal",
                        "type": "float",
                        "limits": (1, 100),
                        "value": 2,
                        "step": 1,
                        "default": 2,
                        "siPrefix": False,
                    },
                    {
                        "name": "maximal",
                        "type": "float",
                        "limits": (50, 1000),
                        "value": 50,
                        "step": 5,
                        "default": 50,
                        "siPrefix": False,
                    },
                ],
            },
        ],
    },
    {
        "name": "Peak parameters for raw trace",
        "type": "group",
        "children": [
            {
                "name": "Peak Threshold (stdev)",
                "type": "float",
                "value": 3,
                "limits": (1, 10),
                "default": 3,
                "step": 0.5,
                "siPrefix": False,
            },
            {
                "name": "Onset Threshold (stdev)",
                "type": "float",
                "value": 2,
                "limits": (1, 10),
                "default": 2,
                "step": 0.5,
                "siPrefix": False,
            },
            {
                "name": "Prominence (mV or pA)",
                "type": "float",
                "limits": (0.1, 10),
                "value": 1,
                "step": 1,
                "default": 1,
                "siPrefix": False,
            },
            # width: Peak width
            {
                "name": "Peak width",
                "type": "float",
                "limits": (1, 200),
                "value": 1,
                "step": 1,
                "default": 1,
                "siPrefix": True,
                "suffix": "mS",
            },
            # wlen: range to look for prominences
            {
                "name": "Waveform post peak duration",
                "type": "float",
                "limits": (0.001, 0.5),
                "value": 0.03,
                "step": 0.005,
                "default": 0.03,
                "siPrefix": True,
                "suffix": "S",
            },
            # Rel_height:  Used for calculation of the peaks width5
            {
                "name": "Rel_height",
                "type": "float",
                "limits": (0.2, 0.8),
                "value": 0.5,
                "step": 0.05,
                "default": 0.5,
                "siPrefix": False,
            },
            {"name": "Linear fit during onset", "type": "bool", "value": False},
            {
                "name": "viewing window",
                "type": "group",
                "children": [
                    {"name": "Limiting view", "type": "bool", "value": False},
                    {
                        "name": "Window length",
                        "type": "float",
                        "limits": (0.01, 60),
                        "value": 1,
                        "step": 0.01,
                        "default": 1,
                        "siPrefix": True,
                        "suffix": "S",
                    },
                ],
            },
        ],
    },
    {
        "name": "Visualization choices",
        "type": "group",
        "children": [
            ## this is for visulaziation purpose
            {
                "name": "Pre event time (before peak)",
                "type": "float",
                "limits": (0.01, 1.0),
                "value": 0.05,
                "step": 0.01,
                "default": 0.05,
                "siPrefix": True,
                "suffix": "S",
            },
            {
                "name": "Post event time (after peak)",
                "type": "float",
                "limits": (0.01, 1.0),
                "value": 0.4,
                "step": 0.01,
                "default": 0.4,
                "siPrefix": True,
                "suffix": "S",
            },
        ],
    },
]


## Template tab for event detection
event_post_processing = [
    {
        "name": "Events",
        "type": "group",
        "children": [
            {"name": "fit each event", "type": "action"},
            {"name": "show averaged event waveform", "type": "action"},
            {"name": "show event amplitude statistics", "type": "action"},
            {"name": "export event waveforms", "type": "action"},
        ],
    },
]


event_detection_helps = [
    {
        "name": "Make templates",
        "type": "group",
        "children": [
            {"name": "Add", "type": "str", "value": "@Ctrl + left click"},
            {
                "name": "Delete",
                "type": "str",
                "value": "@Alt + left click on the red line",
            },
            {
                "name": "Move",
                "type": "str",
                "value": "@Left click red/green line then Drag",
            },
        ],
    },
    {
        "name": "View individual event",
        "type": "group",
        "children": [
            {
                "name": "Arrow key Up or Down",
                "type": "str",
                "value": "select any cell in the <InUse> column. Then navigate up or down",
            },
        ],
    },
    {
        "name": "Edit events",
        "type": "group",
        "children": [
            {
                "name": "Select/Deselect",
                "type": "str",
                "value": "@Left click on the green marker OR press enter in the <Event list> table",
            },
            {"name": "Add new", "type": "str", "value": "@Shift + left click"},
        ],
    },
    {
        "name": "Fitting events problems",
        "type": "group",
        "children": [
            {
                "name": "solution0",
                "type": "str",
                "value": "limit the fitting window to the most significant part",
            },
            {
                "name": "solution1",
                "type": "str",
                "value": "Add more high quality events and try again",
            },
            {"name": "solution2", "type": "str", "value": "Try different parameters"},
        ],
    },
]

Morphor_analysis = [
    {
        "name": "Options",
        "type": "group",
        "children": [
            {
                "name": "Scale bar length",
                "type": "float",
                "limits": (10, 50),
                "value": 50,
                "step": 5,
                "default": 50,
                "siPrefix": False,
            },
            {
                "name": "Rotate tree (degree)",
                "type": "float",
                "limits": (-180, 180),
                "value": 0,
                "step": 5,
                "default": 0,
                "siPrefix": False,
            },
            {
                "name": "Bin size (um)",
                "type": "float",
                "limits": (1, 50),
                "value": 5,
                "step": 1,
                "default": 0,
                "siPrefix": False,
            },
            {
                "name": "Gaussian window size (num. bins)",
                "type": "int",
                "limits": (6, 18),
                "value": 10,
                "step": 2,
                "default": 10,
                "siPrefix": False,
            },
            {
                "name": "Std of Gaussian kernel (num. bins)",
                "type": "int",
                "limits": (2, 4),
                "value": 2,
                "step": 1,
                "default": 2,
                "siPrefix": False,
            },
            {
                "name": "Angle bin (degree)",
                "type": "float",
                "limits": (1, 30),
                "value": 15,
                "step": 5,
                "default": 15,
                "siPrefix": False,
            },
            {"name": "Use full range for density plot", "type": "bool", "value": True},
            {"name": "Show color bar for density plot", "type": "bool", "value": True},
            {"name": "Show axis for density plot", "type": "bool", "value": True},
            {"name": "Ignore diameters", "type": "bool", "value": False},
            {"name": "Draw contour", "type": "bool", "value": True},
            {"name": "Update cell names", "type": "action"},
        ],
    },
    {
        "name": "Measurment",
        "type": "group",
        "children": [
            {"name": "Sholl analysis", "type": "action"},
            {"name": "X axis density", "type": "action"},
            {"name": "Y axis density", "type": "action"},
            {"name": "XY plane density", "type": "action"},
            {"name": "XY polar density", "type": "action"},
            {"name": "Distance to Pia", "type": "action"},
        ],
    },
    {
        "name": "Save options",
        "type": "group",
        "children": [
            {
                "name": "DPI",
                "type": "float",
                "limits": (100, 1000),
                "value": 300,
                "step": 50,
                "default": 300,
                "siPrefix": False,
            },
            {
                "name": "Size",
                "type": "int",
                "limits": (10, 100),
                "value": 60,
                "step": 5,
                "default": 60,
                "siPrefix": False,
            },
            {
                "name": "Format",
                "type": "list",
                "values": ["svg", "eps", "png"],
                "value": "svg",
            },
            {"name": "Export High resolution figure", "type": "action"},
        ],
    },
]

Morphor_legend = [
    {
        "name": "Color scheme",
        "type": "group",
        "children": [
            {"name": "Basal dendtrite", "type": "str", "value": "Purple"},
            {"name": "Apical dendtrite", "type": "str", "value": "Blue"},
            {"name": "Axon", "type": "str", "value": "Red"},
            {"name": "Soma", "type": "str", "value": "Black"},
        ],
    },
    {
        "name": "Help for options",
        "type": "group",
        "children": [
            {
                "name": "ignore diamters",
                "type": "str",
                "value": "When set True, will render each section's diamter accordingly",
            },
        ],
    },
]
