# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 20:37:12 2020

@author: mhu
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.io as io
import scipy.signal as signal
from scipy.optimize import curve_fit
from scipy import integrate
import logging
import importlib
import os
from patchview.ephys import ephys_extractor as efex
from patchview.ephys import ephys_features as ft
import pyqtgraph as pg

sns.set_style()
import matplotlib as mpl

mpl.rcParams["pdf.fonttype"] = 42

from sklearn import linear_model

ransac = linear_model.RANSACRegressor()

class extraEphys(object):
    def __init__(
        self,
        time,
        data,
        datIndex=None,
        stimInfo=None,
        title=None,
        filterHighCutFreq=None,
        dv_cutoff=20.0,
        max_interval=0.005,
        min_height=2.0,
        min_peak=-30.0,
        thresh_frac=0.05,
        baseline_interval=0.1,
        baseline_detect_thresh=0.3,
    ):
        self.time = time  ## time in second
        self.data = data * 1000  ## memebrane voltage in unit of mV
        self.stim_start = (
            stimInfo[0][1]["start"] * stimInfo[0][1]["sampleInteval"]
        )  ## stimuli start
        self.stim_end = (
            stimInfo[0][1]["end"] * stimInfo[0][1]["sampleInteval"]
        )  ## stimuli end
        self.title = title  ## contain the file name + lable of current node
        current_amp = []  ## stimuli amplitude
        for x in stimInfo:
            current_amp.append(
                x[1]["amplitude"]
            )  ## we are using the relative current. As the ephys_extrator assume baseline current is zero pA)  ## pA
        self.current = current_amp
        print(current_amp)
        self.current0_index = np.argsort(np.abs(current_amp))[0]
        self.current_step = np.diff(current_amp)[0]  ## assume constant step
        self.Vholding = x[1]["Vholding"]  ## holding current
        self.sampleRate = 1 / stimInfo[0][1]["sampleInteval"]
        if filterHighCutFreq == None:
            if self.sampleRate > 2.1e4:  ## sampling rate larger than 20K hz
                self.filterHighCutFreq = 10  ## 10 KHz
            else:
                self.filterHighCutFreq = (
                    np.floor((1 / stimInfo[0][1]["sampleInteval"] - 100.0) / 1e3) / 2
                )  ## KHz
        else:
            self.filterHighCutFreq = filterHighCutFreq

        self.spikeList = {}
        self.stimInfo = stimInfo  ## store information about stimuli
        self.datIndex = (
            datIndex  ## store coordinates within the dat file: [group, series]
        )
        ## self.df has the features of each spikes
        ## self.df_related_features, has sweep wise spike features
        self.df, self.df_related_features = self.extract_spike_features(
            dv_cutoff=dv_cutoff,
            max_interval=max_interval,
            min_height=min_height,
            min_peak=min_peak,
            thresh_frac=thresh_frac,
            baseline_interval=baseline_interval,
            baseline_detect_thresh=baseline_detect_thresh,
        )
        self.Cell_Features = self.get_cell_features()

    def data_preparation_pickle(pickleObject):
        """data_preparation_pickled object import the data exported from patchviewer and returns the voltage
        traces, stimulus current magnitudes for all traces and the time trace.

        Parameters
        ----------
        pickle_Object : Pathname to a hdf5 object with keys:
            "data": 3D array: time X chan X sweep
            "time": 1D array
            "stimTime": 1D array
            "stimData": 2D array: time X sweep
            "stimInfo": list with stimulus paramters for each sweep

        Returns
        -------
        time : numpy 1D array of time points (s)
        voltage : numpy 2D array of voltage traces (mV)
        current : numpy 1D array of current stimulus magnitudes
        """
        import pickle

        with open(pickleObject, "rb") as f:
            Ldata = pickle.load(f)
        time = Ldata["time"]  ## seconds
        voltage = Ldata["data"] * 1000  ## to mv
        stim_start = (
            Ldata["stimInfo"][0][1]["start"] * Ldata["stimInfo"][0][1]["sampleInteval"]
        )
        stim_end = (
            Ldata["stimInfo"][0][1]["end"] * Ldata["stimInfo"][0][1]["sampleInteval"]
        )
        # current = Ldata['stimData'][current_start:current_end,:]

        current_amp = []
        for x in Ldata["stimInfo"]:
            current_amp.append(
                x[1]["ampitude"] - x[1]["Vholding"]
            )  ## we are using the relative current. As the ephys_extrator assume baseline current is zero pA
        current0_index = np.argsort(np.abs(current_amp))[0]
        current_step = np.diff(current_amp)[0]  ## assume constant step
        return (
            time,
            voltage,
            stim_start,
            stim_end,
            np.array(current_amp),
            current0_index,
            current_step,
        )

    def extract_spike_features(
        self,
        dv_cutoff=20.0,
        max_interval=0.005,
        min_height=2.0,
        min_peak=-30.0,
        thresh_frac=0.05,
        baseline_interval=0.1,
        baseline_detect_thresh=0.3,
    ):
        """Analyse the voltage traces and extract information for every spike (returned in df), and information for all the spikes
        per current stimulus magnitude.
        Returns
        -------
        df : DataFrame with information for every detected spike (peak_v, peak_index, threshold_v, ...)
        df_related_features : DataFrame with information for every possible used current stimulation magnitude
        """
        current = self.current
        time = self.time
        start = self.stim_start
        end = self.stim_end
        voltage = self.data

        df = pd.DataFrame()
        df_related_features = pd.DataFrame()

        with pg.ProgressDialog(
            "Extracting spikes ", maximum=len(current), busyCursor=True, nested=False
        ) as dlg:
            for c, curr in enumerate(current):
                #            print("sweep: %g, current: %g"%(c, curr))
                current_array = curr * np.ones_like(time)
                start_index = (
                    np.abs(time - start)
                ).argmin()  # Find closest index where the injection current starts
                end_index = (
                    np.abs(time - end)
                ).argmin()  # Find closest index where the injection current ends
                current_array[:start_index] = 0
                current_array[end_index : len(current_array)] = 0
                EphysObject = efex.EphysSweepFeatureExtractor(
                    t=time,
                    v=voltage[:, c],
                    i=current_array,
                    start=start,
                    end=end,
                    filter=self.filterHighCutFreq,
                    dv_cutoff=dv_cutoff,
                    max_interval=max_interval,
                    min_height=min_height,
                    min_peak=min_peak,
                    thresh_frac=thresh_frac,
                    baseline_interval=baseline_interval,
                    baseline_detect_thresh=baseline_detect_thresh,
                )
                EphysObject.process_spikes()

                # Adding peak_height (mV) + code for maximum frequency determination (see further)
                spike_count = 0
                if EphysObject._spikes_df.size:
                    EphysObject._spikes_df["peak_height"] = (
                        EphysObject._spikes_df["peak_v"].values
                        - EphysObject._spikes_df["threshold_v"].values
                    )
                    spike_count = EphysObject._spikes_df["threshold_i"].values.size
                    EphysObject._spikes_df.insert(
                        0, "sweepCurrent", curr
                    )  ## add sweep current for this spike if any
                    EphysObject._spikes_df.insert(
                        0, "sweepCount", c
                    )  ## add sweep id for this spike if any

                df = pd.concat([df, EphysObject._spikes_df], sort=True)

                # Some easily found extra features
                df_features = EphysObject._sweep_features

                # Adding spike count
                df_features.update({"spike_count": spike_count})

                # Adding spike frequency adaptation (ratio of spike frequency of second half to first half)
                SFA = np.nan
                half_stim_index = ft.find_time_index(
                    time, np.float(start + (end - start) / 2)
                )
                if (
                    spike_count > 5
                ):  # We only consider traces with more than 8.333 Hz = 5/600 ms spikes here
                    # but in the end we only take the trace with the max amount of spikes

                    if (
                        np.sum(
                            df.loc[df["threshold_i"] == curr, :]["threshold_index"]
                            < half_stim_index
                        )
                        != 0
                    ):
                        SFA = np.sum(
                            df.loc[df["threshold_i"] == curr, :]["threshold_index"]
                            > half_stim_index
                        ) / np.sum(
                            df.loc[df["threshold_i"] == curr, :]["threshold_index"]
                            < half_stim_index
                        )

                df_features.update({"SFA": SFA})

                # Adding current (pA)
                df_features.update({"current": curr})

                # Adding membrane voltage (mV)
                df_features.update(
                    {"resting_membrane_potential": EphysObject._get_baseline_voltage()}
                )

                # Adding voltage deflection to steady state (mV)
                # voltage_deflection_SS = ft.average_voltage(voltage[:, c], time, start = end - 0.1, end = end)
                # df_features.update({'voltage_deflection': voltage_deflection_SS}) ## this won't work if tau is very large

                (
                    voltage_deflection_SS,
                    voltage_deflection_i,
                ) = EphysObject.voltage_deflection()  #
                df_features.update({"voltage_deflection": voltage_deflection_SS})

                # Adding input resistance (MOhm)
                input_resistance = np.nan
                if (
                    not ("peak_i" in EphysObject._spikes_df.keys()) and not curr == 0
                ):  # We only calculate input resistances
                    #                    try:
                    #                        input_resistance, xx, yy = self.getInputResistance()  ## get the slope of a linear fit between I-V
                    #                    except:
                    #                        input_resistance = np.nan                                                           # from traces without APs
                    input_resistance = (
                        np.abs(
                            voltage_deflection_SS - EphysObject._get_baseline_voltage()
                        )
                        * 1000
                    ) / np.abs(curr)
                    if input_resistance == np.inf:
                        input_resistance = np.nan
                df_features.update({"input_resistance": input_resistance})

                # Adding membrane time constant (s) and voltage plateau level for hyperpolarisation paradigms
                # after stimulus onset
                tau = np.nan
                E_plat = np.nan
                sag_ratio = np.nan
                if (
                    voltage_deflection_SS - df_features["resting_membrane_potential"]
                    < -5
                    or curr < 0
                ):  # We use hyperpolarising steps as required in the object function to estimate the
                    # membrane time constant and E_plateau
                    while True:
                        try:
                            tau = (
                                EphysObject.estimate_time_constant()
                            )  # Result in seconds!
                            break
                        except TypeError:  # Probably a noisy bump for this trace, just keep it to be np.nan
                            break
                    E_plat = ft.average_voltage(
                        voltage[:, c], time, start=end - 0.1, end=end
                    )
                    sag, sag_ratio = EphysObject.estimate_sag()
                df_features.update({"tau": tau})
                df_features.update({"E_plat": E_plat})
                df_features.update({"sag_ratio": sag_ratio})

                # For the rebound and sag time we only are interested in the lowest (-200 pA (usually)) hyperpolarisation trace
                rebound = np.nan
                sag_time = np.nan
                sag_area = np.nan

                if curr < -19:
                    baseline_interval = start  # To calculate the SS voltage
                    v_baseline = EphysObject._get_baseline_voltage()

                    start_index = ft.find_time_index(time, start)
                    end_index = ft.find_time_index(time, end)
                    if (
                        np.flatnonzero(voltage[end_index:, c] > v_baseline).size == 0
                    ):  # So perfectly zero here means
                        # it did not reach it
                        rebound = 0
                    else:
                        index_rebound = (
                            end_index
                            + np.flatnonzero(voltage[end_index:, c] > v_baseline)[0]
                        )
                        if not (
                            time[index_rebound] > (end + 0.15)
                        ):  # We definitely have 150 ms left to calculate the rebound
                            rebound = (
                                ft.average_voltage(
                                    voltage[
                                        index_rebound : index_rebound
                                        + ft.find_time_index(time, 0.15),
                                        c,
                                    ],
                                    time[
                                        index_rebound : index_rebound
                                        + ft.find_time_index(time, 0.15)
                                    ],
                                )
                                - v_baseline
                            )
                        else:  # Work with whatever time is left
                            if time[-1] == time[index_rebound]:
                                rebound = 0
                            else:
                                rebound = (
                                    ft.average_voltage(
                                        voltage[index_rebound:, c], time[index_rebound:]
                                    )
                                    - v_baseline
                                )

                    peak_index = start_index + np.argmin(
                        voltage[start_index:end_index, c]
                    )
                    v_peakMax, peak_index2 = EphysObject.voltage_deflection("max")
                    v_steady = ft.average_voltage(
                        voltage[:, c], time, start=end - baseline_interval, end=end
                    )
                    try:
                        first_index = (
                            start_index
                            + np.flatnonzero(
                                voltage[start_index:peak_index, c] < v_steady
                            )[0]
                        )
                        # First time SS is reached after the max voltage deflection downwards in the sag
                        if (
                            np.flatnonzero(
                                voltage[peak_index:end_index, c] > v_steady
                            ).size
                            == 0
                        ):
                            second_index = end_index
                        else:
                            second_index = (
                                peak_index
                                + np.flatnonzero(
                                    voltage[peak_index:end_index, c] > v_steady
                                )[0]
                            )
                        sag_time = time[second_index] - time[first_index]
                        sag_area = -integrate.cumtrapz(
                            voltage[first_index:second_index, c],
                            time[first_index:second_index],
                        )[-1]

                    except:
                        pass  # print('sag calculation error')

                burst_metric = np.nan
                # print(c)
                if spike_count > 5:
                    burst = EphysObject._process_bursts()
                    if len(burst) != 0:
                        burst_metric = burst[0][0]

                df_features.update({"rebound": rebound})
                df_features.update({"sag_time": sag_time})
                df_features.update({"sag_area": sag_area})
                df_features.update({"burstiness": burst_metric})

                df_related_features = pd.concat(
                    [df_related_features, pd.DataFrame([df_features])], sort=True
                )
                dlg += 1
                if dlg.wasCanceled():
                    print("Canceled stage %s" % c)
                    break
        return df, df_related_features

    def get_cell_features(self, axis=None):
        """Analyse all the features available for the cell per spike and per current stimulation magnitude. Extract typical
        features includig the resting membrane potential (Vm, mV), the input resistance (R_input, MOhm), the membrane time constant
        (tau, ms), the action potential threshold (AP threshold, mV), the action potential amplitude (AP amplitude, mV),
        the action potential width (AP width, ms), the afterhyperpolarisation (AHP, mV), the afterdepolarisation
        (ADP, mV), the adaptation index (AI, %) and the maximum firing frequency (max freq, Hz).

        Parameters
        ----------
        axis : figure axis object (optional, None by default)
        Returns
        ----------
        Cell_Features : DataFrame with values for all required features mentioned above
        """
        current = self.current
        time = self.time
        start = self.stim_start
        end = self.stim_end
        voltage = self.data
        df = self.df
        df_related_features = self.df_related_features
        curr_index_0 = self.current0_index

        current_step = np.diff(current)[0]
        tau_array = df_related_features["tau"].dropna().values
        tau = np.nanmean(tau_array) * 1000
        #        Rm_array = df_related_features['resting_membrane_potential'][:curr_index_0].dropna().values
        #        Ri_array = df_related_features['input_resistance'][:curr_index_0].dropna().values
        Rm_array = df_related_features["resting_membrane_potential"].dropna().values
        Rm = np.nanmean(Rm_array)

        #        Ri_array = df_related_features['input_resistance'].dropna().values
        #        Ri = np.nanmean(Ri_array)
        #        print(Rm_array)
        #        print(Ri_array)
        try:
            (
                Ri,
                xx,
                yy,
                v_rest,
            ) = self.getInputResistance()  ## get the slope of a linear fit between I-V
        except:
            Ri = np.nan
            v_rest = np.nan
        sag_ratio = df_related_features["sag_ratio"].values[
            0
        ]  # Steepest hyperpolarising trace used
        rebound = df_related_features["rebound"].values[
            0
        ]  # Steepest hyperpolarising trace used
        sag_time = df_related_features["sag_time"].values[
            0
        ]  # Steepest hyperpolarising trace used
        sag_area = df_related_features["sag_area"].values[
            0
        ]  # Steepest hyperpolarising trace used

        if not df.empty:
            # The max amount of spikes in 600 ms of the trace showing the max amount of spikes in 600 ms
            max_freq = np.max(df_related_features["spike_count"].values)
            # Take the first trace showing this many spikes if there are many
            current_max_freq = np.flatnonzero(
                df_related_features["spike_count"].values >= max_freq
            )[0]

            # Rebound firing
            rebound_spikes = 0
            EphysObject_rebound = efex.EphysSweepFeatureExtractor(
                t=time[ft.find_time_index(time, end) :],
                v=voltage[ft.find_time_index(time, end) :, 0],
                i=np.zeros_like(time[ft.find_time_index(time, end) :]),
                start=end,
                end=time[-1],
                filter=self.filterHighCutFreq,
            )
            EphysObject_rebound.process_spikes()
            if EphysObject_rebound._spikes_df.size:
                rebound_spikes = EphysObject_rebound._spikes_df[
                    "threshold_i"
                ].values.size

            # Check if there are APs outside the stimilation interval for the highest firing trace. When true, continue looking
            # for lower firing traces untill it shows none anymore. We don't want the neuron to have gone 'wild' (i.e. dying).
            current_max_freq_initial = current_max_freq  # to see for which cells there is going to be much of a difference
            artifact_occurence = False
            EphysObject_end = efex.EphysSweepFeatureExtractor(
                t=time[ft.find_time_index(time, end) :],
                v=voltage[ft.find_time_index(time, end) :, current_max_freq],
                i=np.zeros_like(time[ft.find_time_index(time, end) :]),
                start=end,
                end=time[-1],
                filter=self.filterHighCutFreq,
            )
            EphysObject_start = efex.EphysSweepFeatureExtractor(
                t=time[: ft.find_time_index(time, start) + 1],
                v=voltage[: ft.find_time_index(time, start) + 1, current_max_freq],
                i=np.zeros_like(time[: ft.find_time_index(time, start)]),
                start=0,
                end=start,
                filter=self.filterHighCutFreq,
            )

            EphysObject_end.process_spikes()
            EphysObject_start.process_spikes()
            if EphysObject_end._spikes_df.size or EphysObject_start._spikes_df.size:
                artifact_occurence = True
            while artifact_occurence:
                current_max_freq -= 1

                EphysObject_end = efex.EphysSweepFeatureExtractor(
                    t=time[ft.find_time_index(time, end) :],
                    v=voltage[ft.find_time_index(time, end) :, current_max_freq],
                    i=np.zeros_like(time[ft.find_time_index(time, end) :]),
                    start=end,
                    end=time[-1],
                    filter=self.filterHighCutFreq,
                )
                EphysObject_start = efex.EphysSweepFeatureExtractor(
                    t=time[: ft.find_time_index(time, start) + 1],
                    v=voltage[: ft.find_time_index(time, start) + 1, current_max_freq],
                    i=np.zeros_like(time[: ft.find_time_index(time, start)]),
                    start=0,
                    end=start,
                    filter=self.filterHighCutFreq,
                )
                EphysObject_end.process_spikes()
                EphysObject_start.process_spikes()
                if (
                    not EphysObject_end._spikes_df.size
                    and not EphysObject_start._spikes_df.size
                ):
                    artifact_occurence = False

            # Adding wildness: the feature that for neurogliaform cells specifically can describe whether for highest firing
            # traces the cell sometimes shows APs before and/or after the current stimulation window.
            wildness = (
                df_related_features.iloc[current_max_freq_initial, :].loc["spike_count"]
                - df_related_features.iloc[current_max_freq, :].loc["spike_count"]
            )

            # Adding spike frequency adaptation (ratio of spike frequency of second half to first half for the highest
            # frequency count trace)
            if (
                df_related_features.iloc[current_max_freq, :].loc["spike_count"] < 5
            ):  # If less than 5 spikes we choose not
                # to calculate the SFA ==> np.nan
                SFA = np.nan
            else:
                SFA = df_related_features.iloc[current_max_freq, :].loc["SFA"]
            # Note: we are trying to make sure that if SFA is 0, that it is actually 0 in the way it is defined

            # Adding the Fano factor, a measure of the dispersion of a probability distribution (std^2/mean of the isis)
            # Adding the coefficient of variation. Time intervals between Poisson events should follow an exponential distribution
            # for which the cv should be 1. So if the neuron fires like a Poisson process a cv = 1 should capture that.
            fano_factor = df_related_features.iloc[current_max_freq, :].loc[
                "fano_factor"
            ]
            cv = df_related_features.iloc[current_max_freq, :].loc["cv"]
            AP_fano_factor = df_related_features.iloc[current_max_freq, :].loc[
                "AP_fano_factor"
            ]
            AP_cv = df_related_features.iloc[current_max_freq, :].loc["AP_cv"]

            # Do we have non-Nan values for the burstiness feature?
            non_nan_indexes_BI = ~np.isnan(df_related_features["burstiness"].values)
            if non_nan_indexes_BI.any():

                # Consider only the first and first 5 after threshold reached if possible
                # np.sum will consider a True as a 1 here and a False as a 0 (so you count the True's effectively)

                if np.sum(non_nan_indexes_BI) >= 5:
                    BI_array = df_related_features["burstiness"].values[
                        non_nan_indexes_BI
                    ][0:5]
                    burstiness = np.median(BI_array)
                else:  # Take everything you have
                    BI_array = df_related_features["burstiness"].values[
                        non_nan_indexes_BI
                    ]
                    burstiness = np.median(BI_array)
                if burstiness < 0:
                    burstiness = 0
            else:
                burstiness = 0

            # Do we have non-Nan values for the adaptation index (i.e. traces with more than 1 spike)? We can use this to
            # calculate AP amplitude changes too
            non_nan_indexes_AI = ~np.isnan(df_related_features["isi_adapt"].values)

            if non_nan_indexes_AI.any():

                # Consider only the first 5 after threshold reached if possible
                # np.sum will consider a True as a 1 here and a False as a 0 (so you count the True's effectively)

                if np.sum(non_nan_indexes_AI) >= 5:
                    ISI_adapt_array = df_related_features["isi_adapt"].values[
                        non_nan_indexes_AI
                    ][0:5]
                    ISI_adapt = np.median(ISI_adapt_array)
                    ISI_adapt_average_array = df_related_features[
                        "isi_adapt_average"
                    ].values[non_nan_indexes_AI][0:5]
                    ISI_adapt_average = np.median(ISI_adapt_average_array)
                    AP_amp_adapt_array = df_related_features["AP_amp_adapt"].values[
                        non_nan_indexes_AI
                    ][0:5]
                    AP_amp_adapt = np.median(AP_amp_adapt_array)
                    AP_amp_adapt_average_array = df_related_features[
                        "AP_amp_adapt_average"
                    ].values[non_nan_indexes_AI][0:5]
                    AP_amp_adapt_average = np.median(AP_amp_adapt_average_array)

                else:  # Take everything you have
                    ISI_adapt_array = df_related_features["isi_adapt"].values[
                        non_nan_indexes_AI
                    ]
                    ISI_adapt = np.median(ISI_adapt_array)
                    ISI_adapt_average_array = df_related_features[
                        "isi_adapt_average"
                    ].values[non_nan_indexes_AI]
                    ISI_adapt_average = np.median(ISI_adapt_average_array)
                    AP_amp_adapt_array = df_related_features["AP_amp_adapt"].values[
                        non_nan_indexes_AI
                    ]
                    AP_amp_adapt = np.median(AP_amp_adapt_array)
                    AP_amp_adapt_average_array = df_related_features[
                        "AP_amp_adapt_average"
                    ].values[non_nan_indexes_AI]
                    AP_amp_adapt_average = np.median(AP_amp_adapt_average_array)

            else:
                ISI_adapt = np.nan
                ISI_adapt_average = np.nan
                AP_amp_adapt = np.nan
                AP_amp_adapt_average = np.nan

            # We calculate the latency: the time it takes to elicit the first AP
            df_latency = df_related_features[df_related_features["spike_count"] > 0]

            non_nan_indexes_latency = ~np.isnan(df_latency["latency"].values)
            # Only the first AP is considered for the first trace for which the current clamp stimulation is higher than the
            # current hold
            latency = df_latency["latency"].values[non_nan_indexes_latency][0] * 1000
            try:
                latency_2 = (
                    df_latency["latency"].values[non_nan_indexes_latency][1] * 1000
                )
            except:
                latency_2 = np.nan
            # pdb.set_trace()
            # First index where there is an AP and the current stimulation magnitude is positive
            try:
                if (
                    len(
                        np.where(
                            df.loc[0]["fast_trough_i"].values[
                                ~np.isnan(df.loc[0]["fast_trough_i"].values)
                            ]
                            >= 0
                        )[0]
                    )
                    > 0
                ):
                    index_df = np.where(
                        df.loc[0]["fast_trough_i"].values[
                            ~np.isnan(df.loc[0]["fast_trough_i"].values)
                        ]
                        >= 0
                    )[0][0]

                    non_nan_indexes_ahp = ~np.isnan(df.loc[0]["fast_trough_v"])
                    if non_nan_indexes_ahp.any():
                        AHP = (
                            df.loc[0]["fast_trough_v"].values[index_df]
                            - df.loc[0]["threshold_v"].values[index_df]
                        )
                        # Only first AP is considered
                    else:
                        AHP = 0

                    # ADP (alculated w.r.t. AHP)
                    non_nan_indexes_adp = ~np.isnan(df.loc[0]["adp_v"])
                    if non_nan_indexes_adp.any():
                        ADP = (
                            df.loc[0]["adp_v"].values[index_df]
                            - df.loc[0]["fast_trough_v"].values[index_df]
                        )
                        # Only first AP is considered
                    else:
                        ADP = 0
                    non_nan_indexes_thresh = ~np.isnan(df.loc[0]["threshold_v"])
                    if non_nan_indexes_thresh.any():
                        if df.loc[0]["threshold_v"].size > 1:
                            AP_threshold = df.loc[0]["threshold_v"].values[
                                index_df
                            ]  # Only first AP is considered
                            AP_amplitude = df.loc[0]["peak_height"].values[
                                index_df
                            ]  # Only first AP is considered
                            AP_width = (
                                1000 * df.loc[0]["width"].values[index_df]
                            )  # Only first AP is considered
                            UDR = df.loc[0]["upstroke_downstroke_ratio"].values[
                                index_df
                            ]  # Only first AP is considered
                        else:
                            AP_threshold = df.loc[0][
                                "threshold_v"
                            ]  # Only first AP is considered
                            AP_amplitude = df.loc[0][
                                "peak_height"
                            ]  # Only first AP is considered
                            AP_width = (
                                1000 * df.loc[0]["width"]
                            )  # Only first AP is considered
                            UDR = df.loc[0][
                                "upstroke_downstroke_ratio"
                            ]  # Only first AP is considered
                    else:
                        AP_threshold = 0
                        AP_amplitude = 0
                        AP_width = 0
                else:
                    AP_threshold = 0
                    AP_amplitude = 0
                    AP_width = 0
                    ADP = 0
                    AHP = 0
                    UDR = 0

            except:
                AP_threshold = 0
                AP_amplitude = 0
                AP_width = 0
                ADP = 0
                AHP = 0
                UDR = 0

            # We estimate the rheobase based on a few (i.e. 5)
            # suprathreshold currents steps. A linear fit of the spike frequency w.r.t. to the current injection values of
            # these steps should give the rheobase as the crossing with the x-axis. (Method almost in agreement with Alexandra
            # Naka et al.: "Complementary networks of cortical somatostatin interneurons enforce layer specific control.",
            # they additionally use the subthreshold current step closest to the first suprathreshold one. We think that biases
            # the regression analysis somewhat). We take the min of the first current step showing spikes and the regression line crossing
            # with the x-axis, unless the crossing is at an intercept lower than the last current step still showing no spikes (in
            # that case the first current step showing spikes is simply taken as value for the rheoabse). This method should
            # approximate the rheobase well provided the stimulus interval is long (in our case 600 ms).
            # If a regression cannot be performed, then also simply take the first current step for which spikes have been observed
            # (i.e. not the subthreshold current step!)

            # Only positive currents for this experimental paradigm (i.e. 'spikes' observed for negative currents should not
            # be there)
            # pdb.set_trace()
            # df_rheobase = df_related_features[['current', 'spike_count']][df_related_features['current'] >= 0]
            df_rheobase = df_related_features[["current", "spike_count"]][
                df_related_features["spike_count"] >= 0
            ]
            if len(np.nonzero(df_rheobase["spike_count"].values)[0]) > 4:
                indices = np.nonzero(df_rheobase["spike_count"].values)[0][:5]
                counts = [
                    list(df_rheobase["spike_count"].values[indices]).count(x)
                    for x in df_rheobase["spike_count"].values[indices]
                ]
                if np.max(np.array(counts)) < 3:
                    ransac.fit(
                        df_rheobase["current"].values[indices].reshape(-1, 1),
                        df_rheobase["spike_count"].values[indices].reshape(-1, 1)
                        / (end - start),
                    )
                    line_X = np.concatenate(
                        (df_rheobase["current"].values[indices], np.array([0]))
                    ).reshape(-1, 1)
                    slope = ransac.estimator_.coef_[0][0]
                    sub_thresh_curr = df_rheobase["current"].values[
                        np.nonzero(df_rheobase["spike_count"].values)[0][0] - 1
                    ]  # Last current step with no observed spikes
                    first_supra_thresh_curr = df_rheobase["current"].values[
                        np.nonzero(df_rheobase["spike_count"].values)[0][0]
                    ]  # First current step with observed spikes
                    rheobase = (
                        -ransac.predict(np.array([0]).reshape(-1, 1))[0][0] / slope
                    )
                    if rheobase < sub_thresh_curr:
                        rheobase = first_supra_thresh_curr  # Take the first current step for which spikes are observed
                    elif rheobase > first_supra_thresh_curr:
                        rheobase = first_supra_thresh_curr  # Take the first current step for which spikes are observed
                else:  # A subset can probably not be found by RANSAC to do a regression
                    # Just take the first current step for which spikes have been observed
                    rheobase = df_rheobase["current"].values[
                        np.nonzero(df_rheobase["spike_count"].values)[0][0]
                    ]

            else:  # Not enough datapoints to calculate a rheobase
                # Just take the first current step for which spikes have been observed

                rheobase = df_rheobase["current"].values[
                    np.nonzero(df_rheobase["spike_count"].values)[0][0]
                ]

        else:
            AHP = 0
            ADP = 0
            max_freq = 0
            rebound_spikes = 0
            SFA = np.nan
            fano_factor = np.nan
            cv = np.nan
            norm_sq_isis = np.nan
            ISI_adapt = np.nan
            ISI_adapt_average = np.nan
            AP_amp_adapt = np.nan
            AP_amp_adapt_average = np.nan
            AP_fano_factor = np.nan
            AP_cv = np.nan
            AP_threshold = 0
            AP_amplitude = 0
            AP_width = 0
            UDR = 0
            rheobase = 0
            latency = 0
            latency_2 = np.nan
            burstiness = 0
            wildness = 0
            rebound_spikes = 0

        if np.isnan(ADP):
            ADP = 0

        if rebound < 0:
            rebound = 0

        name_features = [
            "Holding MP (mV)",
            "Fitted MP (mV)",
            "Input resistance (MOhm)",
            "Membrane time constant (ms)",
            "AP threshold (mV)",
            "AP amplitude (mV)",
            "AP width (ms)",
            "Upstroke-to-downstroke ratio",
            "Afterhyperpolarization (mV)",
            "Afterdepolarization (mV)",
            "ISI adaptation index",
            "Max number of APs",
            "Rheobase (pA)",
            "Sag ratio",
            "Latency (ms)",
            "Latency @ +20pA current (ms)",
            "Spike frequency adaptation",
            "ISI Fano factor",
            "ISI coefficient of variation",
            "ISI average adaptation index",
            "Rebound (mV)",
            "Sag time (s)",
            "Sag area (mV*s)",
            "AP amplitude adaptation index",
            "AP amplitude average adaptation index",
            "AP Fano factor",
            "AP coefficient of variation",
            "Burstiness",
            "Wildness",
            "Rebound number of APs",
        ]
        features = [
            Rm,
            v_rest,
            Ri,
            tau,
            AP_threshold,
            AP_amplitude,
            AP_width,
            UDR,
            AHP,
            ADP,
            ISI_adapt,
            max_freq,
            rheobase,
            sag_ratio,
            latency,
            latency_2,
            SFA,
            fano_factor,
            cv,
            ISI_adapt_average,
            rebound,
            sag_time,
            sag_area,
            AP_amp_adapt,
            AP_amp_adapt_average,
            AP_fano_factor,
            AP_cv,
            burstiness,
            wildness,
            rebound_spikes,
        ]
        cell_features = dict(zip(name_features, features))
        Cell_Features = pd.DataFrame([cell_features])
        Cell_Features = Cell_Features.reindex(columns=name_features)
        return Cell_Features

    def alignSpikes(self, feature2Use="peak_index"):
        """get all detected spikes and align with certain feature.

        -------------------
        Parameters:
        feature2Use: should amoung the availabe spike feature:
            ['threshold', 'peak', 'trough', 'upstroke', 'downstroke']
        -------------------
        Return:
        dictionary: with sweep index as key. value is ND array of aligned spikes, time x number of speek
        """
        if "sweepCount" not in self.df.columns:
            print("No spikes detected")
            return
        spikeDF = self.df
        # sweep_feature = self.df_related_features
        feature2Use = "peak_index"
        self.spikeList = {}  ## if not already exist
        self.spikeList_dvdt = {}  ## if not already exist
        filt_coeff = (10000.0) / (
            self.sampleRate / 2.0
        )  # filter kHz -> Hz, then get fraction of Nyquist frequency
        if filt_coeff < 0 or filt_coeff >= 1:
            # raise ValueError("bessel coeff ({:f}) is outside of valid range [0,1); cannot filter sampling frequency {:.1f} kHz with cutoff frequency {:.1f} kHz.".format(filt_coeff, self.sampleRate / 1e3, filter))
            print("bessel coeff  is outside of valid range [0,1)")
            print("sammpling rate: ", self.sampleRate / 1e3, " kHz")
            filt_coeff = 0.95
        b, a = signal.bessel(4, filt_coeff, "low")
        v_filt = signal.filtfilt(b, a, self.data, axis=0)
        dvdt_filt = np.diff(v_filt, axis=0) * self.sampleRate
        for _, spikes_feature in spikeDF.groupby("sweepCount"):
            sweepCount = spikes_feature.sweepCount.iloc[0]

            peaks = spikes_feature[feature2Use].to_list()  ## lis tof peak index
            winBefore = 0.002  ## To do! couple it with paramter tree!
            winAfter = 0.006
            before_ = int(winBefore * self.sampleRate)  # minISI//2
            after_ = int(winAfter * self.sampleRate)  # minISI - minISI//20
            v = np.zeros((after_ + before_, len(peaks)))
            dvdt = np.zeros((after_ + before_, len(peaks)))

            for j, p in enumerate(peaks):

                start_index = p - before_  # before_
                end_index = p + after_  # after_
                # print(start_index, end_index, sweepCount)
                # print(start_index, end_index, sweepCount)
                v[:, j] = v_filt[start_index:end_index, sweepCount]
                dvdt[:, j] = dvdt_filt[start_index:end_index, sweepCount]
            self.spikeList[str(sweepCount)] = v
            self.spikeList_dvdt[str(sweepCount)] = dvdt

    def updateSpikesforCurrentSweep(self, sweepNumber, plotWidget):
        if not self.spikeList:
            self.alignSpikes()
        plotWidget.clear()  ## refresh the layout
        if str(sweepNumber) not in list(self.spikeList.keys()):
            return

        #        plt.setLabels(bottom=('Time', traces[0].XUnit), left=(traces[0].YUnit))
        data = self.spikeList[str(sweepNumber)]
        dvdt_filter = self.spikeList_dvdt[str(sweepNumber)]
        plt = plotWidget.addPlot()
        plt.showGrid(x=True, y=True)
        #        colors = sns.light_palette("red", data.shape[1])
        time = np.arange(data.shape[0]) / self.sampleRate * 1000  ## milisecond
        for j in range(data.shape[1]):

            myPen = pg.mkPen(
                color=pg.intColor(j, hues=data.shape[1])
            )  ## the pen to draw this
            # dvdt = ft.calculate_dvdt(data[:,j], time, self.filter)
            # dvdt = np.diff(data[:,j])*self.sampleRate/1000
            plt.plot(
                data[:, j],
                dvdt_filter[:, j],
                pen=myPen,
                name="Sweep" + str(sweepNumber),
            )
        plt.setTitle(
            "Sweep " + str(sweepNumber + 1) + " Spike count:" + str(data.shape[1])
        )
        plt.setLabel("bottom", "Voltage", units="mV")
        plt.setLabel("left", "dtdt", units="mV/mS")
        return plt

    def plot_alignedSpikes(self, nCol=4, axis=None, mplWidget=None):

        if not self.spikeList:
            self.alignSpikes()
        nSpikes = len(self.spikeList)
        if nSpikes > 0:
            nRow = np.int(np.ceil(nSpikes / nCol))
            if mplWidget == None:  ## pure matplotlib
                fig, axes = plt.subplots(nRow, nCol, figsize=(12, 12))
                for c, spk in enumerate(list(self.spikeList.items())[:16]):
                    axes[c // nCol, np.mod(c, nCol)].plot(spk[1])
                    axes[c // nCol, np.mod(c, nCol)].title.set_text(
                        "Sweep"
                        + str(np.int(spk[0]) + 1)
                        + " n="
                        + str(spk[1].shape[1]),
                        fontsize=10,
                        fontweight="bold",
                    )
                    axes[c // nCol, np.mod(c, nCol)].axis("off")
                    c = c + 1
            else:  ## pyqtgraph MatplotlibWidget.MatplotlibWidget objects defined in PatchViewer 3
                mplWidget.show()
                mplWidget.clf()
                #                axes = mplWidget.subplots(nRow, nCol)
                gs = mplWidget.figure.add_gridspec(4, 4)
                axes = []
                for c, spk in enumerate(list(self.spikeList.items())[:16]):
                    axes.append(
                        mplWidget.figure.add_subplot(gs[c // nCol, np.mod(c, nCol)])
                    )
                    axes[c].plot(spk[1])
                    title = f"sw{np.int(spk[0])+1} {self.current[np.int(spk[0])]:.0f} pA, n={spk[1].shape[1]}"
                    axes[c].title.set_text(title)
                    axes[c].axis("off")
                #                mplWidget.draw()
                #                    print(c)
                mplWidget.figure.tight_layout()

    def plot_nTh_spike(self):
        """plot shape of the nth (up to the 16th) spike across sweep"""
        n = len(self.spikeList)
        nRow = np.int(np.ceil(n / 4))
        fig, ax = plt.subplots(nRow, 4, figsize=(12, 12))
        c = 0
        for c in range(16):
            for key, value in self.spikeList.items():
                if value.shape[1] > c:
                    ax[c // 4, np.mod(c, 4)].plot(value[:, c])
                    ax[c // 4, np.mod(c, 4)].title.set_text(str(c) + " th spike ")
                ax[c // 4, np.mod(c, 4)].axis("off")

    def getInputResistance(self):

        ivs = self.df_related_features[["current", "voltage_deflection"]].loc[
            (self.df_related_features["spike_count"] == 0)
        ]

        if len(ivs) > 6:
            ivs = ivs.iloc[-6:]
        try:
            ransac.fit(
                ivs["current"].values.reshape(-1, 1),
                ivs["voltage_deflection"].values.reshape(-1, 1),
            )
            line_X = ivs["current"].values.reshape(-1, 1)
            v0 = float(
                ransac.estimator_.intercept_
            )  ## membrane potention at zero current level
            prv = ransac.predict(line_X)
            slope = ransac.estimator_.coef_[0][0] * 1000  ## unit of Mohm
        except:
            slope = np.nan
            prv = np.nan
            v0 = np.nan
            line_X = np.nan

        return slope, line_X, prv, v0

    def plot_InputResistance(self, ax=None, font_size=None):
        ivs = self.df_related_features[["current", "voltage_deflection"]].loc[
            (self.df_related_features["spike_count"] == 0)
        ]
        if len(ivs) > 6:
            ivs = ivs.iloc[-6:]
        if font_size == None:
            font_size = 12

        slope, lineX, prv, v_rest = self.getInputResistance()

        ax.plot(
            ivs["current"].values.reshape(-1, 1),
            ivs["voltage_deflection"].values.reshape(-1, 1),
            ".k",
            markersize=12,
        )
        #        ax.plot(ivs['current'].values.reshape(-1, 1), \
        #                ivs['voltage_deflection'].values.reshape(-1, 1)/(end-start), 'b')

        ax.plot(lineX, prv, "--", color="r", linewidth=1.5)
        ax.set_xlabel("Current (pA)", fontsize=font_size)
        ax.set_ylabel("membrane potential (mV)", fontsize=font_size)
        ax.set_title("I vs V below spiking threshold", fontsize=font_size)
        if ~np.isnan(slope):

            ax.annotate(
                "RI="
                + "{0:.1f}".format(slope)
                + "M"
                + "\u03A9"
                + " Vrest="
                + "{0:.1f}".format(v_rest),
                xy=(lineX[1], prv[0] + 5),
                xycoords="data",
                xytext=(-15, 2),
                textcoords="offset points",
                color="b",
                fontsize=14,
            )
        else:
            print("slope fit is not availabe!")

    def plot_currentSpikecount(self, ax=None, font_size=None):
        end = self.stim_end
        start = self.stim_start
        spikes = self.df_related_features[["current", "spike_count"]][
            self.df_related_features["current"] >= 0
        ]
        if font_size == None:
            font_size = 12
        if len(np.nonzero(spikes["spike_count"].values)[0]) > 0:
            indices = np.nonzero(spikes["spike_count"].values)[0]

            ax.plot(
                spikes["current"].values[indices].reshape(-1, 1),
                spikes["spike_count"].values[indices].reshape(-1, 1) / (end - start),
                ".k",
                markersize=5,
            )
            ax.plot(
                spikes["current"].values[indices].reshape(-1, 1),
                spikes["spike_count"].values[indices].reshape(-1, 1) / (end - start),
                "r",
                linewidth=2.0,
            )
            ax.set_xlabel("Current (pA)", fontsize=font_size)
            ax.set_ylabel("Spike frequency (Hz)", fontsize=font_size)
            ax.set_title("Current vs Spike frequency", fontsize=font_size)

    def plot_currentSpikecount_rheobase(self, ax=None, font_size=None):
        end = self.stim_end
        start = self.stim_start
        spikes = self.df_related_features[["current", "spike_count"]][
            self.df_related_features["current"] >= 0
        ]
        if font_size == None:
            font_size = 12
        if len(np.nonzero(spikes["spike_count"].values)[0]) > 0:
            indices = np.nonzero(spikes["spike_count"].values)[0][:4]

            ax.plot(
                spikes["current"].values[indices].reshape(-1, 1),
                spikes["spike_count"].values[indices].reshape(-1, 1) / (end - start),
                ".k",
                markersize=10,
            )
            #                ax.plot(spikes['current'].values[indices].reshape(-1, 1), \
            #                        spikes['spike_count'].values[indices].reshape(-1, 1)/(end-start), '--r')
            ax.set_xlabel("Current (pA)", fontsize=font_size)
            ax.set_ylabel("Spike frequency (Hz)", fontsize=font_size)
            ax.set_title("Rhoeobase", fontsize=font_size)

        df_rheobase = spikes
        if len(np.nonzero(df_rheobase["spike_count"].values)[0]) > 4:

            indices = np.nonzero(df_rheobase["spike_count"].values)[0][:4]
            counts = [
                list(df_rheobase["spike_count"].values[indices]).count(x)
                for x in df_rheobase["spike_count"].values[indices]
            ]

            if np.max(np.array(counts)) < 3:
                ransac.fit(
                    df_rheobase["current"].values[indices].reshape(-1, 1),
                    df_rheobase["spike_count"].values[indices].reshape(-1, 1)
                    / (end - start),
                )
                indices_copy = indices.copy()
                if indices[0] > 0:
                    indices_copy = np.insert(indices_copy, 0, indices[0] - 1)
                line_X = np.concatenate(
                    (df_rheobase["current"].values[indices_copy], np.array([0]))
                ).reshape(-1, 1)
                slope = ransac.estimator_.coef_[0][0]
                sub_thresh_curr = df_rheobase["current"].values[
                    np.nonzero(df_rheobase["spike_count"].values)[0][0] - 1
                ]  # Last current step with no observed spikes
                first_supra_thresh_curr = df_rheobase["current"].values[
                    np.nonzero(df_rheobase["spike_count"].values)[0][0]
                ]  # First current step with observed spikes
                rheobase = -ransac.predict(np.array([0]).reshape(-1, 1))[0][0] / slope
                if rheobase < sub_thresh_curr:
                    rheobase = first_supra_thresh_curr  # Take the first current step for which spikes are observed
                elif rheobase > first_supra_thresh_curr:
                    rheobase = first_supra_thresh_curr  # Take the first current step for which spikes are observed
                prv = ransac.predict(line_X)
                ax.plot(line_X, prv, "--", color="r", linewidth=1.5)
                ax.set_title("Rheobase estimation", fontsize=font_size)
                ax.spines["top"].set_visible(False)
                ax.tick_params(axis="both", which="major", labelsize=font_size)
                ax.set_ylim(
                    [
                        -2.5,
                        np.max(
                            df_rheobase["spike_count"].values[indices].reshape(-1, 1)
                            / (end - start)
                        )
                        + 3,
                    ]
                )
                ax.set_xlim(
                    [
                        np.min(
                            df_rheobase["current"].values[indices_copy].reshape(-1, 1)
                        )
                        - self.current_step,
                        np.max(
                            df_rheobase["current"].values[indices_copy].reshape(-1, 1)
                        )
                        + self.current_step,
                    ]
                )

                # Let's denote the membrane voltage clearly on the plot
                ax.plot(rheobase, 0, marker="|", color="black", ms=50)
                ax.annotate(
                    "",
                    xy=(rheobase - 15, 2),
                    xycoords="data",
                    xytext=(rheobase, 2),
                    textcoords="data",
                    arrowprops={
                        "arrowstyle": "<-",
                        "connectionstyle": "arc3",
                        "lw": 2,
                        "ec": "grey",
                        "shrinkB": 0,
                    },
                )
                ax.annotate(
                    f"Rheobase:{rheobase:.1f} pA",
                    xy=(rheobase, prv[1]),
                    xycoords="data",
                    xytext=(-15, 2),
                    textcoords="offset points",
                    color="b",
                    fontsize=12,
                )

    def plot_rhobose(self, ax=None):
        end = self.stim_end
        start = self.stim_start
        df_rheobase = self.df_related_features[["current", "spike_count"]][
            self.df_related_features["current"] >= 0
        ]
        if len(np.nonzero(df_rheobase["spike_count"].values)[0]) > 4:
            indices = np.nonzero(df_rheobase["spike_count"].values)[0][:5]

            counts = [
                list(df_rheobase["spike_count"].values[indices]).count(x)
                for x in df_rheobase["spike_count"].values[indices]
            ]
            if np.max(np.array(counts)) < 3:
                ransac.fit(
                    df_rheobase["current"].values[indices].reshape(-1, 1),
                    df_rheobase["spike_count"].values[indices].reshape(-1, 1)
                    / (end - start),
                )
                if indices[0] > 0:
                    indices_copy = indices.copy()
                    indices_copy = np.insert(indices_copy, 0, indices[0] - 1)
                line_X = np.concatenate(
                    (df_rheobase["current"].values[indices_copy], np.array([0]))
                ).reshape(-1, 1)
                slope = ransac.estimator_.coef_[0][0]
                sub_thresh_curr = df_rheobase["current"].values[
                    np.nonzero(df_rheobase["spike_count"].values)[0][0] - 1
                ]  # Last current step with no observed spikes
                first_supra_thresh_curr = df_rheobase["current"].values[
                    np.nonzero(df_rheobase["spike_count"].values)[0][0]
                ]  # First current step with observed spikes
                rheobase = -ransac.predict(np.array([0]).reshape(-1, 1))[0][0] / slope
                if rheobase < sub_thresh_curr:
                    rheobase = first_supra_thresh_curr  # Take the first current step for which spikes are observed
                elif rheobase > first_supra_thresh_curr:
                    rheobase = first_supra_thresh_curr  # Take the first current step for which spikes are observed

                ax.plot(
                    df_rheobase["current"].values[indices].reshape(-1, 1),
                    df_rheobase["spike_count"].values[indices].reshape(-1, 1)
                    / (end - start),
                    ".k",
                    markersize=15,
                )
                ax.plot(line_X, ransac.predict(line_X), color="k")
                ax.set_xlabel("Current (pA)", fontsize=12)
                ax.set_ylabel("Spike frequency (Hz)", fontsize=12)
                ax.set_title("Rheobase estimation", fontsize=12)
                ax.spines["top"].set_visible(False)
                ax.tick_params(axis="both", which="major", labelsize=16)
                ax.set_ylim(
                    [
                        -2.5,
                        np.max(
                            df_rheobase["spike_count"].values[indices].reshape(-1, 1)
                            / (end - start)
                        )
                        + 3,
                    ]
                )
                ax.set_xlim(
                    [
                        np.min(
                            df_rheobase["current"].values[indices_copy].reshape(-1, 1)
                        )
                        - 10,
                        np.max(
                            df_rheobase["current"].values[indices_copy].reshape(-1, 1)
                        )
                        + 10,
                    ]
                )

                # Let's denote the membrane voltage clearly on the plot
                ax.plot(rheobase, 0, marker="|", color="black", ms=50)
                ax.annotate(
                    "",
                    xy=(rheobase - 15, 2),
                    xycoords="data",
                    xytext=(rheobase, 2),
                    textcoords="data",
                    arrowprops={
                        "arrowstyle": "<-",
                        "connectionstyle": "arc3",
                        "lw": 2,
                        "ec": "grey",
                        "shrinkB": 0,
                    },
                )
                ax.annotate(
                    f"Rheobase:{rheobase:.1f} pA",
                    xy=(rheobase - 15, 2),
                    xycoords="data",
                    xytext=(-15, 2),
                    textcoords="offset points",
                    color="k",
                    fontsize=14,
                )

    def plot_info(self, el_num=2, axis=None, font_size=None):
        """Analyses a specific data dictionary corresponding to a cell and returns a figure object with annotations on a particular
        trace of how certain features have been calculated.
        Parameters
        ----------
        el_num : integer, from which electrode number has been measured (optional, 2 by default)
        Returns
        -------
        ax : figure object
        """
        current = self.current
        time = self.time
        start = self.stim_start
        end = self.stim_end
        voltage = self.data
        df = self.df
        df_related_features = self.df_related_features
        curr_index_0 = self.current0_index
        df = self.df
        if font_size == None:
            font_size = font_size
        if not np.any(df_related_features["spike_count"].values > 3):
            return
        # First current magnitude (starting from the lowest absolute current stimulation magnitude value) where we find more than
        # three spikes
        current_first = np.flatnonzero(df_related_features["spike_count"].values >= 3)[
            0
        ]
        if axis:
            ax = axis
        else:
            f, ax = plt.subplots(figsize=(10, 10))

        if self.time[-1] < 0.9:
            ax.plot(time, self.voltage[:, current_first], color=np.array([0, 0, 0]))
        else:
            ax.plot(
                time[: ft.find_time_index(time, 0.9)],
                voltage[: ft.find_time_index(time, 0.9), current_first],
                color=np.array([0, 0, 0]),
            )
        # Actual current there
        current_first_magn = self.current[current_first]
        # Amount of spikes there
        spike_count = df_related_features["spike_count"].values[current_first]

        # Find start and end indexes in df for all the spikes in that particular train
        index_start_df = np.flatnonzero(df["threshold_i"].values >= current_first_magn)[
            0
        ]
        index_end_df = (
            np.flatnonzero(df["threshold_i"].values >= current_first_magn)[0]
            + spike_count
        )
        # When last spike is clipped
        if df["clipped"].values[index_end_df - 1]:
            index_end_df = np.flatnonzero(
                df["threshold_i"].values >= current_first_magn
            )[0] + (spike_count - 1)
        # When the spike before the last spike is also clipped
        if df["clipped"].values[index_end_df - 1]:
            index_end_df = np.flatnonzero(
                df["threshold_i"].values >= current_first_magn
            )[0] + (spike_count - 2)

        # Get threshold, peak, upstroke and downstroke indexes
        thresh_index = np.array(
            df["threshold_index"].values[index_start_df:index_end_df], dtype=int
        )
        upstroke_index = np.array(
            df["upstroke_index"].values[index_start_df:index_end_df], dtype=int
        )
        peak_index = np.array(
            df["peak_index"].values[index_start_df:index_end_df], dtype=int
        )
        downstroke_index = np.array(
            df["downstroke_index"].values[index_start_df:index_end_df], dtype=int
        )
        fast_trough_index = np.array(
            df["fast_trough_index"].values[index_start_df:index_end_df], dtype=int
        )
        fast_trough_index = np.array(
            fast_trough_index[~np.isnan(fast_trough_index)], dtype=int
        )
        slow_trough_index = df["slow_trough_index"].values[index_start_df:index_end_df]
        slow_trough_index = np.array(
            slow_trough_index[~np.isnan(slow_trough_index)], dtype=int
        )

        start_index = (
            np.abs(time - start / 10)
        ).argmin()  # Find closest index where the injection current starts (quite a bit before now)
        end_index = (
            np.abs(time - end)
        ).argmin()  # Find closest index where the injection current ends

        # ax.plot(time[np.abs(time - start/2).argmin()], Cell_Features['Rm (mV)'], 'C0.', ms = 15, label = None)
        ax.plot(
            time[thresh_index],
            voltage[thresh_index, current_first],
            "b.",
            ms=10,
            label="AP threshold",
        )
        ax.plot(
            time[upstroke_index],
            voltage[upstroke_index, current_first],
            "r.",
            ms=10,
            label="AP upstroke",
        )
        ax.plot(
            time[peak_index],
            voltage[peak_index, current_first],
            "g.",
            ms=10,
            label="AP peak",
        )
        ax.plot(
            time[downstroke_index],
            voltage[downstroke_index, current_first],
            "k.",
            ms=10,
            label="AP downstroke",
        )
        ax.plot(
            time[fast_trough_index],
            voltage[fast_trough_index, current_first],
            "y.",
            ms=10,
            label="AP fast trough",
        )
        ax.plot(
            time[slow_trough_index],
            voltage[slow_trough_index, current_first],
            "m.",
            ms=10,
            label="AP slow trough\n(if applicable)",
        )
        if time[thresh_index[-1]] + 0.1 > 0.9:
            endtime = 0.9
        else:
            endtime = time[thresh_index[-1]] + 0.1
        ax.set_xlim(time[thresh_index[0]] - 0.05, endtime)
        ax.set_xlabel("Time (s)", fontsize=font_size)
        ax.set_ylabel("Membrane voltage (mV)", fontsize=font_size)
        ax.set_title("First trace showing at least three APs", fontsize=font_size)
        ax.tick_params(axis="both", which="major", labelsize=font_size)
        # ax.legend(['Resting Vm', 'Threshold', 'Upstroke', 'Peak', 'Downstroke', 'Fast trough', \
        #           'Slow trough (if applicable)'], fontsize = 15, loc='upper left', bbox_to_anchor=(1, 1))
        ax.legend(fontsize=font_size, loc="upper right", bbox_to_anchor=(1.1, 1.0))

        # Nice annotations
        if len(thresh_index) > 2:
            ax.plot(
                time[thresh_index[0:3]],
                voltage[thresh_index[0:3], current_first],
                "|",
                color="black",
                ms=200,
            )
            ax.annotate(
                "",
                xy=(
                    time[thresh_index[0]],
                    voltage[thresh_index[0], current_first] - 10,
                ),
                xycoords="data",
                xytext=(
                    time[thresh_index[1]],
                    voltage[thresh_index[0], current_first] - 10,
                ),
                textcoords="data",
                arrowprops={
                    "arrowstyle": "<->",
                    "connectionstyle": "arc3",
                    "lw": 2,
                    "ec": "grey",
                    "shrinkA": 0,
                },
            )
            ax.annotate(
                "",
                xy=(
                    time[thresh_index[1]],
                    voltage[thresh_index[0], current_first] - 10,
                ),
                xycoords="data",
                xytext=(
                    time[thresh_index[2]],
                    voltage[thresh_index[0], current_first] - 10,
                ),
                textcoords="data",
                arrowprops={
                    "arrowstyle": "<->",
                    "connectionstyle": "arc3",
                    "lw": 2,
                    "ec": "grey",
                    "shrinkA": 0,
                },
            )
            ax.annotate(
                "ISI adapt. index = 2nd ISI / 1st ISI",
                xy=(
                    time[thresh_index[1]],
                    voltage[thresh_index[1], current_first] - 30,
                ),
                xycoords="data",
                xytext=(5, 0),
                textcoords="offset points",
                fontsize=font_size,
            )
        return ax

    def plot_info_first_peak(self, el_num=2, axis=None, font_size=None):
        """Analyses a specific data dictionary corresponding to a cell and returns a figure object with annotations on the first
        peak of a particular trace of how certain features have been calculated. Works for data extracted from .mat files, not for
        .asc files.
        Parameters
        ----------
        data : data full of voltage (V) and time (s) for a particular cell
        el_num : integer, from which electrode number has been measured (optional, 2 by default)
        current_step : float, which current step (pA) has been used between consecutive experiments (optional, 20 by default)
        start : start of the stimulation (s) in the voltage trace (optional, default 0.1)
        end : end of the stimulation (s) in the voltage trace (optional, default 0.7)
        axis : axis you'd like to plot information on (optional, None by default)

        Returns
        -------
        ax : figure object
        """

        current = self.current
        time = self.time
        voltage = self.data
        df = self.df
        df_related_features = self.df_related_features

        # First current magnitude (starting from the lowest absolute current stimulation magnitude value) where we find more than
        # three spikes
        if font_size == None:
            font_size = 12
        current_first = np.flatnonzero(
            np.logical_and(
                df_related_features["spike_count"].values >= 1,
                df_related_features["current"].values > 0,
            )
        )[0]

        # Actual current there
        current_first_magn = current[current_first]
        # Amount of spikes there
        spike_count = df_related_features["spike_count"].values[current_first]

        # Find start and end indexes in df for all the spikes in that particular train
        index_start_df = np.flatnonzero(df["threshold_i"].values >= current_first_magn)[
            0
        ]

        # Get threshold, peak, upstroke and downstroke indexes
        thresh_index = np.array(df["threshold_index"].values[index_start_df], dtype=int)
        upstroke_index = np.array(
            df["upstroke_index"].values[index_start_df], dtype=int
        )
        peak_index = np.array(df["peak_index"].values[index_start_df], dtype=int)
        downstroke_index = np.array(
            df["downstroke_index"].values[index_start_df], dtype=int
        )
        fast_trough_index = np.array(
            df["fast_trough_index"].values[index_start_df], dtype=int
        )
        slow_trough_index = np.array(
            df["slow_trough_index"].values[index_start_df], dtype=int
        )
        adp_index = np.array(df["adp_index"].values[index_start_df], dtype=int)
        # slow_trough_index = np.array(slow_trough_index[~np.isnan(slow_trough_index)], dtype = int)

        start_index = thresh_index - 50
        if slow_trough_index.size & (adp_index > 0):
            end_index = adp_index + 50  # Plot 50 time indices after the adp index
        else:
            end_index = (
                fast_trough_index + 50
            )  # Plot 50 time indices after the fast trough index (should always exist)

        # When first spike is clipped
        if df["clipped"].values[index_start_df]:
            return
        if axis:
            ax = axis
        else:
            f, ax = plt.subplots(figsize=(10, 5))

        #        sns.set_context(rc={'lines.markeredgewidth': 1})

        ax.plot(
            time[start_index:end_index],
            voltage[start_index:end_index, current_first],
            color=np.array([0, 0, 0]),
            label=None,
        )

        ax.plot(
            time[thresh_index],
            voltage[thresh_index, current_first],
            "b.",
            ms=12,
            label="AP threshold",
        )
        ax.plot(
            time[upstroke_index],
            voltage[upstroke_index, current_first],
            "r.",
            ms=12,
            label="AP upstroke",
        )
        ax.plot(
            time[peak_index],
            voltage[peak_index, current_first],
            "g.",
            ms=12,
            label="AP peak",
        )
        ax.plot(
            time[downstroke_index],
            voltage[downstroke_index, current_first],
            "k.",
            ms=12,
            label="AP downstroke",
        )
        ax.plot(
            time[fast_trough_index],
            voltage[fast_trough_index, current_first],
            "y.",
            ms=12,
            label="AP fast trough",
        )
        if slow_trough_index.size & (adp_index > 0):
            ax.plot(
                time[adp_index],
                voltage[adp_index, current_first],
                "c.",
                ms=12,
                label="ADP",
            )
            # ax.plot(time[slow_trough_index], voltage[slow_trough_index, current_first], 'm.', ms = 15, label = \
            #       'AP slow trough\n(if applicable)')
        #        ax.legend(fontsize = 10, loc = 'upper right')
        ax.set_xlabel("Time (s)", fontsize=font_size)
        ax.set_ylabel("Membrane voltage (mV)", fontsize=font_size)
        ax.tick_params(axis="both", which="major", labelsize=font_size)
        ax.set_title("First AP", fontsize=font_size)

        # Nice annotations

        # For the AP amplitude
        ax.annotate(
            "",
            xy=(time[peak_index], voltage[peak_index, current_first]),
            xycoords="data",
            xytext=(time[peak_index], voltage[thresh_index, current_first]),
            textcoords="data",
            arrowprops={
                "arrowstyle": "<->",
                "ec": "grey",
                "connectionstyle": "arc3",
                "lw": 2,
                "shrinkB": 0,
            },
        )

        ax.plot(
            time[peak_index],
            voltage[thresh_index, current_first],
            marker="_",
            color="black",
            ms=100,
        )
        ax.plot(
            time[peak_index],
            voltage[peak_index, current_first],
            marker="_",
            color="black",
            ms=100,
        )

        # For the AP width
        width_level = (
            voltage[peak_index, current_first] - voltage[thresh_index, current_first]
        ) / 2 + voltage[thresh_index, current_first]
        width_start_index = (
            peak_index
            - np.flatnonzero(
                voltage[peak_index:thresh_index:-1, current_first] <= width_level
            )[0]
        )
        width_end_index = (
            peak_index
            + np.flatnonzero(
                voltage[peak_index:fast_trough_index, current_first] <= width_level
            )[0]
        )
        ax.plot(
            time[width_start_index],
            voltage[width_end_index, current_first],
            "|",
            color="black",
            ms=100,
        )
        ax.plot(
            time[width_end_index],
            voltage[width_end_index, current_first],
            "|",
            color="black",
            ms=100,
        )

        # The width itself is calculated based on t[width_end_index] - t[width_start_index], but the voltages might be different
        # at the respective indices, thus we choose the arrow to go from v[width_end_index] to v[width_end_index] to make
        # it horizontal (interpretability of the figure improves)
        ax.annotate(
            "",
            xy=(time[width_start_index], voltage[width_end_index, current_first]),
            xycoords="data",
            xytext=(time[width_end_index], voltage[width_end_index, current_first]),
            textcoords="data",
            arrowprops={
                "arrowstyle": "<->",
                "connectionstyle": "arc3",
                "lw": 2,
                "ec": "grey",
                "shrinkA": 0,
            },
        )
        ax.annotate(
            "AP width",
            xy=(time[width_start_index], width_level - 5),
            xycoords="data",
            xytext=(5, -15),
            textcoords="offset points",
            fontsize=font_size,
        )

        # We still need to annotate the AP amplitude based on the width_level!
        ax.annotate(
            "AP amplitude",
            xy=(time[peak_index], width_level + 30),
            xycoords="data",
            xytext=(5, 0),
            textcoords="offset points",
            fontsize=font_size,
        )

        # For the AHP
        ax.plot(
            time[fast_trough_index],
            voltage[thresh_index, current_first],
            marker="_",
            color="black",
            ms=100,
        )
        ax.plot(
            time[fast_trough_index],
            voltage[fast_trough_index, current_first],
            marker="_",
            color="black",
            ms=100,
        )
        ax.annotate(
            "",
            xy=(time[fast_trough_index], voltage[thresh_index, current_first]),
            xycoords="data",
            xytext=(time[fast_trough_index], voltage[fast_trough_index, current_first]),
            textcoords="data",
            arrowprops={
                "arrowstyle": "<->",
                "connectionstyle": "arc3",
                "lw": 2,
                "ec": "grey",
                "shrinkB": 0,
            },
        )
        fast_trough_level = (
            voltage[thresh_index, current_first]
            - voltage[fast_trough_index, current_first]
        ) / 2 + voltage[fast_trough_index, current_first]
        ax.annotate(
            "AHP",
            xy=(time[fast_trough_index], fast_trough_level),
            xycoords="data",
            xytext=(10, -5),
            textcoords="offset points",
            fontsize=font_size,
        )

        # For a possible ADP
        if slow_trough_index.size & (adp_index > 0):
            ax.plot(
                time[adp_index],
                voltage[adp_index, current_first],
                marker="_",
                color="black",
                ms=100,
            )
            ax.plot(
                time[adp_index],
                voltage[fast_trough_index, current_first],
                marker="_",
                color="black",
                ms=100,
            )
            ax.annotate(
                "",
                xy=(time[adp_index], voltage[fast_trough_index, current_first]),
                xycoords="data",
                xytext=(time[adp_index], voltage[adp_index, current_first]),
                textcoords="data",
                arrowprops={
                    "arrowstyle": "<->",
                    "connectionstyle": "arc3",
                    "lw": 2,
                    "ec": "b",
                    "shrinkB": 0,
                },
            )
            adp_level = (
                voltage[adp_index, current_first]
                - voltage[fast_trough_index, current_first]
            ) / 2 + voltage[fast_trough_index, current_first]
            ax.annotate(
                "ADP",
                xy=(time[adp_index], adp_level),
                xycoords="data",
                xytext=(10, -5),
                textcoords="offset points",
                fontsize=font_size,
            )

        return ax

    def plot_max_spikes_trace(self, el_num=2, axis=None):
        """Analyses a specific data dictionary corresponding to a cell and returns a figure object with the trace for which
        the features ISI FF, ISI CV, AP amp FF, AP amp CV, max # number of spikes in 600 ms and the Spike frequency adaptation
        are normally calculated (only in the very specific case of more traces showing the same max # of spikes in 600 ms can
        the trace for which the feature is being calculated differ among the features). Works for data extracted from .mat files,
        not for .asc files.

        Parameters
        ----------
        data : data full of voltage (V) and time (s) for a particular cell
        el_num : integer, from which electrode number has been measured (optional, 2 by default)
        current_step : float, which current step (pA) has been used between consecutive experiments (optional, 20 by default)
        start : start of the stimulation (s) in the voltage trace (optional, default 0.1)
        end : end of the stimulation (s) in the voltage trace (optional, default 0.7)
        axis : axis you'd like to plot information on (optional, None by default)

        Returns
        -------
        ax : figure object
        """

        time = self.time
        start = self.stim_start
        end = self.stim_end
        voltage = self.data
        df_related_features = self.df_related_features

        # The max amount of spikes in 600 ms of the trace showing the max amount of spikes in 600 ms
        max_freq = np.max(df_related_features["spike_count"].values)
        # Take the first trace showing this many spikes if there are many
        current_max_freq = np.flatnonzero(
            df_related_features["spike_count"].values >= max_freq
        )[0]
        # Check if there are peaks outside the stimilation interval for the highest firing trace. When true, continue looking
        # for lower firing traces untill it shows none anymore.
        artifact_occurence = False
        EphysObject_end = efex.EphysSweepFeatureExtractor(
            t=time[ft.find_time_index(time, end) :],
            v=voltage[ft.find_time_index(time, end) :, current_max_freq],
            i=np.zeros_like(time[ft.find_time_index(time, end) :]),
            start=end,
            end=time[-1],
        )
        EphysObject_start = efex.EphysSweepFeatureExtractor(
            t=time[: ft.find_time_index(time, start)],
            v=voltage[: ft.find_time_index(time, start), current_max_freq],
            i=np.zeros_like(time[: ft.find_time_index(time, start)]),
            start=0,
            end=start - (time[1] - time[0]),
        )
        EphysObject_end.process_spikes()
        EphysObject_start.process_spikes()
        if EphysObject_end._spikes_df.size or EphysObject_start._spikes_df.size:
            artifact_occurence = True
        while artifact_occurence:
            current_max_freq -= 1

            EphysObject_end = efex.EphysSweepFeatureExtractor(
                t=time[ft.find_time_index(time, end) :],
                v=voltage[ft.find_time_index(time, end) :, current_max_freq],
                i=np.zeros_like(time[ft.find_time_index(time, end) :]),
                start=end,
                end=time[-1],
            )
            EphysObject_start = efex.EphysSweepFeatureExtractor(
                t=time[: ft.find_time_index(time, start)],
                v=voltage[: ft.find_time_index(time, start), current_max_freq],
                i=np.zeros_like(time[: ft.find_time_index(time, start)]),
                start=0,
                end=start - (time[1] - time[0]),
            )
            EphysObject_end.process_spikes()
            EphysObject_start.process_spikes()
            if (
                not EphysObject_end._spikes_df.size
                and not EphysObject_start._spikes_df.size
            ):
                artifact_occurence = False
        if axis:
            ax = axis
        else:
            f, ax = plt.subplots(figsize=(15, 3))
        sns.set_context(rc={"lines.markeredgewidth": 1})

        if time[-1] < 0.9:
            ax.plot(time, voltage[:, current_max_freq], color="black")
        else:
            ax.plot(
                time[: ft.find_time_index(time, 0.9)],
                voltage[: ft.find_time_index(time, 0.9), current_max_freq],
                color="black",
            )
        ax.set_xlabel("Time (s)", fontsize=font_size)
        ax.set_ylabel("Membrane voltage (mV)", fontsize=font_size)
        ax.tick_params(axis="both", which="major", labelsize=font_size)
        ax.set_title("Highest frequency trace", fontsize=font_size)
        return ax

    def plot_lowest_trace(self, el_num=2, axis=None, font_size=None):
        """Analyses a specific data dictionary corresponding to a cell and returns a figure object with annotations on a particular
        trace of how certain features have been calculated.
        Parameters
        ----------
        data : data full of voltage (V) and time (s) for a particular cell
        el_num : integer, from which electrode number has been measured (optional, 2 by default)
        current_step : float, which current step (pA) has been used between consecutive experiments (optional, 20 by default)
        start : start of the stimulation (s) in the voltage trace (optional, default 0.1)
        end : end of the stimulation (s) in the voltage trace (optional, default 0.7)
        fil : cutoff frequency for 4-pole low-pass Bessel filter in kHz (optional, default 10)
        axis : axis you'd like to plot information on (optional, None by default)

        Returns
        -------
        ax : figure object
        """

        from matplotlib.patches import Polygon

        time = self.time
        start = self.stim_start
        end = self.stim_end
        voltage = self.data
        for j in range(
            10
        ):  ## among the first 10 traces, find the first one that has negative deflection voltage
            if not np.isnan(self.df_related_features.iloc[j].tau):
                ltrace = j
                break
        if font_size == None:
            font_size = 12
        if axis:
            ax = axis
        else:
            f, ax = plt.subplots(figsize=(10, 10))
        if time[-1] < 0.9:
            ax.plot(time, voltage[:, ltrace], color=np.array([0, 0, 0]), label=None)
        else:
            ax.plot(
                time[: ft.find_time_index(time, 0.9)],
                voltage[: ft.find_time_index(time, 0.9), ltrace],
                color=np.array([0, 0, 0]),
                label=None,
            )
        ax.set_xlabel("Time (s)", fontsize=font_size)
        ax.set_ylabel("Membrane voltage (mV)", fontsize=font_size)
        ax.set_title("Lowest trace", fontsize=font_size)
        ax.tick_params(axis="both", which="major", labelsize=font_size)

        baseline_interval = 0.1  # To calculate the SS voltage

        # v_baseline = EphysObject._get_baseline_voltage()
        v_baseline = ft.average_voltage(
            voltage[:, ltrace], time, start=start - 0.1, end=start
        )

        start_index = ft.find_time_index(time, start)
        end_index = ft.find_time_index(time, end)
        index_rebound = end_index + 0.05  # Just a first initial value

        if (
            np.flatnonzero(voltage[end_index:, ltrace] > v_baseline).size == 0
        ):  # So perfectly zero here means
            # it did not reach it
            rebound = 0
        else:
            index_rebound = (
                end_index + np.flatnonzero(voltage[end_index:, ltrace] > v_baseline)[0]
            )
            if not (
                time[index_rebound] > (end + 0.2)
            ):  # We definitely have 100 ms left to calculate the rebound
                rebound = (
                    ft.average_voltage(
                        voltage[
                            index_rebound : index_rebound
                            + ft.find_time_index(time, 0.15),
                            ltrace,
                        ],
                        time[
                            index_rebound : index_rebound
                            + ft.find_time_index(time, 0.15)
                        ],
                    )
                    - v_baseline
                )
            else:  # Work with whatever time is left
                rebound = (
                    ft.average_voltage(
                        voltage[index_rebound:, ltrace], time[index_rebound:]
                    )
                    - v_baseline
                )

        # v_peak, peak_index = EphysObject.voltage_deflection("min")

        peak_index = start_index + np.argmin(voltage[start_index:end_index, ltrace])
        ax.plot(
            time[peak_index],
            voltage[peak_index, ltrace],
            ".",
            c=[0, 0, 0],
            markersize=15,
            label="Sag trough",
        )
        v_steady = ft.average_voltage(
            voltage[:, ltrace], time, start=end - baseline_interval, end=end
        )

        # First time SS is reached after stimulus onset
        first_index = (
            start_index
            + np.flatnonzero(voltage[start_index:peak_index, ltrace] < v_steady)[0]
        )
        # First time SS is reached after the max voltage deflection downwards in the sag
        if np.flatnonzero(voltage[peak_index:end_index, ltrace] > v_steady).size == 0:
            second_index = end_index
        else:
            second_index = (
                peak_index
                + np.flatnonzero(voltage[peak_index:end_index, ltrace] > v_steady)[0]
            )
        # Time_to_SS is the time difference between these two time points
        time_to_SS = time[second_index] - time[first_index]
        # Integration_to_SS is the integration area of the voltage between these two time points
        integration_to_SS = -integrate.cumtrapz(
            voltage[first_index:second_index, ltrace], time[first_index:second_index]
        )[-1]

        # Now let's add nice annotations
        # First up the time and integration to reach the SS
        ax.plot(
            time[first_index],
            voltage[first_index, ltrace],
            "|",
            markersize=30,
            color=[0, 0, 0],
        )
        ax.plot(
            time[second_index],
            voltage[second_index, ltrace],
            "|",
            markersize=30,
            color=[0, 0, 0],
        )
        ax.annotate(
            "",
            xy=(time[first_index], voltage[first_index, ltrace]),
            xycoords="data",
            xytext=(time[second_index], voltage[second_index, ltrace]),
            textcoords="data",
            arrowprops={
                "arrowstyle": "<->",
                "connectionstyle": "arc3",
                "lw": 2,
                "ec": "grey",
                "shrinkA": 0,
            },
        )
        ax.annotate(
            "sag time",
            xy=(time[first_index] + time_to_SS / 2, v_steady),
            xycoords="data",
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=font_size,
        )
        a, b = time[first_index], time[second_index]
        verts = [
            (a, voltage[first_index, ltrace]),
            *zip(
                time[first_index:second_index],
                voltage[first_index:second_index, ltrace],
            ),
            (b, voltage[second_index, ltrace]),
        ]
        poly = Polygon(verts, facecolor="0.9", edgecolor="0.5")
        ax.add_patch(poly)

        # Now the rebound
        if rebound != 0:
            end_index_for_rebound = index_rebound + ft.find_time_index(time, 0.15)
            if time[index_rebound] > (0.9 - 0.15):
                end_index_for_rebound = ft.find_time_index(
                    time, np.max(time)
                )  # Plot till the end (simply what you have left)
            ax.plot(
                time[index_rebound],
                voltage[index_rebound, ltrace],
                "|",
                markersize=10,
                color=[0, 0, 0],
            )
            ax.plot(
                time[end_index_for_rebound],
                voltage[end_index_for_rebound, ltrace],
                "|",
                markersize=10,
                color=[0, 0, 0],
            )
            if time[index_rebound] == time[end_index_for_rebound]:
                return ax
            ax.plot(
                [time[index_rebound], time[end_index_for_rebound]],
                [
                    ft.average_voltage(
                        voltage[index_rebound:end_index_for_rebound, ltrace],
                        time[index_rebound:end_index_for_rebound],
                    ),
                    ft.average_voltage(
                        voltage[index_rebound:end_index_for_rebound, ltrace],
                        time[index_rebound:end_index_for_rebound],
                    ),
                ],
                "-",
                color=[0, 0, 0],
            )
            ax.plot(
                [time[index_rebound], time[end_index_for_rebound]],
                [v_baseline, v_baseline],
                "-",
                color=[0, 0, 0],
            )
            ax.annotate(
                "",
                xy=(
                    time[index_rebound]
                    + (time[end_index_for_rebound] - time[index_rebound]) / 2,
                    v_baseline,
                ),
                xycoords="data",
                xytext=(
                    time[index_rebound]
                    + (time[end_index_for_rebound] - time[index_rebound]) / 2,
                    v_baseline + rebound,
                ),
                textcoords="data",
                arrowprops={
                    "arrowstyle": "<->",
                    "connectionstyle": "arc3",
                    "lw": 2,
                    "ec": "grey",
                    "shrinkB": 0,
                },
            )
            ax.annotate(
                "rebound",
                xy=(
                    time[index_rebound]
                    + (time[end_index_for_rebound] - time[index_rebound]) / 2,
                    v_baseline + rebound / 2,
                ),
                xycoords="data",
                xytext=(0, 2),
                textcoords="offset points",
                fontsize=font_size,
            )
        #        plt.legend(fontsize = font_size)
        return ax

    def plot_firstLast(self, axis=None, font_size=None):
        ## maybe add more annotations
        time = self.time
        start = self.stim_start
        end = self.stim_end
        voltage = self.data

        if axis:
            ax = axis
        else:
            f, ax = plt.subplots(figsize=(8, 8))
        if font_size == None:
            font_size = 12
        ax.plot(time, voltage[:, 0], color="m", label=None)
        ax.plot(time, voltage[:, 5], color="b", label=None)
        ax.plot(time, voltage[:, 10], color="k", label=None)
        ax.plot(time, voltage[:, -1], color="g", label=None)
        ax.plot(time, voltage[:, -10], color="r", label=None)
        ax.set_xlabel("Time (s)", fontsize=font_size)
        ax.set_ylabel("Membrane voltage (mV)", fontsize=font_size)
        ax.set_title(self.title, fontsize=font_size)
        ax.tick_params(axis="both", which="major", labelsize=font_size)

    def plot_w_style(self, el_num=2, filename="fig.pdf", axes=None, mplWidget=None):
        """Analyses a specific data dictionary corresponding to a cell and returns axes with traces and annotations.
        Works for data extracted from .mat files, not for.asc files.
        Parameters
        ----------
        data : data full of voltage (V) and time (s) for a particular cell
        el_num : integer, from which electrode number has been measured (optional, 2 by default)
        start : start of the stimulation (s) in the voltage trace (optional, default 0.1)
        end : end of the stimulation (s) in the voltage trace (optional, default 0.7)
        current_step : float, which current step (pA) has been used between consecutive experiments (optional, 20 by default)

        Returns
        -------
        axes : figure objects

        """

        if mplWidget:
            mplWidget.show()
            mplWidget.clf()
            #            axes = mplWidget.subplots(2,1)
            gs = mplWidget.figure.add_gridspec(2, 2)
            #            self.plot_firstLast(axis = axes[0])
            axes1 = mplWidget.figure.add_subplot(gs[0, 0])
            self.plot_InputResistance(ax=axes1, font_size=8)

            axes2 = mplWidget.figure.add_subplot(gs[0, 1])
            self.plot_currentSpikecount(
                ax=axes2, font_size=8
            )  #            self.plot_info_first_peak( el_num = el_num,  axis = axes[1], font_size = 8)

            axes3 = mplWidget.figure.add_subplot(gs[1, 0])
            self.plot_currentSpikecount_rheobase(ax=axes3, font_size=8)
            #            self.plot_max_spikes_trace(axis = axes[4])
            axes4 = mplWidget.figure.add_subplot(gs[1, 1])
            self.plot_lowest_trace(axis=axes4, font_size=8)
            mplWidget.draw()
            try:
                mplWidget.figure.tight_layout()
            except Exception as e:
                pass

        else:
            sns.set_style("ticks")
            sns.set_context("paper", font_scale=0.85)
            fig = plt.figure(figsize=(12, 10))
            ax1 = plt.subplot2grid((3, 2), (0, 0))
            ax2 = plt.subplot2grid((3, 2), (0, 1))
            ax3 = plt.subplot2grid((3, 2), (1, 0))
            ax4 = plt.subplot2grid((3, 2), (1, 1))
            ax5 = plt.subplot2grid((3, 2), (2, 0))
            ax6 = plt.subplot2grid((3, 2), (2, 1))
            axes = [ax1, ax2, ax3, ax4, ax5, ax6]

            self.plot_firstLast(axis=axes[0])
            self.plot_info(el_num=el_num, axis=axes[1])
            self.plot_info_first_peak(el_num=el_num, axis=axes[2])
            self.plot_currentSpikecount_rheobase(ax=axes[3])
            self.plot_max_spikes_trace(axis=axes[4])
            self.plot_lowest_trace(axis=axes[5])

            plt.tight_layout()
            for axis in axes:
                sns.despine(ax=axis)
            sns.set_context(rc={"lines.markeredgewidth": 3})
            plt.text(
                -0.05,
                1.04,
                "A",
                transform=axes[0].transAxes,
                fontsize=16,
                fontweight="bold",
            )
            plt.text(
                -0.05,
                1.04,
                "B",
                transform=axes[1].transAxes,
                fontsize=16,
                fontweight="bold",
            )
            plt.text(
                -0.05,
                1.04,
                "C",
                transform=axes[2].transAxes,
                fontsize=16,
                fontweight="bold",
            )
            plt.text(
                -0.05,
                1.04,
                "D",
                transform=axes[3].transAxes,
                fontsize=16,
                fontweight="bold",
            )
            plt.text(
                -0.05,
                1.04,
                "E",
                transform=axes[4].transAxes,
                fontsize=16,
                fontweight="bold",
            )
            plt.text(
                -0.05,
                1.04,
                "F",
                transform=axes[5].transAxes,
                fontsize=16,
                fontweight="bold",
            )
            plt.text(
                -0.05,
                1.04,
                "F",
                transform=plt.gcf().get_axes()[5].transAxes,
                fontsize=16,
                fontweight="bold",
            )
            # fig.savefig(filename)
        return axes


def extract_sweep_feature(time, voltage, curr, start, end, EphysObject):

    spikes_df = pd.DataFrame()
    sweep_df = pd.DataFrame()
    start_index = ft.find_time_index(
        time, start
    )  # Find closest index where the injection current starts
    end_index = ft.find_time_index(
        time, end
    )  # Find closest index where the injection current ends
    # Adding peak_height (mV) + code for maximum frequency determination (see further)
    spike_count = 0
    if EphysObject._spikes_df.size:
        EphysObject._spikes_df["peak_height"] = (
            EphysObject._spikes_df["peak_v"].values
            - EphysObject._spikes_df["threshold_v"].values
        )
        spike_count = EphysObject._spikes_df["threshold_i"].values.size

    spikes_df = pd.concat([spikes_df, EphysObject._spikes_df], sort=True)

    # Some easily found extra features
    sweep_df = EphysObject._sweep_features

    # Adding spike count
    sweep_df.update({"spike_count": spike_count})

    # Adding spike frequency adaptation (ratio of spike frequency of second half to first half)
    SFA = np.nan
    half_stim_index = ft.find_time_index(time, np.float(start + (end - start) / 2))
    if (
        spike_count > 5
    ):  # We only consider traces with more than 8.333 Hz = 5/600 ms spikes here
        # but in the end we only take the trace with the max amount of spikes

        if (
            np.sum(
                spikes_df.loc[spikes_df["threshold_i"] == curr, :]["threshold_index"]
                < half_stim_index
            )
            != 0
        ):
            SFA = np.sum(
                spikes_df.loc[spikes_df["threshold_i"] == curr, :]["threshold_index"]
                > half_stim_index
            ) / np.sum(
                spikes_df.loc[spikes_df["threshold_i"] == curr, :]["threshold_index"]
                < half_stim_index
            )

    sweep_df.update({"SFA": SFA})

    # Adding current (pA)
    sweep_df.update({"current": curr})

    # Adding membrane voltage (mV)
    sweep_df.update({"resting_membrane_potential": EphysObject._get_baseline_voltage()})

    # Adding voltage deflection to steady state (mV)
    # voltage_deflection_SS = ft.average_voltage(voltage[:, c], time, start = end - 0.1, end = end)
    # sweep_df.update({'voltage_deflection': voltage_deflection_SS}) ## this won't work if tau is very large

    voltage_deflection_SS, voltage_deflection_i = EphysObject.voltage_deflection()  #
    sweep_df.update({"voltage_deflection": voltage_deflection_SS})

    # Adding input resistance (MOhm)
    input_resistance = np.nan
    if (
        not ("peak_i" in EphysObject._spikes_df.keys()) and not curr == 0
    ):  # We only calculate input resistances
        #                    try:
        #                        input_resistance, xx, yy, v0 = self.getInputResistance()  ## get the slope of a linear fit between I-V
        #                    except:
        #                        input_resistance = np.nan                                                           # from traces without APs
        input_resistance = (
            np.abs(voltage_deflection_SS - EphysObject._get_baseline_voltage()) * 1000
        ) / np.abs(curr)
        if input_resistance == np.inf:
            input_resistance = np.nan
    sweep_df.update({"input_resistance": input_resistance})

    # Adding membrane time constant (s) and voltage plateau level for hyperpolarisation paradigms
    # after stimulus onset
    tau = np.nan
    E_plat = np.nan
    sag_ratio = np.nan
    if (
        voltage_deflection_SS - sweep_df["resting_membrane_potential"] < -5
        or curr < -19
    ):  # We use hyperpolarising steps as required in the object function to estimate the
        # membrane time constant and E_plateau
        while True:
            try:
                tau = EphysObject.estimate_time_constant()  # Result in seconds!
                break
            except TypeError:  # Probably a noisy bump for this trace, just keep it to be np.nan
                break
        E_plat = ft.average_voltage(voltage, time, start=end - 0.1, end=end)
        sag, sag_ratio = EphysObject.estimate_sag()
    sweep_df.update({"tau": tau})
    sweep_df.update({"E_plat": E_plat})
    sweep_df.update({"sag_ratio": sag_ratio})

    # For the rebound and sag time we only are interested in the lowest (-200 pA (usually)) hyperpolarisation trace
    rebound = np.nan
    sag_time = np.nan
    sag_area = np.nan

    if curr < -19:
        baseline_interval = start  # To calculate the SS voltage
        v_baseline = EphysObject._get_baseline_voltage()

        end_index = ft.find_time_index(time, end)
        if (
            np.flatnonzero(voltage[end_index:] > v_baseline).size == 0
        ):  # So perfectly zero here means
            # it did not reach it
            rebound = 0
        else:
            index_rebound = (
                end_index + np.flatnonzero(voltage[end_index:] > v_baseline)[0]
            )
            if not (
                time[index_rebound] > (end + 0.15)
            ):  # We definitely have 150 ms left to calculate the rebound
                rebound = (
                    ft.average_voltage(
                        voltage[
                            index_rebound : index_rebound
                            + ft.find_time_index(time, 0.15)
                        ],
                        time[
                            index_rebound : index_rebound
                            + ft.find_time_index(time, 0.15)
                        ],
                    )
                    - v_baseline
                )
            else:  # Work with whatever time is left
                if time[-1] == time[index_rebound]:
                    rebound = 0
                else:
                    rebound = (
                        ft.average_voltage(
                            voltage[index_rebound:], time[index_rebound:]
                        )
                        - v_baseline
                    )

        v_peak, peak_index = EphysObject.voltage_deflection("min")
        v_peakMax, peak_indexMax = EphysObject.voltage_deflection("max")
        # v_steady = ft.average_voltage(voltage, time, start = end - baseline_interval, end=end)
        v_steady = ft.average_voltage(
            voltage, time, start=0.005, end=baseline_interval - 0.005
        )
        if v_peakMax > v_steady or v_steady - v_peak < -4:
            # The sag should have a minimum depth of 4 mV
            # otherwise we set sag time and sag area to 0
            sag_time = 0
            sag_area = 0
        else:
            # First time SS is reached after stimulus onset
            try:
                first_index = (
                    start_index
                    + np.flatnonzero(voltage[start_index:peak_index] < v_steady)[0]
                )
                # First time SS is reached after the max voltage deflection downwards in the sag
                if np.flatnonzero(voltage[peak_index:end_index] > v_steady).size == 0:
                    second_index = end_index
                else:
                    second_index = (
                        peak_index
                        + np.flatnonzero(voltage[peak_index:end_index] > v_steady)[0]
                    )
                sag_time = time[second_index] - time[first_index]
                sag_area = -integrate.cumtrapz(
                    voltage[first_index:second_index], time[first_index:second_index]
                )[-1]
            except:
                print("line1569")
                print("sag calcuation error")
                sag_time = 0
                sag_area = 0
                pdb.set_trace()

    burst_metric = np.nan
    # print(c)
    if spike_count > 5:
        burst = EphysObject._process_bursts()
        if len(burst) != 0:
            burst_metric = burst[0][0]

    sweep_df.update({"rebound": rebound})
    sweep_df.update({"sag_time": sag_time})
    sweep_df.update({"sag_area": sag_area})
    sweep_df.update({"burstiness": burst_metric})
    return spikes_df, pd.DataFrame([sweep_df])
