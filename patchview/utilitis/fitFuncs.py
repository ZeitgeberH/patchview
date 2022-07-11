# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 17:05:45 2021

@author: MHu
"""
import numpy as np
from scipy.fft import fft, ifft
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_prominences
from scipy.signal import bessel, sosfiltfilt
from scipy.ndimage import gaussian_filter


def singleExp(x, a, b, c):
    return a * np.exp(-x / b) + c


def biExp(x, tau_rise, tau_decay, r=1.0, d=1.0, baseline=0.01):
    return -r * np.exp(-x / tau_rise) + d * np.exp(-x / tau_decay) + baseline


def biExp0(x, tau_rise, tau_decay, delay=0.01, r=1.0, d=1.0):
    return -r * np.exp(-(x - delay) / tau_rise) + d * np.exp(-(x - delay) / tau_decay)


## Bandpass filter
def butter_bandpass_filter(data, lowcut, highcut, fs, order=4):
    filt_coeff = [lowcut / (fs / 2), highcut / (fs / 2)]
    wn = bessel(order, filt_coeff, btype="bandpass", output="sos")
    y = sosfiltfilt(wn, data)
    return y


def signal_deconvolution(rTrace, tp, fs, lowF=1, highF=200):
    """
    deconvolving raw data with a given template

    Parameters
    ----------
    rTrace : numpy array
        filtered traces
    tp : numpy array
        template
    fs: sampling rate
    Returns
    -------
    dTrace : numpy array
        detection trace
    Reference(s):
    [1] A. Pernía-Andrade, S.P. Goswami, Y. Stickler, U. Fröbe, A. Schlögl, and P. Jonas (2012)
     A deconvolution-based method with high sensitivity and temporal resolution for  detection of spontaneous synaptic currents in vitro and in vivo.
    Biophysical Journal Volume 103 October 2012 1–11.
    """
    # pdb.set_trace()
    kpad = 2 * int(fs)
    tp1 = np.pad(tp, (0, 2 * kpad), "constant", constant_values=0)  # post padding
    rTrace1 = np.pad(rTrace, (kpad, kpad), "symmetric")  ## symmetric paddding
    tp_f = fft(tp1)
    rTrace_f = fft(rTrace1)
    D = rTrace_f / tp_f  ## division in the frequency domain
    D = ifft(D).real  ## back to time domain and get real part
    D = butter_bandpass_filter(D, lowF, highF, fs, order=4)
    D = D[kpad:-kpad]  ## removing padding
    return D


def event_hist(D, alg="fd"):
    """
    Parameters
    ----------
    D : TYPE
        DESCRIPTION.
    alg : Str, optional
        DESCRIPTION. Algrithom to calculte optimal bin edges:
       {‘auto’,'fd','doane','scott','stone','rice',‘sturges’}
        The default is 'fd'(Freedman Diaconis Estimator)

    Returns
    -------
    histgram, bin_edges

    """
    eventHist, bin_edges = np.histogram(D, bins=alg)


def getRawPeaks2(rawData, D, Dpeaks, peak_threh, peaksign=1, wlen=500):
    ## peakendIdx
    baselineD = np.mean(D)
    idx = []
    selected_index = []
    if peaksign > 0:  ## positive peak
        for sel, p in enumerate(Dpeaks):
            data0 = D[p : p + wlen]
            initW = 1
            baselineIdx = np.where(data0 <= baselineD)[0]
            downStrokeIdx = 0
            while len(baselineIdx) < 1:
                initW -= 0.01
                if initW < 0.9:
                    downStrokeIdx = 50
                    print("Baseline not reached at 50%")
                    break
                baselineIdx = np.where(data0 <= baselineD * initW)[0]
                if len(baselineIdx) < 1:
                    print("Inloop peak index", p, " ", initW)
            if initW >= 0.9:
                if initW != 1:
                    print("Outloop", baselineIdx)
                pidx = p + baselineIdx[0]
            else:
                pidx = p + downStrokeIdx
            if rawData[pidx] >= peak_threh:
                idx.append(pidx)
                selected_index.append(sel)
    else:  # there is overshoot for negative peak
        for sel, p in enumerate(Dpeaks):
            data0 = D[p : p + wlen]
            initW = 1
            baselineIdx = np.where(data0 < baselineD)[0]
            downStrokeIdx = 0
            while len(baselineIdx) < 1:
                initW -= 0.01
                if initW < 0.9:
                    downStrokeIdx = 50

                    break
                baselineIdx = np.where(data0 < baselineD * initW)[0]

            if initW >= 0.9:
                pidx = p + baselineIdx[0]  ##
            else:
                pidx = p + downStrokeIdx
            data0 = D[pidx : pidx + wlen]
            baselineIdx = np.where(data0 >= baselineD)[0]
            if len(baselineIdx) < 1:
                pidx = pidx + 50
                print("no return to baseline")
            else:
                pidx = pidx + baselineIdx[0]
            if rawData[pidx] <= peak_threh:
                idx.append(pidx)
                selected_index.append(sel)
    return idx, selected_index


def getRawPeaks3(rawData, D, Dpeaks, peak_threh, peaksign=1, wlen=500, avgWlen=10):
    idx = []
    selected_index = []
    # pdb.set_trace()
    data_sm = (
        np.convolve(rawData, np.ones(avgWlen), "same") / avgWlen
    )  ## moving average of dvdt for denoising
    dvdt = np.diff(data_sm)
    dvdt = (
        np.convolve(dvdt, np.ones(avgWlen), "same") / avgWlen
    )  ## moving average of dvdt for denoising
    # dvdt2 = np.convolve(np.diff(dvdt), np.ones(avgWlen), 'same') / avgWlen ## moving average of dvdt for denoising
    if peaksign > 0:  ## positive peak
        for sel, p in enumerate(Dpeaks):
            data0 = rawData[p : p + wlen]
            dvdt0 = dvdt[p : p + wlen]
            for j, dvdt_ in enumerate(dvdt0[:-1]):
                if dvdt0[j] > 0 and dvdt0[j + 1] <= 0 and data0[j] >= peak_threh:
                    idx.append(p + j)
                    selected_index.append(sel)
                    break
    else:  # negative peaks
        for sel, p in enumerate(Dpeaks):
            data0 = rawData[p : p + wlen]
            dvdt0 = dvdt[p : p + wlen]
            for j, dvdt_ in enumerate(dvdt0[:-1]):
                if dvdt0[j] < 0 and dvdt0[j + 1] >= 0 and data0[j] <= peak_threh:
                    idx.append(p + j)
                    selected_index.append(sel)
                    break

    return idx, selected_index


def filterPeakOnset(peakonset_hight, th, sign):
    """
    Filtering detected peaks based on onset voltage.
    peak onset height at opposite end (beyong certain therhold) would be filter out.
    """
    peak_index = []
    if sign > 0:
        for idx, ph in enumerate(peakonset_hight):
            if ph <= th:
                continue
            else:
                peak_index.append(idx)

    else:
        for idx, ph in enumerate(peakonset_hight):
            if ph >= th:
                continue
            else:
                peak_index.append(idx)
    return peak_index


def getRawPeaks(data, Dpeaks, wlen=500):
    """
    The raw data peaks are lagging behind the deconvovled one.
    We can looking for the first local maximal after the deconvoved peak. which should be the actual peak!

    Parameters
    ----------
    data : raw data
    Dpeaks : index for peaks from deconvovled traces
    wlen: number of samples to look ahead
    Returns
    -------
    index for raw peaks

    """
    N = len(data) - 1
    idx = []
    for p in Dpeaks:
        endIdx = np.min([N, p + wlen])
        data0 = data[p:endIdx]
        # pdb.set_trace()
        data0 = gaussian_filter(data0, sigma=5)
        data_ = np.abs(np.diff(np.sign(np.diff(data0))))
        try:
            pidx = p + np.where(data_ > 1)[0][0]
            idx.append(pidx)
        except LookupError:
            print("No peak found!")

    peakHeights = [data[j] for j in idx]
    # pdb.set_trace()
    return idx, peakHeights


def getPeaks(
    D,
    height,
    distance=5,
    width=(2, 50),
    wlen=50,
    prominence=None,
    threshold=None,
    rel_height=1,
):  # wlen=50, width=(2,), distance = 5
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
        if supplied as the maximal required plateau size.:
    """

    peak_idx, peak_prop = find_peaks(
        D,
        height=height,
        wlen=wlen,
        width=width,
        distance=distance,
        rel_height=rel_height,
    )

    # peak_idx, peak_prop = find_peaks(D, height = height, wlen=50, width=(2,), distance = 5)
    peak_heights = D[peak_idx]  ## global peak hight
    # contour_heights = peak_heights - peak_prop['prominences'] ## local peak height
    return peak_idx, peak_heights, peak_prop


curveFitFuncs = {"SingleExponential": singleExp, "BiExponential": biExp}

curveFitPars = {
    "SingleExponential": {"a": 0.0, "b": 0.0, "c": 0.0},
    "BiExponential": {
        "tau1(s)": 0.0,
        "tau2(s)": 0.0,
        "amp1(Normed)": 0.0,
        "amp2(Normed)": 0.0,
        "baseline(Normed)": 0.0,
    },
}
