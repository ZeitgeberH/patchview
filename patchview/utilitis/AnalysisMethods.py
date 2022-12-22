# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 10:44:39 2021

@author: M.H.
"""
import numpy as np
from scipy import signal
import yaml
from patchview.HekaIO.HekaHelpers import HekaBundleInfo
from scipy.ndimage import gaussian_filter

def loadYAML(filename):
    """open yaml file for parameters"""
    with open(filename) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
        return data

def saveYAML(filename, data):
    """open yaml file for parameters"""
    with open(filename) as f:
        yaml.dump(data, filename, default_flow_style=True)

def bandPass_signal(data, fs, highCutOff=None, filter_option=0):
    if filter_option == 0:  ## default fourth order Bessel-Thomson filter
        bessel_filter_sos = Bessel_lowpass(fs, highCutOff)
        f = signal.sosfiltfilt(bessel_filter_sos, data)
    else:  ## Butter
        butter_sos = butter_lowpass(fs)
        f = signal.sosfiltfilt(butter_sos, data)
    return f


def Bessel_lowpass(fs, highCutOff=None):
    """

    Lowpass Bessel/Thomson filter

    Parameters
    ----------
    fs : sampling rate
    highCutOff :

    Returns
    -------
    filter paramters

    """

    if highCutOff == None:
        filt_coeff = 0.95
    elif highCutOff > 1500:
        filt_coeff = (highCutOff - 50) / (
            fs / 2.0
        )  # filter kHz -> Hz, then get fraction of Nyquist frequency
    else:
        filt_coeff = highCutOff / (fs / 2.0)
    if filt_coeff < 0 or filt_coeff >= 1:

        print(
            "bessel coeff ({:f}) is outside of valid range [0,1); \
                        cannot filter sampling frequency {:.1f} kHz with \
                        cutoff frequency {:.1f} kHz.".format(
                filt_coeff, fs / 1e3, highCutOff / 1e3
            )
        )
        print("Using Nyqst frequency for high cutoff")
        filt_coeff = 0.95
    bessel_filter_sos = signal.bessel(4, filt_coeff, btype="low", output="sos")
    return bessel_filter_sos


def butter_lowpass(fs, highCutOff=None, order=2):
    """
    Low pass Butter filter.

    Args:
        - fs       (float) : the sampling rate.
        - high_cut (float) : the high cutoff frequency of the filter.
        - order      (int) : order of the filter, by default defined to 5.
    """
    # calculate the Nyquist frequency
    order = 2
    nyq = 0.5 * fs
    # design filter
    if highCutOff != None:
        high = highCutOff / nyq
    else:
        high = 0.8
    if high >= 1:
        raise ValueError(
            "butter filter high cutoff frequency is higher than Nyquist frequency"
        )
    # returns the filter coefficients: numerator and denominator
    butter_sos = signal.butter(order, high, btype="low", output="sos")
    return butter_sos


def filterDatSeries(data, fs, hcutFreq=None, filter_option=0):
    """
        filter 3D array data

    Parameters
    ----------
    data : np.ndarray
        assume time is along axis 0, channel at axis 1, sweep at axis 2.

    fs : float / int
        sampling rate

    hcutFreq : float, optional
        DESCRIPTION. higtcutoff frequency The default is None.

    filter_option : TYPE, optional
        DESCRIPTION. The default is 0. see `bandPass_signal` for more detains

    Returns
    -------
    data_ : filted data with the same shape
    """
    nTraces = data.shape[1]
    nSweeps = data.shape[2]
    data_ = np.zeros_like(data)
    ## loops through all the sweeps
    for sweep in range(nSweeps):
        for t in range(nTraces):
            data_[:, t, sweep] = bandPass_signal(
                data[:, t, sweep], fs, hcutFreq, filter_option
            )
    return data_


def calculateConnectionTraces(time, fs, data, stimStartSample):
    """
    Clculate connections between traces

    Parameters
    ----------
    time : TYPE
        DESCRIPTION.
    fs :  sampling rate

    data : time X channel X sweeps

    stimStartSample :  sample number where stimuli start

    Returns
    -------
    delaySamples : TYPE
        DESCRIPTION.
    responseSamples : TYPE
        DESCRIPTION.
    avg2 : TYPE
        DESCRIPTION.
    peak : TYPE
        DESCRIPTION.
    peak_idx : TYPE
        DESCRIPTION.
    significantTest : TYPE
        DESCRIPTION.
    baselineMean : TYPE
        DESCRIPTION.
    meanCorrelation : TYPE
        DESCRIPTION.

    """
    avg = np.squeeze(np.mean(data, axis=2))
    nTrace = avg.shape[1]
    responseTime = (
        0.040  ## 100 ms for calcuating baseline And to compare after baseline
    )
    responseOnsetDelay = 0.006
    responseSamples = np.int(fs * responseTime)
    delaySamples = np.int(fs * responseOnsetDelay)
    nStd = 4
    avg2 = np.zeros(avg.shape)
    significantTest = []
    peak = []
    peak_idx = []
    for j in range(nTrace):
        data_ = np.squeeze(data[:, j, :])
        # data_ = self.bandPass_signal(data_.transpose(), hcutFreq)
        avg2[:, j] = np.mean(data_, axis=1)
        t0 = stimStartSample - delaySamples
        baselineData = np.mean(data_[:t0, :], axis=1)
        baselineMean = np.mean(baselineData)
        baselineStd = np.std(baselineData)

        ## calculate mean of response across sweeps
        t1 = stimStartSample + delaySamples
        t2 = stimStartSample + delaySamples + responseSamples
        responseData = np.squeeze(np.mean(data_[t1:t2, :], axis=1))
        trialConsistency = np.corrcoef(data_[t1:t2, :].transpose(), rowvar=True)
        selIdx = np.triu_indices(data_.shape[1], 1)
        meanCorrelation = np.mean(trialConsistency[selIdx])

        peak_idx.append(np.argmax(np.abs(responseData - baselineMean)))
        peak.append(responseData[peak_idx[j]])  ## get peak during this window
        if (
            np.abs(peak[j] - baselineMean) >= nStd * baselineStd
            and meanCorrelation >= 0.1
        ):
            significantTest.append(1)
        else:
            significantTest.append(0)
    return (
        delaySamples,
        responseSamples,
        avg2,
        peak,
        peak_idx,
        significantTest,
        baselineMean,
        meanCorrelation,
    )


def cleanASCfile(filename):
    f_1 = open(filename, "r")
    tempName = filename[:-4] + "_mod.ASC"
    f_2 = open(tempName, "w")
    ASCline = f_1.readline()
    skipOn = False
    while ASCline:
        if ASCline == '("NewContour"\n':
            f_2.writelines(ASCline)
            skipOn = True
        if not skipOn:
            f_2.writelines(ASCline)

        if ASCline.strip()[:11] == "(Resolution":
            skipOn = False

            f_2.writelines("(NewContour)\n")
        ASCline = f_1.readline()
    f_1.close()
    f_2.close()
    return tempName

def padding_(ndiff,   bin_array, padVal=0):
    ''' ndiff: total number to pad
        bin_array: original array
        padVal: padded value
    '''
    leftPadd = ndiff//2
    rightPadd = ndiff -leftPadd
    step = np.diff(bin_array)[0]
    if padVal!=0:
        leftArray = np.linspace(bin_array[0]- leftPadd*step, bin_array[0], leftPadd, endpoint=False)
        rightArray = np.linspace(bin_array[-1]+step, bin_array[-1]+step+ rightPadd*step, rightPadd, endpoint=False)
        pad_array = np.concatenate([leftArray,bin_array, rightArray ])
    else:
        pad_array = np.pad(bin_array, (leftPadd, rightPadd),'constant', constant_values=(padVal, padVal))
    return pad_array
    
def smooth2D(data, window_len=10, std=2):
    ''' smooth 2d image with gaussian kernel
    '''
    truncatN = window_len//2
    return gaussian_filter(data, sigma=std, order=0, output=None, mode='constant', cval=0.0, truncate=truncatN)

def smooth(x, window_len=10, std=2, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
    """
    if len(x) < window_len:
        window_len = len(x)//2
    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]] ## padded
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    elif window == 'gaussian':
        w = signal.gaussian(window_len, std=std)
    else:
        w =eval('np.'+window+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='valid')
    pn = window_len//2
    return y[(pn-1):-pn]

if __name__ == "__main__":
    n = 4
    x = np.array([1,5,9])
    print(padding_(n, x, padVal=0))