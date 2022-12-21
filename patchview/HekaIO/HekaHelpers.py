# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 16:33:46 2021

@author: MHu
"""
from patchview.HekaIO import HEKA_Reader_MAIN as HEKA
import os
import numpy as np


class HekaBundleInfo(object):
    """
    A helpr class wrap aound HEKA reader
    """

    def __init__(self, filePath):
        self.bundle = None
        self.readDatFile(filePath)

    def readDatFile(self, filePath):
        """

        filepath : full path of a .dat file
        """
        path, f = os.path.split(filePath)
        # check extension
        _, ext = os.path.splitext(filePath)
        assert ext == ".dat", "input needs to be a .dat file"

        self.bundle = HEKA.Bundle(filePath)
        self.bundle.closeDatFile()

        self.fileName = f
        self.filePath = path

    def countGroups(self):
        return len(self.bundle.pul.children)

    def getGroupRecord(self, idx):
        assert isinstance(idx, list), "group index need to be a list integer"
        ng = self.countGroups()
        if idx[0] >= ng:
            raise ValueError(
                f"Requested index outside availabe number {ng} groups in this bundle"
            )
        else:
            return self.bundle.pul[idx[0]]

    def countSeries(self, idx):
        ## count number of series in a group
        assert isinstance(idx, list), "input needs to be a list"
        for j in idx:
            assert isinstance(j, int), "element of  the list should be all integers"
        assert idx[0] < self.countGroups() and idx[0] >= 0, "group index out of bounds"
        return len(self.bundle.pul[idx[0]].children)

    def getSeriesRecord(self, idx):
        assert len(idx) >= 2, "series index need to at least length of 2"
        assert isinstance(idx[1], int), "Series index need to be a list of two integers"
        return self.bundle.pul[idx[0]][idx[1]]

    def countSweeps(self, idx):
        """
        Count number of sweeps in this series
        Parameters
        ----------
        idx : [int, int]
        Returns
        -------
        int: number of sweeps in this series
        """
        assert isinstance(idx, list), "input needs to be a list"
        for j in idx:
            assert isinstance(j, int), "element of  list should be all integers"
        assert len(idx) >= 2, "sweeps index need to at least length of 2"
        return len(self.bundle.pul[idx[0]][idx[1]].children)

    def getSweepRecord(self, idx):
        assert len(idx) >= 3 and isinstance(
            idx[2], int
        ), "sweep index need to be a list of integers"
        return self.bundle.pul[idx[0]][idx[1]][idx[2]]

    def getTraceRecord(self, idx):
        assert isinstance(idx, list), "trace index needs to be a list"
        assert len(idx) == 4 and isinstance(
            idx[3], int
        ), "trace index need to be a list of 4 integers"
        return self.bundle.pul[idx[0]][idx[1]][idx[2]][idx[3]]

    def countTraces(self, idx):
        assert isinstance(idx, list), "trace index needs to be a list"
        assert len(idx) >= 3 and isinstance(
            idx[2], int
        ), "trace index need to be a list of integers"
        assert idx[0] < self.countGroups(), "group list index out of range"
        assert idx[1] < self.countSeries(idx), "series list index out of range"
        assert idx[2] < self.countSweeps(idx), "sweep list index out of range"
        return len(self.bundle.pul[idx[0]][idx[1]][idx[2]].children)

    def getSeriesSamplingRate(self, idx):
        idx_ = idx.copy()
        assert len(idx_) >= 2 and isinstance(
            idx_[1], int
        ), "input as a list of 2 integers"
        if len(idx_) < 3:
            idx_.extend([0, 0])
        trace = self.getTraceRecord(idx_)
        return 1.0 / trace.XInterval

    def getSeriesLabel(self, idx):
        assert len(idx) >= 2, "input as a list of 2 integers"
        return self.getSeriesRecord(idx).Label

    def getSeriesData(self, idx):
        assert len(idx) >= 2, "series index needs to be at least length of 2"
        idx_ = idx.copy()
        if len(idx_) < 3:
            idx_.extend([0, 0])
        nSweep = self.countSweeps(idx_)  ## number of sweeps for this series
        nTraces = self.countTraces(idx_)
        nSamples = self.getNumberOfSamplesPerSweep(idx_)
        data = np.zeros((nSamples, nTraces, nSweep))  ## ndrarry, time X traces X sweep
        ## loops through all the sweeps
        for sweep in range(nSweep):
            idx_[2] = sweep  ## change sweeep level index
            for t in range(nTraces):
                idx_[3] = t  ## change sweeep level index
                data_ = self.getSingleTraceData(idx_)
                if len(data_) != nSamples:
                    print(
                        "Sweep:", sweep, "Trace:", t + 1, "# of samples: ", len(data_)
                    )
                    data_ = 0  # np.nan
                data[:, t, sweep] = data_

        return data

    def getChanCount(self, idx):
        """count number of channels in current sweep"""

    def getStim(self, idx):
        time, stim, stimInfo = self.bundle.stim(idx)
        return time, stim, stimInfo

    def getSingleTraceData(self, idx):
        # assert isinstance(idx, list) and len(idx) == 4 , 'series index need to be a list with length of 4'
        if (
            self.countTraces(idx) == 1
        ):  ## it happens that experimenter label is not consistent with actual trace count
            print(
                f"single trace series. Sweep label is not consistent with trace index {idx[-1]}"
            )
            idx[-1] = 0
        return self.bundle.data[idx]

    def getTraceUnit(self, idx):
        print("Not implemented!")

    def getNumberOfSamplesPerSweep(self, idx):
        # assert isinstance(idx, list) and len(idx) == 4 , 'series index need to be a list with length of 4'
        return self.bundle.data[idx].shape[0]

    def getSweepTimeStamps(self, idx):
        """
        get time stamps as 1-D array

        Parameters
        ----------
        idx : [int, int, int, int]

        Returns
        -------
        time : np.array((int,))
        """
        idx_ = idx.copy()
        if len(idx) <= 2:
            idx_.extend([0, 0])
        trace = self.getTraceRecord(idx_)
        nSamples = self.getNumberOfSamplesPerSweep(idx_)
        time = np.linspace(
            trace.XStart, trace.XStart + trace.XInterval * (nSamples - 1), nSamples
        )
        return time
