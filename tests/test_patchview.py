# -*- coding: utf-8 -*-
"""Tests for `patchview` package."""

import random
import unittest
import os
import sys
import numpy

def makePath(path, file):
    return os.path.join(path, file)
from patchview.HekaIO import HEKA_Reader_MAIN as HEKA
from patchview.HekaIO.HekaHelpers import HekaBundleInfo

testFile = "210702s5c08.dat" 

testFile_noExit = "180514.dat"
testFile_wrongExt = "180514c1r1.date"

try:
    DataPath= ".\data"
    bundleTester = HekaBundleInfo(makePath(DataPath, testFile))
except Exception:
    print("faild to load test file!")

class Test_HEKA_Reader(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def test_read_hekaFile(self):
        bundle = HekaBundleInfo(makePath(DataPath, testFile))
        self.assertIsInstance(bundle, HekaBundleInfo)

    def test_wrong_ext(self):
        self.assertRaises(
            AssertionError, HekaBundleInfo, makePath(DataPath, testFile_wrongExt)
        )

    def test_file_not_exist(self):
        self.assertRaises(
            FileNotFoundError, HekaBundleInfo, makePath(DataPath, testFile_noExit)
        )

    def test_counts(self):
        ## group count
        self.assertEqual(bundleTester.countGroups(), 1)
        ## series count
        self.assertRaises(AssertionError, bundleTester.countSeries, 0)
        self.assertRaises(AssertionError, bundleTester.countSeries, [5])
        self.assertEqual(bundleTester.countSeries([0]), 1)
        ## sweep count
        self.assertEqual(bundleTester.countSweeps([0, 0]), 50)
        self.assertRaises(AssertionError, bundleTester.countSweeps, 0)
        self.assertRaises(AssertionError, bundleTester.countSweeps, [0])
        ## trace count
        self.assertRaises(AssertionError, bundleTester.countTraces, [0, 2, 5])
        self.assertRaises(AssertionError, bundleTester.countTraces, [0, 5, 5])
        self.assertEqual(bundleTester.countTraces([0, 0, 4]), 1)

    def test_sampleRate(self):
        self.assertAlmostEqual(
            bundleTester.getSeriesSamplingRate([0, 0, 0, 0]), 25_000
        )

    def test_getGroupRecord(self):
        group = bundleTester.getGroupRecord([0, 0, 0, 0])
        self.assertTrue(
            isinstance(group, HEKA.GroupRecord), "failed to get group record"
        )
        self.assertRaises(ValueError, bundleTester.getGroupRecord, [5, 0, 0, 0])

    def test_getSeriesRecord(self):
        series = bundleTester.getSeriesRecord([0, 0, 0, 0])
        self.assertTrue(
            isinstance(series, HEKA.SeriesRecord), "failed to get trace record"
        )

    def test_sweepRecord(self):
        sweep = bundleTester.getSweepRecord([0, 0, 0, 0])
        self.assertTrue(
            isinstance(sweep, HEKA.SweepRecord), "failed to get trace record"
        )

    def test_traceRecord(self):
        trace = bundleTester.getTraceRecord([0, 0, 0, 0])
        self.assertTrue(
            isinstance(trace, HEKA.TraceRecord), "failed to get trace record"
        )

    def test_stimuli(self):
        time, stim, stimInfo = bundleTester.getStim([0, 0, 0, 0])
        self.assertTrue(time.shape == (24998,), "stimuli time stamp")
        self.assertTrue(isinstance(stimInfo, list), "stimuli information as a list")
        self.assertTrue(stimInfo[1]["duration"] == 0.6, "stimuli duration per sweep")
        self.assertAlmostEqual(stimInfo[1]["amplitude"],-990.9489694417275)

    def test_data(self):
        data = bundleTester.getSingleTraceData([0, 0, 0, 0])
        self.assertTrue(isinstance(data, numpy.ndarray))
        self.assertTrue(data.shape[0] == 25_000)
        self.assertTrue(bundleTester.getNumberOfSamplesPerSweep([0, 0, 0, 0]) == 25_000)

    def test_getSweepTimeStamps(self):
        time = bundleTester.getSweepTimeStamps([0, 0, 0, 0])
        self.assertTrue(isinstance(time, numpy.ndarray))
        self.assertTrue(time[0] == 0)

if __name__ == '__main__':
    unittest.main()

