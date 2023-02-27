# -*- coding: utf-8 -*-
"""
Ref.
https://pynwb.readthedocs.io/en/stable/tutorials/general/extensions.html
https://pynwb.readthedocs.io/en/stable/tutorials/domain/icephys.html#sphx-glr-tutorials-domain-icephys-py

ming.hu@bcm.edu
"""
import os
from pynwb import NWBFile, NWBHDF5IO
from pynwb.icephys import CurrentClampStimulusSeries, CurrentClampSeries
from pynwb.file import Subject
from pynwb.file import LabMetaData
import numpy as np
from patchview.HekaIO.HekaHelpers import HekaBundleInfo
import matplotlib.pyplot as plt
from patchview.utilitis.pvEphy import *
from collections import namedtuple
from patchview.utilitis.pvNDX_class import PatchviewSweep, PatchviewSweepGroup
from datetime import datetime
from dateutil.tz import tzlocal
#TODO:
'''
ADD abstract layer for sweep (regardless of clamp mode)
Use this layer to add all three types: current clamp, voltage clamp, and I=0 clamp
ADD voltageclamp/voltageclampstimulusseries
ADD IZeroClampSeries.
numerically, the only difference between all three are stimulus amplitude and value units.
'''
class dat2NWB(pvNWB):
    def __init__(self, datFile, sessionIdx=None, description='intracellular ephys',
    device_name="Heka EPC10_USB Quadro", protocol='patchview', **kw):
        """
        parameters:sessionIdx
        dataFile: dat filename
        sessionIdx: [group, series]  ## one file for one series
        """
        super(dat2NWB, self).__init__(fileName=None, # not a nwb file
         protocol=protocol,
            **kw
        )
        self.sessionIdx = sessionIdx
        self.description = description
        self.device_name = device_name
        self._HekaBundle(datFile)
        self.idenifier = datFile
        if self.sessionIdx is not None:
            self.session_id = self.idenifier
            swps = self.pvdat.countSweeps(sessionIdx)
            sweepIdx = [[sessionIdx[0], sessionIdx[1], c] for c in range(swps)]
            self.makeNWBcontainer(sweepIdx)
            self.setCellData(self.sessionIdx)

    def _HekaBundle(self, datFile):
        try:
            self.pvdat = HekaBundleInfo(datFile)
            self._bundle = self.pvdat.bundle
        except:
            raise HekaFileError("Could not read .dat file format")
        self.pul = self._bundle.pul  ## get HEKA file

    def _gatherNWBFields(self):
        """
        Gather all the fields of NWBFile class
        """
        sb = Subject(None, 'Mus musculus recordings') # add more fields
        setFields = {'session_description':self.description,
        'identifier':self.idenifier,
        'session_start_time':self._bundle.session_start_time,
        'file_create_date':datetime.now(tzlocal()),
        'subject': sb,
        }
        return setFields

    def _initNWB_(self, nwbFields=None):
        '''Initialize NWB file'''
        if nwbFields is None:
            nwbFields = self._gatherNWBFields()
        self.nwbfile = NWBFile(**nwbFields)
        self.device = self.nwbfile.create_device(name=self.device_name)

    def _setStimulus(self, sweepIndex, stimElectrodeID = 0):
        """Stimulating electrode id is usually in the Pulse tree's @Series' label, such as
        "Sti 1S-cc", or "Cell 5 fp". Or in the stimuli tree's @Stimulation's label
        We assume there's only one channel is in use for stimulating. And we use the metadata
        from the first channel.
        Ref: https://pynwb.readthedocs.io/en/stable/pynwb.icephys.html
        IZeroClampSeries: spontaneous current clamp series should use this class
        """
        self.stimInfo = []
        for j, idx in enumerate(sweepIndex):
            ## TODO: need to work with Ampl file in the .dat file to get the "gain"
            # datSerie = self.pul[idx[0]][idx[1]]
            # seriesName = datSerie.Label
            if j<10:
                sName = '00'+str(j)
            elif j<100:
                sName = '0'+str(j)
            else:
                sName = str(j)
            
            stimTime, stimData, stimInfo = self._bundle.stim(idx)
            self.stimInfo.append(stimInfo)
            ccs = CurrentClampStimulusSeries(
                    name= "sw" + sName,
                    data=np.array(stimData),
                    offset = stimInfo[0]['Vholding'],
                    stimulus_description= stimInfo[0]['EntryName'],
                    starting_time=stimInfo[0]['duration'],
                    rate=1.0 / stimInfo[0]["sampleInteval"],
                    electrode=self.ccrElectrodes[stimElectrodeID],
                    gain=0.02,  # starting_time = stimInfo[stimEpoch]['start'] * stimInfo[0][1]['sampleInteval']
                    sweep_number=j,
                    conversion=1e-12,
                    unit="amperes",
                )
            ccs.bias_current = stimInfo[0]['Vholding']
            self.nwbfile.add_stimulus(ccs, use_sweep_table=True)

    def _setRecordingElectrodes(self, sweepIndex):
        ## get all electrodes
        ## selected session's index
        idx_ = sweepIndex[-1] # assume all sweeps have same setup.
        assert len(idx_)>=2, f'index should have length at least 2. Got {len(idx_)}'
        idx = idx_.copy()
        idx.append(0)
        self.ccrElectrodes = []  ## current clamp recording electectrodes
        for trace in self.pul[idx[0]][idx[1]][idx[2]].children:
            #            print(trace)
            self.ccrElectrodes.append(
                self.nwbfile.create_icephys_electrode(
                    name=trace.Label, description="", device=self.device
                )
            )
        return idx

    def makeCurrentClampSweep(self, p:namedtuple) -> CurrentClampSeries:
        ''' basic unit should be sweep
        sweep_number: sweep number added to the sweep table. This may be different than original sweep number
        in the dat file
        sweepName: name of the sweep, inherited from the parent series. contain original sweep number
        data: 2D array of data, each row is a trace
        eletrodeObj: the electrode object
        srate: sampling rate
        gain: gain of the amplifier
        stimulus description: is it a firing pattern or a connection test etc.
        return a currentclamp sweep object
        '''
        return CurrentClampSeries(
                    name = p.sweepName,
                    data = p.data,
                    rate = p.srate,
                    electrode = p.electrodeObj,
                    conversion = 1.0,
                    unit ="volts",
                    gain = p.gain,
                    sweep_number= p.sweep_number,
                    bias_current = p.bias_current,
                    bridge_balance = p.bridge_balance,
                    capacitance_compensation = p.capacitance_compensation,
                    stimulus_description = p.label,
                    starting_time = p.starting_time,
                )

    def makeNWBcontainer(self, sweepindex: list):
        self._initNWB_() # create NWB file object
        self._setRecordingElectrodes(sweepindex) 
        self._setStimulus(sweepindex) # stimulus data last

    def AddSweeps(self, sweepIndex: list[list[int, int, int]], sweepGroupName=None, traceIndex: list=[]):
        ''' Add sweep to NWB file
        TODO: dispatch different sweep type to different function according to sweep group specification.
        By default, respect the original group-series-sweep structure in dat file.
        Reorganizing the sweeps should be done in separate function.
        We are building a completely new sweep table here. So sweep table index may be different than the
        original one. Save a copy of the source sweep table index.

        Parameters
        -----------
        `sweepIndex`:  list of sweep index (vec3 int) from original dat file
        `sweepGroupName`: name for a new or existing sweep group. TODO: check if a sweepgroup with the same name exist
        `traceIndex`: list of trace index to add. Default, [], all traces available
        '''
        if not hasattr(self, 'nwbfile') or (self.nwbfile is None): # if nwbfile is not created yet
            self.makeNWBcontainer(sweepIndex)
        # parameters for a typical currentclam series
        currentClamp_init_namedTuple = namedtuple('currentClamp_init_namedTuple',
             'sweep_number, sweepName, data, electrodeObj, srate, gain, label,\
                 bias_current, bridge_balance, capacitance_compensation,starting_time')
        # if sweepGroupName is None: 
        #     hashVal = hash(tuple([s2 for s1 in sweepIndex for s2 in s1]))
        #     sweepGroupName = '#'+str(hashVal)
        swIdx = sweepIndex[0]
        if sweepGroupName is None:
            label = self.pul[swIdx[0]][swIdx[1]].Label # parent of current sweep
            label = label.replace(' ', '_')
        else:
            assert isinstance(sweepGroupName, str), f'sweepGroupName should be a string. Got {type(sweepGroupName)}'
            label = sweepGroupName
        pv_sweep_group = PatchviewSweepGroup(name=label) # sweep group object for patchview. annotation purpose
        sweepIndex = sorted(sweepIndex, key=lambda z: (z[1],z[2])) # sort by sweep index of dat file  
        for idx, swIdx in enumerate(sweepIndex):
            datSerie = self.pul[swIdx[0]][swIdx[1]] # parent of current sweep
            serieLabel = 'Series'+str(swIdx[1]+1)+' '+ datSerie.Label
            stimulus_type  = self.querySessionProtocolByLabel(datSerie.Label) # parent class method
            stimulus_name  = serieLabel #+ '_'+stimulus_type
            sweep = datSerie.children[swIdx[2]] # current sweep object
            traces = sweep.children # number of traces in this sweep
            if len(traceIndex) == 0: # default, all traces
                nTrace = len(traces)
                traceIndex = range(nTrace) 
            else: # selected 
                nTrace = len(traceIndex)
            baseIdx = swIdx.copy()
            baseIdx.append(0)
            nTimePoint = len(self._bundle.data[baseIdx])
            data = np.zeros((nTimePoint, nTrace))
            trace_indexes = []
            trace_labels = []
            newTrace_index = 0
            for t, trace in enumerate(sweep.children):
                if not (t in traceIndex):
                    continue
                baseIdx[3] = t
                trace_indexes.append(baseIdx[3])
                trace_labels.append(trace.Label)
                data[:, newTrace_index] = self._bundle.data[baseIdx]
                newTrace_index +=1
            # design is wired here. PatchviewSweep sweep_index is not the same as Sweep table index
            # The latter is global indexed, the former is relative to each sweep group.
            stim_sampleInteval = self.stimInfo[idx][1]["sampleInteval"]
            stim_onset = self.stimInfo[idx][1]['start']*stim_sampleInteval
            stim_offset = self.stimInfo[idx][1]['end']*stim_sampleInteval
            stim_amplitude = self.stimInfo[idx][1]['amplitude']
            stim_holding = self.stimInfo[idx][1]['Vholding']
            pv_sweep = PatchviewSweep(name=f'{idx:03}', stimulus_type=stimulus_type,  # stimulus_type is protocol-aware. stimulus_name are arbitury
            stimulus_name=stimulus_name, local_sweep_index = idx, sweep_table_index= idx, #swIdx[2],
            stimulus_onset = stim_onset, stimulus_offset = stim_offset, stimulus_amplitude = stim_amplitude, 
            stimulus_holding_amplitude = stim_holding,
            trace_indexes = np.array(trace_indexes), trace_labels=trace_labels)
            pv_sweep_group.add_patchview_sweep(pv_sweep)
            eletrodeObj = self.nwbfile.get_ic_electrode(trace_labels[0]) # get the first trace's electrode
            sweepName = serieLabel + '_sw_' + str(swIdx[2])+'_'+ str(t+1)+' traces'
            sw_label = '_'.join([str(d) for d in swIdx.copy()])
            currentClamp_init = currentClamp_init_namedTuple(idx, sweepName,
             data, eletrodeObj, 1.0 / trace.XInterval, gain=0.02, label=sw_label,
                bias_current = 0.0, bridge_balance = 0.0, capacitance_compensation = 0.0,
                starting_time=0.0)
            self.nwbfile.add_acquisition(self.makeCurrentClampSweep(currentClamp_init),
            use_sweep_table=True)
        gpName = serieLabel.replace(' ', '_')
        self.addPachviewSweepGroup({gpName: pv_sweep_group}) # list of PatchviewSweepGroup

    def getSweepIdx(self, sessionIdx:list[int, int]) -> list[list[int, int, int]]:
        datSerie = self.pul[sessionIdx[0]][sessionIdx[1]]
        nSweep = len(datSerie.children)
        return [[sessionIdx[0], sessionIdx[1], s] for s in range(nSweep)]

    def setCellData(self, sessionIdx, use_protocal=True) -> None:
        '''quick way to add all sweeps from a series to NWB file'''
        sweepIdx = self.getSweepIdx(sessionIdx) # get all sweeps from a series
        self.AddSweeps(sweepIdx)
        
class HekaFileError(Exception):
    """Generic Python-exception-derived object raised by any function."""
    pass

if __name__ == "__main__":
    _FILEPATH  = "D:\\DataDemo"
    # fn = os.path.join(_FILEPATH, "test_singleFP.dat")
    fn = os.path.join(_FILEPATH, "190628", "190628s1c6fp.dat")
    datnwb = dat2NWB(fn, [0,0])
    print('dat file is loaded')
    dfs = sorted(datnwb.nwbfile.fields.keys())
    print(dfs)
    # test save nwb
    datnwb.saveNWB(os.path.join(_FILEPATH, 'test.nwb'))
    print('nwb file saved')
    nSweep1 = datnwb.totalNumOfSweeps
    # test load saved nwb with pvNWB Class
    datnwb2 = pvNWB(os.path.join(_FILEPATH,"test.nwb"))
    print('data loaded back by pvNWB')
    dfs2 = sorted(datnwb.nwbfile.fields.keys())
    # extra keys when read back by pvNWB. TODO: need to check why
    # assert  dfs == dfs2, 'keys are not the same'
    # test get data from saved nwb
    nSweep2 = datnwb2.totalNumOfSweeps
    assert nSweep1 == nSweep2, 'number of sweeps are not the same'
    print('sweep number:', nSweep2)
    # read a sweep
    stim, data, meta = datnwb2.getSweepData(0)
    # print('data shape:', data.shape)
    print('meta:', meta)
    # test plot a sweep
    sweepIdx = [0, nSweep2//2, nSweep2-1]
    fig = datnwb2.plotSweep(sweepIdx)
    plt.show()

