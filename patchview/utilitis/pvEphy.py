# -*- coding: utf-8 -*-
"""
A class for Patchview wrap on NWBFile.
providing a unified interface for patchview to handle epysiology data
ming.hu@bcm.edu
"""
import pynwb
from pynwb import NWBFile, NWBHDF5IO
from pynwb.file import Subject
from pynwb.file import LabMetaData
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
import functools
from patchview.utilitis.pvNDX_class import PatchviewSweep, PatchviewSweepGroup
from typing import Union
import os
import yaml

class pvNWB():
    def __init__(self, fileName=None, protocol='patchview',  **kw):
        """ base class for Patchview to handle NWB file
        TODO: rewrite this part to use class that follows Neurodata extension (NDX) protocol
        https://pynwb.readthedocs.io/en/stable/pynwb.file.html#pynwb.file.NWBFile.add_analysis
        the currentclamp series object' stimulus_description is used to
        define sweep groups in patchview. it is a string that can be parsed. think of json format.
        """
        self.nwbfile = None
        self.metaInfo = None
        self.protocol = protocol
        if fileName is not None: # most generic way to load NWB file
            readNWBStatus = self.loadNWB(fileName)
            self.source_filename = fileName
            assert readNWBStatus == 0, 'failed to load NWB file'
            self.loadProtocols(protocol) 
            self.getNWB_metaInfo() # get meta information from NWB file

    def getSweepGroups(self):
        if not hasattr(self, 'sweepGroups'):
            self.sweepGroups = self.setPachviewSweepGroup()
        assert len(self.sweepGroups) > 0, 'no sweepgroup found'
        return self.sweepGroups

    def getSweepGroup(self, idx=0):
        if not hasattr(self, 'sweepGroups'):
            self.sweepGroups = self.setPachviewSweepGroup()
        assert len(self.sweepGroups) > idx, f'index error'
        return self.sweepGroups[idx]        

    def setPachviewSweepGroup(self) -> list[PatchviewSweep]:
        '''Get list of `PatchviewSweep` from NWB file. Create default one if non-exist.

        return
        ------
        return a list of PatchviewSweepGroup
        '''
        promods = self.nwbfile.processing
        sweepGroups = []
        for modName in promods: # assume only one sweep group exist
            for c in promods[modName].children:
                if c.neurodata_type == 'PatchviewSweepGroup':
                    sweepGroups.append(c)
        if len(sweepGroups) == 0:
            print('No patchview sweep group found in NWB file. Creating default sweep group.')
            sweepGroups = self.defineSweepGroups() # overwrite this function to define custom sweep groups
            self.addPachviewSweepGroup(sweepGroups)
        return sweepGroups

    def addPachviewSweepGroup(self, sweepGroup: Union[PatchviewSweepGroup, list[PatchviewSweepGroup]]):
        ''' Add sweep group to NWB file's processing module `patchview_sweep_group`.
        Create it if not exist.

        Parameter
        ---------
        PatchviewSweepGroup: pynwb NDX extension for patchview to handle sweep group
        '''
        assert hasattr(self, 'nwbfile'), "NWB file is not created yet. Please run makeNWBcontainer() first."

        if 'patchView_sweep_group' in self.nwbfile.processing:
            promod_pvgrp = self.nwbfile.get_processing_module('patchView_sweep_group')
        else:
            promod_pvgrp = self.nwbfile.create_processing_module('patchView_sweep_group', 'sweep grouping scheme') # add processing module
        if not isinstance(sweepGroup, list):
            assert isinstance(sweepGroup, PatchviewSweepGroup), f'type of sweepgroup should be PatchviewSweepGroup,Got {type(sweepGroup)}'
            sweepGroup = [sweepGroup]
        else:
            for swp in sweepGroup:
                assert isinstance(swp, PatchviewSweepGroup), f'type of sweepgroup should be PatchviewSweepGroup. Got {type(swp)}'                
        promod_pvgrp.add(sweepGroup)
        
    @property
    def subject(self)-> pynwb.file.Subject:
        ''' return subject info
        fields: age, dob, genotype, sex, species, subject_id
        '''
        return self.nwbfile.fields['subject']

    def saveNWB(self, filename):
        ''' Save as Nb format
        '''
        if self.nwbfile is None:
            print('no data to save')
            return
        with NWBHDF5IO(filename, "w") as io:
            io.write(self.nwbfile)

    def loadNWB(self, filename):
        ''' load from NWB format
        '''
        try:
            self.io = NWBHDF5IO(filename, "r", load_namespaces=True)
            self.nwbfile = self.io.read()
            return 0
        except Exception as e:
            print(f'failed to load file:{filename}')
            print(e)
            return -1

    def getNWB_metaInfo(self):
        ''' gather meta information from NWB file into a lightweight named tuple
        this is where Patchviw will look for meta information about experiment.
        subject info is also copied for convenience. 
        '''
        assert self.nwbfile is not None, 'no NWB file loaded yet'
        fields = ['session_description', 'identifier', 'session_start_time', 
        'experimenter', 'institution'] # selected properties from NWB file
        nwbFields = self.nwbfile.fields.keys()
        fieldsVal = [self.nwbfile.fields[p] if p in nwbFields else None for p in fields]
        if hasattr(self, 'subject'): # subject info if available
            subjetFields = self.subject.fields.keys()
            subjetFields_val = [self.subject.fields[p] for p in subjetFields]
            fields +=self.subject.fields.keys() # add subject fields
            fieldsVal += subjetFields_val # add subject fields values
        # make a namedtuple type
        baseMetaInfo_nt = namedtuple('NWB_metaInfo', fields,
         defaults=(None,)*len(fields), module='Patchview')
        self.metaInfo = baseMetaInfo_nt(*fieldsVal) # ._asdict() to get a dictionary

    @functools.cached_property
    def totalNumOfSweeps(self) -> int:
        assert self.nwbfile is not None, 'no NWB file loaded yet'
        assert hasattr(self.nwbfile, 'sweep_table'), 'no sweep table in NWB file'
        return len(self.nwbfile.stimulus) # what about IZeroSeries

    def querySessionProtocolByLabel(self, label: Union[str, None]) -> str:
            ''' Sort session by its label. Heka saves protocol label in session label field.
            A good practice is to set a informative protocol name in your Heka protocol file.
            Heka file unique feature. Not generlizable to other file
            Parameters
            ----------- 
            lable: read from dat file's Label or None

            Return
            -------
            Stimulus type. Str type for now. Should standardarized as object.
            '''
            if label is None or (not hasattr(self, 'protocol_dict') or (self.protocol_dict is None)):
                return "Custom"
            label = label.replace("-", " ")
            for x in label.split():
                if x in self.protocol_dict["Firing pattern"]:
                    stimType = "Firing pattern"
                    break
                elif x in self.protocol_dict["Connection"]:
                    stimType ="Connection"
                    break
                elif x in self.protocol_dict["Spontaneous"]:  # IZero time-series
                    stimType ="Spontaneous"
                    break
            else:
                stimType = "Custom"
            return stimType

    def loadProtocols(self, protocol):
        '''Load protocol file. 
        If not default, expect a yaml file with a mapping between user's protocol label to standard clamp series
        '''
        if protocol == None:
            self.protocol_dict = None # no mapping needed

        elif protocol== 'patchview': # patchview default. Temporary hack; make here more generic!
            try:
                from appdirs import user_config_dir ## if this is not installed, install it with pip install appdirs
                # assume that the config file is in the Patchview config folder
                useConfigDir = user_config_dir("Patchview")
                confilePath  = os.path.join(useConfigDir,"patchview.yaml")
                with open(confilePath) as f:
                    pars = yaml.load(f, Loader=yaml.FullLoader)
                self.protocol_dict= pars["Protocols"]
                print(f'Protocols loaded from {confilePath}.')
            except Exception as e: 
                print(e)
                self.protocol_dict = None
                print("No protocol file found in Patchview config folder.")
        else:
            try: # assume self.protocol is a file path to the user defined protocol file
                with open(self.protocol) as f:
                    pars = yaml.load(f, Loader=yaml.FullLoader)
                self.protocol_dict= pars["Protocols"]
                print(f'Protocols loaded from {self.protocol}.')
            except Exception as e: 
                print(e)
                self.protocol_dict = None
                print(f"Failing to load protocols from {self.protocols}.")
                print('Please check the file path or the file format to be Yaml.')
                print('Check Patchviw configuration file for example format')  

        if self.protocol_dict is None:
            self.use_protocol = False
        else:
            self.use_protocol = True
        # print(f'loaded protocol: {protocol}\n')
        # print(self.protocol_dict)

    def defineSweepGroups(self) -> list[PatchviewSweep]:
        ''' Create a default intracelluar ephy. sweep group.
        Overwrite this func in a inherit class. See reference implementation from dandi2pvNWB.py
        '''
        pv_sweep_group_v = PatchviewSweepGroup(name = 'voltage_clamp') # sweep group object for patchview. annotation purpose
        pv_sweep_group_c = PatchviewSweepGroup(name = 'current_clamp') 
        pv_sweep_group_u = PatchviewSweepGroup(name = 'custom_clamp') 
        v, c, u = 0, 0, 0
        sweepIndies = range(self.totalNumOfSweeps)
        stimulus = self.nwbfile.stimulus
        stim_id = list(stimulus.keys())
        for idx, swIdx in enumerate(sweepIndies):
            trace_indexes = [0] # single channel case
            trace_labels =  [''] # number of elements need to be match the number of traces in each sweep
            cstim = stimulus[stim_id[idx]]
            neurodata_type = cstim.neurodata_type
            pv_sweep = PatchviewSweep(name=f'sweep {idx:03}', stimulus_type = neurodata_type,
             stimulus_name='stim', sweep_index = idx, source_sweep_table_index = swIdx,
                            trace_indexes = np.array(trace_indexes), # this is a list of list, inner is a list of 4 integers
                            trace_labels=trace_labels)
            if neurodata_type =='VoltageClampStimulusSeries':
                pv_sweep_group_v.add_patchview_sweep(pv_sweep)
                v+=1
            elif neurodata_type =='CurrentClampStimulusSeries':
                pv_sweep_group_c.add_patchview_sweep(pv_sweep)
                c+=1
            else:
                pv_sweep_group_u.add_patchview_sweep(pv_sweep)
                u+=1
        sweepGroups = []
        if c >0:
            sweepGroups.append(pv_sweep_group_c) # pv_sweep_group_v.patchview_sweeps
        if v >0:
            sweepGroups.append(pv_sweep_group_v)
        if u >0:
            sweepGroups.append(pv_sweep_group_u)
        return sweepGroups

    # def defineSweepGroups_namedtuple(self):
    #     ''' Depreted.sort sweep into groups based on stimulus type.
    #     patchview would use this to make sweep groups in the viewer.
    #     A NWB file can contain multiple stimulus types.
    #     User should overwrite this function to make sweep groups.
    #     fields: a list of group names. Default is sweeps (all sweeps, No group)
    #     fields_data: a list of namedtuple: `sweepgroup`: `sweep indices`. Default is all sweeps
    #     return a (nested) namedtupple as sweep group name and value as a list of sweep indices
    #     check example in dandi2pvNWB.py
    #     '''
    #     fields = ['Session'] # default group
    #     self.sweepGroup_def = namedtuple('pvSwgDef', fields, module='Patchview') ## default group
    #     # for each field, create a namedtuple with a list of sweep indices
    #     # By default, one field with all sweeps
    #     pvclampData = {"Sweep_table":range(self.totalNumOfSweeps)} # voltage clamp. to be extended...
    #     sweepGrp_namedtuple = namedtuple('Data', pvclampData.keys())
    #     sweepGrpData = sweepGrp_namedtuple(*pvclampData.values())

    #     # gather all fields' namedtuples into a list and make the final namedtupple with all fields
    #     fields_data = [sweepGrpData] # default group
    #     self.sweepGroups = self.sweepGroup_def(*fields_data)

    def _getSweepTable_structureInfo(self, seriesData):
        ''' information about sweep table series output. stim first or responses first etc.
        Order in the list are arbituary. need to check to dispatch
        '''
        stimIndex = []
        for idx, s in enumerate(seriesData):
            if 'Stimulus' in s.neurodata_type:
                stimIndex.append(idx)
        return stimIndex

    def getSweepData(self, i):
        ''' Get sweep with index i
        return stimulus and response as ndarray
        metaInfo is a dictionary containing meta information of the sweep

        Parameter
        ---------
        i: sweep table index

        Return
        ------
        stim: ndarray of stimulus
        responses: list of numpy ndarray. one for each channel
        metainfo: dictionary of stimulus and response metainformation
        '''
        # self.nSweep = self.getNumberOfSweeps
        if i >=0 and i<self.totalNumOfSweeps:
            seriesData = self.nwbfile.sweep_table.get_series(i)
            stimIndex = self._getSweepTable_structureInfo(seriesData)
            assert len(stimIndex) <= 1, 'more than two stim records found!'
            recordings = [seriesData[i] for i in range(len(seriesData)) if not (i in stimIndex)] # recording is a list of number of recording channels
            responses = [r.data*r.conversion for r in recordings]# extract data and convert to ndarray
            metaInfo = {}
            if stimIndex:
                sti = seriesData[stimIndex[0]] # last entry is stimulus
                stimulus = sti.data*sti.conversion
                metaInfo['stim_rate'] = sti.rate
                metaInfo['stim_time_unit'] = sti.time_unit
                metaInfo['stim_unit'] = sti.unit
            else:
                metaInfo['stim_rate'] = np.nan
                metaInfo['stim_time_unit'] = np.nan
                metaInfo['stim_unit'] = np.nan                
            metaInfo['recording_rate'] = recordings[0].rate
            metaInfo['recording_unit'] = recordings[0].unit
            metaInfo['recording_time_unit'] = recordings[0].time_unit
            for p in ['resistance_comp_bandwidth', 'resistance_comp_correction',
            'capacitance_fast', 'capacitance_slow', 'resistance_comp_prediction',
            'whole_cell_capacitance_comp', 'whole_cell_series_resistance_comp']:
                if hasattr(recordings, p):
                    metaInfo[p] = eval('recording[0].'+p)
                else:
                    metaInfo[p] = ''
            return stimulus, responses, metaInfo
        else:
            print('sweep index is not valid. Number of sweep should be less than ', self.totalNumOfSweeps)
            
    def _getSweepData(self, sweepObj):
        ''' convert sweep to array
        '''
        return sweepObj.data*sweepObj.conversion 
    
    def plotSweep(self, i, axes=None):
        ''' plot sweep of response and stimulus
            if provide an plotting axes, make it contain at least two subplots panel.
            By default, stimulus will be in the second panel.
        '''
        if axes is None:
            fig, axes = plt.subplots(2,1, sharex=True)
        if type(i) is int:
            i = [i]
        assert type(i) is list, "index to be integer or list of integers"
        for s in i: ## plot all series
            stimulus, responses,metaInfo = self.getSweepData(s)
            time_stim = np.arange(len(stimulus))/metaInfo['stim_rate']
            time_data = np.arange(len(responses[0]))/metaInfo['recording_rate']
            for response in responses:
                axes[0].plot(time_data, response)
            axes[1].plot(time_stim, stimulus)
            axes[0].set_xlabel(metaInfo['stim_time_unit'])
            axes[0].set_ylabel(metaInfo['recording_unit'])
            axes[1].set_ylabel(metaInfo['stim_unit'])
            axes[1].set_xlabel(metaInfo['recording_time_unit'])
        plt.show()
        return fig

if __name__ == "__main__":
    import os
    _cd_ = os.path.dirname(os.path.abspath(__file__))
    os.chdir(_cd_)
    testFile = "D:\\DataDemo\\nwb\\sub-626194732_ses-638205339_icephys.nwb"
    # testFile = "sub-626194732_ses-638205339_icephys.nwb" 
    # _FILEPATH = os.path.join('..', '..','tests','data')
    # nwbfn = os.path.join(_FILEPATH, testFile)
    nwbfn = testFile
    x = pvNWB(nwbfn)
    print(f'num sweeps: {x.totalNumOfSweeps}')
    nwb = x.nwbfile
    swGrp = x.defineSweepGroups() #func = range(x.totalNumOfSweeps))
    stim, data = nwb.sweep_table.get_series(10)
    stim, data, meta = x.getSweepData(10)
    x.plotSweep([0, 4, 8, 12, 20, 25])
