# -*- coding: utf-8 -*-
"""
An interface for Patchview to access(read/write) on NWBFile.
ming.hu@bcm.edu
"""
import os, sys
import pynwb
from pynwb import NWBFile, NWBHDF5IO
from pynwb.file import Subject
from pynwb.file import LabMetaData
from pynwb.icephys import CurrentClampStimulusSeries, VoltageClampStimulusSeries
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
import functools
from patchview.utilitis.pvNDX_class import PatchviewSweep, PatchviewSweepGroup
from typing import Union
import yaml
from datetime import datetime
from dateutil.tz import tzlocal

class pvNWB():
    def __init__(self, fileName=None, protocol=None,  **kw): #protocol='patchview'
        """ base class for Patchview to handle NWB file
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
    
    @classmethod
    def make_new_NWBcontainer(cls, src: Union[NWBFile, None]) -> NWBFile:
        selected_fields = [
        'experimenter',
        'experiment_description', 'session_id', 'institution',
        'keywords', 'notes', 'pharmacology', 'protocol', 'related_publications', 'slices',
        'source_script', 'source_script_file_name', 'data_collection', 'surgery', 'virus',
        'stimulus_notes', 'lab', 'epochs', 'trials', 'units']
        required_fields =['session_description', 'identifier',
         'session_start_time','file_create_date','subject']
        if src is not None:
            reqFields = {}
            for f in required_fields:
                if hasattr(src, f):
                    reqFields[f] = getattr(src, f)
            deviceNames = list(src.devices.keys()) 
            devices = [src.devices[dn] for dn in deviceNames]
            for d in devices:
                d.reset_parent()
                d.set_modified()
                d._in_construct_mode = True
                d._AbstractContainer__container_source = None
            nwbfile = NWBFile(**reqFields, devices=devices)
            cls.reset_container(nwbfile.subject, nwbfile)
            return nwbfile           
        else:
            rqFields = {'session_description': 'no description',
            'identifier': 'no identifier',
            'session_start_time': datetime.now(tzlocal()),
            'file_create_date': datetime.now(tzlocal()),
            'subject': Subject('no subject', 'no species')}
            # copy a selected attributes from src and add to rqFields
            reqFields = add_attributes(src, selected_fields, rqFields)
            return NWBFile(**reqFields)

    @classmethod
    def getNumOfSweeps(cls, nwb: NWBFile) -> int:
        assert isinstance(nwb, NWBFile), 'input should be NWBFile'
        return len(nwb.stimulus)

    def set_subject(self, **kw):
        '''Set subject information in NWB file
        kw = {
        'age',
        'description',
        'genotype',
        'sex',
        'species',
        'subject_id',
        'weight',
        'date_of_birth',
        'strain'
        }
        '''
        return Subject(**kw)

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

    @functools.cached_property
    def sweepGroupNames(self):
        return list(self.sweepGroups.keys())

    def getSweepGroups(self):
        if not hasattr(self, 'sweepGroups'):
            self.sweepGroups = self.setPachviewSweepGroup()
        assert len(self.sweepGroups) > 0, 'no sweep groups found'
        return self.sweepGroups

    def getSweepGroup(self, sweep_group_name):
        if not hasattr(self, 'sweepGroups'):
            self.sweepGroups = self.setPachviewSweepGroup()
        assert sweep_group_name in self.sweepGroups.keys(), f'sweep group {sweep_group_name} not found'
        return self.sweepGroups[sweep_group_name]        

    def setPachviewSweepGroup(self) -> dict[str, PatchviewSweepGroup]:
        '''Create a list of `PatchviewSweepGroup` from NWB file. Create default one if non-exist.

        return
        ------
        return a list of PatchviewSweepGroup
        '''
        promods = self.nwbfile.processing
        sweepGroups = {}
        for modName in promods: # assume only one sweep group exist
            for c in promods[modName].children:
                if c.neurodata_type == 'PatchviewSweepGroup':
                    sweepGroups.update({c.name: c})
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
        assert isinstance(sweepGroup, dict), f'type of sweepgroup should be dictionary, Got {type(sweepGroup)}'
        for swp in sweepGroup.keys():
            assert isinstance(sweepGroup[swp], PatchviewSweepGroup), f'type of sweepgroup should be PatchviewSweepGroup. Got {type(swp)}'                
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
        assert self.nwbfile is not None, 'no NWB file loaded'
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

    def defineSweepGroups(self) -> dict:
        ''' Create a default intracelluar ephy. sweep group.
        Overwrite this func in a inherit class. See reference implementation from dandi2pvNWB.py
        '''
        # SweepGroup = {'voltage_clamp':[],
        # 'current_clamp':[],
        # 'custom_clamp':[]}
        # for k in SweepGroup.keys():
        #     SweepGroup[k] = self.nwbfile.sweep_table.where(f'stimulus_type="{k}"')
        pv_sweep_group_v = PatchviewSweepGroup(name = 'voltage_clamp') # sweep group object for patchview. annotation purpose
        pv_sweep_group_c = PatchviewSweepGroup(name = 'current_clamp') 
        pv_sweep_group_u = PatchviewSweepGroup(name = 'custom_clamp') 
        v, c, u = 0, 0, 0
        sweepIndies = list(set(list(self.nwbfile.sweep_table.sweep_number))) #range(self.totalNumOfSweeps)
        stimulus = self.nwbfile.stimulus
        stim_id = list(stimulus.keys())
        for idx, swIdx in enumerate(sweepIndies):
            trace_indexes = [0] # single channel case
            trace_labels =  [''] # number of elements need to be match the number of traces in each sweep
            cstim = stimulus[stim_id[idx]]
            neurodata_type = cstim.neurodata_type

            if neurodata_type =='VoltageClampStimulusSeries':
                v+=1
                sweep_idx = v
                sweep_group = pv_sweep_group_v
            elif neurodata_type =='CurrentClampStimulusSeries':
                c+=1
                sweep_idx = c
                sweep_group = pv_sweep_group_c
            else:
                u+=1
                sweep_idx = u
                sweep_group = pv_sweep_group_u
            pv_sweep = PatchviewSweep(name=f'sweep {sweep_idx:03}', stimulus_type = neurodata_type,
             stimulus_name='Series', local_sweep_index = sweep_idx, sweep_table_index = int(swIdx),
                            trace_indexes = np.array(trace_indexes), # for Heka, only one trace per sweep
                            trace_labels=trace_labels)
            sweep_group.add_patchview_sweep(pv_sweep)
        sweepGroups = {}
        if c >0:
            sweepGroups.update({pv_sweep_group_c.name: pv_sweep_group_c})
        if v >0:
            sweepGroups.update({pv_sweep_group_v.name:pv_sweep_group_v})
        if u >0:
            sweepGroups.update({pv_sweep_group_u.name:pv_sweep_group_u})
        return sweepGroups

    def _getSweepTable_structureInfo(self, seriesData):
        ''' information about sweep table series output. stim first or responses first etc.
        Order in the list are arbituary. need to check to dispatch
        '''
        stimIndex = []
        for idx, s in enumerate(seriesData):
            if 'Stimulus' in s.neurodata_type:
                stimIndex.append(idx)
        return stimIndex

    def getMetaInfoFromSweepGroups(self, sweepGroupKey, sweepIdx, skipChecking=False) -> dict[str, float]:
        ''' Get extra meta information not saved in sweep table from sweep groups

        Parameters
        ----------
        sweepGroupIndex : str. Name of the sweep group
        sweepIdx : int. Sweep index in the sweep group
        
        Returns
        -------
        metaInfo : dict
        '''
        if not skipChecking:
            self.boundCheck(sweepGroupKey, sweepIdx)
        swGrp = self.sweepGroups[sweepGroupKey] # default the last group
        metaInfo = {}
        metaInfo['stimulus_onset'] = swGrp.get_stimulus_onset(sweepIdx)
        metaInfo['stimulus_offset'] = swGrp.get_stimulus_offset(sweepIdx)
        metaInfo['stimulus_amplitude'] = swGrp.get_stimulus_amplitude(sweepIdx)
        metaInfo['stimulus_holding_amplitude'] = swGrp.get_stimulus_holding_amplitude(sweepIdx)
        return metaInfo

    def boundCheck(self, sweepGroupKey, sweepIndex):
        '''Check if the key and the sweep index is valid
        '''
        assert sweepGroupKey in self.sweepGroups.keys(), f'No sweep group named {sweepGroupKey}'
        assert sweepIndex < len(self.sweepGroups[sweepGroupKey].patchview_sweeps), f'No sweep\
             with index {sweepIndex} in sweep group {sweepGroupKey}'
 
    def getSweepData(self, sweep_table_Index):
        ''' Get sweep with index i
        return stimulus and response as ndarray
        metaInfo is a dictionary containing meta information of the sweep

        Parameter
        ---------

        sweep_table_Index: sweep_index in the sweep table 

        Return
        ------
        stim: ndarray of stimulus
        responses: list of numpy ndarray. one for each channel
        metainfo: dictionary of stimulus and response metainformation
        '''
        # self.nSweep = self.getNumberOfSweeps
        # should we create another seperate dictionay to map between sweep group and sweep table index?
        # self.boundCheck(sweepGroupKey, sweepIndex)
        # i = self.sweepGroups[sweepGroupKey].get_sweep_table_index(sweep_index)
        seriesData = self.nwbfile.sweep_table.get_series(sweep_table_Index)
        stimIndex = self._getSweepTable_structureInfo(seriesData)
        assert len(stimIndex) <= 1, 'more than two stim records found!'
        recordings = [seriesData[i] for i in range(len(seriesData)) if not (i in stimIndex)] # recording is a list of number of recording channels
        responses = [r.data*r.conversion for r in recordings]# extract data and convert to ndarray
        metaInfo = {}
        if stimIndex:
            sti = seriesData[stimIndex[0]] # last entry is stimulus
            stimulus = sti.data*sti.conversion
            for p in ['rate','time_unit','unit','game','stimulus_description','comments']:
                if hasattr(sti, p):
                    metaInfo['stim_'+p] = eval('sti.'+p)
                else:
                    metaInfo['stim_'+p] = ''                    
        else:
            metaInfo['stim_rate'] = np.nan
            metaInfo['stim_time_unit'] = np.nan
            metaInfo['stim_unit'] = np.nan
            stimulus = np.nan                
        for p in ['rate','unit','time_unit','bridge_balance', 'capacitance_compensation', 'gain', 'starting_time','conversion', \
            'resistance_comp_bandwidth', 'resistance_comp_correction',\
        'capacitance_fast', 'capacitance_slow', 'resistance_comp_prediction',
        'whole_cell_capacitance_comp', 'whole_cell_series_resistance_comp', 'comments']:
            if hasattr(recordings[0], p):
                metaInfo['recording_'+p] = eval('recordings[0].'+p)
            else:
                metaInfo['recording_'+p] = ''
        return stimulus, responses, metaInfo
    
    @classmethod
    def add_nwb_sweep(cls, nwbfile:NWBFile, NWB_sweepObjList:list, use_sweep_table:bool=True,
     keep_original_sweep_number:bool=False):
        ''' Add sweep to NWB file
        Since we are at pvEphys level, we need some abstractions to add sweeps to NWB file.
        Shall we use ABC or Protocol to define the interface?
        Current function is too convoluted. Need to be refactored into smaller functions.
        TODO: dispatch different sweep type to different function according to sweep group specification.
        By default, respect the original group-series-sweep structure in dat file.
        Reorganizing the sweeps should be done in separate function.
        We are building a completely new sweep table here. So sweep table index may be different than the
        original one. Save a copy of the source sweep table index.

        Parameters
        -----------
        `nwbfile`:  NWB file object

        `NWB_sweepObjList`: list of valid NWB sweep objects: CurrentClampStimulusSeries | VoltageClampStimulusSeries
        `use_sweep_table`: bool
        `keep_original_sweep_number`: bool. If true, use the original sweep number. Otherwise, increment sweep number by 1
        '''       
        for swb in NWB_sweepObjList:
            cls.reset_container(swb, nwbfile)           
            cls.reset_container(swb.electrode.device, nwbfile)
            if swb.electrode.name not in nwbfile.ic_electrodes.keys():
                if swb.electrode.device.name not in nwbfile.devices.keys():
                    device = nwbfile.create_device(name=swb.electrode.device.name,
                    description=swb.electrode.device.description, manufacturer=swb.electrode.device.manufacturer)
                else:
                    device = nwbfile.devices[swb.electrode.device.name]
                swb.fields['electrode'] = nwbfile.create_icephys_electrode(
                        name=swb.electrode.name, description=swb.electrode.description, device=device
                    )
            else:
                swb.fields['electrode'] = nwbfile.ic_electrodes[swb.electrode.name]
            if isinstance(swb, CurrentClampStimulusSeries) or isinstance(swb, VoltageClampStimulusSeries):
                nwbfile.add_stimulus(swb, use_sweep_table=use_sweep_table)
            else:
                nwbfile.add_acquisition(swb, use_sweep_table=use_sweep_table)
        return nwbfile

    @classmethod
    def reset_container(cls, hdmf_obj, parent:NWBFile)->None:
        ''' reset container fields to make sure the object can be added to a new container
        Parameters
        ----------
        hdmf_obj: hdmf object.
        parent: hdmf container type. parent container
        '''
        hdmf_obj._AbstractContainer__container_source = None
        hdmf_obj._in_construct_mode = True
        hdmf_obj.set_modified()
        hdmf_obj.reset_parent()
        hdmf_obj.parent = parent
     
    @classmethod
    def copy_sweeps(cls, src:NWBFile, dst: Union[NWBFile, None], sweep_table_index: list=[]):
        ''' Copy sweeps from source NWB file to destination NWB file.
        This can be done, but it can not be saved. the object ID is 
        '''
        assert isinstance(src, NWBFile), 'src must be an NWBFile object'
        if dst is None:
            dst = pvNWB.make_new_NWBcontainer(src=src)
        dst._check_sweep_table() ## make sure sweep table is created
        swTab_ = src.sweep_table.copy() ## make a copy of the scr sweep table
        cls.reset_container(swTab_, dst) # reset container source and parent
        for s in sweep_table_index: 
            sweep_series = swTab_.get_series(s)
            dst = pvNWB.add_nwb_sweep(dst, sweep_series) 
        return dst


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
        return axes

def add_attributes(src_object, attrNames, init_dict):
    ''' add attributes from src_object to init_dict'''
    for a in attrNames:
        if hasattr(src_object, a):
            init_dict.update({a:getattr(src_object, a)})
    return init_dict

if __name__ == "__main__":
    testFile = "D:\\DataDemo\\nwb\\sub-626194732_ses-638205339_icephys.nwb"
    x = pvNWB(testFile)  # read NWB file

    print(f'num sweeps: {x.totalNumOfSweeps}')
    nwbfile = x.nwbfile # get NWB file object

    create_nwb3 = True
    if create_nwb3: # copy a set of sweeps to a new NWB file
        sweep_table_index = [0,1,5]
        nwb2 = pvNWB.make_new_NWBcontainer(src=x.nwbfile)
        nwb2 = pvNWB.copy_sweeps(x.nwbfile, nwb2, sweep_table_index)
        print('nwb2 sweeps: ', pvNWB.getNumOfSweeps(nwb2))

        with NWBHDF5IO('copy_nwb2.nwb', 'w') as w_io:
            w_io.write(nwb2)

        print('\nready to read file')
        with NWBHDF5IO('copy_nwb2.nwb', 'r', load_namespaces=True) as io:
            print('io created')
            nwb3 = io.read() # read NWB file
            print('nwb2 sweeps: ', pvNWB.getNumOfSweeps(nwb2))
