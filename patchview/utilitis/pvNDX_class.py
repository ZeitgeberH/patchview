from pynwb import register_class, load_namespaces
from pynwb.file import MultiContainerInterface, NWBContainer
import os
import numpy as np
from hdmf.utils import docval,get_docval,popargs_to_dict
import functools
## load namespace
utilities_dir, this_filename = os.path.split(__file__)

labName = "XiaolongJiangLab"
ns_path = os.path.join(utilities_dir,"NDX_files", labName+".namespace.yaml") # check patchview_ndx.py for details
load_namespaces(ns_path)

@register_class('PatchviewSweep', labName)
class PatchviewSweep(NWBContainer):
    __nwbfields__ = ('name', 'stimulating_electrode_index', 'stimulus_type', 'stimulus_name', 'sweep_indexes')
    @docval({'name': 'name', 'type': str, 'doc': 'name of sweep in this sweep'},
            {'name':'sweep_table_index', 'type': int, 'doc': 'sweep number index in NWB source sweep table'},
            {'name':'local_sweep_index','type': int, 'doc': 'local index relatively in this sweep group'},
            {'name': 'stimulating_electrode_index', 'type': int, 'doc': 'index for electrode that apply stimuli', 'default':0},
            {'name': 'stimulus_type', 'type': str, 'doc': 'stimulus type used by Patchview', 'default':''}, 
            {'name': 'stimulus_onset', 'type': float, 'doc': 'stimulus onset time in seconds', 'default': 0.1},
            {'name': 'stimulus_offset', 'type': float, 'doc': 'stimulus offset time in seconds', 'default': 0.7},
            {'name': 'stimulus_amplitude', 'type': float, 'doc': 'stimulus amplitude. Relevent if constant ', 'default': 0.0},
            {'name': 'stimulus_holding_amplitude', 'type': float, 'doc': 'stimulus amplitude during baseline', 'default': 0.0},
            {'name': 'stimulus_name', 'type': str, 'doc': 'stimulus name labelled by experimenter','default':''},
            {'name':'trace_indexes', 'type': None, 'doc': 'indexes for all traces'},
            {'name':'trace_labels', 'type': None, 'doc': 'label for all traces', 'default':[]})

    def __init__(self, **kwargs):
        super().__init__(name=kwargs['name'])
        self.stimulating_electrode_index = kwargs['stimulating_electrode_index']
        self.stimulus_type = kwargs['stimulus_type']
        self.stimulus_name = kwargs['stimulus_name']
        self.stimulus_onset = kwargs['stimulus_onset']
        self.stimulus_offset = kwargs['stimulus_offset']
        self.stimulus_amplitude = kwargs['stimulus_amplitude']
        self.stimulus_holding_amplitude = kwargs['stimulus_holding_amplitude']
        self.trace_indexes = kwargs['trace_indexes'] # original index from whatever format (DAT, abf. etc.)
        self.sweep_table_index = kwargs['sweep_table_index'] # from  NWB sweep table 
        self.local_sweep_index = kwargs['local_sweep_index'] # local index in the current sweep group.
        if len(kwargs['trace_labels']) == 0:
            self.trace_labels = ['trace_'+str(i) for i in range(self.trace_indexes.shape[0])]
        else:
            assert len(kwargs['trace_labels']) == self.trace_indexes.shape[0],\
                 'trace_labels should have the same length as trace_indexes. Do not set it if you do not know what it is.'
            self.trace_labels = kwargs['trace_labels']

@register_class('PatchviewSweepGroup', labName)
class PatchviewSweepGroup(MultiContainerInterface):
    __clsconf__ = {
        'attr': 'patchview_sweeps',# this attr name is tricky!
        'type': PatchviewSweep,
        'add': 'add_patchview_sweep',
        'get': 'get_patchview_sweep',
        'create': 'create_patchview_sweep',
    }

    def __len__(self):
        return len(self.patchview_sweeps)
    
    @functools.cached_property
    def sweepNames(self):
        return list(self.patchview_sweeps.keys())

    def get_sweep(self, sweep_index:int) -> PatchviewSweep:
        '''get sweep by index. patchviewSweep need to be ordered dict'''
        return self.get_patchview_sweep(self.sweepNames[sweep_index])
    
    def get_sweep_table_index(self, sweep_index):
        return self.get_sweep(sweep_index).sweep_table_index

    def get_stimulus_type(self, sweep_index):
        return self.get_sweep(sweep_index).stimulus_type

    def get_stimulus_name(self, sweep_index):
        return self.get_sweep(sweep_index).stimulus_name

    def get_stimulus_onset(self, sweep_index):
        return self.get_sweep(sweep_index).stimulus_onset

    def get_stimulus_offset(self, sweep_index):
        return self.get_sweep(sweep_index).stimulus_offset

    def get_stimulus_amplitude(self, sweep_index):
        return self.get_sweep(sweep_index).stimulus_amplitude
 
    def get_stimulus_holding_amplitude(self, sweep_index):
        return self.get_sweep(sweep_index).stimulus_holding_amplitude
 
    def get_trace_indexes(self, sweep_index):
        return self.get_sweep(sweep_index).trace_indexes

    def get_trace_labels(self, sweep_index):
        return self.get_sweep(sweep_index).trace_labels


if __name__=='__main__':
    from pynwb import NWBHDF5IO, NWBFile
    from datetime import datetime
    from dateutil.tz import tzlocal
    ns_path = os.path.join(utilities_dir,"NDX_files", labName+".namespace.yaml") # check patchview_ndx.py for details
    load_namespaces(ns_path)
    # series_group = PatchviewSeriesGroup(name='testing session')

    sweep_indexes1 = np.array([[0,0,0,1], [0,0,0,2],[0,0,0,3]]) # [group, session, sweep, trace]
    sweep1 = PatchviewSweep(name='sweep_sine', stimulus_type='sine', local_sweep_index=0, sweep_table_index = 4,
     stimulus_name='sine1', trace_indexes=sweep_indexes1, stimulating_electrode_index=0)

    sweep_indexes2 = np.array([[0,0,1,1], [0,0,1,2]]) 
    sweep2 = PatchviewSweep(name='sweep_square', stimulus_type='square',stimulus_name='square', stimulus_onset = 0.1, 
    stimulus_offset=0.7, stimulus_amplitude = 1.0e-12,
     trace_labels = ['recording', 'stimulus'],stimulating_electrode_index=1, local_sweep_index=1, sweep_table_index = 6,
    trace_indexes = sweep_indexes2)
    sweep_group = PatchviewSweepGroup(name='testing sweep group')
    sweep_group.add_patchview_sweep(sweep1)
    sweep_group.add_patchview_sweep(sweep2)
    print(sweep_group)
    print('number of sweeps', len(sweep_group))
    print(sweep_group.get_sweep(0))
    print(sweep_group.get_sweep_table_index(0))
    print(sweep_group.get_stimulus_amplitude(0))
    sweep_group2 = PatchviewSweepGroup(name='testing sweep group2')
    sweep_indexes2b = np.array([[0,0,0,1], [0,0,0,2],[0,0,0,3]]) # [group, session, sweep, trace]
    sweep2b = PatchviewSweep(name='sweep_sine', stimulus_type='sine', local_sweep_index = 0, sweep_table_index = 0,
     stimulus_name='sine1', trace_indexes=sweep_indexes2b, stimulating_electrode_index=0)
    sweep_group2.add_patchview_sweep(sweep2b)

    nwbfile = NWBFile("Pachview_NWB", "101", datetime(2023, 1, 27, tzinfo=tzlocal()))

    pmod = nwbfile.create_processing_module('patchview series group', 'series group 1') # name, description
    pmod.add([sweep_group,sweep_group2])


    with NWBHDF5IO('test_patchviewGroup.nwb', 'w') as io:
        io.write(nwbfile)

    # testing missing trace_labels
    # sweep_indexes3 = np.array([[0,0,0,1], [0,0,0,2]]) # [group, session, sweep, trace]
    # sweep3 = SweepIndexes(name='sweep_sine', stimulus_type='sine',
    #  stimulus_name='sine1', trace_indexes=sweep_indexes1,
    #  trace_labels = ['recording'], # test missing trace_labels
    # stimulating_electrode_index=0)