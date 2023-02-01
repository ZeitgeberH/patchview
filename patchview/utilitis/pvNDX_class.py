from pynwb import register_class, load_namespaces
from pynwb.file import MultiContainerInterface, NWBContainer
import os
import numpy as np
from hdmf.utils import docval,get_docval,popargs_to_dict
## load namespace
utilities_dir, this_filename = os.path.split(__file__)

labName = "XiaolongJiangLab"
ns_path = os.path.join(utilities_dir,"NDX_files", labName+".namespace.yaml") # check patchview_ndx.py for details
load_namespaces(ns_path)

@register_class('PatchviewSweep', labName)
class PatchviewSweep(NWBContainer):
    __nwbfields__ = ('name', 'stimulating_electrode_index', 'stimulus_type', 'stimulus_name', 'sweep_indexes')
    @docval({'name': 'name', 'type': str, 'doc': 'name of sweep in this sweep'},
            {'name':'source_sweep_table_index', 'type': int, 'doc': 'sweep number index in original source sweep table'},
            {'name':'sweep_index','type': int, 'doc': 'local index in sweep'},
            {'name': 'stimulating_electrode_index', 'type': int, 'doc': 'index for electrode that apply stimuli', 'default':0},
            {'name': 'stimulus_type', 'type': str, 'doc': 'stimulus type used by Patchview', 'default':''}, 
            {'name': 'stimulus_name', 'type': str, 'doc': 'stimulus name labelled by experimenter','default':''},
            {'name':'trace_indexes', 'type': None, 'doc': 'indexes for all traces'},
            {'name':'trace_labels', 'type': None, 'doc': 'label for all traces', 'default':[]})

    def __init__(self, **kwargs):
        super().__init__(name=kwargs['name'])
        self.stimulating_electrode_index = kwargs['stimulating_electrode_index']
        self.stimulus_type = kwargs['stimulus_type']
        self.stimulus_name = kwargs['stimulus_name']
        self.trace_indexes = kwargs['trace_indexes'] # original index from whatever format (DAT, abf. etc.)
        self.source_sweep_table_index = kwargs['source_sweep_table_index'] # from original source NWB table if any
        self.sweep_index = kwargs['sweep_index'] # current sweep table. This is what matters if history is not concern
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

# @register_class('PatchviewSeriesGroup', labName)
# class PatchviewSeriesGroup(MultiContainerInterface):
#     __clsconf__ = {
#         'attr': 'patchview_sweep_groups',# this attr name is tricky!
#         'type': PatchviewSweepGroup,
#         'add': 'add_patchview_sweep_group',
#         'get': 'get_patchview_sweep_group',
#         'create': 'create_patchview_sweep_group',
#     }

if __name__=='__main__':
    from pynwb import NWBHDF5IO, NWBFile
    from datetime import datetime
    from dateutil.tz import tzlocal
    ns_path = os.path.join(utilities_dir,"NDX_files", labName+".namespace.yaml") # check patchview_ndx.py for details
    load_namespaces(ns_path)
    # series_group = PatchviewSeriesGroup(name='testing session')

    sweep_indexes1 = np.array([[0,0,0,1], [0,0,0,2],[0,0,0,3]]) # [group, session, sweep, trace]
    sweep1 = PatchviewSweep(name='sweep_sine', stimulus_type='sine', sweep_table_index = 0,
     stimulus_name='sine1', trace_indexes=sweep_indexes1, stimulating_electrode_index=0)

    sweep_indexes2 = np.array([[0,0,1,1], [0,0,1,2]]) 
    sweep2 = PatchviewSweep(name='sweep_square', stimulus_type='square',stimulus_name='square',
     trace_labels = ['recording', 'stimulus'],stimulating_electrode_index=1,sweep_table_index = 1,
    trace_indexes = sweep_indexes2)
    sweep_group = PatchviewSweepGroup(name='testing sweep group')
    sweep_group.add_patchview_sweep(sweep1)
    sweep_group.add_patchview_sweep(sweep2)
    print(sweep_group)

    sweep_group2 = PatchviewSweepGroup(name='testing sweep group2')
    sweep_indexes2b = np.array([[0,0,0,1], [0,0,0,2],[0,0,0,3]]) # [group, session, sweep, trace]
    sweep2b = PatchviewSweep(name='sweep_sine', stimulus_type='sine', sweep_table_index = 0,
     stimulus_name='sine1', trace_indexes=sweep_indexes2b, stimulating_electrode_index=0)
    sweep_group2.add_patchview_sweep(sweep2b)

    # series_group.add_patchview_sweep_group(sweep_group)
    # series_group.add_patchview_sweep_group(sweep_group2)

    nwbfile = NWBFile("Pachview_NWB", "101", datetime(2023, 1, 27, tzinfo=tzlocal()))

    pmod = nwbfile.create_processing_module('patchview series group', 'series group 1') # name, description
    pmod.add([sweep_group,sweep_group2])

    # pmod2 = nwbfile.create_processing_module('patchview sweep group2', 'sweep group 2')
    # pmod2.add_container(sweep_group2)

    with NWBHDF5IO('test_multicontainerinterface4.nwb', 'w') as io:
        io.write(nwbfile)

    # test missing trace_labels
    # sweep_indexes3 = np.array([[0,0,0,1], [0,0,0,2]]) # [group, session, sweep, trace]
    # sweep3 = SweepIndexes(name='sweep_sine', stimulus_type='sine',
    #  stimulus_name='sine1', trace_indexes=sweep_indexes1,
    #  trace_labels = ['recording'], # test missing trace_labels
    # stimulating_electrode_index=0)