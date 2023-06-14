from pynwb.spec import NWBNamespaceBuilder, NWBGroupSpec
from pynwb.spec import NWBDatasetSpec, NWBAttributeSpec
from collections import namedtuple
import os
utilitis_dir, this_filename = os.path.split(__file__)

labName = 'XiaolongJiangLab'
ns_path = labName +".namespace.yaml"
ext_source = labName +".extensions.yaml"

ns_builder = NWBNamespaceBuilder('Extension for use in Xiaolong Jiang Lab', labName, version='0.1.0')
ns_builder.include_type('NWBDataInterface', namespace='core')

pvSweep = NWBGroupSpec(neurodata_type_def='PatchviewSweep',
                      neurodata_type_inc='NWBDataInterface',
                      doc='A sweep with one or multiple traces', quantity='*',
                      datasets=[
                        NWBDatasetSpec(doc='index for all traces/channels ', shape=(None, 4), # may have traces from multiple electrodes
                                        name='trace_indexes', dtype='int', dims=('trace number', 'trace index')),
                        NWBDatasetSpec(doc='label for all traces ', shape=(None, ), # may have traces from multiple electrodes
                        name='trace_labels', dtype='text'),
                      ],
                      attributes=[ #
                        NWBAttributeSpec(name='sweep_table_index',
                                        doc='index in the NWB sweep table. Use this to extract the sweep from NWB file', 
                                        dtype='int',
                                        required=False),
                        NWBAttributeSpec(name='local_sweep_index', 
                                        doc='local relative index in current sweep group', 
                                        dtype='int',
                                        required=False),
                        NWBAttributeSpec(name='stimulating_electrode_index', # read from the electrode table
                                        doc='index for the electrode that applies stimulus', 
                                        dtype='int',
                                        required=False),
                        NWBAttributeSpec(name='stimulus_type', # stimulus type, e.g. 'sine', 'square', 'ramp', 'custom'
                                        doc='stimulus type. Facilitate auto analysis',
                                        dtype='text',
                                        required=False),
                        NWBAttributeSpec(name='stimulus_onset',
                                        doc='stimulus onset time in seconds',
                                        dtype='float',
                                        required=False),
                        NWBAttributeSpec(name='stimulus_offset',
                                        doc='stimulus offset time in seconds',
                                        dtype='float',
                                        required=False),
                        NWBAttributeSpec(name='stimulus_amplitude',
                                        doc='stimulus amplitude',
                                        dtype='float',
                                        required=False),
                        NWBAttributeSpec(name='stimulus_holding_amplitude',
                                        doc='holding current/voltage amplitude',
                                        dtype='float', # 0=constant, 1=variable
                                        required=False),                                        
                        NWBAttributeSpec(name='stimulus_name',
                                        doc='custom human-readable stimulus name', # usually the label of the stimulus by experimenter
                                        dtype='text',
                                        required=False),
                      ])
#  a series (sweep group) contain multiple sweeps
# Note that if youspecify name, quantity cannot be '*', '+', or an integer greater that 1,
#  because you cannot have more than one group of the same name in the same parent group.
#https://pynwb.readthedocs.io/en/stable/pynwb.spec.html#pynwb.spec.NWBGroupSpec.add_group
pvSweepGroup = NWBGroupSpec(neurodata_type_def='PatchviewSweepGroup',
                      neurodata_type_inc='NWBDataInterface',
                      doc='Sweep group for Patchview, parent of pvSweep',
                      quantity = '?',
                      groups = [pvSweep])

#  a series group contain multiple sweep groups. this is useful for connection series
# pvSeriesGroup = NWBGroupSpec(neurodata_type_def='PatchviewSeriesGroup',
#                       neurodata_type_inc='NWBDataInterface',
#                       doc='Series group for Patchview. Parent of pvSweepGroup',
#                       quantity='?',groups = [pvSweepGroup])

ns_builder.add_spec(ext_source, pvSweepGroup)
ns_builder.export(ns_path)


''' NOTES
a FP series can be a `inc` type after `CurrentClamp` type.
A connection series can be a another inc type `CurrentClamp` type,
with at least two simultaneous recording electrodes.

patchviewGroup defined above behaviour similar to `region references` in NWB.
"region references are similar to object references,
but in contrast point to regions (i.e., select subsets) of
datasets with an assigned type." It can be thought of as a more concrete type of 
any sweep-based data type, e.g. `CurrentClampSeries`, `VoltageClampSeries` with 
augmented attributes pointing to the sweep table.

The connection type may be the most useful. Analysis based on this type can be
done in a more general way, e.g. `connection strength`, `connection latency`.

Alternatively, each connection test is still saved in a seperate container or NWB file.
But each of them maintains a reference to all other connection tests. This is similar
to the `region references` in NWB. The advantage is that the connection test can be
saved in a more compact way, e.g. only the connection test data is saved, not the
whole sweep data. The disadvantage is that the connection test data is not as
convenient to access as the `region references` in NWB. 
The above features seems to in dynamic table, especially in `SimultaneousRecordingsTable`,
` SequentialRecordingsTable`, `RepetitionsTable`,`ExperimentalConditionsTable`
https://nwb-schema.readthedocs.io/en/latest/format.html#sec-simultaneousrecordingstable
Need to dig into the NWB dynamic table more.


SweepTable is being replaced by IntracellularRecordingsTable and SimultaneousRecordingsTable tables.
Additional SequentialRecordingsTable, RepetitionsTable, and ExperimentalConditions tables
provide enhanced support for experiment metadata.

'''