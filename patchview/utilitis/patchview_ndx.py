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
                        NWBDatasetSpec(doc='index for all traces ', shape=(None, 4), # may have traces from multiple electrodes
                                        name='trace_indexes', dtype='int', dims=('trace number', 'trace index')),
                        NWBDatasetSpec(doc='label for all traces ', shape=(None, ), # may have traces from multiple electrodes
                        name='trace_labels', dtype='text'),
                      ],
                      attributes=[ #
                        NWBAttributeSpec(name='source_sweep_table_index', # read from the electrode table
                                        doc='index in the source NWB sweep table', 
                                        dtype='int',
                                        required=False),
                        NWBAttributeSpec(name='sweep_index', # read from the electrode table
                                        doc='local relative index in current sweep group', 
                                        dtype='int',
                                        required=False),
                        NWBAttributeSpec(name='stimulating_electrode_index', # read from the electrode table
                                        doc='index for the electrode that applies stimulus', 
                                        dtype='int',
                                        required=False),
                        NWBAttributeSpec(name='stimulus_type', # stimulus type, e.g. 'sine', 'square', 'ramp', 'custom'
                                        doc='stimulus type',
                                        dtype='text',
                                        required=False),
                        NWBAttributeSpec(name='stimulus_name',
                                        doc='custom stimulus name', # usually the label of the stimulus by experimenter
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

#  a series group contain multiple sweep groups
# pvSeriesGroup = NWBGroupSpec(neurodata_type_def='PatchviewSeriesGroup',
#                       neurodata_type_inc='NWBDataInterface',
#                       doc='Series group for Patchview. Parent of pvSweepGroup',
#                       quantity='?',groups = [pvSweepGroup])

ns_builder.add_spec(ext_source, pvSweepGroup)
ns_builder.export(ns_path)