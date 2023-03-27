from patchview.utilitis.pvEphy import pvNWB
from patchview.utilitis.PVdat2NWB import dat2NWB
import os
_FILEPATH  = "D:\\DataDemo"
# datFn = "D:\\DataDemo\\190628s1c6fp.dat"
# fname= "D:\\DataDemo\\190628s1c6fp_sel99.nwb"
fname = "D:\\DataDemo\\nwb\\sub-626194732_ses-638205339_icephys.nwb"
# fname = "D:\\DataDemo\\test_multicontainerinterface.nwb"
# fname = os.path.join(_FILEPATH, '190628s1c6fp.nwb')
# fname = os.path.join(_FILEPATH, 'test_singleFP.nwb')
# fname = "D:\\DataDemo\\nwb\\sub-626194732_ses-638205339_icephys.nwb"

# dat_nwb_obj = dat2NWB(datFn, [0,0], fname, protocol='Patchview')
# nwbfile = dat_nwb_obj.nwbfile
# dat_nwb_obj.saveNWB(fname)

# fname2 = "D:\\DataDemo\\test_2.nwb"
x = pvNWB(fname)
nwbfile = x.nwbfile
# promods = nwbfile.get_processing_module() # need a name if more than one processing module
promods = nwbfile.processing

pv_sweepGroups = x.getPachviewSweepGroup()
# pv_sweepGroups = []
# for modName in promods: # assume only one sweep group exist
#     for c in promods[modName].children:
#         if c.neurodata_type == 'PatchviewSweepGroup':
#             pv_sweepGroups.append(c)

for pvg in pv_sweepGroups: ## may have mulitple groups
    print(pvg.name) ## have type of clamp
    pvGroup = pvg.patchview_sweeps # ->'hdmf.utils.LabelledDict'

# seriesData = nwbfile.sweep_table.get_series(0)
stim, data, metaInfo = x.getSweepData(0)
# len(seriesData)
# stimulus = nwbfile.stimulus
# sti_keys = list(stimulus.keys())
# stimulus[sti_keys[0]]
# # nChildren = len(promods.children)
# # pv_sweepGroup = promods.children[0].patchview_sweeps
# # swp_names = sorted(pv_sweepGroup.keys())
# # sweep0 = pv_sweepGroup[swp_names[0]]
# # trace_index = sweep0.trace_indexes
# # seriesData = nwbfile.sweep_table.get_series(0)
# stim, data, metaInfo = x.getSweepData(0)
# x.io.close()
# print('meta info:', x.metaInfo)

# stim, data, metaInfo = x.getSweepData(0)

# fname =  os.path.join(_FILEPATH, '190628','190628s1c6fp.dat')
# pvEphy = dat2NWB(fname, [0,0])


