from patchview.utilitis.pvEphy import pvNWB
from patchview.utilitis.pvNDX_class import PatchviewSweep, PatchviewSweepGroup
import numpy as np

class dandiNWB(pvNWB):
    def __init__(self, nwbFile, **kw):
        """
        A custom class for Patch-seq data by Gouwens et al.
        https://dandiarchive.org/dandiset/000020/ 
        Notes: their NWB files' sweep table is not compatible with current version of pynwb
        For Heka dat files. Change existing procotol dict may be an easier way.
        For other file, inherite the pvNWB class and write custom loader and sweep group parser.
        """
        super(dandiNWB, self).__init__(nwbFile,
            **kw
        )
    def __post_init__(self):
        self.defineSweepGroups()
        self.addPachviewSweepGroup(self.sweepGroups) # list of PatchviewSweepGroup

    def defineSweepGroups(self) -> list[PatchviewSweep]:
        ''' Dandi Allen Patchseq NWB file sweep groups
        '''
        fields = ['VoltageClampSeries', 'CurrentClampSeries', 'Custom'] # make three fields in sweep group metainformation
        clampData = {}  # current clamp
        current_stimType = {'C1LSCOARSE150216':'LongSquareWave',
            'C1LSCOARSEMICRO':'LongSquareWave_pre',
        'C1RP25PR1S141203':'Ramp', 'C1SSCOARSE150112':'ShortPulse',
            'C1LSFINEST150112':'LongSquare_post'}
        for t in current_stimType.keys():
            clampData[current_stimType[t]]  = [] # add key for each current clamp stim type
        vclampData = {'vclamp':[]} # voltage clamp. to be extended...
        customData = {}
        # parse sweep table
        for s in range(self.totalNumOfSweeps):
            recording, sti = self.nwbfile.sweep_table.get_series(s)
            rtype = recording.neurodata_type
            if  rtype== 'VoltageClampSeries':
                vclampData['vclamp'].append(s)
            elif rtype == 'CurrentClampSeries':
                try:
                    stiDec = sti.stimulus_description.split('_')[0]            
                except:
                    stiDec = 'Custom'
                if stiDec in current_stimType.keys():
                    stimType  = current_stimType[stiDec]
                clampData[stimType].append(s)               
            else:
                print(f'unknown data type:{rtype}')
                customData['custom'].append(s)
        swGroups = []
        for container in [clampData,vclampData, customData]:
            for stimType in  container:
                pv_sweep_group = PatchviewSweepGroup(name=stimType)
                for idx, swIdx in enumerate(container[stimType]): # per sweep
                    pv_sweep = PatchviewSweep(name=f'{idx:03}', stimulus_type=stimType, 
                    stimulus_name=stimType, sweep_index = swIdx, source_sweep_table_index = swIdx,
                        trace_indexes = np.array([0]), trace_labels=['trace1'])
                    pv_sweep_group.add_patchview_sweep(pv_sweep)
                swGroups.append(pv_sweep_group)
        self.sweepGroups = swGroups
        return swGroups

if __name__ == "__main__":
    nwbfn = "D:\\DataDemo\\nwb\\sub-626194732_ses-638205339_icephys.nwb"
    x = dandiNWB(nwbfn)
    print(f'num sweeps: {x.totalNumOfSweeps}')
    print(x.sweepGroups)
    # x.plotSweep([10])