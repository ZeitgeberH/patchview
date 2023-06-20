from patchview.utilitis.pvEphy import pvNWB
from patchview.utilitis.pvNDX_class import PatchviewSweep, PatchviewSweepGroup
import numpy as np

class dandiNWB(pvNWB):
    def __init__(self, nwbFile, **kw):
        """
        A custom class for Patch-seq data by Gouwens et al.
        https://dandiarchive.org/dandiset/000020/ 
        Notes: their NWB files require NDX-MIES extension. If installed (https://github.com/t-b/ndx-MIES)
        this class may not be necessary.
        """
        super(dandiNWB, self).__init__(nwbFile,
            **kw
        )
    def __post_init__(self):
        self.defineSweepGroups()
        self.addPachviewSweepGroup(self.sweepGroups) # list of PatchviewSweepGroup

    def defineSweepGroups(self) -> dict[str, PatchviewSweepGroup]:
        ''' Dandi Allen Patchseq NWB file sweep groups
        '''
        fields = ['VoltageClampSeries', 'CurrentClampSeries', 'Custom'] # make three fields in sweep group metainformation
        clampData = {}  # current clamp
        current_stimType = {'C1LSCOARSE150216':'LongSquareWave',
                            'X1PS':'LongSquareWave',
                            # 'X2LP':'LongSquareWave',
                            'X3LP':'LongSquareWave',
                            'C1LSCOARSE150216_DA_0':'LongSquareWave',
                            'subthreshold_171_DA_0':'LongSquareWave',
                            'C1LSFINEST150112_DA_0':'LongSquareWave',
            'C1LSCOARSEMICRO':'LongSquareWave_pre',
        'C1RP25PR1S141203':'Ramp','X7Ramp':'Ramp',
          'C1SSCOARSE150112':'ShortPulse','X5SP':'ShortPulse',
            'C1LSFINEST150112':'LongSquare_post', 'C2SSTRIPLE171103':'PulseSequence'}
        for t in current_stimType.keys():
            clampData[current_stimType[t]]  = [] # add key for each current clamp stim type
        vclampData = {'vclamp':[]} # voltage clamp. to be extended...
        customData = {}
        # parse sweep table
        minCurrent = 1e10
        maxCurrent = -1e10
        print('total sweep:', self.totalNumOfSweeps)
        for s in range(self.totalNumOfSweeps):
            recording, sti = self.nwbfile.sweep_table.get_series(s)
            rtype = recording.neurodata_type      
            stimulus_amplitude = 0.0
            if  rtype== 'VoltageClampSeries':
                vclampData['vclamp'].append(s)
            elif rtype == 'CurrentClampSeries':
                try:
                    stiDec = sti.stimulus_description.split('_')[0]
                    # print(s, sti.stimulus_description.split('_'))       
                except:
                    stiDec = 'Custom'
                if stiDec in current_stimType.keys():
                    stimType  = current_stimType[stiDec]
        
                    if stimType == 'LongSquareWave': ## assume constant step current
                        minCurrent = np.min([minCurrent, np.min(sti.data[2000:-2000])])
                        maxCurrent = np.max([maxCurrent, np.max(sti.data[2000:-2000])])
     
                    
                    clampData[stimType].append(s)               
            else:
                print(f'unknown data type:{rtype}')
                customData['custom'].append(s)

        print(f'long square wave min current: {minCurrent}, max current: {maxCurrent}')
        swGroups = {}
        for container in [clampData,vclampData, customData]:
            for stimType in  container:
                pv_sweep_group = PatchviewSweepGroup(name=stimType)
                for idx, swIdx in enumerate(container[stimType]): # per sweep
                    stimulus_onset, stimulus_offset, stimulus_amplitude = 0.0, 0.0, 0.0
                    if stimType == 'LongSquareWave' and len(container[stimType]) > 0:
                        stimulus_amplitude = minCurrent+ idx*(maxCurrent-minCurrent)/len(container[stimType])
                        stimulus_onset, stimulus_offset = 1.0, 2.0
                    pv_sweep = PatchviewSweep(name=f'{idx:02}', stimulus_type=stimType, 
                    stimulus_name=stimType, local_sweep_index = idx, sweep_table_index = swIdx,
                    stimulus_offset = stimulus_offset, stimulus_onset = stimulus_onset,
                    stimulus_amplitude = stimulus_amplitude, stimulus_holding_amplitude = 0.0,
                        trace_indexes = np.array([0]), trace_labels=['trace1'])
                    pv_sweep_group.add_patchview_sweep(pv_sweep)
                swGroups.update({stimType:pv_sweep_group})
        self.sweepGroups = swGroups
        return swGroups

if __name__ == "__main__":
    nwbfn = "D:\\DataDemo\\nwb\\sub-626194732_ses-638205339_icephys.nwb"
    x = dandiNWB(nwbfn)
    print(f'num sweeps: {x.totalNumOfSweeps}')
    swpGrps = x.getSweepGroups()
    print(f'num sweep groups: {len(swpGrps)}')
    # x.plotSweep([10])