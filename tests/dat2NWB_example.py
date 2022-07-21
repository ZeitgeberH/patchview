# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 20:47:51 2020

@author: Ming Hu (ming.hu@bcm.edu)
"""
from pynwb import NWBHDF5IO
from patchview.utilitis.PVdat2NWB import dat2NWB
import matplotlib.pyplot as plt
import numpy as np
import os
import pdb
#%% Write NWB file from dat file
testFile = "210917s3c04.dat" 
nwbData = dat2NWB(os.path.join('Data',testFile), [0,0])

saveNwbFileName = 'test2.nwb'
nwbData.saveNBF(os.path.join('Data',saveNwbFileName))

print(f'Number of sweeps: {nwbData.getNumberOfSweeps()}')

stim, resp = nwbData.getSweep(0)
nwbData.plotSweep([0, 5, 20])

#%% Check if saved NWB file has the right cotents.
#with NWBHDF5IO('200225s5c22.nwb', 'r') as io:
io = NWBHDF5IO(os.path.join('Data',saveNwbFileName), 'r') 
nwbfile = io.read()

print(nwbfile.fields.keys())

#%% electrode name
electrodes = list(nwbfile.icephys_electrodes.keys())
stimElectrodeName = [x for x in electrodes if 'ccs' in x]
voltageElectrodeName = [y for y in electrodes if 'Vmon-' in y]
print(electrodes, stimElectrodeName,voltageElectrodeName)
### get more information about an electrode by: nwbfile.get_icephys_electrode('ccsElectrode0')
#%%
#%%
def plot_sweep(sweep,ax=None):
    if ax is None:
        _, ax = plt.subplots()
    dat = sweep.data[:]
    yy = dat *sweep.conversion
    xx = np.arange(len(dat))/sweep.rate    
    ax.plot(xx, yy)   

nSweeps = len(nwbfile.sweep_table)//2
fig, axs = plt.subplots(2,1, sharex=True)
for i in range(nSweeps):
    stimulus, response = nwbfile.sweep_table.get_series(i)
    plot_sweep(stimulus,  ax=axs[1])
    plot_sweep(response, ax=axs[0])
axs[0].set_xlabel('')
axs[0].set_ylabel(response.unit)
axs[1].set_ylabel(stimulus.unit)
axs[1].set_xlabel('time (s)')
plt.show()
#%%
io.close()