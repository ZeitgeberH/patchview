'''
Example of how to load the spike list output by Patchview.
And plot the spike waveforms of the first sweep of each neuron
'''
import pickle
import matplotlib.pyplot as plt

dataName = 'Aligned_spike_list3.pickle'
# load the spike list
with open(dataName, 'rb') as f:
    spike_list = pickle.load(f)

# get the names of the neurons
neuronNames = list(spike_list.keys())

# an example of how to plot the spike waveforms of the first neuron
firstSpikeOnly = True
plt.figure()
for neuron in neuronNames: # get the name of neuron
    sweep_names = sorted(list(spike_list[neuron].keys())) # get the names of the sweeps
    print(f'Cell {neuron} has {len(sweep_names)} sweeps')
    if len(sweep_names) > 0:
        
        ## plot the spike waveforms of this neuron's first sweep above rehbose
        sw = sweep_names[0]
        if firstSpikeOnly:
            plt.plot(spike_list[neuron][sw][:,0])
        else:
            plt.plot(spike_list[neuron][sw])
        
        ## all sweeps
        # for sw in sweep_names:
        #     if firstSpikeOnly:
        #         plt.plot(spike_list[neuron][sw][:,0])
        #     else:
        #         plt.plot(spike_list[neuron][sw])

        plt.xlabel('Time (ms)')
        plt.ylabel('Voltage (mV)')
        plt.title(f'Spike waveforms of {neuron}')
plt.show()