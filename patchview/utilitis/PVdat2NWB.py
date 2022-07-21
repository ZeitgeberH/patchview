# -*- coding: utf-8 -*-
"""
https://pynwb.readthedocs.io/en/stable/tutorials/domain/icephys.html#sphx-glr-tutorials-domain-icephys-py
ming.hu@bcm.edu
"""
from datetime import datetime
from dateutil.tz import tzlocal
from pynwb import NWBFile, NWBHDF5IO
from pynwb.icephys import CurrentClampStimulusSeries, CurrentClampSeries
import numpy as np
from patchview.HekaIO import HEKA_Reader_MAIN as HEKA
import matplotlib.pyplot as plt

class dat2NWB(NWBFile):
    def __init__(self, datFile, sessionIdx, **kw):
        """
        parameters:sessionIdx
        dataFile: dat filename
        sessionIdx: [group, series]  ## one file for one series
        """
        self._HekaBundle(datFile)
        super(dat2NWB, self).__init__(
            "",
            datFile,
            self._bundle.session_start_time,
            file_create_date=datetime.now(tzlocal()),
            institution="Baylor College of Medicine",
            lab="Xiaolong Jiang lab",
            **kw
        )
        """
        (session_description, identifier, session_start_time, 
        file_create_date=None, timestamps_reference_time=None,
        experimenter=None, experiment_description=None, session_id=None,
        institution=None, keywords=None, notes=None, pharmacology=None,
        protocol=None, related_publications=None, slices=None,
        source_script=None, source_script_file_name=None, 
        data_collection=None, surgery=None, virus=None,
        stimulus_notes=None, lab=None, acquisition=None,
        analysis=None, stimulus=None, stimulus_template=None, 
        epochs=None, epoch_tags=set(), trials=None, invalid_times=None,
        intervals=None, units=None, processing=None, lab_meta_data=None,
        electrodes=None, electrode_groups=None, ic_electrodes=None,
        sweep_table=None, imaging_planes=None, ogen_sites=None, 
        devices=None, subject=None, scratch=None, icephys_electrodes=None
        """
        self.sessionIdx = sessionIdx
        self.device = self.create_device(name="Heka EPC10_USB Quadro")  #
        self._setStimulusElectrodes()
        self._setRecordingElectrodes()
        self._setCellData()

    def _HekaBundle(self, datFile):
        try:
            self._bundle = HEKA.Bundle(datFile)
        except:
            raise HekaFileError("Could not read .dat file format")
        self.pul = self._bundle.pul  ## get pulse tree if HEKA file is read

    def _setStimulusElectrodes(self):
        """Stimulating electrode id is usually in the Pulse tree's @Series' label, such as
        "Sti 1S-cc", or "Cell 5 fp". Or in the stimuli tree's @Stimulation's label
        We assume there's only one channel is in use for stimulating. And we use the metadata
        from the first channel.
        Ref: https://pynwb.readthedocs.io/en/stable/pynwb.icephys.html
        """
        idx = list(self.sessionIdx.copy())
        idx.append(0)  ## assume only 1 channel is stimulating
        stimTime, stimData, stimInfo = self._bundle.stim(idx)

        self.SeriesName = (
            stimInfo[0]["EntryName"] + "_G" + str(idx[0]) + "_S" + str(idx[1])
        )
        self.ccsElectrodes = []  ## current clamp stimulating electrode
        self.ccsElectrodes.append(
            self.create_icephys_electrode(
                name="ccsElectrode0",  # trace.Label,
                description=self.SeriesName,  # trace.Label,
                device=self.device,
            )
        )
        datSerie = self.pul[idx[0]][idx[1]]
        #        stimData = np.zeros((len(stimData),len(datSerie.children)))
        for j, sweep in enumerate(datSerie.children):
            #            ## TODO: need to work with Ampl file in the .dat file to get the "gain"
            idx = list(self.sessionIdx.copy())
            idx.append(j)  ## append this sweep
            if j<10:
                sName = '00'+str(j)
            elif j<100:
                sName = '0'+str(j)
            else:
                sName = str(j)
            stimTime, stimData, stimInfo = self._bundle.stim(idx)
            ##??? which version support "ValueError: Specifying rate and timestamps is not supported."
            self.add_stimulus(
                CurrentClampStimulusSeries(
                    name=self.SeriesName + "_sw" + sName,
                    data=np.array(stimData),
                    #                        starting_time = stimTime[-1]*j,
                    rate=1.0 / stimInfo[0]["sampleInteval"],
                    electrode=self.ccsElectrodes[0],
                    gain=0.02,  # starting_time = stimInfo[stimEpoch]['start'] * stimInfo[0][1]['sampleInteval']
                    sweep_number=j,
                    conversion=1e-12,
                    unit="amperes",
                ),use_sweep_table=True
            )

    def _setRecordingElectrodes(self):
        ## get all electrodes
        idx = list(self.sessionIdx.copy())
        idx = list([idx[0], idx[1]])
        idx.append(0)
        self.ccrElectrodes = []  ## current clamp recording electectrodes
        for trace in self.pul[idx[0]][idx[1]][idx[2]].children:
            #            print(trace)
            self.ccrElectrodes.append(
                self.create_icephys_electrode(
                    name=trace.Label, description="", device=self.device
                )
            )

    def _setCellData(self):
        idx = list(self.sessionIdx.copy())
        datSerie = self.pul[idx[0]][idx[1]]

        #        self.get_ic_electrode
        for s, sweep in enumerate(datSerie.children):
            idx = list(self.sessionIdx.copy())  ## make sure base idx is the same
            baseIdx = list([idx[0], idx[1]])
            baseIdx.append(s)
            baseIdx.append(0)
            if s<10:
                sName = '00'+str(s)
            elif s <100:
                sName = '0'+str(s)
            else:
                sName = str(s)
            for t, trace in enumerate(sweep.children):
                baseIdx[3] = t
                self.add_acquisition(
                    CurrentClampSeries(
                        name=self.SeriesName + "_sw" + sName +'_'+ trace.Label,
                        data=self._bundle.data[baseIdx],
                        rate=1.0 / trace.XInterval,
                        electrode=self.get_ic_electrode(trace.Label),
                        conversion=1.0,
                        unit="volts",
                        gain=0.02,
                        sweep_number=s,
                    ), use_sweep_table=True
                )

    def saveNBF(self, filename):
        ''' Save as Nb format
        '''
        with NWBHDF5IO(filename, "w") as io:
            io.write(self)

    def getNumberOfSweeps(self):
        if ~ hasattr(self, 'nSweep'):
            self.nSweep = len(self.sweep_table)//2
        return self.nSweep

    def getSweep(self, i):
        ''' Get series with index i
        return stimulus and response as sweep structure
        '''
        self.getNumberOfSweeps()
        if i >=0 and i<self.nSweep:
            stimulus, response = self.sweep_table.get_series(i)
            return stimulus, response
        else:
            print('sweep index is not valid. Number of sweep should be less than ', self.nSweep)
            
    
    def _sweepPlot(self, sweep,ax=None):
        ''' helper visualizer
        '''
        if ax is None:
            _, ax = plt.subplots()
        dat = sweep.data[:]
        yy = dat *sweep.conversion
        xx = np.arange(len(dat))/sweep.rate    
        ax.plot(xx, yy)   
    
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
            stimulus, response = self.getSweep(s)
            self._sweepPlot(response, ax=axes[0])
            self._sweepPlot(stimulus,  ax=axes[1])
            axes[0].set_xlabel('')
            axes[0].set_ylabel(response.unit)
            axes[1].set_ylabel(stimulus.unit)
            axes[1].set_xlabel('time (s)')
        plt.show()


class HekaFileError(Exception):
    """Generic Python-exception-derived object raised by any function."""
    pass
