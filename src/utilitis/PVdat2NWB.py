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

from HekaIO import HEKA_Reader_MAIN as HEKA

class XLJdat2NWB(NWBFile):
    def __init__(self, datFile, sessionIdx, **kw):
        """
        parameters:sessionIdx
        dataFile: dat filename
        sessionIdx: [group, series]  ## one file for one series
        """
        self._HekaBundle(datFile)
        super(XLJdat2NWB, self).__init__(
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
        self.device = self.create_device(name="Heka EPC10_USB Quadro")  ##?
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
        We assume there's only one channel is in use for stimululating. And we use the metadata
        from the first channel.
        Ref: https://pynwb.readthedocs.io/en/stable/pynwb.icephys.html
        """
        idx = list(self.sessionIdx.copy())
        idx.append(0)  ## assume only 1 channel is stimulating
        stimTime, stimData, stimInfo = self._bundle.stim(idx)

        #        for elec in self.pul[idx[0]][idx[1]][idx[2]].children:
        #            self.StimElecs.append(self.create_icephys_electrode(name=stimInfo[0]['EntryName'], #trace.Label,
        #                                    description='Great data', #trace.Label,
        #                                   device=self.device))

        #                    stim_ = {'EntryName': stimulatingElectrode, 'start': cumSamples,
        #                     'end': cumSamples + segSamples, 'duration': seg.Duration,
        #                     'amplitude': seg_V,  'sampleInteval': sampleInteval,
        #             'Vholding': Holding, 'Stim2Dac':stim2dacType,'sealresistance':sealResistance}
        ## grab the epoc when the real stimulu starts.
        #        stimEpoch = 0
        ##        print(stimInfo)
        #        for stimIdx, stim in enumerate(stimInfo):
        #            if abs(stim['amplitude'] - stim['Vholding']) > 1:  ## first non-zero amplitude
        #                stimEpoch = stimIdx
        ##                print('stimEpoch: %g' % stimEpoch)
        #                break
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
            stimTime, stimData, stimInfo = self._bundle.stim(idx)
            ##??? which version support "ValueError: Specifying rate and timestamps is not supported."
            self.add_stimulus(
                CurrentClampStimulusSeries(
                    name=self.SeriesName + "_sw" + str(j),
                    data=np.array(stimData),
                    #                        starting_time = stimTime[-1]*j,
                    rate=1.0 / stimInfo[0]["sampleInteval"],
                    electrode=self.ccsElectrodes[0],
                    gain=0.02,  # starting_time = stimInfo[stimEpoch]['start'] * stimInfo[0][1]['sampleInteval']
                    sweep_number=j,
                    conversion=1e-12,
                    unit="amperes",
                )
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
            for t, trace in enumerate(sweep.children):
                baseIdx[3] = t
                self.add_acquisition(
                    CurrentClampSeries(
                        name=self.SeriesName + "sw" + str(s) + trace.Label,
                        data=self._bundle.data[baseIdx],
                        rate=1.0 / trace.XInterval,
                        electrode=self.get_ic_electrode(trace.Label),
                        conversion=1.0,
                        unit="volts",
                        gain=0.02,
                        sweep_number=s,
                    )
                )

    def saveNBF(self, filename):
        with NWBHDF5IO(filename, "w") as io:
            io.write(self)


#        print('writing done!')


class HekaFileError(Exception):
    """Generic Python-exception-derived object raised by any function."""

    pass
