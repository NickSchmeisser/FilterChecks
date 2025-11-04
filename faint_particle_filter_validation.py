#!/usr/bin/env python3
"""Validation script for the faint particle filter."""
#ruff: noqa: F401 PLW0603 PTH118 PTH207
import os
import argparse
from icecube import dataio, dataclasses, icetray
from icecube.icetray import I3Tray
from glob import glob
import numpy as np
import h5py

#Global counts
q_count = 0
fpf_count = 0

#Variable labels
labels = np.array([ "#Launches", "#Doubles",
                        "#Azimuth", "#Zenith", "SLC_fraction",
                        "Trigger_length_ns", "LineFit_zenith_rad",
                        "LineFit_azimuth_rad", "MPE_fit_zenith_rad",
                        "MPE_fit_azimuth_rad" ])
labels = [np.bytes_(label) for label in labels]

#Load lookup tables for distances and directions
distance_lookup = np.load("/data/user/nschmeisser/TFT/lookup_files/distance_lookup.npy")
direction_lookup_zen = np.load("/data/user/nschmeisser/TFT/lookup_files/direction_lookup_zen.npy")
direction_lookup_azi = np.load("/data/user/nschmeisser/TFT/lookup_files/direction_lookup_azi.npy")

class Launch:
    """Launch class."""

    def __init__(self, OMKey, Time):
        """Initialize the Launch object with an OMKey and Time."""
        self.OMKey = OMKey
        self.Time = Time

class FPFVariables(icetray.I3Module):
    """A module for calculating and storing FPF variables for a given event."""

    def __init__(self, context):
        """Initialize the FPFVariables module."""
        icetray.I3Module.__init__(self, context)

    def Configure(self):
        """Initialize the global 'var' array."""
        global var
        var = [np.array([]) for _ in range(len(labels))]

    def DAQ(self, frame):
        """Process the q-frames and count those filtered by the FPF."""
        global q_count,fpf_count
        q_count += 1
        if frame.Has("QOfflineFilterMask"):
            filters = frame.Get("QOfflineFilterMask")
            if filters["FaintParticleFilter_24"].condition_passed:
                fpf_count += 1
        self.PushFrame(frame)

    def Physics(self, frame):
        """Process p-frames."""
        self.Check(frame)
        self.PushFrame(frame)

    def Check(self, frame):
        """Check the conditions for offline filters and call GetVariables for FPF events."""
        if frame.Has("OfflineFilterMask") and (frame.Has("FPF_MaskedInIcePulses") or frame.Has("SplitFaintPulses")):
            filters = frame.Get("OfflineFilterMask")
            if filters["FaintParticleFilter_24"].condition_passed:
                self.GetVariables(frame)

    def GetVariables(self, frame):
        """Extract and calculate various variables from the event frame."""
        global cut_doublet_max, cut_doublet_min, binwidth

        dc_strings = {79, 80, 81, 82, 83, 84, 85, 86}
        ic_dc_strings = {25, 26, 27, 34, 35, 36, 37, 44, 45, 46, 47, 54}
        dc_bounds = (11, 60)
        ic_dc_bounds = (39, 60)

        binwidth = 20
        cut_doublet_max = 310000
        cut_doublet_min = 10000

        triggers = frame["I3TriggerHierarchy"]
        triggers_unique = {k.key.config_id for k in triggers if k.key.config_id is not None}
        trigger_id = 33001 if 33001 in triggers_unique else 1011 if 1011 in triggers_unique else None

        for trigger in triggers:
            if trigger.key.config_id == trigger_id:
                pulse_time_delta = 150
                event_start = trigger.time - pulse_time_delta
                event_stop = event_start + trigger.length + pulse_time_delta

                pulses_window = []
                slcfraction_window = []

                if frame.Has("FPF_MaskedInIcePulses"):
                    pulses_frame = "FPF_MaskedInIcePulses"
                    direction_mpe = "FPF_MPEFit_hive"
                    direction_line = "FPF_LineFit_hive"
                else:
                    pulses_frame = "SplitFaintPulses"
                    direction_mpe = "FPF_MPEFit"
                    direction_line = "FPF_LineFit"

                pulse_map = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulses_frame)

                for omkey, pulse_series in pulse_map.items():
                    is_dc = omkey[0] in dc_strings and dc_bounds[0] <= omkey[1] <= dc_bounds[1]
                    is_ic_dc = omkey[0] in ic_dc_strings and ic_dc_bounds[0] <= omkey[1] <= ic_dc_bounds[1]

                    if not (is_dc or is_ic_dc):
                        continue

                    for pulse in pulse_series:
                        if event_start <= pulse.time <= event_stop:
                            pulses_window.append(Launch(omkey, pulse.time))
                            slcfraction_window.append(1 if pulse.flags in {3, 5, 7} else 0)

                if pulses_window:
                    self.process_pulses(frame, pulses_window, np.array(slcfraction_window), trigger, direction_mpe, direction_line)

    def process_pulses(self, frame, pulses_window, slcfraction_window, trigger, direction_mpe, direction_line):
        """Process the collected pulses and calculate variables."""
        pulses_window.sort(key=lambda launch: launch.Time)

        slcfrac = len(slcfraction_window[slcfraction_window==0])/len(pulses_window)

        surv_doublets, doub_ind = self.CalcDoubles(pulses_window)
        zen, azi = self.CalcDirection(doub_ind, pulses_window)

        azihist = np.histogram(azi, bins=np.arange(0, 360 + binwidth, binwidth))[0]
        zenhist = np.histogram(zen, bins=np.arange(0, 180 + binwidth, binwidth))[0]

        var[0] = np.append(var[0], len(pulses_window))
        var[1] = np.append(var[1], len(surv_doublets))
        var[2] = np.append(var[2], np.max(azihist))
        var[3] = np.append(var[3], np.max(zenhist))
        var[4] = np.append(var[4], slcfrac)
        var[5] = np.append(var[5], trigger.length)
        var[6] = np.append(var[6], frame[direction_line].dir.zenith)
        var[7] = np.append(var[7], frame[direction_line].dir.azimuth)
        var[8] = np.append(var[8], frame[direction_mpe].dir.zenith)
        var[9] = np.append(var[9], frame[direction_mpe].dir.azimuth)

    def CalcDirection(self, indices_doub, selectedpulsemap):
        """Calculate zenith and azimuth directions for Doubles."""
        num_pairs = len(indices_doub) // 2
        zenith = np.zeros(num_pairs)
        azimuth = np.zeros(num_pairs)

        for i in range(num_pairs):
            hit1 = selectedpulsemap[int(indices_doub[2 * i])]
            hit2 = selectedpulsemap[int(indices_doub[2 * i + 1])]
            dir_lookup_zen = direction_lookup_zen[hit1.OMKey[0], hit1.OMKey[1], hit2.OMKey[0], hit2.OMKey[1]]
            dir_lookup_azi = direction_lookup_azi[hit1.OMKey[0], hit1.OMKey[1], hit2.OMKey[0], hit2.OMKey[1]]

            zenith[i] = np.degrees(dir_lookup_zen)
            azimuth[i] = np.degrees(dir_lookup_azi)

        return zenith, azimuth

    def CalcDoubles(self, selectedpulsemap):
        """Calculate Double velocities and return the Doubles that meet the criteria."""
        cut_doublet_min = 10000
        cut_doublet_max = 310000

        double_velocity = []
        index_list = []

        for indlist1, pulse1 in enumerate(selectedpulsemap):
            for indlist2, pulse2 in enumerate(selectedpulsemap[indlist1 + 1:], start=indlist1 + 1):
                if pulse1.OMKey != pulse2.OMKey and pulse1.Time != pulse2.Time:
                    distance12 = self.GetDistance(pulse1.OMKey, pulse2.OMKey)
                    timediff = abs(pulse1.Time - pulse2.Time)
                    vel_doublet = distance12 / timediff * 10**6

                    if cut_doublet_min < vel_doublet < cut_doublet_max:
                        double_velocity.append(vel_doublet)
                        index_list.extend([indlist1, indlist2])

        return np.array(double_velocity), np.array(index_list)

    def GetDistance(self, hit1, hit2):
        """Get the distance between two hits."""
        if hit2 < hit1:
            hit1, hit2 = hit2, hit1
        return distance_lookup[hit1[0], hit1[1], hit2[0], hit2[1]]

def main(in_dir, out_dir, run_num):
    """Execute main function."""
    '''
    # Argument parsing for input and output directories
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--infolder", dest="INPUT", required=True, help="Input folder containing .i3.zst files")
    parser.add_argument("-o", "--outfolder", dest="OUTPUT", required=True, help="Output folder for HDF5 files")
    parser.add_argument("-rn", "--runnumber",nargs='*', dest="RUNNUMBER", required=False, help="Runnumber on the respective date")
    args = parser.parse_args()
    '''


    # Get the list of files to process
    runnumber=run_num #args.RUNNUMBER
    print(runnumber)
    infiles=[]
    pattern = str(in_dir)+ "/*"+str(runnumber[0])+ "*.zst"
    files = glob(pattern)
    print(files)
    for j in files:
        infiles.append(j)

    # Process each file
    for count, input_file in enumerate(infiles, start=1):
        tray = I3Tray()
        tray.AddModule("I3Reader", "Reader", Filename=input_file)
        tray.AddModule(FPFVariables, "FPFVariables")
        tray.AddModule("TrashCan", "Trash")

        tray.Execute()
        tray.Finish()
        del tray

        # Write the data to HDF5 file
        with h5py.File(f"{out_dir}/variable_{count}.h5", "w") as hf:
            for i in range(10):
                dataset = hf.create_dataset(f"var{i}", data=var[i])
                dataset.attrs["label"] = labels[i]

    # Save summary of the counts
    counts = np.array([q_count, fpf_count])
    np.save(f"{out_dir}/summary.npy", counts)

if __name__ == "__main__":
    main()
