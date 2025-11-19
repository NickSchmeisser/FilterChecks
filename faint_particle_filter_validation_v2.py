#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/icetray-start
#METAPROJECT /data/user/nschmeisser/Software/IceTrayExoticGenerator/build
"""Validation script for the faint particle filter."""
#ruff: noqa: F401 PLW0603 PTH118 PTH207
import os
import argparse
from icecube import dataio, dataclasses, icetray
from icecube.icetray import I3Tray
from glob import glob
import numpy as np
import h5py
from pathlib import Path
import unittest
from scipy.stats import kstest, poisson

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

class Combine_variables:
    num_datasets = 10
    output_file_path = "/data/user/nschmeisser/TFT/Pass3_automation/processed/2018/fpf/130755"

    def combine_h5_files_in_subfolder(self, output_file_path, num_datasets = 10):
        """Combine all h5 files from subruns in output_file_path and generate one combined file per subfolder."""
        combined_data = []
        labels = [None] * num_datasets
    
        for i in range(num_datasets):
            temp_list = []
            dataset_label = None
    
            for file_name in os.listdir(output_file_path):
                if file_name.endswith(".h5"):
                    file_path = Path(output_file_path) / file_name
                    with h5py.File(file_path, "r") as hf:
                        data = hf[f"var{i}"][:]
                        temp_list.append(data)
    
                        if "label" in hf[f"var{i}"].attrs and dataset_label is None:
                            attr = hf[f"var{i}"].attrs["label"]
                            if isinstance(attr, bytes):
                                dataset_label = attr.decode("utf-8")
                            else:
                                dataset_label = attr
            if temp_list:
                combined_data.append(np.concatenate(temp_list))
    
                if dataset_label:
                    labels[i] = dataset_label
                else:
                    labels[i] = f"No label found for dataset {i}"

        output_file = Path(output_file_path) / "combined_data.h5"
        with h5py.File(output_file, "w") as hf:
            for i, data in enumerate(combined_data):
                dataset_name = f"var{i}"
                hf.create_dataset(dataset_name, data=data)
                hf[dataset_name].attrs["label"] = labels[i].encode("utf-8")

    def delete_subrun_files(self, output_file_path):
        for file_name in os.listdir(output_file_path):
            if(file_name.endswith(".h5") and not file_name.endswith("combined_data.h5")):
                output_file = Path(output_file_path) / str(file_name)
                os.remove(output_file)

class FPFTest(): #unittest.TestCase):
    num_datasets = 10
    output_file_path = "/data/user/nschmeisser/TFT/Pass3_automation/processed/2018/fpf/130755"
    ref_path = '/data/user/nschmeisser/TFT/Pass3_automation/reference2018/2018/fpf/130764'
    ref_livetime = 28767.12
    pvalue_threshold = 1e-4
    livetime = 0
                     
    def fpf_rate_check(self, ref_livetime = ref_livetime):
        ref_rate_file = np.load(self.ref_path + "/summary.npy")
        ref_rate = ref_rate_file[1] / ref_livetime #reference rate/reference livetime
        summary = np.load(self.output_file_path+"/summary.npy") #read in counts of run
        counts = summary[1]
        expected_counts = ref_rate * self.livetime
        pvalue = poisson.pmf(counts, expected_counts)
        if (pvalue < self.pvalue_threshold):
            print("p_value for rate is lower than expected.")
        np.save(str(self.output_file_path)+"/fpf_rate_p_value.npy", np.array([pvalue]))
        
    def fpf_variable_check(self, num_datasets = num_datasets):
        ref_file = h5py.File(os.path.join(self.ref_path, "combined_data.h5"), 'r')
        num_vars = num_datasets
        labels = [None] * num_vars
        p_values = [None] * num_datasets
        for i in range(num_vars):
            ref_value = ref_file[f'var{i}'][:]
            file_path = os.path.join(self.output_file_path, "combined_data.h5")
            with h5py.File(file_path, 'r') as hf:
                data = hf[f'var{i}'][:]
                if(len(data>=1)):
                    pvalue = kstest(ref_value, data).pvalue
                    p_values[i] = pvalue
                    if (pvalue < self.pvalue_threshold):
                        print(f"p_value for Variable {i} is lower than expected.")
                    #message = f"The test for variable {i} failed"
                    #self.assertGreater(pvalue, self.pvalue_threshold, message)
        np.save(str(self.output_file_path)+"/fpf_p_values.npy", np.array(p_values))



def main(in_dir, out_dir, run_num, livetime):
    """Execute main function."""
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
        with h5py.File(str(out_dir)+f"/variable_{count}.h5", "w") as hf:
            for i in range(10):
                dataset = hf.create_dataset(f"var{i}", data=var[i])
                dataset.attrs["label"] = labels[i]
    print("Read i3-Files!")
    
    # Save summary of the counts
    counts = np.array([q_count, fpf_count])
    print(str(out_dir)+"/summary.npy")
    np.save(str(out_dir)+"/summary.npy", counts)
    cv = Combine_variables()
    cv.combine_h5_files_in_subfolder(output_file_path = out_dir)
    cv.delete_subrun_files(output_file_path = out_dir)
    print("Combined subruns!")
    
    # Performing statistical tests
    fpf_tests = FPFTest()
    fpf_tests.output_file_path = str(out_dir)
    fpf_tests.livetime = float(livetime[0])
    fpf_tests.fpf_rate_check()
    print("Rate checks for Faint Particle Filter done!")
    fpf_tests.fpf_variable_check()
    print("Variable checks for Faint Particle Filter done!")
    

if __name__ == "__main__":
    # Argument parsing for input and output directories
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--infolder", dest="INPUT", required=True, help="Input folder containing .i3.zst files")
    parser.add_argument("-o", "--outfolder", dest="OUTPUT", required=True, help="Output folder for HDF5 files")
    parser.add_argument("-rn", "--runnumber",nargs='*', dest="RUNNUMBER", required=True, help="Runnumber on the respective date")
    parser.add_argument("-lt", "--livetime",nargs='*', dest="LIVETIME", required=True, help="Livetime of the respective run")
    args = parser.parse_args()
    main(args.INPUT, args.OUTPUT, args.RUNNUMBER, args.LIVETIME)
