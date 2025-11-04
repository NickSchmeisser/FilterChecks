#!/usr/bin/env python3
"""Validation script for the monopole filter."""
#ruff: noqa: F401 PLW0603 PTH118 PTH207
import argparse
from icecube import dataio, dataclasses, recclasses, icetray
from icecube.icetray import I3Tray
from glob import glob
import numpy as np
import h5py
import os

#Global counts
monopolefilter_count = 0
q_count = 0

# Labels for extracted variables
labels = [
    "MonopoleFilter_IC_LineFitI_zenith", "MonopoleFilter_IC_LineFitI_azimuth",
    "MonopoleFilter_IC_Pulses_First_Fit_zenith", "MonopoleFilter_IC_Pulses_First_Fit_azimuth",
    "MonopoleFilter_IC_Pulses_Second_Fit_zenith", "MonopoleFilter_IC_Pulses_Second_Fit_azimuth",
    "n_hit_doms", "n_hit_strings", "n_hit_doms_one_pulse", "n_pulses",
    "speed", "timelength_fwhm", "timelength_last_first", "timelength_maxgap",
    "zpattern", "avg_dom_dist_q_tot_dom", "empty_hits_track_length",
    "track_hits_separation_length", "track_hits_distribution_smoothness",
    "MonopoleFilter_IC_Pulses_First_Fit_speed", "MonopoleFilter_IC_Pulses_Second_Fit_speed",
    "MonopoleFilter_IC_LineFitIParams_vel", "MonopoleFilter_IC_LineFitIParams_vel_x",
    "MonopoleFilter_IC_LineFitIParams_vel_y", "MonopoleFilter_IC_LineFitIParams_vel_z",
    "MonopoleFilter_IC_LineFitIParams_vel_nhits",
]
labels = [np.bytes_(label) for label in labels]

p_frame_keys = np.array([
    ["MonopoleFilter_IC_LineFitI", "dir"], ["MonopoleFilter_IC_Pulses_First_Fit", "dir"],
    ["MonopoleFilter_IC_Pulses_Second_Fit", "dir"], ["MonopoleFilter_IC_HitMultiplicityValues", "n_hit_doms"],
    ["MonopoleFilter_IC_HitMultiplicityValues", "n_hit_strings"],
    ["MonopoleFilter_IC_HitMultiplicityValues", "n_hit_doms_one_pulse"],
    ["MonopoleFilter_IC_HitMultiplicityValues", "n_pulses"], ["MonopoleFilter_IC_LineFitI", "speed"],
    ["MonopoleFilter_IC_TimeCharacteristicsValues", "timelength_fwhm"],
    ["MonopoleFilter_IC_TimeCharacteristicsValues", "timelength_last_first"],
    ["MonopoleFilter_IC_TimeCharacteristicsValues", "timelength_maxgap"],
    ["MonopoleFilter_IC_TimeCharacteristicsValues", "zpattern"],
    ["MonopoleFilter_IC_TrackCharacteristicsValues", "avg_dom_dist_q_tot_dom"],
    ["MonopoleFilter_IC_TrackCharacteristicsValues", "empty_hits_track_length"],
    ["MonopoleFilter_IC_TrackCharacteristicsValues", "track_hits_separation_length"],
    ["MonopoleFilter_IC_TrackCharacteristicsValues", "track_hits_distribution_smoothness"],
    ["MonopoleFilter_IC_Pulses_First_Fit", "speed"], ["MonopoleFilter_IC_Pulses_Second_Fit", "speed"],
    ["MonopoleFilter_IC_LineFitIParams", "LFVel"], ["MonopoleFilter_IC_LineFitIParams", "LFVelX"],
    ["MonopoleFilter_IC_LineFitIParams", "LFVelY"], ["MonopoleFilter_IC_LineFitIParams", "LFVelZ"],
    ["MonopoleFilter_IC_LineFitIParams", "NHits"],
])

class MonopoleFilterVariables(icetray.I3Module):
    """Extracts variables from p-frames."""

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)

    def Configure(self):
        """Initialize the global 'var' array."""
        global var
        var = [np.array([]) for _ in range(len(labels))]  # Initialize empty arrays for all variables

    def DAQ(self, frame):
        """Process the q-frames and count those filtered by the MonopoleFilter."""
        global q_count, monopolefilter_count
        q_count += 1

        if frame.Has("QOfflineFilterMask"):
            filters = frame.Get("QOfflineFilterMask")
            if filters["MonopoleFilter_24"].condition_passed:
                monopolefilter_count += 1

        self.PushFrame(frame)

    def Physics(self, frame):
        """Process p-frames."""
        self.Check(frame)
        self.PushFrame(frame)

    def Check(self, frame):
        """Check offline filters and call GetVariables for MonopoleFilter events."""
        if frame.Has("OfflineFilterMask") and frame.Has("SplitInIcePulses"):
            filters = frame.Get("OfflineFilterMask")
            if filters["MonopoleFilter_24"].condition_passed:
                self.GetVariables(frame)

    def GetVariables(self, frame):
        """Extract variables from the frame."""
        var_count = -1

        for k in np.arange(0, 8, 2):
            var_count += 1
            if frame.Has(p_frame_keys[var_count][0]) and p_frame_keys[var_count][1] == "dir" and p_frame_keys[var_count][0] in \
        ["MonopoleFilter_IC_LineFitI", "MonopoleFilter_IC_Pulses_First_Fit", "MonopoleFilter_IC_Pulses_Second_Fit"]:
                    key = frame.Get(p_frame_keys[var_count][0])
                    v = getattr(key, p_frame_keys[var_count][1])
                    var[k] = np.append(var[k], v.zenith)
                    var[k + 1] = np.append(var[k + 1], v.azimuth)

        for k in np.arange(6, 26, 1):
            if frame.Has(p_frame_keys[k-3][0]):
                key = frame.Get(p_frame_keys[k-3][0])
                v = getattr(key, p_frame_keys[k-3][1])
                var[k] = np.append(var[k], v)

def main(in_dir, out_dir, run_num):
    """Execute main function."""
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--infolder", dest="INPUT", required=True, help="Input folder containing .i3.zst files")
    parser.add_argument("-o", "--outfolder", dest="OUTPUT", required=True, help="Output folder for HDF5 files")
    parser.add_argument("-rn", "--runnumbers",nargs='*', dest="RUNNUMBERS", required=False, help="Runnumbers on the respective date")
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
    
    for count, input_file in enumerate(infiles, start=1):
        tray = I3Tray()
        tray.AddModule("I3Reader", "Reader", Filename=input_file)
        tray.AddModule(MonopoleFilterVariables, "MonopoleFilterVariables")
        tray.AddModule("TrashCan", "Trash")
        tray.Execute()
        tray.Finish()
        del tray

        # Save extracted variables in an HDF5 file
        with h5py.File(f"{out_dir}/variable_{count}.h5", "w") as hf:
            for i in range(10):
                dataset = hf.create_dataset(f"var{i}", data=var[i])
                dataset.attrs["label"] = labels[i]

    # Save summary counts
    counts = np.array([q_count, monopolefilter_count])
    np.save(f"{out_dir}/summary.npy", counts)

if __name__ == "__main__":
    main()
