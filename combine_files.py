#!/usr/bin/env python3
"""Combine subfiles of faint and monpole filter validation scripts."""

import os
import h5py
import numpy as np
import argparse
from pathlib import Path

def combine_h5_files_in_subfolder(subfolder_path, output_file_path, num_datasets):
    """Combine all h5 files from subfolders and generate one combined file per subfolder."""
    combined_data = []
    labels = [None] * num_datasets

    for i in range(num_datasets):
        temp_list = []
        dataset_label = None

        for file_name in os.listdir(subfolder_path):
            if file_name.endswith(".h5"):
                file_path = Path(subfolder_path) / file_name
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

    with h5py.File(output_file_path, "w") as hf:
        for i, data in enumerate(combined_data):
            dataset_name = f"var{i}"
            hf.create_dataset(dataset_name, data=data)
            hf[dataset_name].attrs["label"] = labels[i].encode("utf-8")

def process_all_subfolders(root_folder, num_datasets):
    """Loop through all subfolders and call combine function."""
    for subfolder in os.listdir(root_folder):
        subfolder_path = Path(root_folder) / subfolder
        if Path(subfolder_path).is_dir():
            output_file_path = Path(subfolder_path) / "combined_data.h5"
            try:
                combine_h5_files_in_subfolder(subfolder_path, output_file_path, num_datasets)
            except:
                print("Error in directory: " + str(subfolder_path))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine HDF5 files in subfolders.")
    parser.add_argument("-f", "--root_folder", type=str, required=True, help="Root directory containing subfolders with .h5 files")
    parser.add_argument("-n","--num_datasets", type=int, required=True, help="Number of keys in h5 file (26 for MonopoleFilter and 10 for FaintParticleFilter")

    args = parser.parse_args()
    process_all_subfolders(args.root_folder, args.num_datasets)
