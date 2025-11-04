#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/icetray-start
#METAPROJECT /data/user/nschmeisser/Software/IceTrayExoticGenerator/build
import os
import csv
import numpy as np
from glob import glob
import monopole_filter_validation
import combine_files
import argparse

out_dir = "/data/user/nschmeisser/TFT/Pass3_automation/processed/"

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-rundir", "--run_dir", dest="RUNDIR", required=False, help="Location of run")
parser.add_argument("-runnr", "--run_nr", dest="RUNNR", required=False, help="Number of run")
parser.add_argument("-y", "--runyear", dest="RUNYEAR", required=False, help="Process only specific year", type=int)
args = parser.parse_args()


run_numbers = [args.RUNNR]
run_dirs = [args.RUNDIR]

years = [args.RUNYEAR]
run_years = [args.RUNYEAR]

print(years)

for i in years:
    directory_name = str(i)
    try:
        os.mkdir(out_dir + directory_name)
        print(f"Directory '{directory_name}' created successfully.")
    except FileExistsError:
        print(f"Directory '{directory_name}' already exists.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{directory_name}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

for i in years:
    directory_name = str(i)
    try:
        os.mkdir(out_dir + directory_name + "/mono")
        print(f"Directory '{directory_name}' created successfully.")
    except FileExistsError:
        print(f"Directory '{directory_name}' already exists.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{directory_name}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

for i in range(len(run_dirs)):
    if(run_years[i] in years):
        directory_name = out_dir + str(run_years[i])+"/mono/" + str(run_numbers[i])
        try:
            os.mkdir(directory_name)
            print(f"Directory '{directory_name}' created successfully.")
        except FileExistsError:
            print(f"Directory '{directory_name}' already exists.")
        except PermissionError:
            print(f"Permission denied: Unable to create '{directory_name}'.")
        except Exception as e:
            print(f"An error occurred: {e}")

for i in range(len(run_dirs)):
    if(run_years[i] in years):
        directory_name = out_dir + str(run_years[i])+"/mono/" + str(run_numbers[i])
        print(directory_name)
        print(run_dirs[i])
        print(run_numbers[i])
        #!python "/data/user/nschmeisser/TFT/Pass3_automation/0_faint_particle_filter_validation.py" -i run_dirs[i] -o directory_name -rn run_numbers[i]
        monopole_filter_validation.main(in_dir=run_dirs[i], out_dir = directory_name, run_num=run_numbers[i])
        print("Run: "+str(run_numbers[i])+" done!")
'''
for i in years:
    directory_name = out_dir + str(i)+"/mono/"
    combine_files.process_all_subfolders(directory_name, 10)
'''

