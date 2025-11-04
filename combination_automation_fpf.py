import os
import csv
import numpy as np
from glob import glob
import combine_files

goodrun_dir ='/data/ana/IceCube/2024/filtered/OfflinePass3.1/IC86_3024_GoodRunInfo.txt'
out_dir = "/data/user/nschmeisser/TFT/Pass3_automation/processed/"

goodrun = np.genfromtxt(goodrun_dir,skip_header=1, dtype='str')
num_datasets=10

run_numbers = []
livetimes = []
run_dirs = []

for i in goodrun:
    run_numbers.append(i[0])
    livetimes.append(float(i[3]))
    run_dirs.append(i[7])

years = []
run_years = []

for i in run_dirs:
    run_year = i[18:22]
    run_years.append(run_year)
    if ( run_year not in years):
        years.append(run_year)

for i in years:
    print(i)
    combination_dir = out_dir + str(i) + "/fpf/"
    combine_files.process_all_subfolders(combination_dir, num_datasets)
    print("Year done!")