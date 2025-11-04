import os
import csv
import numpy as np
from glob import glob

goodrun_dir ='/data/ana/IceCube/2024/filtered/OfflinePass3.1/IC86_3024_GoodRunInfo.txt'
goodrun = np.genfromtxt(goodrun_dir,skip_header=1, dtype='str')

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

d=open('/home/nschmeisser/TFT/automation/fpf_automation_perrun.dag', 'w')

for i in range(len(run_numbers)):
        fname="/home/nschmeisser/TFT/automation/fpfscan_"+str(run_numbers[i])+"_jobscript.sub"
        f = open(fname, 'w')
        d.write("JOB " +"FPFJob_"+str(run_numbers[i])+ ' ' + fname + '\n')

        f.write('executable=/data/user/nschmeisser/TFT/Pass3_automation/pass3_automation_fpf_perrun.py\n')

        f.write("output = /home/nschmeisser/TFT/automation/logs/fpfscan_"+str(run_numbers[i])+ ".out\n")
        f.write("error = /home/nschmeisser/TFT/automation/logs/fpfscan_"+str(run_numbers[i])+ ".err\n")
        f.write("log = /home/nschmeisser/TFT/automation/logs/fpfscan_"+str(run_numbers[i])+ ".log\n")

        f.write('request_memory = 1.5 GB\n')
        f.write("Arguments=-rundir "+str(run_dirs[i]) + " -runnr "+str(run_numbers[i]) + " -y "+ str(run_years[i]) + "\n")
        f.write('queue\n')
        f.close()
d.close()