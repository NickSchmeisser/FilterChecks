import os
import csv
import numpy as np
from glob import glob
import argparse

def main(goodrun_dir = '/data/ana/IceCube/2024/filtered/OfflinePass3.1/IC86_3024_GoodRunInfo.txt', out_dir = "/data/user/nschmeisser/TFT/Pass3_automation/processed2/"):
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
            os.mkdir(out_dir + directory_name + "/fpf")
            print(f"Directory '{directory_name}' created successfully.")
        except FileExistsError:
            print(f"Directory '{directory_name}' already exists.")
        except PermissionError:
            print(f"Permission denied: Unable to create '{directory_name}'.")
        except Exception as e:
            print(f"An error occurred: {e}")
    
    for i in range(len(run_dirs)):
        if(run_years[i] in years):
            directory_name = out_dir + str(run_years[i])+"/fpf/" + str(run_numbers[i])
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
            directory_name = out_dir + str(run_years[i])+"/fpf/" + str(run_numbers[i])
            print(directory_name)
            print(run_dirs[i])
            print(run_numbers[i])
            
    
    d=open('/home/nschmeisser/TFT/automation2/fpf_automation_final.dag', 'w')
    
    for i in range(len(run_dirs)):
        directory_name = out_dir + str(run_years[i])+"/fpf/" + str(run_numbers[i])
        fname="/home/nschmeisser/TFT/automation2/fpfscan_"+str(run_numbers[i])+"_jobscript.sub"
        f = open(fname, 'w')
        d.write("JOB " +"FPFJob_"+str(run_numbers[i])+ ' ' + fname + '\n')

        f.write('executable=/data/user/nschmeisser/TFT/Pass3_automation/faint_particle_filter_validation_v2.py\n')

        f.write("output = /home/nschmeisser/TFT/automation2/logs/fpfscan_"+str(run_numbers[i])+ ".out\n")
        f.write("error = /home/nschmeisser/TFT/automation2/logs/fpfscan_"+str(run_numbers[i])+ ".err\n")
        f.write("log = /home/nschmeisser/TFT/automation2/logs/fpfscan_"+str(run_numbers[i])+ ".log\n")

        f.write('request_memory = 1.5 GB\n')
        f.write(f'Arguments = -i {run_dirs[i]} -o {directory_name} -rn {run_numbers[i]} -lt {livetimes[i]} \n')
        f.write('queue\n')
        f.close()
    d.close()

if __name__ == "__main__":
    # Argument parsing for input and output directories
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-grl", "--goodrunlist", dest="GOODRUNS", required=True, help="Input GoodRunInfo.txt file")
    parser.add_argument("-o", "--outfolder", dest="OUTPUT", required=True, help="Output folder for all runs in GoodRunInfo")
    args = parser.parse_args()
    main(args.GOODRUNS,  args.OUTPUT)