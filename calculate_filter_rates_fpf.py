import os
import csv
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
plt.style.use('tableau-colorblind10')
from scipy.stats import kstest, poisson

goodrun_dir ='/data/ana/IceCube/2024/filtered/OfflinePass3.1/IC86_3024_GoodRunInfo.txt'
goodrun = np.genfromtxt(goodrun_dir,skip_header=1, dtype='str')

#reference_path = "/data/user/nschmeisser/TFT/25/2025/fpf/0613"
#ref_file = np.load(reference_path+"/summary.npy")
ref_rate = 0#ref_file[1]/(28811.29+28669.57)

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

out_dir = "/data/user/nschmeisser/TFT/Pass3_automation/processed/"

rates = [[], [], []]#years, rates, p_values
for i in range(len(run_numbers)):
    summary_dir = out_dir+str(run_years[i])+"/fpf/"+str(run_numbers[i])
    if(os.path.isdir(summary_dir)):
        summary = np.load(summary_dir+"/summary.npy")
        if(i>100 and ref_rate==0):
            ref_rate = summary[1]/livetimes[i]
        rates[0].append(run_years[i])
        rates[1].append(summary[1]/livetimes[i])
        if(rates[1][-1]>10 or rates[1][-1]<0.8):
           print(run_numbers[i])

j=0
for i in range(len(run_numbers)):
    summary_dir = out_dir+str(run_years[i])+"/fpf/"+str(run_numbers[i])
    if(os.path.isdir(summary_dir)):
        expected_counts = ref_rate * livetimes[i]
        p_value = poisson.pmf(rates[1][j]*livetimes[i], expected_counts)
        rates[2].append(p_value)
        j+=1

fig = plt.figure(figsize=(7, 5))
ax0=plt.subplot2grid((1,1), (0, 0))
ax0.plot([int(i) for i in rates[0]], rates[1], linestyle='None',markersize=7, marker='.')
ax0.set_xlabel("Years", fontsize=15)
ax0.set_ylabel("Rate in Hz", fontsize=15)
ax0.grid(True, linewidth=0.5)
#ax0.set_yscale("log")
fig.tight_layout()
fig.savefig("./plots_fpf/fpf_rates.png", bbox_inches='tight')
plt.close(fig)

fig = plt.figure(figsize=(7, 5))
ax0=plt.subplot2grid((1,1), (0, 0))
ax0.plot([int(i) for i in rates[0]], rates[2], linestyle='None',markersize=7, marker='.')
ax0.set_xlabel("Years", fontsize=15)
ax0.set_ylabel("p_values", fontsize=15)
ax0.set_title("Filter rates", fontsize=15)
ax0.grid(True, linewidth=0.5)
ax0.set_yscale("log")
fig.tight_layout()
fig.savefig("./plots_fpf/fpf_rates_pvalues.png", bbox_inches='tight')
plt.close(fig)

fig = plt.figure(figsize=(7, 5))
ax0=plt.subplot2grid((1,1), (0, 0))
ax0.plot([int(i) for i in rates[0]], rates[1], linestyle='None',markersize=7, marker='.')
ax0.set_xlabel("Years", fontsize=15)
ax0.set_ylabel("Rate in Hz", fontsize=15)
ax0.set_ylim(0.9,1.25)
ax0.grid(True, linewidth=0.5)
#ax0.set_yscale("log")
fig.tight_layout()
fig.savefig("./plots_fpf/fpf_rates_zoom.png", bbox_inches='tight')
plt.close(fig)

fig = plt.figure(figsize=(7, 5))
ax0=plt.subplot2grid((1,1), (0, 0))
hist_main, binedges = np.histogram(rates[1], bins=100)
yerr_main = np.sqrt(hist_main)
X = np.array([binedges[:-1], binedges[1:]]).T.flatten()
bincenters = (binedges[:-1] + binedges[1:]) / 2
ax0.plot(X, np.repeat(hist_main, 2), color=f"C0", linewidth=2)
ax0.errorbar(bincenters, hist_main, yerr=yerr_main, fmt=".", color=f"C0")
ax0.set_xlabel("Rates in Hz", fontsize=15)
ax0.set_ylabel("Counts", fontsize=15)
ax0.grid(True, linewidth=0.5)
#ax0.set_yscale("log")
fig.tight_layout()
fig.savefig("./plots_fpf/fpf_rates_hist.png", bbox_inches='tight')
plt.close(fig)

fig = plt.figure(figsize=(7, 5))
ax0=plt.subplot2grid((1,1), (0, 0))
hist_main, binedges = np.histogram(rates[1], bins=1000)
yerr_main = np.sqrt(hist_main)
X = np.array([binedges[:-1], binedges[1:]]).T.flatten()
bincenters = (binedges[:-1] + binedges[1:]) / 2
ax0.plot(X, np.repeat(hist_main, 2), color=f"C0", linewidth=2)
ax0.errorbar(bincenters, hist_main, yerr=yerr_main, fmt=".", color=f"C0")
ax0.set_xlabel("Rates in Hz", fontsize=15)
ax0.set_ylabel("Counts", fontsize=15)
ax0.grid(True, linewidth=0.5)
ax0.set_xlim(0,2)
#ax0.set_yscale("log")
fig.tight_layout()
fig.savefig("./plots_fpf/fpf_rates_hist_zoom.png", bbox_inches='tight')
plt.close(fig)


