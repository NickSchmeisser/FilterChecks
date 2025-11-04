import argparse
import os
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('tableau-colorblind10')
import warnings
from scipy.stats import kstest, poisson
"""
%run plotting_fpf_mf.py \
  -i "/data/user/tstuerwald/TFT/24/testrun/fpf_dev9" \
  -ir "/data/user/tstuerwald/TFT/24/testrun/reference/fpf" \
  -o"/data/user/tstuerwald/TFT/24/testrun/test"\
  -n 
"""
labels_fpf = np.array([
    "#Launches", "#Doubles", "#Azimuth", "#Zenith", 
    "SLC fraction", "Trigger length [ns]", "LineFit zenith [rad]",
    "LineFit azimuth [rad]", "MPE fit zenith[rad]", "MPE fit azimuth[rad]"
])

def load_combined_data(base_path):
    years_combined = []
    years_labels = []
    years = [ f.name for f in os.scandir(base_path) if (f.is_dir() and (f.name)[0]!='.') ]
    num_vars = 10

    for j in range(len(years)):
        combined = []
        labels = [None] * num_vars
        for i in range(num_vars):
            runs = [ f.path for f in os.scandir(str(base_path)+str(years[j])+"/fpf/") if f.is_dir() ]
            temp_list = []
            for run in runs:
                file_path = os.path.join(run, "combined_data.h5")
                with h5py.File(file_path, 'r') as hf:
                    data = hf[f'var{i}'][:]
                    temp_list.append(data)  
    
                    if labels[i] is None and 'label' in hf[f'var{i}'].attrs:
                        label = hf[f'var{i}'].attrs['label']
                        if isinstance(label, bytes):
                            label = label.decode('utf-8')
                        labels[i] = label
    
            combined_data = np.concatenate(temp_list)
            combined.append(combined_data)
        years_combined.append(combined)
        years_labels.append(labels)

    return years_combined, years_labels, years

def KS_tests(combined_years, combined_ref, labels, years, output_dir):
    for i, var_name in enumerate(labels):
        p_values = [[],[],[]] #Label/variable, Year, p_value
        var_ref = combined_ref[i]
        for j in range(len(combined_years)):
            var_main = combined_years[j][i]
            p_values[0].append(var_name)
            p_values[1].append(years[j])
            pvalue = kstest(var_ref,var_main).pvalue
            p_values[2].append(pvalue)
        #Plotting p_values for variable
        fig = plt.figure(figsize=(7, 3.5))
        ax0=plt.subplot2grid((1,1), (0, 0))
        print(p_values[1])
        print(p_values[0])
        ax0.plot([int(i) for i in p_values[1]],p_values[2],linestyle='None',markersize=7, marker='.')
        ax0.set_xlabel("Years", fontsize=15)
        ax0.set_ylabel("p_value", fontsize=15)
        ax0.set_title(var_name, fontsize=15)
        ax0.grid(True, linewidth=0.5)
        ax0.set_yscale("log")
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, f"{var_name}_pvalues.png"), bbox_inches='tight')
        plt.close(fig)
            

def plot_distributions(combined_years, combined_ref, labels, years, output_dir, normalize=False):
    os.makedirs(output_dir, exist_ok=True)

    # Special bin settings and x-axis options
    bin_settings = {
        "n_hit_doms": {"bins": 30, "xlim": (0, 550)},
        "speed": {"bins": 50, "xlim": (1e-2, 1.5), "xscale": "log"},
        "#Zenith": {"bins": 100},
        "#Azimuth": {"bins": 100},
    }

    for i, var_name in enumerate(labels):
        #var_main = combined_main[i]
        var_ref = combined_ref[i]
        fig = plt.figure(figsize=(7, 8))
        ax0 = plt.subplot2grid((3, 2), (0, 0), rowspan=2, colspan=2)
        ax1 = plt.subplot2grid((3, 2), (2, 0), rowspan=1, colspan=2)

        ax0.set_yscale("log")

        # Define binning
        settings = bin_settings.get(var_name, {})
        bins = settings.get("bins", 50)
        borders = [min([np.nanmin(j[i]) for j in combined_years]), max([np.nanmax(j[i]) for j in combined_years])]
        binning = np.linspace(borders[0], borders[1], bins + 1)

        for j in range(len(combined_years)):
            combined_main = combined_years[j]
            var_main = combined_main[i]
            main_label = years[j]
            '''
            # Define binning
            settings = bin_settings.get(var_name, {})
            bins = settings.get("bins", 50)
            borders = [np.nanmin(var_main), np.nanmax(var_main)]
            binning = np.linspace(borders[0], borders[1], bins + 1)
            '''
            # Histogram
            hist_main, binedges = np.histogram(var_main, bins=binning)
            hist_ref, _ = np.histogram(var_ref, bins=binning)
            bincenters = (binedges[:-1] + binedges[1:]) / 2
            binwidth = binedges[1] - binedges[0]
            X = np.array([binedges[:-1], binedges[1:]]).T.flatten()
    
            yerr_main = np.sqrt(hist_main)
            yerr_ref = np.sqrt(hist_ref)
    
            # Normalize histograms if requested
            if normalize:
                nevents_main = float(sum(hist_main))
                hist_main = hist_main / nevents_main / binwidth
                yerr_main = yerr_main / nevents_main / binwidth
    
                nevents_ref = float(sum(hist_ref))
                hist_ref = hist_ref / nevents_ref / binwidth
                yerr_ref = yerr_ref / nevents_ref / binwidth
    
                y_label = "Density"
            else:
                y_label = "Count"
    
            # Plot histograms
            
            ax0.plot(X, np.repeat(hist_main, 2), label=str(main_label)+f"_{var_name}", color=f"C{j}", linewidth=2)
            ax0.errorbar(bincenters, hist_main, yerr=yerr_main, fmt=".", color=f"C{j}")

            # Ratio plot
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                ratio = np.divide(hist_ref, hist_main, out=np.zeros_like(hist_ref, dtype=float), where=hist_main != 0)
            ax1.plot(X, np.repeat(ratio, 2), color=f"C{j}", linestyle='none', marker='x', label=f"reference/{main_label}")
            ax1.set_ylim(0.6, 1.4)
            ax1.set_xlabel(var_name, fontsize=10)
            ax1.set_ylabel(f"reference/{main_label}", fontsize=14)
            
        ax0.plot(X, np.repeat(hist_ref, 2), label="reference", color=f"C{j+1}", linewidth=2)
        ax0.errorbar(bincenters, hist_ref, yerr=yerr_ref, fmt=".", color=f"C{j+1}")
        ax0.set_ylabel(y_label, fontsize=14)
        ax0.legend(loc="best")
        ax0.grid(True, linewidth=0.5)

        # Apply optional settings
        if "xlim" in settings:
            ax0.set_xlim(*settings["xlim"])
        if "xscale" in settings:
            ax0.set_xscale(settings["xscale"])

        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, f"{var_name}.png"), bbox_inches='tight')
        plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description="Combine HDF5 data and plot variable distributions.")
    parser.add_argument('-i', '--main_base', required=True, help='Path to main data folder')
    parser.add_argument('-ir', '--ref_base', required=True, help='Path to reference data folder')
    #parser.add_argument('-d', '--main_dates', nargs='+', required=True, help='Dates for main data')
    #parser.add_argument('-dr', '--ref_dates', nargs='+', required=True, help='Dates for reference data')
    #parser.add_argument('-v', '--num_vars', type=int, help='Number of variables (26 for MF, 10 for FPF)')
    parser.add_argument('-o', '--plot_dir', type=str, help='Output directory for plots')
    #parser.add_argument('-l', '--main_label', type=str, help='Legend label for main dataset')
    #parser.add_argument('-lr', '--ref_label', type=str, help='Legend label for reference dataset')
    parser.add_argument('-n', '--normalize', action='store_true', help='Normalize histograms')

    args = parser.parse_args()
    print("Loading main dataset...")
    combined_years, labels_years, years = load_combined_data(args.main_base)
    print("Loading reference dataset...")
    combined_years_ref, labels_years_ref, years_ref = load_combined_data(args.ref_base)
    '''
    if labels_main != labels_ref:
        print("Warning: Labels mismatch between datasets!")
    #Use FPF labels in case FPF data are plotted
    '''
    print("Plotting distributions...")
    plot_distributions(
        combined_years,
        combined_years_ref[0],
        labels_years[0],
        years,
        args.plot_dir,
        normalize=args.normalize
    )
    print("Performing KS_tests...")
    KS_tests(
        combined_years,
        combined_years_ref[0],
        labels_years[0],
        years,
        args.plot_dir
    )

    print(f"Done! Plots saved to: {args.plot_dir}")

if __name__ == "__main__":
    main()
