import numpy as np
import scipy as sp

import os, pickle
import matplotlib.pyplot as plt
from math import sqrt
from fc_syst_tools import *
from pprint import pprint
from collections import defaultdict

import pdb


def main():
    Strue = 350.
    #Btrue = 350.
    #ratio = 8.

    # Create cases, here a 5x5 grid of bkgd amount and syst/stat err
    cases = {}
    titles = {}
    positions = {}
    ratios = {}
    for x, Btrue in enumerate([50, 150, 250, 350, 450]):
        for y, errratio in enumerate([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45]):
            B0 = Btrue # For these tests, get nuissance right

            staterrAbs = sqrt(Strue + 2*Btrue)
            systerrAbs = errratio*staterrAbs
            biasratio = (B0 - Btrue) / systerrAbs
            systerr = systerrAbs / B0

            biasstr = "correct"
            if biasratio > 0 or biasratio < 0:
                intratio = int(biasratio*10)
                if biasratio > 0:
                    signstr = "p"
                else:
                    signstr = "n"
                biasstr = f"{signstr}{intratio:02d}off"

            # construct the short name
            case = f"{biasstr}_{Strue:g}_{Btrue:03g}_{errratio*100:g}"

            # construct the arguments
            cases[case] = {"systerr":systerr, "Strue": Strue, "Btrue":Btrue, "B0":B0}

            # Construct the title
            biasstr = "$ correct"
            if biasratio > 0 or biasratio < 0:
                biasstr = f" = {biasratio:+.1g}\sigma$"
            titles[case] = "$\sigma_{syst}/\sigma_{stat} = " + f"{errratio:.2}" + ", B_{true}/S_{true} = "+f"{B0/Strue:.2}$"


            ratios[case] = systerrAbs/staterrAbs

            # Set positions
            positions[case] = (y,x)


    labels = { "FC": "Profiled FC",
               "HC": "Highland-Cousins",
               "Wilks": "Wilks' Theorem" }

    # Define colors and line styles for visual distinction
    method_colors = {
        'Wilks': 'red',
        'HC': 'green',
        'FC': 'blue'
    }


    all_pvals = {}

    for case, cargs in cases.items():
        #print(case)
        filename = "pvals_"+case+".pkl"
        # Check if the file exists
        if os.path.exists(filename):
            with open(filename, 'rb') as file:
                all_pvals[case] = pickle.load(file)
            print("Loaded ",case," pvals from disk.")
        else:
            all_pvals[case] = generate_pvals(**cargs, NTests=40000)
            with open(filename, 'wb') as file:
                pickle.dump(all_pvals[case], file)
            print("Generated ",case," pvals and saved to disk.")

    print("")
    print("Done generating/loading. Time to plot.")


    #
    # Calculate the summary statistics
    #
    specific_zs            = [1, 2, 3]
    specific_significances = [2 * sp.stats.norm.cdf(-z) for z in specific_zs]
    specific_accuracies    = defaultdict(lambda: defaultdict(dict))  

    significance_levels = np.logspace(-3, -0.1, 75)

    specific_accuracies, coverage_accuracies, specific_errors, errors = calculate_statistics(
        all_pvals, significance_levels, specific_zs, specific_significances
    )


    #
    # Individual Plots
    #
    fig, axs = plt.subplots(9, 5, figsize=(13,15), sharey=True)

    for case in cases:
        ax = axs[positions[case]]
        ax.set_title(titles[case])

        for method, color in method_colors.items():
            accuracies = coverage_accuracies[case][method]
            ax.plot(significance_levels, accuracies, label=labels[method], color=color)

            lower_bound = [a - e for a, e in zip(accuracies, errors[case][method])]
            upper_bound = [a + e for a, e in zip(accuracies, errors[case][method])]
            ax.fill_between(significance_levels, lower_bound, upper_bound, alpha=0.2, color=color)

        ax.set_xlabel('Significance Level' if positions[case][0] == 1 else "")
        ax.set_ylabel('Coverage Accuracy' if positions[case][1] == 0 else "")
        if positions[case] == (0,0):
            ax.legend()
        apply_coverage_style(ax)

    plt.tight_layout()
    fig.savefig("coverage_bkgd_syst.png")
    fig.savefig("coverage_bkgd_syst.pdf")



    # Set global font size using rcParams
    plt.rcParams.update({
        'font.size': 12,        # Sets the base default font size
        'axes.titlesize': 14,   # Specific font size for subplot titles
        'axes.labelsize': 12,   # Font size for x and y labels
        'xtick.labelsize': 12,  # Font size for x-axis tick labels
        'ytick.labelsize': 12,  # Font size for y-axis tick labels
        'legend.fontsize': 12   # Font size for legend
    })



    #
    # Summary Plots
    #

    # Create a figure with 3 subplots (one for each z value)
    fig, axs = plt.subplots(1, 3, figsize=(10, 5), sharey=True)  # Adjust size as needed, with 3 horizontal subplots

    # Iterate over z_values to create each subplot
    for i, z in enumerate([1, 2, 3]):
        ax = axs[i]  # Select the corresponding subplot axis
        ax.tick_params(axis='both', direction='in', which='both', top='true', right='true')
        for method, color in sorted(method_colors.items()):
            # Use the calculation function to get means and std devs
            x_vals, means, std_devs = calculate_mean_and_std(cases, method, z, ratios, specific_accuracies)
            
            # Plotting mean values with error band
            ax.plot(x_vals, means, linestyle='-', color=color, label=labels[method] if i == 0 else "_nolegend_")
            ax.fill_between(x_vals, [m - s for m, s in zip(means, std_devs)], [m + s for m, s in zip(means, std_devs)],
                            color=color, alpha=0.2)

        # Add a horizontal line at y=0
        ax.axhline(0, color='gray', linewidth=1, linestyle='-')

        # Customize each subplot
        ax.set_title(f"Coverage Accuracy at {z}$\sigma$ C.L.")
        ax.set_xlabel(r"Ratio $\sigma_{syst}/\sigma_{stat}$")
        
        if i == 0:  # Only add the legend to the first subplot
            ax.set_ylabel('Coverage Accuracy')
            ax.legend()

    # Adjust layout to prevent overlap
    fig.tight_layout()
    fig.savefig("accuracies_summary_systbkgd.png")
    fig.savefig("accuracies_summary_systbkgd.pdf")



if __name__ == "__main__":
    main()

