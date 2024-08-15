import numpy as np
import scipy as sp

import os, pickle
import matplotlib.pyplot as plt
from math import sqrt
from fc_syst_tools import *
from pprint import pprint
from collections import defaultdict
import itertools

import pdb


def main():


    Strue = 350.
    Btrue = 150.

    cases = {}
    titles = {}
    positions = {}
    ratios = {}
    biases = {}
    for y, biasratio in enumerate([-1, -0.5, 0, 0.25, 0.5, 1, 2]):
        for x, errratio in enumerate([0.1, 0.2, 0.45]):
            B0 = Btrue # For these tests, get nuissance right

            staterrAbs = sqrt(Strue + 2*Btrue)
            systerrAbs = errratio*staterrAbs
            B0 = biasratio * systerrAbs + Btrue
            systerr = systerrAbs / Btrue

            biasstr = "correct"
            if biasratio > 0 or biasratio < 0:
                intratio = int(biasratio*10)
                if biasratio > 0:
                    signstr = "p"
                else:
                    signstr = "n"
                biasstr = f"{signstr}{abs(intratio):02d}off"

            # construct the short name
            case = f"{biasstr}_{Strue:g}_{Btrue:03g}_{errratio*100:g}"

            # construct the arguments
            cases[case] = {"systerr":systerr, "Strue": Strue, "Btrue":Btrue, "B0":B0}

            # Construct the title
            biasstr = "$ correct"
            if biasratio > 0 or biasratio < 0:
                biasstr = f" = {biasratio:+.1g}\sigma$"
            titles[case] = "$\sigma_{syst}/\sigma_{stat} = " + f"{errratio:.2}" + ", B_{0}" + biasstr
            ratios[case] = systerrAbs/staterrAbs
            biases[case] = biasratio

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
            all_pvals[case] = generate_pvals(**cargs, sigmarng=5.5)
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
    maxy, maxx = map(max, zip(*positions.values()))
    fig, axs = plt.subplots(maxy+1, maxx+1, figsize=(13,15), sharey=True)
    print("Shape",axs.shape)

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
    fig.savefig("coverage_syst_bias.png")
    fig.savefig("coverage_syst_bias.pdf")

    

    #
    # Summary Plots
    #

    # Step 1: Organize the data for plotting
    plot_data = {1: {}, 2: {}, 3: {}}
    for case, method_data in specific_accuracies.items():
        for method, z_data in method_data.items():
            for z, accuracy in z_data.items():
                if method not in plot_data[z]:
                    plot_data[z][method] = {}
                ratio = ratios[case]
                if ratio not in plot_data[z][method]:
                    plot_data[z][method][ratio] = {'biases': [], 'accuracies': []}
                plot_data[z][method][ratio]['biases'].append(biases[case])
                plot_data[z][method][ratio]['accuracies'].append(accuracy)


    # Step 1: Prepare a 3x3 grid of subplots
    fig, axs = plt.subplots(3, 3, figsize=(10, 7), sharex='col', sharey='row')
    unique_ratios = sorted(set(ratios.values()))  # Assuming ratios[case] gives the ratio for that case

    # Iterate over each subplot to populate it with appropriate data
    for j, ratio in enumerate(unique_ratios):
        for i, z in enumerate([1, 2, 3]):
            ax = axs[j][i]
            ax.tick_params(axis='both', direction='in', which='both', top=True, right=True)
            ax.set_ylim(-1, 1)  # Set the y-axis range for all plots

            # Check if data exists for this z-score and ratio
            if z in plot_data:
                for method in sorted(plot_data[z]):
                    if ratio in plot_data[z][method]:
                        data = plot_data[z][method][ratio]
                        sorted_data = sorted(zip(data['biases'], data['accuracies']))
                        biases_sorted, accuracies_sorted = zip(*sorted_data)
                        # Extracting errors corresponding to each bias point
                        errors = [specific_errors[case][method][z] for case in cases if ratios[case] == ratio]

                        # Plotting the data with error bands
                        ax.plot(biases_sorted, accuracies_sorted, label=labels[method],
                                color=method_colors[method], linestyle='-')  # Standard line style
                        ax.fill_between(biases_sorted,
                                        [a - e for a, e in zip(accuracies_sorted, errors)],
                                        [a + e for a, e in zip(accuracies_sorted, errors)],
                                        color=method_colors[method], alpha=0.2)  # Error band

            # Customize the subplot
            ax.axhline(0, color='gray', linewidth=1, linestyle='-')
            ax.axvline(0, color='gray', linewidth=1, linestyle='-')
            if j == 0:
                ax.set_title(f"Coverage Accuracy at {z}$\sigma$ C.L.")
            if i == 0:
                ax.set_ylabel("$\sigma_{syst}/\sigma_{stat} = "+f' {ratio}''$\nCoverage Accuracy')
            if j == len(unique_ratios) - 1:
                ax.set_xlabel(r"Bias in Estimated $B_{0}$ ($\sigma_{syst}$)")

    # Position the legend in the upper left subplot
    axs[0][0].legend(loc='upper left')

    fig.tight_layout()
    fig.savefig("accuracies_summary_systbias.png")
    fig.savefig("accuracies_summary_systbias.pdf")



if __name__ == "__main__":
    main()

