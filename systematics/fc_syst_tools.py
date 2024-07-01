import numpy as np
import scipy as sp
from math import sqrt
from tqdm import tqdm
import matplotlib.pyplot as plt

#
# Optimal value of B given observation, syst
#
def Bhat(S, N, B0, sigsq):
    return (N*B0 + (N-S)*sigsq)/(N + sigsq)

#
# Chisquared for S given observation, syst
#
def chisq(S, N, B0, sigsq):
    B = Bhat(S, N, B0, sigsq)
    if sigsq>0:
        return (N - S - B)**2/N + (B-B0)**2/sigsq
    else:
        return (N - S - B)**2/N

#
# Throw Profiled FC Pseudo experiments to get a chisq distribution
#
def ProfileFC(S1, Ndata, B0, sigsq, Nexp):
    chisquares = []

    # Calculate the profile value of the nuisance parameter
    B1 = Bhat(S1, Ndata, B0, sigsq)

    for i in range(Nexp):
        # Create a pseudo experiment based on
        # test point S1, profile value B1
        N = np.random.poisson(S1+B1)

        # Perform FC calculation
        chiTrue = chisq(S1, N, B0, sigsq)
        chiBest = 0  # With only 1 data point, a perfect fit is always possible
        chisquares.append(chiTrue - chiBest)
        
    chisquares.sort()
    return chisquares

#
# Throw Highland-Cousins Pseudo experiments to get a chisq distribution
#
def HighlandCousins(S1, B0, sigsq, Nexp):
    chisquares = []

    sigma = sqrt(sigsq)

    for i in range(Nexp):
        # In HC, draw nuissance parameter from it's prior distribution
        B1 = np.random.normal(B0, sigma)

        # Create a pseudo experiment based on
        # test point S1, profile value B1
        N = np.random.poisson(S1+B1)

        # Perform FC calculation
        chiTrue = chisq(S1, N, B0, sigsq)
        chiBest = 0  # With only 1 data point, a perfect fit is always possible
        chisquares.append(chiTrue - chiBest)
        
    chisquares.sort()
    return chisquares

#
# Get the CDF at point x for given distribution chisq
#
def cdf(chisqs, x, err=True):
    N = len(chisqs)
    count = sum(value <= x for value in chisqs)
    if count == 0:
        if err:
            return 0.0, 0.
        else:
            return 0.0
    elif count == N:
        if err:
            return 1.0, 0.
        else:
            return 1.0
    else:
        lower_value = chisqs[count - 1]
        upper_value = chisqs[count]
        lower_cdf = count / N
        upper_cdf = (count + 1) / N
        p = lower_cdf + (upper_cdf - lower_cdf) * (x - lower_value) / (upper_value - lower_value)

    if not err:
        return p
    else:
        q = 1.-p
        err = sqrt(p*q/N)
        return p, err

#
# Calculate the inverse of the cumulative distribution function (CDF) of a list of data.
#
def critval(chisqs, cdf_value):
    """
    Calculate the inverse of the cumulative distribution function (CDF) of a list of data.

    Parameters:
        chisqs (list): List of positive floating point values.
        cdf_value (float): Target cumulative probability value.

    Returns:
        float: Input value corresponding to the given CDF value.
    """
    N = len(chisqs)   

    # Handle edge cases
    if cdf_value <= 0.0:
        return chisqs[0]
    elif cdf_value >= 1.0:
        return chisqs[-1]

    # Find the two consecutive CDF values that bracket the target CDF value
    for i in range(N - 1):
        lower_cdf = i / N
        upper_cdf = (i + 1) / N
        if lower_cdf <= cdf_value <= upper_cdf:
            #print(f"DEBUG: {lower_cdf}, {cdf_value}, {upper_cdf}")
            lower_value = chisqs[i]
            upper_value = chisqs[i + 1]
            interpolated_x = lower_value + (upper_value - lower_value) * (cdf_value - lower_cdf) / (upper_cdf - lower_cdf)
            #print(f"{lower_value} < {interpolated_x} < {upper_value}")
            return interpolated_x

#
# Generate p-values with Wilks, FC, and HC for given 
# parameters of the toy experiments.
#
def generate_pvals( systerr = 0.04,    # Systematic uncertainty
                    sigmarng = 4,
                    Nexp    = 100000,      # Number of pseudoexperiments to use
                    Strue   = 350.,
                    Btrue   = 150.,
                    B0      = 150.
                   ):

    # Test cases where nuissance parameter is correct
    pvals = {}
    pvals["FC"] = []
    pvals["Wilks"] = []
    pvals["HC"] = []

    # Test coverage at the true values
    S1 = Strue
    sigsq = (systerr*B0)**2

    HCchisquares = HighlandCousins(S1, B0, sigsq, Nexp)

    # Test every possible integer in +/- sigma range on count
    Ntrue = Strue+Btrue
    sigma = sqrt(Ntrue)
    bottom = int(Ntrue - sigmarng*sigma + 0.5)
    top    = int(Ntrue + sigmarng*sigma + 1.5)
    minprob = sp.stats.norm.pdf(bottom, Ntrue, sigma)

    
    for Ndata in tqdm(range(bottom, top)):
        chisqData = chisq(S1, Ndata, B0, sigsq)

        # These values are all the same for the same Ndata
        FCchisquares = ProfileFC(S1, Ndata, B0, sigsq, Nexp)
        pWilks = 1 - sp.stats.chi2.cdf(chisqData, df=1) 
        pFC    = 1 - cdf(FCchisquares, chisqData, err=False)
        pHC    =  1 - cdf(HCchisquares, chisqData, err=False)

        # Recereate full distro by adding entries in proportion to probability
        prob = sp.stats.norm.pdf(Ndata, Ntrue, sigma)
        deweight = int(prob / minprob + 0.5)
        pvals["Wilks"].append( (pWilks,deweight))
        pvals["FC"].append(    (pFC,   deweight))
        pvals["HC"].append(    (pHC,   deweight))

    return pvals

#
# Calculate how accurate the coverage of the true value is for
# the given sets of toy experiment p-values
# 
def calculate_accuracies(weighted_values, significance_levels):
    # Separate values and weights
    values, weights = zip(*weighted_values)
    values = np.array(values)
    weights = np.array(weights)

    coverage_accuracies = []    
    errors = []

    for level in significance_levels:
        intended_coverage = 1-level
        actual_coverage = np.average(np.array(values) > level, weights=weights)
        coverage_accuracy = (actual_coverage - intended_coverage)/level

        # Calculate the binomial error considering weights
        total = np.sum(weights)
        err = sqrt(actual_coverage * (1 - actual_coverage) / total)/level
        
        coverage_accuracies.append(coverage_accuracy)
        errors.append(err)
        #print(f'Method: {method:5}, Intended: {intended_coverage:.3f}, Actual Coverage: {actual_coverage:.3f}, Coverage Accuracy: {coverage_accuracy:.3f}')

    return coverage_accuracies, errors


def calculate_statistics(all_pvals, significance_levels, specific_zs, specific_significances):
    specific_accuracies = {}
    coverage_accuracies = {}
    specific_errors = {}
    errors = {}

    for case, pvals in all_pvals.items():
        specific_accuracies[case] = {}
        specific_errors[case] = {}
        coverage_accuracies[case] = {method: [] for method in pvals}
        errors[case] = {method: [] for method in pvals}

        for method, values in pvals.items():
            coverage_accuracies[case][method], errors[case][method] = calculate_accuracies(values, significance_levels)

            specific_accuracies[case][method] = {}
            specific_errors[case][method] = {}
            for z, x in zip(specific_zs, specific_significances):
                idx = (np.abs(significance_levels - x)).argmin()
                specific_accuracies[case][method][z] = coverage_accuracies[case][method][idx]
                specific_errors[case][method][z]     = errors[case][method][idx]

    return specific_accuracies, coverage_accuracies, specific_errors, errors

def calculate_mean_and_std(cases, method, z, ratios, specific_accuracies):
    from collections import defaultdict
    import numpy as np

    # Group accuracies by unique ratios
    grouped_accuracies = defaultdict(list)
    for case in cases:
        grouped_accuracies[ratios[case]].append(specific_accuracies[case][method][z])
    
    # Calculate mean and standard deviation for each unique ratio
    sorted_ratios = sorted(grouped_accuracies.keys())
    means = []
    std_devs = []
    for ratio in sorted_ratios:
        accuracies = grouped_accuracies[ratio]
        means.append(np.mean(accuracies))
        std_devs.append(np.std(accuracies))
    
    return sorted_ratios, means, std_devs


def apply_coverage_style(ax):
    # Root-style ticks
    ax.tick_params(axis='both', direction='in', which='both', top='true', right='true')
    ax.set_ylim(-1, 1)
    ax.set_xscale('log')  # Set the x-axis to a logarithmic scale
    ax.axhline(0, color='grey', lw=0.5)  # Adds a line at zero for reference
    for z in [1, 2, 3]:
        ax.axvline(2*sp.stats.norm.cdf(-1*z), color='grey', lw=0.5)  # Adds a line at zero for reference 


