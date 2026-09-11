import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import argparse
import json
import csv
from scipy import stats
from scipy.special import gammaln
from scipy.ndimage import gaussian_filter
import emcee

# import custom libs
sys.path.append("../python")
from cuts import *
from plotting import *
from libs import *
from name_index_association import *


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-j', type=str, help='the json file with the parameters')
args = parser.parse_args()

# Load the parameters from the json file
with open(args.j) as f:
    parameters = json.load(f)

output_folder_base = parameters["folders"]["output_folder_base"]
if output_folder_base[-1] == "/":
    output_folder_base = output_folder_base[:-1]

if not os.path.exists(output_folder_base):
    os.makedirs(output_folder_base)

# ------------------ Load candidate counts from a data_data_comparison.py run -----------------------
data_data_output_folder_base = parameters["folders"]["data_data_output_folder_base"]
data_data_folder = data_data_output_folder_base + "_" + parameters["analysis"]["TP_RATE"] + "_cuts_" + str(parameters["analysis"]["apply_cuts"])

events_csv_path = os.path.join(data_data_folder, "events_passing_all_cuts.csv")
with open(events_csv_path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    is_signal_column = [row["IsSignal"] for row in reader]

N_on = is_signal_column.count("Signal")
N_off = is_signal_column.count("Background")

print(f"Read candidate counts from {events_csv_path}")
print(f"N_on (signal candidates) = {N_on}")
print(f"N_off (background candidates) = {N_off}")

# ------------------ Livetimes from the JSON config -----------------------
T_on = parameters["signal"]["run_parameters"]["spill_on"]["total_time"]
T_off = parameters["bkg"]["run_parameters"]["spill_off"]["total_time"]

print(f"T_on (hours) = {T_on}")
print(f"T_off (hours) = {T_off}")

# ------------------ Simple formula significance (Wilks, 1 dof) -----------------------
if N_off <= 0 or T_off <= 0:
    print(f"Cannot estimate expected background: N_off={N_off}, T_off={T_off}. Exiting.")
    sys.exit(1)

mu = N_off * (T_on / T_off)

if N_on > 0:
    test_statistic = max(0.0, 2 * (N_on * np.log(N_on / mu) + mu - N_on))
else:
    test_statistic = 0.0

# Wilks' theorem: this test statistic is asymptotically chi2 with 1 dof, i.e.
# the square of a standard normal. Working with sqrt(test_statistic) directly
# (instead of round-tripping through chi2.cdf/norm.ppf) avoids catastrophic
# loss of precision for large excesses, where 1 - p_value rounds to 1.0.
sign = 1.0 if N_on >= mu else -1.0
sigma_equivalent = sign * np.sqrt(test_statistic)
p_value = stats.norm.sf(sigma_equivalent)

print(f"Expected background in on-window (mu) = {mu}")
print(f"Test statistic = {test_statistic}")
print(f"P-value (one-sided, 1 dof) = {p_value}")
print(f"Significance (sigma equivalent) = {sigma_equivalent}")

# ------------------ Save results to file -----------------------
output_file = os.path.join(output_folder_base, "significance_results.txt")
with open(output_file, "w") as f:
    f.write(f"N_on (signal candidates): {N_on}\n")
    f.write(f"N_off (background candidates): {N_off}\n")
    f.write(f"T_on (hours): {T_on}\n")
    f.write(f"T_off (hours): {T_off}\n")
    f.write(f"Expected background in on-window (mu): {mu}\n")
    f.write(f"Test statistic: {test_statistic}\n")
    f.write(f"P-value (one-sided, 1 dof): {p_value}\n")
    f.write(f"Significance (sigma equivalent): {sigma_equivalent}\n")

print(f"Results saved to {output_file}")

# ------------------ MCMC sampling of the on/off Poisson likelihood -----------------------
mcmc_settings = parameters["settings"]
nwalkers = mcmc_settings["nwalkers"]
nsteps = mcmc_settings["nsteps"]
nburn = mcmc_settings["nburn"]
b_rate_limits = mcmc_settings["b_rate_limits"]
n_rate_limits = mcmc_settings["n_rate_limits"]


def loglike(theta, N_on, N_off, T_on, T_off):
    # theta: [R_beamOFF, R_nu], both in events/hour
    R_beamOFF, R_nu = theta
    expected_on = T_on * (R_beamOFF + R_nu)
    expected_off = T_off * R_beamOFF
    logl_on = N_on * np.log(expected_on + 1e-15) - expected_on - gammaln(N_on + 1)
    logl_off = N_off * np.log(expected_off + 1e-15) - expected_off - gammaln(N_off + 1)
    return logl_on + logl_off


def logprior(theta, b_rate_limits, n_rate_limits):
    R_beamOFF, R_nu = theta
    if b_rate_limits[0] < R_beamOFF < b_rate_limits[1] and n_rate_limits[0] < R_nu < n_rate_limits[1]:
        return 0.0
    return -np.inf


def logpost(theta, N_on, N_off, T_on, T_off, b_rate_limits, n_rate_limits):
    lp = logprior(theta, b_rate_limits, n_rate_limits)
    if not np.isfinite(lp):
        return -np.inf
    return lp + loglike(theta, N_on, N_off, T_on, T_off)


ndim = 2
initial_positions = np.array([b_rate_limits[0], n_rate_limits[0]]) + np.random.rand(nwalkers, ndim) * np.array(
    [b_rate_limits[1] - b_rate_limits[0], n_rate_limits[1] - n_rate_limits[0]]
)

sampler = emcee.EnsembleSampler(
    nwalkers, ndim, logpost, args=(N_on, N_off, T_on, T_off, b_rate_limits, n_rate_limits)
)
sampler.run_mcmc(initial_positions, nsteps, progress=True)

samples = sampler.get_chain()
flat_samples = sampler.get_chain(discard=nburn, flat=True)

R_beamOFF_samples = flat_samples[:, 0]
R_nu_samples = flat_samples[:, 1]

R_beamOFF_median = np.median(R_beamOFF_samples)
R_nu_median = np.median(R_nu_samples)
R_nu_std = np.std(R_nu_samples)
R_nu_16, R_nu_84 = np.percentile(R_nu_samples, [16, 84])
R_nu_2p5, R_nu_97p5 = np.percentile(R_nu_samples, [2.5, 97.5])

# Approximate "sigma" from the posterior, valid when the posterior is
# roughly Gaussian; the Phase-1 Wilks-theorem sigma above is the more
# rigorous number, this is a cross-check from the full sampled posterior.
mcmc_sigma_approx = R_nu_median / R_nu_std if R_nu_std > 0 else float("nan")

print(f"MCMC: R_beamOFF posterior median = {R_beamOFF_median}")
print(f"MCMC: R_nu posterior median = {R_nu_median} "
      f"(68% CI: [{R_nu_16}, {R_nu_84}], 95% CI: [{R_nu_2p5}, {R_nu_97p5}])")
print(f"MCMC: significance approx (R_nu_median / R_nu_std) = {mcmc_sigma_approx}")

# ---- walker trace diagnostic plot ----
param_names = ["R_beamOFF", "R_nu"]
fig, axes = plt.subplots(ndim, 1, figsize=(8, 6), sharex=True)
for i in range(ndim):
    for k in range(nwalkers):
        axes[i].plot(samples[:, k, i], alpha=0.3, color="k")
    axes[i].set_ylabel(param_names[i])
axes[-1].set_xlabel("Step")
plt.tight_layout()
plt.savefig(os.path.join(output_folder_base, "mcmc_chains.png"), dpi=150)
plt.clf()
plt.close(fig)

# ---- 2D posterior density with credible-region contours ----
H, xedges, yedges = np.histogram2d(
    R_beamOFF_samples, R_nu_samples, bins=(100, 100), range=[b_rate_limits, n_rate_limits]
)
H = H.T
# Smoothing before contouring removes the jaggedness of the finite-sample
# histogram; kept as a separate density so both the raw and smoothed contours
# can be produced and compared side by side.
H_smooth = gaussian_filter(H, sigma=1.5)

xcenters = 0.5 * (xedges[:-1] + xedges[1:])
ycenters = 0.5 * (yedges[:-1] + yedges[1:])
X, Y = np.meshgrid(xcenters, ycenters)
max_idx = np.unravel_index(np.argmax(H), H.shape)
R_beamOFF_mode = xcenters[max_idx[1]]
R_nu_mode = ycenters[max_idx[0]]

# Credible regions: for a k-dimensional Gaussian the probability enclosed
# within radius n (in units of standard deviations) is the CDF of a chi2
# distribution with k=ndim degrees of freedom evaluated at n^2. The radii
# below (1, 3, 5) are only used to pick a spread of enclosed fractions to
# draw contours at; the plot itself reports the enclosed probability, not
# the radius, to keep the vocabulary purely Bayesian (credible regions,
# not sigma significances).
contour_radii = [1, 3, 5]
enclosed_fractions = [stats.chi2.cdf(r ** 2, df=ndim) for r in contour_radii]


def plot_posterior_contours(contour_density, output_path, title_suffix=""):
    density_sorted = np.sort(contour_density.flatten())[::-1]
    density_cumsum = np.cumsum(density_sorted)
    density_cumsum /= density_cumsum[-1]

    thresholds = []
    for frac in enclosed_fractions:
        idx = np.where(density_cumsum >= frac)[0]
        thresholds.append(density_sorted[idx[0]] if len(idx) > 0 else density_sorted[-1])

    # matplotlib needs strictly increasing levels; larger sigma -> larger
    # enclosed fraction -> lower density threshold, so sort ascending and drop
    # duplicates (can happen at high sigma if the binning is too coarse).
    order = np.argsort(thresholds)
    contour_levels = []
    contour_labels = {}
    for i in order:
        level = thresholds[i]
        if contour_levels and level <= contour_levels[-1]:
            continue
        contour_levels.append(level)
        pct = enclosed_fractions[i] * 100
        pct_str = f"{pct:.4f}%" if pct >= 99.99 else f"{pct:.1f}%"
        contour_labels[level] = f"{pct_str} CR"

    plt.figure(figsize=(7, 6))
    plt.hist2d(R_beamOFF_samples, R_nu_samples, bins=(100, 100), range=[b_rate_limits, n_rate_limits])
    plt.colorbar(label="Counts")
    if len(contour_levels) > 0:
        CS = plt.contour(X, Y, contour_density, levels=contour_levels, colors="white", linewidths=1.0)
        # Place one label per level manually, on the midpoint of its longest
        # segment: clabel's automatic placement silently skips contours (e.g.
        # the innermost 1-sigma ring) that are too small to fit the text inline.
        manual_positions = []
        for level, segs in zip(CS.levels, CS.allsegs):
            if segs:
                longest = max(segs, key=len)
                manual_positions.append(tuple(longest[len(longest) // 2]))
        plt.clabel(CS, fmt=lambda v: contour_labels.get(v, f"{v:.0f}"), colors="white", fontsize=8,
                   manual=manual_positions if manual_positions else False)
    plt.scatter(R_beamOFF_mode, R_nu_mode, color="red", marker="x", s=80,
                label=f"Mode: R_beamOFF={R_beamOFF_mode:.3f}, R_nu={R_nu_mode:.3f}")
    plt.xlabel("R_beamOFF (events / hour)")
    plt.ylabel("R_nu (events / hour)")
    plt.title(f"MCMC posterior of (R_beamOFF, R_nu){title_suffix}")
    plt.legend(loc="upper right")
    plt.savefig(output_path, dpi=200)
    plt.clf()


plot_posterior_contours(
    H, os.path.join(output_folder_base, "mcmc_posterior_contours_raw.png"), title_suffix=" [raw]"
)
plot_posterior_contours(
    H_smooth, os.path.join(output_folder_base, "mcmc_posterior_contours_smoothed.png"), title_suffix=" [smoothed]"
)

# ---- save MCMC results ----
mcmc_output_file = os.path.join(output_folder_base, "mcmc_significance_results.txt")
with open(mcmc_output_file, "w") as f:
    f.write(f"nwalkers: {nwalkers}\n")
    f.write(f"nsteps: {nsteps}\n")
    f.write(f"nburn (discarded): {nburn}\n")
    f.write(f"R_beamOFF prior range (events/hour): {b_rate_limits}\n")
    f.write(f"R_nu prior range (events/hour): {n_rate_limits}\n")
    f.write(f"R_beamOFF posterior median: {R_beamOFF_median}\n")
    f.write(f"R_beamOFF posterior mode: {R_beamOFF_mode}\n")
    f.write(f"R_nu posterior median: {R_nu_median}\n")
    f.write(f"R_nu posterior mode: {R_nu_mode}\n")
    f.write(f"R_nu 68% credible interval: [{R_nu_16}, {R_nu_84}]\n")
    f.write(f"R_nu 95% credible interval: [{R_nu_2p5}, {R_nu_97p5}]\n")
    f.write(f"MCMC significance approx (R_nu_median / R_nu_std): {mcmc_sigma_approx}\n")
    f.write(f"Simple-formula significance (sigma equivalent, for comparison): {sigma_equivalent}\n")

print(f"MCMC results saved to {mcmc_output_file}")
