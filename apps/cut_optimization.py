import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import uproot
import argparse
import json
import csv
from scipy import stats

import emcee

sys.path.append("../python")
from cuts import cut_thresholds, cuts
from plotting import *
from libs import *
from name_index_association import *



parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-j', type=str, help='the json file with the parameters')
args = parser.parse_args()

# Load the parameters from the json file
with open(args.j) as f:
    parameters = json.load(f)

nsteps = parameters["settings"]["nsteps"] 
nwalkers = parameters["settings"]["nwalkers"]
nburn = parameters["settings"]["nburn"]

events_target = parameters["events_target"] if "events_target" in parameters else None

output_folder_base = parameters["folders"]["output_folder_base"] + f"_nsteps_{nsteps}"
output_folder_opt = output_folder_base + "/optimized"
output_folder_def = output_folder_base + "/default"
if output_folder_base[-1] == "/":
    output_folder_base = output_folder_base[:-1]

if not os.path.exists(output_folder_base):
    os.makedirs(output_folder_base)
if not os.path.exists(output_folder_opt):
    os.makedirs(output_folder_opt)
if not os.path.exists(output_folder_def):
    os.makedirs(output_folder_def)

print(f"Output folder: {output_folder_base}")
data_events, data_true_events, data_full_events, data_weights, data_labels, data_filenames = load_from_folder(  parameters["folders"]["mother_folder_data"], 
                                                                                                    get_truth=parameters["data"]["GET_TRUTH"], 
                                                                                                    get_full_array=False, 
                                                                                                    use_combined=parameters["data"]["USE_COMBINED"], 
                                                                                                    weights_mode=parameters["data"]["WEIGHT_MODE"], 
                                                                                                    parameters=parameters["data"], 
                                                                                                    spill_status=parameters["data"]["SPILL_STATUS"])
  



mc_events, mc_true_events, mc_full_events, mc_weights, mc_labels, mc_filenames = load_from_folder(  parameters["folders"]["mother_folder_mc"], 
                                                                                                    get_truth=parameters["mc"]["GET_TRUTH"], 
                                                                                                    get_full_array=False, 
                                                                                                    use_combined=parameters["mc"]["USE_COMBINED"], 
                                                                                                    weights_mode=parameters["mc"]["WEIGHT_MODE"], 
                                                                                                    parameters=parameters["mc"], 
                                                                                                    spill_status=parameters["mc"]["SPILL_STATUS"])

print('------------')
print("TA cut:")
# Apply the trigger activity cut
print(len(mc_events), len(mc_weights), len(mc_true_events), len(mc_labels))
print(f"MC events before TA filter: {np.sum(mc_weights)}")

index_mc = np.where(mc_true_events[:, true_dict["triggerActivityFlag"]] == 1)[0]
mc_events = mc_events[index_mc]
mc_weights = mc_weights[index_mc]
mc_true_events = mc_true_events[index_mc]
mc_labels = mc_labels[index_mc]
mc_filenames = [mc_filenames[i] for i in index_mc]

print(f"MC events after TA filter: {np.sum(mc_weights)}")



print(f"Loaded {len(data_events)} off spill data events")
print(f"Loaded {len(mc_events)} MC events")
# ---------------------------------------
# ------------------- Data loaded, small parenthesis: make a simple csv file with the variables used for the selection -----------------------

interesting_variables_data = {
    "neutrino_pfp_in_slice": data_events[:, aggregate_dict["numberOfPFParticles"]],
    "vertexZ": data_events[:, aggregate_dict["vertexZ"]],
    "vertexY": data_events[:, aggregate_dict["vertexY"]],
    "num_pf_particles": data_events[:, aggregate_dict["numberOfPFParticles"]],
    "directionZ": np.maximum(data_events[:, aggregate_dict["directionZ"]], data_events[:, aggregate_dict["directionZ2"]]),
    "energyDepositedInFirst10cm": data_events[:, aggregate_dict["energyDepositedInFirst10cm"]],
    "energyDepositedInFifth10cm": data_events[:, aggregate_dict["energyDepositedInFifth10cm"]],
    "energyDepositedInFifteenth10cm": data_events[:, aggregate_dict["energyDepositedInFifteenth10cm"]],
    "ROI_Z_size": data_events[:, aggregate_dict["zROIEnd"]] - data_events[:, aggregate_dict["zROIStart"]],
    "ROI_Z_starting_point_close_to_vertexZ": np.abs(data_events[:, aggregate_dict["vertexZ"]] - data_events[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"]),
    "Neutrino_Tail_Length_Density": data_events[:, aggregate_dict["lengthOfMuonTrack"]]/(data_events[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"])
}

interesting_variables_mc = {
    "neutrino_pfp_in_slice": mc_events[:, aggregate_dict["numberOfPFParticles"]],
    "vertexZ": mc_events[:, aggregate_dict["vertexZ"]],
    "vertexY": mc_events[:, aggregate_dict["vertexY"]],
    "num_pf_particles": mc_events[:, aggregate_dict["numberOfPFParticles"]],
    "directionZ": np.maximum(mc_events[:, aggregate_dict["directionZ"]], mc_events[:, aggregate_dict["directionZ2"]]),
    "energyDepositedInFirst10cm": mc_events[:, aggregate_dict["energyDepositedInFirst10cm"]],
    "energyDepositedInFifth10cm": mc_events[:, aggregate_dict["energyDepositedInFifth10cm"]],
    "energyDepositedInFifteenth10cm": mc_events[:, aggregate_dict["energyDepositedInFifteenth10cm"]],
    "ROI_Z_size": mc_events[:, aggregate_dict["zROIEnd"]] - mc_events[:, aggregate_dict["zROIStart"]],
    "ROI_Z_starting_point_close_to_vertexZ": np.abs(mc_events[:, aggregate_dict["vertexZ"]] - mc_events[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"]),
    "Neutrino_Tail_Length_Density": mc_events[:, aggregate_dict["lengthOfMuonTrack"]]/(mc_events[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"])
}

# # save these variables in a csv file for later use in the optimization
# with open(os.path.join(output_folder_base, "interesting_variables_data.csv"), "w", newline="") as f:
#     writer = csv.writer(f)
#     writer.writerow(interesting_variables_data.keys())
#     for i in range(len(data_events)):
#         writer.writerow([interesting_variables_data[key][i] for key in interesting_variables_data.keys()])

# with open(os.path.join(output_folder_base, "interesting_variables_mc.csv"), "w", newline="") as f:
#     writer = csv.writer(f)
#     writer.writerow(interesting_variables_mc.keys())
#     for i in range(len(mc_events)):
#         writer.writerow([interesting_variables_mc[key][i] for key in interesting_variables_mc.keys()])



# ---------------------------------------
# ------------------- Data loaded, define functions -----------------------

# def apply_all_cuts( sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
#                     bkg_events, bkg_weights, bkg_labels, bkg_filenames, cut_thresholds, skip_cut=None):
#     # Helper to apply all cuts except skip_cut (by name)
#     # Returns: sig_events, sig_weights, sig_true_events, sig_labels, bkg_events, bkg_weights, bkg_labels
#     # NOTE: cut_thresholds is captured via default argument to avoid late-binding closure issues
#     def _make_cuts(ct):
#         return [
#             ("neutrino_pfp_in_slice", lambda s, b, ct=ct: (
#                 np.where(s[:, aggregate_dict["numberOfPFParticles"]] >= ct["neutrino_pfp_in_slice"])[0],
#                 np.where(b[:, aggregate_dict["numberOfPFParticles"]] >= ct["neutrino_pfp_in_slice"])[0]
#             )),
#             ("vertex_fiducial_volume", lambda s, b, ct=ct: (
#                 np.where(
#                     (s[:, aggregate_dict["vertexZ"]] >= ct["vertex_z_min"]) &
#                     (s[:, aggregate_dict["vertexY"]] <= ct["vertex_y_max"])
#                 )[0],
#                 np.where(
#                     (b[:, aggregate_dict["vertexZ"]] >= ct["vertex_z_min"]) &
#                     (b[:, aggregate_dict["vertexY"]] <= ct["vertex_y_max"])
#                 )[0]
#             )),
#             ("daughter_particles", lambda s, b, ct=ct: (
#                 np.where(s[:, aggregate_dict["numberOfPFParticles"]] >= ct["num_pf_particles_min"])[0],
#                 np.where(b[:, aggregate_dict["numberOfPFParticles"]] >= ct["num_pf_particles_min"])[0]
#             )),
#             ("max_directionZ", lambda s, b, ct=ct: (
#                 np.where(np.maximum(
#                     s[:, aggregate_dict["directionZ"]],
#                     s[:, aggregate_dict["directionZ2"]]
#                 ) > ct["direction_z_max"])[0],
#                 np.where(np.maximum(
#                     b[:, aggregate_dict["directionZ"]],
#                     b[:, aggregate_dict["directionZ2"]]
#                 ) > ct["direction_z_max"])[0]
#             )),
#             ("energy_first10cm", lambda s, b, ct=ct: (
#                 np.where(
#                     (s[:, aggregate_dict["vertexZ"]] >= 450) |
#                     (s[:, aggregate_dict["energyDepositedInFirst10cm"]] > ct["energy_first10cm_min"])
#                 )[0],
#                 np.where(
#                     (b[:, aggregate_dict["vertexZ"]] >= 450) |
#                     (b[:, aggregate_dict["energyDepositedInFirst10cm"]] > ct["energy_first10cm_min"])
#                 )[0]
#             )),
#             ("energy_fifth10cm", lambda s, b, ct=ct: (
#                 np.where(
#                     (s[:, aggregate_dict["vertexZ"]] >= 410) |
#                     (s[:, aggregate_dict["energyDepositedInFifth10cm"]] > ct["energy_fifth10cm_min"])
#                 )[0],
#                 np.where(
#                     (b[:, aggregate_dict["vertexZ"]] >= 410) |
#                     (b[:, aggregate_dict["energyDepositedInFifth10cm"]] > ct["energy_fifth10cm_min"])
#                 )[0]
#             )),
#             ("energy_fifteenth10cm", lambda s, b, ct=ct: (
#                 np.where(
#                     (s[:, aggregate_dict["vertexZ"]] >= 310) |
#                     (s[:, aggregate_dict["energyDepositedInFifteenth10cm"]] > ct["energy_fifteenth10cm_min"])
#                 )[0],
#                 np.where(
#                     (b[:, aggregate_dict["vertexZ"]] >= 310) |
#                     (b[:, aggregate_dict["energyDepositedInFifteenth10cm"]] > ct["energy_fifteenth10cm_min"])
#                 )[0]
#             )),
#             ("ROI_Z_size", lambda s, b, ct=ct: (
#                 np.where(s[:, aggregate_dict["zROIEnd"]] - s[:, aggregate_dict["zROIStart"]] > ct["roi_z_size_min"])[0],
#                 np.where(b[:, aggregate_dict["zROIEnd"]] - b[:, aggregate_dict["zROIStart"]] > ct["roi_z_size_min"])[0]
#             )),
#             ("ROI_Z_starting_point_close_to_vertexZ", lambda s, b, ct=ct: (
#                 np.where(np.abs(s[:, aggregate_dict["vertexZ"]] - s[:, aggregate_dict["zROIStart"]]*ct["roi_z_scale_factor"]) < ct["roi_z_vertex_distance_max"])[0],
#                 np.where(np.abs(b[:, aggregate_dict["vertexZ"]] - b[:, aggregate_dict["zROIStart"]]*ct["roi_z_scale_factor"]) < ct["roi_z_vertex_distance_max"])[0]
#             )),
#             ("Neutrino_Tail_Length_Density", lambda s, b, ct=ct: (
#                 np.where(s[:, aggregate_dict["lengthOfMuonTrack"]]/(s[:, aggregate_dict["zROIStart"]]*ct["roi_z_scale_factor"]) < ct["tail_length_density_max"])[0],
#                 np.where(b[:, aggregate_dict["lengthOfMuonTrack"]]/(b[:, aggregate_dict["zROIStart"]]*ct["roi_z_scale_factor"]) < ct["tail_length_density_max"])[0]
#             ))
#         ]
#     local_cuts = _make_cuts(cut_thresholds)

#     # Work on copies
#     s_ev, s_w, s_true, s_lab, s_fnames = sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames
#     b_ev, b_w, b_lab, b_fnames = bkg_events, bkg_weights, bkg_labels, bkg_filenames

#     for name, cutfunc in local_cuts:
#         if skip_cut is not None and name == skip_cut:
#             continue
#         idx_s, idx_b = cutfunc(s_ev, b_ev)
#         s_ev = s_ev[idx_s]
#         s_w = s_w[idx_s]
#         if len(s_true) > 0:
#             s_true = s_true[idx_s]
#         s_lab = s_lab[idx_s]
#         s_fnames = [s_fnames[i] for i in idx_s]
#         b_ev = b_ev[idx_b]
#         b_w = b_w[idx_b]
#         b_lab = b_lab[idx_b]
#         b_fnames = [b_fnames[i] for i in idx_b]
#         # print(f"After cut {name}: N signal events = {len(s_ev)}, N background events = {len(b_ev)}")

#     return s_ev, s_w, s_true, s_lab, s_fnames, b_ev, b_w, b_lab, b_fnames

def apply_all_cuts( sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
                    bkg_events, bkg_weights, bkg_labels, bkg_filenames, cut_thresholds, skip_cut=None):
    # Helper to apply all cuts except skip_cut (by name)
    # Returns: sig_events, sig_weights, sig_true_events, sig_labels, bkg_events, bkg_weights, bkg_labels
    # NOTE: cut_thresholds is captured via default argument to avoid late-binding closure issues
    def _make_cuts(ct):
        return [
            ("neutrino_pfp_in_slice", lambda s, b, ct=ct: (
                np.where(s[:, aggregate_dict["numberOfPFParticles"]] >= ct["neutrino_pfp_in_slice"])[0],
                np.where(b[:, aggregate_dict["numberOfPFParticles"]] >= ct["neutrino_pfp_in_slice"])[0]
            )),
            ("vertex_fiducial_volume", lambda s, b, ct=ct: (
                np.where(
                    (s[:, aggregate_dict["vertexZ"]] >= ct["vertex_z_min"]) &
                    (s[:, aggregate_dict["vertexY"]] <= ct["vertex_y_max"])
                )[0],
                np.where(
                    (b[:, aggregate_dict["vertexZ"]] >= ct["vertex_z_min"]) &
                    (b[:, aggregate_dict["vertexY"]] <= ct["vertex_y_max"])
                )[0]
            )),
            ("daughter_particles", lambda s, b, ct=ct: (
                np.where(s[:, aggregate_dict["numberOfPFParticles"]] >= ct["num_pf_particles_min"])[0],
                np.where(b[:, aggregate_dict["numberOfPFParticles"]] >= ct["num_pf_particles_min"])[0]
            )),
            ("max_directionZ", lambda s, b, ct=ct: (
                np.where(np.maximum(
                    s[:, aggregate_dict["directionZ"]],
                    s[:, aggregate_dict["directionZ2"]]
                ) > ct["direction_z_max"])[0],
                np.where(np.maximum(
                    b[:, aggregate_dict["directionZ"]],
                    b[:, aggregate_dict["directionZ2"]]
                ) > ct["direction_z_max"])[0]
            )),
            ("energy_first10cm", lambda s, b, ct=ct: (
                np.where(
                    (s[:, aggregate_dict["energyDepositedInFirst10cm"]] > ct["energy_first10cm_min"])
                )[0],
                np.where(
                    (b[:, aggregate_dict["energyDepositedInFirst10cm"]] > ct["energy_first10cm_min"])
                )[0]
            )),
            ("energy_fifth10cm", lambda s, b, ct=ct: (
                np.where(
                    (s[:, aggregate_dict["energyDepositedInFifth10cm"]] > ct["energy_fifth10cm_min"])
                )[0],
                np.where(
                    (b[:, aggregate_dict["energyDepositedInFifth10cm"]] > ct["energy_fifth10cm_min"])
                )[0]
            )),
            ("energy_fifteenth10cm", lambda s, b, ct=ct: (
                np.where(
                    (s[:, aggregate_dict["energyDepositedInFifteenth10cm"]] > ct["energy_fifteenth10cm_min"])
                )[0],
                np.where(
                    (b[:, aggregate_dict["energyDepositedInFifteenth10cm"]] > ct["energy_fifteenth10cm_min"])
                )[0]
            )),
            ("ROI_Z_size", lambda s, b, ct=ct: (
                np.where(s[:, aggregate_dict["zROIEnd"]] - s[:, aggregate_dict["zROIStart"]] > ct["roi_z_size_min"])[0],
                np.where(b[:, aggregate_dict["zROIEnd"]] - b[:, aggregate_dict["zROIStart"]] > ct["roi_z_size_min"])[0]
            )),
            ("ROI_Z_starting_point_close_to_vertexZ", lambda s, b, ct=ct: (
                np.where(np.abs(s[:, aggregate_dict["vertexZ"]] - s[:, aggregate_dict["zROIStart"]]*ct["roi_z_scale_factor"]) < ct["roi_z_vertex_distance_max"])[0],
                np.where(np.abs(b[:, aggregate_dict["vertexZ"]] - b[:, aggregate_dict["zROIStart"]]*ct["roi_z_scale_factor"]) < ct["roi_z_vertex_distance_max"])[0]
            )),
            ("Neutrino_Tail_Length_Density", lambda s, b, ct=ct: (
                np.where(s[:, aggregate_dict["lengthOfMuonTrack"]]/(s[:, aggregate_dict["zROIStart"]]*ct["roi_z_scale_factor"]) < ct["tail_length_density_max"])[0],
                np.where(b[:, aggregate_dict["lengthOfMuonTrack"]]/(b[:, aggregate_dict["zROIStart"]]*ct["roi_z_scale_factor"]) < ct["tail_length_density_max"])[0]
            ))
        ]
    local_cuts = _make_cuts(cut_thresholds)

    # Work on copies
    s_ev, s_w, s_true, s_lab, s_fnames = sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames
    b_ev, b_w, b_lab, b_fnames = bkg_events, bkg_weights, bkg_labels, bkg_filenames

    for name, cutfunc in local_cuts:
        if skip_cut is not None and name == skip_cut:
            continue
        idx_s, idx_b = cutfunc(s_ev, b_ev)
        s_ev = s_ev[idx_s]
        s_w = s_w[idx_s]
        if len(s_true) > 0:
            s_true = s_true[idx_s]
        s_lab = s_lab[idx_s]
        s_fnames = [s_fnames[i] for i in idx_s]
        b_ev = b_ev[idx_b]
        b_w = b_w[idx_b]
        b_lab = b_lab[idx_b]
        b_fnames = [b_fnames[i] for i in idx_b]
        # print(f"After cut {name}: N signal events = {len(s_ev)}, N background events = {len(b_ev)}")

    return s_ev, s_w, s_true, s_lab, s_fnames, b_ev, b_w, b_lab, b_fnames


def log_likelihood(theta, mc_events, mc_weights, data_events, data_weights):
    # Unpack the parameters
    vertex_z_min, vertex_y_max, num_pf_particles_min, direction_z_max, energy_first10cm_min, energy_fifth10cm_min, energy_fifteenth10cm_min, roi_z_size_min, roi_z_vertex_distance_max, tail_length_density_max = theta    
    # Apply the cuts to the MC and off spill data events
    cut_thresholds = {
        "neutrino_pfp_in_slice": 0,
        "vertex_z_min": vertex_z_min,
        "vertex_y_max": vertex_y_max,
        "num_pf_particles_min": num_pf_particles_min,
        "direction_z_max": direction_z_max,
        "energy_first10cm_min": energy_first10cm_min,
        "energy_fifth10cm_min": energy_fifth10cm_min,
        "energy_fifteenth10cm_min": energy_fifteenth10cm_min,
        "roi_z_size_min": roi_z_size_min, 
        "roi_z_vertex_distance_max": roi_z_vertex_distance_max,
        "tail_length_density_max": tail_length_density_max,
        "roi_z_scale_factor": 460/100
    }
    
    data_before_cuts = np.sum(data_weights)
    mc_before_cuts = np.sum(mc_weights)
    
    mc_events_cut, mc_weights_cut, mc_true_events_cut, mc_labels_cut, mc_filenames_cut, data_events_cut, data_weights_cut, data_labels_cut, data_filenames_cut = apply_all_cuts(
        mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames,
        data_events, data_weights, data_labels, data_filenames,
        cut_thresholds=cut_thresholds
    )

    if len(data_weights_cut) < 5 or np.sum(mc_weights_cut) < 3:
        return -1e10 # if too few events pass the cuts, return a very low likelihood to discourage the optimizer from going in that direction

    efficiency = np.sum(mc_weights_cut) / mc_before_cuts if mc_before_cuts > 0 else 0
    purity = np.sum(mc_weights_cut) / (np.sum(mc_weights_cut) + np.sum(data_weights_cut)) if (np.sum(mc_weights_cut) + np.sum(data_weights_cut)) > 0 else 0

    T_on = 7.65 + 8.88
    T_off = 44.688
    n_on = (np.sum(mc_weights_cut)+np.sum(data_weights_cut))*T_on
    n_off = np.sum(data_weights_cut)*T_off
    mu = n_off * (T_on / T_off)
    significance = 2 * (n_on * np.log(n_on / mu) + mu - n_on)

    if not np.isfinite(significance):
        significance = 0

    return significance


def log_likelihood_and_values(theta, mc_events, mc_weights, data_events, data_weights):
    # Unpack the parameters
    vertex_z_min, vertex_y_max, num_pf_particles_min, direction_z_max, energy_first10cm_min, energy_fifth10cm_min, energy_fifteenth10cm_min, roi_z_size_min, roi_z_vertex_distance_max, tail_length_density_max = theta    
    # Apply the cuts to the MC and off spill data events
    cut_thresholds = {
        "neutrino_pfp_in_slice": 0,
        "vertex_z_min": vertex_z_min,
        "vertex_y_max": vertex_y_max,
        "num_pf_particles_min": num_pf_particles_min,
        "direction_z_max": direction_z_max,
        "energy_first10cm_min": energy_first10cm_min,
        "energy_fifth10cm_min": energy_fifth10cm_min,
        "energy_fifteenth10cm_min": energy_fifteenth10cm_min,
        "roi_z_size_min": roi_z_size_min, 
        "roi_z_vertex_distance_max": roi_z_vertex_distance_max,
        "tail_length_density_max": tail_length_density_max,
        "roi_z_scale_factor": 460/100
    }
    

    data_before_cuts = np.sum(data_weights)
    mc_before_cuts = np.sum(mc_weights)
    
    mc_events_cut, mc_weights_cut, mc_true_events_cut, mc_labels_cut, mc_filenames_cut, data_events_cut, data_weights_cut, data_labels_cut, data_filenames_cut = apply_all_cuts(
        mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames,
        data_events, data_weights, data_labels, data_filenames,
        cut_thresholds=cut_thresholds
    )

    efficiency = np.sum(mc_weights_cut) / mc_before_cuts if mc_before_cuts > 0 else 0
    purity = np.sum(mc_weights_cut) / (np.sum(mc_weights_cut) + np.sum(data_weights_cut)) if (np.sum(mc_weights_cut) + np.sum(data_weights_cut)) > 0 else 0
    f1_score = 2 * (efficiency * purity) / (efficiency + purity) if (efficiency + purity) > 0 else 0

    # compute significance as well
    T_on = 7.65 + 8.88
    T_off = 44.688
    n_on = (np.sum(mc_weights_cut)+np.sum(data_weights_cut))*T_on
    n_off = np.sum(data_weights_cut)*T_off
    mu = n_off * (T_on / T_off)
    significance = 2 * (n_on * np.log(n_on / mu) + mu - n_on)

    if not np.isfinite(significance):
        significance = 0



    total_mc_events_after_cuts = np.sum(mc_weights_cut)
    total_data_events_after_cuts = np.sum(data_weights_cut)
    return f1_score, efficiency, purity, total_mc_events_after_cuts, total_data_events_after_cuts, significance


# threshold_boundaries = {
#     "vertex_z_min": (5, 40),
#     "vertex_y_max": (500, 600),
#     "num_pf_particles_min": (2, 10),
#     "direction_z_max": (0.93, 0.99),
#     "energy_first10cm_min": (500, 5000),
#     "energy_fifth10cm_min": (5000, 12000),
#     "energy_fifteenth10cm_min": (500, 5000),
#     "roi_z_size_min": (8, 20), 
#     "roi_z_vertex_distance_max": (15, 35),
#     "tail_length_density_max": (0.05, 0.35),
# }

threshold_boundaries = {
    "vertex_z_min": (0, 460),
    "vertex_y_max": (0, 600),
    "num_pf_particles_min": (0, 15),
    "direction_z_max": (0.95, 0.99),
    "energy_first10cm_min": (0, 10000),
    "energy_fifth10cm_min": (0, 10000),
    "energy_fifteenth10cm_min": (0, 10000),
    "roi_z_size_min": (0, 50), 
    "roi_z_vertex_distance_max": (0, 50),
    "tail_length_density_max": (0, 1),
}


def logprior(args):
    vertex_z_min, vertex_y_max, num_pf_particles_min, direction_z_max, energy_first10cm_min, energy_fifth10cm_min, energy_fifteenth10cm_min, roi_z_size_min, roi_z_vertex_distance_max, tail_length_density_max = args
    if not (threshold_boundaries["vertex_z_min"][0] <= vertex_z_min <= threshold_boundaries["vertex_z_min"][1]):
        return -np.inf
    if not (threshold_boundaries["vertex_y_max"][0] <= vertex_y_max <= threshold_boundaries["vertex_y_max"][1]):
        return -np.inf
    if not (threshold_boundaries["num_pf_particles_min"][0] <= num_pf_particles_min <= threshold_boundaries["num_pf_particles_min"][1]):
        return -np.inf
    if not (threshold_boundaries["direction_z_max"][0] <= direction_z_max <= threshold_boundaries["direction_z_max"][1]):
        return -np.inf
    if not (threshold_boundaries["energy_first10cm_min"][0] <= energy_first10cm_min <= threshold_boundaries["energy_first10cm_min"][1]):
        return -np.inf
    if not (threshold_boundaries["energy_fifth10cm_min"][0] <= energy_fifth10cm_min <= threshold_boundaries["energy_fifth10cm_min"][1]):
        return -np.inf
    if not (threshold_boundaries["energy_fifteenth10cm_min"][0] <= energy_fifteenth10cm_min <= threshold_boundaries["energy_fifteenth10cm_min"][1]):
        return -np.inf
    if not (threshold_boundaries["roi_z_size_min"][0] <= roi_z_size_min <= threshold_boundaries["roi_z_size_min"][1]):
        return -np.inf
    if not (threshold_boundaries["roi_z_vertex_distance_max"][0] <= roi_z_vertex_distance_max <= threshold_boundaries["roi_z_vertex_distance_max"][1]):
        return -np.inf
    if not (threshold_boundaries["tail_length_density_max"][0] <= tail_length_density_max <= threshold_boundaries["tail_length_density_max"][1]):
        return -np.inf
    return 0


def logpost(args, mc_events, mc_weights, data_events, data_weights):
    if not np.isfinite(logprior(args)):
        return -np.inf
    return logprior(args) + log_likelihood(args, mc_events, mc_weights, data_events, data_weights)

default_cut_thresholds = {
    "neutrino_pfp_in_slice": 0,
    "vertex_z_min": 20,
    "vertex_y_max": 550,
    "num_pf_particles_min": 6,
    "direction_z_max": 0.97,
    "energy_first10cm_min": 3000,
    "energy_fifth10cm_min": 10000,
    "energy_fifteenth10cm_min": 3000,
    "roi_z_size_min": 15,
    "roi_z_vertex_distance_max": 20,
    "tail_length_density_max": 0.2,
    "roi_z_scale_factor": 460/100
}

# Initial guess for the parameters
manual_guess_parameters = [
    cut_thresholds["vertex_z_min"],
    cut_thresholds["vertex_y_max"],
    cut_thresholds["num_pf_particles_min"],
    cut_thresholds["direction_z_max"],
    cut_thresholds["energy_first10cm_min"],
    cut_thresholds["energy_fifth10cm_min"],
    cut_thresholds["energy_fifteenth10cm_min"],
    cut_thresholds["roi_z_size_min"],
    cut_thresholds["roi_z_vertex_distance_max"],
    cut_thresholds["tail_length_density_max"]
]
print("Initial parameters:", manual_guess_parameters)


# ---------------------------------------
# ------------------- Start optimization -----------------------


# Run the MCMC sampler
ndim = len(manual_guess_parameters)
initial_params = np.random.rand(nwalkers, len(manual_guess_parameters)) * [threshold_boundaries[key][1] - threshold_boundaries[key][0] for key in threshold_boundaries.keys()] + [threshold_boundaries[key][0] for key in threshold_boundaries.keys()]
sampler = emcee.EnsembleSampler(nwalkers, ndim, logpost, args=(mc_events, mc_weights, data_events, data_weights))
sampler.run_mcmc(initial_params, nsteps, progress=True)

samples = sampler.get_chain(discard=nburn, flat=False)

fig, axes = plt.subplots(figsize = (10, 20), ncols = 1, nrows = ndim)

# get the histogram of the log likelihood values and plot it
scores = sampler.get_log_prob(flat=True)
plt.figure()
plt.hist(scores, bins=50)
plt.xlabel("Log likelihood")
plt.ylabel("Number of samples")
plt.title("Histogram of log likelihood values")
plt.savefig(os.path.join(output_folder_base, "log_likelihood_histogram.png"))
plt.close()
# plot the sequence of likelihood values for each walker
plt.figure()
plt.plot(sampler.get_log_prob(), alpha=0.7)
plt.xlabel("Step number")
plt.ylabel("Log likelihood")
plt.title("Log likelihood values for each walker")
plt.savefig(os.path.join(output_folder_base, "log_likelihood_chains.png"))
plt.close()





for i in range(ndim):
    for k in range(nwalkers):
        axes[i].plot(samples[: , k, i], alpha=0.7)
    axes[i].set_xlabel("Step number")

plt.savefig(os.path.join(output_folder_base, "mcmc_chains.png"))
plt.close()

# Get the best parameters
best_params = sampler.get_chain(flat=True)[np.argmax(sampler.get_log_prob(flat=True))]
print("Best parameters:", best_params)

# Save results to a file
with open(os.path.join(output_folder_base, "cut_optimization_results.txt"), "w") as f:
    f.write("Best parameters:\n")
    for i, key in enumerate(threshold_boundaries.keys()):
        f.write(f"{key}: {best_params[i]}\n")

    # write total efficiency, purity and score
    cut_thresholds_optimized = {
        "neutrino_pfp_in_slice": 0,
        "vertex_z_min": best_params[0],
        "vertex_y_max": best_params[1],
        "num_pf_particles_min": best_params[2],
        "direction_z_max": best_params[3],
        "energy_first10cm_min": best_params[4],
        "energy_fifth10cm_min": best_params[5],
        "energy_fifteenth10cm_min": best_params[6],
        "roi_z_size_min": best_params[7], 
        "roi_z_vertex_distance_max": best_params[8],
        "tail_length_density_max": best_params[9],
        "roi_z_scale_factor": 460/100
    }
    mc_events_cut, mc_weights_cut, mc_true_events_cut, mc_labels_cut, mc_filenames_cut, data_events_cut, data_weights_cut, data_labels_cut, data_filenames_cut = apply_all_cuts(
        mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames,
        data_events, data_weights, data_labels, data_filenames,
        cut_thresholds=cut_thresholds_optimized
    )

    f1_score, efficiency, purity, total_mc_events_after_cuts, total_data_events_after_cuts, significance = log_likelihood_and_values(best_params, mc_events, mc_weights, data_events, data_weights)

    f.write("\n")
    f.write("Number of MC events after cuts normalized: " + str(np.sum(mc_weights_cut)) + "\n")
    f.write("Number of off spill data events after cuts normalized: " + str(np.sum(data_weights_cut)) + "\n")
    f.write(f"MC/off spill data ratio after cuts: " + str(np.sum(mc_weights_cut)/np.sum(data_weights_cut)) + "\n")
    f.write(f"Efficiency: {efficiency}\n")
    f.write(f"Purity: {purity}\n")
    f.write(f"F1 Score: {f1_score}\n")
    f.write(f"Significance: {significance}\n")  
    f.write("Number of MC events after cuts absolute: " + str(len(mc_weights_cut)) + "\n")
    f.write("Number of off spill data events after cuts absolute: " + str(len(data_weights_cut)) + "\n")

    # compare with the default cuts
    cut_thresholds = {
        "neutrino_pfp_in_slice": 0,
        "vertex_z_min": 20,
        "vertex_y_max": 550,
        "num_pf_particles_min": 6,
        "direction_z_max": 0.97,
        "energy_first10cm_min": 3000,
        "energy_fifth10cm_min": 10000,
        "energy_fifteenth10cm_min": 3000,
        "roi_z_size_min": 15,
        "roi_z_vertex_distance_max": 20,
        "tail_length_density_max": 0.2,
        "roi_z_scale_factor": 460/100
    }
    mc_events_cut, mc_weights_cut, mc_true_events_cut, mc_labels_cut, mc_filenames_cut, data_events_cut, data_weights_cut, data_labels_cut, data_filenames_cut = apply_all_cuts(
        mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames,
        data_events, data_weights, data_labels, data_filenames,
        cut_thresholds=cut_thresholds
    )
    # efficiency = np.sum(mc_weights_cut) / np.sum(mc_weights) if np.sum(mc_weights) > 0 else 0
    # purity = np.sum(mc_weights_cut) / (np.sum(mc_weights_cut) + np.sum(data_weights_cut)) if (np.sum(mc_weights_cut) + np.sum(data_weights_cut)) > 0 else 0
    # f1_score = 2 * (efficiency * purity) / (efficiency + purity) if (efficiency + purity) > 0 else 0
    f1_score, efficiency, purity, total_mc_events_after_cuts, total_data_events_after_cuts, significance = log_likelihood_and_values(
        [cut_thresholds["vertex_z_min"], cut_thresholds["vertex_y_max"], cut_thresholds["num_pf_particles_min"], cut_thresholds["direction_z_max"], cut_thresholds["energy_first10cm_min"], cut_thresholds["energy_fifth10cm_min"], cut_thresholds["energy_fifteenth10cm_min"], cut_thresholds["roi_z_size_min"], cut_thresholds["roi_z_vertex_distance_max"], cut_thresholds["tail_length_density_max"]],
        mc_events, mc_weights, data_events, data_weights
    )
    f.write("\n")
    f.write("Default cuts:\n")
    f.write("Number of MC events after cuts normalized: " + str(np.sum(mc_weights_cut)) + "\n")
    f.write("Number of off spill data events after cuts normalized: " + str(np.sum(data_weights_cut)) + "\n")
    f.write("MC/off spill data ratio after cuts: " + str(np.sum(mc_weights_cut)/np.sum(data_weights_cut)) + "\n")
    f.write(f"Efficiency: {efficiency}\n")
    f.write(f"Purity: {purity}\n")
    f.write(f"F1 Score: {f1_score}\n")
    f.write(f"Significance: {significance}\n")
    f.write("Number of MC events after cuts absolute: " + str(len(mc_weights_cut)) + "\n")
    f.write("Number of off spill data events after cuts absolute: " + str(len(data_weights_cut)) + "\n")

    # extra plot, show the distribution of the vertex Z variable for the MC and data events before and after cuts, using the optimized cuts
    plt.figure(figsize=(10, 6))
    plt.hist(mc_events[:, aggregate_dict["vertexZ"]], bins=50, alpha=0.5, label="MC before cuts", color='blue', weights=mc_weights/np.sum(mc_weights))
    plt.hist(data_events[:, aggregate_dict["vertexZ"]], bins=50, alpha=0.5, label="Data before cuts", color='orange', weights=data_weights/np.sum(data_weights))
    plt.hist(mc_events_cut[:, aggregate_dict["vertexZ"]], bins=50, alpha=0.5, label="MC after optimized cuts", color='green', weights=mc_weights_cut/np.sum(mc_weights_cut))
    plt.hist(data_events_cut[:, aggregate_dict["vertexZ"]], bins=50, alpha=0.5, label="Data after optimized cuts", color='red', weights=data_weights_cut/np.sum(data_weights_cut))
    plt.xlabel("Vertex Z")
    plt.ylabel("Number of events")
    plt.title("Distribution of vertex Z before and after cuts")
    plt.legend()
    plt.savefig(os.path.join(output_folder_opt, "vertexZ_distribution_before_after_cuts.png"))
    plt.close() 

# for each parameter, fix the others to the best value and vary it to see how the score changes, and plot the results
for i, key in enumerate(threshold_boundaries.keys()):
    param_values = np.linspace(max(threshold_boundaries[key][0], best_params[i]*0.6), min(threshold_boundaries[key][1], best_params[i]*1.4), 20)
    scores = []
    efficiencies = []
    purities = []
    significances = []
    total_mc = []
    total_data = []
    for value in param_values:
        params = best_params.copy()
        params[i] = value
        f1_score, efficiency, purity, total_mc_events, total_data_events, significance = log_likelihood_and_values(params, mc_events, mc_weights, data_events, data_weights)
        scores.append(f1_score)
        efficiencies.append(efficiency)
        purities.append(purity)
        significances.append(significance)
        total_mc.append(total_mc_events)
        total_data.append(total_data_events)

    fig, axes = plt.subplots(1, 4, figsize=(15, 4))
    axes[0].plot(param_values, scores, color='blue', linewidth=2)
    axes[0].axvline(x=best_params[i], color='r', linestyle='--', label='Optimized value')
    axes[0].set_xlabel(key)
    axes[0].set_ylabel("Score")
    axes[0].set_title(f"Score vs {key}")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].plot(param_values, efficiencies, color='green', linewidth=2)
    axes[1].axvline(x=best_params[i], color='r', linestyle='--', label='Optimized value')
    axes[1].set_xlabel(key)
    axes[1].set_ylabel("Efficiency")
    axes[1].set_title(f"Efficiency vs {key}")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    axes[2].plot(param_values, purities, color='orange', linewidth=2)
    axes[2].axvline(x=best_params[i], color='r', linestyle='--', label='Optimized value')
    axes[2].set_xlabel(key)
    axes[2].set_ylabel("Purity")
    axes[2].set_title(f"Purity vs {key}")
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    axes[3].plot(param_values, significances, color='magenta', linewidth=2)
    axes[3].axvline(x=best_params[i], color='r', linestyle='--', label='Optimized value')
    axes[3].set_xlabel(key)
    axes[3].set_ylabel("Significance")
    axes[3].set_title(f"Significance vs {key}")
    axes[3].legend()
    axes[3].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder_opt, f"{key}_score_efficiency_purity_significance.png"))
    plt.close()
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].plot(param_values, total_mc, color='purple', linewidth=2)
    axes[0].axvline(x=best_params[i], color='r', linestyle='--', label='Optimized value')
    axes[0].set_xlabel(key)
    axes[0].set_ylabel("Number of events")
    axes[0].set_title(f"Total MC events after cuts vs {key}")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].plot(param_values, total_data, color='cyan', linewidth=2)
    axes[1].axvline(x=best_params[i], color='r', linestyle='--', label='Optimized value')
    axes[1].set_xlabel(key)
    axes[1].set_ylabel("Number of events")
    axes[1].set_title(f"Total off spill data events after cuts vs {key}")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder_opt, f"{key}_total_events.png"))
    plt.close()


# do the same for the default cuts, varying each parameter around the default value

for i, key in enumerate(threshold_boundaries.keys()):
    param_values = np.linspace(max(threshold_boundaries[key][0], default_cut_thresholds[key]*0.6), min(threshold_boundaries[key][1], default_cut_thresholds[key]*1.4), 20)
    scores = []
    efficiencies = []
    purities = []
    significances = []
    total_mc = []
    total_data = []
    for value in param_values:
        cut_thresholds = default_cut_thresholds.copy()
        cut_thresholds[key] = value
        mc_events_cut, mc_weights_cut, mc_true_events_cut, mc_labels_cut, mc_filenames_cut, data_events_cut, data_weights_cut, data_labels_cut, data_filenames_cut = apply_all_cuts(
            mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames,
            data_events, data_weights, data_labels, data_filenames,
            cut_thresholds=cut_thresholds
        )
        efficiency = np.sum(mc_weights_cut) / np.sum(mc_weights) if np.sum(mc_weights) > 0 else 0
        purity = np.sum(mc_weights_cut) / (np.sum(mc_weights_cut) + np.sum(data_weights_cut)) if (np.sum(mc_weights_cut) + np.sum(data_weights_cut)) > 0 else 0
        f1_score = 2 * (efficiency * purity) / (efficiency + purity) if (efficiency + purity) > 0 else 0
        T_on = 7.65 + 8.88
        T_off = 44.688
        n_on = (np.sum(mc_weights_cut) + np.sum(data_weights_cut)) * T_on
        n_off = np.sum(data_weights_cut) * T_off
        mu = n_off * (T_on / T_off) if n_off > 0 else 1e-10
        sig_val = 2 * (n_on * np.log(n_on / mu) + mu - n_on) if n_on > 0 and mu > 0 else 0
        if not np.isfinite(sig_val):
            sig_val = 0
        scores.append(f1_score)
        efficiencies.append(efficiency)
        purities.append(purity)
        significances.append(sig_val)
        total_mc.append(np.sum(mc_weights_cut))
        total_data.append(np.sum(data_weights_cut))

    fig, axes = plt.subplots(1, 4, figsize=(15, 4))
    axes[0].plot(param_values, scores, color='blue', linewidth=2)
    axes[0].axvline(x=default_cut_thresholds[key], color='r', linestyle='--', label='Default value')
    axes[0].set_xlabel(key)
    axes[0].set_ylabel("Score")
    axes[0].set_title(f"Score vs {key}")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].plot(param_values, efficiencies, color='green', linewidth=2)
    axes[1].axvline(x=default_cut_thresholds[key], color='r', linestyle='--', label='Default value')
    axes[1].set_xlabel(key)
    axes[1].set_ylabel("Efficiency")
    axes[1].set_title(f"Efficiency vs {key}")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    axes[2].plot(param_values, purities, color='orange', linewidth=2)
    axes[2].axvline(x=default_cut_thresholds[key], color='r', linestyle='--', label='Default value')
    axes[2].set_xlabel(key)
    axes[2].set_ylabel("Purity")
    axes[2].set_title(f"Purity vs {key}")
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    axes[3].plot(param_values, significances, color='magenta', linewidth=2)
    axes[3].axvline(x=default_cut_thresholds[key], color='r', linestyle='--', label='Default value')
    axes[3].set_xlabel(key)
    axes[3].set_ylabel("Significance")
    axes[3].set_title(f"Significance vs {key}")
    axes[3].legend()
    axes[3].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_folder_def, f"{key}_score_efficiency_purity_significance.png"))
    plt.close()
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].plot(param_values, total_mc, color='purple', linewidth=2)
    axes[0].axvline(x=default_cut_thresholds[key], color='r', linestyle='--', label='Default value')
    axes[0].set_xlabel(key)
    axes[0].set_ylabel("Number of events")
    axes[0].set_title(f"Total MC events after cuts vs {key}")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].plot(param_values, total_data, color='cyan', linewidth=2)
    axes[1].axvline(x=default_cut_thresholds[key], color='r', linestyle='--', label='Default value')
    axes[1].set_xlabel(key)
    axes[1].set_ylabel("Number of events")
    axes[1].set_title(f"Total off spill data events after cuts vs {key}")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder_def, f"{key}_total_events.png"))
    plt.close()



