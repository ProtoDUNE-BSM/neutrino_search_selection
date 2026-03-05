import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import uproot
import argparse
import json
import csv
from scipy import stats


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
output_folder_base = output_folder_base + "_"+parameters["analysis"]["TP_RATE"]+"_cuts_"+str(parameters["analysis"]["apply_cuts"])

output_folder_normalized = output_folder_base+"/output_test_normalized/"
output_folder_absolute = output_folder_base+"/output_test_absolute/"
output_folder_1hour = output_folder_base+"/output_test_1hour/"
output_folder_full_time = output_folder_base+"/output_test_full_time/"

if not os.path.exists(output_folder_base):
    os.makedirs(output_folder_base)
if not os.path.exists(output_folder_normalized):
    os.makedirs(output_folder_normalized)
if not os.path.exists(output_folder_absolute):
    os.makedirs(output_folder_absolute)
if not os.path.exists(output_folder_1hour):
    os.makedirs(output_folder_1hour)
if not os.path.exists(output_folder_full_time):
    os.makedirs(output_folder_full_time)


mc_events, mc_true_events, mc_full_events, mc_weights, mc_labels, mc_filenames = load_from_folder(  parameters["folders"]["mother_folder_mc"], 
                                                                                                    get_truth=parameters["mc"]["GET_TRUTH"], 
                                                                                                    get_full_array=True, 
                                                                                                    use_combined=parameters["mc"]["USE_COMBINED"], 
                                                                                                    weights_mode=parameters["mc"]["WEIGHT_MODE"], 
                                                                                                    parameters=parameters["mc"], 
                                                                                                    spill_status=parameters["mc"]["SPILL_STATUS"])

print('------------')
print("TA cut:")
# Apply the trigger activity cut
print(len(mc_events), len(mc_weights), len(mc_true_events), len(mc_labels))
print(f"MC events: {np.sum(mc_weights)}")

index_mc = np.where(mc_true_events[:, true_dict["triggerActivityFlag"]] == 1)[0]
mc_events = mc_events[index_mc]
mc_weights = mc_weights[index_mc]
mc_true_events = mc_true_events[index_mc]
mc_labels = mc_labels[index_mc]
mc_filenames = [mc_filenames[i] for i in index_mc]

print(f"MC events: {np.sum(mc_weights)}")

sig_events, sig_true_events, sig_full_events, sig_weights, sig_labels, sig_filenames = load_from_folder(  parameters["folders"]["mother_folder_signal"], 
                                                                                                    get_truth=parameters["signal"]["GET_TRUTH"], 
                                                                                                    get_full_array=True, 
                                                                                                    use_combined=parameters["signal"]["USE_COMBINED"], 
                                                                                                    weights_mode=parameters["signal"]["WEIGHT_MODE"], 
                                                                                                    parameters=parameters["signal"], 
                                                                                                    spill_status=parameters["signal"]["SPILL_STATUS"])
                                                                            
bkg_events, bkg_true_events, bkg_full_events, bkg_weights, bkg_labels, bkg_filenames = load_from_folder(    parameters["folders"]["mother_folder_bkg"],
                                                                                                            get_truth=parameters["bkg"]["GET_TRUTH"],
                                                                                                            get_full_array=True, 
                                                                                                            use_combined=parameters["bkg"]["USE_COMBINED"], 
                                                                                                            weights_mode=parameters["bkg"]["WEIGHT_MODE"], 
                                                                                                            parameters=parameters["bkg"], 
                                                                                                            spill_status=parameters["bkg"]["SPILL_STATUS"])



print(f"Loaded {len(sig_events)} signal events")
print(f"Loaded {len(bkg_events)} background events")
print(f"Loaded {len(mc_events)} MC events")

# ------------------ PLOTS -----------------------

run_all_triple_hist_stack(
    sig_events, sig_weights, bkg_events, bkg_weights, mc_events, mc_weights,
    output_folder_base= output_folder_base,
    label_sig = parameters["signal"]["label"],
    label_bkg = parameters["bkg"]["label"],
    nbins = parameters["analysis"]["nbins"],
    spill_status = f"{parameters['signal']['SPILL_STATUS']}_{parameters['bkg']['SPILL_STATUS']}",
    apply_cuts=parameters["analysis"]["apply_cuts"],
    cuts=cuts_with_mc, 
    total_time_with_correct_tps_on=parameters["signal"]["run_parameters"]["spill_on"]["total_time"]
    )

# -----------------------------------------

# ------------------ Save hist value to file -----------------------
sig_events_orig = sig_events.copy()
sig_weights_orig = sig_weights.copy()
sig_true_events_orig = sig_true_events.copy()
sig_filenames_orig = sig_filenames.copy()
sig_labels_orig = sig_labels.copy()
bkg_events_orig = bkg_events.copy()
bkg_weights_orig = bkg_weights.copy()
bkg_labels_orig = bkg_labels.copy()
bkg_filenames_orig = bkg_filenames.copy()
mc_events_orig = mc_events.copy()
mc_weights_orig = mc_weights.copy()
mc_true_events_orig = mc_true_events.copy()
mc_filenames_orig = mc_filenames.copy()
mc_labels_orig = mc_labels.copy()


signal_remaining_events = []
background_remaining_events = []
mc_remaining_events = []

s_ev_mc, s_w_mc, s_true_mc, s_lab_mc, s_fnames_mc, b_ev_mc, b_w_mc, b_lab_mc, b_fnames_mc, mc_ev_mc, mc_w_mc, mc_true_mc, mc_lab_mc, mc_fnames_mc = apply_all_cuts_MC(
    sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
    bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
    mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
    cuts=cuts_with_mc, skip_cut="max_directionZ"
)
this_hist_n_bins = 40
hist, bin = np.histogram(np.maximum(
            b_ev_mc[:, aggregate_dict["directionZ"]],
            b_ev_mc[:, aggregate_dict["directionZ2"]]
        ), bins=this_hist_n_bins, range=(-1,1), weights=b_w_mc)

# save the only the hist to file
np.savetxt(os.path.join(output_folder_base, "bkg_max_directionZ_hist.txt"), hist)
hist, bin = np.histogram(np.maximum(
            b_ev_mc[:, aggregate_dict["directionZ"]],
            b_ev_mc[:, aggregate_dict["directionZ2"]]
        ), bins=this_hist_n_bins, range=(-1,1))

# save the only the hist to file
np.savetxt(os.path.join(output_folder_base, "bkg_max_directionZ_hist_absolute.txt"), hist)

hist, bin = np.histogram(np.maximum(
            s_ev_mc[:, aggregate_dict["directionZ"]],
            s_ev_mc[:, aggregate_dict["directionZ2"]]
        ), bins=this_hist_n_bins, range=(-1,1), weights=s_w_mc)
np.savetxt(os.path.join(output_folder_base, "sig_max_directionZ_hist.txt"), hist)

hist, bin = np.histogram(np.maximum(
            s_ev_mc[:, aggregate_dict["directionZ"]],
            s_ev_mc[:, aggregate_dict["directionZ2"]]
        ), bins=this_hist_n_bins, range=(-1,1))
np.savetxt(os.path.join(output_folder_base, "sig_max_directionZ_hist_absolute.txt"), hist)

# do the same for MC
hist, bin = np.histogram(np.maximum(
            mc_ev_mc[:, aggregate_dict["directionZ"]],
            mc_ev_mc[:, aggregate_dict["directionZ2"]]
        ), bins=this_hist_n_bins, range=(-1,1), weights=mc_w_mc)
np.savetxt(os.path.join(output_folder_base, "mc_max_directionZ_hist.txt"), hist)
hist, bin = np.histogram(np.maximum(
            mc_ev_mc[:, aggregate_dict["directionZ"]],
            mc_ev_mc[:, aggregate_dict["directionZ2"]]
        ), bins=this_hist_n_bins, range=(-1,1))
np.savetxt(os.path.join(output_folder_base, "mc_max_directionZ_hist_absolute.txt"), hist)



# ----------------------------------------------
# Apply cuts in order, save the results in a latex table
create_table_cuts(sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
                    bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
                    mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
                    cuts_names=cuts_names, skip_cut=None, output_folder=output_folder_base)


# get more information about the events passing all cuts
s_ev, s_w, s_true, s_lab, s_fnames, b_ev, b_w, b_lab, b_fnames = apply_all_cuts(
    sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
    bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
    cuts=cuts
)

# save to file the list of events passing all cuts
csv_path = os.path.join(output_folder_base, "events_passing_all_cuts.csv")
with open(csv_path, "w", newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["EventID", "Filename", "IsSignal"])
    for i in range(len(s_ev)):
        writer.writerow([
            int(s_ev[i, aggregate_dict["eventID"]]),
            s_fnames[i],
            "Signal"
        ])
    for i in range(len(b_ev)):
        writer.writerow([
            int(b_ev[i, aggregate_dict["eventID"]]),
            b_fnames[i],
            "Background"
        ])
print(f"List of events passing all cuts saved to {csv_path}\n")



