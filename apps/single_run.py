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


events_target = parameters["events_target"] if "events_target" in parameters else None

output_folder_base = parameters["folders"]["output_folder_base"]
if output_folder_base[-1] == "/":
    output_folder_base = output_folder_base[:-1]
output_folder_base = output_folder_base + "_TP_"+parameters["analysis"]["TP_RATE"]+"_cuts_"+str(parameters["analysis"]["apply_cuts"])+"_spill_"+parameters["run"]["SPILL_STATUS"]

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

run_events, run_true_events, run_full_events, run_weights, run_labels, run_filenames = load_from_folder(  parameters["folders"]["mother_folder_run"], 
                                                                                                    get_truth=parameters["run"]["GET_TRUTH"], 
                                                                                                    get_full_array=False, 
                                                                                                    use_combined=parameters["run"]["USE_COMBINED"], 
                                                                                                    weights_mode=parameters["run"]["WEIGHT_MODE"], 
                                                                                                    parameters=parameters["run"], 
                                                                                                    spill_status=parameters["run"]["SPILL_STATUS"])

if parameters["run"]["TP_RATE"] == "mc":
    print('------------')
    print("TA cut:")
    # Apply the trigger activity cut
    print(len(run_events), len(run_weights), len(run_true_events), len(run_labels))
    print(f"run events: {np.sum(run_weights)}")

    index_run = np.where(run_true_events[:, true_dict["triggerActivityFlag"]] == 1)[0]
    run_events = run_events[index_run]
    run_weights = run_weights[index_run]
    run_true_events = run_true_events[index_run]
    run_labels = run_labels[index_run]
    run_filenames = [run_filenames[i] for i in index_run]


print(f"Events: {len(run_events)}, Weights: {len(run_weights)}, True events: {len(run_true_events)}, Labels: {len(run_labels)}, Filenames: {len(run_filenames)}")

# ------------------ PLOTS -----------------------
run_all_triple_hist_single(
    run_events, run_weights,
    output_folder_base=output_folder_base, 
    label_sig=parameters["run"]["label"], 
    nbins=parameters["analysis"]["nbins"], 
    spill_status = f"{parameters['run']['SPILL_STATUS']}",
    apply_cuts=parameters["analysis"]["apply_cuts"],
    cuts=cuts_single, 
    total_time_with_correct_tps_on=parameters["run"]["run_parameters"]["spill_on"]["total_time"]
    )


# -----------------------------------------
# ------------------ Save hist value to file -----------------------
run_events_orig = run_events.copy()
run_weights_orig = run_weights.copy()
run_true_events_orig = run_true_events.copy()
run_filenames_orig = run_filenames.copy()
run_labels_orig = run_labels.copy()

run_remaining_events = []

run_ev, run_w, run_true, run_lab, run_fnames = apply_all_cuts_single(
    run_events_orig, run_weights_orig, run_true_events_orig, run_labels_orig, run_filenames_orig,
    cuts=cuts_single, skip_cut="max_directionZ"
)

hist, bin = np.histogram(np.maximum(
            run_ev[:, aggregate_dict["directionZ"]],
            run_ev[:, aggregate_dict["directionZ2"]]
        ), bins=parameters["analysis"]["nbins"], range=(-1,1), weights=run_w)
np.savetxt(os.path.join(output_folder_base, "run_max_directionZ_hist.txt"), hist)

hist, bin = np.histogram(np.maximum(
            run_ev[:, aggregate_dict["directionZ"]],
            run_ev[:, aggregate_dict["directionZ2"]]
        ), bins=parameters["analysis"]["nbins"], range=(-1,1))
np.savetxt(os.path.join(output_folder_base, "run_max_directionZ_hist_absolute.txt"), hist)



# ----------------------------------------------
# Apply cuts in order, save the results in a latex table
create_table_cuts_single(  run_events_orig, run_weights_orig, run_true_events_orig, run_labels_orig, run_filenames_orig,
                    cuts_names=cuts_names, skip_cut=None, output_folder=output_folder_base, events_target = events_target
)

# ----------------------------------------------
# Apply cuts in order, save the results in a latex table
dump_information_events_single(  run_events_orig, run_weights_orig, run_true_events_orig, run_labels_orig, run_filenames_orig,
                    cuts_names=cuts_names, skip_cut=None, output_folder=output_folder_base, events_target = events_target
)


