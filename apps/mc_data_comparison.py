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

data_events, data_true_events, data_full_events, data_weights, data_labels, data_filenames = load_from_folder(  parameters["folders"]["mother_folder_data"], 
                                                                                                    get_truth=parameters["data"]["GET_TRUTH"], 
                                                                                                    get_full_array=True, 
                                                                                                    use_combined=parameters["data"]["USE_COMBINED"], 
                                                                                                    weights_mode=parameters["data"]["WEIGHT_MODE"], 
                                                                                                    parameters=parameters["data"], 
                                                                                                    spill_status=parameters["data"]["SPILL_STATUS"])
  




print(f"Loaded {len(data_events)} data events")
print(f"Loaded {len(mc_events)} MC events")

# ------------------ PLOTS -----------------------
run_all_triple_hist(
    mc_events, mc_weights, data_events, data_weights, 
    output_folder_base=output_folder_base, 
    label_sig=parameters["mc"]["label"], 
    label_bkg=parameters["data"]["label"], 
    nbins=parameters["analysis"]["nbins"], 
    spill_status = f"{parameters['mc']['SPILL_STATUS']}_{parameters['data']['SPILL_STATUS']}",
    apply_cuts=parameters["analysis"]["apply_cuts"],
    cuts=cuts, 
    total_time_with_correct_tps_on=parameters["mc"]["run_parameters"]["spill_on"]["total_time"]
    )
    
# -----------------------------------------
# ------------------ Save hist value to file -----------------------
mc_events_orig = mc_events.copy()
mc_weights_orig = mc_weights.copy()
mc_true_events_orig = mc_true_events.copy()
mc_filenames_orig = mc_filenames.copy()
mc_labels_orig = mc_labels.copy()
data_events_orig = data_events.copy()
data_weights_orig = data_weights.copy()
data_labels_orig = data_labels.copy()
data_filenames_orig = data_filenames.copy()

mc_remaining_events = []
signal_remaining_events = []

mc_ev, mc_w, mc_true, mc_lab, mc_fnames, data_ev, data_w, data_lab, data_fnames,  = apply_all_cuts(
    mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
    data_events_orig, data_weights_orig, data_labels_orig, data_filenames_orig,
    cuts=cuts, skip_cut="max_directionZ"
)

hist, bin = np.histogram(np.maximum(
            mc_ev[:, aggregate_dict["directionZ"]],
            mc_ev[:, aggregate_dict["directionZ2"]]
        ), bins=parameters["analysis"]["nbins"], range=(-1,1), weights=mc_w)
np.savetxt(os.path.join(output_folder_base, "mc_max_directionZ_hist.txt"), hist)

hist, bin = np.histogram(np.maximum(
            mc_ev[:, aggregate_dict["directionZ"]],
            mc_ev[:, aggregate_dict["directionZ2"]]
        ), bins=parameters["analysis"]["nbins"], range=(-1,1))
np.savetxt(os.path.join(output_folder_base, "mc_max_directionZ_hist_absolute.txt"), hist)

# do the same for data
hist, bin = np.histogram(np.maximum(
            data_ev[:, aggregate_dict["directionZ"]],
            data_ev[:, aggregate_dict["directionZ2"]]
        ), bins=parameters["analysis"]["nbins"], range=(-1,1), weights=data_w)
np.savetxt(os.path.join(output_folder_base, "data_max_directionZ_hist.txt"), hist)
hist, bin = np.histogram(np.maximum(
            data_ev[:, aggregate_dict["directionZ"]],
            data_ev[:, aggregate_dict["directionZ2"]]
        ), bins=parameters["analysis"]["nbins"], range=(-1,1))
np.savetxt(os.path.join(output_folder_base, "data_max_directionZ_hist_absolute.txt"), hist)



# ----------------------------------------------
# Apply cuts in order, save the results in a latex table
create_table_cuts(  mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
                    data_events_orig, data_weights_orig, data_labels_orig, data_filenames_orig,
                    cuts_names=cuts_names, skip_cut=None, output_folder=output_folder_base)







