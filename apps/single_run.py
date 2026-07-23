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
                                                                                                    spill_status=parameters["run"]["SPILL_STATUS"],
                                                                                                    timing=parameters["timing"])

if parameters["run"]["TP_RATE"] == "mc":
    print('------------')
    print(f"Before TA cut we have: {np.sum(run_weights)} per hour")
    print("TA cut:")
    # Apply the trigger activity cut
    print(len(run_events), len(run_weights), len(run_true_events), len(run_labels))
    print(f"run events: {np.sum(run_weights)}")

    index_run = np.where(run_true_events[:, true_dict["triggerActivityFlag"]] == 1)[0]
    print(f"Activity flag: {run_true_events[:, true_dict['triggerActivityFlag']].sum()} / {len(run_true_events)}")
    print(f"Index run: {len(index_run)}")

    run_events = run_events[index_run]
    run_weights = run_weights[index_run]
    run_true_events = run_true_events[index_run]
    run_labels = run_labels[index_run]
    run_filenames = [run_filenames[i] for i in index_run]
    print(f"After TA cut we have: {np.sum(run_weights)} per hour")


print(f"Events: {len(run_events)}, Weights: {len(run_weights)}, True events: {len(run_true_events)}, Labels: {len(run_labels)}, Filenames: {len(run_filenames)}")

# spillStatusFlag
print((f"Spill on files: {len(run_events[:, aggregate_dict['spillStatusFlag'] == 1])}"))
print((f"Spill off files: {len(run_events[:, aggregate_dict['spillStatusFlag'] == 0])}"))

# ------------------ PLOTS -----------------------
run_tp_rate = parameters["run"]["TP_RATE"]
if run_tp_rate == "all":
    total_time_with_correct_tps_on = parameters["timing"]["total_time_all_tp_rate"]
elif run_tp_rate == "high":
    total_time_with_correct_tps_on = parameters["timing"]["total_time_high_tp_rate"]
elif run_tp_rate == "low":
    total_time_with_correct_tps_on = parameters["timing"]["total_time_low_tp_rate"]
elif run_tp_rate == "mc":
    total_time_with_correct_tps_on = parameters["run"]["run_parameters"]["spill_on"]["total_time"]
else:
    print(f"Unknown TP_RATE: {run_tp_rate}. Exiting.")
    sys.exit(1)

run_all_triple_hist_single(
    run_events, run_weights,
    output_folder_base=output_folder_base,
    label_sig=parameters["run"]["label"],
    nbins=parameters["analysis"]["nbins"],
    spill_status = f"{parameters['run']['SPILL_STATUS']}",
    apply_cuts=parameters["analysis"]["apply_cuts"],
    cuts=cuts_single,
    total_time_with_correct_tps_on=total_time_with_correct_tps_on
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

# ----------------------------------------------
# MC vertex resolution plots (only when truth info is available)
if len(run_true_events_orig) > 0:
    run_ev_cut, run_w_cut, run_true_cut, run_lab_cut, run_fnames_cut = apply_all_cuts_single(
        run_events_orig, run_weights_orig, run_true_events_orig, run_labels_orig, run_filenames_orig,
        cuts=cuts_single, skip_cut=None
    )

    if len(run_ev_cut) > 0:
        label = parameters["run"]["label"]
        nbins = parameters["analysis"]["nbins"]

        dx = run_ev_cut[:, aggregate_dict["vertexX"]] - run_true_cut[:, true_dict["vertexX"]]
        dy = run_ev_cut[:, aggregate_dict["vertexY"]] - run_true_cut[:, true_dict["vertexY"]]
        dz = run_ev_cut[:, aggregate_dict["vertexZ"]] - run_true_cut[:, true_dict["vertexZ"]]
        dist3d = np.sqrt(dx**2 + dy**2 + dz**2)

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        axes[0, 0].hist(dx, bins=nbins, range=(-15, 15), weights=run_w_cut/np.sum(run_w_cut), color="orange", alpha=0.7)
        axes[0, 0].set_xlabel("Reco - True Vertex X [cm]")
        axes[0, 0].set_ylabel("Counts")
        axes[0, 0].set_title("Vertex X Residual")
        axes[0, 0].grid(alpha=0.5)

        axes[0, 1].hist(dy, bins=nbins, range=(-15, 15), weights=run_w_cut/np.sum(run_w_cut), color="orange", alpha=0.7)
        axes[0, 1].set_xlabel("Reco - True Vertex Y [cm]")
        axes[0, 1].set_ylabel("Counts")
        axes[0, 1].set_title("Vertex Y Residual")
        axes[0, 1].grid(alpha=0.5)

        axes[1, 0].hist(dz, bins=nbins, range=(-15, 15), weights=run_w_cut/np.sum(run_w_cut), color="orange", alpha=0.7)
        axes[1, 0].set_xlabel("Reco - True Vertex Z [cm]")
        axes[1, 0].set_ylabel("Counts")
        axes[1, 0].set_title("Vertex Z Residual")
        axes[1, 0].grid(alpha=0.5)

        axes[1, 1].hist(dist3d, bins=nbins, range=(0, 25), weights=run_w_cut/np.sum(run_w_cut), color="orange", alpha=0.7)
        axes[1, 1].set_xlabel("3D Distance [cm]")
        axes[1, 1].set_ylabel("Counts")
        axes[1, 1].set_title("3D Reco-True Vertex Distance")
        axes[1, 1].grid(alpha=0.5)

        plt.suptitle(f"{label}: Vertex Resolution (all cuts applied)", fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder_base, "mc_vertex_resolution.png"))
        plt.close()
        print(f"Saved mc_vertex_resolution.png with {len(run_ev_cut)} events passing all cuts.")
