import numpy as np
import os
import sys
import argparse
import json
from collections import defaultdict

# import custom libs
sys.path.append("../python")
from cuts import *
from libs import *
from name_index_association import *


def strip_job_folder(filename, job_folders):
    # eventID-only equality is not enough to spot a job-submission overlap:
    # eventID is a local counter that legitimately repeats across the many
    # leaf (subrun-like) directories inside a single job folder. The real
    # signature of two overlapping job submissions is the same eventID
    # coming from the same relative sub-path (e.g. "1/001/combined_45228.0")
    # under two different job-ID folders, so duplicates are keyed on that
    # relative path rather than on the raw eventID.
    #
    # get_files() resolves USE_COMBINED files on the /pnfs/.../persistent/
    # mirror when available (falling back to /pnfs/.../scratch/ otherwise),
    # so the recorded filename doesn't always share the exact prefix given
    # in mother_folder_signal. Match on the job-ID directory name (the last
    # path component of each job folder) instead, which is present in the
    # filename either way.
    job_ids = {job_folder.rstrip("/").split("/")[-1] for job_folder in job_folders}
    parts = filename.split("/")
    for i, part in enumerate(parts):
        if part in job_ids:
            return part, "/".join(parts[i + 1:])
    return None, filename


def find_duplicate_files(events, filenames, job_folders):
    groups = defaultdict(list)
    for ev, fname in zip(events, filenames):
        job_folder, rel_path = strip_job_folder(fname, job_folders)
        eid = int(ev[aggregate_dict["eventID"]])
        groups[rel_path].append((eid, job_folder, fname))
    duplicates = {rel_path: entries for rel_path, entries in groups.items() if len(entries) > 1}
    return duplicates


def write_section(f, title, events, filenames, job_folders):
    n_total = len(events)
    n_unique = len(set(int(ev[aggregate_dict["eventID"]]) for ev in events))
    duplicates = find_duplicate_files(events, filenames, job_folders)

    f.write(f"--- {title} ---\n")
    f.write(f"Total events: {n_total}\n")
    f.write(f"Unique eventIDs: {n_unique}\n")
    f.write(f"Overlapping eventIDs (same eventID, same file, found under >1 job folder): {len(duplicates)}\n")
    for rel_path in sorted(duplicates):
        eid = duplicates[rel_path][0][0]
        f.write(f"  eventID {eid} (file: {rel_path}):\n")
        for _, job_folder, fname in duplicates[rel_path]:
            f.write(f"    {fname}\n")
    f.write("\n")
    return duplicates


parser = argparse.ArgumentParser(description='Check for duplicate/overlapping events within a run, before and after cuts.')
parser.add_argument('-j', type=str, help='the json file with the parameters')
args = parser.parse_args()

with open(args.j) as f:
    parameters = json.load(f)

output_folder_base = parameters["folders"]["output_folder_base"]
if output_folder_base[-1] == "/":
    output_folder_base = output_folder_base[:-1]

if not os.path.exists(output_folder_base):
    os.makedirs(output_folder_base)

run_parameters = parameters["run"]
run_label = run_parameters.get("label", "run")

sig_events, sig_true_events, sig_full_events, sig_weights, sig_labels, sig_filenames = load_from_folder(
    parameters["folders"]["mother_folder_signal"],
    get_truth=run_parameters["GET_TRUTH"],
    get_full_array=False,
    use_combined=run_parameters["USE_COMBINED"],
    weights_mode=run_parameters["WEIGHT_MODE"],
    parameters=run_parameters,
    spill_status=run_parameters["SPILL_STATUS"],
    timing=parameters.get("timing")
)

cut_events, cut_weights, cut_true_events, cut_labels, cut_filenames = apply_all_cuts_single(
    sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
    cuts=cuts_single, skip_cut=None
)

job_folders = parameters["folders"]["mother_folder_signal"]

output_file = os.path.join(output_folder_base, "overlap_summary.txt")
with open(output_file, "w") as f:
    f.write(f"Overlap check for {run_label}\n")
    f.write("=" * (len(run_label) + 18) + "\n\n")
    f.write("Input job folders:\n")
    for folder in job_folders:
        f.write(f"  {folder}\n")
    f.write("\n")

    duplicates_before = write_section(f, "BEFORE CUTS", sig_events, sig_filenames, job_folders)
    duplicates_after = write_section(f, "AFTER CUTS", cut_events, cut_filenames, job_folders)

print(f"Before cuts: {len(sig_events)} events, {len(duplicates_before)} overlapping eventIDs")
print(f"After cuts: {len(cut_events)} events, {len(duplicates_after)} overlapping eventIDs")
print(f"Summary written to {output_file}")
