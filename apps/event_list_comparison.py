import os
import re
import csv
import json
import argparse
from collections import defaultdict


def load_events(cfg):
    # Run number isn't always literally in the Filename: some job outputs only
    # carry a job ID, which job_id_to_run translates to the real run number.
    run_regex = re.compile(cfg["run_regex"])
    job_id_to_run = {int(k): int(v) for k, v in cfg.get("job_id_to_run", {}).items()}

    events = {}
    with open(cfg["csv"]) as f:
        reader = csv.DictReader(f)
        for row in reader:
            m = run_regex.search(row["Filename"])
            if not m:
                raise ValueError(f"run_regex did not match filename: {row['Filename']}")
            run = int(m.group(1))
            if job_id_to_run:
                run = job_id_to_run[run]
            events[(run, int(row["EventID"]))] = row["IsSignal"]
    return events


def format_by_run(keys):
    grouped = defaultdict(list)
    for run, eid in keys:
        grouped[run].append(eid)
    lines = []
    for run in sorted(grouped):
        eids = sorted(grouped[run])
        lines.append(f"Run {run} ({len(eids)}): {', '.join(str(e) for e in eids)}")
    return lines


def write_group_file(path, is_signal, common, lost, new_only, old_events, new_events, old_label, new_label):
    common_keys = [k for k in common if new_events[k] == is_signal]
    lost_keys = [k for k in lost if old_events[k] == is_signal]
    new_keys = [k for k in new_only if new_events[k] == is_signal]

    with open(path, "w") as f:
        f.write(f"{is_signal} comparison: {old_label} vs {new_label}\n\n")
        f.write(f"Number of common {is_signal.lower()}s: {len(common_keys)}\n")
        f.write(f"Number of unique old {is_signal.lower()}s: {len(lost_keys)}\n")
        f.write(f"Number of unique new {is_signal.lower()}s: {len(new_keys)}\n\n")

        sections = [
            (f"Common {is_signal.lower()} event IDs", common_keys),
            (f"Unique old {is_signal.lower()} event IDs", lost_keys),
            (f"Unique new {is_signal.lower()} event IDs", new_keys),
        ]
        for title, keys in sections:
            f.write(f"--- {title} ---\n")
            for line in format_by_run(keys):
                f.write(line + "\n")
            f.write("\n")


parser = argparse.ArgumentParser(description="Diff two events_passing_all_cuts.csv files, matching on (run, eventID).")
parser.add_argument('-j', type=str, help='the json file with the parameters')
args = parser.parse_args()

with open(args.j) as f:
    parameters = json.load(f)

old_cfg = parameters["old"]
new_cfg = parameters["new"]

output_folder_base = parameters["output_folder_base"]
if output_folder_base[-1] == "/":
    output_folder_base = output_folder_base[:-1]
if not os.path.exists(output_folder_base):
    os.makedirs(output_folder_base)

old_events = load_events(old_cfg)
new_events = load_events(new_cfg)

old_keys = set(old_events)
new_keys = set(new_events)

common = old_keys & new_keys
lost = old_keys - new_keys       # in old, not in new
new_only = new_keys - old_keys   # in new, not in old

signal_path = os.path.join(output_folder_base, "signal_comparison.txt")
write_group_file(signal_path, "Signal", common, lost, new_only, old_events, new_events, old_cfg["label"], new_cfg["label"])

background_path = os.path.join(output_folder_base, "background_comparison.txt")
write_group_file(background_path, "Background", common, lost, new_only, old_events, new_events, old_cfg["label"], new_cfg["label"])

print(f"{old_cfg['label']}: {len(old_events)} events, {new_cfg['label']}: {len(new_events)} events")
print(f"Common: {len(common)}  Lost: {len(lost)}  New: {len(new_only)}")
print(f"Signal comparison written to {signal_path}")
print(f"Background comparison written to {background_path}")
