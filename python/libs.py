import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import uproot
from scipy import stats
from scipy.special import gammaln
from name_index_association import *
from cuts import *

# set the default parameters for matplotlib
plt.rcParams.update({
    "figure.figsize": (8, 6),
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "lines.linewidth": 2,
    "grid.alpha": 0.5
})

def create_passing_matrix(all_events, passing_index, weights=None, parameters=None):
    true_positives = 0
    false_positives = 0
    true_negatives = 0
    false_negatives = 0
    if weights is None:
        for idx, ev in enumerate(all_events):
            if idx in passing_index:
                if ev[aggregate_dict["trueOriginID"]] >= 1:
                    true_positives += 1
                else:
                    false_positives += 1
            else:
                if ev[aggregate_dict["trueOriginID"]] >= 1:
                    false_negatives += 1
                else:
                    true_negatives += 1
    else:
        for idx, ev in enumerate(all_events):
            if idx in passing_index:
                if ev[aggregate_dict["trueOriginID"]] >= 1:
                    true_positives += weights[idx]
                else:
                    false_positives += weights[idx]
            else:
                if ev[aggregate_dict["trueOriginID"]] >= 1:
                    false_negatives += weights[idx]
                else:
                    true_negatives += weights[idx]

    return true_positives, false_positives, true_negatives, false_negatives

def plot_matrix(tp, fp, tn, fn, save_name, title, parameters=None):
    plt.hist2d([0, 1], [0, 1], bins=(2, 2), range=((0, 2), (0, 2)), cmap='Blues')
    plt.text(1.5, 1.5, f"{tp:.2f}", ha='center', va='center', color='white', fontsize=16)
    plt.text(1.5, 0.5, f"{fp:.2f}", ha='center', va='center', color='black', fontsize=16)
    plt.text(0.5, 1.5, f"{fn:.2f}", ha='center', va='center', color='black', fontsize=16)
    plt.text(0.5, 0.5, f"{tn:.2f}", ha='center', va='center', color='white', fontsize=16)
    plt.xticks([0.5, 1.5], ["Not passing", "Passing"], fontsize=12)
    plt.yticks([0.5, 1.5], ["True not neutrino", "True neutrino"], fontsize=12)
    plt.subplots_adjust(left=0.26)
    plt.title(title, fontsize=16)
    plt.savefig(save_name)
    plt.clf()



def get_files(folder, get_truth=False, get_full_array=False, use_combined=False, parameters=None):
    print(f"Opening {folder}")
    events = []
    true_events = []
    full_reco_events = []
    file_info = []
    filenames = []
    counter = 0
    for root, dirs, files in os.walk(folder):
        counter = 0
        print(f"Found {len(files)} files")
        # if "combined.root" in files and use_combined:
        #     files = ["combined.root"]
        # remove combined.root
        # if "combined.root" in files:
        #     files.remove("combined.root")
        if use_combined:
            if "combined.root" not in files:
                print("Combining files")
                print(f"source /exp/dune/app/users/dpullia/neutrino_search_selection/scripts/hadd_script.sh {root}")
                os.system(f"source /exp/dune/app/users/dpullia/neutrino_search_selection/scripts/hadd_script.sh {root}")
            files = ["combined.root"]

        for file in files:
            if file == "combined.root" and len(files) > 1:
                continue
            sys.stdout.write(f"\rPercentage: {counter/len(files):.2%}, File: {os.path.join(root, file)}")
            sys.stdout.flush()
            counter += 1
            if file.endswith(".root"):
                # use the bash timeout command to see if it takes more than 5 seconds to do the head command on the file
                # if it does, then skip the file
                try:
                    os.system(f"timeout 5 head -n 1 {os.path.join(root, file)} > /dev/null")
                except Exception as e:
                    print(f"Error with file {file}: {e}")
                    continue

                try:
                    f = uproot.open(os.path.join(root, file), timeout=1)
                except Exception as e:
                    print(f"Error opening file {file}: {e}")
                    continue
                # Get the aggregate information
                try:
                    tree = f["ana/tree_aggregate"]            
                except KeyError:
                    print(f"KeyError: {file} does not contain ana/tree_aggregate")
                    continue
                eventIDs = tree["eventID"].array()
                vertexX = tree["vertexX"].array()
                vertexY = tree["vertexY"].array()
                vertexZ = tree["vertexZ"].array()
                reconstructedEnergy = tree["reconstructedEnergy"].array()
                directionX = tree["directionX"].array()
                directionY = tree["directionY"].array()
                directionZ = tree["directionZ"].array()
                numberOfHits = tree["numberOfHits"].array()
                numberOfPFParticles = tree["numberOfPFParticles"].array()
                trueOriginID = tree["trueOriginID"].array()
                eventSequenceNumber = tree["eventSequenceNumber"].array()
                passSelectionCriterion = tree["passSelectionCriterion"].array()
                spillStatusFlag = tree["spillStatusFlag"].array()
                eventTimestamp = tree["eventTimestamp"].array()
                totalNumberOfHits = tree["totalNumberOfHits"].array()
                triggerCandidateCount = tree["triggerCandidateCount"].array()
                groundShakeCount = tree["groundShakeCount"].array()
                sumOfLastADCTicks = tree["sumOfLastADCTicks"].array()
                sumOfTriggeredADCTicks = tree["sumOfTriggeredADCTicks"].array()
                passSecondSelectionCriterion = tree["passSecondSelectionCriterion"].array()
                energyDepositedInFirst10cm = tree["energyDepositedInFirst10cm"].array()
                energyDepositedInSecond10cm = tree["energyDepositedInSecond10cm"].array()
                energyDepositedInThird10cm = tree["energyDepositedInThird10cm"].array()
                energyDepositedInFourth10cm = tree["energyDepositedInFourth10cm"].array()
                energyDepositedInFifth10cm = tree["energyDepositedInFifth10cm"].array()
                energyDepositedInSixth10cm = tree["energyDepositedInSixth10cm"].array()
                energyDepositedInSeventh10cm = tree["energyDepositedInSeventh10cm"].array()
                energyDepositedInEighth10cm = tree["energyDepositedInEighth10cm"].array()
                energyDepositedInNinth10cm = tree["energyDepositedInNinth10cm"].array()
                energyDepositedInTenth10cm = tree["energyDepositedInTenth10cm"].array()
                energyDepositedInEleventh10cm = tree["energyDepositedInEleventh10cm"].array()
                energyDepositedInTwelfth10cm = tree["energyDepositedInTwelfth10cm"].array()
                energyDepositedInThirteenth10cm = tree["energyDepositedInThirteenth10cm"].array()
                energyDepositedInFourteenth10cm = tree["energyDepositedInFourteenth10cm"].array()
                energyDepositedInFifteenth10cm = tree["energyDepositedInFifteenth10cm"].array()
                energyDepositedInFirst10cmBefore = tree["energyDepositedInFirst10cmBefore"].array()
                energyDepositedInSecond10cmBefore = tree["energyDepositedInSecond10cmBefore"].array()
                zROIStart = tree["zROIStart"].array()
                zROIEnd = tree["zROIEnd"].array()
                timeROIStart = tree["timeROIStart"].array()
                timeROIEnd = tree["timeROIEnd"].array()
                timeFitMean = tree["timeFitMean"].array()
                timeFitSigma = tree["timeFitSigma"].array()
                directionX2 = tree["directionX2"].array()
                directionY2 = tree["directionY2"].array()
                directionZ2 = tree["directionZ2"].array()
                lengthOfMuonTrack = tree["lengthOfMuonTrack"].array()
                # create a numpy array of events
                event = np.array([
                    eventIDs, vertexX, vertexY, vertexZ, reconstructedEnergy, directionX, directionY, directionZ,
                    numberOfHits, numberOfPFParticles, trueOriginID, eventSequenceNumber, passSelectionCriterion,
                    spillStatusFlag, eventTimestamp, totalNumberOfHits, triggerCandidateCount, groundShakeCount,
                    sumOfLastADCTicks, sumOfTriggeredADCTicks, passSecondSelectionCriterion,
                    energyDepositedInFirst10cm, energyDepositedInSecond10cm, energyDepositedInThird10cm,
                    energyDepositedInFourth10cm, energyDepositedInFifth10cm, energyDepositedInSixth10cm,
                    energyDepositedInSeventh10cm, energyDepositedInEighth10cm, energyDepositedInNinth10cm,
                    energyDepositedInTenth10cm, energyDepositedInEleventh10cm, energyDepositedInTwelfth10cm,
                    energyDepositedInThirteenth10cm, energyDepositedInFourteenth10cm, energyDepositedInFifteenth10cm,
                    energyDepositedInFirst10cmBefore, energyDepositedInSecond10cmBefore, zROIStart, zROIEnd,
                    timeROIStart, timeROIEnd, timeFitMean, timeFitSigma,
                    directionX2, directionY2, directionZ2,
                    lengthOfMuonTrack
                ])
                if (len(eventIDs)==0):
                    print(f"Warning: empty event in file {file}, skipping")
                    os.remove(os.path.join(root, file))
                    continue
                events.append(event)
                # Get file information
                if get_truth:
                    tree = f["ana/file_info"]
                    filePOT = tree["totalAccumulatedPOT"].array()
                    file_info.append(filePOT)

                # Get the truth information
                if get_truth:
                    tree = f["ana/tree_truth"]
                    true_eventIDs = tree["eventID"].array()
                    true_vertexX = tree["vertexX"].array()
                    true_vertexY = tree["vertexY"].array()
                    true_vertexZ = tree["vertexZ"].array()
                    true_energy = tree["energy"].array()
                    true_momentumX = tree["momentumX"].array()
                    true_momentumY = tree["momentumY"].array()
                    true_momentumZ = tree["momentumZ"].array()
                    true_pdgCode = tree["pdgCode"].array()
                    true_motherPdgCode = tree["motherPdgCode"].array()
                    true_subrunTotalPOT = tree["subrunTotalPOT"].array()
                    true_subrunGoodPOT = tree["subrunGoodPOT"].array()
                    true_accumulatedPOT = tree["accumulatedPOT"].array()
                    true_triggerActivityFlag = tree["triggerActivityFlag"].array()
                    true_eventSequenceNumber = tree["eventSequenceNumber"].array()
                    # create a numpy array of events
                    true_event = np.array([
                        true_eventIDs, true_vertexX, true_vertexY, true_vertexZ, true_energy,
                        true_momentumX, true_momentumY, true_momentumZ, true_pdgCode, true_motherPdgCode,
                        true_subrunTotalPOT, true_subrunGoodPOT, true_accumulatedPOT,
                        true_triggerActivityFlag, true_eventSequenceNumber
                    ])
                    true_events.append(true_event)
                    if (true_event.shape[1]) != event.shape[1]:
                        os.remove(os.path.join(root, file))
                        continue
                if get_full_array:
                    tree = f["ana/tree_reco"]
                    # Updated variable names and keys to match the rest of the code
                    eventIDs = tree["eventID"].array()
                    vertexX = tree["vertexX"].array()
                    vertexY = tree["vertexY"].array()
                    vertexZ = tree["vertexZ"].array()
                    pdgCode = tree["pdgCode"].array()
                    energy = tree["energy"].array()
                    directionX = tree["directionX"].array()
                    directionY = tree["directionY"].array()
                    directionZ = tree["directionZ"].array()
                    nHits = tree["numberOfHits"].array()
                    nPFPs = tree["numberOfPFParticles"].array()
                    trueOriginID = tree["trueOriginID"].array()
                    eventNumber = tree["eventSequenceNumber"].array()
                    passCut = tree["passSelectionCriterion"].array()
                    
                    # create a numpy array of events    
                    full_event = np.array([eventIDs, vertexX, vertexY, vertexZ, pdgCode, energy, directionX, directionY, directionZ, nHits, nPFPs, trueOriginID, eventNumber, passCut])
                    full_event = np.transpose(full_event)
                    event_ids= event[aggregate_dict["eventID"], :]
                    split_indices = []
                    for unique_event_id in event_ids[:-1]:
                        indices = np.where(full_event[:, full_reco_dict["eventID"]] == unique_event_id)[0]
                        if len(indices) > 0:
                            split_indices.append(indices[-1] + 1)
                        else:
                            if len(split_indices) == 0:
                                split_indices.append(0)
                            else:
                                split_indices.append(split_indices[-1])

                    # print("Split indices: ", split_indices)
                    sublist = np.split(full_event, split_indices, axis=0)
                    # print("Sublist: ", [len(sublist[i]) for i in range(len(sublist))])
                    if len(sublist) != event.shape[1]:
                        if event.shape[1] == 0 and sublist[0].shape[0] == 0:
                            print("Warning: empty event, skipping")
                            os.remove(os.path.join(root, file))
                            continue
                        print("Error: sublist length does not match events length")
                        print(f"sublist length: {len(sublist)}, events length: {event.shape}, true event length: {true_event.shape}")
                        print(f"event ids: {event[aggregate_dict['eventID'], :]}")
                        sys.exit(1)

                    full_reco_events.extend(sublist)

                # filenames.extend([os.path.join(root, file)] * event.shape[1])
                base_name = file.split(".root")[0]
                base_name = os.path.join(root, base_name)
                filenames.extend([base_name+f"_{ev[aggregate_dict['eventID']]}" for ev in event.T])


    print("\n")
    events = np.concatenate(events, axis=1)
    events = np.transpose(events)
    if get_truth:
        # true_events = np.array(true_events)
        # true_events = np.squeeze(true_events)
        # if len(files) == 1:
        #     true_events = np.transpose(true_events)
        true_events = np.concatenate(true_events, axis=1).T

    if get_full_array:
        full_reco_events = np.array(full_reco_events)

    file_info = np.array(file_info)

    filenames = np.array(filenames)
    return events, true_events, full_reco_events, filenames, file_info



def load_from_folder(folders, get_truth=False, get_full_array=False, use_combined=False, weights_mode="POT", parameters=None, spill_status="off"):
    all_events = []
    all_true_events = []
    all_full_events = []
    all_weights = []
    all_labels = []
    filenames = []
    total_time = 0
    total_POT_MC = 0

    all_folders = []
    for root, dirs, files in os.walk(folders):
        for d in dirs:
            all_folders.append(os.path.join(root, d))
    if all_folders == []:
        all_folders = [folders]
    max_length = max([len(f.split("/")) for f in all_folders])
    all_folders = [f for f in all_folders if len(f.split("/")) == max_length]


    for f in all_folders:
        ev, true_ev, ev_full, ev_filenames, ev_file_info = get_files(f, get_truth=get_truth, get_full_array=get_full_array, use_combined=use_combined, parameters=parameters)
        print(len(ev), len(true_ev), len(ev_full), len(ev_filenames))

        if spill_status == "off" and parameters["TP_RATE"] != "mc":
            # Get the events that are not in the spill
            spill_index = np.where(ev[:, aggregate_dict["spillStatusFlag"]] == 0)[0]
            ev = ev[spill_index]
            if get_truth:
                true_ev = true_ev[spill_index]
            if get_full_array:
                ev_full = [ev_full[i] for i in spill_index] 
            ev_filenames = [ev_filenames[i] for i in spill_index]
        elif spill_status == "on" and parameters["TP_RATE"] != "mc":
            # Get the events that are in the spill
            spill_index = np.where(ev[:, aggregate_dict["spillStatusFlag"]] == 1)[0]
            ev = ev[spill_index]
            if get_truth:
                true_ev = true_ev[spill_index]
            if get_full_array:
                ev_full = [ev_full[i] for i in spill_index]
            ev_filenames = [ev_filenames[i] for i in spill_index]
        elif spill_status == "all" or parameters["TP_RATE"] == "mc":
            pass
        else:
            print("Spill status not recognized")
            sys.exit(1)     
        # Get the flavour and decay from the file name
        flavour = f.split("/")[-2]

        try:
            if weights_mode == "POT_1hour":
                total_POT_MC += np.sum(ev_file_info)
            elif weights_mode == "none":
                weight = 1
            elif weights_mode == "1hour":
                pass
            else:
                print("Weights mode not recognized")
                sys.exit(1)
        except KeyError as e:
            print(f"KeyError: Run not found in parameters: {e}")
            print(f"Using weight = 1")
            weight = 1
            total_time = 1  # Set total time to 1 to avoid division by zero
        all_events.append(ev)
        all_true_events.append(true_ev)
        all_full_events.extend(ev_full)
        all_labels.append(np.array([f"{flavour}"] * ev.shape[0]))
        filenames.extend(ev_filenames)

    all_events = np.concatenate(all_events)
    all_true_events = np.concatenate(all_true_events)
    all_labels = np.concatenate(all_labels)
    if weights_mode == "1hour":
        job = f.split("/")[-3]
        all_weights = np.ones(all_events.shape[0]) * (1/parameters["run_parameters"][f"spill_{spill_status}"]["total_time"])
    elif weights_mode == "POT_1hour":
        all_weights = np.array([parameters["run_parameters"][f"spill_{spill_status}"]["total_POT"] / total_POT_MC] * all_events.shape[0])
        all_weights = all_weights * (1/parameters["run_parameters"][f"spill_{spill_status}"]["total_time"])
    else:
        all_weights = np.concatenate(all_weights)

    if spill_status == "on" and parameters["TP_RATE"] != "mc":
        files = parameters["tp_monitor_files"]
        monitor_tp_time = []
        monitor_tp_rate = []
        for file in files:
            data = np.loadtxt(file, delimiter=",", usecols=(1,2))
            monitor_tp_time.extend(data[:, 1])
            monitor_tp_rate.extend(data[:, 0])
        monitor_tp_time = np.array(monitor_tp_time) * 16E-9
        monitor_tp_rate = np.array(monitor_tp_rate)


        if parameters["TP_RATE"] == "low":
            index_all = np.array([tp_rate_below_threshold(ts, 2.1E6, monitor_tp_time, monitor_tp_rate, parameters=parameters) for ts in all_events[:, aggregate_dict["eventTimestamp"]]])
            print(np.sum(index_all))
            all_events = all_events[index_all]
            all_true_events = all_true_events[index_all]
            all_labels = all_labels[index_all]
            all_weights = all_weights[index_all]
            filenames = [filenames[i] for i in range(len(filenames)) if index_all[i]]

        elif parameters["TP_RATE"] == "high":
            index_all = np.array([tp_rate_above_threshold(ts, 2.1E6, monitor_tp_time, monitor_tp_rate, parameters=parameters) for ts in all_events[:, aggregate_dict["eventTimestamp"]]])
            print(np.sum(index_all))
            all_events = all_events[index_all]
            all_true_events = all_true_events[index_all]
            all_labels = all_labels[index_all]
            all_weights = all_weights[index_all]
            filenames = [filenames[i] for i in range(len(filenames)) if index_all[i]]

    print(f"Total events: {len(all_events)}, Total true events: {len(all_true_events)}, Total full events: {len(all_full_events)}, Total weights: {len(all_weights)}, Total labels: {len(all_labels)}, Total filenames: {len(filenames)}")

    return all_events, all_true_events, all_full_events, all_weights, all_labels, filenames


def triple_hist(hist_sig, sig_weights, hist_bkg, bkg_weights, parameters=None, output_folder=None, spill_status="off", plot_name="sig_bkg_energy_spectrum", title="Background vs Signal Energy Spectrum", range=(0, 100), bins=50, xlabel="Energy (GeV)", ylabel="Counts", label_sig = "Signal", label_bkg = "Background", total_time_with_correct_tps_on=(2.93733+5.5266)):
    color_sig = "orange"
    # color_sig = "green"
    color_bkg = "blue"
   
    # Step 1: Compare energy spectrums between signal off spill and background
    output_folder_normalized = output_folder+"/output_test_normalized/"
    output_folder_absolute = output_folder+"/output_test_absolute/"
    output_folder_1hour = output_folder+"/output_test_1hour/"
    output_folder_full_time = output_folder+"/output_test_full_time/"

    hist_bkg_pure, hist_bkg_pure_bins = np.histogram(hist_bkg, bins=bins, range=range)
    # Compute the errors
    hist_bkg_errors_absolute = np.sqrt(hist_bkg_pure)
    hist_bkg_errors_normalized = hist_bkg_errors_absolute / hist_bkg_pure.sum() if hist_bkg_pure.sum() > 0 else 0
    hist_bkg_errors_1hour = hist_bkg_errors_absolute * bkg_weights[0]




    bin_width = (range[1] - range[0]) / bins

    plt.figure(figsize=(9, 6))


    normalization_sig = sig_weights.sum() if sig_weights.sum() > 0 else 1
    normalization_bkg = bkg_weights.sum() if bkg_weights.sum() > 0 else 1

    plt.hist(hist_bkg, bins=bins, label=label_bkg, density=False, range=range, color=color_bkg, alpha = 0.5, weights=bkg_weights/normalization_bkg)
    plt.hist(hist_sig, bins=bins, label=label_sig, density=False, range=range, color=color_sig, alpha = 0.5, weights=sig_weights/normalization_sig)
    plt.hist(hist_bkg, bins=bins, density=False, range=range, histtype='step', color=color_bkg, weights=bkg_weights/normalization_bkg)
    plt.hist(hist_sig, bins=bins, density=False, range=range, histtype='step', color=color_sig, weights=sig_weights/normalization_sig)
    # add error bars
    # plt.errorbar(hist_bkg_pure_bins[:-1]+bin_width*0.5, hist_bkg_pure / hist_bkg_pure.sum() if hist_bkg_pure.sum() > 0 else 0, yerr=hist_bkg_errors_normalized, fmt='none', ecolor='black', capsize=5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"{title}: Normalized.")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(output_folder_normalized, f"{plot_name}_normalized.png"))
    plt.clf()

    plt.hist(hist_bkg, bins=bins, label=label_bkg, density=False, range=range, color=color_bkg, alpha=0.5)
    plt.hist(hist_sig, bins=bins, label=label_sig, density=False, range=range, color=color_sig, alpha=0.5)
    plt.hist(hist_bkg, bins=bins, density=False, range=range, histtype='step', color=color_bkg)
    plt.hist(hist_sig, bins=bins, density=False, range=range, histtype='step', color=color_sig)
    # add error bars
    plt.errorbar(hist_bkg_pure_bins[:-1]+bin_width*0.5, hist_bkg_pure, yerr=hist_bkg_errors_absolute, fmt='none', ecolor='black', capsize=5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"{title}: Absolute.")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(output_folder_absolute, f"{plot_name}_absolute.png"))
    plt.clf()

    plt.hist(hist_bkg, bins=bins, label=label_bkg, density=False, range=range, weights=bkg_weights, color=color_bkg, alpha=0.5)
    plt.hist(hist_sig, bins=bins, label=label_sig, density=False, range=range, weights=sig_weights, color=color_sig, alpha=0.5)
    plt.hist(hist_bkg, bins=bins, density=False, range=range, histtype='step', color=color_bkg, weights=bkg_weights)
    plt.hist(hist_sig, bins=bins, density=False, range=range, histtype='step', color=color_sig, weights=sig_weights)
    # add error bars
    plt.errorbar(hist_bkg_pure_bins[:-1]+bin_width*0.5, hist_bkg_pure * bkg_weights[0], yerr=hist_bkg_errors_1hour, fmt='none', ecolor='black', capsize=5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"{title}: 1 hour.")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(output_folder_1hour, f"{plot_name}_1hour.png"))
    plt.clf()

    
    bkg_weights_full = bkg_weights * total_time_with_correct_tps_on
    sig_weights_full = sig_weights * total_time_with_correct_tps_on
    hist_bkg_errors_full = hist_bkg_errors_absolute * bkg_weights_full[0]
    hist_sig_pure, hist_sig_pure_bins = np.histogram(hist_sig, bins=bins, range=range)
    hist_sig_errors_absolute = np.sqrt(hist_sig_pure)
    hist_sig_errors_full = hist_sig_errors_absolute * sig_weights_full[0]






    plt.hist(hist_bkg, bins=bins, label=label_bkg, density=False, range=range, weights=bkg_weights_full, color=color_bkg, alpha=0.5)
    plt.hist(hist_sig, bins=bins, label=label_sig, density=False, range=range, weights=sig_weights_full, color=color_sig, alpha=0.5)
    plt.hist(hist_bkg, bins=bins, density=False, range=range, histtype='step', color=color_bkg, weights=bkg_weights_full)
    plt.hist(hist_sig, bins=bins, density=False, range=range, histtype='step', color=color_sig, weights=sig_weights_full)
    # add error bars
    plt.errorbar(hist_bkg_pure_bins[:-1]+bin_width*0.5, hist_bkg_pure * bkg_weights_full[0], yerr=hist_bkg_errors_full, fmt='none', ecolor='black', capsize=2)
    # plt.errorbar(hist_sig_pure_bins[:-1]+bin_width*0.5, hist_sig_pure * sig_weights_full[0], yerr=hist_sig_errors_full, fmt='none', ecolor='black', capsize=2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"{title}: Full data time.")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(output_folder_full_time, f"{plot_name}_full_sig_time.png"))
    plt.clf()
    plt.close()


def triple_hist_stack(hist_sig, sig_weights, hist_bkg, bkg_weights, hist_mc, mc_weights, parameters=None, output_folder=None, spill_status="off", plot_name="sig_bkg_energy_spectrum", title="Background vs Signal Energy Spectrum", range=(0, 100), bins=50, xlabel="Energy (GeV)", ylabel="Counts", label_sig = "Signal", label_bkg = "Background",total_time_with_correct_tps_on=(2.93733+5.5266)):
    # color_sig = "orange"
    color_sig = "black"
    color_bkg = "blue"
    color_mc = "orange"
    # Step 1: Compare energy spectrums between signal off spill and background
    output_folder_normalized = output_folder+"/output_test_normalized/"
    output_folder_absolute = output_folder+"/output_test_absolute/"
    output_folder_1hour = output_folder+"/output_test_1hour/"
    output_folder_full_time = output_folder+"/output_test_full_time/"

    hist_bkg_pure, hist_bkg_pure_bins = np.histogram(hist_bkg, bins=bins, range=range)
    # Compute the errors
    hist_bkg_errors_absolute = np.sqrt(hist_bkg_pure)
    hist_bkg_errors_normalized = hist_bkg_errors_absolute / hist_bkg_pure.sum() if hist_bkg_pure.sum() > 0 else 0
    hist_bkg_errors_1hour = hist_bkg_errors_absolute * bkg_weights[0]




    bin_width = (range[1] - range[0]) / bins

    plt.figure(figsize=(10, 6))


    normalization_sig = sig_weights.sum() if sig_weights.sum() > 0 else 1
    normalization_bkg = bkg_weights.sum() if bkg_weights.sum() > 0 else 1
    normalization_mc = mc_weights.sum() if mc_weights.sum() > 0 else 1

    plt.hist([hist_bkg, hist_mc], bins=bins, label=[label_bkg, "Neutrino MC"], density=False, range=range, color=[color_bkg, color_mc], alpha=0.5, weights=[bkg_weights/normalization_bkg, mc_weights/normalization_mc], stacked=True)
    # plt.hist(hist_sig, bins=bins, label=label_sig, density=False, range=range, color=color_sig, alpha = 0.5, weights=sig_weights/normalization_sig)
    plt.hist([hist_bkg, hist_mc], bins=bins, density=False, range=range, histtype='step', color=[color_bkg, color_mc], weights=[bkg_weights/normalization_bkg, mc_weights/normalization_mc], stacked=True)
    plt.hist(hist_sig, bins=bins, density=False, range=range, histtype='step', color=color_sig, weights=sig_weights/normalization_sig, label=label_sig, linewidth=3)
    # add error bars
    # plt.errorbar(hist_bkg_pure_bins[:-1]+bin_width*0.5, hist_bkg_pure / hist_bkg_pure.sum() if hist_bkg_pure.sum() > 0 else 0, yerr=hist_bkg_errors_normalized, fmt='none', ecolor='black', capsize=5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"{title}: Normalized.")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(output_folder_normalized, f"{plot_name}_normalized.png"))
    plt.clf()

    plt.hist([hist_bkg, hist_mc], bins=bins, label=[label_bkg, "Neutrino MC"], density=False, range=range, color=[color_bkg, color_mc], alpha=0.5, weights=[bkg_weights, mc_weights], stacked=True)
    # plt.hist(hist_sig, bins=bins, label=label_sig, density=False, range=range, color=color_sig, alpha=0.5)
    plt.hist([hist_bkg, hist_mc], bins=bins, density=False, range=range, histtype='step', color=[color_bkg, color_mc], weights=[bkg_weights, mc_weights], stacked=True)
    plt.hist(hist_sig, bins=bins, density=False, range=range, histtype='step', color=color_sig, label=label_sig, linewidth=3)
    # add error bars
    plt.errorbar(hist_bkg_pure_bins[:-1]+bin_width*0.5, hist_bkg_pure, yerr=hist_bkg_errors_absolute, fmt='none', ecolor='black', capsize=5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"{title}: Absolute.")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(output_folder_absolute, f"{plot_name}_absolute.png"))
    plt.clf()

    plt.hist([hist_bkg, hist_mc], bins=bins, label=[label_bkg, "Neutrino MC"], density=False, range=range, weights=[bkg_weights, mc_weights], color=[color_bkg, color_mc], alpha=0.5, stacked=True)
    # plt.hist(hist_sig, bins=bins, label=label_sig, density=False, range=range, weights=sig_weights, color=color_sig, alpha=0.5)
    plt.hist([hist_bkg, hist_mc], bins=bins, density=False, range=range, histtype='step', color=[color_bkg, color_mc], weights=[bkg_weights, mc_weights], stacked=True)
    plt.hist(hist_sig, bins=bins, density=False, range=range, histtype='step', color=color_sig, weights=sig_weights, label=label_sig, linewidth=3)
    # add error bars
    plt.errorbar(hist_bkg_pure_bins[:-1]+bin_width*0.5, hist_bkg_pure * bkg_weights[0], yerr=hist_bkg_errors_1hour, fmt='none', ecolor='black', capsize=5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"{title}: 1 hour.")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(output_folder_1hour, f"{plot_name}_1hour.png"))
    plt.clf()

    
    bkg_weights_full = bkg_weights * total_time_with_correct_tps_on
    sig_weights_full = sig_weights * total_time_with_correct_tps_on
    mc_weights_full = mc_weights * total_time_with_correct_tps_on
    hist_bkg_errors_full = hist_bkg_errors_absolute * bkg_weights_full[0]
    hist_sig_pure, hist_sig_pure_bins = np.histogram(hist_sig, bins=bins, range=range)
    hist_sig_errors_absolute = np.sqrt(hist_sig_pure)
    hist_sig_errors_full = hist_sig_errors_absolute * sig_weights_full[0]






    plt.hist([hist_bkg, hist_mc], bins=bins, label=[label_bkg, "Neutrino MC"], density=False, range=range, weights=[bkg_weights_full, mc_weights_full], color=[color_bkg, color_mc], alpha=0.5, stacked=True)
    # plt.hist(hist_sig, bins=bins, label=label_sig, density=False, range=range, weights=sig_weights_full, color=color_sig, alpha=0.5)
    plt.hist([hist_bkg, hist_mc], bins=bins, density=False, range=range, histtype='step', color=[color_bkg, color_mc], weights=[bkg_weights_full, mc_weights_full], stacked=True)
    plt.hist(hist_sig, bins=bins, density=False, range=range, histtype='step', color=color_sig, weights=sig_weights_full, label=label_sig, linewidth=3)
    # add error bars
    plt.errorbar(hist_bkg_pure_bins[:-1]+bin_width*0.5, hist_bkg_pure * bkg_weights_full[0], yerr=hist_bkg_errors_full, fmt='none', ecolor='black', capsize=2)
    plt.errorbar(hist_sig_pure_bins[:-1]+bin_width*0.5, hist_sig_pure * sig_weights_full[0], yerr=hist_sig_errors_full, fmt='none', ecolor='black', capsize=2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"{title}: Full data time.")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(output_folder_full_time, f"{plot_name}_full_sig_time.png"))
    plt.clf()
    plt.close()



def apply_all_cuts( sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
                    bkg_events, bkg_weights, bkg_labels, bkg_filenames, cuts, skip_cut=None):
    # Helper to apply all cuts except skip_cut (by name)
    # Returns: sig_events, sig_weights, sig_true_events, sig_labels, bkg_events, bkg_weights, bkg_labels

    # Work on copies
    s_ev, s_w, s_true, s_lab, s_fnames = sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames
    b_ev, b_w, b_lab, b_fnames = bkg_events, bkg_weights, bkg_labels, bkg_filenames

    for name, cutfunc in cuts:
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


def apply_all_cuts_MC( sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
                    bkg_events, bkg_weights, bkg_labels, bkg_filenames, mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames, cuts, skip_cut=None):
    # Helper to apply all cuts except skip_cut (by name)
    # Returns: sig_events, sig_weights, sig_true_events, sig_labels, bkg_events, bkg_weights, bkg_labels

    # Work on copies
    s_ev, s_w, s_true, s_lab, s_fnames = sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames
    b_ev, b_w, b_lab, b_fnames = bkg_events, bkg_weights, bkg_labels, bkg_filenames
    mc_ev, mc_w, mc_true, mc_lab, mc_fnames = mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames

    for name, cutfunc in cuts:
        if skip_cut is not None and name == skip_cut:
            continue
        idx_s, idx_b, idx_mc = cutfunc(s_ev, b_ev, mc_ev)
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

        mc_ev = mc_ev[idx_mc]
        mc_w = mc_w[idx_mc]
        if len(mc_true) > 0:
            mc_true = mc_true[idx_mc]
        mc_lab = mc_lab[idx_mc]
        mc_fnames = [mc_fnames[i] for i in idx_mc]
        # print(f"After cut {name}: N signal events = {len(s_ev)}, N background events = {len(b_ev)}, N MC events = {len(mc_ev)}")

    return s_ev, s_w, s_true, s_lab, s_fnames, b_ev, b_w, b_lab, b_fnames, mc_ev, mc_w, mc_true, mc_lab, mc_fnames



def create_table_cuts_with_mc( sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
                    bkg_events, bkg_weights, bkg_labels, bkg_filenames, mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames, 
                    cuts_names, skip_cut=None, output_folder=None):
    # Work on copies
    s_ev, s_w, s_true, s_lab, s_fnames = sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames
    b_ev, b_w, b_lab, b_fnames = bkg_events, bkg_weights, bkg_labels, bkg_filenames
    mc_ev, mc_w, mc_true, mc_lab, mc_fnames = mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames

    absolute_sig, absolute_bkg, absolute_mc = [len(s_ev)], [len(b_ev)], [len(mc_ev)]
    normalized_sig, normalized_bkg, normalized_mc = [s_w.sum()], [b_w.sum()], [mc_w.sum()]

    for index, name in enumerate(cuts_names):
        cutfunc = next((func for n, func in cuts_with_mc if n == name), None)
        short_name = shorter_cut_names[index] if index < len(shorter_cut_names) else name
        if cutfunc is None:
            continue
        if skip_cut is not None and name == skip_cut:
            continue
        idx_s, idx_b, idx_mc = cutfunc(s_ev, b_ev, mc_ev)
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

        mc_ev = mc_ev[idx_mc]
        mc_w = mc_w[idx_mc]
        if len(mc_true) > 0:
            mc_true = mc_true[idx_mc]
        mc_lab = mc_lab[idx_mc]
        mc_fnames = [mc_fnames[i] for i in idx_mc]
        absolute_sig.append(len(s_ev))
        absolute_bkg.append(len(b_ev))
        absolute_mc.append(len(mc_ev))
        normalized_sig.append(s_w.sum())
        normalized_bkg.append(b_w.sum())
        normalized_mc.append(mc_w.sum())
        # print(f"After cut {name}: N signal events = {len(s_ev)}, N background events = {len(b_ev)}, N MC events = {len(mc_ev)}")
    with open(os.path.join(output_folder, "cut_table_absolute.txt"), "w") as f:
        # setup a latex table
        f.write("\\begin{tabular}{|c|c|c|c|}\n")
        f.write("\\hline\n")
        f.write("Cut & Signal Events & Background Events & MC Events \\\\\n")
        f.write("\\hline\n")
        f.write(f"Initial & {absolute_sig[0]} & {absolute_bkg[0]} & {absolute_mc[0]} \\\\\n")
        for index, name in enumerate(cuts_names):
            short_name = shorter_cut_names[index] if index < len(shorter_cut_names) else name
            f.write(f"{short_name} & {absolute_sig[index+1]} & {absolute_bkg[index+1]} & {absolute_mc[index+1]} \\\\\n")
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
    with open(os.path.join(output_folder, "cut_table_normalized.txt"), "w") as f:
        # setup a latex table
        f.write("\\begin{tabular}{|c|c|c|c|}\n")
        f.write("\\hline\n")
        f.write("Cut & Signal Events & Background Events & MC Events \\\\\n")
        f.write("\\hline\n")
        f.write(f"Initial & {normalized_sig[0]:.2f} & {normalized_bkg[0]:.2f} & {normalized_mc[0]:.2f} \\\\\n")
        for index, name in enumerate(cuts_names):
            short_name = shorter_cut_names[index] if index < len(shorter_cut_names) else name
            f.write(f"{short_name} & {normalized_sig[index+1]:.2f} & {normalized_bkg[index+1]:.2f} & {normalized_mc[index+1]:.2f} \\\\\n")
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")


    return s_ev, s_w, s_true, s_lab, s_fnames, b_ev, b_w, b_lab, b_fnames, mc_ev, mc_w, mc_true, mc_lab, mc_fnames



def create_table_cuts( sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
                    bkg_events, bkg_weights, bkg_labels, bkg_filenames, 
                    cuts_names, skip_cut=None, output_folder=None):
    # Work on copies
    s_ev, s_w, s_true, s_lab, s_fnames = sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames
    b_ev, b_w, b_lab, b_fnames = bkg_events, bkg_weights, bkg_labels, bkg_filenames

    absolute_sig, absolute_bkg = [len(s_ev)], [len(b_ev)]
    normalized_sig, normalized_bkg = [s_w.sum()], [b_w.sum()]

    for index, name in enumerate(cuts_names):
        cutfunc = next((func for n, func in cuts if n == name), None)
        short_name = shorter_cut_names[index] if index < len(shorter_cut_names) else name
        if cutfunc is None:
            continue
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

        absolute_sig.append(len(s_ev))
        absolute_bkg.append(len(b_ev))
        normalized_sig.append(s_w.sum())
        normalized_bkg.append(b_w.sum())
        # print(f"After cut {name}: N signal events = {len(s_ev)}, N background events = {len(b_ev)}, N MC events = {len(mc_ev)}")
    with open(os.path.join(output_folder, "cut_table_absolute.txt"), "w") as f:
        # setup a latex table
        f.write("\\begin{tabular}{|c|c|c|}\n")
        f.write("\\hline\n")
        f.write("Cut & Signal Events & Background Events \\\\\n")
        f.write("\\hline\n")
        f.write(f"Initial & {absolute_sig[0]} & {absolute_bkg[0]} \\\\\n")
        for index, name in enumerate(cuts_names):
            short_name = shorter_cut_names[index] if index < len(shorter_cut_names) else name
            f.write(f"{short_name} & {absolute_sig[index+1]} & {absolute_bkg[index+1]} \\\\\n")
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
    with open(os.path.join(output_folder, "cut_table_normalized.txt"), "w") as f:
        # setup a latex table
        f.write("\\begin{tabular}{|c|c|c|}\n")
        f.write("\\hline\n")
        f.write("Cut & Signal Events & Background Events \\\\\n")
        f.write("\\hline\n")
        f.write(f"Initial & {normalized_sig[0]:.2f} & {normalized_bkg[0]:.2f} \\\\\n")
        for index, name in enumerate(cuts_names):
            short_name = shorter_cut_names[index] if index < len(shorter_cut_names) else name
            f.write(f"{short_name} & {normalized_sig[index+1]:.2f} & {normalized_bkg[index+1]:.2f} \\\\\n")
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")


    return s_ev, s_w, s_true, s_lab, s_fnames, b_ev, b_w, b_lab, b_fnames



def tp_rate_below_threshold(timestamp, threshold, monitor_tp_time, monitor_tp_rate, parameters=None):
    # find closest timestamp to element in monitor_tp_time
    timestamp = timestamp * 1E-9
    closest_index = np.argmin(np.abs(monitor_tp_time - timestamp))
    if max(monitor_tp_rate[closest_index-3:closest_index+3]) < threshold:
        return True
    return False


def tp_rate_above_threshold(timestamp, threshold, monitor_tp_time, monitor_tp_rate, parameters=None):
    # find closest timestamp to element in monitor_tp_time
    timestamp = timestamp * 1E-9
    closest_index = np.argmin(np.abs(monitor_tp_time - timestamp))
    if max(monitor_tp_rate[closest_index-3:closest_index+3]) > threshold:
        return True
    return False






