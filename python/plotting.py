import os
import sys

sys.path.append("../python")
from libs import *
from name_index_association import *
from cuts import *

def run_all_triple_hist(
    sig_events, sig_weights, bkg_events, bkg_weights, output_folder_base, label_sig, label_bkg, nbins, spill_status = f"spill_status", apply_cuts=False, cuts=cuts,  total_time_with_correct_tps_on=(2.93733+5.5266)
    ):
    
    sig_events_orig = sig_events.copy()
    sig_weights_orig = sig_weights.copy()
    sig_true_events_orig = np.ones((sig_events.shape[0], 1)) * -1  # Placeholder, as true events are not provided
    sig_labels_orig = np.zeros((sig_events.shape[0], 1))  # Placeholder, as labels are not provided
    sig_filenames_orig = np.zeros((sig_events.shape[0], 1))  # Placeholder, as filenames are not provided   
    bkg_events_orig = bkg_events.copy()
    bkg_weights_orig = bkg_weights.copy()
    bkg_labels_orig = np.zeros((bkg_events.shape[0], 1))  # Placeholder, as labels are not provided
    bkg_filenames_orig = np.zeros((bkg_events.shape[0], 1))  # Placeholder, as filenames are not provided
    # ----------------------------------------------
    # Step 1: Compare energy spectrums between {label_sig} off spill and {label_bkg}
    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut=None
        )

    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["reconstructedEnergy"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["reconstructedEnergy"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_energy_spectrum",
        title=f"{label_bkg} vs {label_sig}. Energy Spectrum",
        range=(0, 100),
        bins=nbins,
        xlabel="Energy (GeV)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # ----------------------------------------------
    # # Step 2: Compare the direction spectrums between {label_sig} off spill and {label_bkg}
    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut="max_directionZ"
        )


    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["directionX"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["directionX"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_directionX_spectrum",
        title=f"{label_bkg} vs {label_sig}. Direction X Spectrum",
        range=(-1, 1),
        bins=nbins,
        xlabel="Direction X",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["directionY"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["directionY"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_directionY_spectrum",
        title=f"{label_bkg} vs {label_sig}. Direction Y Spectrum",
        range=(-1, 1),
        bins=nbins,
        xlabel="Direction Y",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["directionZ"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["directionZ"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_directionZ_spectrum",
        title=f"{label_bkg} vs {label_sig}. Direction Z Spectrum",
        range=(-1, 1),
        bins=nbins,
        xlabel="Direction Z",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    # ----------------------------------------------
    # # Step 2.5: Compare the direction spectrums between {label_sig} off spill and {label_bkg}
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["directionX2"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["directionX2"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_directionX2_spectrum",
        title=f"{label_bkg} vs {label_sig}. Direction X2 Spectrum",
        range=(-1, 1),
        bins=nbins,
        xlabel="Direction X2",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["directionY2"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["directionY2"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_directionY2_spectrum",
        title=f"{label_bkg} vs {label_sig}. Direction Y2 Spectrum",
        range=(-1, 1),
        bins=nbins,
        xlabel="Direction Y2",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["directionZ2"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["directionZ2"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_directionZ2_spectrum",
        title=f"{label_bkg} vs {label_sig}. Direction Z2 Spectrum",
        range=(-1, 1),
        bins=nbins,
        xlabel="Direction Z2",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    # ----------------------------------------------
    # # Step 2.6: Do an array with the max value of directionZ and directionZ2

    triple_hist(
        hist_sig=np.maximum(
            sig_events[:, aggregate_dict["directionZ"]],
            sig_events[:, aggregate_dict["directionZ2"]]
        ),
        sig_weights=sig_weights,
        hist_bkg=np.maximum(
            bkg_events[:, aggregate_dict["directionZ"]],
            bkg_events[:, aggregate_dict["directionZ2"]]
        ),
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_max_directionZ_directionZ2_spectrum",
        title=f"{label_bkg} vs {label_sig}. $\\cos(\\theta_{{\\mathrm{{beam}}}})_{{\\mathrm{{max}}}}$ Spectrum",
        range=(-1, 1),
        bins=nbins,
        xlabel="$\\cos(\\theta_{\\mathrm{beam}})_{{\\mathrm{{max}}}}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # print to text the number of events with directionZ > 0.97 or directionZ2 > 0.97
    hist_sig_max, bin_sig_max = np.histogram(
        np.maximum(
            sig_events[:, aggregate_dict["directionZ"]],
            sig_events[:, aggregate_dict["directionZ2"]]
        ),
        bins=80,
        range=(-1, 1),
        weights=sig_weights
    )
    hist_bkg_max, bin_bkg_max = np.histogram(
        np.maximum(
            bkg_events[:, aggregate_dict["directionZ"]],
            bkg_events[:, aggregate_dict["directionZ2"]]
        ),
        bins=80,
        range=(-1, 1),
        weights=bkg_weights
    )
    print("Hist signal:", hist_sig_max)
    print("Bin signal:", bin_sig_max)
    print("Hist background:", hist_bkg_max)
    print("Bin background:", bin_bkg_max)
    # save to file
    with open(os.path.join(output_folder_base, "max_directionZ_directionZ2.txt"), "w") as f:
        f.write("Bin edges:\n")
        f.write(", ".join([str(b) for b in bin_sig_max]) + "\n")
        f.write("Signal histogram:\n")
        f.write(", ".join([str(h) for h in hist_sig_max]) + "\n")
        f.write("Background histogram:\n")
        f.write(", ".join([str(h) for h in hist_bkg_max]) + "\n")
        

    # 2d plot, Z vs Z2
    plt.figure(figsize=(8, 6))
    plt.hist2d(
        sig_events[:, aggregate_dict["directionZ"]],
        sig_events[:, aggregate_dict["directionZ2"]],
        weights=sig_weights,
        bins=nbins,
        range=[[-1, 1], [-1, 1]],
        cmap='viridis',
        label=label_sig,
        cmin=0.001
    )
    plt.colorbar(label='Counts')
    plt.xlabel('$\\cos(\\theta_{{\\mathrm{{beam}}}})$ Hits')
    plt.ylabel('$\\cos(\\theta_{{\\mathrm{{beam}}}})$ PFP')
    plt.title(f"{label_sig} $\\cos(\\theta_{{\\mathrm{{beam}}}})$ Hits vs $\\cos(\\theta_{{\\mathrm{{beam}}}})$ PFP")
    plt.savefig(os.path.join(output_folder_base, "sig_directionZ_vs_directionZ2.png"))
    plt.close()
    plt.figure(figsize=(8, 6))
    plt.hist2d(
        bkg_events[:, aggregate_dict["directionZ"]],
        bkg_events[:, aggregate_dict["directionZ2"]],
        weights=bkg_weights,
        bins=nbins,
        range=[[-1, 1], [-1, 1]],
        cmap='plasma',
        label=label_bkg,
    )
    plt.colorbar(label='Counts')
    plt.xlabel('$\\cos(\\theta_{{\\mathrm{{beam}}}})$ Hits')
    plt.ylabel('$\\cos(\\theta_{{\\mathrm{{beam}}}})$ PFP')
    plt.title(f"{label_bkg} $\\cos(\\theta_{{\\mathrm{{beam}}}})$ Hits vs $\\cos(\\theta_{{\\mathrm{{beam}}}})$ PFP")
    plt.savefig(os.path.join(output_folder_base, "bkg_directionZ_vs_directionZ2.png"))
    plt.close()

    # on the same histogram, directionZ vs directionZ2 vs max direction
    plt.figure(figsize=(8, 6))
    plt.hist(
        bkg_events[:, aggregate_dict["directionZ"]],
        bins=nbins,
        range=(-1, 1),
        weights=bkg_weights/np.sum(bkg_weights),
        alpha=0.4,
        label=f"{label_bkg} $\\cos(\\theta_{{\\mathrm{{beam}}}})$ Hits",
        color='blue',
    )
    plt.hist(
        bkg_events[:, aggregate_dict["directionZ2"]],
        bins=nbins,
        range=(-1, 1),
        weights=bkg_weights/np.sum(bkg_weights),
        alpha=0.4,
        label=f"{label_bkg} $\\cos(\\theta_{{\\mathrm{{beam}}}})$ PFP",
        color='green',
    )
    plt.hist(
        np.maximum(
            bkg_events[:, aggregate_dict["directionZ"]],
            bkg_events[:, aggregate_dict["directionZ2"]]
        ),
        bins=nbins,
        range=(-1, 1),
        weights=bkg_weights/np.sum(bkg_weights),
        alpha=0.4,
        label=f"{label_bkg} $\\cos(\\theta_{{\\mathrm{{beam}}}})_{{\\mathrm{{max}}}}$",
        color='red',
    )
    plt.xlabel('$\\cos(\\theta_{\\mathrm{beam}})$')
    plt.ylabel('Counts')
    plt.title(f"{label_bkg} $\\cos(\\theta_{{\\mathrm{{beam}}}})$ Hits, PFP and Max")
    plt.legend()
    plt.savefig(os.path.join(output_folder_base, "bkg_direction_histograms.png"))
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.hist(
        sig_events[:, aggregate_dict["directionZ"]],
        bins=nbins,
        range=(-1, 1),
        weights=sig_weights/np.sum(sig_weights),
        alpha=0.4,
        label=f"{label_sig} $\\cos(\\theta_{{\\mathrm{{beam}}}})$ Hits",
        color='blue',
    )
    plt.hist(
        sig_events[:, aggregate_dict["directionZ2"]],
        bins=nbins,
        range=(-1, 1),
        weights=sig_weights/np.sum(sig_weights),
        alpha=0.4,
        label=f"{label_sig} $\\cos(\\theta_{{\\mathrm{{beam}}}})$ PFP",
        color='green',
    )
    plt.hist(
        np.maximum(
            sig_events[:, aggregate_dict["directionZ"]],
            sig_events[:, aggregate_dict["directionZ2"]]
        ),
        bins=nbins,
        range=(-1, 1),
        weights=sig_weights/np.sum(sig_weights),
        alpha=0.4,
        label=f"{label_sig} $\\cos(\\theta_{{\\mathrm{{beam}}}})_{{\\mathrm{{max}}}}$",
        color='red',
    )
    plt.xlabel('$\\cos(\\theta_{\\mathrm{beam}})$')
    plt.ylabel('Counts')
    plt.title(f"{label_sig} $\\cos(\\theta_{{\\mathrm{{beam}}}})_{{\\mathrm{{max}}}}$ Hits, PFP and Max")
    plt.legend()
    plt.savefig(os.path.join(output_folder_base, "sig_direction_histograms.png"))
    plt.close()

    # ----------------------------------------------
    # # Step 2: Compare the energy deposited in the first 10 cm between {label_sig} off spill and {label_bkg}
    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut="energy_first10cm"
        )

    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInFirst10cm"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInFirst10cm"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_energyDepositedInFirst10cm_spectrum",
        title=f"{label_bkg} vs {label_sig}." + "  $C_{0-10}$",
        range=(100, 100000),
        bins=nbins,
        xlabel="$C_{0-10}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )


    # ----------------------------------------------
    # Step 6: Compare the energy deposited in the fifth 10 cm between {label_sig} off spill and {label_bkg}
    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut="energy_fifth10cm"
        )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInFifth10cm"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInFifth10cm"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_energyDepositedInFifth10cm_spectrum",
        title=f"{label_bkg} vs {label_sig}. " + "  $C_{40-50}$",
        range=(100, 100000),
        bins=nbins,
        xlabel="$C_{40-50}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # ----------------------------------------------
    # Step 7: Compare the energy deposited in the sixth 10 cm between {label_sig} off spill and {label_bkg}
    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut="energy_sixth10cm"
        )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInSixth10cm"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInSixth10cm"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_energyDepositedInSixth10cm_spectrum",
        title=f"{label_bkg} vs {label_sig}. " + "  $C_{50-60}$",
        range=(100, 100000),
        bins=nbins,
        xlabel="$C_{50-60}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    # ----------------------------------------------
    # Step 9: Compare the energy deposited in the eighth 10 cm between {label_sig} off spill and {label_bkg}
    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut="energy_eighth10cm"
        )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInEighth10cm"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInEighth10cm"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_energyDepositedInEighth10cm_spectrum",
        title=f"{label_bkg} vs {label_sig}. " + "  $C_{70-80}$",
        range=(100, 100000),
        bins=nbins,
        xlabel="$C_{70-80}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    # ----------------------------------------------
    # Step 12: Compare the energy deposited in the eleventh 10 cm between {label_sig} off spill and {label_bkg}
    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut="energy_eleventh10cm"
        )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInEleventh10cm"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInEleventh10cm"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_energyDepositedInEleventh10cm_spectrum",
        title=f"{label_bkg} vs {label_sig}. " + "  $C_{100-110}$",
        range=(100, 100000),
        bins=nbins,
        xlabel="$C_{100-110}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    # ----------------------------------------------
    # Step 14: Compare the energy deposited in the thirteenth 10 cm between {label_sig} off spill and {label_bkg}
    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
        cuts=cuts, skip_cut="energy_thirteenth10cm"
        )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInThirteenth10cm"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInThirteenth10cm"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_energyDepositedInThirteenth10cm_spectrum",
        title=f"{label_bkg} vs {label_sig}. " + "  $C_{120-130}$",
        range=(100, 100000),
        bins=nbins,
        xlabel="$C_{120-130}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    # ----------------------------------------------
    # Step 16: Compare the energy deposited in the fifteenth 10 cm between {label_sig} off spill and {label_bkg}
    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut="energy_fifteenth10cm"
        )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInFifteenth10cm"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInFifteenth10cm"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_energyDepositedInFifteenth10cm_spectrum",
        title=f"{label_bkg} vs {label_sig}. " + "  $C_{140-150}$",
        range=(100, 100000),
        bins=nbins,
        xlabel="$C_{140-150}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    # ----------------------------------------------
    # Step 19: Compare ratio of energy deposited in the first 10 cm before the event to the energy deposited in the fifth 10 cm between {label_sig} off spill and {label_bkg}
    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut="ratio_before_after"
        )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInFirst10cmBefore"]]/sig_events[:, aggregate_dict["energyDepositedInFifth10cm"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInFirst10cmBefore"]]/bkg_events[:, aggregate_dict["energyDepositedInFifth10cm"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_energyDepositedInFirst10cmBefore_spectrum_ratio",
        title=f"{label_bkg} vs {label_sig}. " + "  $C_{-10-0}/C_{40-50}$",
        range=(0, 20),
        bins=nbins,
        xlabel="$C_{-10-0}/C_{40-50}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # ----------------------------------------------
    # Step 20: Compare the sigma of the time fit between {label_sig} off spill and {label_bkg}
    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut="timeFitSigma"
        )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["timeFitSigma"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["timeFitSigma"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_timeFitSigma_spectrum",
        title=f"{label_bkg} vs {label_sig}. Time Fit Sigma Spectrum",
        range=(0, 100),
        bins=nbins,
        xlabel="Time Fit Sigma (ns)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    # ----------------------------------------------
    # Step 20: Compare the ROI size between {label_sig} off spill and {label_bkg}
    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut="ROI_Z_size"
        )
    triple_hist(
        hist_sig=(sig_events[:, aggregate_dict["zROIEnd"]]-sig_events[:, aggregate_dict["zROIStart"]])*460/100,
        sig_weights=sig_weights,
        hist_bkg=(bkg_events[:, aggregate_dict["zROIEnd"]]-bkg_events[:, aggregate_dict["zROIStart"]])*460/100,
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_ROI_Z_size_spectrum",
        title=f"{label_bkg} vs {label_sig}. ROI Size Spectrum",
        range=(0, 460),
        bins=nbins,
        xlabel="ROI Size Z",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    # ----------------------------------------------
    # Step 22: compare the ROI starting position
    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut="ROI_Z_starting_point"
        )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["zROIStart"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["zROIStart"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_ROI_Z_start_spectrum",
        title=f"{label_bkg} vs {label_sig}. ROI Z Start Spectrum",
        range=(0, 100),
        bins=nbins,
        xlabel="ROI Z Start (cm)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

        # ----------------------------------------------
    # Step 22: compare the ROI starting position
    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut="ROI_Z_starting_point_close_to_vertexZ"
        )
    triple_hist(
        hist_sig=np.abs(sig_events[:, aggregate_dict["vertexZ"]] - sig_events[:, aggregate_dict["zROIStart"]]*460/100),
        sig_weights=sig_weights,
        hist_bkg=np.abs(bkg_events[:, aggregate_dict["vertexZ"]] - bkg_events[:, aggregate_dict["zROIStart"]]*460/100),
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_ROI_Z_start_spectrum",
        title=f"{label_bkg} vs {label_sig}. ROI Z Start Spectrum",
        range=(0, 100),
        bins=nbins,
        xlabel="ROI Z Start (cm)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    
    # ----------------------------------------------
    # Step 23: compare the lenght of neutrino tail and compare the density of neutrino tail

    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut="Neutrino_Tail_Length_Density"
        )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["numberOfHitsInMuonRegion"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["numberOfHitsInMuonRegion"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_length_of_neutrino_tail",
        title=f"{label_bkg} vs {label_sig}. Length of neutrino tail",
        range=(-2, 460),
        bins=nbins,
        xlabel="Length of neutrino tail (cm)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["numberOfHitsInMuonRegion"]]/(sig_events[:, aggregate_dict["zROIStart"]]*460/100),
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["numberOfHitsInMuonRegion"]]/(bkg_events[:, aggregate_dict["zROIStart"]]*460/100),
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_tail_density",
        title=f"{label_bkg} vs {label_sig}. Density of neutrino tail",
        range=(0,1),
        bins=nbins,
        xlabel="Density ",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # ----------------------------------------------
    # Step 23: compare vertex position Z vs ROI starting position Z
    if apply_cuts:

        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut=None
        )   
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["vertexX"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["vertexX"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_vertexX_spectrum",
        title=f"{label_bkg} vs {label_sig}. Vertex X Spectrum",
        range=(0, 600),
        bins=nbins,
        xlabel="Vertex X (cm)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["vertexY"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["vertexY"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_vertexY_spectrum",
        title=f"{label_bkg} vs {label_sig}. Vertex Y Spectrum",
        range=(0, 600),
        bins=nbins,
        xlabel="Vertex Y (cm)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["vertexZ"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["vertexZ"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_vertexZ_spectrum",
        title=f"{label_bkg} vs {label_sig}. Vertex Z Spectrum",
        range=(0, 460),
        bins=nbins,
        xlabel="Vertex Z (cm)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # PFParticles number
    if apply_cuts:
        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut="daughter_particles"
        )  

    triple_hist(
        hist_sig=sig_events[:, aggregate_dict["numberOfPFParticles"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["numberOfPFParticles"]],
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_numberOfPFParticles_spectrum",
        title=f"{label_bkg} vs {label_sig}. Number of PFParticles Spectrum",
        range=(0, 25),
        bins=25,
        xlabel="Number of PFParticles",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # ROI Z start vs vertex Z
    if apply_cuts:
        sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames, bkg_events, bkg_weights, bkg_labels, bkg_filenames = apply_all_cuts(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            cuts=cuts, skip_cut="ROI_Z_starting_point_close_to_vertexZ"
        )  

    triple_hist(
        hist_sig=np.abs(sig_events[:, aggregate_dict["vertexZ"]] - sig_events[:, aggregate_dict["zROIStart"]]*460/100),
        sig_weights=sig_weights,
        hist_bkg=np.abs(bkg_events[:, aggregate_dict["vertexZ"]] - bkg_events[:, aggregate_dict["zROIStart"]]*460/100),
        bkg_weights=bkg_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_distance_vertexZ_to_ROI_start",
        title=f"{label_bkg} vs {label_sig}. Distance Vertex Z to ROI Start",
        range=(0, 200),
        bins=nbins,
        xlabel="Distance from Vertex Z to ROI Start (cm)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

def run_all_triple_hist_stack(
    sig_events, sig_weights, bkg_events, bkg_weights, mc_events, mc_weights,
    output_folder_base, label_sig, label_bkg, nbins, spill_status=f"spill_status",
    apply_cuts=False, cuts=cuts, total_time_with_correct_tps_on=(2.93733+5.5266)
):
    sig_events_orig = sig_events.copy()
    sig_weights_orig = sig_weights.copy()
    sig_true_events_orig = np.ones((sig_events.shape[0], 1)) * -1
    sig_labels_orig = np.zeros((sig_events.shape[0], 1))
    sig_filenames_orig = np.zeros((sig_events.shape[0], 1))

    bkg_events_orig = bkg_events.copy()
    bkg_weights_orig = bkg_weights.copy()
    bkg_labels_orig = np.zeros((bkg_events.shape[0], 1))
    bkg_filenames_orig = np.zeros((bkg_events.shape[0], 1))

    mc_events_orig = mc_events.copy()
    mc_weights_orig = mc_weights.copy()
    mc_labels_orig = np.zeros((mc_events.shape[0], 1))
    mc_filenames_orig = np.zeros((mc_events.shape[0], 1))
    mc_true_events_orig = np.ones((mc_events.shape[0], 1))

    # Step 1: Energy spectrum
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig, cuts=cuts, skip_cut=None
        )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["reconstructedEnergy"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["reconstructedEnergy"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["reconstructedEnergy"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_energy_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC. Energy Spectrum",
        range=(0, 100),
        bins=nbins,
        xlabel="Energy (GeV)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # Step 2: Direction spectra (skip max_directionZ)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut="max_directionZ"
        )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["directionX"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["directionX"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["directionX"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_directionX_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC. Direction X Spectrum",
        range=(-1, 1),
        bins=nbins,
        xlabel="Direction X",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["directionY"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["directionY"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["directionY"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_directionY_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC. Direction Y Spectrum",
        range=(-1, 1),
        bins=nbins,
        xlabel="Direction Y",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["directionZ"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["directionZ"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["directionZ"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_directionZ_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC. Direction Z Spectrum",
        range=(-1, 1),
        bins=nbins,
        xlabel="Direction Z",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # Step 2.5: Direction X2/Y2/Z2
    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["directionX2"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["directionX2"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["directionX2"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_directionX2_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC. Direction X2 Spectrum",
        range=(-1, 1),
        bins=nbins,
        xlabel="Direction X2",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["directionY2"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["directionY2"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["directionY2"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_directionY2_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC. Direction Y2 Spectrum",
        range=(-1, 1),
        bins=nbins,
        xlabel="Direction Y2",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["directionZ2"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["directionZ2"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["directionZ2"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_directionZ2_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC. Direction Z2 Spectrum",
        range=(-1, 1),
        bins=nbins,
        xlabel="Direction Z2",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # Step 2.6: Max(directionZ, directionZ2)
    triple_hist_stack(
        hist_sig=np.maximum(
            sig_events[:, aggregate_dict["directionZ"]],
            sig_events[:, aggregate_dict["directionZ2"]],
        ),
        sig_weights=sig_weights,
        hist_bkg=np.maximum(
            bkg_events[:, aggregate_dict["directionZ"]],
            bkg_events[:, aggregate_dict["directionZ2"]],
        ),
        bkg_weights=bkg_weights,
        hist_mc=np.maximum(
            mc_events[:, aggregate_dict["directionZ"]],
            mc_events[:, aggregate_dict["directionZ2"]],
        ),
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_max_directionZ_directionZ2_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC."+" $\\cos(\\theta_{{\\mathrm{{beam}}}})_{{\\mathrm{{max}}}}$ Spectrum",
        range=(-1, 1),
        bins=nbins,
        xlabel="$\\cos(\\theta_{{\\mathrm{{beam}}}})_{{\\mathrm{{max}}}}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # Step: Energy in first 10 cm (skip energy_first10cm)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut="energy_first10cm"
        )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInFirst10cm"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInFirst10cm"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["energyDepositedInFirst10cm"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_energyDepositedInFirst10cm_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC." + "  $C_{0-10}$",
        range=(100, 100000),
        bins=nbins,
        xlabel="$C_{0-10}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # Fifth 10 cm (skip energy_fifth10cm)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut="energy_fifth10cm"
        )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInFifth10cm"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInFifth10cm"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["energyDepositedInFifth10cm"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_energyDepositedInFifth10cm_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC." + "  $C_{40-50}$",
        range=(100, 100000),
        bins=nbins,
        xlabel="$C_{40-50}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # Sixth 10 cm (skip energy_sixth10cm)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut="energy_sixth10cm"
        )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInSixth10cm"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInSixth10cm"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["energyDepositedInSixth10cm"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_energyDepositedInSixth10cm_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC." + "  $C_{50-60}$",
        range=(100, 100000),
        bins=nbins,
        xlabel="$C_{50-60}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # Eighth 10 cm (skip energy_eighth10cm)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut="energy_eighth10cm"
        )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInEighth10cm"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInEighth10cm"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["energyDepositedInEighth10cm"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_energyDepositedInEighth10cm_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC." + "  $C_{70-80}$",
        range=(100, 100000),
        bins=nbins,
        xlabel="$C_{70-80}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # Eleventh 10 cm (skip energy_eleventh10cm)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut="energy_eleventh10cm"
        )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInEleventh10cm"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInEleventh10cm"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["energyDepositedInEleventh10cm"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_energyDepositedInEleventh10cm_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC." + " $C_{100-110}$",
        range=(100, 100000),
        bins=nbins,
        xlabel="$C_{100-110}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # Thirteenth 10 cm (skip energy_thirteenth10cm)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut="energy_thirteenth10cm"
        )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInThirteenth10cm"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInThirteenth10cm"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["energyDepositedInThirteenth10cm"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_energyDepositedInThirteenth10cm_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC." + " $C_{120-130}$",
        range=(100, 100000),
        bins=nbins,
        xlabel="$C_{120-130}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # Fifteenth 10 cm (skip energy_fifteenth10cm)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut="energy_fifteenth10cm"
        )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInFifteenth10cm"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInFifteenth10cm"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["energyDepositedInFifteenth10cm"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_energyDepositedInFifteenth10cm_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC." + " $C_{140-150}$",
        range=(100, 100000),
        bins=nbins,
        xlabel="$C_{140-150}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # Ratio before/after (skip ratio_before_after)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut="ratio_before_after"
        )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["energyDepositedInFirst10cmBefore"]] / np.maximum(
            sig_events[:, aggregate_dict["energyDepositedInFifth10cm"]], 1e-12),
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["energyDepositedInFirst10cmBefore"]] / np.maximum(
            bkg_events[:, aggregate_dict["energyDepositedInFifth10cm"]], 1e-12),
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["energyDepositedInFirst10cmBefore"]] / np.maximum(
            mc_events[:, aggregate_dict["energyDepositedInFifth10cm"]], 1e-12),
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_energyDepositedInFirst10cmBefore_ratio",
        title=f"{label_sig} vs {label_bkg} + MC. Ratio $C_{0-10} / C_{40-50}$",
        range=(0, 20),
        bins=nbins,
        xlabel="$C_{0-10} / C_{40-50}$",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # Time fit sigma (skip timeFitSigma)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut="timeFitSigma"
        )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["timeFitSigma"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["timeFitSigma"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["timeFitSigma"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_timeFitSigma_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC. Time Fit Sigma",
        range=(0, 100),
        bins=nbins,
        xlabel="Time Fit Sigma (ns)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # ROI Z size (skip ROI_Z_size)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut="ROI_Z_size"
        )

    triple_hist_stack(
        hist_sig=(sig_events[:, aggregate_dict["zROIEnd"]] - sig_events[:, aggregate_dict["zROIStart"]]) * 460 / 100,
        sig_weights=sig_weights,
        hist_bkg=(bkg_events[:, aggregate_dict["zROIEnd"]] - bkg_events[:, aggregate_dict["zROIStart"]]) * 460 / 100,
        bkg_weights=bkg_weights,
        hist_mc=(mc_events[:, aggregate_dict["zROIEnd"]] - mc_events[:, aggregate_dict["zROIStart"]]) * 460 / 100,
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_ROI_Z_size_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC. ROI Size Z",
        range=(0, 460),
        bins=nbins,
        xlabel="ROI Size Z (cm)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # ROI Z start (skip ROI_Z_starting_point)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut="ROI_Z_starting_point"
        )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["zROIStart"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["zROIStart"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["zROIStart"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_ROI_Z_start_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC. ROI Z Start",
        range=(0, 100),
        bins=nbins,
        xlabel="ROI Z Start (cm)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # Distance vertexZ to ROI start (skip ROI_Z_starting_point_close_to_vertexZ)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut="ROI_Z_starting_point_close_to_vertexZ"
        )

    triple_hist_stack(
        hist_sig=np.abs(sig_events[:, aggregate_dict["vertexZ"]] - sig_events[:, aggregate_dict["zROIStart"]] * 460 / 100),
        sig_weights=sig_weights,
        hist_bkg=np.abs(bkg_events[:, aggregate_dict["vertexZ"]] - bkg_events[:, aggregate_dict["zROIStart"]] * 460 / 100),
        bkg_weights=bkg_weights,
        hist_mc=np.abs(mc_events[:, aggregate_dict["vertexZ"]] - mc_events[:, aggregate_dict["zROIStart"]] * 460 / 100),
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_distance_vertexZ_to_ROI_start",
        title=f"{label_sig} vs {label_bkg} + MC. Distance Vertex Z to ROI Start",
        range=(0, 200),
        bins=nbins,
        xlabel="Distance Vertex Z to ROI Start (cm)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # Neutrino tail length and density (skip Neutrino_Tail_Length_Density)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut="Neutrino_Tail_Length_Density"
        )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["numberOfHitsInMuonRegion"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["numberOfHitsInMuonRegion"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["numberOfHitsInMuonRegion"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_length_of_neutrino_tail",
        title=f"{label_sig} vs {label_bkg} + MC. Length of neutrino tail",
        range=(-2, 460),
        bins=nbins,
        xlabel="Length of neutrino tail (cm)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["numberOfHitsInMuonRegion"]] / np.maximum(sig_events[:, aggregate_dict["zROIStart"]] * 460 / 100, 1e-12),
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["numberOfHitsInMuonRegion"]] / np.maximum(bkg_events[:, aggregate_dict["zROIStart"]] * 460 / 100, 1e-12),
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["numberOfHitsInMuonRegion"]] / np.maximum(mc_events[:, aggregate_dict["zROIStart"]] * 460 / 100, 1e-12),
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_tail_density",
        title=f"{label_sig} vs {label_bkg} + MC. Density of neutrino tail",
        range=(0, 1),
        bins=nbins,
        xlabel="Density",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # Vertex positions (no skip)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut=None
        )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["vertexX"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["vertexX"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["vertexX"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_vertexX_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC. Vertex X Spectrum",
        range=(-400, 400),
        bins=nbins,
        xlabel="Vertex X (cm)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["vertexY"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["vertexY"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["vertexY"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_vertexY_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC. Vertex Y Spectrum",
        range=(-1, 600),
        bins=nbins,
        xlabel="Vertex Y (cm)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )
    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["vertexZ"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["vertexZ"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["vertexZ"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_vertexZ_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC. Vertex Z Spectrum",
        range=(0, 460),
        bins=nbins,
        xlabel="Vertex Z (cm)",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    # make a scatter plot of vertexZ vs vertexY
    plt.figure(figsize=(8, 6))
    plt.plot(
        sig_events[:, aggregate_dict["vertexZ"]],
        sig_events[:, aggregate_dict["vertexY"]],
        'o', label=label_sig
    )
    plt.plot(
        bkg_events[:, aggregate_dict["vertexZ"]],
        bkg_events[:, aggregate_dict["vertexY"]],
        'o', label=label_bkg
    )
    # set axis range, x from 0 to 460, y from 0 to 600
    plt.xlim(0,460)
    plt.ylim(0,600)

    plt.xlabel("Vertex Z (cm)")
    plt.ylabel("Vertex Y (cm)")
    plt.title(f"{label_sig} vs {label_bkg}. Vertex Z vs Vertex Y")
    plt.legend()
    plt.savefig(f"{output_folder_base}/sig_bkg_vertexZ_vs_vertexY_scatter.png")
    plt.close()

    # Number of PFParticles (skip daughter_particles)
    if apply_cuts:
        (sig_events, sig_weights, sig_true_events, sig_labels, sig_filenames,
         bkg_events, bkg_weights, bkg_labels, bkg_filenames,
         mc_events, mc_weights, mc_true_events, mc_labels, mc_filenames) = apply_all_cuts_MC(
            sig_events_orig, sig_weights_orig, sig_true_events_orig, sig_labels_orig, sig_filenames_orig,
            bkg_events_orig, bkg_weights_orig, bkg_labels_orig, bkg_filenames_orig,
            mc_events_orig, mc_weights_orig, mc_true_events_orig, mc_labels_orig, mc_filenames_orig,
            cuts=cuts, skip_cut="daughter_particles"
        )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["numberOfPFParticles"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["numberOfPFParticles"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["numberOfPFParticles"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,
        plot_name="sig_bkg_mc_numberOfPFParticles_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC. Number of PFParticles",
        range=(0, 25),
        bins=25,
        xlabel="Number of PFParticles",
        ylabel="Counts",
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )

    triple_hist_stack(
        hist_sig=sig_events[:, aggregate_dict["totalNumberOfHits"]],
        sig_weights=sig_weights,
        hist_bkg=bkg_events[:, aggregate_dict["totalNumberOfHits"]],
        bkg_weights=bkg_weights,
        hist_mc=mc_events[:, aggregate_dict["totalNumberOfHits"]],
        mc_weights=mc_weights,
        
        output_folder=output_folder_base,
        spill_status=spill_status,  
        plot_name="sig_bkg_mc_totalNumberOfHits_spectrum",
        title=f"{label_sig} vs {label_bkg} + MC. Total Number of Hits",
        range=(0, 200000),
        bins=nbins,
        xlabel="Total Number of Hits",
        ylabel="Counts",    
        label_sig=label_sig,
        label_bkg=label_bkg,
        total_time_with_correct_tps_on=total_time_with_correct_tps_on
    )




