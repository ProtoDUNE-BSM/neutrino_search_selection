import numpy as np
from name_index_association import *


# Cut position thresholds dictionary
cut_thresholds = {
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


# List of cuts as (name, function)
cuts = [
    ("vertex_fiducial_volume", lambda s, b: (
        np.where(
            (s[:, aggregate_dict["vertexZ"]] >= cut_thresholds["vertex_z_min"]) &
            (s[:, aggregate_dict["vertexY"]] <= cut_thresholds["vertex_y_max"])
        )[0],
        np.where(
            (b[:, aggregate_dict["vertexZ"]] >= cut_thresholds["vertex_z_min"]) &
            (b[:, aggregate_dict["vertexY"]] <= cut_thresholds["vertex_y_max"])
        )[0]
    )),
    ("daughter_particles", lambda s, b: (
        np.where(s[:, aggregate_dict["numberOfPFParticles"]] >= cut_thresholds["num_pf_particles_min"])[0],
        np.where(b[:, aggregate_dict["numberOfPFParticles"]] >= cut_thresholds["num_pf_particles_min"])[0]
    )),
    ("max_directionZ", lambda s, b: (
        np.where(np.maximum(
            s[:, aggregate_dict["directionZ"]],
            s[:, aggregate_dict["directionZ2"]]
        ) > cut_thresholds["direction_z_max"])[0],
        np.where(np.maximum(
            b[:, aggregate_dict["directionZ"]],
            b[:, aggregate_dict["directionZ2"]]
        ) > cut_thresholds["direction_z_max"])[0]
    )),
    ("energy_first10cm", lambda s, b: (
        np.where(
            (s[:, aggregate_dict["energyDepositedInFirst10cm"]] > cut_thresholds["energy_first10cm_min"])
        )[0],
        np.where(
            (b[:, aggregate_dict["energyDepositedInFirst10cm"]] > cut_thresholds["energy_first10cm_min"])
        )[0]
    )),
    ("energy_fifth10cm", lambda s, b: (
        np.where(
            (s[:, aggregate_dict["energyDepositedInFifth10cm"]] > cut_thresholds["energy_fifth10cm_min"])
        )[0],
        np.where(
            (b[:, aggregate_dict["energyDepositedInFifth10cm"]] > cut_thresholds["energy_fifth10cm_min"])
        )[0]
    )),
    ("energy_fifteenth10cm", lambda s, b: (
        np.where(
            (s[:, aggregate_dict["energyDepositedInFifteenth10cm"]] > cut_thresholds["energy_fifteenth10cm_min"])
        )[0],
        np.where(
            (b[:, aggregate_dict["energyDepositedInFifteenth10cm"]] > cut_thresholds["energy_fifteenth10cm_min"])
        )[0]
    )),    
    ("ROI_Z_size", lambda s, b: (
        np.where(s[:, aggregate_dict["zROIEnd"]] - s[:, aggregate_dict["zROIStart"]] > cut_thresholds["roi_z_size_min"])[0],
        np.where(b[:, aggregate_dict["zROIEnd"]] - b[:, aggregate_dict["zROIStart"]] > cut_thresholds["roi_z_size_min"])[0]
    )),
    ("ROI_Z_starting_point_close_to_vertexZ", lambda s, b: (
        np.where(np.abs(s[:, aggregate_dict["vertexZ"]] - s[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"]) < cut_thresholds["roi_z_vertex_distance_max"])[0],
        np.where(np.abs(b[:, aggregate_dict["vertexZ"]] - b[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"]) < cut_thresholds["roi_z_vertex_distance_max"])[0]
    )),
    ("Neutrino_Tail_Length_Density", lambda s, b: (
        np.where(s[:, aggregate_dict["numberOfHitsInMuonRegion"]]/(s[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"]) < cut_thresholds["tail_length_density_max"])[0],
        np.where(b[:, aggregate_dict["numberOfHitsInMuonRegion"]]/(b[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"]) < cut_thresholds["tail_length_density_max"])[0]
    ))
]

cuts_with_mc = [
    ("vertex_fiducial_volume", lambda s, b, mc: (
        np.where(
            (s[:, aggregate_dict["vertexZ"]] >= cut_thresholds["vertex_z_min"]) &
            (s[:, aggregate_dict["vertexY"]] <= cut_thresholds["vertex_y_max"])
        )[0],
        np.where(
            (b[:, aggregate_dict["vertexZ"]] >= cut_thresholds["vertex_z_min"]) &
            (b[:, aggregate_dict["vertexY"]] <= cut_thresholds["vertex_y_max"])
        )[0],
        np.where(
            (mc[:, aggregate_dict["vertexZ"]] >= cut_thresholds["vertex_z_min"]) &
            (mc[:, aggregate_dict["vertexY"]] <= cut_thresholds["vertex_y_max"])
        )[0]
    )),
    ("daughter_particles", lambda s, b, mc: (
        np.where(s[:, aggregate_dict["numberOfPFParticles"]] >= cut_thresholds["num_pf_particles_min"])[0],
        np.where(b[:, aggregate_dict["numberOfPFParticles"]] >= cut_thresholds["num_pf_particles_min"])[0],
        np.where(mc[:, aggregate_dict["numberOfPFParticles"]] >= cut_thresholds["num_pf_particles_min"])[0]
    )),
    ("max_directionZ", lambda s, b, mc: (
        np.where(np.maximum(
            s[:, aggregate_dict["directionZ"]],
            s[:, aggregate_dict["directionZ2"]]
        ) > cut_thresholds["direction_z_max"])[0],
        np.where(np.maximum(
            b[:, aggregate_dict["directionZ"]],
            b[:, aggregate_dict["directionZ2"]]
        ) > cut_thresholds["direction_z_max"])[0],
        np.where(np.maximum(
            mc[:, aggregate_dict["directionZ"]],
            mc[:, aggregate_dict["directionZ2"]]
        ) > cut_thresholds["direction_z_max"])[0]
    )),
    ("energy_first10cm", lambda s, b, mc: (
        np.where(
            (s[:, aggregate_dict["energyDepositedInFirst10cm"]] > cut_thresholds["energy_first10cm_min"])
        )[0],
        np.where(
            (b[:, aggregate_dict["energyDepositedInFirst10cm"]] > cut_thresholds["energy_first10cm_min"])
        )[0],
        np.where(
            (mc[:, aggregate_dict["energyDepositedInFirst10cm"]] > cut_thresholds["energy_first10cm_min"])
        )[0]
    )),
    ("energy_fifth10cm", lambda s, b, mc: (
        np.where(
            (s[:, aggregate_dict["energyDepositedInFifth10cm"]] > cut_thresholds["energy_fifth10cm_min"])
        )[0],
        np.where(
            (b[:, aggregate_dict["energyDepositedInFifth10cm"]] > cut_thresholds["energy_fifth10cm_min"])
        )[0],
        np.where(
            (mc[:, aggregate_dict["energyDepositedInFifth10cm"]] > cut_thresholds["energy_fifth10cm_min"])
        )[0]
    )),
    ("energy_fifteenth10cm", lambda s, b, mc: (
        np.where(
            (s[:, aggregate_dict["energyDepositedInFifteenth10cm"]] > cut_thresholds["energy_fifteenth10cm_min"])
        )[0],
        np.where(
            (b[:, aggregate_dict["energyDepositedInFifteenth10cm"]] > cut_thresholds["energy_fifteenth10cm_min"])
        )[0],
        np.where(
            (mc[:, aggregate_dict["energyDepositedInFifteenth10cm"]] > cut_thresholds["energy_fifteenth10cm_min"])
        )[0]
    )),
    ("ROI_Z_size", lambda s, b, mc: (
        np.where(s[:, aggregate_dict["zROIEnd"]] - s[:, aggregate_dict["zROIStart"]] > cut_thresholds["roi_z_size_min"])[0],
        np.where(b[:, aggregate_dict["zROIEnd"]] - b[:, aggregate_dict["zROIStart"]] > cut_thresholds["roi_z_size_min"])[0],
        np.where(mc[:, aggregate_dict["zROIEnd"]] - mc[:, aggregate_dict["zROIStart"]] > cut_thresholds["roi_z_size_min"])[0]
    )),
    ("ROI_Z_starting_point_close_to_vertexZ", lambda s, b, mc: (
        np.where(np.abs(s[:, aggregate_dict["vertexZ"]] - s[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"]) < cut_thresholds["roi_z_vertex_distance_max"])[0],
        np.where(np.abs(b[:, aggregate_dict["vertexZ"]] - b[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"]) < cut_thresholds["roi_z_vertex_distance_max"])[0],
        np.where(np.abs(mc[:, aggregate_dict["vertexZ"]] - mc[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"]) < cut_thresholds["roi_z_vertex_distance_max"])[0]
    )),
    ("Neutrino_Tail_Length_Density", lambda s, b, mc: (
        np.where(s[:, aggregate_dict["numberOfHitsInMuonRegion"]]/(s[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"]) < cut_thresholds["tail_length_density_max"])[0],
        np.where(b[:, aggregate_dict["numberOfHitsInMuonRegion"]]/(b[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"]) < cut_thresholds["tail_length_density_max"])[0],
        np.where(mc[:, aggregate_dict["numberOfHitsInMuonRegion"]]/(mc[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"]) < cut_thresholds["tail_length_density_max"])[0]
    ))
]
cuts_single = [
    ("vertex_fiducial_volume", lambda s: (
        np.where(
            (s[:, aggregate_dict["vertexZ"]] >= cut_thresholds["vertex_z_min"]) &
            (s[:, aggregate_dict["vertexY"]] <= cut_thresholds["vertex_y_max"])
        )[0]
    )),
    ("daughter_particles", lambda s: (
        np.where(s[:, aggregate_dict["numberOfPFParticles"]] >= cut_thresholds["num_pf_particles_min"])[0]
    )),
    ("max_directionZ", lambda s: (
        np.where(np.maximum(
            s[:, aggregate_dict["directionZ"]],
            s[:, aggregate_dict["directionZ2"]]
        ) > cut_thresholds["direction_z_max"])[0]
    )),
    ("energy_first10cm", lambda s: (
        np.where(
            (s[:, aggregate_dict["energyDepositedInFirst10cm"]] > cut_thresholds["energy_first10cm_min"])
        )[0]
    )),
    ("energy_fifth10cm", lambda s: (
        np.where(
            (s[:, aggregate_dict["energyDepositedInFifth10cm"]] > cut_thresholds["energy_fifth10cm_min"])
        )[0]
    )),
    ("energy_fifteenth10cm", lambda s: (
        np.where(
            (s[:, aggregate_dict["energyDepositedInFifteenth10cm"]] > cut_thresholds["energy_fifteenth10cm_min"])
        )[0]
    )),    
    ("ROI_Z_size", lambda s: (
        np.where(s[:, aggregate_dict["zROIEnd"]] - s[:, aggregate_dict["zROIStart"]] > cut_thresholds["roi_z_size_min"])[0]
    )),
    ("ROI_Z_starting_point_close_to_vertexZ", lambda s: (
        np.where(np.abs(s[:, aggregate_dict["vertexZ"]] - s[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"]) < cut_thresholds["roi_z_vertex_distance_max"])[0]
    )),
    ("Neutrino_Tail_Length_Density", lambda s: (
        np.where(s[:, aggregate_dict["numberOfHitsInMuonRegion"]]/(s[:, aggregate_dict["zROIStart"]]*cut_thresholds["roi_z_scale_factor"]) < cut_thresholds["tail_length_density_max"])[0]
    ))
]



cuts_names = [
    "vertex_fiducial_volume",
    "daughter_particles",
    "max_directionZ",
    "energy_first10cm",
    "energy_fifth10cm",
    "energy_fifteenth10cm",
    "ROI_Z_size",
    "ROI_Z_starting_point_close_to_vertexZ",
    "Neutrino_Tail_Length_Density"
]
shorter_cut_names = [
    'Vertex fiducial volume',
    'N daughter particles',
    'Max $\\cos(\\theta_{{\\mathrm{{beam}}}})$',
    'Cut on $C_{0-10}$',
    'Cut on $C_{40-50}$',
    'Cut on $C_{140-150}$',
    'ROI Z size',
    'ROI Z close to vertex',
    'Tail Length Density',
]


