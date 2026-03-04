import numpy as np
from name_index_association import *


# List of cuts as (name, function)
cuts = [
    ("vertex_fiducial_volume", lambda s, b: (
        np.where(
            (s[:, aggregate_dict["vertexZ"]] >= 20) &
            (s[:, aggregate_dict["vertexY"]] <= 550)
        )[0],
        np.where(
            (b[:, aggregate_dict["vertexZ"]] >= 20) &
            (b[:, aggregate_dict["vertexY"]] <= 550)
        )[0]
    )),
    ("daughter_particles", lambda s, b: (
        np.where(s[:, aggregate_dict["numberOfPFParticles"]] >= 6)[0],
        np.where(b[:, aggregate_dict["numberOfPFParticles"]] >= 6)[0]
    )),
    ("max_directionZ", lambda s, b: (
        np.where(np.maximum(
            s[:, aggregate_dict["directionZ"]],
            s[:, aggregate_dict["directionZ2"]]
        ) > 0.97)[0],
        np.where(np.maximum(
            b[:, aggregate_dict["directionZ"]],
            b[:, aggregate_dict["directionZ2"]]
        ) > 0.97)[0]
    )),
    ("energy_first10cm", lambda s, b: (
        np.where(
            (s[:, aggregate_dict["energyDepositedInFirst10cm"]] > 3000))[0],
        np.where(
            (b[:, aggregate_dict["energyDepositedInFirst10cm"]] > 3000))[0]
    )),    
    ("energy_fifth10cm", lambda s, b: (
        np.where(
            (s[:, aggregate_dict["energyDepositedInFifth10cm"]] > 10000))[0],
        np.where(
            (b[:, aggregate_dict["energyDepositedInFifth10cm"]] > 10000))[0]
    )),
    ("energy_fifteenth10cm", lambda s, b: (
        np.where(
            (s[:, aggregate_dict["energyDepositedInFifteenth10cm"]] > 3000))[0],
        np.where(
            (b[:, aggregate_dict["energyDepositedInFifteenth10cm"]] > 3000))[0]
    )),    
    ("ROI_Z_size", lambda s, b: (
        np.where(s[:, aggregate_dict["zROIEnd"]] - s[:, aggregate_dict["zROIStart"]] > 15)[0],
        np.where(b[:, aggregate_dict["zROIEnd"]] - b[:, aggregate_dict["zROIStart"]] > 15)[0]
    )),
    ("ROI_Z_starting_point_close_to_vertexZ", lambda s, b: (
        np.where(np.abs(s[:, aggregate_dict["vertexZ"]] - s[:, aggregate_dict["zROIStart"]]*460/100) < 20)[0],
        np.where(np.abs(b[:, aggregate_dict["vertexZ"]] - b[:, aggregate_dict["zROIStart"]]*460/100) < 20)[0]
    )),
    ("Neutrino_Tail_Length_Density", lambda s, b: (
        np.where(s[:, aggregate_dict["numberOfHitsInMuonRegion"]]/(s[:, aggregate_dict["zROIStart"]]*460/100) < 0.2)[0],
        np.where(b[:, aggregate_dict["numberOfHitsInMuonRegion"]]/(b[:, aggregate_dict["zROIStart"]]*460/100) < 0.2)[0]
    ))
]

cuts_with_mc = [
    ("vertex_fiducial_volume", lambda s, b, mc: (
        np.where(
            (s[:, aggregate_dict["vertexZ"]] >= 20) &
            (s[:, aggregate_dict["vertexY"]] <= 550)
        )[0],
        np.where(
            (b[:, aggregate_dict["vertexZ"]] >= 20) &
            (b[:, aggregate_dict["vertexY"]] <= 550)
        )[0],
        np.where(
            (mc[:, aggregate_dict["vertexZ"]] >= 20) &
            (mc[:, aggregate_dict["vertexY"]] <= 550)
        )[0]
    )),
    ("daughter_particles", lambda s, b, mc: (
        np.where(s[:, aggregate_dict["numberOfPFParticles"]] >= 6)[0],
        np.where(b[:, aggregate_dict["numberOfPFParticles"]] >= 6)[0],
        np.where(mc[:, aggregate_dict["numberOfPFParticles"]] >= 6)[0]
    )),
    ("max_directionZ", lambda s, b, mc: (
        np.where(np.maximum(
            s[:, aggregate_dict["directionZ"]],
            s[:, aggregate_dict["directionZ2"]]
        ) > 0.97)[0],
        np.where(np.maximum(
            b[:, aggregate_dict["directionZ"]],
            b[:, aggregate_dict["directionZ2"]]
        ) > 0.97)[0],
        np.where(np.maximum(
            mc[:, aggregate_dict["directionZ"]],
            mc[:, aggregate_dict["directionZ2"]]
        ) > 0.97)[0]
    )),
    ("energy_first10cm", lambda s, b, mc: (
        np.where(
            (s[:, aggregate_dict["energyDepositedInFirst10cm"]] > 3000))[0],
        np.where(
            (b[:, aggregate_dict["energyDepositedInFirst10cm"]] > 3000))[0],
        np.where(
            (mc[:, aggregate_dict["energyDepositedInFirst10cm"]] > 3000))[0]
    )),
    ("energy_fifth10cm", lambda s, b, mc: (
        np.where(
            (s[:, aggregate_dict["energyDepositedInFifth10cm"]] > 10000))[0],
        np.where(
            (b[:, aggregate_dict["energyDepositedInFifth10cm"]] > 10000))[0],
        np.where(
            (mc[:, aggregate_dict["energyDepositedInFifth10cm"]] > 10000))[0]
    )),
    ("energy_fifteenth10cm", lambda s, b, mc: (
        np.where(
            (s[:, aggregate_dict["energyDepositedInFifteenth10cm"]] > 3000))[0],
        np.where(
            (b[:, aggregate_dict["energyDepositedInFifteenth10cm"]] > 3000))[0],
        np.where(
            (mc[:, aggregate_dict["energyDepositedInFifteenth10cm"]] > 3000))[0]
    )),
    ("ROI_Z_size", lambda s, b, mc: (
        np.where(s[:, aggregate_dict["zROIEnd"]] - s[:, aggregate_dict["zROIStart"]] > 15)[0],
        np.where(b[:, aggregate_dict["zROIEnd"]] - b[:, aggregate_dict["zROIStart"]] > 15)[0],
        np.where(mc[:, aggregate_dict["zROIEnd"]] - mc[:, aggregate_dict["zROIStart"]] > 15)[0]
    )),
    ("ROI_Z_starting_point_close_to_vertexZ", lambda s, b, mc: (
        np.where(np.abs(s[:, aggregate_dict["vertexZ"]] - s[:, aggregate_dict["zROIStart"]]*460/100) < 20)[0],
        np.where(np.abs(b[:, aggregate_dict["vertexZ"]] - b[:, aggregate_dict["zROIStart"]]*460/100) < 20)[0],
        np.where(np.abs(mc[:, aggregate_dict["vertexZ"]] - mc[:, aggregate_dict["zROIStart"]]*460/100) < 20)[0]
    )),
    ("Neutrino_Tail_Length_Density", lambda s, b, mc: (
        np.where(s[:, aggregate_dict["numberOfHitsInMuonRegion"]]/(s[:, aggregate_dict["zROIStart"]]*460/100) < 0.2)[0],
        np.where(b[:, aggregate_dict["numberOfHitsInMuonRegion"]]/(b[:, aggregate_dict["zROIStart"]]*460/100) < 0.2)[0],
        np.where(mc[:, aggregate_dict["numberOfHitsInMuonRegion"]]/(mc[:, aggregate_dict["zROIStart"]]*460/100) < 0.2)[0]
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


