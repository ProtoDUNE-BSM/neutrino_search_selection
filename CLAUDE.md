# neutrino_search_selection

DUNE ProtoDUNE-SP neutrino search and selection analysis framework. Reads ROOT files produced by the LArSoft analysis module, applies selection cuts, and produces diagnostic plots and LaTeX cut-flow tables.

## Environment setup

**Always source the LCG environment before running any Python or shell scripts:**

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc15-opt/setup.sh
```

This provides Python 3, uproot, numpy, matplotlib, scipy, emcee, and ROOT.

## Repository layout

```
apps/               Entry-point scripts (run with python3 <script> -j <json>)
  single_run.py         Analyse a single data run
  mc_data_comparison.py MC vs data comparison
  data_data_comparison.py Two-dataset data/data comparison
  cut_optimization.py   MCMC-based cut optimisation (uses emcee)

python/             Shared library modules (added to sys.path by apps)
  cuts.py               Cut thresholds dict + cuts/cuts_with_mc/cuts_single lists
  libs.py               get_files(), load_from_folder(), apply_all_cuts*(), triple_hist*(), create_table_cuts*()
  plotting.py           run_all_triple_hist_single() and other plot helpers
  name_index_association.py  aggregate_dict, true_dict, full_reco_dict column-index maps

analysis_settings/  JSON parameter files, one per app mode
  single_run/standard.json
  cut_optimization/standard.json
  mc_data/standard.json
  data_data/standard.json

scripts/            Bash wrappers that invoke the apps with a chosen JSON
  run_standard_single_run.sh
  run_standard_cut_optimization.sh
  run_standard_mc_data.sh
  run_standard_data_data.sh
  run_standard_many_single_run.sh
  hadd_script.sh        Combines ROOT files into combined.root via hadd

out/                Analysis output (plots, cut tables) — not committed
```

## Data model

Events are stored as NumPy 2-D arrays (rows = events, columns = variables). Column indices are defined in `name_index_association.py`:

- `aggregate_dict` — 49-column reconstructed aggregate tree (`ana/tree_aggregate`)
- `true_dict` — 15-column truth tree (`ana/tree_truth`)
- `full_reco_dict` — 14-column per-PFP reco tree (`ana/tree_reco`)

ROOT files live on `/pnfs/dune/scratch/users/dpullia/` (FNAL dCache).

## Cut framework (`python/cuts.py`)

`cut_thresholds` dict holds the current optimised values. Three parallel cut lists exist:

| List | Signature | Used by |
|------|-----------|---------|
| `cuts` | `(sig, bkg) → (idx_s, idx_b)` | mc_data_comparison |
| `cuts_with_mc` | `(sig, bkg, mc) → (idx_s, idx_b, idx_mc)` | cut_optimization, data_data |
| `cuts_single` | `(sig,) → idx_s` | single_run |

Current cuts (in order):
1. `neutrino_pfp_in_slice` — numberOfPFParticles ≥ 0
2. `vertex_fiducial_volume` — vertexZ ≥ 24.15, vertexY ≤ 415.24
3. `daughter_particles` — numberOfPFParticles ≥ 7.6
4. `max_directionZ` — max(directionZ, directionZ2) > 0.965
5. `energy_first10cm` — energyDepositedInFirst10cm > 774 (if vertexZ < 450)
6. `energy_fifth10cm` — energyDepositedInFifth10cm > 1735 (if vertexZ < 410)
7. `energy_fifteenth10cm` — energyDepositedInFifteenth10cm > 3207 (if vertexZ < 310)
8. `ROI_Z_size` — zROIEnd − zROIStart > 23.76
9. `ROI_Z_starting_point_close_to_vertexZ` — |vertexZ − zROIStart×4.6| < 34.21
10. `Neutrino_Tail_Length_Density` — lengthOfMuonTrack / (zROIStart×4.6) < 0.43

## Running analyses
First cd into the scripts/ directory, then run the desired script. For example:

```bash
cd scripts/
source run_standard_single_run.sh
source run_standard_cut_optimization.sh
source run_standard_mc_data.sh
source run_standard_data_data.sh
source run_standard_many_single_run.sh
```


## JSON parameter file structure

Each JSON has a `"folders"` block (input ROOT dirs, output dir) and one or more mode blocks (`"run"`, `"data"`, `"mc"`) with keys:
- `SPILL_STATUS`: `"on"` | `"off"` | `"all"`
- `GET_TRUTH`: bool
- `WEIGHT_MODE`: `"1hour"` | `"POT_1hour"` | `"none"`
- `USE_COMBINED`: bool (triggers hadd if combined.root missing)
- `TP_RATE`: `"all"` | `"low"` | `"high"` | `"mc"`
- `run_parameters.spill_{on|off|all}.{total_POT, total_time}`

Cut optimisation adds a `"settings"` block: `nsteps`, `nwalkers`, `nburn` (emcee MCMC).

## Key paths

- Input data: `/pnfs/dune/scratch/users/dpullia/`
- Output: `/exp/dune/app/users/dpullia/neutrino_search_selection/out/`
- TP monitor CSVs: `/exp/dune/app/users/dpullia/unblinding/all_beam_monitor_*.csv`
