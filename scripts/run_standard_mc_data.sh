#!bin/bash

json_file=../analysis_settings/mc_data/official_4_cuts_true.json
python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/mc_data_comparison.py -j $json_file

