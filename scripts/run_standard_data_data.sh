#!bin/bash

# json_file=../analysis_settings/data_data/standard.json
json_file=../analysis_settings/data_data/official_4.json
python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/data_data_comparison.py -j $json_file
