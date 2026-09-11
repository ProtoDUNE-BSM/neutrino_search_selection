#!bin/bash

json_file=../analysis_settings/event_list_comparison/official_3_vs_old.json
python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/event_list_comparison.py -j $json_file

json_file=../analysis_settings/event_list_comparison/official_4_vs_old.json
python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/event_list_comparison.py -j $json_file
