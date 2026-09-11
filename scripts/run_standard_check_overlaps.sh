#!bin/bash

json_file=../analysis_settings/check_overlaps/official_4_29424.json
python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/check_overlaps.py -j $json_file

json_file=../analysis_settings/check_overlaps/official_4_29425.json
python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/check_overlaps.py -j $json_file
