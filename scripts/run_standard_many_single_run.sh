#!bin/bash

# ----
json_file=../analysis_settings/single_run/official_29424.json

python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/single_run.py -j $json_file

json_file=../analysis_settings/single_run/official_29425.json

python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/single_run.py -j $json_file

# json_file=../analysis_settings/single_run/official_31036.json

# python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/single_run.py -j $json_file
