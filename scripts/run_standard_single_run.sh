#!bin/bash

# json_file=../analysis_settings/single_run/new_29424.json
# json_file=../analysis_settings/single_run/first_100_31036.json
# json_file=../analysis_settings/single_run/old_wnp04.json
# json_file=../analysis_settings/single_run/hits_test_120k.json
json_file=../analysis_settings/single_run/old_MC.json


python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/single_run.py -j $json_file
