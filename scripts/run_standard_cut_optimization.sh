#!bin/bash

json_file=../analysis_settings/cut_optimization/standard.json

# /exp/dune/app/users/dpullia/neutrino_search_selection/extra_libs/lib/python3.13/site-packages
# append to PYTHONPATH
export PYTHONPATH=$PYTHONPATH:/exp/dune/app/users/dpullia/neutrino_search_selection/extra_libs/lib/python3.13/site-packages

python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/cut_optimization.py -j $json_file
