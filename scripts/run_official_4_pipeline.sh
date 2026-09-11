#!/bin/bash

json_file=../analysis_settings/data_data/official_4.json
python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/data_data_comparison.py -j $json_file
json_file=../analysis_settings/data_data/official_4_high.json
python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/data_data_comparison.py -j $json_file
json_file=../analysis_settings/data_data/official_4_low.json
python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/data_data_comparison.py -j $json_file


# append to PYTHONPATH so emcee (installed under extra_libs) is importable
export PYTHONPATH=$PYTHONPATH:/exp/dune/app/users/dpullia/neutrino_search_selection/extra_libs/lib/python3.13/site-packages

json_file=../analysis_settings/significance_computation/official_4.json
python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/significance_computation.py -j $json_file

json_file=../analysis_settings/significance_computation/official_4_low.json
python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/significance_computation.py -j $json_file
