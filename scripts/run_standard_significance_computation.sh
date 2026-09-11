#!bin/bash

json_file=../analysis_settings/significance_computation/official_3.json

# append to PYTHONPATH so emcee (installed under extra_libs) is importable
export PYTHONPATH=$PYTHONPATH:/exp/dune/app/users/dpullia/neutrino_search_selection/extra_libs/lib/python3.13/site-packages

python3 /exp/dune/app/users/dpullia/neutrino_search_selection/apps/significance_computation.py -j $json_file
