#!/bin/bash

find /pnfs/dune/persistent/users/dpullia/2026/official_4_29424/fnal/18946 -type f -exec ls -la --block-size=G {} \;
find /pnfs/dune/persistent/users/dpullia/2026/official_4_29424/fnal/18947 -type f -exec ls -la --block-size=G {} \;
find /pnfs/dune/persistent/users/dpullia/2026/official_4_29424/fnal/18948 -type f -exec ls -la --block-size=G {} \;
find /pnfs/dune/persistent/users/dpullia/2026/official_4_29425/fnal/18942 -type f -exec ls -la --block-size=G {} \;
find /pnfs/dune/persistent/users/dpullia/2026/official_4_29425/fnal/18943 -type f -exec ls -la --block-size=G {} \;
find /pnfs/dune/persistent/users/dpullia/2026/official_4_29425/fnal/18944 -type f -exec ls -la --block-size=G {} \;
find /pnfs/dune/persistent/users/dpullia/2026/official_4_29425/fnal/18945 -type f -exec ls -la --block-size=G {} \;

scp /pnfs/dune/persistent/users/dpullia/official_4_reco_data.tar.gz dapullia@lxplus.cern.ch:/eos/user/d/dapullia/dune/pau/
