#!/bin/bash

OUTPUT_ARCHIVE="/pnfs/dune/persistent/users/dpullia/official_4_reco_data.tar.gz"

tar -czvf "$OUTPUT_ARCHIVE" \
    /pnfs/dune/persistent/users/dpullia/2026/official_4_29424/fnal/18946 \
    /pnfs/dune/persistent/users/dpullia/2026/official_4_29424/fnal/18947 \
    /pnfs/dune/persistent/users/dpullia/2026/official_4_29424/fnal/18948 \
    /pnfs/dune/persistent/users/dpullia/2026/official_4_29425/fnal/18942 \
    /pnfs/dune/persistent/users/dpullia/2026/official_4_29425/fnal/18943 \
    /pnfs/dune/persistent/users/dpullia/2026/official_4_29425/fnal/18944 \
    /pnfs/dune/persistent/users/dpullia/2026/official_4_29425/fnal/18945
