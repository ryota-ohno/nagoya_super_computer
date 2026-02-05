#!/bin/bash
#PJM -L "rscunit=fx"
#PJM -L "rscgrp=fx-small"
#PJM -L "node=1"
#PJM -L "elapse=1:00:00"
#PJM -j
#PJM -S
#PJM "--norestart"

python /data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/ditBu_BTBT/src/step1_8_xyz_0.py --auto-dir test --monomer-name ditBu_BTBT --num-nodes 7
