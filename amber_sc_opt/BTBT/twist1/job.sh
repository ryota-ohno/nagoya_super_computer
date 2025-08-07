#!/bin/bash
#PJM -L "rscunit=fx"
#PJM -L "rscgrp=fx-small"
#PJM -L "node=1"
#PJM -L "elapse=24:00:00"
#PJM -j
#PJM -S
#PJM "--norestart"

python /data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/BTBT/src/step1_8_xyz_new.py --auto-dir twist1 --monomer-name BTBT --num-nodes 5
