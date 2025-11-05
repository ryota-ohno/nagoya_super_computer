#!/bin/bash
#PJM -L "rscunit=fx"
#PJM -L "rscgrp=fx-small"
#PJM -L "node=1"
#PJM -L "elapse=24:00:00"
#PJM -j
#PJM -S
#PJM "--norestart"

nohup python /data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/pentacene/src/step1_8_xyz_z.py --auto-dir test --monomer-name pentacene --num-nodes 5 &
