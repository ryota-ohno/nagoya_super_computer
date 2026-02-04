##tetracene層内計算
import os
os.environ['HOME'] ='/data/group1/z40145w'
import pandas as pd
import argparse
import subprocess
import numpy as np

def init_process(args):
    auto_dir = f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/ditBu_BTBT/{args.auto_dir}'
    monomer_name=args.monomer_name
    df_init=pd.read_csv(os.path.join(auto_dir,'step1_init_params.csv'))
    A2_list=[-10,-8,-6,-4,-2,0,2,4,6,8,10]
    for A2 in A2_list:
        dir_name = f'{A2}'
        os.makedirs(os.path.join(auto_dir,f'{dir_name}'), exist_ok=True)
        df_init_=df_init[df_init['A2']==A2]
        df_init_.to_csv(os.path.join(auto_dir,f'{dir_name}/step1_init_params.csv'),index=False)
        os.chdir(os.path.join(auto_dir,f'{dir_name}'))
        job_lines=[
        '#!/bin/bash \n',
        '#PJM -L "rscunit=fx"\n',
        '#PJM -L "rscgrp=fx-small"\n',
        '#PJM -L "node=1"\n',
        '#PJM -L "elapse=24:00:00"\n',
        '#PJM -j\n',
        '#PJM -S\n',
        '#PJM "--norestart"\n',
        '\n',
        f'python /data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/ditBu_BTBT/src/step1_8_xyz.py --auto-dir {args.auto_dir}/{dir_name} --monomer-name {monomer_name} --num-nodes 5\n',
        '\n',
        '#sleep 12 \n'
            ]
        with open(os.path.join(auto_dir,f'{dir_name}/job.sh'),'w')as f:
            f.writelines(job_lines)
        subprocess.run(['pjsub',os.path.join(auto_dir,f'{dir_name}/job.sh')])

def update_value_in_df(df,index,key,value):
    df.loc[index,key]=value
    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--isTest',action='store_true')
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes gaussian, gaussview and csv')
    parser.add_argument('--monomer-name',type=str,help='name of monomer to be calculated')
    ##maxnum-machine2 がない
    args = parser.parse_args()

    print("----main process----")
    init_process(args)
    print("----finish process----")    