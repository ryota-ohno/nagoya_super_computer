##tetracene層内計算
import os
os.environ['HOME'] ='/data/group1/z40145w'
import pandas as pd
import argparse
import subprocess

def init_process(args):
    auto_dir = f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/pentacene/{args.auto_dir}'
    df_init=pd.read_csv(os.path.join(auto_dir,'step1_init_params.csv'))
    theta_list=[5,10,15,20,25,30,35,40,45]
    for theta in theta_list:
        dir_name = f'{theta}'
        os.makedirs(os.path.join(auto_dir,f'{dir_name}'), exist_ok=True)
        df_init_=df_init[df_init['theta']==theta]
        df_init_.to_csv(os.path.join(auto_dir,f'{dir_name}/step1_init_params.csv'),index=False)
        
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
        f'python /data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/pentacene/src/step1_8_xyz_z_new.py --auto-dir {args.auto_dir}/{dir_name} --monomer-name pentacene --num-nodes 10\n',
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
    ##maxnum-machine2 がない
    args = parser.parse_args()

    print("----main process----")
    init_process(args)
    print("----finish process----")    