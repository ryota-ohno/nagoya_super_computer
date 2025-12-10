##tetracene層内計算
import os
os.environ['HOME'] ='/data/group1/z40145w'
import pandas as pd
import argparse
import subprocess

def init_process(args):
    auto_dir = f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/pentacene/step3_para/{args.auto_dir}'
    df_init=pd.read_csv(os.path.join(auto_dir,'step3_init_params.csv'))
    for i in range(1):
        cx=0.0
        dir_name = f'{cx}'
        os.makedirs(os.path.join(auto_dir,f'{dir_name}'), exist_ok=True)
        df_init_=df_init[df_init['cx']==cx]
        df_init_.to_csv(os.path.join(auto_dir,f'{dir_name}/step3_init_params.csv'),index=False)
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
        f'python /data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/pentacene/step3_para/src/step3_xyz_new_new_0.py --auto-dir {args.auto_dir}/{dir_name} --monomer-name pentacene --num-nodes 5\n',
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