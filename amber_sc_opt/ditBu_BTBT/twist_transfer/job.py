##tetracene層内計算
import os
os.environ['HOME'] ='/data/group1/z40145w'
import pandas as pd
import argparse
import subprocess
import numpy as np

def init_process(args):
    auto_dir = f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/BTBT/{args.auto_dir}'
    df_init=pd.read_csv(os.path.join(auto_dir,'step1_init_params.csv'))
    z_list=[np.round(z,1) for z in np.linspace(-4.0,4.0,41)]  
    for z in z_list:
        dir_name = f'{z}'
        os.makedirs(os.path.join(auto_dir,f'{dir_name}'), exist_ok=True)
        df_init_=df_init[df_init['z']==z]
        df_init_.to_csv(os.path.join(auto_dir,f'{dir_name}/step1_init_params.csv'),index=False)
        
        job_lines=[
        '#!/bin/bash \n',
        '#PJM -L "rscunit=fx"\n',
        '#PJM -L "rscgrp=fx-small"\n',
        '#PJM -L "node=1"\n',
        '#PJM -L "elapse=3:00:00"\n',
        '#PJM -j\n',
        '#PJM -S\n',
        '#PJM "--norestart"\n',
        '\n',
        f'python /data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/BTBT/src/step1_8_xyz_new.py --auto-dir {args.auto_dir}/{dir_name} --monomer-name BTBT --num-nodes 5\n',
        '\n',
        '#sleep 5 \n'
            ]
        with open(os.path.join(auto_dir,f'{dir_name}/job.sh'),'w')as f:
            f.writelines(job_lines)
        subprocess.run(['pjsub',os.path.join(auto_dir,f'{dir_name}/job.sh')])

def main_process(args):
    auto_dir = f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/BTBT/{args.auto_dir}'
    df_init=pd.read_csv(os.path.join(auto_dir,'step1_init_params.csv'))
    inprogress = True
    z_list=[np.round(z,1) for z in np.linspace(-4.0,4.0,41)]  
    while inprogress:
        for z in z_list:
            dir_name = f'{z}'
            path_dir=os.path.join(auto_dir,f'{dir_name}')
            if os.path.exists(os.path.join(path_dir,'done.txt')):
                df_init.loc[df_init['z']==z, 'status'] = 'Done'
                df_init.to_csv(os.path.join(auto_dir,'step1_init_params.csv'),index=False)
        df_init_done=df_init[df_init['status']=='Done']
        if len(df_init)==len(df_init_done):
            inprogress = False
    print('Done')

def result_process(args):
    subprocess.run(['rm','*.sh.*'])
    auto_dir = f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/BTBT/{args.auto_dir}'
    df_tot=[]
    z_list=[np.round(z,1) for z in np.linspace(-4.0,4.0,41)]  
    for z in z_list:
        dir_name = f'{z}'
        path_dir=os.path.join(auto_dir,f'{dir_name}')
        df=pd.read_csv(os.path.join(path_dir,'step1.csv'))
        df_tot.append(df)
    df_=pd.concat(df_tot, ignore_index=True)
    df_.to_csv(os.path.join(auto_dir,'step1.csv'),index=False)

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
    main_process(args)
    result_process(args)
    print("----finish process----")
    