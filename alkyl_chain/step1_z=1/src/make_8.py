import os
import numpy as np
import pandas as pd
import subprocess
from utils import Rod, R2atom

MONOMER_LIST = ["C5",'C6']
############################汎用関数###########################
def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A2,A3):
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv('/data/group1/z40145w/Working/nagoya_super_computer/alkyl_chain/monomer/{}_re.csv'.format(monomer_name))
    atoms_array_xyzR=df_mono[['X','Y','Z','R']].values
    
    ex = np.array([1.,0.,0.]); ey = np.array([0.,1.,0.]); ez = np.array([0.,0.,1.])

    xyz_array = atoms_array_xyzR[:,:3]
    xyz_array = np.matmul(xyz_array,Rod(-ex,A2).T)
    xyz_array = np.matmul(xyz_array,Rod(ez,A3).T)
    xyz_array = xyz_array + T_vec
    R_array = atoms_array_xyzR[:,3].reshape((-1,1))
    
    if monomer_name in MONOMER_LIST:
        return np.concatenate([xyz_array,R_array],axis=1)
    else:
        raise RuntimeError('invalid monomer_name={}'.format(monomer_name))
        
def get_xyzR_lines(xyzR_array,file_description,machine_type):
    lines = [     
        '%mem=24GB\n',
        '%nproc=48\n',
        '#P TEST B3LYP/6-311G** EmpiricalDispersion=GD3BJ counterpoise=2\n',###汎関数や基底関数系は適宜変更する
        '\n',
        file_description+'\n',
        '\n',
        '0 1 0 1 0 1\n'
    ]
    mol_len = len(xyzR_array)//2
    atom_index = 0
    mol_index = 0
    for x,y,z,R in xyzR_array:
        atom = R2atom(R)
        mol_index = atom_index//mol_len + 1
        line = '{}(Fragment={}) {} {} {}\n'.format(atom,mol_index,x,y,z)     
        lines.append(line)
        atom_index += 1
    return lines

# 実行ファイル作成
def get_one_exe(file_name,machine_type):
    file_basename = os.path.splitext(file_name)[0]
    cc_list=[
        '#!/bin/bash \n',
        '#PJM -L "rscunit=fx"\n',
        '#PJM -L "rscgrp=fx-small"\n',
        '#PJM -L "node=1"\n',
        '#PJM -L "elapse=24:00:00"\n',
        '#PJM -j\n',
        '#PJM -S\n',
        '#PJM "--norestart"\n',
        '\n',
        'export g16root=/home/center/opt/aarch64/apps/gaussian16/c01 \n',
        'source $g16root/g16/bsd/g16.profile \n',
        '\n',
        'export GAUSS_SCRDIR=/data/tmp/Gtmp/$PJM_JOBID \n',
        'mkdir  /data/tmp/Gtmp/$PJM_JOBID \n',
        '\n',
        'g16 < {}.inp > {}.log \n'.format(file_basename,file_basename),
        '\n',
        'rm -rf  /data/tmp/Gtmp/$PJM_JOBID \n',
        '\n',
        '\n',
        '#sleep 5 \n'
#          '#sleep 500 \n'
            ]

    return cc_list

######################################## 特化関数 ########################################

##################gaussview##################
def make_xyzfile(monomer_name,params_dict,structure_type):
    a1 = params_dict.get('a1',0.0); z1 = params_dict.get('z1',0.0)
    a2 = params_dict.get('a2',0.0); b2 = params_dict.get('b2',0.0); z2 = params_dict.get('z2',0.0)
    A2 = params_dict.get('A2',0.0); A3 = params_dict.get('theta',0.0)

    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3)
    
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a1,0,z1,A2,A3)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,a2,b2,z2,A2,A3)##1,2がb方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a1+a2,b2,z1+z2,A2,A3)##1,2がb方向
    monomer_array_t2 = get_monomer_xyzR(monomer_name,a1-a2,-b2,z1-z2,A2,A3)##1,2がb方向
    
    xyz_list=['400 \n','polyacene9 \n']##4分子のxyzファイルを作成
    
    if structure_type == 1:##隣接8分子について対称性より3分子でエネルギー計算
        monomers_array_4 = np.concatenate([monomer_array_i,monomer_array_p1],axis=0)
    elif structure_type == 2:##隣接8分子について対称性より3分子でエネルギー計算
        monomers_array_4 = np.concatenate([monomer_array_i,monomer_array_p2],axis=0)
    elif structure_type == 3:##隣接8分子について対称性より3分子でエネルギー計算
        monomers_array_4 = np.concatenate([monomer_array_i,monomer_array_p1,monomer_array_p2,monomer_array_t1,monomer_array_t2],axis=0)
    
    for x,y,z,R in monomers_array_4:
        atom = R2atom(R)
        line = '{} {} {} {}\n'.format(atom,x,y,z)     
        xyz_list.append(line)
    
    return xyz_list

def make_xyz(monomer_name,params_dict,structure_type):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    for key,val in params_dict.items():
        if key in ['a1','z1','a2','b2','z2']:
            val = np.round(val,2)
        elif key in ['A1','A2','theta']:
            val = int(val)
        xyzfile_name += '_{}={}'.format(key,val)
    return xyzfile_name + f'_{structure_type}.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict,machine_type,structure_type):
    a1 = params_dict.get('a1',0.0); z1 = params_dict.get('z1',0.0)
    a2 = params_dict.get('a2',0.0); b2 = params_dict.get('b2',0.0); z2 = params_dict.get('z2',0.0)
    A2 = params_dict.get('A2',0.0); A3 = params_dict.get('theta',0.0)

    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3)
    
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a1,0,z1,A2,A3)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,a2,b2,z2,A2,A3)##1,2がb方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a1+a2,b2,z1+z2,A2,A3)##1,2がb方向
    monomer_array_t2 = get_monomer_xyzR(monomer_name,a1-a2,-b2,z1-z2,A2,A3)##1,2がb方向
    
    dimer_array_p1 = np.concatenate([monomer_array_i,monomer_array_p1])
    dimer_array_p2 = np.concatenate([monomer_array_i,monomer_array_p2])
    dimer_array_t1 = np.concatenate([monomer_array_i,monomer_array_t1])
    dimer_array_t2 = np.concatenate([monomer_array_i,monomer_array_t2])
    
    file_description = f'{monomer_name}_a1={a1}_z1={z1}_a2={a2}_b2={b2}_z2={z2}'
    line_list_dimer_p1 = get_xyzR_lines(dimer_array_p1,file_description+'_p1',machine_type)
    line_list_dimer_p2 = get_xyzR_lines(dimer_array_p2,file_description+'_p2',machine_type)
    line_list_dimer_t1 = get_xyzR_lines(dimer_array_t1,file_description+'_t1',machine_type)
    line_list_dimer_t2 = get_xyzR_lines(dimer_array_t2,file_description+'_t2',machine_type)
    
    if structure_type == 1:##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = ['$ RunGauss\n'] + line_list_dimer_p1 + ['\n\n\n']
    elif structure_type == 2:##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = ['$ RunGauss\n'] + line_list_dimer_p2 + ['\n\n\n']
    elif structure_type == 3:##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = ['$ RunGauss\n'] + line_list_dimer_t1 + ['\n\n--Link1--\n'] + line_list_dimer_t2 + ['\n\n\n']
    
    
    file_name = get_file_name_from_dict(monomer_name,params_dict,structure_type)
    os.makedirs(os.path.join(auto_dir,'gaussian'),exist_ok=True)
    gij_xyz_path = os.path.join(auto_dir,'gaussian',file_name)
    with open(gij_xyz_path,'w') as f:
        f.writelines(gij_xyz_lines)
    
    return file_name

def get_file_name_from_dict(monomer_name,params_dict,structure_type):
    file_name = ''
    file_name += monomer_name
    for key,val in params_dict.items():
        if key in ['a1','z1','a2','b2','z2']:
            val = val
        elif key in ['A1','A2','theta']:
            val = int(val)
        file_name += '_{}={}'.format(key,val)
    return file_name + f'_{structure_type}.inp'
    
def exec_gjf(auto_dir, monomer_name, params_dict, machine_type,structure_type,isTest=True):
    inp_dir = os.path.join(auto_dir,'gaussian')
    xyz_dir = os.path.join(auto_dir,'gaussview')
    print(params_dict)
    
    xyzfile_name = make_xyz(monomer_name, params_dict,structure_type)
    xyz_path = os.path.join(xyz_dir,xyzfile_name)
    xyz_list = make_xyzfile(monomer_name,params_dict,structure_type)
    with open(xyz_path,'w') as f:
        f.writelines(xyz_list)
    
    file_name = make_gjf_xyz(auto_dir, monomer_name, params_dict,machine_type,structure_type)
    cc_list = get_one_exe(file_name,machine_type)
    sh_filename = os.path.splitext(file_name)[0]+'.r1'
    sh_path = os.path.join(inp_dir,sh_filename)
    with open(sh_path,'w') as f:
        f.writelines(cc_list)
    if not(isTest):
        subprocess.run(['pjsub',sh_path])
    log_file_name = os.path.splitext(file_name)[0]+'.log'
    return log_file_name
    
############################################################################################