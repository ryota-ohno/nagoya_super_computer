import os
import numpy as np
import pandas as pd
import subprocess
from utils import Rod, R2atom

MONOMER_LIST = ['naphthalene','anthracene','tetracene','pentacene','hexacene']
############################汎用関数###########################
def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A2,A3):
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv('/data/group1/z40145w/Working/nagoya_super_computer/polyacene/monomer/{}.csv'.format(monomer_name))
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
        
def get_xyzR_lines(xyza_array,file_description):
    lines = [     
        '%mem=24GB\n',
        f'%nproc=48\n',
        '#P TEST b3lyp/6-311G** EmpiricalDispersion=GD3 counterpoise=2\n',###汎関数や基底関数系は適宜変更する
        '\n',
        file_description+'\n',
        '\n',
        '0 1 0 1 0 1\n'
    ]
    mol_len = len(xyza_array)//2
    atom_index = 0
    mol_index = 0
    for x,y,z,R in xyza_array:
        atom = R2atom(R)
        mol_index = atom_index//mol_len + 1
        line = '{}(Fragment={}) {} {} {}\n'.format(atom,mol_index,x,y,z)     
        lines.append(line)
        atom_index += 1
    return lines

# 実行ファイル作成
def get_one_exe(file_name):
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
def make_xyzfile(monomer_name,params_dict):
    x = params_dict.get('x',0.0); y = params_dict.get('y',0.0); z = params_dict.get('z',0.0)
    
    monomer_array_1 = get_monomer_xyzR(monomer_name,0,0,0,0,0)
    monomer_array_2 = get_monomer_xyzR(monomer_name,x,y,z,0,0)##1,2がb方向

    xyz_list=['400 \n','polyacene9 \n']##4分子のxyzファイルを作成
    
    monomers_array_4 = np.concatenate([monomer_array_1,monomer_array_2],axis=0)
    for x,y,z,R in monomers_array_4:
        atom = R2atom(R)
        line = '{} {} {} {}\n'.format(atom,x,y,z)     
        xyz_list.append(line)
    
    return xyz_list

def make_xyz(monomer_name,params_dict):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    for key,val in params_dict.items():
        if key in ['x','y','z']:
            val = np.round(val,2)
        elif key in ['A1','A2','theta']:
            val = int(val)
        xyzfile_name += '_{}={}'.format(key,val)
    return xyzfile_name + f'.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict):
    x = params_dict.get('x',0.0); y = params_dict.get('y',0.0); z = params_dict.get('z',0.0)
    
    monomer_array_1 = get_monomer_xyzR(monomer_name,0,0,0,0,0)
    monomer_array_2 = get_monomer_xyzR(monomer_name,x,y,z,0,0)##1,2がb方向

    dimer_array = np.concatenate([monomer_array_1,monomer_array_2])
    
    file_description = f'{monomer_name}_x={x}_x={y}_z={z}'
    line_list_dimer = get_xyzR_lines(dimer_array,file_description)
    
    gij_xyz_lines = ['$ RunGauss\n'] + line_list_dimer + ['\n\n\n']
    
    file_name = get_file_name_from_dict(monomer_name,params_dict)
    os.makedirs(os.path.join(auto_dir,'gaussian'),exist_ok=True)
    gij_xyz_path = os.path.join(auto_dir,'gaussian',file_name)
    with open(gij_xyz_path,'w') as f:
        f.writelines(gij_xyz_lines)
    
    return file_name

def get_file_name_from_dict(monomer_name,params_dict):
    file_name = ''
    file_name += monomer_name
    for key,val in params_dict.items():
        if key in ['a','b','z']:
            val = val
        elif key in ['A2','theta']:
            val = int(val)
        file_name += '_{}={}'.format(key,val)
    return file_name + '.inp'
    
def exec_gjf(auto_dir, monomer_name, params_dict,isTest=True):
    inp_dir = os.path.join(auto_dir,'gaussian')
    xyz_dir = os.path.join(auto_dir,'gaussview')
    print(params_dict)
    
    xyzfile_name = make_xyz(monomer_name, params_dict)
    xyz_path = os.path.join(xyz_dir,xyzfile_name)
    xyz_list = make_xyzfile(monomer_name,params_dict)
    with open(xyz_path,'w') as f:
        f.writelines(xyz_list)
    
    file_name = make_gjf_xyz(auto_dir, monomer_name, params_dict)
    cc_list = get_one_exe(file_name)
    sh_filename = os.path.splitext(file_name)[0]+'.r1'
    sh_path = os.path.join(inp_dir,sh_filename)
    with open(sh_path,'w') as f:
        f.writelines(cc_list)
    if not(isTest):
        subprocess.run(['pjsub',sh_path])
    log_file_name = os.path.splitext(file_name)[0]+'.log'
    return log_file_name
    
############################################################################################