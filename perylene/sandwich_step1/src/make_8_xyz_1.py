import os
import numpy as np
import pandas as pd
import subprocess
from utils import Rod, R2atom

MONOMER_LIST = ["perylene"]
############################汎用関数###########################
def get_monomer_xyzR1(monomer_name,Ga,Gb,Gc,Ta,Tb,Tc,A2,A3):
    T_vec = np.array([Ta,Tb,Tc]);G_vec = np.array([Ga,Gb,Gc])
    df_mono=pd.read_csv('~/Working/nagoya_super_computer/perylene/monomer/{}1.csv'.format(monomer_name))
    
    ex = np.array([1.,0.,0.]); ey = np.array([0.,1.,0.]); ez = np.array([0.,0.,1.])

    atoms_array_xyzR1=df_mono[['X','Y','Z','R']].values
    xyz_array1 = atoms_array_xyzR1[:,:3]
    xyz_array1 = xyz_array1 + G_vec/2
    xyz_array1 = np.matmul(xyz_array1,Rod(-ex,A2).T)
    xyz_array1 = np.matmul(xyz_array1,Rod(ez,A3).T)
    xyz_array1 = xyz_array1 + T_vec
    R_array1 = atoms_array_xyzR1[:,3].reshape((-1,1))
    xyzR_array1=np.concatenate([xyz_array1,R_array1],axis=1)
    return xyzR_array1
        
def get_monomer_xyzR2(monomer_name,Ga,Gb,Gc,Ta,Tb,Tc,A2,A3):
    T_vec = np.array([Ta,Tb,Tc]);G_vec = np.array([Ga,Gb,Gc])
    df_mono=pd.read_csv('~/Working/nagoya_super_computer/perylene/monomer/{}1.csv'.format(monomer_name))
    
    ex = np.array([1.,0.,0.]); ey = np.array([0.,1.,0.]); ez = np.array([0.,0.,1.])

    atoms_array_xyzR1=df_mono[['X','Y','Z','R']].values
    xyz_array1 = atoms_array_xyzR1[:,:3]
    xyz_array1 = xyz_array1 - G_vec/2
    xyz_array1 = np.matmul(xyz_array1,Rod(-ex,A2).T)
    xyz_array1 = np.matmul(xyz_array1,Rod(ez,A3).T)
    xyz_array1 = xyz_array1 + T_vec
    R_array1 = atoms_array_xyzR1[:,3].reshape((-1,1))
    xyzR_array1=np.concatenate([xyz_array1,R_array1],axis=1)
    return xyzR_array1
        
def get_xyzR_lines(xyzR_array,file_description):
    lines = [     
        '%mem=24GB\n',
        '%nproc=45\n',
        '#P TEST b3lyp/6-311G** EmpiricalDispersion=GD3 counterpoise=2\n',
        '# symmetry = none\n',
        '#integral=NoXCTest\n',
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
def get_one_exe(file_name):
    file_basename = os.path.splitext(file_name)[0]
    #mkdir
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
def make_xyzfile(monomer_name,params_dict,isInterlayer=False):
    a_ = params_dict['a']; b_ = params_dict['b']
    Ga = params_dict['Ga']; Gb = params_dict['Gb']; Gc = params_dict['Gc']
    A2 = params_dict.get('A2',0.0); A3 = params_dict['theta']
    
    monomer_array_i1 = get_monomer_xyzR1(monomer_name,Ga,Gb,Gc,0,0,0,A2,A3)
    monomer_array_i2 = get_monomer_xyzR2(monomer_name,Ga,Gb,Gc,0,0,0,A2,A3)
    
    monomer_array_b_1 = get_monomer_xyzR1(monomer_name,Ga,Gb,Gc,0,-b_,0,A2,A3)##1,2がb方向
    monomer_array_b_2 = get_monomer_xyzR2(monomer_name,Ga,Gb,Gc,0,b_,0,A2,A3)##1,2がb方向
    monomer_array_a_1 = get_monomer_xyzR1(monomer_name,Ga,Gb,Gc,-a_,0,0,A2,A3)##3,4がa方向
    monomer_array_a_2 = get_monomer_xyzR2(monomer_name,Ga,Gb,Gc,a_,0,0,A2,A3)##3,4がa方向
    monomer_array_t1_1 = get_monomer_xyzR1(monomer_name,-Ga,Gb,Gc,a_/2,b_/2,0,A2,-A3)
    monomer_array_t2_1 = get_monomer_xyzR1(monomer_name,-Ga,Gb,Gc,a_/2,-b_/2,0,A2,-A3)
    monomer_array_t3_1 = get_monomer_xyzR1(monomer_name,-Ga,Gb,Gc,-a_/2,-b_/2,0,A2,-A3)
    monomer_array_t4_1 = get_monomer_xyzR1(monomer_name,-Ga,Gb,Gc,-a_/2,b_/2,0,A2,-A3)
    monomer_array_t1_2 = get_monomer_xyzR2(monomer_name,-Ga,Gb,Gc,a_/2,b_/2,0,A2,-A3)
    monomer_array_t2_2 = get_monomer_xyzR2(monomer_name,-Ga,Gb,Gc,a_/2,-b_/2,0,A2,-A3)
    monomer_array_t3_2 = get_monomer_xyzR2(monomer_name,-Ga,Gb,Gc,-a_/2,-b_/2,0,A2,-A3)
    monomer_array_t4_2 = get_monomer_xyzR2(monomer_name,-Ga,Gb,Gc,-a_/2,b_/2,0,A2,-A3)
    xyz_list=['500 \n','polyacene9 \n']##4分子のxyzファイルを作成
    monomers_array_4 = np.concatenate([monomer_array_i1,monomer_array_i2,monomer_array_b_1,monomer_array_b_2,monomer_array_a_1,monomer_array_a_2,
                                       monomer_array_t1_1,monomer_array_t1_2,monomer_array_t2_1,monomer_array_t2_2,monomer_array_t3_1,monomer_array_t3_2,monomer_array_t4_1,monomer_array_t4_2],axis=0)
    
    for x,y,z,R in monomers_array_4:
        atom = R2atom(R)
        line = '{} {} {} {}\n'.format(atom,x,y,z)     
        xyz_list.append(line)
    
    return xyz_list

def make_xyz(monomer_name,params_dict):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    for key,val in params_dict.items():
        if key in ['a','b','cx','cy','cz','theta']:
            val = np.round(val,2)
        elif key in ['A2']:#,'theta']:
            val = int(val)
        xyzfile_name += '_{}={}'.format(key,val)
    return xyzfile_name + '.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict,isInterlayer):
    a_ = params_dict['a']; b_ = params_dict['b']
    Ga = params_dict['Ga']; Gb = params_dict['Gb']; Gc = params_dict['Gc']
    A2 = params_dict.get('A2',0.0); A3 = params_dict['theta']
    Rt = params_dict.get('Rt',0.0); Rp = params_dict.get('Rp',0.0)
    
    monomer_array_i1 = get_monomer_xyzR1(monomer_name,Ga,Gb,Gc,0.0,0.0,0.0,A2,A3)
    if a_ > b_:
        monomer_array_p1 = get_monomer_xyzR2(monomer_name,Ga,Gb,Gc,0.0,b_,Rp,A2,A3)##p1がb方向
    else:
        monomer_array_p1 = get_monomer_xyzR2(monomer_name,Ga,Gb,Gc,a_,0.0,2*Rt-Rp,A2,A3)
    monomer_array_t1 = get_monomer_xyzR1(monomer_name,-Ga,Gb,Gc,a_/2,b_/2,Rt,A2,-A3)
    monomer_array_t2 = get_monomer_xyzR2(monomer_name,-Ga,Gb,Gc,a_/2,b_/2,Rt,A2,-A3)
    
    dimer_array_t1 = np.concatenate([monomer_array_i1,monomer_array_t1])
    dimer_array_t2 = np.concatenate([monomer_array_i1,monomer_array_t2])
    dimer_array_p1 = np.concatenate([monomer_array_i1,monomer_array_p1])
    
    file_description = '{}_a={}_b={}_A3={}'.format(monomer_name,round(a_,1),round(b_,1),round(A3,2))
    line_list_dimer_p1 = get_xyzR_lines(dimer_array_p1,file_description+'_p1')
    line_list_dimer_t1 = get_xyzR_lines(dimer_array_t1,file_description+'_t1')
    line_list_dimer_t2 = get_xyzR_lines(dimer_array_t2,file_description+'_t2')
    
    gij_xyz_lines = ['$ RunGauss\n'] + line_list_dimer_t1 + ['\n\n--Link1--\n'] + line_list_dimer_t2 + ['\n\n--Link1--\n'] + line_list_dimer_p1 + ['\n\n\n']
    
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
        if key in ['a','b','Rt','Rp','cx','cy','cz','Ga','Gb','Gc','theta']:
            val = np.round(val,2)
        elif key in ['A2']:#,'theta']:
            val = int(val)
        file_name += '_{}={}'.format(key,val)
    return file_name + '.inp'
    
def exec_gjf(auto_dir, monomer_name, params_dict,isInterlayer,isTest=True):
    inp_dir = os.path.join(auto_dir,'gaussian')
    xyz_dir = os.path.join(auto_dir,'gaussview')
    print(params_dict)
    
    xyzfile_name = make_xyz(monomer_name, params_dict)
    xyz_path = os.path.join(xyz_dir,xyzfile_name)
    xyz_list = make_xyzfile(monomer_name,params_dict,isInterlayer=False)
    with open(xyz_path,'w') as f:
        f.writelines(xyz_list)
    
    file_name = make_gjf_xyz(auto_dir, monomer_name, params_dict, isInterlayer)
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