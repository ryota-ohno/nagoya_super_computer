import os
import numpy as np
import pandas as pd
import subprocess
from utils import Rod, R2atom

MONOMER_LIST = ['C5','C6','C9','C3','C7']
############################汎用関数###########################
def get_monomer_xyza(monomer_name,Ta,Tb,Tc,A2,A3):
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv('/data/group1/z40145w/Working/nagoya_super_computer/alkyl_chain/monomer/{}_.csv'.format(monomer_name))
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
        '%nproc=48\n',
        '#P TEST pbepbe/6-311G** EmpiricalDispersion=GD3BJ counterpoise=2\n',###汎関数や基底関数系は適宜変更する
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
def make_xyzfile(monomer_name,params_dict):
    a = params_dict.get('a',0.0); b = params_dict.get('b',0.0); z = params_dict.get('z',0.0)
    cx = params_dict.get('cy',0.0); cy = params_dict.get('cy',0.0); cz = params_dict.get('cz',0.0)
    A3 = params_dict['theta']; A2 = 0
    
    monomer_array_c = get_monomer_xyza(monomer_name,cx,cy,cz,0,A3)
    monomer_array_i = get_monomer_xyza(monomer_name,0,0,0,0,A3)
    monomer_array_p1 = get_monomer_xyza(monomer_name,0,b,2*z,A2,A3)##1,2がb方向
    monomer_array_p2 = get_monomer_xyza(monomer_name,0,b,-2*z,A2,A3)##1,2がb方向
    monomer_array_p3 = get_monomer_xyza(monomer_name,a,0,0,A2,A3)##3,4がa方向
    monomer_array_p4 = get_monomer_xyza(monomer_name,-a,0,0,A2,A3)##3,4がa方向
    monomer_array_t1 = get_monomer_xyza(monomer_name,a/2,b/2,z,A2,-A3)
    monomer_array_t2 = get_monomer_xyza(monomer_name,a/2,-b/2,-z,A2,-A3)
    monomer_array_t3 = get_monomer_xyza(monomer_name,-a/2,-b/2,-z,A2,-A3)
    monomer_array_t4 = get_monomer_xyza(monomer_name,-a/2,b/2,z,A2,-A3)
    xyz_list=['400 \n','polyacene9 \n']##4分子のxyzファイルを作成
    monomers_array_4 = np.concatenate([monomer_array_c,monomer_array_i,monomer_array_p1,monomer_array_p3,monomer_array_p2,monomer_array_p4,monomer_array_t1,monomer_array_t2,monomer_array_t3,monomer_array_t4],axis=0)
    
    for x,y,z,R in monomers_array_4:
        atom = R2atom(R)
        line = '{} {} {} {}\n'.format(atom,x,y,z)     
        xyz_list.append(line)
    
    return xyz_list

def make_xyz(monomer_name,params_dict):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    for key,val in params_dict.items():
        if key in ['a','b','z','cx','cy','cz']:
            val = np.round(val,2)
        elif key in ['A1','A2','theta']:#]:
            val = int(val)
        xyzfile_name += '_{}={}'.format(key,val)
    return xyzfile_name + '.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict,isInterlayer):
    a = params_dict.get('a',0.0); b = params_dict.get('b',0.0); z = params_dict.get('z',0.0)
    cx = params_dict.get('cy',0.0); cy = params_dict.get('cy',0.0); cz = params_dict.get('cz',0.0)
    A3 = params_dict['theta']; A2 = 0
    
    monomer_array_c = get_monomer_xyza(monomer_name,cx,cy,cz,0,A3)
    monomer_array_i = get_monomer_xyza(monomer_name,0,0,0,0,A3)
    monomer_array_p1 = get_monomer_xyza(monomer_name,0,b,2*z,A2,A3)##1,2がb方向
    monomer_array_p2 = get_monomer_xyza(monomer_name,0,b,-2*z,A2,A3)##1,2がb方向
    monomer_array_p3 = get_monomer_xyza(monomer_name,a,0,0,A2,A3)##3,4がa方向
    monomer_array_p4 = get_monomer_xyza(monomer_name,-a,0,0,A2,A3)##3,4がa方向
    
    monomer_array_t1 = get_monomer_xyza(monomer_name,a/2,b/2,z,A2,-A3)
    monomer_array_t2 = get_monomer_xyza(monomer_name,a/2,-b/2,-z,A2,-A3)
    monomer_array_t3 = get_monomer_xyza(monomer_name,-a/2,-b/2,-z,A2,-A3)
    monomer_array_t4 = get_monomer_xyza(monomer_name,-a/2,b/2,z,A2,-A3)

    monomer_array_c_ = get_monomer_xyza(monomer_name,cx,cy,cz,0,-A3)
    monomer_array_i_ = get_monomer_xyza(monomer_name,0,0,0,0,-A3)
    monomer_array_p1_ = get_monomer_xyza(monomer_name,0,b,2*z,A2,-A3)##1,2がb方向
    monomer_array_p2_ = get_monomer_xyza(monomer_name,0,b,-2*z,A2,-A3)##1,2がb方向
    monomer_array_p3_ = get_monomer_xyza(monomer_name,a,0,0,A2,-A3)##3,4がa方向
    monomer_array_p4_ = get_monomer_xyza(monomer_name,-a,0,0,A2,-A3)##3,4がa方向
    
    dimer_array_i0 = np.concatenate([monomer_array_c,monomer_array_i])
    dimer_array_ip1 = np.concatenate([monomer_array_c,monomer_array_p1])
    dimer_array_ip2 = np.concatenate([monomer_array_c,monomer_array_p2])
    dimer_array_ip3 = np.concatenate([monomer_array_c,monomer_array_p3])
    dimer_array_ip4 = np.concatenate([monomer_array_c,monomer_array_p4])
    
    dimer_array_i0_ = np.concatenate([monomer_array_c_,monomer_array_i_])
    dimer_array_ip1_ = np.concatenate([monomer_array_c_,monomer_array_p1_])
    dimer_array_ip2_ = np.concatenate([monomer_array_c_,monomer_array_p2_])
    dimer_array_ip3_ = np.concatenate([monomer_array_c_,monomer_array_p3_])
    dimer_array_ip4_ = np.concatenate([monomer_array_c_,monomer_array_p4_])
    
    dimer_array_it1 = np.concatenate([monomer_array_c,monomer_array_t1])
    dimer_array_it2 = np.concatenate([monomer_array_c,monomer_array_t2])
    dimer_array_it3 = np.concatenate([monomer_array_c,monomer_array_t3])
    dimer_array_it4 = np.concatenate([monomer_array_c,monomer_array_t4])
    

    file_description = '{}_A2={}_A3={}'.format(monomer_name,int(A2),round(A3,2))
    lld_i0 = get_xyzR_lines(dimer_array_i0,file_description+'_i0')
    lld_ip1 = get_xyzR_lines(dimer_array_ip1,file_description+'_ip1')
    lld_ip2 = get_xyzR_lines(dimer_array_ip2,file_description+'_ip2')
    lld_ip3 = get_xyzR_lines(dimer_array_ip3,file_description+'_ip3')
    lld_ip4 = get_xyzR_lines(dimer_array_ip4,file_description+'_ip4')
    lld_i0_ = get_xyzR_lines(dimer_array_i0_,file_description+'_i0_')
    lld_ip1_ = get_xyzR_lines(dimer_array_ip1_,file_description+'_ip1_')
    lld_ip2_ = get_xyzR_lines(dimer_array_ip2_,file_description+'_ip2_')
    lld_ip3_ = get_xyzR_lines(dimer_array_ip3_,file_description+'_ip3_')
    lld_ip4_ = get_xyzR_lines(dimer_array_ip4_,file_description+'_ip4_')
    lld_it1 = get_xyzR_lines(dimer_array_it1,file_description+'_it1')
    lld_it2 = get_xyzR_lines(dimer_array_it2,file_description+'_it1')
    lld_it3 = get_xyzR_lines(dimer_array_it3,file_description+'_it1')
    lld_it4 = get_xyzR_lines(dimer_array_it4,file_description+'_it4')
    
    if monomer_name in MONOMER_LIST and not(isInterlayer):##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = ['$ RunGauss\n']  + lld_i0 + ['\n\n--Link1--\n'] + lld_ip1 + ['\n\n--Link1--\n'] + lld_ip2 + ['\n\n--Link1--\n'] + lld_ip3 + ['\n\n--Link1--\n'] + lld_ip4 + ['\n\n--Link1--\n']
        + lld_i0_ + ['\n\n--Link1--\n'] + lld_ip1_ + ['\n\n--Link1--\n'] + lld_ip2_ + ['\n\n--Link1--\n'] + lld_ip3_ + ['\n\n--Link1--\n'] + lld_ip4_ + ['\n\n--Link1--\n'] 
        + lld_it1 + ['\n\n--Link1--\n'] + lld_it2 + ['\n\n--Link1--\n'] + lld_it3 + ['\n\n--Link1--\n'] + lld_it4 + ['\n\n\n']
    
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
        if key in ['a','b','z','cx','cy','cz']:
            val = np.round(val,2)
        elif key in ['A1','A2','theta']:#]:
            val = int(val)
        file_name += '_{}={}'.format(key,val)
    return file_name + '.inp'
    
def exec_gjf(auto_dir, monomer_name, params_dict,isInterlayer,isTest=True):
    inp_dir = os.path.join(auto_dir,'gaussian')
    xyz_dir = os.path.join(auto_dir,'gaussview')
    print(params_dict)
    
    xyzfile_name = make_xyz(monomer_name, params_dict)
    xyz_path = os.path.join(xyz_dir,xyzfile_name)
    xyz_list = make_xyzfile(monomer_name,params_dict)
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