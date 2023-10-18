import os
import numpy as np
import pandas as pd
import subprocess
from utils import Rod, R2atom

MONOMER_LIST = ["BTBT","BBBT","benzene","naphthalene","anthracene","tetracene","pentacene","hexacene","coronene","PTA"]
############################汎用関数###########################
def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A1,A2,A3,phi=0.0,isFF=False):
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv('~/Working/nagoya_super_computer/coronene/monomer/{}1.csv'.format(monomer_name))
    atoms_array_xyzR=df_mono[['X','Y','Z','R']].values
    
    ex = np.array([1.,0.,0.]); ey = np.array([0.,1.,0.]); ez = np.array([0.,0.,1.])

    xyz_array = atoms_array_xyzR[:,:3]
    xyz_array = np.matmul(xyz_array,Rod(ez,A3).T)
    xyz_array = np.matmul(xyz_array,Rod(-ex,A2).T)
    xyz_array = np.matmul(xyz_array,Rod(ey,A1).T)
    xyz_array = xyz_array + T_vec
    R_array = atoms_array_xyzR[:,3].reshape((-1,1))
    
    if monomer_name in MONOMER_LIST:
        return np.concatenate([xyz_array,R_array],axis=1)
    
    else:
        raise RuntimeError('invalid monomer_name={}'.format(monomer_name))
        
def get_xyzR_lines(xyzR_array,file_description):
    lines = [     
        '%mem=24GB\n',
        '%nproc=42\n',
        '#P TEST b3lyp/6-311G** EmpiricalDispersion=GD3 counterpoise=2\n',
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
    A1 = params_dict.get('A1',0.0); A2 = params_dict.get('A2',0.0); A3 = params_dict['theta']
    phi1 = params_dict.get('phi1',0.0); phi2 = params_dict.get('phi2',0.0)
    print(phi1, phi2)

    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A1,A2,A3, phi1)
    
    monomer_array_p1 = get_monomer_xyzR(monomer_name,0,b_,0,A1,A2,A3, phi1)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,0,-b_,0,A1,A2,A3, phi1)##1,2がb方向
    monomer_array_p3 = get_monomer_xyzR(monomer_name,a_,0,0,A1,A2,A3, phi1)##3,4がa方向
    monomer_array_p4 = get_monomer_xyzR(monomer_name,-a_,0,0,A1,A2,A3, phi1)##3,4がa方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a_/2,b_/2,0,-A1,A2,-A3, phi2)
    monomer_array_t2 = get_monomer_xyzR(monomer_name,a_/2,-b_/2,0,-A1,A2,-A3, phi2)
    monomer_array_t3 = get_monomer_xyzR(monomer_name,-a_/2,-b_/2,0,-A1,A2,-A3, phi2)
    monomer_array_t4 = get_monomer_xyzR(monomer_name,-a_/2,b_/2,0,-A1,A2,-A3, phi2)
    xyz_list=['400 \n','polyacene9 \n']##4分子のxyzファイルを作成
    monomers_array_4 = np.concatenate([monomer_array_i,monomer_array_p1,monomer_array_p3,monomer_array_p2,monomer_array_p4,monomer_array_t1,monomer_array_t2,monomer_array_t3,monomer_array_t4],axis=0)
    
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
        elif key in ['A1','A2']:#,'theta']:
            val = int(val)
        xyzfile_name += '_{}={}'.format(key,val)
    return xyzfile_name + '.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict,isInterlayer):
    a_ = params_dict['a']; b_ = params_dict['b']; c = np.array([params_dict.get('cx',0.0),params_dict.get('cy',0.0),params_dict.get('cz',0.0)])
    A1 = params_dict.get('A1',0.0); A2 = params_dict.get('A2',0.0); A3 = params_dict['theta']
    Rt = params_dict.get('Rt',0.0); Rp = params_dict.get('Rp',0.0)
    phi1 = params_dict.get('phi1',0.0); phi2 = params_dict.get('phi2',0.0)
    print(phi1, phi2)
    
    monomer_array_i = get_monomer_xyzR(monomer_name,0.01,0.01,0.01,A1,A2,A3, phi1)
    
    monomer_array_p1 = get_monomer_xyzR(monomer_name,0.01,b_+0.01,Rp+0.01,A1,A2,A3, phi1)##p1がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,a_+0.01,0+0.01,2*Rt-Rp+0.01,A1,A2,A3, phi1)##p2がa方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a_/2+0.01,b_/2+0.01,Rt+0.01,-A1,A2,-A3, phi2)
    monomer_array_t2 = get_monomer_xyzR(monomer_name,a_/2+0.01,-b_/2+0.01,Rt-Rp+0.01,-A1,A2,-A3, phi2)
    monomer_array_t3 = get_monomer_xyzR(monomer_name,-a_/2+0.01,-b_/2+0.01,-Rt+0.01,-A1,A2,-A3, phi2)
    monomer_array_t4 = get_monomer_xyzR(monomer_name,-a_/2+0.01,b_/2+0.01,-Rt+Rp+0.01,-A1,A2,-A3, phi2)
    
    dimer_array_t1 = np.concatenate([monomer_array_i,monomer_array_t1])
    dimer_array_t2 = np.concatenate([monomer_array_i,monomer_array_t2])
    dimer_array_t3 = np.concatenate([monomer_array_i,monomer_array_t3])
    dimer_array_t4 = np.concatenate([monomer_array_i,monomer_array_t4])
    dimer_array_p1 = np.concatenate([monomer_array_i,monomer_array_p1])
    dimer_array_p2 = np.concatenate([monomer_array_i,monomer_array_p2])
    
    file_description = '{}_A1={}_A2={}_A3={}'.format(monomer_name,int(A1),int(A2),round(A3,2))
    line_list_dimer_p1 = get_xyzR_lines(dimer_array_p1,file_description+'_p1')
    line_list_dimer_p2 = get_xyzR_lines(dimer_array_p2,file_description+'_p2')
    line_list_dimer_t1 = get_xyzR_lines(dimer_array_t1,file_description+'_t1')
    line_list_dimer_t2 = get_xyzR_lines(dimer_array_t2,file_description+'_t2')
    line_list_dimer_t3 = get_xyzR_lines(dimer_array_t3,file_description+'_t3')
    line_list_dimer_t4 = get_xyzR_lines(dimer_array_t4,file_description+'_t4')
    
    if monomer_name in MONOMER_LIST and not(isInterlayer):##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = ['$ RunGauss\n'] + line_list_dimer_t1 + ['\n\n--Link1--\n'] + line_list_dimer_p1 + ['\n\n--Link1--\n'] + line_list_dimer_p2 + ['\n\n\n']
    #elif monomer_name in MONOMER_LIST and isInterlayer:
        #gij_xyz_lines = ['$ RunGauss\n'] + line_list_dimer_i0 + ['\n\n--Link1--\n'] + line_list_dimer_ip1+ ['\n\n--Link1--\n'] + line_list_dimer_ip2+ ['\n\n--Link1--\n'] + line_list_dimer_it1 + ['\n\n--Link1--\n'] + line_list_dimer_it2 + ['\n\n\n']#['\n\n--Link1--\n'] + line_list_dimer_it3 + ['\n\n--Link1--\n'] + line_list_dimer_it4 + ['\n\n\n']
    elif monomer_name=='mono-C9-BTBT':##tshaped ４分子を全て計算
        gij_xyz_lines = ['$ RunGauss\n'] + line_list_dimer_p1 + ['\n\n--Link1--\n'] + line_list_dimer_t1 + ['\n\n--Link1--\n'] + line_list_dimer_t2 + ['\n\n--Link1--\n'] + line_list_dimer_t3 + ['\n\n--Link1--\n'] + line_list_dimer_t4 + ['\n\n\n']
    
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
        if key in ['a','b','Rt','Rp','cx','cy','cz','theta']:
            val = np.round(val,2)
        elif key in ['A1','A2']:#,'theta']:
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