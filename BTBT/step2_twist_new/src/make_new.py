##pbepbe+d3bjで計算
import os
import numpy as np
import pandas as pd
import subprocess
from utils import Rod, R2atom

MONOMER_LIST = ["BTBT"];MONOMER_LIST2 = ["mono-C4-BTBT","mono-C9-BTBT"]
############################汎用関数###########################
def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A2,A3,phi):
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv('~/Working/nagoya_super_computer/BTBT/step2_twist_new/monomer/{}.csv'.format(monomer_name))
    atoms_array_xyzR=df_mono[['X','Y','Z','R']].values
    
    ex = np.array([1.,0.,0.]); ey = np.array([0.,1.,0.]); ez = np.array([0.,0.,1.])

    xyz_array = atoms_array_xyzR[:,:3]
    xyz_array = np.matmul(xyz_array,Rod(-ex,A2).T)
    xyz_array = np.matmul(xyz_array,Rod(ez,A3).T)
    xyz_array = xyz_array + T_vec
    R_array = atoms_array_xyzR[:,3].reshape((-1,1))
    
    if monomer_name in MONOMER_LIST:
        return np.concatenate([xyz_array,R_array],axis=1)
    elif monomer_name in MONOMER_LIST2:
        C0_index = 5;C1_index = 23 #アルキルの根本
        C0=xyz_array[C0_index];C1=xyz_array[C1_index]
        n1=C1-C0;n1/=np.linalg.norm(n1)
        xyz_array[C1_index:] = np.matmul((xyz_array[C1_index:]-C0),Rod(n1,phi).T) + C0
        return np.concatenate([xyz_array,R_array],axis=1)
    else:
        raise RuntimeError('invalid monomer_name={}'.format(monomer_name))
        
def get_xyzR_lines(xyzR_array,file_description):
    lines = [     
        '%mem=24GB\n',
        '%nproc=48\n',
        '#P TEST pbepbe/6-311G** EmpiricalDispersion=GD3BJ counterpoise=2\n',
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
        '#PJM -L "elapse=2:00:00"\n',
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
def make_gaussview_xyz(auto_dir,monomer_name,params_dict,isInterlayer=False):
    a_ = params_dict['a']; b_ = params_dict['b']
    Rt = params_dict['Rt']; A2 = params_dict['A2']; A3 = params_dict['theta']
    phi1 = params_dict.get('phi1',0.0); phi2 = params_dict.get('phi2',0.0)
    print(phi1, phi2)
    a =np.array([a_,0,0])
    b =np.array([0,b_,0])
    
    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3, phi1)
    if a_>b_:
        monomer_array_p1 = get_monomer_xyzR(monomer_name,0,b_,2*Rt,A2,A3, phi1)
        monomer_array_p2 = get_monomer_xyzR(monomer_name,0,-b_,-2*Rt,A2,A3, phi1)
    else:
        monomer_array_p1 = get_monomer_xyzR(monomer_name,a_,0,0,A2,A3, phi1)
        monomer_array_p2 = get_monomer_xyzR(monomer_name,-a_,0,0,A2,A3, phi1)
    
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a_/2,b_/2,Rt,A2,-A3, phi2)
    monomer_array_t2 = get_monomer_xyzR(monomer_name,a_/2,-b_/2,Rt,A2,-A3, phi2)
    monomer_array_t3 = get_monomer_xyzR(monomer_name,-a_/2,-b_/2,Rt,A2,-A3, phi2)
    monomer_array_t4 = get_monomer_xyzR(monomer_name,-a_/2,b_/2,Rt,A2,-A3, phi2)

    monomers_array = np.concatenate([monomer_array_i,monomer_array_p1,monomer_array_t1,monomer_array_t2,monomer_array_p2,monomer_array_t3,monomer_array_t4],axis=0)
    
    file_description = 'Rt={}_A2={}_A3={}'.format(np.round(Rt,1),round(A2),round(A3))
    lines = get_xyzR_lines(monomers_array,file_description)
    lines.append('Tv {} {} {}\n'.format(a[0],a[1],a[2]))
    lines.append('Tv {} {} {}\n'.format(b[0],b[1],b[2]))
    #lines.append('Tv {} {} {}\n\n\n'.format(c[0],c[1],c[2]))
    
    os.makedirs(os.path.join(auto_dir,'gaussview'),exist_ok=True)
    output_path = os.path.join(
        auto_dir,
        'gaussview/{}_Rt={}_A2={}_A3={}_a={}_b={}.gjf'.format(monomer_name,round(Rt),round(A2),round(A3),np.round(a_,2),np.round(b_,2))
    )
            
    with open(output_path,'w') as f:
        f.writelines(lines)

def make_gjf_xyz(auto_dir,monomer_name,params_dict,isInterlayer):
    a_ = params_dict['a']; b_ = params_dict['b']
    Rt = params_dict['Rt']; A2 = params_dict['A2']; A3 = params_dict['theta']
    phi1 = params_dict.get('phi1',0.0); phi2 = params_dict.get('phi2',0.0)
    print(phi1, phi2)
    
    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3, phi1)
    if b_ > a_:
        monomer_array_p1 = get_monomer_xyzR(monomer_name,a_,0,0,A2,A3, phi1)
    else:
        monomer_array_p1 = get_monomer_xyzR(monomer_name,0,b_,2*Rt,A2,A3, phi1)
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a_/2,b_/2,Rt,A2,-A3, phi2)
    monomer_array_t4 = get_monomer_xyzR(monomer_name,-a_/2,b_/2,Rt,A2,-A3, phi2)
    
    dimer_array_t1 = np.concatenate([monomer_array_i,monomer_array_t1])
    dimer_array_t4 = np.concatenate([monomer_array_i,monomer_array_t4])
    dimer_array_p1 = np.concatenate([monomer_array_i,monomer_array_p1])
    
    file_description = '{}_Rt={}_A2={}_A3={}'.format(monomer_name,round(Rt,1),int(A2),round(A3,2))
    line_list_dimer_p1 = get_xyzR_lines(dimer_array_p1,file_description+'_p1')
    line_list_dimer_t1 = get_xyzR_lines(dimer_array_t1,file_description+'_t1')
    line_list_dimer_t4 = get_xyzR_lines(dimer_array_t4,file_description+'_t4')
    
    gij_xyz_lines = ['$ RunGauss\n'] + line_list_dimer_t1 + ['\n\n--Link1--\n'] + line_list_dimer_t4 + ['\n\n--Link1--\n'] + line_list_dimer_p1 + ['\n\n\n']#+ ['\n\n--Link1--\n'] + line_list_dimer_p2 + ['\n\n\n']
    
    file_name = get_file_name_from_dict(monomer_name,params_dict)
    inp_dir = os.path.join(auto_dir,'gaussian')
    gij_xyz_path = os.path.join(inp_dir,file_name)
    with open(gij_xyz_path,'w') as f:
        f.writelines(gij_xyz_lines)
    
    return file_name

def get_file_name_from_dict(monomer_name,paras_dict):
    file_name = ''
    file_name += monomer_name
    for key,val in paras_dict.items():
        if key in ['a','b','Rt','cx','cy','cz','theta']:
            val = np.round(val,2)
        elif key in ['A2']:#,'theta']:
            val = int(val)
        file_name += '_{}={}'.format(key,val)
    return file_name + '.inp'
    
def exec_gjf(auto_dir, monomer_name, params_dict,isInterlayer,isTest=True):
    inp_dir = os.path.join(auto_dir,'gaussian')
    print(params_dict)
    
    file_name = make_gjf_xyz(auto_dir, monomer_name, params_dict,isInterlayer)
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