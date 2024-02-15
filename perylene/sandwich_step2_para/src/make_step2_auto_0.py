import numpy as np
from scipy.ndimage.filters import maximum_filter
import os
import numpy as np
import pandas as pd
import sys
from tqdm import tqdm
import numpy as np
from scipy import signal
import scipy.spatial.distance as distance
from utils import Rod, R2atom
import subprocess

############################汎用関数###########################
def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A1,A2,A3,phi=0.0,isFF=False):
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv('/data/group1/z40145w/Working/nagoya_super_computer/perylene/monomer/{}1.csv'.format(monomer_name))
    atoms_array_xyzR=df_mono[['X','Y','Z','R']].values
    
    ex = np.array([1.,0.,0.]); ey = np.array([0.,1.,0.]); ez = np.array([0.,0.,1.])

    xyz_array = atoms_array_xyzR[:,:3]
    xyz_array = np.matmul(xyz_array,Rod(ez,A3).T)
    xyz_array = np.matmul(xyz_array,Rod(-ex,A2).T)
    xyz_array = np.matmul(xyz_array,Rod(ey,A1).T)
    xyz_array = xyz_array + T_vec
    R_array = atoms_array_xyzR[:,3].reshape((-1,1))
    
    return np.concatenate([xyz_array,R_array],axis=1)        

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
        '#P TEST b3lyp/6-311G** Empiricaldispersion=GD3 counterpoise=2\n',
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

def make_gjf_xyz(auto_dir,monomer_name,params_dict):##計算する際のジョブファイル作成
    a_ = params_dict['a']; b_ = params_dict['b']
    Ga = params_dict['Ga']; Gb = params_dict['Gb']; Gc = params_dict['Gc']
    A2 = params_dict.get('A2',0.0); A3 = params_dict['theta']
    ##get_monomer_xyzRで分子の座標データの作成
    monomer_array_i1 = get_monomer_xyzR1(monomer_name,Ga,Gb,Gc,0,0,0,A2,A3)##対称性よりs-pの面内は片方で十分
    z_list1=[np.round(z,1) for z in np.linspace(np.round(-2.5,1),np.round(2.5,1),int(np.round(np.round(5,1)/0.1))+1)]
    file_base_name = monomer_name + '_step2_'
    file_base_name += 'a={}_b={}_theta={}_'.format(a_,b_,A3)
    os.makedirs(os.path.join(auto_dir,'gaussian'),exist_ok=True)##auto_dir\gaussianに格納　そこへのパスを指定している
    for z in z_list1:
        gij_xyz_lines1 = ['$ RunGauss\n'];gij_xyz_lines2 = ['$ RunGauss\n'];gij_xyz_lines3 = ['$ RunGauss\n']
        if a_ > b_:
            monomer_array_p1 = get_monomer_xyzR2(monomer_name,Ga,Gb,Gc,0.0,b_,z,A2,A3)##p1がb方向
        else:
            monomer_array_p1 = get_monomer_xyzR2(monomer_name,Ga,Gb,Gc,a_,0.0,z,A2,A3)
        monomer_array_t1 = get_monomer_xyzR1(monomer_name,-Ga,Gb,Gc,a_/2,b_/2,z,A2,-A3)
        monomer_array_t2 = get_monomer_xyzR2(monomer_name,-Ga,Gb,Gc,a_/2,b_/2,z,A2,-A3)##誘導体はtが等価でないから4つつくる
        dimer_array_t1 = np.concatenate([monomer_array_i1,monomer_array_t1])##2分子の座標データを結合
        dimer_array_t2 = np.concatenate([monomer_array_i1,monomer_array_t2])##2分子の座標データを結合
        dimer_array_p1 = np.concatenate([monomer_array_i1,monomer_array_p1])
        file_description = '{}_z={}_A2={}_A3={}'.format(monomer_name,round(z,1),int(A2),round(A3,2))##ファイル名の角度部分　位置情報はそれぞれ後で加える
        line_list_dimer_p1 = get_xyzR_lines(dimer_array_p1,file_description+'_p1')##2分子の計算ファイルの文章部分の作成　位置情報をファイル名に加えた
        line_list_dimer_t1 = get_xyzR_lines(dimer_array_t1,file_description+'_t1')
        line_list_dimer_t2 = get_xyzR_lines(dimer_array_t2,file_description+'_t2')
        gij_xyz_lines1 = gij_xyz_lines1 + line_list_dimer_t1 + ['\n\n\n']##p*2+t*4でいける
        gij_xyz_lines2 = gij_xyz_lines2 + line_list_dimer_t2 + ['\n\n\n']##p*2+t*4でいける
        gij_xyz_lines3 = gij_xyz_lines3 + line_list_dimer_p1 + ['\n\n\n']
        file_name1 = file_base_name;file_name2 = file_base_name;file_name3 = file_base_name
        file_name1 +='z={}_1.inp'.format(z)
        file_name2 +='z={}_2.inp'.format(z)
        file_name3 +='z={}_3.inp'.format(z)
        gij_xyz_path1 = os.path.join(auto_dir,'gaussian',file_name1)##ファイルへのパス
        gij_xyz_path2 = os.path.join(auto_dir,'gaussian',file_name2)##ファイルへのパス
        gij_xyz_path3 = os.path.join(auto_dir,'gaussian',file_name3)##ファイルへのパス
        with open(gij_xyz_path1,'w') as f1: ##.inpファイルを作成してそこに上で作った計算をまとめたものを書き込む
            f1.writelines(gij_xyz_lines1)##.inpファイルの作成完了
        with open(gij_xyz_path2,'w') as f2: ##.inpファイルを作成してそこに上で作った計算をまとめたものを書き込む
            f2.writelines(gij_xyz_lines2)##.inpファイルの作成完了
        with open(gij_xyz_path3,'w') as f3: ##.inpファイルを作成してそこに上で作った計算をまとめたものを書き込む
            f3.writelines(gij_xyz_lines3)##.inpファイルの作成完了
    return file_base_name

def get_one_exe(file_basename):
    
    #mkdir
    cc_list=[
        '#!/bin/bash \n',
         '#PJM -L "rscunit=fx"\n',
         '#PJM -L "rscgrp=fx-small"\n',
         '#PJM -L "node=1"\n',
         '#PJM -L "elapse=4:00:00"\n',
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

auto_dir='/data/group1/z40145w/Working/nagoya_super_computer/perylene/sandwich_step2_para/1'
monomer_name='perylene'
params_dict={'a':10.6,'b':11.2,'theta':55,'Ga':3.5,'Gb':-1.4,'Gc':1.0}
isTest = False

inp_dir = os.path.join(auto_dir,'gaussian')

#print(params_dict)
z_list1=[np.round(z,1) for z in np.linspace(np.round(-2.5,1),np.round(2.5,1),int(np.round(np.round(5,1)/0.1))+1)]
file_base_name = make_gjf_xyz(auto_dir, monomer_name, params_dict)
for z in z_list1:
    file_basename1 = file_base_name;file_basename2 = file_base_name;file_basename3 = file_base_name
    file_basename1 +='z={}_1'.format(z);file_basename2 +='z={}_2'.format(z);file_basename3 +='z={}_3'.format(z)
    cc_list1 = get_one_exe(file_basename1)
    sh_filename1 = file_basename1 + '.r1'
    sh_path1 = os.path.join(inp_dir,sh_filename1)
    with open(sh_path1,'w') as f:
        f.writelines(cc_list1)
    if not(isTest):
        subprocess.run(['cd',inp_dir])
        subprocess.run(['pjsub',sh_path1])
    cc_list2 = get_one_exe(file_basename2)
    sh_filename2 = file_basename2 + '.r1'
    sh_path2 = os.path.join(inp_dir,sh_filename2)
    with open(sh_path2,'w') as f:
        f.writelines(cc_list2)
    if not(isTest):
        subprocess.run(['cd',inp_dir])
        subprocess.run(['pjsub',sh_path2])
    cc_list3 = get_one_exe(file_basename3)
    sh_filename3 = file_basename3 + '.r1'
    sh_path3 = os.path.join(inp_dir,sh_filename3)
    with open(sh_path3,'w') as f:
        f.writelines(cc_list3)
    if not(isTest):
        subprocess.run(['cd',inp_dir])
        subprocess.run(['pjsub',sh_path3])
    
############################################################################################