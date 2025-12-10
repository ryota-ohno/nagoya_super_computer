import os
import subprocess
from utils import Rod, R2atom
import csv
import time

def concatenate(array_list):
    total_array=[]
    for arr in array_list:
        total_array.extend(arr)
    return total_array

def matmul(mat1, mat2):
    n = len(mat1)  # 行数
    result = [[0.0 for _ in range(3)] for _ in range(n)]
    for i in range(n):        # mat1 の各行
        for j in range(3):    # mat2 の各列
            for k in range(3):
                result[i][j] += mat1[i][k] * mat2[k][j]
    return result

def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A2,A3):
    csv_path=f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/{monomer_name}/monomer/{monomer_name}.csv';cols=['X','Y','Z','R'];atoms_array_xyzR = []
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # 指定列の値をfloatに変換してリスト化
            values = [float(row[col]) for col in cols]
            atoms_array_xyzR.append(values)
    xyz_array = [];R_array=[]
    for x,y,z,r in atoms_array_xyzR:
        xyz_array.append([x,y,z]);R_array.append(r)
    xyz_array = matmul(xyz_array,Rod([-1,0,0],A2))
    xyz_array = matmul(xyz_array,Rod([0,0,1],A3))
    xyz_array_ = []
    for x,y,z in xyz_array:
        xyz_array_.append([x+Ta,y+Tb,z+Tc])
    xyzR_array=[]
    for i in range(len(xyz_array_)):
        xyzR_array.append([xyz_array_[i][0],xyz_array_[i][1],xyz_array_[i][2],R_array[i]])
    return xyzR_array
        
line1='@<TRIPOS>MOLECULE\npentacene\n   72    80     2     0     0\nSMALL\nbcc\n\n\n@<TRIPOS>ATOM\n'
line2='@<TRIPOS>BOND\n'

bond_lines=[[1, 1, 2, 'ar'], [2, 1, 3, 'ar'], [3, 1, 7, 'ar'], [4, 2, 4, 'ar'], [5, 2, 9, 'ar'], [6, 3, 5, 'ar'], [7, 3, 11, '1'], [8, 4, 6, 'ar'],
            [9, 4, 14, '1'], [10, 5, 6, 'ar'], [11, 5, 16, '1'], [12, 6, 15, '1'], [13, 7, 8, 'ar'], [14, 7, 12, '1'], [15, 8, 10, 'ar'], [16, 8, 17, 'ar'],
            [17, 9, 10, 'ar'], [18, 9, 13, '1'], [19, 10, 19, 'ar'], [20, 17, 18, 'ar'], [21, 17, 21, '1'], [22, 18, 20, 'ar'], [23, 18, 23, 'ar'], [24, 19, 20, 'ar'],
            [25, 19, 22, '1'], [26, 20, 25, 'ar'], [27, 23, 24, 'ar'], [28, 23, 27, '1'], [29, 24, 26, 'ar'], [30, 24, 29, 'ar'], [31, 25, 26, 'ar'], [32, 25, 28, '1'],
            [33, 26, 31, 'ar'], [34, 29, 30, 'ar'], [35, 29, 33, '1'], [36, 30, 32, 'ar'], [37, 30, 34, '1'], [38, 31, 32, 'ar'], [39, 31, 36, '1'], [40, 32, 35, '1']]

line3='@<TRIPOS>SUBSTRUCTURE\n     1 RES1        1 GROUP             0 ****  ****    0  \n     2 RES2       37 GROUP             0 ****  ****    0 \n\n'

para_list=[]
with open(r'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/pentacene/monomer/pentacene_mono.mol2')as f:
    for line in f:
        #print(line)
        s=line.split()
        if len(s)==9:
            para_list.append([s[5],float(s[8])])
        if (line.find('BOND')>-1):
            break

def make_mono_mol2(xyzr_list_1):
    num_mol1=int(len(xyzr_list_1)/36)
    lines=[f'@<TRIPOS>MOLECULE\npentacene\n   {int(36*num_mol1)}    {int(40*num_mol1)}     1     0     0\nSMALL\nbcc\n\n\n@<TRIPOS>ATOM\n']
    for i in range(num_mol1):
        for j in range(36):
            x,y,z,r=xyzr_list_1[36*i+j]
            atom_type,charge=para_list[j]
            lines.append(f'  {36*i+j+1} {R2atom(r)} {x} {y} {z} {atom_type} 1 RES1 {charge}\n')
    lines.append('@<TRIPOS>BOND\n')
    for i in range(num_mol1):
        for j in range(40):
            a,b,c,type=bond_lines[j]
            a+=40*i;b+=36*i;c+=36*i
            lines.append(f'{a} {b} {c} {type} \n')
    line3=f'@<TRIPOS>SUBSTRUCTURE\n     1 RES1        1 GROUP             0 ****  ****    0  \n\n\n'
    lines.append(line3)
    return lines

def make_dimer_mol2(xyzr_list_1,xyzr_list_2):
    num_mol1=int(len(xyzr_list_1)/36);num_mol2=int(len(xyzr_list_2)/36)
    lines=[f'@<TRIPOS>MOLECULE\npentacene\n   {int(36*(num_mol1+num_mol2))}    {int(40*(num_mol1+num_mol2))}     1     0     0\nSMALL\nbcc\n\n\n@<TRIPOS>ATOM\n']
    for i in range(num_mol1):
        for j in range(36):
            x,y,z,r=xyzr_list_1[36*i+j]
            atom_type,charge=para_list[j]
            lines.append(f'  {36*i+j+1} {R2atom(r)} {x} {y} {z} {atom_type} 1 RES1 {charge}\n')
    for i in range(num_mol2):
        for j in range(36):
            x,y,z,r=xyzr_list_2[36*i+j]
            atom_type,charge=para_list[j]
            lines.append(f'  {36*(i+num_mol1)+j+1} {R2atom(r)} {x} {y} {z} {atom_type} 2 RES2 {charge}\n')
    lines.append('@<TRIPOS>BOND\n')
    for i in range(num_mol1+num_mol2):
        for j in range(40):
            a,b,c,type=bond_lines[j]
            a+=40*i;b+=36*i;c+=36*i
            lines.append(f'{a} {b} {c} {type} \n')
    line3=f'@<TRIPOS>SUBSTRUCTURE\n     1 RES1        1 GROUP             0 ****  ****    0  \n     2 RES2       {1+36*num_mol1} GROUP             0 ****  ****    0 \n\n'
    lines.append(line3)
    return lines

# 実行ファイル作成
def get_one_exe(auto_dir,file_name):
    file_basename = file_name
    lines_job=['#!/bin/bash\n','\n',
    'module load amber\n','\n',]
    lines_job.append(f'tleap -f {file_basename}_tleap.in\n')
    for i in range(1,4):
        file_basename_=file_basename+f'_{i}'
        lines_job.append(f'sander -O -i FF_calc.in -o {file_basename_}.out -p {file_basename_}.prmtop -c {file_basename_}.inpcrd -r min.rst -ref {file_basename_}.inpcrd\n')
        #lines_job.append(f'rm {file_basename_}.frcmod\n')
        lines_job.append(f'rm {file_basename_}.inpcrd\n')
        lines_job.append(f'rm {file_basename_}.prmtop\n')
    
    lines_tleap=['source /home/center/opt/aarch64/apps/amber/19.0/dat/leap/cmd/leaprc.gaff\n']
    for i in range(1,4):
        file_basename_=file_basename+f'_{i}'
        lines_tleap.append(f'MOL = loadmol2 {file_basename_}.mol2\n')
        lines_tleap.append(f'loadamberparams pentacene.frcmod\n')
        lines_tleap.append(f'saveamberparm MOL {file_basename_}.prmtop {file_basename_}.inpcrd\n')
    lines_tleap.append('quit\n\n')
        
    file_tleap = os.path.join(auto_dir,f'amber/{file_basename}_tleap.in')
    with open(file_tleap,'w')as f:
        f.writelines(lines_tleap)
    file_job = os.path.join(auto_dir,f'amber/job_{file_basename}.sh')
    with open(file_job,'w')as f:
        f.writelines(lines_job)
    file_job_base = os.path.join(auto_dir,f'amber/job_{file_basename}')
    return file_job_base,f'{file_basename}'

######################################## 特化関数 ########################################

##################gaussview##################
def make_xyzfile(monomer_name,params_dict):
    a = float(params_dict.get('a',0.0));b = float(params_dict.get('b',0.0)); z1 = float(params_dict.get('z1',0.0)); z2 = float(params_dict.get('z2',0.0))
    A2 = float(params_dict.get('A2',0.0)); A3 = float(params_dict.get('theta',0.0))
    cx = float(params_dict.get('cx',0.0));cy = float(params_dict.get('cy',0.0)); cz = float(params_dict.get('cz',0.0))
    
    monomer_array_c = get_monomer_xyzR(monomer_name,cx,cy,cz,A2,A3)
    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3)
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a,0,z1-z2,A2,A3)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,-a,0,z2-z1,A2,A3)##1,2がb方向
    monomer_array_p3 = get_monomer_xyzR(monomer_name,0,b,z1+z2,A2,A3)##1,2がb方向
    monomer_array_p4 = get_monomer_xyzR(monomer_name,0,-b,-z1-z2,A2,A3)##1,2がb方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a/2,b/2,z1,A2,-A3)##1,2がb方向
    monomer_array_t2 = get_monomer_xyzR(monomer_name,-a/2,b/2,z2,A2,-A3)##1,2がb方向
    monomer_array_t3 = get_monomer_xyzR(monomer_name,a/2,-b/2,-z2,A2,-A3)##1,2がb方向
    monomer_array_t4 = get_monomer_xyzR(monomer_name,-a/2,-b/2,-z1,A2,-A3)##1,2がb方向
    
    xyz_list=['400 \n','polyacene9 \n']##4分子のxyzファイルを作成
    
    monomers_array_4 = concatenate([monomer_array_c,monomer_array_i,monomer_array_p1,monomer_array_p2,monomer_array_p3,monomer_array_p4,
                                    monomer_array_t1,monomer_array_t2,monomer_array_t3,monomer_array_t4])
    
    for x,y,z,R in monomers_array_4:
        atom = R2atom(R)
        line = '{} {} {} {}\n'.format(atom,x,y,z)     
        xyz_list.append(line)
    
    return xyz_list

def make_xyz(monomer_name,params_dict):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    for key,val in params_dict.items():
        val=float(val)
        if key in ['a','b','z1','z2','cx','cy','cz']:
            val = round(val,1)
        elif key in ['A1','A2','theta']:
            val = int(val)
        xyzfile_name += '_{}_{}'.format(key,val)
    return xyzfile_name + '.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict):
    a = float(params_dict.get('a',0.0));b = float(params_dict.get('b',0.0)); z1 = float(params_dict.get('z1',0.0)); z2 = float(params_dict.get('z2',0.0))
    A2 = float(params_dict.get('A2',0.0)); A3 = float(params_dict.get('theta',0.0))
    cx = float(params_dict.get('cx',0.0));cy = float(params_dict.get('cy',0.0)); cz = float(params_dict.get('cz',0.0))
    
    monomer_array_c = get_monomer_xyzR(monomer_name,cx,cy,cz,A2,A3)
    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3)
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a,0,z1-z2,A2,A3)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,-a,0,z2-z1,A2,A3)##1,2がb方向
    monomer_array_p3 = get_monomer_xyzR(monomer_name,0,b,z1+z2,A2,A3)##1,2がb方向
    monomer_array_p4 = get_monomer_xyzR(monomer_name,0,-b,-z1-z2,A2,A3)##1,2がb方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a/2,b/2,z1,A2,-A3)##1,2がb方向
    monomer_array_t2 = get_monomer_xyzR(monomer_name,-a/2,b/2,z2,A2,-A3)##1,2がb方向
    monomer_array_t3 = get_monomer_xyzR(monomer_name,a/2,-b/2,-z2,A2,-A3)##1,2がb方向
    monomer_array_t4 = get_monomer_xyzR(monomer_name,-a/2,-b/2,-z1,A2,-A3)##1,2がb方向

    monomer_array_c_ = get_monomer_xyzR(monomer_name,cx,cy,cz,A2,-A3)
    monomer_array_i_ = get_monomer_xyzR(monomer_name,0,0,0,A2,-A3)
    monomer_array_p1_ = get_monomer_xyzR(monomer_name,a,0,z1-z2,A2,-A3)##1,2がb方向
    monomer_array_p2_ = get_monomer_xyzR(monomer_name,-a,0,z2-z1,A2,-A3)##1,2がb方向
    monomer_array_p3_ = get_monomer_xyzR(monomer_name,0,b,z1+z2,A2,-A3)##1,2がb方向
    monomer_array_p4_ = get_monomer_xyzR(monomer_name,0,-b,-z1-z2,A2,-A3)##1,2がb方向
    
    layer_1=concatenate([monomer_array_i,monomer_array_p1,monomer_array_p2,monomer_array_p3,monomer_array_p4])
    layer_2=concatenate([monomer_array_t1,monomer_array_t2,monomer_array_t3,monomer_array_t4])
    layer_3=concatenate([monomer_array_i_,monomer_array_p1_,monomer_array_p2_,monomer_array_p3_,monomer_array_p4_])

    dimer_1=make_dimer_mol2(monomer_array_c,layer_1);dimer_2=make_dimer_mol2(monomer_array_c,layer_2);dimer_3=make_dimer_mol2(monomer_array_c_,layer_3)
    dimer_list=[dimer_1,dimer_2,dimer_3]
    file_name = get_file_name_from_dict(monomer_name,params_dict)
    os.makedirs(os.path.join(auto_dir,'amber'),exist_ok=True)
    for i in range(3):
        file_mol2 = file_name+ f'_{int(i+1)}.mol2'
        os.makedirs(os.path.join(auto_dir,'amber'),exist_ok=True)
        gij_xyz_path = os.path.join(auto_dir,'amber',file_mol2)
        line=dimer_list[i]
        with open(gij_xyz_path,'w') as f:
            f.writelines(line)
        i+=1
    return file_name

def make_monomer(auto_dir,monomer_name,params_dict):
    a = float(params_dict.get('a',0.0));b = float(params_dict.get('b',0.0)); z1 = float(params_dict.get('z1',0.0)); z2 = float(params_dict.get('z2',0.0))
    A2 = float(params_dict.get('A2',0.0)); A3 = float(params_dict.get('theta',0.0))
    
    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3)
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a,0,z1-z2,A2,A3)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,-a,0,z2-z1,A2,A3)##1,2がb方向
    monomer_array_p3 = get_monomer_xyzR(monomer_name,0,b,z1+z2,A2,A3)##1,2がb方向
    monomer_array_p4 = get_monomer_xyzR(monomer_name,0,-b,-z1-z2,A2,A3)##1,2がb方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a/2,b/2,z1,A2,-A3)##1,2がb方向
    monomer_array_t2 = get_monomer_xyzR(monomer_name,-a/2,b/2,z2,A2,-A3)##1,2がb方向
    monomer_array_t3 = get_monomer_xyzR(monomer_name,a/2,-b/2,-z2,A2,-A3)##1,2がb方向
    monomer_array_t4 = get_monomer_xyzR(monomer_name,-a/2,-b/2,-z1,A2,-A3)##1,2がb方向

    monomer_array_i_ = get_monomer_xyzR(monomer_name,0,0,0,A2,-A3)
    monomer_array_p1_ = get_monomer_xyzR(monomer_name,a,0,z1-z2,A2,-A3)##1,2がb方向
    monomer_array_p2_ = get_monomer_xyzR(monomer_name,-a,0,z2-z1,A2,-A3)##1,2がb方向
    monomer_array_p3_ = get_monomer_xyzR(monomer_name,0,b,z1+z2,A2,-A3)##1,2がb方向
    monomer_array_p4_ = get_monomer_xyzR(monomer_name,0,-b,-z1-z2,A2,-A3)##1,2がb方向
    
    layer_1=concatenate([monomer_array_i,monomer_array_p1,monomer_array_p2,monomer_array_p3,monomer_array_p4])
    layer_2=concatenate([monomer_array_t1,monomer_array_t2,monomer_array_t3,monomer_array_t4])
    layer_3=concatenate([monomer_array_i_,monomer_array_p1_,monomer_array_p2_,monomer_array_p3_,monomer_array_p4_])

    mono_1=make_mono_mol2(layer_1);mono_2=make_dimer_mol2(layer_2);mono_3=make_dimer_mol2(layer_3)
    dimer_list=[mono_1,mono_2,mono_3]
    file_name = 'layer'
    os.makedirs(os.path.join(auto_dir,'amber'),exist_ok=True)
    for i in range(3):
        file_mol2 = file_name+ f'_{int(i+1)}.mol2'
        os.makedirs(os.path.join(auto_dir,'amber'),exist_ok=True)
        gij_xyz_path = os.path.join(auto_dir,'amber',file_mol2)
        line=dimer_list[i]
        with open(gij_xyz_path,'w') as f:
            f.writelines(line)
        i+=1
    return file_name

def get_file_name_from_dict(monomer_name,params_dict):
    file_name = ''
    file_name += monomer_name
    for key,val in params_dict.items():
        val=float(val)
        if key in ['a','b','z1','z2','cx','cy','cz']:
            val = round(val,1)
        elif key in ['A2','theta']:
            val = int(val)
        file_name += '_{}_{}'.format(key,val)
    return file_name
    
def exec_gjf(auto_dir, monomer_name, params_dict,isTest=True):
    xyz_dir = os.path.join(auto_dir,'gaussview')
    xyzfile_name = make_xyz(monomer_name, params_dict)
    xyz_path = os.path.join(xyz_dir,xyzfile_name)
    xyz_list = make_xyzfile(monomer_name,params_dict)
    with open(xyz_path,'w') as f:
        f.writelines(xyz_list)
    
    file_name = make_gjf_xyz(auto_dir, monomer_name, params_dict)
    file_job_base,log_file_name = get_one_exe(auto_dir,file_name)
    file_job=file_job_base+'.sh'
    if not(isTest):
        subprocess.run(['chmod','+x',file_job])
        subprocess.run([file_job])
    time.sleep(0.1)
    return log_file_name
    
def exec_mono_gjf(auto_dir, monomer_name, params_dict,isTest=True):
    xyz_dir = os.path.join(auto_dir,'gaussview')
    xyzfile_name = make_xyz(monomer_name, params_dict)
    xyz_path = os.path.join(xyz_dir,xyzfile_name)
    file_name = make_monomer(auto_dir, monomer_name, params_dict)
    file_job_base,log_file_name = get_one_exe(auto_dir,file_name)
    file_job=file_job_base+'.sh'
    if not(isTest):
        subprocess.run(['chmod','+x',file_job])
        subprocess.run([file_job])
    time.sleep(0.1)
    return log_file_name
    
############################################################################################