import os
import subprocess
from utils import Rod, R2atom
import csv

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
bond_lines=[' 1 1 2 ar\n', ' 2 1 3 ar\n', ' 3 1 7 ar\n', ' 4 2 4 ar\n', ' 5 2 9 ar\n', ' 6 3 5 ar\n', ' 7 3 11 1\n', ' 8 4 6 ar\n', 
            ' 9 4 14 1\n', ' 10 5 6 ar\n', ' 11 5 16 1\n', ' 12 6 15 1\n', ' 13 7 8 ar\n', ' 14 7 12 1\n', ' 15 8 10 ar\n', ' 16 8 17 ar\n', 
            ' 17 9 10 ar\n', ' 18 9 13 1\n', ' 19 10 19 ar\n', ' 20 17 18 ar\n', ' 21 17 21 1\n', ' 22 18 20 ar\n', ' 23 18 23 ar\n', ' 24 19 20 ar\n', 
            ' 25 19 22 1\n', ' 26 20 25 ar\n', ' 27 23 24 ar\n', ' 28 23 27 1\n', ' 29 24 26 ar\n', ' 30 24 29 ar\n', ' 31 25 26 ar\n', ' 32 25 28 1\n', 
            ' 33 26 31 ar\n', ' 34 29 30 ar\n', ' 35 29 33 1\n', ' 36 30 32 ar\n', ' 37 30 34 1\n', ' 38 31 32 ar\n', ' 39 31 36 1\n', ' 40 32 35 1\n', 
            ' 41 37 38 ar\n', ' 42 37 39 ar\n', ' 43 37 43 ar\n', ' 44 38 40 ar\n', ' 45 38 45 ar\n', ' 46 39 41 ar\n', ' 47 39 47 1\n', ' 48 40 42 ar\n', 
            ' 49 40 50 1\n', ' 50 41 42 ar\n', ' 51 41 52 1\n', ' 52 42 51 1\n', ' 53 43 44 ar\n', ' 54 43 48 1\n', ' 55 44 46 ar\n', ' 56 44 53 ar\n', 
            ' 57 45 46 ar\n', ' 58 45 49 1\n', ' 59 46 55 ar\n', ' 60 53 54 ar\n', ' 61 53 57 1\n', ' 62 54 56 ar\n', ' 63 54 59 ar\n', ' 64 55 56 ar\n', 
            ' 65 55 58 1\n', ' 66 56 61 ar\n', ' 67 59 60 ar\n', ' 68 59 63 1\n', ' 69 60 62 ar\n', ' 70 60 65 ar\n', ' 71 61 62 ar\n', ' 72 61 64 1\n', 
            ' 73 62 67 ar\n', ' 74 65 66 ar\n', ' 75 65 69 1\n', ' 76 66 68 ar\n', ' 77 66 70 1\n', ' 78 67 68 ar\n', ' 79 67 72 1\n', ' 80 68 71 1\n']
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

def get_xyzR_lines(xyzr_array):
    lines=[]
    lines.append(line1)
    mol=int(len(xyzr_array)/2)
    for i in range(mol):
        x,y,z,r=xyzr_array[i]
        atom_type,charge=para_list[i]
        lines.append(f'  {i+1} {R2atom(r)} {x} {y} {z} {atom_type} 1 RES1 {charge}\n')
    for i in range(mol):
        x,y,z,r=xyzr_array[i+mol]
        atom_type,charge=para_list[i]
        lines.append(f'  {i+1+mol} {R2atom(r)} {x} {y} {z} {atom_type} 2 RES2 {charge}\n')   
    lines.append(line2)
    for line in bond_lines:
        lines.append(line)
    lines.append(line3)
    return lines

# 実行ファイル作成
def get_one_exe(auto_dir,file_name):
    file_basename = os.path.splitext(file_name)[0]
    lines_job=[
'#!/bin/bash\n','\n',
'module load amber\n','\n',
f'parmchk2 -i {file_basename}.mol2 -f mol2 -o {file_basename}.frcmod\n',
f'tleap -f {file_basename}_tleap.in\n',
f'sander -O -i FF_calc.in -o {file_basename}.out -p {file_basename}.prmtop -c {file_basename}.inpcrd -r min.rst -ref {file_basename}.inpcrd\n',
f'rm {file_basename}.frcmod\n',
f'rm {file_basename}.inpcrd\n',
f'rm {file_basename}.prmtop\n',
]
    
    lines_tleap=['source /home/center/opt/aarch64/apps/amber/19.0/dat/leap/cmd/leaprc.gaff\n',
f'MOL = loadmol2 {file_basename}.mol2\n',
f'loadamberparams {file_basename}.frcmod\n',
f'saveamberparm MOL {file_basename}.prmtop {file_basename}.inpcrd\n',
'quit\n']
    file_job = os.path.join(auto_dir,f'amber/job_{file_basename}.sh')
    file_tleap = os.path.join(auto_dir,f'amber/{file_basename}_tleap.in')
    
    with open(file_job,'w')as f:
        f.writelines(lines_job)
    with open(file_tleap,'w')as f:
        f.writelines(lines_tleap)

    return file_job,f'{file_basename}.out'

######################################## 特化関数 ########################################

##################gaussview##################
def make_xyzfile(monomer_name,params_dict,structure_type):
    a = float(params_dict.get('a',0.0));b = float(params_dict.get('b',0.0))
    z1 = float(params_dict.get('z1',0.0)); z2 = float(params_dict.get('z2',0.0)); za = float(params_dict.get('za',0.0)); zb = float(params_dict.get('zb',0.0))
    A2 = float(params_dict.get('A2',0.0)); A3 = float(params_dict.get('theta',0.0))

    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3)
    
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a,0,za,A2,A3)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,0,b,zb,A2,A3)##1,2がb方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a/2,b/2,z1,A2,-A3)##1,2がb方向
    monomer_array_t2 = get_monomer_xyzR(monomer_name,-a/2,b/2,z2,A2,-A3)##1,2がb方向
    
    xyz_list=['400 \n','polyacene9 \n']##4分子のxyzファイルを作成
    
    if structure_type == 1:##隣接8分子について対称性より3分子でエネルギー計算
        monomers_array_4 = concatenate([monomer_array_i,monomer_array_p1])
    elif structure_type == 2:##隣接8分子について対称性より3分子でエネルギー計算
        monomers_array_4 = concatenate([monomer_array_i,monomer_array_p2])
    elif structure_type == 3:##隣接8分子について対称性より3分子でエネルギー計算
        monomers_array_4 = concatenate([monomer_array_i,monomer_array_t1])
    elif structure_type == 4:##隣接8分子について対称性より3分子でエネルギー計算
        monomers_array_4 = concatenate([monomer_array_i,monomer_array_t2])
    
    for x,y,z,R in monomers_array_4:
        atom = R2atom(R)
        line = '{} {} {} {}\n'.format(atom,x,y,z)     
        xyz_list.append(line)
    
    return xyz_list

def make_xyz(monomer_name,params_dict,structure_type):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    for key,val in params_dict.items():
        val=float(val)
        if key in ['a','b','z1','z2','za','zb']:
            val = round(val,2)
        elif key in ['A1','A2','theta']:
            val = int(val)
        xyzfile_name += '_{}_{}'.format(key,val)
    return xyzfile_name + f'_{structure_type}.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict,structure_type):
    a = float(params_dict.get('a',0.0));b = float(params_dict.get('b',0.0))
    z1 = float(params_dict.get('z1',0.0)); z2 = float(params_dict.get('z2',0.0)); za = float(params_dict.get('za',0.0)); zb = float(params_dict.get('zb',0.0))
    A2 = float(params_dict.get('A2',0.0)); A3 = float(params_dict.get('theta',0.0))

    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3)
    
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a,0,za,A2,A3)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,0,b,zb,A2,A3)##1,2がb方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a/2,b/2,z1,A2,-A3)##1,2がb方向
    monomer_array_t2 = get_monomer_xyzR(monomer_name,-a/2,b/2,z2,A2,-A3)##1,2がb方向
    
    dimer_array_p1 = concatenate([monomer_array_i,monomer_array_p1]);dimer_array_p2 = concatenate([monomer_array_i,monomer_array_p2])
    dimer_array_t1 = concatenate([monomer_array_i,monomer_array_t1]);dimer_array_t2 = concatenate([monomer_array_i,monomer_array_t2])
    
    line_list_dimer_p1 = get_xyzR_lines(dimer_array_p1);line_list_dimer_p2 = get_xyzR_lines(dimer_array_p2)
    line_list_dimer_t1 = get_xyzR_lines(dimer_array_t1);line_list_dimer_t2 = get_xyzR_lines(dimer_array_t2)
    
    if structure_type == 1:##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = line_list_dimer_p1 + ['\n\n\n']
    elif structure_type == 2:##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = line_list_dimer_p2 + ['\n\n\n']
    elif structure_type == 3:##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = line_list_dimer_t1 + ['\n\n\n']
    elif structure_type == 4:##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = line_list_dimer_t2 + ['\n\n\n']
    
    file_name = get_file_name_from_dict(monomer_name,params_dict,structure_type)
    os.makedirs(os.path.join(auto_dir,'amber'),exist_ok=True)
    gij_xyz_path = os.path.join(auto_dir,'amber',file_name)
    with open(gij_xyz_path,'w') as f:
        f.writelines(gij_xyz_lines)
    
    return file_name

def get_file_name_from_dict(monomer_name,params_dict,structure_type):
    file_name = ''
    file_name += monomer_name
    for key,val in params_dict.items():
        val=float(val)
        if key in ['a','b','z1','z2','za','zb']:
            val = round(val,2)
        elif key in ['A1','A2','theta']:
            val = int(val)
        file_name += '_{}_{}'.format(key,val)
    return file_name + f'_{structure_type}.mol2'
    
def exec_gjf(auto_dir, monomer_name, params_dict,structure_type,isTest=True):
    xyz_dir = os.path.join(auto_dir,'gaussview')
    xyzfile_name = make_xyz(monomer_name, params_dict,structure_type)
    xyz_path = os.path.join(xyz_dir,xyzfile_name)
    xyz_list = make_xyzfile(monomer_name,params_dict,structure_type)
    with open(xyz_path,'w') as f:
        f.writelines(xyz_list)
    
    file_name = make_gjf_xyz(auto_dir, monomer_name, params_dict,structure_type)
    file_job,log_file_name = get_one_exe(auto_dir,file_name)
    if not(isTest):
        subprocess.run(['chmod','+x',file_job])
        subprocess.run([file_job])
    return log_file_name
    
############################################################################################