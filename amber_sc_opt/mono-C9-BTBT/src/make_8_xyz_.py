import os
import subprocess
from utils import Rod, R2atom
import csv
import math

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

def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A2,A3,phi):
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
    
    C0_index = 5;C1_index = 23
    C0=xyz_array[C0_index];C1=xyz_array[C1_index]
    nx=C1[0]-C0[0];ny=C1[1]-C0[1];nz=C1[2]-C0[2]
    n1_=math.sqrt(nx**2+ny**2+nz**2)
    n1=[nx/n1_,ny/n1_,nz/n1_]
    print(n1)
    xyz_array_=[]; i=0
    for x,y,z in xyz_array:
        if i < 23:
            xyz_array_.append([x,y,z])
        else:
            xyz_array_1=[[x-C0[0],y-C0[1],z-C0[2]]]
            xyz_array_2=matmul(xyz_array_1,Rod(n1,phi))
            xyz_array_.append([xyz_array_2[0][0]+C0[0],xyz_array_2[0][1]+C0[1],xyz_array_2[0][2]+C0[2]])
        i+=1

    xyz_array_f = []
    for x,y,z in xyz_array_:
        xyz_array_f.append([x+Ta,y+Tb,z+Tc])
    xyzR_array=[]
    for i in range(len(xyz_array_f)):
        xyzR_array.append([xyz_array_f[i][0],xyz_array_f[i][1],xyz_array_f[i][2],R_array[i]])
    return xyzR_array
        
line1='@<TRIPOS>MOLECULE\npentacene\n   102    108     2     0     0\nSMALL\nrbcc\n\n\n@<TRIPOS>ATOM\n'
line2='@<TRIPOS>BOND\n'
bond_lines=[' 1 1 2 ar\n', ' 2 1 3 ar\n', ' 3 1 7 1\n', ' 4 2 4 ar\n', ' 5 2 9 1\n', ' 6 3 5 ar\n', ' 7 3 21 1\n', ' 8 4 6 ar\n', ' 9 4 23 1\n', 
            ' 10 5 6 ar\n', ' 11 5 22 1\n', ' 12 6 24 1\n', ' 13 7 8 2\n', ' 14 7 10 1\n', ' 15 8 9 1\n', ' 16 8 12 1\n', ' 17 10 11 1\n', ' 18 11 12 ar\n',
            ' 19 11 13 ar\n', ' 20 12 15 ar\n', ' 21 13 14 ar\n', ' 22 13 17 1\n', ' 23 14 16 ar\n', ' 24 14 18 1\n', ' 25 15 16 ar\n', ' 26 15 20 1\n', ' 27 16 19 1\n',
            ' 28 24 25 1\n', ' 29 24 33 1\n', ' 30 24 34 1\n', ' 31 25 26 1\n', ' 32 25 35 1\n', ' 33 25 36 1\n', ' 34 26 27 1\n', ' 35 26 37 1\n', ' 36 26 38 1\n',
            ' 37 27 28 1\n', ' 38 27 39 1\n', ' 39 27 40 1\n', ' 40 28 29 1\n', ' 41 28 41 1\n', ' 42 28 42 1\n', ' 43 29 30 1\n', ' 44 29 43 1\n', ' 45 29 44 1\n',
            ' 46 30 31 1\n', ' 47 30 45 1\n', ' 48 30 46 1\n', ' 49 31 32 1\n', ' 50 31 47 1\n', ' 51 31 48 1\n', ' 52 32 49 1\n', ' 53 32 50 1\n', ' 54 32 51 1\n',
            ' 55 52 53 ar\n', ' 56 52 54 ar\n', ' 57 52 58 1\n', ' 58 53 55 ar\n', ' 59 53 60 1\n', ' 60 54 56 ar\n', ' 61 54 72 1\n', ' 62 55 57 ar\n', ' 63 55 74 1\n',
            ' 64 56 57 ar\n', ' 65 56 73 1\n', ' 66 57 75 1\n', ' 67 58 59 2\n', ' 68 58 61 1\n', ' 69 59 60 1\n', ' 70 59 63 1\n', ' 71 61 62 1\n', ' 72 62 63 ar\n',
            ' 73 62 64 ar\n', ' 74 63 66 ar\n', ' 75 64 65 ar\n', ' 76 64 68 1\n', ' 77 65 67 ar\n', ' 78 65 69 1\n', ' 79 66 67 ar\n', ' 80 66 71 1\n', ' 81 67 70 1\n',
            ' 82 75 76 1\n', ' 83 75 84 1\n', ' 84 75 85 1\n', ' 85 76 77 1\n', ' 86 76 86 1\n', ' 87 76 87 1\n', ' 88 77 78 1\n', ' 89 77 88 1\n', ' 90 77 89 1\n',
            ' 91 78 79 1\n', ' 92 78 90 1\n', ' 93 78 91 1\n', ' 94 79 80 1\n', ' 95 79 92 1\n', ' 96 79 93 1\n', ' 97 80 81 1\n', ' 98 80 94 1\n', ' 99 80 95 1\n',
            ' 100 81 82 1\n', ' 101 81 96 1\n', ' 102 81 97 1\n', ' 103 82 83 1\n', ' 104 82 98 1\n', ' 105 82 99 1\n', ' 106 83 100 1\n', ' 107 83 101 1\n', ' 108 83 102 1\n']
line3='@<TRIPOS>SUBSTRUCTURE\n     1 RES1        1 GROUP             0 ****  ****    0  \n     2 RES2       52 GROUP             0 ****  ****    0 \n\n'

para_list=[]
with open(r'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/mono-C9-BTBT/monomer/mono_esp.mol2')as f:
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
f'sander -O -i FF_calc.in -o {file_basename}.out -p {file_basename}.prmtop -c {file_basename}.inpcrd -r min.rst -ref {file_basename}.inpcrd\n']
    
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
    a = float(params_dict.get('a',0.0));b = float(params_dict.get('b',0.0)); z = float(params_dict.get('z',0.0))
    A2 = float(params_dict.get('A2',0.0)); A3 = float(params_dict.get('theta',0.0)); phi = float(params_dict.get('phi',0.0))

    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3,phi)
    
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a,0,0,A2,A3,phi)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,0,b,2*z,A2,A3,phi)##1,2がb方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a/2,b/2,z,A2,-A3,-phi)##1,2がb方向
    monomer_array_t2 = get_monomer_xyzR(monomer_name,-a/2,-b/2,-z,A2,-A3,-phi)##1,2がb方向
    
    xyz_list=['400 \n','polyacene9 \n']##4分子のxyzファイルを作成
    
    if structure_type == 1:##隣接8分子について対称性より3分子でエネルギー計算
        monomers_array_4 = concatenate([monomer_array_i,monomer_array_p1])
    elif structure_type == 2:##隣接8分子について対称性より3分子でエネルギー計算
        monomers_array_4 = concatenate([monomer_array_i,monomer_array_p2])
    elif structure_type == 3 :##隣接8分子について対称性より3分子でエネルギー計算
        monomers_array_4 = concatenate([monomer_array_i,monomer_array_t1])
    elif structure_type == 4 :##隣接8分子について対称性より3分子でエネルギー計算
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
        if key in ['a','b','z']:
            val = round(val,2)
        elif key in ['A1','A2','theta','phi']:
            val = int(val)
        xyzfile_name += '_{}_{}'.format(key,val)
    return xyzfile_name + f'_{structure_type}.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict,structure_type):
    a = float(params_dict.get('a',0.0));b = float(params_dict.get('b',0.0)); z = float(params_dict.get('z',0.0))
    A2 = float(params_dict.get('A2',0.0)); A3 = float(params_dict.get('theta',0.0)); phi = float(params_dict.get('phi',0.0))

    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3,phi)
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a,0,0,A2,A3,phi)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,0,b,2*z,A2,A3,phi)##1,2がb方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a/2,b/2,z,A2,-A3,-phi)##1,2がb方向
    monomer_array_t2 = get_monomer_xyzR(monomer_name,-a/2,-b/2,-z,A2,-A3,-phi)##1,2がb方向
    
    dimer_array_p1 = concatenate([monomer_array_i,monomer_array_p1]);dimer_array_p2 = concatenate([monomer_array_i,monomer_array_p2])
    dimer_array_t1 = concatenate([monomer_array_i,monomer_array_t1]);dimer_array_t2 = concatenate([monomer_array_i,monomer_array_t2])
    
    line_list_dimer_p1 = get_xyzR_lines(dimer_array_p1);line_list_dimer_p2 = get_xyzR_lines(dimer_array_p2)
    line_list_dimer_t1 = get_xyzR_lines(dimer_array_t1);line_list_dimer_t2 = get_xyzR_lines(dimer_array_t2)
    
    if structure_type == 1:##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = line_list_dimer_p1 
    elif structure_type == 2:##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = line_list_dimer_p2 
    elif structure_type == 3:##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = line_list_dimer_t1 
    elif structure_type == 4:##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = line_list_dimer_t2 
    
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
        if key in ['a','b','z']:
            val = round(val,2)
        elif key in ['A1','A2','theta','phi']:
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